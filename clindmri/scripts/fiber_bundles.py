#! /usr/bin/env python
##########################################################################
# NSAP - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import glob
import numpy
import nibabel
import scipy
import json

# Dipy import
from dipy.segment.clustering import QuickBundles

# Clindmri import
from clindmri.tractography.pydipy import deterministic
from clindmri.tractography.gaussian_processes import FiberGP
from clindmri.tractography.utils import resample
from clindmri.tractography.utils import filter_by_length
from clindmri.clustering.agglomerative_clustering import fdist
from clindmri.clustering.agglomerative_clustering import agglomerative_clustering
from clindmri.registration.utils import extract_image
from clindmri.segmentation.fsl import bet2
from clindmri.plot.slicer import plot_image
import clindmri.plot.pvtk as pvtk
from clindmri.plot.colors import line_colors


# Global parameters
diffusion_file = "/volatile/imagen/dmritest/001/raw/hardi-b1500-1-001.nii.gz"
bvecs_file = "/volatile/imagen/dmritest/001/raw/hardi-b1500-1-001.bvec"
bvals_file = "/volatile/imagen/dmritest/001/raw/hardi-b1500-1-001.bval"
mask_file = "/volatile/imagen/dmritest/001/cortex_mask.nii.gz" #None
outdir = "/volatile/imagen/dmritest/001/processed/tGP"
use_vtk = False
nb_seeds_per_voxel = 1
display_amount = 20000
rfactor = 3
speed_factor = 1
min_length = 10
actor_ang = (-90, 0, -90)


"""
Define first a quality check folder.
"""
qcdir = os.path.join(outdir, "qc")
if not os.path.isdir(qcdir):
    os.makedirs(qcdir)


"""
Test
----
"""
fiber1 = []
deriv = 1
tensor1 = []
eigval = 0.1
for i in range(40):
    fiber1.append([20 * numpy.sin(i * numpy.pi / 40 ), i, 5])
    eigval += 0.1
fiber1= numpy.asarray(fiber1)

#fiber1 = scipy.signal.resample(fiber1, 200)

fiber2 = []
deriv = 1
for i in range(40):
    fiber2.append([ - 20 * numpy.sin(i * numpy.pi / 40 ) + 20, i, 5])
fiber2= numpy.asarray(fiber2)

if 0:
    fiber1 = numpy.array([
        [ 1.7326782,   0.25029223,  5.        ],
        [ 2.70803832,  1.71739357,  5.        ],
        [ 3.67926713,  2.43806875,  5.        ],
        [ 4.0906216,   3.27442217,  5.        ],
        [ 4.27240855,  4.4456303,   5.        ],
        [ 5.2095575,   5.5033897,   5.        ],
        [ 5.53538539,  6.29320943,  5.        ],
        [ 5.59439468,  7.36038572,  5.        ],
        [ 6.06738832,  8.2942982,   5.        ],
        [ 6.92363806,  9.03236326,  5.        ]])
    fiber2 = numpy.array([
        [ 1.1796562,   0.16609008,  5.        ],
        [ 1.19413916,  1.96607455,  5.        ],
        [ 1.85886967,  2.93353053,  5.        ],
        [ 2.79900016,  3.09145943,  5.        ],
        [ 3.26364434,  4.95455423,  5.        ],
        [ 3.46143541,  5.12321982,  5.        ],
        [ 3.70400735,  6.00718094,  5.        ],
        [ 4.64417862,  7.91226521,  5.        ],
        [ 4.86004664,  8.64782218,  5.        ],
        [ 5.01891052,  9.34209871,  5.        ]])


"""
Non-diffusion-weighted mask
---------------------------

For tractography, we need to generate a mask within which we
constrain tractography. We first select the first non-diffusion weighted
volume of the DTI sequence and then use 'bet2' on this image with a fractional
intensity threshold of 0.25 (this is generally a robust threshold to
remove unwanted tissue from a non-diffusion weighted image) and the 'm' option
that creates a binary 'nodif_brain_mask' image.
"""
if mask_file is None:

    # Extract the b0 map
    bvals = numpy.loadtxt(bvals_file).tolist()
    b0_index = bvals.index(0)
    b0_file = os.path.join(outdir, "nodif.nii.gz")
    if not os.path.isfile(b0_file):
        extract_image(diffusion_file, index=b0_index, out_file=b0_file)
    snap_file = os.path.join(qcdir, "nodif.png")
    plot_image(b0_file, snap_file=snap_file, name="nodif")

    # Extract the brain mask from the b0 map
    b0_brain_file = os.path.join(outdir, "nodif_brain")
    bet_files = glob.glob(b0_brain_file + "*")
    if len(bet_files) == 0:
        (output, mask_file, mesh_file, outline_file,
         inskull_mask_file, inskull_mesh_file,
         outskull_mask_file, outskull_mesh_file, outskin_mask_file,
         outskin_mesh_file, skull_mask_file) = bet2(
            b0_file,
            b0_brain_file,
            m=True,
            f=0.25)
    else:
        mask_file = sorted(bet_files)[0]
        if not os.path.isfile(mask_file):
            raise IOError("FileDoesNotExist: '{0}'.".format(mask_file))
    snap_file = os.path.join(qcdir, "bet.png")
    plot_image(b0_file, contour_file=mask_file, snap_file=snap_file,
               name="bet")


"""
Tractography
------------

We use dipy to generate a streamline tractography. Play with the
'nb_seeds_per_voxel' parameter to increase/decrease the fiber density. Then
we count the tracks that start and end at each label pair.
Finally, computed fibers are length filtered.
"""
track_outdir = os.path.join(outdir, "streamline")
track_file = os.path.join(track_outdir, "fibers.trk")
track_filtered_file = os.path.join(track_outdir, "filtered_fibers.trk")
if not os.path.isdir(track_outdir):
    os.mkdir(track_outdir)
if not os.path.isfile(track_file):
    trackvis_fibers, trackvis_header = deterministic(
        diffusion_file,
        bvecs_file,
        bvals_file,
        track_file,
        mask_file=mask_file,
        order=4,
        nb_seeds_per_voxel=nb_seeds_per_voxel,
        step=0.5)
else:
    trackvis_fibers, trackvis_header = nibabel.trackvis.read(
        track_file, as_generator=False, points_space="voxel")
if not os.path.isfile(track_filtered_file):
    fibers = [track_item[0] for track_item in trackvis_fibers]
    fibers = filter_by_length(fibers, min_length)
    streamlines = ((track, None, None) for track in fibers)
    nibabel.trackvis.write(track_filtered_file, streamlines, trackvis_header,
                           points_space="voxel")
else:
    trackvis_fibers, trackvis_header = nibabel.trackvis.read(
        track_filtered_file, as_generator=False, points_space="voxel")
    fibers = [track_item[0] for track_item in trackvis_fibers]

if use_vtk:
    ren = pvtk.ren()
    colors = line_colors(fibers)
    actor = pvtk.tubes(fibers, colors)
    actor.RotateX(actor_ang[0])
    actor.RotateY(actor_ang[1])
    actor.RotateZ(actor_ang[2])
    pvtk.add(ren, actor)
    ren.SetBackground(1, 1, 1)
    pvtk.record(ren, qcdir, "fibers", az_ang=45, n_frames=2)
    pvtk.record(ren, qcdir, "fibers", n_frames=36, az_ang=10, animate=True,
                delay=25)
    pvtk.show(ren)
    pvtk.clear(ren)


"""
Fiber clustering
----------------

Based on an agglomerative clustering, and a geometric distance.
"""

clustering_outdir = os.path.join(outdir, "clustering")
cluster_file = os.path.join(clustering_outdir, "clusters.json")
if not os.path.isdir(clustering_outdir):
    os.mkdir(clustering_outdir)
if not os.path.isfile(cluster_file):
    fibers_18 = [resample(track, nb_pol=18) for track in fibers]
    qb = QuickBundles(threshold=10.)
    clusters_ = qb.cluster(fibers_18)
    clusters = {}
    for cnt, cluster in enumerate(clusters_):
        clusters[str(cnt)] = {"indices": cluster.indices}
    with open(cluster_file, "w") as open_file:
        json.dump(clusters, open_file, indent=4)
else:
    with open(cluster_file) as open_file:
        clusters = json.load(open_file)

if 1: #use_vtk:
    ren = pvtk.ren()
    colors = numpy.ones((len(fibers),))
    nb_clusters = len(clusters)
    for clusterid, item in clusters.items():
        indices = item["indices"]
        colors[indices] = numpy.random.rand()
    actor = pvtk.line(fibers, colors.tolist())
    actor.RotateX(actor_ang[0])
    actor.RotateY(actor_ang[1])
    actor.RotateZ(actor_ang[2])
    pvtk.add(ren, actor)
    ren.SetBackground(1, 1, 1)
    #pvtk.record(ren, qcdir, "clusters", az_ang=45, n_frames=2)
    pvtk.record(ren, qcdir, "clusters", n_frames=36, az_ang=10, animate=True,
                delay=25)
    pvtk.show(ren)
    pvtk.clear(ren)
    print(stop)


"""
Bundle statistic
----------------

Combine the N fiber of a bundle by simply averaging the GPs corresponding to
these fibers, and obtain a GP that corresponds to the indicator function
of the fiber bundle.
"""
stats_outdir = os.path.join(outdir, "statistics")
bundle_id = "50"
mean_file = os.path.join(stats_outdir, "mean_{0}.nii.gz".format(bundle_id))
var_file = os.path.join(stats_outdir, "var_{0}.nii.gz".format(bundle_id))
if not os.path.isdir(stats_outdir):
    os.mkdir(stats_outdir)
if not (os.path.isfile(mean_file) and os.path.isfile(var_file)):
    bundle_mean_array = numpy.zeros(trackvis_header["dim"])
    bundle_var_array = numpy.zeros(trackvis_header["dim"])
    nb_fibers = len(clusters[bundle_id]["indices"])
    for index in clusters[bundle_id]["indices"]:
        track = fibers[index]
        gp = FiberGP(track, rfactor)
        bundle_mean_array += gp.get_mean_field(shape=trackvis_header["dim"])
        bundle_var_array += gp.get_variance_field(shape=trackvis_header["dim"])
    bundle_mean_array /= float(nb_fibers)
    bundle_var_array /= float(nb_fibers)**2
    nifti_image = nibabel.Nifti1Image(bundle_mean_array, numpy.eye(4))
    nibabel.save(nifti_image, mean_file)
    nifti_image = nibabel.Nifti1Image(bundle_var_array, numpy.eye(4))
    nibabel.save(nifti_image, var_file)

    if 1: #use_vtk:
        bundle_fibers = [fibers[index] 
                         for index in clusters[bundle_id]["indices"]]
        ren = pvtk.ren()
        actor = pvtk.line(bundle_fibers, 1)
        pvtk.add(ren, actor)
        ren.SetBackground(1, 1, 1)
        pvtk.record(ren, qcdir, "bundle_{0}".format(bundle_id), az_ang=45,
                    n_frames=2)
        pvtk.show(ren)
        pvtk.clear(ren)


