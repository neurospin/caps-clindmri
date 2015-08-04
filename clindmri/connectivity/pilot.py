#! /usr/bin/env python
##########################################################################
# NSAP - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# Import
import nibabel
import os
import numpy
import matplotlib.pyplot as plt
import clindmri.plot.pvtk as pvtk
from clindmri.segmentation.freesurfer import cortex
from clindmri.tractography.pydipy import deterministic
from clindmri.tractography import loadtxt
from clindmri.extensions.freesurfer import read_cortex_surface_segmentation
from clindmri.connectivity.anatomical import connectivity_matrix


fsdir = "/volatile/imagen/dmritest/000000022453/fs"
t1file = (
    "/volatile/imagen/dmritest/000000022453/ADNI_MPRAGE/000000022453s012a1001.nii.gz")
dfile = (
    "/volatile/imagen/dmritest/000000022453/DTI/000000022453s011a1001.nii.gz")
bvecfile = "/volatile/imagen/dmritest/000000022453/DTI/000000022453s011a1001.bvec"
bvalfile = "/volatile/imagen/dmritest/000000022453/DTI/000000022453s011a1001.bval"
outdir = "/volatile/imagen/dmritest/000000022453/processed"
labelsfile = os.path.join(outdir, "labels.nii.gz")
trffile = os.path.join(outdir, "diff_to_anat.trf")
t2file = os.path.join(outdir, "diff_to_anat.nii.gz")
maskfile = os.path.join(outdir, "mask.nii.gz")

if 0:
    physical_to_index = numpy.linalg.inv(nibabel.load(t1file).get_affine())
    seg = read_cortex_surface_segmentation(fsdir, physical_to_index, None)   
    ren = pvtk.ren()
    for hemi in ["lh", "rh"]:
        surf = seg[hemi]
        ctab = [item["color"] for _, item in surf.metadata.items()]
        actor = pvtk.surface(surf.vertices, surf.triangles, surf.labels, ctab)
        pvtk.add(ren, actor)
        pvtk.show(ren)
        pvtk.clear(ren)
        actor = pvtk.surface(surf.inflated_vertices, surf.triangles,
                             surf.labels, ctab)
        pvtk.add(ren, actor)
        pvtk.show(ren)
        pvtk.clear(ren)
if 0:
    mask, label = cortex(t1file, fsdir, outdir, dfile)
if 0:
    mask = os.path.join(outdir, "cortex000000022453s011a1001-mask.nii.gz")
    tracks = deterministic(dfile, bvecfile, bvalfile, outdir, mask_file=mask,
                           order=4, nb_seeds_per_voxel=1, step=0.5, fmt="%.4f")
if 0:
    tracks = os.path.join(outdir, "000000022453s011a1001.trk")
    streamlines = loadtxt(tracks)
    ren = pvtk.ren()
    actor = pvtk.line(streamlines, 0)
    pvtk.add(ren, actor)
    pvtk.show(ren)

if 1:
    label = os.path.join(outdir, "cortex000000022453s011a1001-labels.nii.gz")
    tracks = os.path.join(outdir, "000000022453s011a1001.trk")
    connectivity = connectivity_matrix(tracks, label, symmetric=True)
    plt.imshow(numpy.log1p(connectivity), interpolation="nearest")
    plt.show()
    
