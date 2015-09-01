#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" A script for cortical connectivity using FreeSurfer Labels and Dipy
determinist tractography.

We will describe a method for performing deterministic tractography in a single
subject using automatically generated cortical labels from Freesurfer resulting
in a deterministic region connectivity matrix.

Prerequisites
=============

This analysis requires the following MR data:

* A T1 image: this should cover the entire brain at high-resolution
  (~1mm3 voxels).
* A DTI sequence with at least 30 directions, though depending on your SnR
  fewer directions can be acceptable in certain cases. This sequence has
  to be preprocessed in order to deal with suceptibility artifacts, spikes
  and Eddy current distortions.
* The b-values & b-vectors for the gradients: these files should be generated
  in the process from converting your raw images (probably in DICOM format) to
  the standard research format (NIFTI).

Software
========

This analysis is using Freesurfer & FSL & nilearn & dipy, all of them beeing
freely available.
"""

from clindmri import signaturedecorator
import os
import shutil
import numpy
import nibabel
import glob
from clindmri.tractography import loadtxt
from clindmri.extensions.freesurfer import read_cortex_surface_segmentation
from clindmri.estimation.fsl import dtifit
from clindmri.segmentation.freesurfer import cortex
from clindmri.registration.utils import extract_image
from clindmri.segmentation.fsl import bet2
from clindmri.tractography.pydipy import deterministic
from clindmri.connectivity.anatomical import diffusion_connectivity_matrix
from clindmri.connectivity.anatomical import anatomical_connectivity_matrix
from clindmri.plot.slicer import plot_image
from clindmri.plot.slicer import  plot_matrix
import clindmri.plot.pvtk as pvtk


t1_file = "/volatile/imagen/dmritest/001/raw/t1-1mm-1-001.nii.gz"
fsdir = "/volatile/imagen/dmritest/001/processed/fs"
diffusion_file = "/volatile/imagen/dmritest/001/raw/hardi-b1500-1-001.nii.gz"
bvecs_file = "/volatile/imagen/dmritest/001/raw/hardi-b1500-1-001.bvec"
bvals_file = "/volatile/imagen/dmritest/001/raw/hardi-b1500-1-001.bval"
outdir = "/volatile/imagen/dmritest/001/processed/detconnect-t1"
use_t1_space = True
nb_seeds_per_voxel = 1
use_vtk = True

"""
Define first a quality check folder.
"""

qcdir = os.path.join(outdir, "qc")
if not os.path.isdir(qcdir):
    os.makedirs(qcdir)

"""
Structural processing
=====================

The T1 image needs to be processed through Freesurfer's standard recon-all
pipeline. There are many resources for how to do this online, namely the
Freesurfer wiki. The proposed analysis run the pipeline this way:

* The 'recon-all' imports the data and creates the standard folder layout in 
  $SUBJECTS_DIR/$SUBJECT_NAME where the SUBJECT_NAME is passed through the '-s'
  option.
* Then we execute all three steps of the Freesurfer pipeline (with the '-all'
  flag).

This process usually takes between 20-40 hours depending on the quality of
data.
"""

if fsdir is None:
    raise NotImplementedError()
    # recon-all -s 0001 -i T1.nii.gz

if use_vtk:
    physical_to_index = numpy.linalg.inv(nibabel.load(t1_file).get_affine())
    hemi_surfaces = read_cortex_surface_segmentation(fsdir, physical_to_index)   
    ren = pvtk.ren()
    for hemi in ["lh", "rh"]:
        surface = hemi_surfaces[hemi]
        ctab = [item["color"] for _, item in surface.metadata.items()]
        actor = pvtk.surface(surface.vertices, surface.triangles, surface.labels,
                             ctab)
        pvtk.add(ren, actor)
        pvtk.record(ren, qcdir, hemi + "_white")
        pvtk.clear(ren)
        actor = pvtk.surface(surface.inflated_vertices, surface.triangles,
                             surface.labels, ctab)
        pvtk.add(ren, actor)
        pvtk.record(ren, qcdir, hemi + "_inflated")
        pvtk.clear(ren)

"""
Diffusion Processing
====================

At this point we have a motion- & artifact-corrected image, the corrected
gradient table and a mask of the non-diffusion-weighted image.

From our DTI data, we need to produce the following information:

* The mask of the non-diffusion-weighted image.
* A Fractional Anisotropy (FA) image for registration to the T1.
* A streamline determinist traxctography representing a model of the white
  matter organization.

Non-diffusion-weighted mask
---------------------------

For probabalistic tractography, we need to generate a mask within which we
constrain tractography. We first select the first non-diffusion weighted
volume of the DTI sequence and then use 'bet2' on this image with a fractional
intensity threshold of 0.25 (this is generally a robust threshold to
remove unwanted tissue from a non-diffusion weighted image) and the 'm' option
that creates a binary 'nodif_brain_mask' image.
"""

bvals = numpy.loadtxt(bvals_file).tolist()
b0_index = bvals.index(0)
b0_file = os.path.join(outdir, "nodif.nii.gz")
if not os.path.isfile(b0_file):
    extract_image(
        diffusion_file,
        index=b0_index,
        out_file=b0_file)

snap_file = os.path.join(qcdir, "nodif.pdf")
plot_image(b0_file, snap_file=snap_file, name="nodif")


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

snap_file = os.path.join(qcdir, "bet.pdf")
plot_image(b0_file, contour_file=mask_file, snap_file=snap_file, name="bet")

"""
Generating the FA image
-----------------------

There are much better ways to create an FA image than this method, but we're
only using this image for registration purposes. It is probably not a good
idear to use this image in a whole-brain FA analysis.
"""

dtifit_outdir = os.path.join(outdir, "dtifit")
if not os.path.isdir(dtifit_outdir):
    os.mkdir(dtifit_outdir)
if len(os.listdir(dtifit_outdir)) == 0:
    (v1_file, v2_file, v3_file, l1_file,
     l2_file, l3_file, md_file, fa_file, 
     s0_file, tensor_file, m0_file) = dtifit(
        diffusion_file, 
        bvecs_file,
        bvals_file,
        mask_file,
        dtifit_outdir)
else:
    fa_file = glob.glob(os.path.join(dtifit_outdir, "dtifit_FA.*"))[0]
    if not os.path.isfile(fa_file):
        raise IOError("FileDoesNotExist: '{0}'.".format(fa_file))

snap_file = os.path.join(qcdir, "fa.pdf")
plot_image(fa_file, snap_file=snap_file, name="fa")

"""
Registration of the FA image to the T1 image
--------------------------------------------

Register the FA image to the T1 image and generate a 'cortex_whitematter-mask'
and 'cortex_gyri-labels' in the FA space. It also generate seed regions
to compute the connectivity matrix in the FA space in the 'gyri' folder.
"""

cortex_outdir = os.path.join(outdir, "cortex")
if not os.path.isdir(cortex_outdir):
    os.mkdir(cortex_outdir)
if len(os.listdir(cortex_outdir)) == 0:
    mask_file, _, _, reg_file, trf_file = cortex(
        t1_file,
        fsdir,
        cortex_outdir,
        dest_file=fa_file,
        generate_mask=True,
        generate_seeds=False)
    if use_t1_space:
        _, label_file, _, _, _ = cortex(
            t1_file,
            fsdir,
            cortex_outdir,
            generate_mask=False,
            generate_seeds=False)

else:
    reg_file = os.path.join(cortex_outdir, "cortex_dest_to_t1.nii.gz")
    trf_file = os.path.join(cortex_outdir, "cortex_dest_to_t1.trf")
    mask_file = os.path.join(cortex_outdir, "cortex_mask.nii.gz")
    label_file = os.path.join(cortex_outdir, "cortex_gyri_labels.nii.gz")
    for restored_file in [reg_file, trf_file, mask_file, label_file]:
        if not os.path.isfile(restored_file):
            raise IOError("FileDoesNotExist: '{0}'.".format(restored_file))

snap_file = os.path.join(qcdir, "affine.pdf")
plot_image(t1_file, edge_file=reg_file, snap_file=snap_file, name="affine")
snap_file = os.path.join(qcdir, "cortex_mask.pdf")
plot_image(b0_file, overlay_file=mask_file, snap_file=snap_file,
           name="cortex mask")
snap_file = os.path.join(qcdir, "gyri.pdf")
if use_t1_space:
    plot_image(t1_file, overlay_file=label_file, snap_file=snap_file,
               name="gyri", overlay_cmap="cold_hot")
else:
    plot_image(b0_file, overlay_file=label_file, snap_file=snap_file,
               name="gyri", overlay_cmap="cold_hot")

"""
Tractography
------------

We use dipy to generate a streamline tractography. Play with the
'nb_seeds_per_voxel' parameter to increase/decrease the fiber density. Then
we count the tracks that start and end at each label pair.
"""

track_outdir = os.path.join(outdir, "streamline")
track_file = os.path.join(track_outdir, "fibers.txt")
if not os.path.isdir(track_outdir):
    os.mkdir(track_outdir)
if not os.path.isfile(track_file):
    deterministic(
        diffusion_file,
        bvecs_file,
        bvals_file,
        track_outdir,
        mask_file=mask_file,
        order=4,
        nb_seeds_per_voxel=nb_seeds_per_voxel,
        step=0.5,
        fmt="%.4f")

if use_vtk:
    tracks = loadtxt(track_file)
    actor = pvtk.line(tracks, scalar=1)
    pvtk.add(ren, actor)
    pvtk.record(ren, qcdir, "fibers", az_ang=45, n_frames=2)
    pvtk.clear(ren)
    

connect_outdir = os.path.join(outdir, "det_region_connectivity")
if not os.path.isdir(connect_outdir):
    os.mkdir(connect_outdir)
if len(os.listdir(connect_outdir)) == 0:
    if use_t1_space:
        proba_file, network_file = anatomical_connectivity_matrix(
            track_file,
            label_file,
            t1_file,
            fa_file,
            trf_file,
            connect_outdir,
            symmetric=True)

    else:
        proba_file, network_file = diffusion_connectivity_matrix(
            track_file,
            label_file,
            connect_outdir,
            symmetric=True)
else:
    proba_file = os.path.join(connect_outdir, "det_paths.nii.gz")
    network_file = os.path.join(connect_outdir, "det_network_matrix")
    for restored_file in [proba_file, network_file]:
        if not os.path.isfile(restored_file):
            raise IOError("FileDoesNotExist: '{0}'.".format(restored_file))

snap_file = os.path.join(qcdir, "fiber_density_map.pdf")
if use_t1_space:
    plot_image(t1_file, overlay_file=proba_file, snap_file=snap_file,
               name="density", overlay_cmap="cold_hot")
else:
    plot_image(b0_file, overlay_file=proba_file, snap_file=snap_file,
               name="density", overlay_cmap="cold_hot")
snap_file = os.path.join(qcdir, "det_gyri_connectogram.pdf")
plot_matrix(network_file, snap_file=snap_file, name="det gyri connectogram",
            transform=numpy.log1p)



