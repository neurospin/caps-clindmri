#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" A script for cortical connectivity using FreeSurfer Labels and FSL
probabilistic tractography.

We will describe a method for performing probabilistic tractography in a single
subject using automatically generated cortical labels from Freesurfer resulting
in a probabilistic region connectivity matrix.

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

This analysis is using Freesurfer & FSL & nilearn, all of them beeing freely
available.
"""

import os
import shutil
import numpy
import nibabel
import glob
from clindmri.extensions.freesurfer import read_cortex_surface_segmentation
from clindmri.estimation.fsl import dtifit
from clindmri.segmentation.freesurfer import cortex
from clindmri.registration.utils import extract_image
from clindmri.segmentation.fsl import bet2
from clindmri.tractography.fsl import bedpostx_datacheck
from clindmri.tractography.fsl import bedpostx
from clindmri.tractography.fsl import probtrackx2
from clindmri.plot.slicer import plot_image
from clindmri.plot.slicer import plot_matrix
import clindmri.plot.pvtk as pvtk


t1_file = "/volatile/imagen/dmritest/001/raw/t1-1mm-1-001.nii.gz"
fsdir = "/volatile/imagen/dmritest/001/processed/fs"
diffusion_file = "/volatile/imagen/dmritest/001/raw/hardi-b1500-1-001.nii.gz"
bvecs_file = "/volatile/imagen/dmritest/001/raw/hardi-b1500-1-001.bvec"
bvals_file = "/volatile/imagen/dmritest/001/raw/hardi-b1500-1-001.bval"
outdir = "/volatile/imagen/dmritest/001/processed/probaconnect"
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
        actor = pvtk.surface(surface.vertices, surface.triangles,
                             surface.labels, ctab)
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
* PDFs characterizing the underlying diffusion process.

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
    mask_file, label_file, seeds, reg_file, trf_file = cortex(
        t1_file,
        fsdir,
        cortex_outdir,
        dest_file=fa_file,
        generate_mask=True,
        generate_seeds=True)
else:
    seeds_outdir = os.path.join(cortex_outdir, "gyri")
    reg_file = os.path.join(cortex_outdir, "cortex_dest_to_t1.nii.gz")
    mask_file = os.path.join(cortex_outdir, "cortex_mask.nii.gz")
    label_file = os.path.join(cortex_outdir, "cortex_gyri_labels.nii.gz")
    for restored_file in [reg_file, mask_file]:
        if not os.path.isfile(restored_file):
            raise IOError("FileDoesNotExist: '{0}'.".format(restored_file))
    if not os.path.isdir(seeds_outdir):
        raise IOError("DirectoryDoesNotExist: '{0}'.".format(seeds_outdir))
    seeds = [os.path.join(seeds_outdir, name) for
             name in os.listdir(seeds_outdir)]

snap_file = os.path.join(qcdir, "affine.pdf")
plot_image(t1_file, edge_file=reg_file, snap_file=snap_file, name="affine")
snap_file = os.path.join(qcdir, "cortex_mask.pdf")
plot_image(b0_file, overlay_file=mask_file, snap_file=snap_file,
           name="cortex mask")
snap_file = os.path.join(qcdir, "gyri.pdf")
plot_image(b0_file, overlay_file=label_file, snap_file=snap_file,
           name="gyri", overlay_cmap="cold_hot")

"""
Generating PDFs
---------------

We use 'bedpostx' to generate PDFs of the diffusion direction and get on
with tractography. 'bedpostx' takes about 5 hours of compute time. This routine
need specific files that are checked with the 'bedpostx_datacheck' command.
"""

bedpostx_indir = os.path.join(outdir, "bedpostx")
bedpostx_outdir = os.path.join(outdir, "bedpostx.bedpostX")
if not os.path.isdir(bedpostx_indir):
    os.mkdir(bedpostx_indir)
if len(os.listdir(bedpostx_outdir)) == 0:
    shutil.copy2(mask_file, bedpostx_indir)
    data_ext = ".".join(diffusion_file.split(".")[1:])
    shutil.copy2(diffusion_file, os.path.join(bedpostx_indir,
                                              "data." + data_ext))
    shutil.copy2(bvecs_file, os.path.join(bedpostx_indir, "bvecs"))
    shutil.copy2(bvals_file, os.path.join(bedpostx_indir, "bvals"))
    if not bedpostx_datacheck(bedpostx_indir):
        raise IOError("'{0}' does not contain the data expected by "
                      "'bedpostx'.".format(bedpostx_indir))

    (bedpostx_outdir, merged_th, merged_ph,
     merged_f, mean_th, mean_ph,
     mean_f, dyads) = bedpostx(
        bedpostx_indir)
else:
    merged_files = glob.glob(os.path.join(bedpostx_outdir, "merged*"))
    if len(merged_files) == 0:
        raise IOError("FilesDoNotExist: in '{0}'.".format(bedpostx_outdir))

"""
Tractography
============

We compute the fibers starting at a specific gyrus, connecting the other
gyri only with the 'network' option. The result is stored in the 'fdt_paths'
image. The 'fdt_network_matrix' is a 2D region-by-region connectivity matrix.
"""

probtrackx2_outdir = os.path.join(outdir, "proba_region_connectivity")
masks_file = os.path.join(probtrackx2_outdir, "masks.txt")
if not os.path.isdir(probtrackx2_outdir):
    os.mkdir(probtrackx2_outdir)
if not os.path.isfile(masks_file):
    numpy.savetxt(masks_file, numpy.asarray(seeds), delimiter=" ", fmt="%s")
    proba_files, network_file = probtrackx2(
        network=True,
        seed=masks_file,
        loopcheck=True,
        onewaycondition=True,
        samples=os.path.join(bedpostx_outdir, "merged"),
        mask=mask_file,
        dir=probtrackx2_outdir)
else:
    proba_file = glob.glob(os.path.join(probtrackx2_outdir, "fdt_paths*"))[0]
    network_file = os.path.join(probtrackx2_outdir, "fdt_network_matrix")
    weights_file = os.path.join(probtrackx2_outdir, "waytotal")
    for restored_file in [network_file, proba_file, weights_file]:
        if not os.path.isfile(restored_file):
            raise IOError("FileDoesNotExist: '{0}'.".format(restored_file))

snap_file = os.path.join(qcdir, "fiber_density_map.pdf")
plot_image(b0_file, overlay_file=proba_file, snap_file=snap_file,
           name="density", overlay_cmap="cold_hot")
snap_file = os.path.join(qcdir, "prob_gyri_connectogram.pdf")
plot_matrix(network_file, snap_file=snap_file, name="prob gyri connectogram",
            transform=numpy.log1p)
