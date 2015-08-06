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
in a probabilistice region connectivity matrix.

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

This analysis is using Freesurfer & FSL, both of which are freely available.
"""

import os
import shutil
import numpy
from clindmri.estimation.fsl import dtifit
from clindmri.segmentation.freesurfer import cortex
from clindmri.registration.utils import extract_image
from clindmri.segmentation.fsl import bet2
from clindmri.tractography.fsl import bedpostx_datacheck
from clindmri.tractography.fsl import bedpostx


t1_file = "/volatile/imagen/dmritest/001/raw/t1-1mm-1-001.nii.gz"
fsdir = "/volatile/imagen/dmritest/001/processed/fs"
diffusion_file = "/volatile/imagen/dmritest/001/raw/hardi-b1500-1-001.nii.gz"
bvecs_file = "/volatile/imagen/dmritest/001/raw/hardi-b1500-1-001.bvec"
bvals_file = "/volatile/imagen/dmritest/001/raw/hardi-b1500-1-001.bval"
outdir = "/volatile/imagen/dmritest/001/processed/probaconnect"


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
    raise NotImplementedError
    # recon-all -s 0001 -i T1.nii.gz


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
extract_image(diffusion_file,
              index=b0_index,
              out_file=b0_file)


b0_brain_file = os.path.join(outdir, "nodif_brain")
(output, mask_file, mesh_file, 
 outline_file, inskull_mask_file,
 inskull_mesh_file, outskull_mask_file,
 outskull_mesh_file, outskin_mask_file,
 outskin_mesh_file, skull_mask_file) = bet2(b0_file,
                                            b0_brain_file,
                                            m=True,
                                            f=0.25)

"""
Generating the FA image
-----------------------

There are much better ways to create an FA image than this method, but we're
only using this image for registration purposes. It is probably not a good
idear to use this image in a whole-brain FA analysis.
"""

dtifit_outdir = os.path.join(outdir, "dtifit")
if len(os.listdir(dtifit_outdir)) == 0:
    (v1_file, v2_file, v3_file, l1_file,
     l2_file, l3_file, md_file, fa_file, 
     s0_file, tensor_file, m0_file) = dtifit(diffusion_file, 
                                             bvecs_file,
                                             bvals_file,
                                             mask_file,
                                             dtifit_outdir)
else:
    fa_file = os.path.join(dtifit_outdir, "dtifit_FA.nii.gz")
    if not os.path.isfile(fa_file):
        raise IOError("FileDoesNotExist: '{0}'.".format(fa_file))


"""
Registration of the FA image to the T1 image
--------------------------------------------

Register the FA image to the T1 image and generate a 'cortex_whitematter-mask'
and 'cortex_gyri-labels' in the FA space. It also generate seed regions
to compute the connectivity matrix in the FA space in the 'gyri' folder.
"""

cortex_outdir = os.path.join(outdir, "cortex")
if len(os.listdir(cortex_outdir)) == 0:
    mask_file, label_file, seeds = cortex(t1_file,
                                          fsdir,
                                          cortex_outdir,
                                          dest_file=fa_file,
                                          generate_mask=True,
                                          generate_seeds=True)
else:
    seeds_outdir = os.path.join(cortex_outdir, "gyri")
    if not os.path.isdir(seeds_outdir):
        raise IOError("DirectoryDoesNotExist: '{0}'.".format(seeds_outdir))
    seeds = [os.path.join(seeds_outdir, name) for
             name in os.listdir(seeds_outdir)]

"""
Generating PDFs
---------------

We use 'bedpostx' to generate PDFs of the diffusion direction and get on
with tractography. 'bedpostx' takes about 5 hours of compute time. This routine
need specific files that are checked with the 'bedpostx_datacheck' command.
"""

bedpostx_indir = os.path.join(outdir, "bedpostx")
if not os.path.isdir(bedpostx_indir):
    os.mkdir(bedpostx_indir)
shutil.copy2(mask_file, bedpostx_indir)
data_ext = ".".join(diffusion_file.split(".")[1:])
shutil.copy2(diffusion_file, os.path.join(bedpostx_indir, "data." + data_ext))
shutil.copy2(bvecs_file, os.path.join(bedpostx_indir, "bvecs"))
shutil.copy2(bvals_file, os.path.join(bedpostx_indir, "bvals"))
if not bedpostx_datacheck(bedpostx_indir):
    raise IOError("'{0}' does not contain the data expected by "
                  "'bedpostx'.".format(bedpostx_indir))

(outdir, merged_th, merged_ph,
 merged_f, mean_th, mean_ph,
 mean_f, dyads) = bedpostx(bedpostx_indir)


""" 
Tractography
============

At this point, the jobs for each region can be run independently.
We compute the fibers starting at a specific gyrus, connecting the other
gyri only with the 'network' option. The result is stored in the 'fdt_paths'
image. The 'matrix_seeds_to_all_targets' is a 2D voxel-by-target matrix.

By collapsing across all seed voxels and dividing by the total number of
streamlines generated during the run, a 1xN array of percentages representing
the proportion of streamlines that reached each target.

By doing this for the N regions and stacking the arrays, we generate a NxN
connectivity matrix.
"""

probtrackx2_outdir = os.path.join(outdir, "connectivity")
masks_file = os.path.join(probtrackx2_outdir, "masks.txt")
numpy.savetxt(masks_file, numpy.asarray(seeds), delimiter=" ", fmt="%s")
probtrackx2(network=True,
            seed=masks_file,
            loopcheck=True,
            onewaycondition=True,
            samples=os.path.join(outdir, "merged"),
            mask=mask_file,
            dir=probtrackx2_outdir)




















