#! /usr/bin/env python
# -*- coding: utf-8 -*-

##########################################################################
# NSAp - Copyright (C) CEA, 2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

import os
import logging
import time
import argparse

from clindmri.connectivity.connectogram import (register_diffusion_to_anatomy,
                                                create_masks_for_tractography,
                                                probtrackx2_connectogram,
                                                qc_dif2anat_registration,
                                                qc_tractography_masks,
                                                qc_connectogram,
                                                TractoMaskTypes)


def create_logger(log_dir):
    """
    Create a logger with console logging (info level) + log file (debug level).
    """
    logger = logging.getLogger(__file__)
    logger.setLevel(logging.INFO)

    # file logger
    log_filename = "probabilist_connectogram_%s.log" % time.strftime("%Y-%m-%d_%H:%M:%S")
    if log_dir:
        log_path = os.path.join(log_dir, log_filename)
    else:
        log_path = log_filename
    file_handler = logging.FileHandler(log_path)
    formatter = logging.Formatter('%(asctime)s :: %(message)s')
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)

    # console logger
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    logger.info("Log path: %s" % log_path)

    return logger


def main(outdir,
         subject_id,
         nodif_brain,
         bedpostx_dir,
         tracto_mask_type,
         cortical_atlas   = "Desikan",
         add_subcortical  = True,
         nsamples         = 5000,
         nsteps           = 2000,
         fs_subjects_dir  = None):
    """

    Required
    --------
    outdir:       Str, path to the directory where to output.
    subject_id:   Str, the subject identifier. The subject should exist in
                  the Freesurfer <SUBJECTS_DIR>.
    nodif_brain:  Str, path to a preprocessed brain-only volume with bvalue=0.
    bedpostx_dir: Str, path to the bedpostx output directory.
    tracto_mask_type: Str, the type of tractography mask to create, allowed types:
                  "nodif_brain" (whole brain), "wm", "wm_dilated_1vox_6conn" or
                  "wm_dilated_1vox_14conn".
                  Two of the proposed white matter masks are dilated because a
                  non-dilated white matter mask does not overlap with the "gray"
                  'subcortical regions, therefore the samples will never get there.
                  Moreover the right and left white matter regions are much less
                  connected without dilation, therefore the connectogram shows
                  few interhemisphere connections with a simple white matter mask.
                  "wm_dilated_1vox_6conn" means white matter dilated by 1 voxel
                  based one a 6-connexity structuring element.

    Optional
    --------
    cortical_atlas:  Str, the cortical atlas to use, "Desikan" (default) or "Destrieux".
    add_subcortical: Bool, if True include subcortical regions in the connectogram
                     (Thalamus, Caudate, Putamen, Pallidum, Hippocampus,
                      Amygdala, Accumbens-area, and VentralDC).
    nsamples:        Int (default 5000), number of samples per voxel to initiate
                     in seed region.
    nsteps:          Int (default 2000), maximum number of steps for a sample.
    fs_subjects_dir: Str, the Freesurfer <SUBJECTS_DIR>. By default the Freesurfer
                     $SUBJECTS_DIR environment variable is used.
    """
    # If <outdir> does not exist, create it
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    logger = create_logger(outdir)

    logger.info("outdir: %s" % outdir)

    if fs_subjects_dir is None:
        if "SUBJECTS_DIR" in os.environ:
            fs_subjects_dir = os.environ["SUBJECTS_DIR"]
        else:
            raise ValueError("Missing <SUBJECTS_DIR>: set the $SUBJECTS_DIR "
                             "environment variable for Freesurfer or pass it "
                             "as an argument.")

    logger.info("SUBJECTS_DIR: %s" % fs_subjects_dir)

    # Use boundary based registration from Freesurfer to infer an affine transformation
    # between diffusion space and anatomy (T1)
    logger.info("Registration of diffusion to T1 using Freesurfer bbregister...")
    dif2anat_dat = register_diffusion_to_anatomy(outdir,
                                                 nodif_brain,
                                                 subject_id,
                                                 fs_subjects_dir)

    logger.info("Creating snapshots to QC registration...")
    qc_dir = qc_dif2anat_registration(outdir,
                                      nodif_brain,
                                      dif2anat_dat,
                                      subject_id,
                                      fs_subjects_dir)
    logger.info("The snapshots are in %s" % qc_dir)

    # Creation of the masks for the tractography (connectogram computation):
    #   - a white matter mask = the tractography mask
    #   - a list of masks to be used as seeds in probtrackx2
    logger.info("Creation of seed and tractography masks...")
    logger.info("Cortical atlas: %s" % cortical_atlas)
    logger.info("Include subcortical regions: %s" % str(add_subcortical))
    seed_masks, tracto_mask, stop_mask, avoid_mask = \
        create_masks_for_tractography(outdir           = outdir,
                                      nodif_brain      = nodif_brain,
                                      dif2anat_dat     = dif2anat_dat,
                                      subject_id       = subject_id,
                                      cortical_atlas   = cortical_atlas,
                                      add_subcortical  = add_subcortical,
                                      tracto_mask_type = tracto_mask_type,
                                      fs_subjects_dir  = fs_subjects_dir)

    logger.info("Creating snapshots to QC masks...")
    qc_dir = qc_tractography_masks(outdir          = outdir,
                                   nodif_brain     = nodif_brain,
                                   tracto_mask     = tracto_mask,
                                   seed_masks      = seed_masks,
                                   subject_id      = subject_id,
                                   cortical_atlas  = cortical_atlas,
                                   fs_subjects_dir = fs_subjects_dir)
    logger.info("The snapshots are in %s" % qc_dir)

    # Compute the connectogram (network_matrix file)
    logger.info("Running probtrackx2, nsamples: %i, nsteps: %i" % (nsamples, nsteps))
    logger.info("Tractography mask: %s" % tracto_mask)
    proba_matrix, network_matrix = \
        probtrackx2_connectogram(outdir       = outdir,
                                 bedpostx_dir = bedpostx_dir,
                                 seed_masks   = seed_masks,
                                 tracto_mask  = tracto_mask,
                                 stop_mask    = stop_mask,
                                 avoid_mask   = avoid_mask,
                                 subdir       = "probtrackx2",
                                 nsamples     = nsamples,
                                 nsteps       = nsteps)

    # Create connectogram PNG file
    logger.info("Creating snapshots of the connectogram...")
    qc_dir = qc_connectogram(outdir         = outdir,
                             tracto_mask    = tracto_mask,
                             proba_matrix   = proba_matrix,
                             network_matrix = network_matrix,
                             seed_masks     = seed_masks)
    logger.info("The napshots are in %s" % qc_dir)


def get_cmd_line_args():
    """
    Create a command line argument parser, run it and return a dict mapping
    <argument name> -> <argument value>.
    """
    usage = "%(prog)s <subject id> <nodif_brain> <bedpostx_dir> <outdir> [options]"
    parser = argparse.ArgumentParser(prog  = "python probabilist_connectogram.py",
                                     usage = usage)

    # Required arguments
    parser.add_argument("subject_id",
                        help="The name of the subject's folder in <SUBJECTS_DIR>.")
    parser.add_argument("nodif_brain",
                        help="A preprocessed brain-only volume with bvalue=0.")
    parser.add_argument("bedpostx_dir",
                        help="The bedpostx output directory for the subject's DWI data.")
    parser.add_argument("outdir", help="Directory where to output.")

    # Optional arguments
    parser.add_argument("--cortical-atlas", default="Desikan",
                        choices=["Desikan", "Destrieux"], metavar="<atlas name>",
                        help="Cortical atlas name, 'Desikan' (default) or 'Destrieux'")
    parser.add_argument("--remove-subcortical", action="store_true",
                        help="Remove subcortical regions from the connectogram "
                             "(Thalamus, Caudate, Putamen, Pallidum, Hippocampus, "
                             "Amygdala, Accumbens-area and VentralDC).")
    parser.add_argument("--tracto-mask-type", default="nodif_brain",
                        choices=TractoMaskTypes.choices, metavar="<tracto mask type>",
                        help='The type of tractography mask to create, allowed types: '
                             '"nodif_brain" (default, whole brain), "wm", '
                             '"wm_dilated_1vox_6conn" or "wm_dilated_1vox_14conn". '
                             'Two of the proposed white matter masks are dilated because a '
                             'non-dilated white matter mask does not overlap with the "gray" '
                             'subcortical regions, therefore the samples will never get there. '
                             'Moreover the right and left white matter regions are much less '
                             'connected without dilation, therefore the connectogram shows '
                             'few interhemisphere connections with a simple white matter mask. '
                             '"wm_dilated_1vox_6conn" means white matter dilated by 1 voxel '
                             'based one a 6-connexity structuring element.')
    parser.add_argument("--nsamples", type=int, default=5000, metavar="<nsamples>",
                        help="Number of samples per voxel to initiate in seed "
                             "region (default 5000).")
    parser.add_argument("--nsteps", type=int, default=2000, metavar="<nsteps>",
                        help="Maximum number of steps for a sample (default 2000).")
    parser.add_argument("--fs-subjects-dir", metavar="<Freesurfer subjects directory>",
                        help="To bypass the $SUBJECTS_DIR environment variable.")

    # Create a dict of arguments to pass to the 'main' function
    args   = parser.parse_args()
    kwargs = vars(args)

    # Adapt one argument to the 'main' interface
    kwargs["add_subcortical"] = not kwargs["remove_subcortical"]
    del kwargs["remove_subcortical"]

    return kwargs


if __name__ == "__main__":

    kwargs = get_cmd_line_args()
    print kwargs
    main(**kwargs)

