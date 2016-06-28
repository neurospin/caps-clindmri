#! /usr/bin/env python
# -*- coding: utf-8 -*-

##########################################################################
# NSAp - Copyright (C) CEA, 2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

import os
import argparse

from clindmri.connectivity.connectogram import CORTICAL_ATLASES, STOP_MASK_TYPES
from clindmri.connectivity.connectogram_pipeline import connectogram_seeding_wm_pipeline


def is_dir(dirpath):
    """ Check direcory's existence - argparse 'type' argument.
    """
    if not os.path.isdir(dirpath):
        raise argparse.ArgumentError("Directory does not exist: %s" % dirpath)
    return dirpath


def is_file(filepath):
    """ Check file's existence - argparse 'type' argument.
    """
    if not os.path.isfile(filepath):
        raise argparse.ArgumentError("File does not exist: %s" % filepath)
    return filepath


def get_cmd_line_args():
    """
    Create a command line argument parser, run it and return a dict mapping
    <argument name> -> <argument value>.
    """
    usage = ("%(prog)s -s <subject id> -i <nodif_brain> -m <nodif_brain_mask> "
             "-b <bedpostx_dir> -o <outdir> [options]")
    parser = argparse.ArgumentParser(prog="python probabilist_connectogram.py",
                                     usage=usage)

    # Required arguments
    parser.add_argument("-s", "--subject-id",
                        required=True,
                        help="Subject id used with Freesurfer.")
    parser.add_argument("-i", "--nodif-brain",
                        type=is_file,
                        required=True,
                        help="A preprocessed brain-only volume with bvalue=0.")
    parser.add_argument("-m", "--nodif-brain-mask",
                        type=is_file,
                        required=True,
                        help="A brain-only mask volume in diffusion space.")
    parser.add_argument("-b", "--bedpostx-dir",
                        type=is_dir,
                        required=True,
                        help="The bedpostx output directory.")
    parser.add_argument("-o", "--outdir",
                        required=True,
                        help="Directory where to output.")

    # Optional arguments
    parser.add_argument("-a", "--cortical-atlas",
                        default="Desikan",
                        choices=CORTICAL_ATLASES,
                        metavar="<atlas name>",
                        help="Cortical atlas name")
    parser.add_argument("-p", "--stop-mask-type",
                        default="target_rois",
                        choices=STOP_MASK_TYPES,
                        metavar="<stop mask type>",
                        help='Which type of stop mask to use:'
                             '\n- "target_rois": stop a sample as soon as it reaches a target region'
                             '\n- "inverse_wm":  stop a sample as soon as it leaves the white matter')
    parser.add_argument("-l", "--nsamples",
                        type=int,
                        default=5000,
                        metavar="<nsamples>",
                        help="Number of samples per voxel to initiate in seed region (default 5000).")
    parser.add_argument("-e", "--nsteps",
                        type=int,
                        default=2000,
                        metavar="<nsteps>",
                        help="Maximum number of steps for a sample (default 2000).")
    parser.add_argument('-g', "--steplength",
                        type=float,
                        default=0.5,
                        metavar="<steplength>",
                        help="Step length in mm.")
    parser.add_argument("-d", "--subjects-dir",
                        metavar="<Freesurfer subjects directory>",
                        help="To bypass the $SUBJECTS_DIR environment variable.")
    parser.add_argument("-f", "--fsl-init",
                        default="/etc/fsl/5.0/fsl.sh",
                        metavar="<FSL init script>",
                        help="To initialize FSL's environment.")

    # Create a dict of arguments to pass to the 'main' function
    args = parser.parse_args()
    kwargs = vars(args)

    return kwargs


kwargs = get_cmd_line_args()
print kwargs
connectogram_seeding_wm_pipeline(**kwargs)
