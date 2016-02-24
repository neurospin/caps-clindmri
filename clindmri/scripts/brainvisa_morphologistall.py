#!/usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
# 
# First activate brainvisa env:
# /neurospin/brainvisa/build/Ubuntu-14.04-x86_64/bug_fix/bin/bv_env
##########################################################################

# System import
from __future__ import print_function
import argparse
import os
import shutil


# Bredala import
try:
    import bredala
    bredala.USE_PROFILER = False
    bredala.register("clindmri.segmentation.morphologist",
                     names=["morphologist_all"])
    bredala.register("clindmri.extensions.morphologist", names=["set_bvenv"])
except:
    pass

# Clindmri import
from clindmri.segmentation.morphologist import morphologist_all
from clindmri.extensions.morphologist import set_bvenv


# Parameters to keep trace
__hopla__ = ["t1file", "subjectdir", "wffile", "wfid", "wfstatus", "spmrun",
             "spmdir"]


# Script documentation
doc = """
Morphologist segmentation
~~~~~~~~~~~~~~~~~~~~~~~~~

Performs all the Morphologist steps.

Steps:

1- Ensure image orientation and reorient it if needed (Prepare Subject for
   Anatomical Pipeline).
2- Computation of a brain mask (Brain Mask Segmentation).
3- Computation of a mask for each hemisphere (Split Brain Mask).
4- A grey/white classification of each hemisphere to perform "Voxel Based
   Morphometry" (Grey White Classification) and spherical triangulation of
   cortical hemispheres (Grey White Surface).
5- Spherical triangulation of the external interface of the cortex of one or
   two hemispheres (Get Spherical Hemi Surface).
6- Computation of a graph representing the cortical fold topography
   (Cortical Fold Graph).
7- Automatic identification of the cortical sulci (Automatic Sulci
   Recognition), located in the "sulci" toolbox.

To check the worklow submission, use the 'soma_workflow_gui' command.

Command:

python $HOME/git/caps-clindmri/clindmri/scripts/brainvisa_morphologistall.py \
    -v 2 \
    -r /i2bm/local/spm8-standalone/run_spm8.sh \
    -p /i2bm/local/spm8-standalone \
    -b /neurospin/brainvisa/build/Ubuntu-14.04-x86_64/bug_fix/bin/bv_env \
    -o /tmp/morphologist \
    -n morphologist \
    -s subject0001 \
    -a /tmp/morphologist/3DT1.nii.gz \
    -e
"""


def is_file(filearg):
    """ Type for argparse - checks that file exists but does not open.
    """
    if not os.path.isfile(filearg):
        raise argparse.ArgumentError(
            "The file '{0}' does not exist!".format(filearg))
    return filearg


def is_directory(dirarg):
    """ Type for argparse - checks that directory exists.
    """
    if not os.path.isdir(dirarg):
        raise argparse.ArgumentError(
            "The directory '{0}' does not exist!".format(dirarg))
    return dirarg


parser = argparse.ArgumentParser(description=doc)
parser.add_argument(
    "-v", "--verbose", dest="verbose", type=int, choices=[0, 1, 2], default=0,
    help="increase the verbosity level: 0 silent, [1, 2] verbose.")
parser.add_argument(
    "-e", "--erase", dest="erase", action="store_true",
    help="if activated, clean the subject folder.")
parser.add_argument(
    "-p", "--spmdir", dest="spmdir", metavar="PATH", required=True,
    help="the standalone SPM directory.", type=is_directory)
parser.add_argument(
    "-r", "--spmrun", dest="spmrun", metavar="FILE", required=True,
    help="the path to the standalone SPM execution file.", type=is_file)
parser.add_argument(
    "-b", "--bvenv", dest="bvenv", metavar="FILE", required=True,
    help="the 'bvenv' binary path.", type=is_file)
parser.add_argument(
    "-o", "--outdir", dest="outdir", required=True, metavar="PATH",
    help="the output directory processing home directory.", type=is_directory)
parser.add_argument(
    "-n", "--name", dest="name", default="morphologist",
    help="a name used to build the result folder.")
parser.add_argument(
    "-s", "--subjectid", dest="subjectid", required=True,
    help="the subject identifier used also to build the result folder.")
parser.add_argument(
    "-a", "--anatfile", dest="anatfile", metavar="FILE", required=True,
    help="the subject anatomical image to be processed.", type=is_file)
args = parser.parse_args()


"""
First check if the subject morphologist directory exists on the file system,
and clean it if requested.
"""
if args.verbose > 0:
    print("[info] Start Morphologist morphologist_all...")
    print("[info] Directory: {0}.".format(args.outdir))
    print("[info] Name: {0}.".format(args.name))
    print("[info] Subject: {0}.".format(args.subjectid))
    print("[info] Anatomy: {0}.".format(args.anatfile))
    print("[info] SPM dir: {0}.".format(args.spmdir))
    print("[info] SPM run: {0}.".format(args.spmrun))
t1file = args.anatfile
subjectdir = os.path.join(args.outdir, args.name, args.subjectid)
spmrun = args.spmrun
spmdir = args.spmdir
if os.path.isdir(subjectdir) and args.erase:
    shutil.rmtree(subjectdir)


"""
Set the brainvisa environment if not found.
"""
set_bvenv(args.bvenv)


"""
Morphologist: all steps
"""
wffile, wfid, wfstatus = morphologist_all(
    t1file, args.subjectid, args.outdir, study=args.name, waittime=10,
    spmexec=spmrun, spmdir=spmdir)
if args.verbose > 1:
    print("[result] In folder: {0}.".format(subjectdir))
    print("[result] Workflow file: {0}.".format(wffile))
    print("[result] Workflow id: {0}.".format(wfid))
    print("[result] Workflow status: {0}.".format(wfstatus))

