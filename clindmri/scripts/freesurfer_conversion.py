#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013-2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
from __future__ import print_function
import argparse
import os
import shutil

# Bredala import
import bredala
bredala.USE_PROFILER = False
bredala.register("clindmri.segmentation.freesurfer", names=["mri_convert",
                 "resample_cortical_surface", "conformed_to_native_space",
                 "surf_convert", "qc"])

# Clindmri import
from clindmri.segmentation.freesurfer import mri_convert
from clindmri.segmentation.freesurfer import resample_cortical_surface
from clindmri.segmentation.freesurfer import conformed_to_native_space
from clindmri.segmentation.freesurfer import surf_convert
from clindmri.segmentation.freesurfer import qc


# Script documentation
doc = """
Freesurfer inspection
~~~~~~~~~~~~~~~~~~~~~

Inspect the results returned by the FreeSurfer cortical reconstruction 
pipeline.

Steps:

1. Nifti conversions: aseg - aparc+aseg - aparc.a2009s+aseg - wm - t1.
   Export FreeSurfer '.mgz' images of interest in Nifti format. These
   images are resliced like the 'rawavg.mgz' file, have a '.native'
   suffix and are stored in a 'convert' folder. 
2. Registration matrix: between the conformed space (orig.mgz)
   and the native anatomical (rawavg.mgz).
3. Surface conversions: resample the white or pial FreeSurfer
   surfaces at different resolutions (impacts the number of vertex)
   with common mesh that can be directly used in a longitudinal
   setting. The results are also stored in a 'convert' folder with
   a '.native' suffix and the considered level in the file name. Vetex
   are expressed in the index coordinate system.
4. Quality controls: segmentations/t1 overlays - pial/white parcellation.
   Results are stored in the 'convert/qc' folder.

Command:

python $HOME/git/clindmri/scripts/freesurfer_conversion.py 
    -v 2
    -c /i2bm/local/freesurfer/SetUpFreeSurfer.sh
    -d /volatile/imagen/dmritest/freesurfer
    -s 000043561374
    -e
    -g

Local multi-processing:

from hopla import hopla
import os
myhome = os.environ["HOME"]
status, exitcodes = hopla(
    os.path.join(myhome, "git", "caps-clindmri", "clindmri", "scripts",
                 "freesurfer_conversion.py"), 
    c="/i2bm/local/freesurfer/SetUpFreeSurfer.sh",
    d="/volatile/imagen/dmritest/freesurfer",
    s=["000043561374", "000085724167", "000052904972"],
    e=True,
    hopla_iterative_kwargs=["s"],
    hopla_cpus=3,
    hopla_logfile="/volatile/imagen/dmritest/freesurfer/conversion.log",
    hopla_verbose=1)

Cluster multi-processing:

from hopla.converter import hopla
import hopla.demo as demo
import os

apath = os.path.abspath(os.path.dirname(demo.__file__))
script = os.path.join(os.path.dirname(demo.__file__),
                      "my_ls_script.py")
status, exitcodes = hopla(
    script, hopla_iterative_kwargs=["d"], d=[apath, apath], v=1,
    hopla_verbose=1, hopla_cpus=2, hopla_cluster=True,
    hopla_cluster_logdir="/home/ag239446/test/log",
    hopla_cluster_python_cmd="/usr/bin/python2.7",
    hopla_cluster_queue="Cati_LowPrio")
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
    "-e", "--errase", dest="errase", action="store_true",
    help="if activated, clean the result folder.")
parser.add_argument(
    "-c", "--config", dest="fsconfig", metavar="FILE", required=True,
    help="the FreeSurfer configuration file.", type=is_file)
parser.add_argument(
    "-d", "--fsdir", dest="fsdir", required=True, metavar="PATH",
    help="the FreeSurfer processing home directory.", type=is_directory)
parser.add_argument(
    "-s", "--subjectid", dest="subjectid", required=True,
    help="the subject identifier.")
parser.add_argument(
    "-g", "--graph", dest="graphics", action="store_true",
    help="if activated compute quality controls on the FreeSurfer outputs.")
args = parser.parse_args()


"""
First construct the subject FreeSurfer directory and check its existance on
the file system. Create a new 'convert' level in the FreeSurfer hierarchy
to store the conversion results.
"""
if args.verbose > 0:
    print("[info] Start FreeSurfer conversions...")
    print("[info] Directory: {0}.".format(args.fsdir))
    print("[info] Subject: {0}.".format(args.subjectid))
subjdir = os.path.join(args.fsdir, args.subjectid)
if not os.path.isdir(subjdir):
    raise ValueError(
        "'{0}' is not a FreeSurfer subject folder.".format(subjdir))
convertdir = os.path.join(subjdir, "convert")
if not os.path.isdir(convertdir):
    os.mkdir(convertdir)
elif args.errase:
    shutil.rmtree(convertdir)
    os.mkdir(convertdir)


"""
Step 1: Nifti conversions.
"""
if args.verbose > 0:
    print("[info] Start Nifti conversions...")
niftifiles = {}
for modality in ["aparc+aseg", "aparc.a2009s+aseg", "aseg", "wm", "rawavg"]:
    regex = os.path.join(args.subjectid, "mri", "{0}.mgz".format(modality))
    niftifiles[modality] = mri_convert(
        args.fsdir, regex, output_directory=None, reslice=True,
        interpolation="nearest", fsconfig=args.fsconfig)
    if args.verbose > 1:
        print("[result] {0}: {1}.".format(modality, niftifiles[modality]))


"""
Step 2: Registration matrix.
"""
if args.verbose > 0:
    print("[info] Start Registration matrix...")
regex = os.path.join(args.subjectid, "mri")
trffile = conformed_to_native_space(
    args.fsdir, regex, output_directory=None, fsconfig=args.fsconfig)
if args.verbose > 1:
    print("[result] trffile: {0}.".format(trffile))


"""
Step 3: Surface conversions.
"""
if args.verbose > 0:
    print("[info] Start surface conversions...")
surfaces = {}
annotations = []
for modality in ["pial", "white"]:
    for hemi in ["lh", "rh"]:
        name = "{0}.{1}".format(hemi, modality)
        regex = os.path.join(args.subjectid, "surf", name)
        resamplefiles, annotfiles = resample_cortical_surface(
            args.fsdir, regex, output_directory=None, orders=[4, 5, 6, 7],
            surface_name=modality,
            fsconfig=args.fsconfig)
        annotations.extend(annotfiles)
        surfaces[name] = surf_convert(
            args.fsdir, niftifiles["rawavg"], resamplefiles,
            output_directory=None,  rm_orig=True, fsconfig=args.fsconfig)
        if args.verbose > 1:
            print("[result] {0}: {1}.".format(name, surfaces[name]))
annotations = list(set(annotations))
if args.verbose > 1:
    print("[result] Annotations: {0}.".format(annotations))


"""
Step 4: Quality controls.
"""
if args.graphics:
    if args.verbose > 0:
        print("[info] Start quality control...")
    whitefiles = surfaces["lh.white"] + surfaces["rh.white"]
    pialfiles = surfaces["lh.pial"] + surfaces["rh.pial"]
    qcfiles = qc(
        niftifiles["rawavg"], niftifiles["wm"], niftifiles["aseg"],  whitefiles,
        pialfiles, annotations, actor_ang=[90., 180., 0.],
        output_directory=None, fsconfig=args.fsconfig)
    if args.verbose > 1:
        print("[result] QC: {0}.".format(qcfiles))
        
    
