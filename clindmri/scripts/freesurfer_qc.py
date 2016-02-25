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
import numpy
import nibabel

# Bredala import
try:
    import bredala
    bredala.USE_PROFILER = False
    bredala.register("clindmri.segmentation.freesurfer",
                     names=["population_statistic", "parse_fs_lut"])
    bredala.register("clindmri.plot.slicer", names=["plot_image"])
    bredala.register("clindmri.plot.polar", names=["polarplot"])
    bredala.register("clindmri.segmentation.fsl", names=["fslreorient2std"])
except:
    pass

# Clindmri import
from clindmri.segmentation.freesurfer import population_statistic
from clindmri.segmentation.freesurfer import parse_fs_lut
from clindmri.plot.slicer import plot_image
from clindmri.plot.polar import polarplot
from clindmri.segmentation.fsl import fslreorient2std


# Parameters to keep trace
__hopla__ = ["subjdir", "fsconfig", "fslconfig", "mosaics", "cutlower",
             "cutupper", "cutcoords"]


# Script documentation
doc = """
Freesurfer QC
~~~~~~~~~~~~~

Inspect the results returned by the FreeSurfer cortical reconstruction
pipeline.

Steps:

1. Create the population statistic
2. Create polar plots
3. Create t1 overlays mosaic

Command:

python $HOME/git/caps-clindmri/clindmri/scripts/freesurfer_qc.py \
    -v 2 \
    -k /etc/fsl/5.0/fsl.sh \
    -c /i2bm/local/freesurfer/SetUpFreeSurfer.sh \
    -d /neurospin/senior/nsap/data/V4/freesurfer \
    -s ag110371 \
    -l 50 \
    -u 20 \
    -p 14 \
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
    help="if activated, clean the result folder.")
parser.add_argument(
    "-c", "--fsconfig", dest="fsconfig", metavar="FILE", required=True,
    help="the FreeSurfer configuration file.", type=is_file)
parser.add_argument(
    "-k", "--fslconfig", dest="fslconfig", metavar="FILE", required=True,
    help="the FSL configuration file.", type=is_file)
parser.add_argument(
    "-d", "--fsdir", dest="fsdir", required=True, metavar="PATH",
    help="the FreeSurfer processing home directory.", type=is_directory)
parser.add_argument(
    "-s", "--subjectid", dest="subjectid", required=True,
    help="the subject identifier.")
parser.add_argument(
    "-l", "--cutlower", dest="cutlower", default=0, type=int,
    help="exclude the 'cutlower' first slices.")
parser.add_argument(
    "-u", "--cutupper", dest="cutupper", default=0, type=int,
    help="exclude the 'cutupper' last slices.")
parser.add_argument(
    "-p", "--cutcoords", dest="cutcoords", nargs="+", required=True, type=int,
    help="the slicing strategy a 1-uplet to slice in the 'display_mode' direction.")
args = parser.parse_args()


"""
First construct the subject FreeSurfer directory and check its existance on
the file system. Create a new 'qc' level in the FreeSurfer hierarchy
to store the qc results.
"""
if args.verbose > 0:
    print("[info] Start FreeSurfer QC...")
    print("[info] Directory: {0}.".format(args.fsdir))
    print("[info] Subject: {0}.".format(args.subjectid))
subjdir = os.path.join(args.fsdir, args.subjectid)
fsconfig = args.fsconfig
fslconfig = args.fslconfig
cutlower = args.cutlower
cutupper = args.cutupper
if len(args.cutcoords) == 1:
    cutcoords = args.cutcoords[0]
elif len(args.cutcoords) == 3:
    cutcoords = args.cutcoords
else:
    raise Exception("Unsupported 'cutcoords' format.")
if not os.path.isdir(subjdir):
    raise ValueError(
        "'{0}' is not a FreeSurfer subject folder.".format(subjdir))
qcdir = os.path.join(subjdir, "qc")
if not os.path.isdir(qcdir):
    os.mkdir(qcdir)
elif args.erase:
    shutil.rmtree(qcdir)
    os.mkdir(qcdir)

"""
Create the population statistic and get the subjects measures
"""
popstats = population_statistic(args.fsdir)
indstats = population_statistic(args.fsdir, args.subjectid)


"""
Create polar plots
"""
for name, cohort_stats in popstats.items():
    individual_stats = indstats[name]
    snapfile = os.path.join(qcdir, "polarplot-{0}.png".format(name))
    polarplot(individual_stats, cohort_stats, snapfile,
              name="polarplot-{0}".format(name))

"""
Compute t1-images overlay mosaics
"""
# Create the FreeSurfer LUT
fs_lut_names, fs_lut_colors = parse_fs_lut(os.path.join(
    os.path.dirname(fsconfig), "FreeSurferColorLUT.txt"))
cmap = []
nb_values = numpy.asarray(fs_lut_colors.keys()).max()
cmap = numpy.zeros((nb_values, 4), dtype=numpy.single)
for key, color in fs_lut_colors.items():
    if key > 0:
        cmap[key - 1, :3] = color
cmap[:, 3] = 160.
cmap /= 255.

# Get file path
data = {
    "aparc.a2009s+aseg.native": None,
    "rawavg.native": None,
    "aseg.native": None,
    "wm.native": None,
    "aparc+aseg.native": None
}
mripath = os.path.join(args.fsdir, args.subjectid, "convert")
for basename in data:
    fpath = os.path.join(mripath, basename + ".nii.gz")
    if not os.path.isfile(fpath):
        raise ValueError("'{0}' is not a valid file name.".format(fpath))
    data[basename] = fpath

# Need to reorient the images to MNI standard space
for basename, fpath in data.items():
    reotfpath = os.path.join(qcdir, "reo." + os.path.basename(fpath))
    fslreorient2std(fpath, reotfpath, shfile=fslconfig)
    data[basename] = reotfpath

# Mosaics: t1 vs rest of the world
t1file = data.pop("rawavg.native")
mosaics = []
for basename, fpath in data.items():

    # Troncate the color map based on the label max
    array = nibabel.load(fpath).get_data()
    order = sorted(set(array.flatten()))
    ccmap = cmap[order[1]: order[-1] + 1]

    # Compute overlay mosaics
    for axis in ["x", "y", "z"]:
        qcname = "{0}t1-{1}".format(axis, basename)
        snap_file = os.path.join(qcdir, qcname + ".pdf")
        plot_image(t1file, overlay_file=fpath, snap_file=snap_file,
                   name=qcname, overlay_cmap=ccmap, cut_coords=cutcoords,
                   cutlower=cutlower, cutupper=cutupper,
                   display_mode=axis)
        mosaics.append(snap_file)

