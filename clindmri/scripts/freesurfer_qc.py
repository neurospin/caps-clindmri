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
import copy

# Bredala import
try:
    import bredala
    bredala.USE_PROFILER = False
    bredala.register("clindmri.segmentation.freesurfer",
        names=["population_statistic", "parse_fs_lut"])
    bredala.register("clindmri.plot.slicer", names=["xyz_mosaics"])
    bredala.register("clindmri.plot.polar", names=["polarplot"])
    bredala.register("clindmri.segmentation.fsl", names=["fslreorient2std"])
except:
    pass

# Clindmri import
from clindmri.segmentation.freesurfer import population_statistic
from clindmri.segmentation.freesurfer import parse_fs_lut
from clindmri.segmentation.freesurfer import mri_binarize
from clindmri.plot.slicer import xyz_mosaics
from clindmri.plot.polar import polarplot
from clindmri.segmentation.fsl import fslreorient2std
from clindmri.extensions.freesurfer.reader import TriSurface


# Parameters to keep trace
__hopla__ = ["subjdir", "fsconfig", "fslconfig", "mosaics", "cutlower",
             "cutupper", "nbslices", "outdir"]


# Script documentation
doc = """
Freesurfer QC
~~~~~~~~~~~~~

Inspect the results returned by the FreeSurfer cortical reconstruction
pipeline.

Steps:

1- Create the population statistic
2- Create polar plots
3- Create t1 overlays mosaic
4- Create white/pial meash overlays mosaics

Command:

python $HOME/git/caps-clindmri/clindmri/scripts/freesurfer_qc.py \
    -v 2 \
    -k /etc/fsl/5.0/fsl.sh \
    -c /i2bm/local/freesurfer/SetUpFreeSurfer.sh \
    -d /neurospin/senior/nsap/data/V4/freesurfer \
    -s ag110371 \
    -o /neurospin/senior/nsap/data/V4/qc/freesurfer \
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
    "-o", "--outdir", dest="outdir", metavar="PATH", type=is_directory,
    help="the FreeSurfer qc home directory, default is 'fsdir'.")
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
    "-p", "--nbslices", dest="nbslices", required=True, type=int,
    help="number of slice in the 'display_mode' direction.")
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
auto_cutupper = False
if cutupper == 0:
    auto_cutupper = True
auto_cutlower = False
if cutlower == 0:
    auto_cutlower = True
nbslices = args.nbslices
if not os.path.isdir(subjdir):
    raise ValueError(
        "'{0}' is not a FreeSurfer subject folder.".format(subjdir))
outdir = args.outdir
if outdir is None:
    outdir = args.fsdir
    qcdir = os.path.join(subjdir, "qc")
else:
    qcdir = os.path.join(outdir, args.subjectid)
if not os.path.isdir(qcdir):
    os.mkdir(qcdir)
elif args.erase:
    shutil.rmtree(qcdir)

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
nb_values = numpy.asarray(fs_lut_colors.keys()).max() + 1
cmap = numpy.zeros((nb_values, 4), dtype=numpy.single)
for key, color in fs_lut_colors.items():
    cmap[key, :3] = color
    if color != (0, 0, 0):
        cmap[key, 3] = 200.
cmap /= 255.

# Set null opacity FreeSurfer 'unknown' labels
for label, label_name in fs_lut_names.items():
    if "unknown" in label_name.lower():
        cmap[label, 3] = 0.

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
origt1file = data["rawavg.native"]
for basename, fpath in data.items():
    reotfpath = os.path.join(qcdir, "reo." + os.path.basename(fpath))
    fslreorient2std(fpath, reotfpath, shfile=fslconfig)
    data[basename] = reotfpath

# Mosaics: t1 vs rest of the world
t1file = data.pop("rawavg.native")
mosaics = []
aseg_thr = 256
for basename, fpath in data.items():

    # Update color map opacities
    ccmap = copy.deepcopy(cmap)
    if "aparc" in basename:
        ccmap[:aseg_thr, 3] = 0.3
    else:
        ccmap[:aseg_thr, 3] = 0.6

    # Create mosaic files
    mosaics.extend(
        xyz_mosaics(t1file, fpath, nbslices, "t1-{0}".format(basename),
                    ccmap, qcdir, cutupper=cutupper, cutlower=cutlower))


"""
Create white/pial mesh overlays mosaics
"""
surf = {
    "lh.white.7.native": None,
    "rh.white.7.native": None,
    "lh.pial.7.native": None,
    "rh.pial.7.native": None
}
colors = {
    "lh.white.7.native": 41,
    "rh.white.7.native": 41,
    "lh.pial.7.native": 123,
    "rh.pial.7.native": 123
}
# Construct the surfaces mesh volume
t1im = nibabel.load(origt1file)
t1affine = t1im.get_affine()
t1shape = t1im.get_data().shape
meshfile = os.path.join(qcdir, "meshs.nii.gz")
mesharray = numpy.zeros(t1shape, dtype=numpy.uint)
for basename in surf:
    # > construc path
    fpath = os.path.join(mripath, basename)
    if not os.path.isfile(fpath):
        raise ValueError("'{0}' is not a valid file name.".format(fpath))
    surf[basename] = fpath

    # > load mesh
    name = basename.split(".")[1]
    annot_basename = basename.replace(
        name, "aparc.annot").replace(".native", "")
    annotfile = os.path.join(mripath, annot_basename)
    if not os.path.isfile(fpath):
        raise ValueError("'{0}' is not a valid file name.".format(fpath))
    surface = TriSurface.load(fpath, annotfile=annotfile)

    # > binarize mesh
    indices = numpy.round(surface.vertices).astype(int).T
    indices[0, numpy.where(indices[0] >= t1shape[0])] = 0
    indices[1, numpy.where(indices[1] >= t1shape[1])] = 0
    indices[2, numpy.where(indices[2] >= t1shape[2])] = 0
    mesharray[indices.tolist()] = colors[basename]

# Save the mesh volume
meshim = nibabel.Nifti1Image(mesharray, t1affine)
nibabel.save(meshim, meshfile)

# Need to reorient the image to MNI standard space
reomeshfile = os.path.join(qcdir, "reo." + os.path.basename(meshfile))
fslreorient2std(meshfile, reomeshfile, shfile=fslconfig)

# Update color map opacities
ccmap = copy.deepcopy(cmap)
ccmap[:aseg_thr, 3] = 1.0

# Create mosaic files
mosaics.extend(
    xyz_mosaics(t1file, reomeshfile, nbslices, "meshs", ccmap, qcdir,
                cutupper=cutupper, cutlower=cutlower))



