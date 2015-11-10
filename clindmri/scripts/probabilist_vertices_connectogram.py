#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013-2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System imports
from __future__ import print_function
import os
import argparse
import shutil
import glob
import numpy

# Bredala import
import bredala
bredala.USE_PROFILER = False
bredala.register("clindmri.connectivity.fsl", names=["get_profile"])
bredala.register("clindmri.tractography.fsl", names=["probtrackx2"])
bredala.register("clindmri.registration.fsl", names=["flirt"])
bredala.register("clindmri.segmentation.freesurfer", names=["mri_vol2surf",
    "conformed_to_native_space"])

# Clindmri imports
from clindmri.connectivity.fsl import get_profile
from clindmri.segmentation.freesurfer import conformed_to_native_space
from clindmri.extensions.freesurfer import read_cortex_surface_segmentation
from clindmri.extensions.freesurfer.reader import TriSurface
from clindmri.registration.fsl import flirt
from clindmri.extensions.fsl import flirt2aff


# Parameters to keep trace
__hopla__ = ["bedpostxdir", "fsdir", "subjectid", "t1_file", "nodif_file",
             "outdir", "trf_file", "dat_file", "ico_order", "hemi"]


# Script documentation
doc = """
Probabilist connectogram
~~~~~~~~~~~~~~~~~~~~~~~~

A script for cortical connectivity using FreeSurfer cortical surface
vertices and FSL probabilistic tractography.

We will describe a method for performing FSL probabilistic tractography in a
single subject using automatically generated cortical surface vertices from
Freesurfer resulting in a probabilistic vertices connectivity matrix.

Prerequisites
-------------

This analysis requires the following MR data:

* A T1 image: this should cover the entire brain at high-resolution
  (~1mm3 voxels).
* The diffusion image(s).
* The subject FSL bedpostx modelization folder
* The subject FreeSurfer folder

Software
--------

This analysis is done using Freesurfer & FSL, all of them beeing freely
available.

List options
------------

>> python probabilist_vertices_connectogram.py --help

Example
-------

>> python $HOME/git/caps-clindmri/clindmri/scripts/probabilist_vertices_connectogram.py \
          -v 2 \
          -e \
          -c /i2bm/local/freesurfer/SetUpFreeSurfer.sh \
          -d /volatile/imagen/dmritest/001/processed \
          -s fs \
          -b /volatile/imagen/dmritest/001/processed/probaconnect/bedpostx.bedpostX \
          -a /volatile/imagen/dmritest/001/processed/probaconnect/t1-1mm-1-001.nii.gz \
          -n /volatile/imagen/dmritest/001/processed/probaconnect/nodif.nii.gz \
          -o /volatile/imagen/dmritest/001/processed/probaconnect/proba_vertices_connectivity \
          -i 7 \
          --hemi rh \
          --indices 1
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
    "-c", "--config", dest="fsconfig", metavar="FILE", required=True,
    help="the FreeSurfer configuration file.", type=is_file)
parser.add_argument(
    "-d", "--fsdir", dest="fsdir", required=True, metavar="PATH",
    help="the FreeSurfer processing home directory.", type=is_directory)
parser.add_argument(
    "-s", "--subjectid", dest="subjectid", required=True,
    help="the subject identifier.")
parser.add_argument(
    "-b", "--beddir", dest="bedpostxdir", required=True, metavar="PATH",
    help="the FSL bedpostx output directory.", type=is_directory)
parser.add_argument(
    "-a", "--anat", dest="anat_file", metavar="FILE", required=True,
    help="a T1 anatomical image file", type=is_file)
parser.add_argument(
    "-n", "--nodif", dest="nodif_file", metavar="FILE", required=True,
    help="the diffsuion sequence b0 image file.", type=is_file)
parser.add_argument(
    "-o", "--outdir", dest="outdir", metavar="PATH", required=True,
    help=("the output directory where the connectogram of each subject will "
          "be saved."), type=is_directory)
parser.add_argument(
    "--trf", dest="trf_file", default=None, metavar="FILE",
    help="a transformation file: diffusion space -> structural space",
    type=is_file)
parser.add_argument(
    "--dat", dest="dat_file", default=None, metavar="FILE",
    help=("a transformation file as computed by tkregister2: structural "
          "space -> FS space)."), type=is_file)
parser.add_argument(
    "-i", "--icoorder", dest="ico_order", default=7, type=int,
    choices=range(8),
    help=("specifies the order of the icosahedral tesselation (in [0, 7]) used "
          "to define the surface resolution."))
parser.add_argument(
    "--hemi", dest="hemi", default="rh", choices=["lh", "rh"],
    help="select the hemisphere to be processed.")
parser.add_argument(
    "--indices", dest="vertices_indices", default=None, nargs="+", type=int,
    help=("the vertices indices used as seeding point to compute the "
          "connectogram. If None all the hemisphere vertices are considered."))
args = parser.parse_args()


# Required arguments
bedpostxdir = args.bedpostxdir
fsdir = args.fsdir
subjectid = args.subjectid
t1_file = args.anat_file
nodif_file = args.nodif_file
outdir = args.outdir
fsconfig = args.fsconfig

# Optional arguments
trf_file = args.trf_file
dat_file = args.dat_file
ico_order = args.ico_order
hemi = args.hemi
vertices_indices = args.vertices_indices


"""
First construct the subject connectogram output directory: check
its existance on the file systema and clean it if requested.
"""
if args.verbose > 0:
    print("[info] Start probabilist connectogram...")
    print("[info] FS dir: {0}.".format(fsdir))
    print("[info] Bedpostx dir: {0}.".format(bedpostxdir))
    print("[info] Output dir: {0}.".format(outdir))
    print("[info] T1 file: {0}.".format(t1_file))
    print("[info] Nodif file: {0}.".format(nodif_file))
    print("[info] Subject: {0}.".format(subjectid))
    print("[info] Trf file: {0}.".format(trf_file))
    print("[info] Dat file: {0}.".format(dat_file))
    print("[info] Ico order: {0}.".format(ico_order))
connectdir = os.path.join(outdir, subjectid)
if not os.path.isdir(connectdir):
    os.makedirs(connectdir)
elif args.erase:
    shutil.rmtree(connectdir)
    os.makedirs(connectdir)


"""
Check that all the data are available to compute the connectogram.
"""
subjdir = os.path.join(fsdir, subjectid)
if not os.path.isdir(subjdir):
    raise ValueError(
        "'{0}' is not a FreeSurfer subject folder.".format(subjdir))
convertdir = os.path.join(subjdir, "convert")
if not os.path.isdir(convertdir):
    raise ValueError(
        "'{0}' has not been generated with the "
        "'clindmri.scripts.freesurfer_conversion' script.".format(convertdir))
whitefile = os.path.join(convertdir,
                         "{0}.white.{1}.native".format(hemi, ico_order))
if not os.path.isfile(whitefile):
    raise ValueError(
        "'{0}' is not a file. Generate it through the "
        "'clindmri.scripts.freesurfer_conversion' script.".format(whitefile))
if not os.path.isdir(bedpostxdir):
    raise ValueError(
        "'{0}' is not a folder.".format(bedpostxdir))
merged_prefix = os.path.join(bedpostxdir, "merged")
merged_files = glob.glob(merged_prefix + "*")
if len(merged_files) == 0:
    raise ValueError(
        "'{0}' is not a FSL bedpostx folder, see "
        "'clindmri.tractography.fsl.bedpostx'.".format(bedpostxdir))
nodifmask_file = os.path.join(bedpostxdir, "nodif_brain_mask.nii.gz")
if not os.path.isfile(nodifmask_file):
    raise ValueError(
        "'{0}' is not a mask file. Generate it through the "
        "'??' script.".format(nodifmask_file))


"""
If no '.dat' registration file is provided, find it in the FreeSurfer subject
directory.
"""
if dat_file is None:
    dat_file = os.path.join(subjdir, "convert", "register.native.dat")
    if not os.path.isfile(dat_file):
        raise ValueError("'{0}' has not been generated with the "
                         "'clindmri.scripts.freesurfer_conversion' "
                         "script.".format(dat_file))


""" 
If no '.trf' file is provided, register the nodif image on the t1 image
to get it.
"""
if trf_file is None:
    trf_file = os.path.join(connectdir, "dmri_to_t1.trf")
    reg_file = os.path.join(connectdir, "nodif_to_t1.nii.gz")
    flirt(nodif_file, t1_file, omat=trf_file, out=reg_file, usesqform=False,
          cost="normmi", dof=6)


""" 
Launch the tractography on the requested point of the cortical surface on the
selected hemisphere
"""
# Load the white mesh in the diffusion space
physical_to_index = numpy.linalg.inv(nibabel.load(t1_file).get_affine())
voxel_diff_to_t1 = flirt2aff(trf_file, nodif_file, t1_file)
voxel_t1_to_diff = numpy.linalg.inv(voxel_diff_to_t1)

# Get deformation between the ras and ras-tkregister spaces
surface = TriSurface.load(whitefile)
asegfile = os.path.join(segfile, "aseg.mgz")
translation = tkregister_translation(asegfile, fsconfig)
deformation = numpy.dot(voxel_t1_to_dest, numpy.dot(physical_to_index, translation))
surface.vertices = apply_affine_on_mesh(surface.vertices, deformation)



# Select the vertices of interest
if vertices_indices is None:
    vertices_indices = range(len(surface.vertices))

# Go through all the hemisphere vertices
textures = {}
for index in vertices_indices:

    # Select the seeding vertex
    point = surface.vertices[index]

    # Create a directory for each seeding vertex in order to avoid collision
    # Raise an exception if the directory has already been created
    vertexdir = os.path.join(connectdir, "{0}_{1}".format(hemi, index))
    os.mkdir(vertexdir)

    # Write seeding vertex coordinates to file
    seed_file = os.path.join(vertexdir, "fdt_coordinates.txt")
    with open(seed_file, "w") as open_file:
        for coordinate in point:
            print(coordinate, file=open_file)

    # Launch the tractography on this seeding vertex
    textures[repr(point.tolist())] = get_profile(
        ico_order, nodif_file, nodifmask_file, seed_file, merged_prefix,
        vertexdir, t1_file, trf_file, dat_file, fsdir, subjectid, fsconfig)

