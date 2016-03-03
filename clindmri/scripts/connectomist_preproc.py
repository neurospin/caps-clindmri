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
try:
    import bredala
    bredala.USE_PROFILER = False
    bredala.register(
        "clindmri.preproc.connectomist.complete_preprocessing",
        names=["complete_preprocessing"])
    bredala.register(
        "clindmri.preproc.connectomist.import_and_qspace_model",
        names=["dwi_data_import_and_qspace_sampling"])
    bredala.register(
        "clindmri.preproc.connectomist.mask",
        names=["dwi_rough_mask_extraction"])
    bredala.register(
        "clindmri.preproc.connectomist.outliers",
        names=["dwi_outlier_detection"])
    bredala.register(
        "clindmri.preproc.connectomist.susceptibility",
        names=["dwi_susceptibility_artifact_correction"])
    bredala.register(
        "clindmri.preproc.connectomist.registration",
        names=["dwi_to_anatomy"])
    bredala.register(
        "clindmri.preproc.connectomist.eddy_current_and_motion",
        names=["dwi_eddy_current_and_motion_correction",
               "export_eddy_motion_results_to_nifti"])
except:
    pass

# Clindmri import
from clindmri.preproc.connectomist.complete_preprocessing import (
    complete_preprocessing)
from clindmri.extensions.connectomist.manufacturers import MANUFACTURERS


# Parameters to keep trace
__hopla__ = ["outdir", "subjectid", "subjdir", "dwi", "bvec", "bval",
             "manufacturer", "delta_te", "partial_fourier_factor",
             "parallel_acceleration_factor", "b0_magnitude", "b0_phase",
             "b0_magnitude", "invertx", "inverty", "invertz", "negative_sign",
             "echo_spacing", "epi_factor", "b0_field", "water_fat_shift",
             "delete_steps", "morphologist_dir"]

# Script documentation
doc = """
Connectomist preproc
~~~~~~~~~~~~~~~~~~~~

Function that runs all preprocessing tabs from Connectomist.

Steps:

1- Create the preprocessing output directory if not existing.
2- Import files to Connectomist and choose q-space model.
3- Create a brain mask.
4- Detect and correct outlying diffusion slices.
5- Susceptibility correction.
6- Eddy current and motion correction.
7- Export result as a Nifti with a .bval and a .bvec.
8- Export outliers.
9- Delete intermediate files and directories if requested.

Command:

python $HOME/git/caps-clindmri/clindmri/scripts/connectomist_preproc.py \
    -v 2 \
    -o /tmp/morphologist/connectomist \
    -s ab130187 \
    -i /tmp/morphologist/ab130187/000009_DTI/DTI.nii.gz \
    -r /tmp/morphologist/ab130187/000009_DTI/DTI.bvec \
    -b /tmp/morphologist/ab130187/000009_DTI/DTI.bval \
    -m Siemens \
    -t 2.46 \
    -f 0.75 \
    -a 2 \
    -u /tmp/morphologist/ab130187/000015_B0MAP/B0MAP.nii.gz \
    -p /tmp/morphologist/ab130187/000016_B0MAP/B0MAP.nii.gz \
    -x \
    -c 0.75 \
    -l 3.0 \
    -g /neurospin/senior/nsap/data/V0/morphologist \
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
    "-o", "--outdir", dest="outdir", required=True, metavar="PATH",
    help="the Connectomist preprocessing home directory.", type=is_directory)
parser.add_argument(
    "-s", "--subjectid", dest="subjectid", required=True,
    help="the subject identifier.")
parser.add_argument(
    "-i", "--dwifile", dest="dwifile", metavar="FILE", required=True,
    help="the subject diffusion image to be processed.", type=is_file)
parser.add_argument(
    "-r", "--bvecfile", dest="bvecfile", metavar="FILE", required=True,
    help="the subject diffusion bvec to be processed.", type=is_file)
parser.add_argument(
    "-b", "--bvalfile", dest="bvalfile", metavar="FILE", required=True,
    help="the subject diffusion bval to be processed.", type=is_file)
parser.add_argument(
    "-m", "--manufacturer", dest="manufacturer", choices=MANUFACTURERS.keys(),
    help="the MRI scanner manufacturer.")
parser.add_argument(
    "-t", "--delta_te", dest="delta_te", type=float, required=True,
    help=("difference in msec between the 2 echoes in B0 magnitude map "
          "acquisition."))
parser.add_argument(
    "-f", "--partial_fourier_factor", dest="partial_fourier_factor",
    type=float, required=True,
    help="percentage of k-space plane acquired (]0;1]).")
parser.add_argument(
    "-a", "--parallel_acceleration_factor",
    dest="parallel_acceleration_factor", type=float, required=True,
    help="nb of parallel acquisition in k-space plane.")
parser.add_argument(
    "-u", "--b0_magnitude", dest="b0_magnitude", metavar="FILE", required=True,
    help="path to B0 magnitude map, also contains phase for GE.", type=is_file)
parser.add_argument(
    "-p", "--b0_phase", dest="b0_phase", metavar="FILE",
    help="not for GE, path to B0 phase map.", type=is_file)
parser.add_argument(
    "-x", "--invertx", dest="invertx", action="store_true",
    help="if True invert x-axis in diffusion model.")
parser.add_argument(
    "-y", "--inverty", dest="inverty", action="store_true",
    help="if True invert y-axis in diffusion model.")
parser.add_argument(
    "-z", "--invertz", dest="invertz", action="store_true",
    help="if True invert z-axis in diffusion model.")
parser.add_argument(
    "-n", "--negative_sign", dest="negative_sign",  action="store_true",
    help=("if True invert direction of unwarping in usceptibility-distortion"
          "correction."))
parser.add_argument(
    "-c", "--echo_spacing", dest="echo_spacing", type=float, required=False,
    help=("not for Philips, acquisition time in ms between 2 centers of 2 "
          "consecutively acquired lines in k-space."))
parser.add_argument(
    "-k", "--epi_factor", dest="epi_factor", type=int, required=False,
    help="nb of echoes after one RF pulse, i.e. echo train length.")
parser.add_argument(
    "-l", "--b0_field", dest="b0_field", type=float, default=3.0,
    help="Philips only, B0 field intensity.")
parser.add_argument(
    "-w", "--water_fat_shift", dest="water_fat_shift", type=float,
    default=4.68, help="Philips only, the spins frequency offsets.")
parser.add_argument(
    "-d", "--delete_steps", dest="delete_steps", action="store_true",
    help="if True remove all intermediate files and directories at the end.")
parser.add_argument(
    "-g", "--morphologist_dir", dest="morphologist_dir", required=True,
    metavar="PATH", help="the path to the morphologist processings.",
    type=is_directory)
args = parser.parse_args()


"""
First check if the subject Connectomist directory exists on the file system,
and clean it if requested.
"""
if args.verbose > 0:
    print("[info] Start Connectomist preproc...")
    print("[info] Directory: {0}.".format(args.outdir))
    print("[info] Subject: {0}.".format(args.subjectid))
    print("[info] DWI: {0}.".format(args.dwifile))
    print("[info] BVEC: {0}.".format(args.bvecfile))
    print("[info] BVAL: {0}.".format(args.bvalfile))
    print("[info] Magnitude: {0}.".format(args.b0_magnitude))
    print("[info] Phase: {0}.".format(args.b0_phase))
outdir = args.outdir
subjectid = args.subjectid
subjdir = os.path.join(args.outdir, subjectid, "preproc")
dwi = args.dwifile
bvec = args.bvecfile
bval = args.bvalfile
manufacturer = args.manufacturer
delta_te = args.delta_te
partial_fourier_factor = args.partial_fourier_factor
parallel_acceleration_factor = args.parallel_acceleration_factor
b0_magnitude = args.b0_magnitude
b0_phase = args.b0_phase
invertx = args.invertx
inverty = args.inverty
invertz = args.invertz
negative_sign = args.negative_sign
echo_spacing = args.echo_spacing
epi_factor = args.epi_factor
b0_field = args.b0_field
water_fat_shift = args.water_fat_shift
delete_steps = args.delete_steps
morphologist_dir = args.morphologist_dir
if not os.path.isdir(subjdir):
    os.makedirs(subjdir)
elif os.path.isdir(subjdir) and args.erase:
    shutil.rmtree(subjdir)
    os.mkdir(subjdir)

"""
Connectomist preproc: all steps
"""
preproc_dwi, preproc_bval, preproc_bvec = complete_preprocessing(
    subjdir,
    subjectid,
    dwi,
    bval,
    bvec,
    manufacturer,
    delta_te,
    partial_fourier_factor,
    parallel_acceleration_factor,
    b0_magnitude,
    b0_phase=b0_phase,
    invertX=invertx,
    invertY=inverty,
    invertZ=invertz,
    negative_sign=negative_sign,
    echo_spacing=echo_spacing,
    EPI_factor=epi_factor,
    b0_field=b0_field,
    water_fat_shift=water_fat_shift,
    delete_steps=delete_steps,
    morphologist_dir=morphologist_dir,
    nb_tries=10,
    path_connectomist=(
        "/i2bm/local/Ubuntu-14.04-x86_64/ptk/bin/connectomist"))
if args.verbose > 1:
    print("[result] In folder: {0}.".format(subjdir))
