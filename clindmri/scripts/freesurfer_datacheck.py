#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013-2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System modules
from __future__ import print_function
import os
import re
import argparse


doc = """
FreeSurfer Data Check
---------------------

Run this command to check if your input FreeSurfer processing home
directory ('fsdir') contains the expected FreeSurfer output files for all
the subjects.
It is possible to check in details one specific subject specifying his
identifier.
This script can also checks the 'freesurfer_conversion' outputs.

python $HOME/git/caps-clindmri/clindmri/scripts/freesurfer_datacheck.py \
    -v 2 \
    -d /neurospin/imagen/BL/processed/freesurfer \
    -r '\d{12}'
"""


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
    "-d", "--fsdir", dest="fsdir", required=True, metavar="PATH",
    help="the FreeSurfer processing home directory.", type=is_directory)
parser.add_argument(
    "-s", "--subjectid", dest="subjectid",
    help=("the subject identifier. If specified the complete subject "
          "structure is exposed."))
parser.add_argument(
    "-c", "--conversion", dest="conversion_count", type=int,
    help="the number of files expected in the converts folder.")
parser.add_argument(
    "-q", "--qc", dest="qc_count", type=int,
    help="the number of files expected in the qc folder.")
parser.add_argument(
    "-r", "--regex", dest="regex", required=True,
    help=("the expression to be applied on the 'fsdir' subfolders in order to"
          "retrieve the subjects only."), type=str)
args = parser.parse_args()


"""
Expected FreeSurfer data structure
----------------------------------

First we define the expected organization of the FreeSurfer output files: we
consider the files count only as a criterion.
"""

fsstruct = {
    "bem": 0,
    "label": 69,
    "mri": 35,
    os.path.join("mri", "orig"): 1,
    os.path.join("mri", "transforms"): 13,
    os.path.join("mri", "transforms", "bak"): 0,
    "scripts": 11,
    "src": 0,
    "stats": 18,
    "surf": 70,
    "tmp": 0,
    "touch": 67,
    "trash": 0,
    "extrapaths": [""]
}
if args.conversion_count is not None:
    fsstruct["convert"] = args.conversion_count
if args.qc_count is not None:
    fsstruct[os.path.join("convert", "qc")] = args.qc_count


"""
Browsing
--------

Detect each subject in the FreeSurfer processing home directory 'fsdir' and
then check the subject processing tree.
"""
status = {}
if args.subjectid is not None:
    fsdir_subfolders = [args.subjectid]
else:
    fsdir_subfolders = os.listdir(args.fsdir)
for sid in fsdir_subfolders:

    # Check if the subject id is valid using the input regex
    match = re.findall(args.regex, sid)
    if len(match) != 1 or match[0] != sid:
        continue

    # Store the subject tree folder-folder files counts
    status[sid] = {"extrapaths": []}
    siddir = os.path.join(args.fsdir, sid)
    for path, dirs, files in os.walk(siddir):
        rpath = path.replace(siddir, "").lstrip(os.sep)
        if rpath in fsstruct:
            status[sid][rpath] = len(files)
        else:
            status[sid]["extrapaths"].append(rpath)


"""
Display a summary
-----------------

Depending on the verbosity, display the number of complete/failed processings,
or the list of subjects that didn't run properly.
"""
failed_sids = []
total_success = 0
total_failed = 0
total = len(status)
for sid, sid_status in status.items():
    if sid_status == fsstruct:
        total_success += 1
    else:
        total_failed += 1
        failed_sids.append(sid)
if args.verbose >= 0:
    print("SUCCESS: {0}/{1}".format(total_success, total))
    print("FAILED: {0}/{1}".format(total_failed, total))
if args.verbose >= 1:
    print("FAILED SIDS: {0}".format(failed_sids))
if args.subjectid is not None:
    print("TREE '{0}': ".format(args.subjectid))
    for key, value in status[args.subjectid].items():
        print("{0}: {1} (observed) - {2} (reference)".format(key, value,
                                                             fsstruct[key]))
