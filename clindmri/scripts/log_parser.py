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
Log parser
~~~~~~~~~~

A script to browse an execution log produced when treating large
datasets. It enables us to detect the failing jobs and an attached identifier
(called also values) from a regular expression.

Command
-------

python $HOME/git/caps-clindmri/clindmri/scripts/log_parser.py \
    -v 2 \
    -l /neurospin/imagen/BL/processed/freesurfer_logs/fsreconall_BL_3.log \
    -e exitcode \
    -s cmd \
    -r '\d{12}'

python $HOME/git/caps-clindmri/clindmri/scripts/log_parser.py \
    -v 2 \
    -l /neurospin/imagen/BL/processed/freesurfer_logs/fsreconall_BL_1.log \
    -e exitcode \
    -s inputs \
    -r '\d{12}'
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
    "-l", "--log", dest="logfile", metavar="FILE", required=True,
    help="a processing log report.", type=is_file)
parser.add_argument(
    "-e", "--exit", dest="exitcode", required=True,
    help="the name of the exitcode field.", type=str)
parser.add_argument(
    "-s", "--selection", dest="selection", required=True,
    help="the name of the field that contain the data of interest.", type=str)
parser.add_argument(
    "-r", "--regex", dest="regex", required=True,
    help=("the expression to be applied on the selected field value in order "
          "retrieve the value(s) of interest."), type=str)
args = parser.parse_args()

"""
Description
-----------

Use a two browsing document strategy:
* locate the exicode field.
* locate the fields of interest and the associated values.

Then display a summary.
"""

# Define the field of interest
exitcode_field = ".{0} = ".format(args.exitcode)
selection_field = ".{0} = ".format(args.selection)

# Parse the log file
with open(args.logfile, "r") as open_file:
    lines = open_file.readlines()

# Create a status parameter that links a job name with its exitcode and desired
# values
status = {}

"""
Step 1
------

Fist browse the job names and associated status (0 or 1).
"""
for line in lines:
    if exitcode_field in line:
        first_part, second_part = line.split(exitcode_field)
        job_name = first_part.split()[-1]
        if job_name in status:
            raise ValueError("A job with name '{0}' has already been "
                             "discovered.".format(job_name))
        exitcode = second_part[0]
        try:
            exitcode = int(exitcode)
            if exitcode not in (0, 1):
                raise
        except:
            raise ValueError("'{0}' is not a valid exit code. Expected values "
                             "are in (0, 1).".format(exitcode))
        status[job_name] = {"exitcode": exitcode, "value": None}


"""
Step 2
------

Then browse the desired jobs information.
"""
for line in lines:
    if selection_field in line:
        first_part, second_part = line.split(selection_field)
        job_name = first_part.split()[-1]
        if job_name not in status:
            raise ValueError("No job with name '{0}' has been  detected "
                             "during the first browsing.".format(job_name))
        values = set(re.findall(args.regex, second_part))
        if len(values) == 1:
            value = values.pop()
        else:
            value = " - ".join(values)
        status[job_name]["value"] = value


"""
Step 3
------

Display a summary.
"""
failed_values = []
failed_job_names = []
total_success = 0
total_failed = 0
total = len(status)
for job_name, job_status in status.items():
    if job_status["exitcode"] == 0:
        total_success += 1
    else:
        total_failed += 1
        failed_job_names.append(job_name)
        failed_values.append(job_status["value"])
if args.verbose >= 0:
    print("SUCCESS: {0}/{1}".format(total_success, total))
    print("FAILED: {0}/{1}".format(total_failed, total))
if args.verbose >= 1:
    print("FAILED VALUES: {0}".format(failed_values))
if args.verbose >= 2:
    print("FAILED JOBS: {0}".format(failed_job_names))
