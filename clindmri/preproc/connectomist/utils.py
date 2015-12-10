#!/usr/bin/env python
# -*- coding: utf-8 -*-

##########################################################################
# NSAp - Copyright (C) CEA, 2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

import os
import time
import pprint
import subprocess
import gzip
import shutil

from .exceptions import ConnectomistError, BadFileError

###############################################################################
# UTILITY FUNCTIONS
###############################################################################


def create_parameter_file(algorithm, parameters_dict, output_directory):
    """
    Writes the .py file that Connectomist uses when working in command line.

    Parameters
    ----------
    algorithm:        Str, name of Connectomist's tab.
    parameters_dict:  Dict, parameter values for the tab.
    output_directory: Str, path to directory where to write the parameter file.
                      If not existing the directory is created.

    Returns
    -------
    parameter_file:   Str, path to the created parameter file.
    """

    # If not existing create output_directory
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    parameter_file = os.path.join(output_directory, "%s.py" % algorithm)

    with open(parameter_file, "w") as f:
        f.write("algorithmName = '%s'\n" % algorithm)

        # Pretty text to write, without the first "{"
        pretty_dict = pprint.pformat(parameters_dict)[1:]
        f.write("parameterValues = {\n " + pretty_dict)

    return parameter_file


def run_connectomist(algorithm, parameter_file, output_directory,
                     path_connectomist=None, nb_tries=10):
    """
    Makes the command line call to run a specific tab from Connectomist.

    Parameters
    ----------
    algorithm:         Str, name of Connectomist's tab.
    paramter_file:     Str, path to the parameter file for the tab.
    output_directory:  Str, path to directory where the algorithm outputs.
    path_connectomist: Str (optional), path to the Connectomist executable.
    nb_tries:          Int, nb of times to try an algorithm if it fails.
                       It often crashes when running in parallel. The reason
                       why it crashes is unknown.

    Raises
    ------
    ConnectomistError: If Connectomist call failed.
    """

    if not path_connectomist:
        path_connectomist = "/i2bm/local/Ubuntu-14.04-x86_64/ptk/bin/connectomist"

    # Map algorithm name to a list of files that should be created.
    # It is meant to check that the call to Connectomist worked, because the
    # exit code is 0 even when it fails.
    files_to_check = {
        "DWI-Data-Import-And-QSpace-Sampling":
            ["t2.ima", "dw.ima", "acquisition_parameters.py"],
        "DWI-Rough-Mask-Extraction": ["mask.ima"],
        "DWI-Outlier-Detection": ["t2_wo_outlier.ima", "dw_wo_outlier.ima"],
        "DWI-Susceptibility-Artifact-Correction":
            ["t2_wo_susceptibility.ima", "dw_wo_susceptibility.ima"],
        "DWI-Eddy-Current-And-Motion-Correction":
            ["t2_wo_eddy_current_and_motion.ima", "dw_wo_eddy_current_and_motion.ima"]
    }

    # List filenames to be checked for this call
    to_check = [os.path.join(output_directory, f) for f in files_to_check.get(algorithm, [])]

    # Command to be run
    cmd = "%s -p %s -f %s" % (path_connectomist, algorithm, parameter_file)

    # Run the command, multiple times if it fails.
    nb_tried = 0
    while nb_tried < nb_tries:
        nb_tried += 1

        try:
            subprocess.check_call(cmd, shell=True)
        except:
            pass

        success = all(map(os.path.isfile, to_check))
        if success:
            return
        else:
            time.sleep(3)  # wait 3 sec before retrying

    # If the command failed nb_tries times, raise an exception
    raise ConnectomistError("Connectomist call failed with cmd:\n%s" % cmd)


###############################################################################
# Wrappers to Ptk command line tools, may disappear in the future.
###############################################################################

def nifti_to_gis(nifti, gis, nb_tries=10):
    """
    Function that wraps the PtkNifti2GisConverter command line tool from
    Connectomist to make it capsul-isable.

    Parameters
    ----------
    nifti:    Str, path to the input Nifti file to be converted.
    gis:      Str, path without extension to the 3 output GIS files.
    nb_tries: Int, number of times to try the conversion. It often fails
              when using parallel processing.

    Raises
    ------
    ConnectomistError: If call to PtkNifti2GisConverter failed nb_tries times.

    <unit>
        <output name="path_gis" type="File" description="Path to output
            .ima file (Gis format)."/>

        <input name="nifti" type="File" description="Path to input Nifti file."/>
        <input name="gis" type="Str" description="Path to the Gis .ima file.
            If no extension, .ima is added."/>
        <input name="nb_tries" type="Int"/>
    </unit>
    """
    # Check input existence
    if not os.path.isfile(nifti):
        raise BadFileError(nifti)

    # Add extension if there is none
    if not gis.endswith("ima"):
        gis += ".ima"

    # Call command line tool
    cmd = ["PtkNifti2GisConverter", "-i", nifti, "-o", gis,
           "-verbose", "False", "-verbosePluginLoading", "False"]
    nb_tried = 0
    while nb_tried < nb_tries:
        subprocess.check_call(cmd)
        nb_tried += 1
        if os.path.isfile(gis):
            return gis
        else:
            time.sleep(3)  # wait 10 sec before retrying

    raise ConnectomistError("Conversion nifti_to_gis failed, cmd:\n%s" % " ".join(cmd))


def gz_compress(file_to_compress, clean=True):
    """
    Compress a file with gzip, the output path is the same but with ".gz" extension.
    If 'clean' is True, the input file is deleted, to keep only the compressed
    version.

    The function raises ValueError if the input file does not exist and
    Exception if the output compressed file is not created.
    """
    if not os.path.isfile(file_to_compress):
        raise ValueError("File to compress does not exist: {}".format(file_to_compress))

    gz_file = file_to_compress + ".gz"
    with open(file_to_compress, 'rb') as f_in, gzip.open(gz_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    if clean:
        os.remove(file_to_compress)

    if not os.path.isfile(gz_file):
        raise Exception("Failed to compress %s to %s" % (file_to_compress, gz_file))
    else:
        return gz_file


def gis_to_nifti(gis, nifti, nb_tries=10):
    """
    Function that wraps the PtkGis2NiftiConverter command line tool from
    Connectomist to make it capsul-isable.

    Parameters
    ----------
    gis:      Str, path to the Gis .ima file.
    nifti:    Str, path to the output Nifti file.
    nb_tries: Int, number of times to try the conversion. It often fails
              when using parallel processing.

    Returns
    -------
    nifti: path to output file.

    Raises
    ------
    ConnectomistError: If call to PtkGis2NiftiConverter failed nb_tries times.

    <unit>
        <output name="nifti" type="File" />

        <input name="gis"      type="Str"  />
        <input name="nifti"    type="File" />
        <input name="nb_tries" type="Int"  />
    </unit>
    """
    # Check input existence
    if not os.path.isfile(gis):
        raise BadFileError(gis)

    # The command line tool does not handle .gz properly
    compress_to_gz = False
    if nifti.endswith(".gz"):
        nifti = nifti[:-3] # remove ".gz"
        compress_to_gz = True

    # Add .nii extension if not the case
    if not nifti.endswith(".nii"):
        nifti += ".nii"

    # Call command line tool:
    # it creates a Nifti + a .minf file (metainformation)
    cmd = ["PtkGis2NiftiConverter", "-i", gis, "-o", nifti,
           "-verbose", "False", "-verbosePluginLoading", "False"]

    nb_tried = 0
    while nb_tried < nb_tries:
        subprocess.check_call(cmd)
        nb_tried += 1
        if os.path.isfile(nifti):
            if compress_to_gz:
                nifti = compress_to_gz(nifti)
            return nifti
        else:
            time.sleep(3)  # wait 3 sec before retrying

    raise ConnectomistError("Conversion gis_to_nifti failed, cmd:\n%s" % " ".join(cmd))


def concatenate_volumes(path_inputs, path_output, axis="t", nb_tries=10):
    """
    Function that wraps the PtkCat command line tool from Connectomist, with
    only the basic arguments. It allows concatenating volumes. In particular to
    concatenante T2 and DW volumes in one file at the end of the preprocessing.

    Parameters
    ----------
    path_inputs: List of str, paths to input volumes.
    path_output: Str, path to output volume.
    axis:        Str, axis along which the concatenation is done.
    nb_tries:    Int, number of times to try the concatenation.
                 The reason why it sometimes fails is not known.

    Raises
    ------
    ConnectomistError: If call to PtkCat failed
    """
    # check input existence
    for path in path_inputs:
        if not os.path.isfile(path):
            raise BadFileError(path)

    # Run command line tool
    cmd = ["PtkCat", "-i"] + path_inputs + ["-o", path_output, "-t", axis,
           "-verbose", "False", "-verbosePluginLoading", "False"]
    nb_tried = 0
    while nb_tried < nb_tries:
        subprocess.check_call(cmd)
        nb_tried += 1
        if os.path.isfile(path_output):
            return path_output
        else:
            time.sleep(3)  # wait 3 sec before retrying

    raise ConnectomistError("Failed to concatenate volumes, cmd:\n%s" % " ".join(cmd))
