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

from .exceptions import ConnectomistError, ConnectomistRuntimeError, BadFileError

###############################################################################
# UTILITY FUNCTIONS
###############################################################################


def create_parameter_file(algorithm_name, parameters_value, output_directory):
    """
    Writes the .py file that Connectomist uses when working in command line.

    Parameters
    ----------
    algorithm_name:   Str, name of Connectomist's tab.
    parameters_value: Dict, parameter values for the tab.
    output_directory: Str, path to directory where to write the parameter file.
                      If not existing the directory is created.

    Returns
    -------
    parameter_file:   Str, path to the created parameter file.
    """

    # If not existing create output_directory
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    parameter_file = os.path.join(output_directory, "%s.py" % algorithm_name)

    with open(parameter_file, "w") as f:
        f.write("algorithmName = '%s'\n" % algorithm_name)

        # Pretty text to write, without the first "{"
        pretty_dict = pprint.pformat(parameters_value)[1:]
        f.write("parameterValues = {\n " + pretty_dict)

    return parameter_file


def run_connectomist(algorithm_name, parameter_file, connectomist_path=None):
    """
    Makes the command line call to run a specific tab from Connectomist.

    Parameters
    ----------
    algorithm_name:    Str, name of Connectomist's tab.
    paramter_file:     Str, path to the parameter file for the tab.
    connectomist_path: Str (optional), path to the Connectomist executable.

    Raises
    ------
    ConnectomistRuntimeError: If Connectomist call failed.
    """

    if not connectomist_path:
        connectomist_path = "/i2bm/local/Ubuntu-14.04-x86_64/ptk/bin/connectomist"

    cmd = "%s -p %s -f %s" % (connectomist_path, algorithm_name, parameter_file)
    try:
        subprocess.check_call(cmd, shell=True)
    except:
        raise ConnectomistRuntimeError(algorithm_name, parameter_file)


###############################################################################
# Wrappers to Ptk command line tools, may disappear in the future.
###############################################################################

def nifti_to_gis(path_nifti, path_gis, nb_tries=3):
    """
    Function that wraps the PtkNifti2GisConverter command line tool from
    Connectomist to make it capsul-isable.

    Parameters
    ----------
    path_nifti: Str, path to the input Nifti file to be converted.
    path_gis:   Str, path without extension to the 3 output GIS files.
    nb_tries:   Int, number of times to try the conversion.
                The reason why it sometimes fails is not known.

    Raises
    ------
    ConnectomistError: If call to PtkNifti2GisConverter failed nb_tries times.

    <unit>
        <output name="path_gis" type="File" description="Path to output
            .ima file (Gis format)."/>

        <input name="path_nifti" type="File" description="Path to input Nifti file."/>
        <input name="path_gis" type="Str" description="Path to the Gis .ima file.
            If no extension, .ima is added."/>
        <input name="nb_tries" type="Int"/>
    </unit>
    """
    # Check input existence
    if not os.path.isfile(path_nifti):
        raise BadFileError(path_nifti)
    
    # Add extension if there is none
    if not path_gis.endswith("ima"):
        path_gis += ".ima"
    
    # Call command line tool
    cmd = ["PtkNifti2GisConverter", "-i", path_nifti, "-o", path_gis,
           "-verbose", "False", "-verbosePluginLoading", "False"]
    nb_tried = 0
    while nb_tried < nb_tries:
        subprocess.check_call(cmd)
        nb_tried += 1
        if os.path.isfile(path_gis):
            return path_gis
        else:
            time.sleep(10)  # wait 10 sec before retrying

    raise ConnectomistError("Conversion nifti_to_gis failed, cmd:\n%s" % " ".join(cmd))


def gis_to_nifti(path_gis, path_nifti, nb_tries=3):
    """
    Function that wraps the PtkGis2NiftiConverter command line tool from
    Connectomist to make it capsul-isable.

    Parameters
    ----------
    path_gis:   Str, path to the Gis .ima file.
    path_nifti: Str, path to the output Nifti file.
    nb_tries:   Int, number of times to try the conversion.
                The reason why it sometimes fails is not known.

    Returns
    -------
    path_nifti: path to output file.

    Raises
    ------
    ConnectomistError: If call to PtkGis2NiftiConverter failed nb_tries times.

    <unit>
        <output name="path_nifti" type="File" />

        <input name="path_gis"   type="Str" description="Path to the Gis .ima file."/>
        <input name="path_nifti" type="File"  />
        <input name="nb_tries"   type="Int"   />
    </unit>
    """
    # Check input existence
    if not os.path.isfile(path_gis):
        raise BadFileError(path_gis)

    # The command line tool does not handle .gz properly
    if path_nifti.endswith(".gz"):
        path_nifti = path_nifti[:-3] # remove ".gz"
    
    # Add .nii extension if not the case
    if not path_nifti.endswith(".nii"):
        path_nifti += ".nii"

    # Call command line tool:
    # it creates a Nifti + a .minf file (metainformation)
    cmd = ["PtkGis2NiftiConverter", "-i", path_gis, "-o", path_nifti,
           "-verbose", "False", "-verbosePluginLoading", "False"]

    nb_tried = 0
    while nb_tried < nb_tries:
        subprocess.check_call(cmd)
        nb_tried += 1
        if os.path.isfile(path_nifti):
            return path_nifti
        else:
            time.sleep(10)  # wait 10 sec before retrying

    raise ConnectomistError("Conversion gis_to_nifti failed, cmd:\n%s" % " ".join(cmd))


def concatenate_volumes(path_inputs, path_output, axis="t", nb_tries=3):
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
            time.sleep(10)  # wait 10 sec before retrying

    raise ConnectomistError("Failed to concatenate volumes, cmd:\n%s" % " ".join(cmd))
