#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
import pprint
import subprocess

from .exceptions import ConnectomistError, ConnectomistRuntimeError

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
    Function that wraps the PtkNiftiToGisConverter command line tool from
    Connectomist to make it capsul-isable.

    Parameters
    ----------
    path_nifti: Str, path to the input Nifti file to be converted.
    path_gis:   Str, path without extension to the 3 output GIS files.
    nb_tries:   Int, number of times to try the conversion.
                The reason why it sometimes fails is not known.

    <process>
        <return name="path_gis_ima"  type="File" desc="Path to output Gis .ima file."/>
        <return name="path_gis_minf" type="File" desc="Path to output Gis .minf file."/>
        <return name="path_gis_dim"  type="File" desc="Path to output Gis .dim file."/>
        <input  name="path_nifti"    type="File" desc="Path to input Nifti file."/>
        <input  name="path_gis"      type="Str"  desc="Path without extension to
                                                       the 3 Gis files."/>
        <input  name="nb_tries"      type="Int"  optional="True"/>
    </process>
    """

    cmd = ["PtkNifti2GisConverter", "-i", path_nifti, "-o", path_gis, 
           "-verbose", "False", "-verbosePluginLoading", "False"]
    
    path_gis_ima  = path_gis + ".ima"
    path_gis_minf = path_gis + ".ima.minf"
    path_gis_dim  = path_gis + ".dim"
    
    nb_tried = 0
    while nb_tried < nb_tries:
        subprocess.check_call(cmd)
        nb_tried += 1
        if os.path.isfile(path_gis_ima) and os.path.isfile(path_gis_minf):
            return path_gis_ima, path_gis_minf, path_gis_dim

    raise ConnectomistError("Conversion nifti_to_gis failed, cmd:\n%s" % " ".join(cmd))


def gis_to_nifti(path_gis, path_nifti, nb_tries=3):
    """
    Function that wraps the PtkGis2NiftiConverter command line tool from
    Connectomist to make it capsul-isable.

    Parameters
    ----------
    path_gis:   Str, path without extension to the 3 input GIS files.
    path_nifti: Str, path to the output Nifti file.
    nb_tries:   Int, number of times to try the conversion.
                The reason why it sometimes fails is not known.

    Returns
    -------
    path_nifti: path to output file.

    <process>
        <return name="path_nifti" type="File"  desc="Path to output Nifti file."/>
        <input  name="path_gis"   type="Str"   desc="Path without extension to
                                                     the 3 input Gis files."/>
        <input  name="path_nifti"  type="File" desc="Path to output Nifti file."/>
        <input  name="nb_tries"    type="Int"  optional="True"/>
    </process>
    """

    # Add .nii extension if it's not the case
    if not (path_nifti.endswith(".nii") or path_nifti.endswith(".nii.gz")):
        path_nifti += ".nii"

    # Call command line tool:
    # It creates a Nifti + a .minf file (metainformation)
    cmd = ["PtkGis2NiftiConverter", "-i", path_gis, "-o", path_nifti, 
           "-verbose", "False", "-verbosePluginLoading", "False"]
    
    nb_tried = 0
    while nb_tried < nb_tries:
        subprocess.check_call(cmd)
        nb_tried += 1
        if os.path.isfile(path_nifti):
            return path_nifti

    raise ConnectomistError("Conversion gis_to_nifti failed for %s" % path_gis)


def concatenate_volumes(path_inputs, path_output, axis="t", nb_tries=3):
    """
    Function that wraps the PtkCat command line tool from Connectomist to make
    it capsul-isable, with only the basic arguments. It allows concatenating
    volumes. In particular to concatenante T2 and DW volumes in one file at the
    end of the preprocessing.

    Parameters
    ----------
    path_input:  List of str, paths to input volumes.
    path_output: Str, path to output volume.
    axis:        Str, axis along which the concatenation is done.
    nb_tries:    Int, number of times to try the concatenation.
                 The reason why it sometimes fails is not known.

    <process>
        <return name="path_output" type="File"/>
        <input  name="path_inputs" type="List"/>
        <input  name="path_output" type="File"/>
        <input  name="axis"        type="Str" optional="True"/>
        <input  name="nb_tries"    type="Int" optional="True"/>
    </process>
    """
    cmd = ["PtkCat", "-i"] + path_inputs + ["-o", path_output, "-t", axis,
           "-verbose", "False", "-verbosePluginLoading", "False"]
    nb_tried = 0
    while nb_tried < nb_tries:
        subprocess.check_call(cmd)
        nb_tried += 1
        if os.path.isfile(path_output + ".ima"):
            return path_output

    raise ConnectomistError("Concatenate_volumes failed, cmd:\n%s" % " ".join(cmd))
