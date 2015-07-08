#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" Module that defines functions to wrap Connectomist's tabs and simplify
    scripting.
"""

import os
import pprint
import numpy as np
import subprocess
import tempfile
import shutil
from multiprocessing import Pool

###############################################################################
# MODULE VARIABLES
###############################################################################

# Note: we use simplified manufacturer names to simplify usage.
# e.g. we use "GE" where Connectomist uses "GE HealthCare",
# but a DICOM header would use "GE Medical Systems"
MANUFACTURERS = {"Bruker" : 0,
                 "GE"     : 1,
                 "Philips": 2,
                 "Siemens": 3}


###############################################################################
# MODULE EXCEPTION
###############################################################################

class ConnectomistRuntimeError(Exception):
    """ New exception type to differentiate Connectomist's runtime error from
        other exceptions. No specific code.
    """
    pass

###############################################################################
# UTILITY FUNCTIONS
###############################################################################


def create_parameter_file(algorithm_name, parameter_values, output_directory):
    """
    Writes the .py file that Connectomist uses when working in command line.

    Parameters
    ----------
    algorithm_name:   Str, name of Connectomist's tab.
    parameter_values: Dict, parameter values for the tab.
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
        pretty_dict = pprint.pformat(parameter_values)[1:]
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
    ConnectomistRuntimeError: if Connectomist's exit status != 0
    """

    if not connectomist_path:
        connectomist_path = "/i2bm/local/Ubuntu-14.04-x86_64/ptk/bin/connectomist"

    cmd = "%s -p %s -f %s" % (connectomist_path, algorithm_name, parameter_file)
    exit_status = os.system(cmd)

    if exit_status:
        raise ConnectomistRuntimeError("Connectomist call for %s failed. "
                                       "Exit status %d" % (algorithm_name, exit_status))



def nifti_to_gis(path_nifti, path_gis, verbose=False):
    """
    Function that wraps the PtkNiftiToGisConverter command line tool from
    Connectomist to make it capsul-isable.

    Parameters
    ----------
    path_nifti: Str, path to the input Nifti file to be converted.
    path_gis:   Str, path without extension to the 3 output GIS files.
    verbose:    Bool, print processing infos or not, default False.

    <process>

        <return name="path_gis_dim"  type="File" desc="Path to output Gis .dim file."/>
        <return name="path_gis_ima"  type="File" desc="Path to output Gis .ima file."/>
        <return name="path_gis_minf" type="File" desc="Path to output Gis .minf file."/>

        <input name="path_nifti" type="File" desc="Path to input Nifti file."/>
        <input name="path_gis"   type="Str"  desc="Path without extension to
                                                   the 3 Gis files."/>
        <input name="verbose"    type="Bool" optional="True"/>

    </process>
    """

    cmd = ["PtkNifti2GisConverter", "-i", path_nifti, "-o", path_gis,
           "-verbose", str(verbose), "-verbosePluginLoading", str(False)]
    subprocess.check_call(cmd)

    path_gis_dim  = path_gis + ".dim"
    path_gis_ima  = path_gis + ".ima"
    path_gis_minf = path_gis + ".ima.minf"

    return path_gis_dim, path_gis_ima, path_gis_minf



def gis_to_nifti(path_gis, path_nifti, verbose=False):
    """
    Function that wraps the PtkGis2NiftiConverter command line tool from
    Connectomist to make it capsul-isable.

    Parameters
    ----------
    path_gis:   Str, path without extension to the 3 input GIS files.
    path_nifti: Str, path to the output Nifti file.
    verbose:    Bool, print processing infos or not, default False.

    Returns
    -------
    path_nifti: path to output file.

    <process>

        <return name="path_nifti" type="File" desc="Path to output Nifti file."/>

        <input name="path_gis" type="Str" desc="Path without extension to the
                                                3 input Gis files."/>
        <input name="path_nifti" type="File" desc="Path to output Nifti file."/>
        <input name="verbose"    type="Bool" optional="True"/>

    </process>
    """

    # Add .nii extension if it's not the case
    if not (path_nifti.endswith(".nii") or path_nifti.endswith(".nii.gz")):
        path_nifti += ".nii"

    # Call command line tool:
    # It creates a Nifti + a .minf file (metainformation)
    cmd = ["PtkGis2NiftiConverter", "-i", path_gis, "-o", path_nifti,
           "-verbose", str(verbose), "-verbosePluginLoading", str(False)]
    subprocess.check_call(cmd)

    return path_nifti



def concatenate_volumes(path_inputs, path_output, axis="t", verbose=False):
    """ Function that wraps the PtkCat command line tool from Connectomist
        to make it capsul-isable, with only the basic arguments. It allows
        concatenating volumes. In particulier concatenante T2 and DW volumes
        in one file at the end of preprocessing.

        Parameters
        ----------
        path_input:  List of str, paths to input volumes.
        path_output: Str, path to output volume.
        axis:        Str, axis along which the concatenation is done.
        verbose:     Bool (optional), default False.

        <process>

            <return name="path_output" type="File"/>

            <input name="path_inputs" type="List"/>
            <input name="path_output" type="File"/>
            <input name="axis"        type="Str"  optional="True"/>
            <input name="verbose"     type="Bool" optional="True"/>

        </process>
    """

    cmd = ["PtkCat", "-i"] + path_inputs + ["-o", path_output, "-t", axis,
           "-verbose", str(verbose), "-verbosePluginLoading", str(False)]
    subprocess.check_call(cmd)

    return path_output


###############################################################################
# PREPROCESSING INITIALIZATION: GATHER INPUT FILES IN GIS FORMAT
###############################################################################


def gather_files_for_preprocessing(path_nifti, path_bval, path_bvec,
                                   working_directory=None, path_b0_maps=None):
    """
    Create a directory, "01-InputData", in working_directory with all files
    in Gis format that are needed to start the preprocessing.
    If working_directory is not given, a temporary directory is automatically
    created.

    Parameters
    ----------
    path_nifti:        Str,
    path_bval:         Str,
    path_bvec:         Str,
    working_directory; Str (optional)
    path_b0_maps:      List of paths (optional)

    Returns
    -------
    working_directory, path_gis, path_bval, path_bvec, path_b0_maps: New paths.
    """

    if working_directory: # If a path is given
        if not os.path.isdir(working_directory):
            # Directory does not exist, create it;
            os.mkdir(working_directory)
    else: # No path given, create a random directory
        working_directory = tempfile.mkdtemp(prefix="connectomist_")

    # Directory where DWI data is gathered
    input_dir = os.path.join(working_directory, "01-InputData")
    if not os.path.isdir(input_dir): # if it does not exist, create it
        os.mkdir(input_dir)

    # Convert nifti to gis
    filename_wo_ext = os.path.basename(path_nifti).rsplit(".gz", 1)[0].rsplit(".nii",1)[0]
    path_gis        = os.path.join(input_dir, filename_wo_ext)
    nifti_to_gis(path_nifti, path_gis)

    # Copy bval and bvec files, with homogeneous names
    new_path_bval = path_gis + ".bval"
    shutil.copy(path_bval, new_path_bval)
    new_path_bvec = path_gis + ".bvec"
    shutil.copy(path_bvec, new_path_bvec)

    new_path_b0_maps = []
    if path_b0_maps:
        for path in path_b0_maps:
            # Convert Nifti to Gis
            filename_wo_ext = os.path.basename(path).rsplit(".gz",1)[0].rsplit(".nii",1)[0]
            path_gis        = os.path.join(input_dir, filename_wo_ext)
            nifti_to_gis(path, path_gis)
            new_path_b0_maps.append(path_gis)

    return working_directory, path_gis, new_path_bval, new_path_bvec, new_path_b0_maps



###############################################################################
# FUNCTIONS THAT WRAPS CONNECTOMIST'S TABS
###############################################################################
# DWI-Data-Import-And-QSpace-Sampling

def dwi_data_import_and_qspace_sampling(path_gis, path_bval, path_bvec,
                                        output_directory, manufacturer,
                                        invert_axes=(False, False, False)):
    """
    Wrapper to Connectomist's "DWI & Q-space" tab.

    Parameters
    ----------
    path_gis:         Str, path to .ima (Gis format) diffusion-weighted data.
    path_bvec:        Str, path to .bval file.
    path_bval:        Str, path to .bvec file.
    manufacturer:     Str, name of the manufacturer (e.g. "Siemens", "GE",
                      "Philips" or "Bruker").
    output_directory: Str, path to Connectomist output work directory.
    invert_axes:      Tuple (optional), 3 boolean indicating whether x,y,z axis
                      from diffusion model should be inverted or not.

    <process>
        <return name="output_directory" type="Directory"/>

        <input name="path_gis"  type="File" desc="Path to .ima (Gis format)
                                                  diffusion-weighted data."/>
        <input name="path_bval" type="File" desc="Path to .bval file."/>
        <input name="path_bvec" type="File" desc="Path to .bvec file."/>
        <input name="output_directory" type="Directory"/>
        <input name="invert_axes"      type="List"
               desc="List of 3 boolean indicating whether x,y,z axis from
                     diffusion model should be inverted or not."/>

    </process>
    """

    algorithm_name = "DWI-Data-Import-And-QSpace-Sampling"

    # Dict with all parameters for connectomist
    parameter_values = {
        # Parameters are ordered as they appear in connectomist's GUI

        # ---------------------------------------------------------------------
        # Field: "Diffusion weighted-images"
        "fileNameDwi":  path_gis, # "DW data"
        "sliceAxis":    2,        # "Slice axis", default "Z-axis"
        "phaseAxis":    1,        # "Phase axis", default "Y-axis"
        "manufacturer": None,

        # Subfield: "Advanced parameters"
        "flipAlongX": 0,          # "Flip data along x"
        "flipAlongY": 0,
        "flipAlongZ": 0,
        "numberOfDiscarded":   0,    # "#discarded images at beginning"
        "numberOfT2":          None, # "#T2"
        "numberOfRepetitions": 1,    # "#repetitions"
        # ---------------------------------------------------------------------
        # Field: "Rotation of field of view", default is identity matrix
        "qSpaceTransform_xx": 1.0,
        "qSpaceTransform_xy": 0.0,
        "qSpaceTransform_xz": 0.0,
        "qSpaceTransform_yx": 0.0,
        "qSpaceTransform_yy": 1.0,
        "qSpaceTransform_yz": 0.0,
        "qSpaceTransform_zx": 0.0,
        "qSpaceTransform_zy": 0.0,
        "qSpaceTransform_zz": 1.0,
        # ---------------------------------------------------------------------
        # Field: "Q-space sampling"
        "qSpaceSamplingType":     4, # default "spherical single-shell custom"
        "qSpaceChoice5BValue": 1300,
        "qSpaceChoice5OrientationFileNames": path_bvec,

        # Apparently Connectomist uses 2 as True, and 0 as False.
        "invertXAxis": 2 if invert_axes[0] else 0,
        "invertYAxis": 2 if invert_axes[1] else 0,
        "invertZAxis": 2 if invert_axes[2] else 0,

        # In this field but not used/handled parameters
        "qSpaceChoice1MaximumBValue": 1000, # case Cartesian
        "qSpaceChoice2BValue":        1000, # case spherical single-shell PTK
        "qSpaceChoice3BValue":        1000, # case spherical single-shell SMS
        "qSpaceChoice4BValue":        1000, # case spherical single-shell GEHC
        "qSpaceChoice6BValues":         "",
        "qSpaceChoice7BValues":         "",
        "qSpaceChoice8BValues":         "",
        "qSpaceChoice9BValues":         "",
        "qSpaceChoice10BValues":        "",
        "qSpaceChoice11BValues":        "",
        "qSpaceChoice12BValues":        "",
        "qSpaceChoice13BValues":        "",
        "qSpaceChoice1NumberOfSteps":         11,
        "qSpaceChoice2NumberOfOrientations":   6,
        "qSpaceChoice3NumberOfOrientations":   6,
        "qSpaceChoice4NumberOfOrientations":   6,
        "qSpaceChoice6NumberOfOrientations":   6,
        "qSpaceChoice7NumberOfOrientations":   6,
        "qSpaceChoice8NumberOfOrientations":   6,
        "qSpaceChoice9OrientationFileNames":  "",
        "qSpaceChoice10NumberOfOrientations": "",
        "qSpaceChoice11NumberOfOrientations": "",
        "qSpaceChoice12NumberOfOrientations": "",
        "qSpaceChoice13OrientationFileNames": "",
        # ---------------------------------------------------------------------
        # Field: "Diffusion time (in ms)"
        "diffusionTime": 1.0,
        # ---------------------------------------------------------------------
        # Field: "Work directory"
        "outputWorkDirectory": output_directory,
        # ---------------------------------------------------------------------
        # unknown parameter
        "_subjectName": "",
    }

    if manufacturer not in MANUFACTURERS:
        raise ValueError("Bad argument: manufacturer should in {}"
                         .format(set(MANUFACTURERS)))
    parameter_values["manufacturer"] = MANUFACTURERS[manufacturer]

    # -------------------------------------------------------------------------
    # Read .bval file to infer nb of T2, nb of shells...

    try:
        bvalues = np.loadtxt(path_bval)
    except Exception as e:
        e.message += "\nFailed to read bval file: %s" % path_bval
        raise

    nb_T2     = np.sum(bvalues == 0)  # nb of volumes where bvalue=0
    bvals_set = set(bvalues) - {0}    # set of non-zero bvalues
    nb_shells = len(bvals_set)

    parameter_values["numberOfT2"] = nb_T2

    if nb_shells == 1:
        parameter_values["qSpaceSamplingType"] = 4 # spherical single-shell custom
    else:
        raise NotImplementedError("Multiple shell models not handled yet.")

    # -------------------------------------------------------------------------
    parameter_file = create_parameter_file(algorithm_name, parameter_values,
                                           output_directory)
    run_connectomist(algorithm_name, parameter_file)

    return output_directory


###############################################################################
# DWI-Rough-Mask-Extraction

def dwi_rough_mask_extraction(raw_dwi_directory, output_directory):
    """
    Wrapper to Connectomist's "Rough mask" tab.

    Parameters
    ----------
    raw_dwi_directory: Str, path to Connectomist Raw DWI folder.
    output_directory:  Str, path to Connectomist output work directory.

    <process>
        <return name="output_directory" type="Directory"/>

        <input name="raw_dwi_directory" type="Directory"/>
        <input name="output_directory"  type="Directory"/>
    </process>
    """
    algorithm_name = "DWI-Rough-Mask-Extraction"

    parameter_values = {
    # -------------------------------------------------------------------------
    # Used parameters
     'outputWorkDirectory': output_directory,
     'rawDwiDirectory':    raw_dwi_directory,
    # -------------------------------------------------------------------------
    # Parameters not used/handled by the code
     '_subjectName': '',
     'anatomy': '',
     'dwToT1RegistrationParameter': {
             'applySmoothing':                   1,
             'floatingLowerThreshold':         0.0,
             'initialParametersRotationX':       0,
             'initialParametersRotationY':       0,
             'initialParametersRotationZ':       0,
             'initialParametersScalingX':      1.0,
             'initialParametersScalingY':      1.0,
             'initialParametersScalingZ':      1.0,
             'initialParametersShearingXY':    0.0,
             'initialParametersShearingXZ':    0.0,
             'initialParametersShearingYZ':    0.0,
             'initialParametersTranslationX':    0,
             'initialParametersTranslationY':    0,
             'initialParametersTranslationZ':    0,
             'initializeCoefficientsUsingCenterOfGravity': True,
             'levelCount':                      32,
             'maximumIterationCount':         1000,
             'maximumTestGradient':         1000.0,
             'maximumTolerance':              0.01,
             'optimizerName':                    0,
             'optimizerParametersRotationX':     5,
             'optimizerParametersRotationY':     5,
             'optimizerParametersRotationZ':     5,
             'optimizerParametersScalingX':   0.05,
             'optimizerParametersScalingY':   0.05,
             'optimizerParametersScalingZ':   0.05,
             'optimizerParametersShearingXY': 0.05,
             'optimizerParametersShearingXZ': 0.05,
             'optimizerParametersShearingYZ': 0.05,
             'optimizerParametersTranslationX': 30,
             'optimizerParametersTranslationY': 30,
             'optimizerParametersTranslationZ': 30,
             'referenceLowerThreshold':        0.0,
             'resamplingOrder':                  1,
             'similarityMeasureName':            1,
             'stepSize':                       0.1,
             'stoppingCriterionError':        0.01,
             'subSamplingMaximumSizes':       '64',
             'transform3DType':                  0},
     'maskClosingRadius':        0.0,
     'maskDilationRadius':       4.0,
     'morphologistBrainMask':     '',
     'noiseThresholdPercentage': 2.0,
     'strategyRoughMaskFromT1':    0,
     'strategyRoughMaskFromT2':    1}


    parameter_file = create_parameter_file(algorithm_name, parameter_values,
                                           output_directory)
    run_connectomist(algorithm_name, parameter_file)

    return output_directory


###############################################################################
# DWI-Outlier-Detection

def dwi_outlier_detection(raw_dwi_directory, rough_mask_directory,
                          output_directory):
    """
    Wrapper to Connectomist's "Outliers" tab.

    Parameters
    ----------
    raw_dwi_directory:    Str, path to Connectomist Raw DWI folder.
    rough_mask_directory: Str, path to Connectomist Rough Mask folder.
    output_directory:     Str, path to Connectomist output work directory.

    <process>
        <return name="output_directory" type="Directory"/>

        <input name="raw_dwi_directory"    type="Directory"/>
        <input name="rough_mask_directory" type="Directory"/>
    </process>
    """

    algorithm_name = "DWI-Outlier-Detection"

    parameter_values = {'rawDwiDirectory':     raw_dwi_directory,
                        'roughMaskDirectory':  rough_mask_directory,
                        'outputWorkDirectory': output_directory,
                        '_subjectName': '',
                        'discardedOrientationList': '',
                        'outlierFactor': 3.0}

    parameter_file = create_parameter_file(algorithm_name, parameter_values,
                                           output_directory)
    run_connectomist(algorithm_name, parameter_file)

    return output_directory


###############################################################################
# DWI-Susceptibility-Artifact-Correction

def dwi_susceptibility_artifact_correction():
    """
    Wrapper to Connectomist's "Susceptibility" tab.
    """
    algorithm_name = "DWI-Susceptibility-Artifact-Correction"
    raise NotImplementedError(algorithm_name)


###############################################################################
# DWI-Eddy-Current-And-Motion-Correction

# TO BE COMPLETED: capsul
def dwi_eddy_current_and_motion_correction(raw_dwi_directory,
                                           rough_mask_directory,
                                           outliers_directory,
                                           output_directory):
    """
    Wrapper to Connectomist's "Eddy current & motion" tab.

    Parameters
    ----------
    raw_dwi_directory:    Str, path to Connectomist Raw DWI directory.
    rough_mask_directory: Str, path to Connectomist Rough Mask directory.
    outliers_directory:   Str, path to Connectomist Outliers directory.
    output_directory:     Str, path to Connectomist output work directory.

    <process>
        <return name="output_directory"/>
    </process>
    """
    algorithm_name = "DWI-Eddy-Current-And-Motion-Correction"

    parameter_values = {
            # -----------------------------------------------------------------
            # Used parameters
            "rawDwiDirectory":       raw_dwi_directory,
            "roughMaskDirectory":    rough_mask_directory,
            "correctedDwiDirectory": outliers_directory,
            "outputWorkDirectory":   output_directory,
            "eddyCurrentCorrection": 2,
            "motionCorrection":      1,
            # -----------------------------------------------------------------
            # Parameters not used/handled by the code
            "_subjectName": "",
            "fileNameMotionTransform": "",
            "eddyCurrentCorrectionOptions": {
                "optimizerParametersTranslationY":  2,
                "optimizerParametersTranslationX":  2,
                "optimizerParametersTranslationZ":  2,
                "maximumTestGradient":         1000.0,
                "subSamplingMaximumSizes":       "64",
                "optimizerName":                    0,
                "optimizerParametersShearingYZ": 0.01,
                "initialParametersTranslationZ":    0,
                "optimizerParametersShearingXZ": 0.01,
                "initialParametersRotationX":       0,
                "initialParametersRotationY":       0,
                "initialParametersRotationZ":       0,
                "maximumIterationCount":         1000,
                "registrationResamplingOrder":      1,
                "optimizerParametersShearingXY": 0.01,
                "optimizerParametersRotationX":     2,
                "optimizerParametersRotationZ":     2,
                "optimizerParametersScalingZ":   0.01,
                "optimizerParametersScalingY":   0.01,
                "backgroundResamplingLevel":        0,
                "initialParametersShearingXZ":    0.0,
                "initialParametersShearingXY":    0.0,
                "outputResamplingOrder":            3,
                "optimizerParametersScalingX":   0.01,
                "lowerThreshold":                 0.0,
                "stoppingCriterionError":        0.01,
                "maximumTolerance":              0.01,
                "levelCount":                      32,
                "initialParametersTranslationX":    0,
                "initialParametersTranslationY":    0,
                "optimizerParametersRotationY":     2,
                "stepSize":                       0.1,
                "initialParametersShearingYZ":    0.0,
                "similarityMeasureName":            1,
                "applySmoothing":                   1,
                "initialParametersScalingX":      1.0,
                "initialParametersScalingY":      1.0,
                "initialParametersScalingZ":      1.0
            },
            "motionCorrectionOptions": {
                "subSamplingMaximumSizes":       "64",
                "optimizerParametersTranslationX":  2,
                "stoppingCriterionError":        0.01,
                "stepSize":                       0.1,
                "optimizerParametersTranslationY":  2,
                "optimizerName":                    0,
                "optimizerParametersTranslationZ":  2,
                "maximumTestGradient":         1000.0,
                "initialParametersRotationX":       0,
                "initialParametersRotationY":       0,
                "initialParametersRotationZ":       0,
                "maximumIterationCount":         1000,
                "applySmoothing":                   1,
                "optimizerParametersRotationX":     2,
                "optimizerParametersRotationZ":     2,
                "backgroundResamplingLevel":        0,
                "outputResamplingOrder":            3,
                "lowerThreshold":                 0.0,
                "maximumTolerance":              0.01,
                "initialParametersTranslationZ":    0,
                "levelCount":                      32,
                "initialParametersTranslationX":    0,
                "initialParametersTranslationY":    0,
                "registrationResamplingOrder":      1,
                "similarityMeasureName":            1,
                "optimizerParametersRotationY":     2
            }
        }

    parameter_file = create_parameter_file(algorithm_name, parameter_values,
                                           output_directory)
    run_connectomist(algorithm_name, parameter_file)

    return output_directory

###############################################################################


def export_corrected_nifti_with_metafiles(eddy_motion_directory, path_nifti,
                                          path_bval=None, path_bvec=None):
    """
    After Eddy current and motion correction, export the result as Nifti with
    metafiles (.bval and .bvec with corrrected directions).

    Parameters
    ----------
    eddy_motion_directory: Str, path to Connectomist "Eddy Current & Motion"
                           directory.
    path_nifti: Str, path to output Nifti file.
    path_bval:  Str (optional), path to output .bval file.
                By default same as path_nifti but with .bval extension.
    path_bvec:  Same a path_bval for .bvec file.
    
    Returns
    -------
    path_nifti, path_bval, path_bvec: (Str, Str, Str) paths to output files.
    
    <process>
    
        <return name="path_nifti" type="File" />
        <return name="path_bval"  type="File" />
        <return name="path_bvec"  type="File" />
    
        <input name="eddy_motion_directory" type="Directory" 
               desc="Path to Connectomist "Eddy Current & Motion" directory."/>
        <input name="path_nifti" type="File" desc="Path to output Nifti file."/>
        <input name="path_bval" optional="True" type="File"
                   desc="Path to the output .bval file. By default it is the
                         same as path_nifti with .bval extension."/>
        <input name="path_bvec" optional="True" type="File"
               desc="path to the output .bval file. By default it is the
                     same as path_nifti with .bvec extension."/>

    </process>
    """

    # Concatenate preprocessed T2 and preprocessed DW volumes
    path_t2 = os.path.join(eddy_motion_directory, "t2_wo_eddy_current_and_motion")
    path_dw = os.path.join(eddy_motion_directory, "dw_wo_eddy_current_and_motion")
    path_t2_dw_gis = os.path.join(eddy_motion_directory, "t2_dw_wo_eddy_current_and_motion")
    concatenate_volumes([path_t2, path_dw], path_t2_dw_gis)

    # Convert to Nifti
    gis_to_nifti(path_t2_dw_gis, path_nifti)
    
    # Load a dict from the metadata file (<nifti path>.minf).
    path_minf = path_dw + ".ima.minf"
    variables = dict() # To store variables created by execfile call.
    execfile(path_minf, variables)

    # If the dict is correct, there should be an "attributes" key.
    if not "attributes" in variables:
        raise Exception("Failed to read metadata from %s" % path_minf)

    # Get bvalues and create .bval file
    bvalues = np.array(variables["attributes"]["bvalues"], dtype=np.int)
    bvalues = np.concatenate(([0], bvalues)) # don't forget b=0 associated to T2 
    if not path_bval:
        path_bval = path_nifti.rsplit(".gz", 1)[0].rsplit(".nii", 1)[0] + ".bval"
    np.savetxt(path_bval, bvalues, newline=" ", fmt="%d")

    # Get gradient directions and create .bvec file
    directions = np.array(variables["attributes"]["diffusion_gradient_orientations"])
    if not path_bvec:
        path_bvec = path_nifti.rsplit(".gz", 1)[0].rsplit(".nii", 1)[0] + ".bvec"
    np.savetxt(path_bvec, directions.T, fmt="%.10f")
    
    return path_nifti, path_bval, path_bvec


###############################################################################

# TO BE COMPLETED Susceptibility + DOC + capsul
def complete_preprocessing(path_nifti, path_bval, path_bvec, manufacturer,
                           working_directory=None, path_b0_maps=None,
                           invert_axes=(False, False, False)):
    """
    Function that runs all preprocessing tabs from Connectomist.

    Parameters
    ----------
    path_nifti         Str, path to input Nifti diffusion volumes.
    path_bval:         Str, path to Nifti's associated .bval file.
    path_bvec:         Str, path to Nifti's associated .bval file.
    manufacturer:      Str, name of manufacturer associated to DWI acquisition.
    working_directory: Str (optional), path to directory where the preprocessing
                       will be done. If not path given, a random temp directory
                       is created.
    path_b0_maps:      List of path(s) (optional), 1 or 2 Nifti path(s) to B0 maps.
    invert_axes:       Tuple (optional), 3 boolean indicating whether x,y,z axis
                       from diffusion model should be inverted or not.


    Returns
    -------
    path_nifti, path_bval, path_bvec: (Str, Str, Str) paths to output files.

    <process>
    </process>
    """

    # Step 1 - Bring all files in Gis format in a directory
    working_directory, path_gis, path_bval, path_bvec, path_b0_maps = \
        gather_files_for_preprocessing(path_nifti, path_bval, path_bvec,
                                       working_directory, path_b0_maps)

    # Step 2 - Import files to Connectomist and choose diffusion model
    raw_dwi_directory = os.path.join(working_directory, "02-Raw_Dwi")
    dwi_data_import_and_qspace_sampling(path_gis, path_bval, path_bvec,
                                        raw_dwi_directory, manufacturer,
                                        invert_axes)

    # Step 3 - Create a brain mask
    rough_mask_directory = os.path.join(working_directory, "03-Rough_Mask")
    dwi_rough_mask_extraction(raw_dwi_directory, rough_mask_directory)

    # Step 4 - Detect and filter outlying diffusion slices
    outlier_directory = os.path.join(working_directory, "04-Outlier")
    dwi_outlier_detection(raw_dwi_directory, rough_mask_directory, outlier_directory)

    # Step 5 - Susceptibility correction
#    susceptibility_directory = os.path.join(working_directory, "05-Susceptiblity")
#    dwi_susceptibility_artifact_correction()

    # Step 6 - Eddy current and motion correction
    eddy_motion_directory = os.path.join(working_directory, "06-Eddy_Current_And_Motion")
    dwi_eddy_current_and_motion_correction(raw_dwi_directory, rough_mask_directory,
                                           outlier_directory, eddy_motion_directory)

    # Step 7 - Export result as a Nifti with a .bval and a .bvec
    nifti_filename = os.path.basename(path_gis).rsplit(".ima",1)[0] + "_preprocessed.nii"
    path_nifti = os.path.join(working_directory, nifti_filename)
    path_nifti, path_bval, path_bvec = \
        export_corrected_nifti_with_metafiles(eddy_motion_directory, path_nifti)
    
    return path_nifti, path_bval, path_bvec



###############################################################################

# TO BE COMPLETED: DOC + capsul + test
def parallel_preprocessing(nb_processes, path_niftis, path_bvals, path_bvecs,
                           manufacturers, working_directories,
                           path_b0_maps, tuples_invert_axes):
    """
    Function to make complete_preprocessing() run in parallel. The arguments
    are the same as complete_preprocessing() but with lists for the data args.
    """
    p = Pool(nb_processes)

    # Create a list of tuples.
    # Each tuple is one "set" of arguments for a complete_preprocessing() call.
    args = zip(path_niftis, path_bvals, path_bvecs, manufacturers,
               working_directories, tuples_invert_axes)

    # The map() method from multiprocessing.Pool() works only with single-
    # argument functions. To allow parallelization of multiple-argument
    # functions we pass all the arguments as a tuple to an intermediate
    # function, i.e. that behaves like single argument function, that unpacks
    # the arguments and properly calls the function to be run in //.
    def intermediate_hack(args):
        return complete_preprocessing(*args)
    p.map(intermediate_hack, args)


