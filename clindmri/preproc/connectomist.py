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
import logging
import nibabel
from multiprocessing import Process, Manager

###############################################################################
# MODULE VARIABLES
###############################################################################

# Note: we use simplified manufacturer names to simplify usage.
# e.g. we use "GE" where Connectomist uses "GE HealthCare",
# but a DICOM header would use "GE Medical Systems"
MANUFACTURERS = {
    "Bruker":  0,
    "GE":      1,
    "Philips": 2,
    "Siemens": 3
}

# Messages for communication between processes in multiprocessing
FLAG_STOP_PROCESS     = "STOP_WORK"
FLAG_PROCESS_FINISHED = "PROCESS_HAS_FINISHED"

###############################################################################
# MODULE EXCEPTION
###############################################################################

class ConnectomistWrapperError(Exception):
    """ New exception type to differentiate errors from modules wrapping
        Connectomist from other exceptions. No specific code.
    """
    def __init__(self, message):
        super(ConnectomistWrapperError, self).__init__(message)


class ConnectomistRuntimeError(ConnectomistWrapperError):
    """ Error thrown when call to the Connectomist software failed.
    """
    def __init__(self, algorithm_name, parameter_file):
        message = "Connectomist call for %s failed, with parameters: %s." \
                  % (algorithm_name, parameter_file)
        super(ConnectomistRuntimeError, self).__init__(message)


class BadManufacturerNameError(ConnectomistWrapperError):
    """ Error thrown when an incorrect manufacturer name is detected.
    """
    def __init__(self, manufacturer):
        message = "Incorrect manufacturer name: %s, should be in {}." \
                  .format(set(MANUFACTURERS))
        super(BadManufacturerNameError, self).__init__(message)


class MissingParametersError(ConnectomistWrapperError):
    """ Error thrown when needed parameters were not all given.
    """
    def __init__(self, algorithm_name, missing_parameters):
        message = "Missing parameters for {}: {}.".format(algorithm_name,
                                                          missing_parameters)
        super(MissingParametersError, self).__init__(message)


class BadFileError(ConnectomistWrapperError):
    """ Error thrown when a file is missing or corrupted.
    """
    def __init__(self, file_path):
        message = "Missing or corrupted file: %s" % file_path
        super(BadFileError, self).__init__(message)


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
        <input name="path_gis"    type="Str"  desc="Path without extension to
                                                    the 3 input Gis files."/>
        <input name="path_nifti"  type="File" desc="Path to output Nifti file."/>
        <input name="verbose"     type="Bool" optional="True"/>
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

# TO BE COMPLETED
def gather_and_format_input_files(path_nifti,
                                  path_bval,
                                  path_bvec,
                                  path_b0_magnitude,
                                  path_b0_phase     = None,
                                  working_directory = None):
    """
    Create a directory, "01-Input_Data", in working_directory with all files
    in Gis format that are needed to start the preprocessing.
    If working_directory is not given, a temporary directory is automatically
    created.

    Parameters
    ----------
    path_nifti:        Str, path to input Nifti file to be preprocessed.
    path_bval:         Str, path to .bval file associated to Nifti.
    path_bvec:         Str, path to .bvec file associated to Nifti.
    path_b0_magnitude: Str, path to B0 magnitude map, also contains phase for GE.
    path_b0_phase:     Str, not for GE, path to B0 phase map.
    working_directory: Str (optional), path to folder where all the preproces-
                       sing will be done. If no path is given, a random tempo-
                       rary directory will be created.

    Returns
    -------
    working_directory, path_gis, path_bval, path_bvec, path_b0_maps: New paths.
    """

    # Check that files exist
    files = [path_nifti, path_bval, path_bvec, path_b0_magnitude]
    if path_b0_phase:
        files.append(path_b0_phase)

    for path in files:
        if not os.path.isfile(path):
            raise BadFileError(path)
            
    if working_directory:  # If a path is given
        if not os.path.isdir(working_directory):
            # Directory does not exist, create it;
            os.mkdir(working_directory)
    else:  # No path given, create a random directory
        working_directory = tempfile.mkdtemp(prefix="connectomist_")
    
    # Directory where DWI data is gathered
    input_directory = os.path.join(working_directory, "01-Input_Data")
    if not os.path.isdir(input_directory):  # if it does not exist, create it
        os.mkdir(input_directory)

    # If there is only one b0 map file and if this file contains 2 volumes,
    # split it in 2 files: magnitude and phase
    if path_b0_magnitude and not path_b0_phase:
        b0_maps = nibabel.load(path_b0_magnitude)
        if b0_maps.shape[-1] == 2:
            b0_maps   = nibabel.load(path_b0_magnitude)
            voxels    = b0_maps.get_data()
            header    = b0_maps.get_header()
            affine    = b0_maps.get_affine()
            magnitude = nibabel.Nifti1Image(voxels[:,:,:,0], affine, header)
            phase     = nibabel.Nifti1Image(voxels[:,:,:,1], affine, header)
            path_b0_magnitude = os.path.join(input_directory, "b0_magnitude.nii.gz")
            path_b0_phase     = os.path.join(input_directory, "b0_phase.nii.gz")
            magnitude.to_filename(path_b0_magnitude)
            phase    .to_filename(path_b0_phase)

    # Convert Nifti to Gis
    filename_wo_ext = os.path.basename(path_nifti).rsplit(".gz", 1)[0].rsplit(".nii",1)[0]
    path_gis        = os.path.join(input_directory, filename_wo_ext)
    nifti_to_gis(path_nifti, path_gis)

    # Copy bval and bvec files, with homogeneous names
    path_bval_copy = path_gis + ".bval"
    path_bvec_copy = path_gis + ".bvec"
    shutil.copy(path_bval, path_bval_copy)
    shutil.copy(path_bvec, path_bvec_copy)

    # Convert and rename B0 map(s)
    path_b0_magnitude_gis = os.path.join(input_directory, "b0_magnitude.gis")
    nifti_to_gis(path_b0_magnitude, path_b0_magnitude_gis)

    if path_b0_phase:
        path_b0_phase_gis = os.path.join(input_directory, "b0_phase.gis")
        nifti_to_gis(path_b0_phase, path_b0_phase_gis)
    else:
        path_b0_phase_gis = None

    return (path_gis, path_bval_copy, path_bvec_copy, path_b0_magnitude_gis,
            path_b0_phase_gis, working_directory)


###############################################################################
# FUNCTIONS THAT WRAPS CONNECTOMIST'S TABS
###############################################################################
# DWI-Data-Import-And-QSpace-Sampling

def dwi_data_import_and_qspace_sampling(path_gis,
                                        path_bval,
                                        path_bvec,
                                        output_directory,
                                        manufacturer,
                                        invertX=True,
                                        invertY=False,
                                        invertZ=False):
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
    invertX:          Bool, if True invert x-axis of diffusion model.
    invertY           Same as invertX for y-axis.
    invertZ:          Same as invertX for z-axis.

    <process>
        <return name="output_directory" type="Directory"/>
        <input name="path_gis"  type="File" desc="Path to .ima (Gis format)
                                                  diffusion-weighted data."/>
        <input name="path_bval" type="File" desc="Path to .bval file."/>
        <input name="path_bvec" type="File" desc="Path to .bvec file."/>
        <input name="output_directory" type="Directory"/>
        <input name="invertX" type="Bool" desc="If True invert x-axis of
                                                diffusion model"/>
        <input name="invertY" type="Bool" desc="Same as invertX but for y-axis"/>
        <input name="invertZ" type="Bool" desc="Same as invertX but for z-axis"/>
    </process>
    """

    algorithm_name = "DWI-Data-Import-And-QSpace-Sampling"

    # Dict with all parameters for connectomist
    parameters_value = {
        # Parameters are ordered as they appear in connectomist's GUI

        # ---------------------------------------------------------------------
        # Field: "Diffusion weighted-images"
        "fileNameDwi":  path_gis,  # "DW data"
        "sliceAxis":    2,         # "Slice axis", default "Z-axis"
        "phaseAxis":    1,         # "Phase axis", default "Y-axis"
        "manufacturer": None,

        # Subfield: "Advanced parameters"
        "flipAlongX": 0,  # "Flip data along x"
        "flipAlongY": 0,
        "flipAlongZ": 0,
        "numberOfDiscarded":   0,     # "#discarded images at beginning"
        "numberOfT2":          None,  # "#T2"
        "numberOfRepetitions": 1,     # "#repetitions"
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
        "qSpaceSamplingType":     4,  # default "spherical single-shell custom"
        "qSpaceChoice5BValue": 1300,
        "qSpaceChoice5OrientationFileNames": path_bvec,

        # Apparently Connectomist uses 2 as True, and 0 as False.
        "invertXAxis": 2 if invertX else 0,
        "invertYAxis": 2 if invertY else 0,
        "invertZAxis": 2 if invertZ else 0,

        # In this field but not used/handled parameters
        "qSpaceChoice1MaximumBValue": 1000,  # case Cartesian
        "qSpaceChoice2BValue":        1000,  # case spherical single-shell PTK
        "qSpaceChoice3BValue":        1000,  # case spherical single-shell SMS
        "qSpaceChoice4BValue":        1000,  # case spherical single-shell GEHC
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
        raise BadManufacturerNameError(manufacturer)

    parameters_value["manufacturer"] = MANUFACTURERS[manufacturer]

    # Read .bval file to infer nb of T2, nb of shells...
    try:
        bvalues = np.loadtxt(path_bval)
        if set(bvalues) == {0}:  # If only 0s raise Exception
            raise Exception
    except:
        raise BadFileError(path_bval)

    nb_T2     = np.sum(bvalues == 0)  # nb of volumes where bvalue=0
    bvals_set = set(bvalues) - {0}    # set of non-zero bvalues
    nb_shells = len(bvals_set)

    parameters_value["numberOfT2"] = nb_T2

    if nb_shells == 1:
        # Spherical single-shell custom
        parameters_value["qSpaceSamplingType"] = 4
    else:
        raise ConnectomistWrapperError("Multiple shell models not handled.")

    # Check validity of .bvec file.
    # If bvec file does not exist or filled with 0s, raise Error
    if (not os.path.isfile(path_bvec)) or np.loadtxt(path_bvec).max() == 0:
        raise BadFileError(path_bvec)

    parameter_file = create_parameter_file(algorithm_name,
                                           parameters_value,
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
        <return name="output_directory"  type="Directory"/>
        <input  name="raw_dwi_directory" type="Directory"/>
        <input  name="output_directory"  type="Directory"/>
    </process>
    """
    algorithm_name = "DWI-Rough-Mask-Extraction"

    parameters_value = {
        # ---------------------------------------------------------------------
        # Used parameters
        'outputWorkDirectory': output_directory,
        'rawDwiDirectory':    raw_dwi_directory,
        # ---------------------------------------------------------------------
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
            'transform3DType':                  0
        },
        'maskClosingRadius':        0.0,
        'maskDilationRadius':       4.0,
        'morphologistBrainMask':     '',
        'noiseThresholdPercentage': 2.0,
        'strategyRoughMaskFromT1':    0,
        'strategyRoughMaskFromT2':    1
    }

    parameter_file = create_parameter_file(algorithm_name,
                                           parameters_value,
                                           output_directory)
    run_connectomist(algorithm_name, parameter_file)

    return output_directory


###############################################################################
# DWI-Outlier-Detection

def dwi_outlier_detection(raw_dwi_directory,
                          rough_mask_directory,
                          output_directory):
    """
    Wrapper to Connectomist's "Outliers" tab.

    Parameters
    ----------
    raw_dwi_directory:    Str, path to Connectomist Raw DWI folder.
    rough_mask_directory: Str, path to Connectomist Rough Mask folder.
    output_directory:     Str, path to Connectomist output work directory.

    <process>
        <return name="output_directory"     type="Directory"/>
        <input  name="raw_dwi_directory"    type="Directory"/>
        <input  name="rough_mask_directory" type="Directory"/>
    </process>
    """

    algorithm_name = "DWI-Outlier-Detection"

    parameters_value = {'rawDwiDirectory':     raw_dwi_directory,
                        'roughMaskDirectory':  rough_mask_directory,
                        'outputWorkDirectory': output_directory,
                        '_subjectName': '',
                        'discardedOrientationList': '',
                        'outlierFactor': 3.0}

    parameter_file = create_parameter_file(algorithm_name,
                                           parameters_value,
                                           output_directory)
    run_connectomist(algorithm_name, parameter_file)

    return output_directory


###############################################################################
# DWI-Susceptibility-Artifact-Correction

# TO BE COMPLETED: capsul
def dwi_susceptibility_artifact_correction(raw_dwi_directory,
                                           rough_mask_directory,
                                           outliers_directory,
                                           output_directory,
                                           manufacturer,
                                           delta_TE,
                                           partial_fourier_factor,
                                           parallel_acceleration_factor,
                                           path_b0_magnitude,
                                           path_b0_phase             = None,
                                           negative_sign             = False,
                                           echo_spacing              = None,
                                           EPI_factor                = None,
                                           b0_field                  = 3.0,
                                           water_fat_shift_per_pixel = 4.68):
    """
    Wrapper to Connectomist's "Susceptibility" tab.

    Parameters
    ----------
    raw_dwi_directory:    Str, path to Connectomist Raw DWI folder.
    rough_mask_directory: Str, path to Connectomist Rough Mask folder.
    outliers_directory:   Str, path to Connectomist Outliers folder.
    output_directory:     Str, path to Connectomist output work directory.
    manufacturer:         Str, manufacturer name (e.g. "Siemens", "GE"...).
    delta_TE:             Float, difference in seconds between the 2 echoes
                          in B0 magnitude map acquisition.
    partial_fourier_factor: Float (]0;1]), percentage of k-space plane acquired.
    parallel_acceleration_factor: Int, nb of parallel acquisition in k-space plane.
    path_b0_magnitude:    Str, path to B0 magnitude map, also contains phase for GE.
    path_b0_phase:        Str, not for GE, path to B0 phase map.
    negative_sign:        Bool, if True invert direction of unwarping in
                          susceptibility-distortion correction.
    echo_spacing:         Float, not for Philips, acquisition time in ms between
                          2 centers of 2 consecutively acquired lines in k-space.
    EPI_factor:           Int, nb of echoes after one excitation (90 degrees),
                          i.e. echo train length.
    b0_field:             Float, Philips only, B0 field intensity, by default 3.0.
    water_fat_shift_per_pixel: Float, Philips only, default 4.68Hz

    <process>
    </process>
    """
    algorithm_name = "DWI-Susceptibility-Artifact-Correction"

    parameters_value = {
        # ---------------------------------------------------------------------
        # Paths parameters
        'rawDwiDirectory':              raw_dwi_directory,
        'roughMaskDirectory':        rough_mask_directory,
        'outlierFilteredDwiDirectory': outliers_directory,
        'outputWorkDirectory':           output_directory,

        # ---------------------------------------------------------------------
        # Fieldmap correction only
        'correctionStrategy':            0,
        'importDwToB0Transformation':    0,
        'generateDwToB0Transformation':  1,
        'fileNameDwToB0Transformation': '',

        # ---------------------------------------------------------------------
        # Bruker parameters
        'brukerDeltaTE':                    2.46,
        'brukerEchoSpacing':                0.75,
        'brukerPhaseNegativeSign':             0,
        'brukerPartialFourierFactor':        1.0,
        'brukerParallelAccelerationFactor':    1,
        'brukerFileNameFirstEchoB0Magnitude': '',
        'brukerFileNameB0PhaseDifference':    '',

        # ---------------------------------------------------------------------
        # GE parameters
        'geDeltaTE':                                       2.46,
        'geEchoSpacing':                                   0.75,
        'gePhaseNegativeSign':                                0,
        'gePartialFourierFactor':                           1.0,
        'geParallelAccelerationFactor':                       1,
        'geFileNameDoubleEchoB0MagnitudePhaseRealImaginary': '',

        # ---------------------------------------------------------------------
        # Philips parameters
        'philipsDeltaTE':                    2.46,
        'philipsEchoSpacing':                0.75,  # Not requested in GUI; ignored
        'philipsPhaseNegativeSign':             0,
        'philipsPartialFourierFactor':        1.0,
        'philipsParallelAccelerationFactor':    1,
        'philipsFileNameFirstEchoB0Magnitude': '',
        'philipsFileNameB0PhaseDifference':    '',
        'philipsEPIFactor':                   128,
        'philipsStaticB0Field':               3.0,
        'philipsWaterFatShiftPerPixel':       0.0,

        # ---------------------------------------------------------------------
        # Siemens parameters
        'siemensDeltaTE':                     2.46,
        'siemensEchoSpacing':                 0.75,
        'siemensPhaseNegativeSign':              2,
        'siemensPartialFourierFactor':         1.0,
        'siemensParallelAccelerationFactor':     1,
        'siemensFileNameDoubleEchoB0Magnitude': "",
        'siemensFileNameB0PhaseDifference':     "",

        # ---------------------------------------------------------------------
        # Parameters not used/handled by the code
        '_subjectName': '',
        'DwToB0RegistrationParameter': {
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
            'optimizerParametersRotationX':    10,
            'optimizerParametersRotationY':    10,
            'optimizerParametersRotationZ':    10,
            'optimizerParametersScalingX':   0.05,
            'optimizerParametersScalingY':   0.05,
            'optimizerParametersScalingZ':   0.05,
            'optimizerParametersShearingXY': 0.05,
            'optimizerParametersShearingXZ': 0.05,
            'optimizerParametersShearingYZ': 0.05,
            'optimizerParametersTranslationX': 10,
            'optimizerParametersTranslationY': 10,
            'optimizerParametersTranslationZ': 10,
            'referenceLowerThreshold':        0.0,
            'resamplingOrder':                  1,
            'similarityMeasureName':            1,
            'stepSize':                       0.1,
            'stoppingCriterionError':        0.01,
            'subSamplingMaximumSizes':       '56',
            'transform3DType':                  0
        },
    }

    if manufacturer not in MANUFACTURERS:
        raise BadManufacturerNameError(manufacturer)

    # Maps required parameters in Connectomist, for each manufacturer, to the
    # arguments of the function.
    args_map = {
        "Bruker": {
            "brukerDeltaTE":                       delta_TE,
            "brukerPartialFourierFactor":          partial_fourier_factor,
            "brukerParallelAccelerationFactor":    parallel_acceleration_factor,
            "brukerFileNameFirstEchoB0Magnitude":  path_b0_magnitude,
            "brukerFileNameB0PhaseDifference":     path_b0_phase,
            "brukerPhaseNegativeSign":             negative_sign,
            'brukerEchoSpacing':                   echo_spacing
        },
        "GE": {
            "geDeltaTE":                           delta_TE,
            "gePartialFourierFactor":              partial_fourier_factor,
            "geParallelAccelerationFactor":        parallel_acceleration_factor,
            "geFileNameDoubleEchoB0MagnitudePhaseRealImaginary": path_b0_magnitude,
            "gePhaseNegativeSign":                 negative_sign,
            "geEchoSpacing":                       echo_spacing

        },
        "Philips": {
            "philipsDeltaTE":                      delta_TE,
            "philipsPartialFourierFactor":         partial_fourier_factor,
            "philipsParallelAccelerationFactor":   parallel_acceleration_factor,
            "philipsFileNameFirstEchoB0Magnitude": path_b0_magnitude,
            "philipsFileNameB0PhaseDifference":    path_b0_phase,
            "philipsPhaseNegativeSign":            negative_sign,
            "philipsEPIFactor":                    EPI_factor,
            "philipsStaticB0Field":                b0_field,
            "philipsWaterFatShiftPerPixel":        water_fat_shift_per_pixel
        },
        "Siemens": {
            "siemensDeltaTE":                       delta_TE,
            "siemensPartialFourierFactor":          partial_fourier_factor,
            "siemensParallelAccelerationFactor":    parallel_acceleration_factor,
            "siemensFileNameDoubleEchoB0Magnitude": path_b0_magnitude,
            "siemensFileNameB0PhaseDifference":     path_b0_phase,
            "siemensPhaseNegativeSign":             negative_sign,
            "siemensEchoSpacing":                   echo_spacing
        }
    }

    required_parameters = set(args_map[manufacturer])

    # Check that all needed parameters have been passed
    missing_parameters = [p for p in required_parameters
                          if args_map[manufacturer][p] == None]

    if len(missing_parameters) > 0:
        raise MissingParametersError(algorithm_name, missing_parameters)

    # Set given parameters
    for p in required_parameters:
        parameters_value[p] = args_map[manufacturer][p]

    parameter_file = create_parameter_file(algorithm_name,
                                           parameters_value,
                                           output_directory)
    run_connectomist(algorithm_name, parameter_file)

    return output_directory


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

    parameters_value = {
        # ---------------------------------------------------------------------
        # Used parameters
        "rawDwiDirectory":       raw_dwi_directory,
        "roughMaskDirectory":    rough_mask_directory,
        "correctedDwiDirectory": outliers_directory,
        "outputWorkDirectory":   output_directory,
        "eddyCurrentCorrection": 2,
        "motionCorrection":      1,
        # ---------------------------------------------------------------------
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

    parameter_file = create_parameter_file(algorithm_name,
                                           parameters_value,
                                           output_directory)
    run_connectomist(algorithm_name, parameter_file)

    return output_directory

###############################################################################


def export_corrected_nifti_with_metafiles(eddy_motion_directory,
                                          path_nifti,
                                          path_bval=None,
                                          path_bvec=None):
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
               desc="Path to the output .bval file. By default it is the same
                     as path_nifti with .bval extension."/>
        <input name="path_bvec" optional="True" type="File"
               desc="path to the output .bval file. By default it is the same
                     as path_nifti with .bvec extension."/>
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
    try:
        path_minf = path_dw + ".ima.minf"
        context_dict = dict()  # To store variables created by execfile call.
        execfile(path_minf, context_dict)

        # Verify dict structure
        if "bvalues" not in context_dict["attributes"]:
            raise Exception

        if "diffusion_gradient_orientations" not in context_dict["attributes"]:
            raise Exception
    except:
        raise BadFileError(path_minf)

    # Get bvalues and create .bval file
    bvalues = np.array(context_dict["attributes"]["bvalues"], dtype=np.int)
    bvalues = np.concatenate(([0], bvalues))  # don't forget b=0 associated to T2
    if not path_bval:
        path_bval = path_nifti.rsplit(".gz", 1)[0].rsplit(".nii", 1)[0] + ".bval"
    np.savetxt(path_bval, bvalues, newline=" ", fmt="%d")

    # Get gradient directions and create .bvec file
    directions = np.array(context_dict["attributes"]["diffusion_gradient_orientations"])
    if not path_bvec:
        path_bvec = path_nifti.rsplit(".gz", 1)[0].rsplit(".nii", 1)[0] + ".bvec"
    np.savetxt(path_bvec, directions.T, fmt="%.10f")

    return path_nifti, path_bval, path_bvec


###############################################################################

# TO BE COMPLETED: capsul
def complete_preprocessing(path_nifti,
                           path_bval,
                           path_bvec,
                           manufacturer,
                           delta_TE,
                           partial_fourier_factor,
                           parallel_acceleration_factor,
                           path_b0_magnitude,
                           path_b0_phase             = None,
                           invertX                   = True,
                           invertY                   = False,
                           invertZ                   = False,
                           negative_sign             = False,
                           echo_spacing              = None,
                           EPI_factor                = None,
                           b0_field                  = 3.0,
                           water_fat_shift_per_pixel = 4.68,
                           working_directory         = None):
    """
    Function that runs all preprocessing tabs from Connectomist.

    Parameters
    ----------
    path_nifti         Str, path to input Nifti DWI data.
    path_bval:         Str, path to Nifti's associated .bval file.
    path_bvec:         Str, path to Nifti's associated .bval file.
    manufacturer:      Str, manufacturer name (e.g. "Siemens", "GE"...).
    delta_TE:          Float, difference in seconds between the 2 echoes in B0
                       magnitude map acquisition.
    partial_fourier_factor: Float (]0;1]), percentage of k-space plane acquired.
    parallel_acceleration_factor: Int, nb of parallel acquisition in k-space plane.
    path_b0_magnitude: Str, path to B0 magnitude map, also contains phase for GE.
    path_b0_phase:     Str, not for GE, path to B0 phase map.
    invertX:           Bool, if True invert x-axis in diffusion model.
    invertY:           Bool, same as invertX but for y-axis.
    invertZ:           Bool, same as invertX but for z-axis.
    negative_sign:     Bool, if True invert direction of unwarping in
                       susceptibility-distortion correction.
    echo_spacing:      Float, not for Philips, acquisition time in ms between
                       2 centers of 2 consecutively acquired lines in k-space.
    EPI_factor:        Int, nb of echoes after one excitation (90 degrees),
                       i.e. echo train length.
    b0_field:          Float, Philips only, B0 field intensity, by default 3.0.
    water_fat_shift_per_pixel: Float, Philips only, default 4.68Hz
    working_directory: Str (optional), path to folder where all the preproces-
                       sing will be done. If no path is given, a random tempo-
                       rary directory will be created.


    Returns
    -------
    path_nifti, path_bval, path_bvec, working_directory: Paths to output dir/files.

    <process>
    </process>
    """

    # Check parameters

    # Step 1 - Bring all files in Gis format in a directory: "01-Input_Data"
    path_gis, path_bval, path_bvec, path_b0_magnitude, path_b0_phase, \
    working_directory = gather_and_format_input_files(path_nifti,
                                                      path_bval,
                                                      path_bvec,
                                                      path_b0_magnitude,
                                                      path_b0_phase,
                                                      working_directory)

    # Step 2 - Import files to Connectomist and choose diffusion model
    raw_dwi_directory  = os.path.join(working_directory, "02-Raw_Dwi")
    dwi_data_import_and_qspace_sampling(path_gis,
                                        path_bval,
                                        path_bvec,
                                        raw_dwi_directory,
                                        manufacturer,
                                        invertX,
                                        invertY,
                                        invertZ)

    # Step 3 - Create a brain mask
    rough_mask_directory = os.path.join(working_directory, "03-Rough_Mask")
    dwi_rough_mask_extraction(raw_dwi_directory, rough_mask_directory)

    # Step 4 - Detect and filter outlying diffusion slices
    outliers_directory = os.path.join(working_directory, "04-Outliers")
    dwi_outlier_detection(raw_dwi_directory,
                          rough_mask_directory,
                          outliers_directory)

    # Step 5 - Susceptibility correction
    susceptibility_directory = os.path.join(working_directory, "05-Susceptiblity")
    dwi_susceptibility_artifact_correction(raw_dwi_directory,
                                           rough_mask_directory,
                                           outliers_directory,
                                           susceptibility_directory,
                                           manufacturer,
                                           delta_TE,
                                           partial_fourier_factor,
                                           parallel_acceleration_factor,
                                           path_b0_magnitude,
                                           path_b0_phase,
                                           negative_sign,
                                           echo_spacing,
                                           EPI_factor,
                                           b0_field,
                                           water_fat_shift_per_pixel)

    # Step 6 - Eddy current and motion correction
    eddy_motion_directory = os.path.join(working_directory,
                                         "06-Eddy_Current_And_Motion")
    dwi_eddy_current_and_motion_correction(raw_dwi_directory,
                                           rough_mask_directory,
                                           outliers_directory,
                                           eddy_motion_directory)

    # Step 7 - Export result as a Nifti with a .bval and a .bvec
    preprocessed_directory = os.path.join(working_directory, "07-Preprocessed")
    if not os.path.isdir(preprocessed_directory):
        os.mkdir(preprocessed_directory)
    nifti_name = os.path.basename(path_gis).rsplit(".ima",1)[0] + ".nii"
    path_nifti = os.path.join(preprocessed_directory, nifti_name)
    path_nifti, path_bval, path_bvec = \
        export_corrected_nifti_with_metafiles(eddy_motion_directory, path_nifti)

    return path_nifti, path_bval, path_bvec



###############################################################################
# To make preprocessing work in parallel

def parallel_worker(work_queue, result_queue):
    """ Function to make complete_preprocessing work in parallel processing.
    """
    while True:
        new_work = work_queue.get()
        if new_work == FLAG_STOP_PROCESS:
            result_queue.put(FLAG_PROCESS_FINISHED)
            break
        kwargs = new_work
        try:
            path_nifti, path_bval, path_bvec = complete_preprocessing(**kwargs)
            result_queue.put((path_nifti, path_bval, path_bvec))
        except ConnectomistWrapperError as e:
            result_queue.put(e.message)
        except Exception as e:
            e.message += "\nUNKNOWN ERROR HAPPENED"
            result_queue.put(e.message)


# TO BE COMPLETED: capsul + test
def parallel_preprocessing(nb_processes, list_kwargs, log_path=None):
    """
    Function to make complete_preprocessing() run in parallel.

    Parameters
    ----------
    list_kwargs: List of dicts, each dict stores the arguments for one call to
            complete_preprocessing(), i.e. one job.
    """

    ###########################################################################
    # SETTING UP LOGGING SYSTEM
    ###########################################################################
    logger = logging.getLogger(__file__)
    logger.handlers = []
    logger.setLevel(logging.DEBUG)

    if len(logger.handlers) == 0:
        if not log_path:
            log_path = tempfile.mkstemp(prefix="preprocessing_", suffix=".log")[1]
        file_handler = logging.FileHandler(log_path, mode="w")
        file_handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        file_handler   .setFormatter(formatter)
        logger.addHandler(file_handler)

        logger.info("Path to log file: %s" % log_path)

    ###########################################################################

    # Data structures for parallel processing
    manager = Manager()  # multiprocessing.Manager()
    work_queue, result_queue = manager.Queue(), manager.Queue()

    # Add jobs in work_queue
    for kwargs in list_kwargs:
        work_queue.put(kwargs)
    
    # Add poison pills to stop the remote workers
    # When a process gets this job, it will stop
    for n in range(nb_processes):
        work_queue.put(FLAG_STOP_PROCESS)

    # Define processes
    workers = []
    for i in range(nb_processes):
        worker = Process(target=parallel_worker, args=(work_queue, result_queue))
        worker.daemon = True
        workers.append(worker)
        worker.start()

    # Process results and log everything
    nb_finished_processes = 0
    try:
        while True:
            new_result = result_queue.get()
            print new_result
            if new_result == FLAG_PROCESS_FINISHED:
                nb_finished_processes += 1
                if nb_finished_processes == nb_processes:
                    break
            elif type(new_result) is str:
                logger.warning(new_result)
            else:
                logger.info("Successfull preprocessing, resulting files:"
                            "\n%s\n%s\n%s" % new_result)
    except KeyboardInterrupt:  # To stop if user uses ctrl+c
        for worker in workers:
            worker.terminate()
            worker.join()


