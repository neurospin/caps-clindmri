#!/usr/bin/env python
# -*- coding: utf-8 -*-

##########################################################################
# NSAp - Copyright (C) CEA, 2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

import os
import shutil
import numpy as np
import nibabel

from .manufacturers import MANUFACTURERS
from .exceptions    import (ConnectomistError, BadManufacturerNameError,
                            BadFileError)
from .utils         import create_parameter_file, run_connectomist, nifti_to_gis


def gather_and_format_input_files(output_directory,
                                  path_dwi,
                                  path_bval,
                                  path_bvec,
                                  path_b0_magnitude,
                                  path_b0_phase = None):
    """
    Gather all files needed to start the preprocessing in the right format
    (Gis format for images and B0 maps).

    Parameters
    ----------
    output_directory:  Str, path to directory where to gather all files.
    path_dwi:          Str, path to input Nifti file to be preprocessed.
    path_bval:         Str, path to .bval file associated to Nifti.
    path_bvec:         Str, path to .bvec file associated to Nifti.
    path_b0_magnitude: Str, path to B0 magnitude map, may also contain phase.
    path_b0_phase:     Str, not required if phase is already contained in
                       path_b0_magnitude.

    Returns
    -------
    raw_dwi_directory: Connectomist's import data directory

    <unit>
        <output name="raw_dwi_directory" type="Directory"     />

        <input name="output_directory"       type="Directory" />
        <input name="path_dwi"               type="File"      />
        <input name="path_bval"              type="File"      />
        <input name="path_bvec"              type="File"      />
        <input name="path_b0_magnitude"      type="File"      />
        <input name="path_b0_phase"          type="File"      />
    </unit>
    """

    # Check that files exist
    files = [path_dwi, path_bval, path_bvec, path_b0_magnitude]
    if path_b0_phase:
        files.append(path_b0_phase)
    for path in files:
        if not os.path.isfile(path):
            raise BadFileError(path)

    # Create the directory if not existing
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # If there is only one b0 map file and if this file contains 2 volumes,
    # split it in 2 files: magnitude and phase, assuming the first one is magnitude
    if path_b0_magnitude and not path_b0_phase:
        b0_maps = nibabel.load(path_b0_magnitude)
        if b0_maps.shape[-1] == 2:
            b0_maps   = nibabel.load(path_b0_magnitude)
            voxels    = b0_maps.get_data()
            header    = b0_maps.get_header()
            affine    = b0_maps.get_affine()
            magnitude = nibabel.Nifti1Image(voxels[:,:,:,0], affine, header)
            phase     = nibabel.Nifti1Image(voxels[:,:,:,1], affine, header)
            path_b0_magnitude = os.path.join(output_directory, "b0_magnitude.nii.gz")
            path_b0_phase     = os.path.join(output_directory, "b0_phase.nii.gz")
            magnitude.to_filename(path_b0_magnitude)
            phase    .to_filename(path_b0_phase)

    # Convert Nifti to Gis
    path_dwi = nifti_to_gis(path_dwi, os.path.join(output_directory, "dwi.ima"))

    # Copy bval and bvec files, with homogeneous names
    path_bval_copy = os.path.join(output_directory, "dwi.bval")
    path_bvec_copy = os.path.join(output_directory, "dwi.bvec")
    shutil.copyfile(path_bval, path_bval_copy)
    shutil.copyfile(path_bvec, path_bvec_copy)

    # Convert and rename B0 map(s)
    path_b0_magnitude = nifti_to_gis(path_b0_magnitude,
                                     os.path.join(output_directory, "b0_magnitude.ima"))

    if path_b0_phase:
        path_b0_phase = nifti_to_gis(path_b0_phase,
                                     os.path.join(output_directory, "b0_phase.ima"))
    else:
        path_b0_phase = None

    return (path_dwi, path_bval_copy, path_bvec_copy, path_b0_magnitude,
            path_b0_phase)


def dwi_data_import_and_qspace_sampling(output_directory,
                                        path_dwi,
                                        path_bval,
                                        path_bvec,
                                        manufacturer,
                                        path_b0_magnitude,
                                        path_b0_phase = None,
                                        invertX       = True,
                                        invertY       = False,
                                        invertZ       = False,
                                        subject_id    = None):
    """
    Wrapper to Connectomist's "DWI & Q-space" tab.

    Parameters
    ----------
    output_directory: Str, path to Connectomist's output directory.
    path_dwi:         Str, path to Nifti diffusion-weighted data.
    path_bvec:        Str, path to .bval file associated to the Nifti.
    path_bval:        Str, path to .bvec file associated to the Nifti.
    manufacturer:     Str, name of the manufacturer (e.g. "Siemens", "GE",
                      "Philips" or "Bruker").
    invertX:          Bool, if True invert x-axis of diffusion model.
    invertY           Same as invertX for y-axis.
    invertZ:          Same as invertX for z-axis.

    <unit>
        <output name="raw_dwi_directory" type="Directory" description="Path to
            Connectomist output directory."/>

        <input name="output_directory" type="Directory" />
        <input name="path_dwi"  type="File"             />
        <input name="path_bval" type="File"             />
        <input name="path_bvec" type="File"             />
        <input name="manufacturer" type="Str" description="Name of the MRI
            manufacturer (e.g. 'Siemens', 'GE', 'Philips' or 'Bruker')." />
        <input name="path_b0_magnitude" type="File"     />
        <input name="path_b0_phase"     type="File"     />
        <input name="invertX" type="Bool" description="If True invert x-axis
            of diffusion model."/>
        <input name="invertY" type="Bool" description="Same as invertX but for
            y-axis."/>
        <input name="invertZ" type="Bool" description="Same as invertX but for
            z-axis."/>
        <input name="subject_id" type="Str" description="Subject's identifier." />
    </unit>
    """

    (path_dwi, path_bval, path_bvec, path_b0_magnitude, path_b0_phase) = \
        gather_and_format_input_files(output_directory,
                                      path_dwi,
                                      path_bval,
                                      path_bvec,
                                      path_b0_magnitude,
                                      path_b0_phase)

    algorithm_name = "DWI-Data-Import-And-QSpace-Sampling"

    # Dict with all parameters for connectomist
    parameters_value = {
        # Parameters are ordered as they appear in connectomist's GUI

        # ---------------------------------------------------------------------
        # Field: "Diffusion weighted-images"
        "fileNameDwi":   path_dwi,  # "DW data"
        "sliceAxis":            2,  # "Slice axis", default "Z-axis"
        "phaseAxis":            1,  # "Phase axis", default "Y-axis"
        "manufacturer":      None,

        # Subfield: "Advanced parameters"
        "flipAlongX":           0,  # "Flip data along x"
        "flipAlongY":           0,
        "flipAlongZ":           0,
        "numberOfDiscarded":    0,  # "#discarded images at beginning"
        "numberOfT2":        None,  # "#T2"
        "numberOfRepetitions":  1,  # "#repetitions"
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
        "qSpaceChoice1MaximumBValue":       1000,  # case Cartesian
        "qSpaceChoice2BValue":              1000,
        "qSpaceChoice3BValue":              1000,
        "qSpaceChoice4BValue":              1000,
        "qSpaceChoice6BValues":               "",
        "qSpaceChoice7BValues":               "",
        "qSpaceChoice8BValues":               "",
        "qSpaceChoice9BValues":               "",
        "qSpaceChoice10BValues":              "",
        "qSpaceChoice11BValues":              "",
        "qSpaceChoice12BValues":              "",
        "qSpaceChoice13BValues":              "",
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
        "_subjectName": subject_id if subject_id else "",
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
        raise ConnectomistError("Multiple shell models not handled. "
                                "Path to .bval file: %s" % path_bval)

    # Check validity of .bvec file.
    # If bvec file does not exist or filled with 0s, raise Error
    if (not os.path.isfile(path_bvec)) or np.loadtxt(path_bvec).max() == 0:
        raise BadFileError(path_bvec)

    parameter_file = create_parameter_file(algorithm_name,
                                           parameters_value,
                                           output_directory)
    run_connectomist(algorithm_name, parameter_file)

    # Capsul needs the output name to be different from input arguments
    raw_dwi_directory = output_directory
    return raw_dwi_directory
