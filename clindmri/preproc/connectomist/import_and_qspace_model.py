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
                                  dwi,
                                  bval,
                                  bvec,
                                  b0_magnitude,
                                  b0_phase = None):
    """
    Gather all files needed to start the preprocessing in the right format
    (Gis format for images and B0 maps).

    Parameters
    ----------
    output_directory: Str, path to directory where to gather all files.
    dwi:              Str, path to input Nifti file to be preprocessed.
    bval:             Str, path to .bval file associated to Nifti.
    bvec:             Str, path to .bvec file associated to Nifti.
    b0_magnitude:     Str, path to B0 magnitude map, may also contain phase.
    b0_phase:         Str, not required if phase is already contained in
                      b0_magnitude.

    Returns
    -------
    raw_dwi_directory: Connectomist's import data directory

    <unit>
        <output name="raw_dwi_directory" type="Directory"     />

        <input name="output_directory"  type="Directory" />
        <input name="dwi"               type="File"      />
        <input name="bval"              type="File"      />
        <input name="bvec"              type="File"      />
        <input name="b0_magnitude"      type="File"      />
        <input name="b0_phase"          type="File"      />
    </unit>
    """

    # Check that files exist
    files = [dwi, bval, bvec, b0_magnitude]
    if b0_phase:
        files.append(b0_phase)
    for path in files:
        if not os.path.isfile(path):
            raise BadFileError(path)

    # Create the directory if not existing
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # If there is only one b0 map file and if this file contains 2 volumes,
    # split it in 2 files: magnitude and phase, assuming the first one is magnitude
    if b0_magnitude and not b0_phase:
        b0_maps = nibabel.load(b0_magnitude)
        if b0_maps.shape[-1] == 2:
            b0_maps   = nibabel.load(b0_magnitude)
            voxels    = b0_maps.get_data()
            header    = b0_maps.get_header()
            affine    = b0_maps.get_affine()
            magnitude = nibabel.Nifti1Image(voxels[:,:,:,0], affine, header)
            phase     = nibabel.Nifti1Image(voxels[:,:,:,1], affine, header)
            b0_magnitude = os.path.join(output_directory, "b0_magnitude.nii.gz")
            b0_phase     = os.path.join(output_directory, "b0_phase.nii.gz")
            magnitude.to_filename(b0_magnitude)
            phase    .to_filename(b0_phase)

    # Convert Nifti to Gis
    dwi = nifti_to_gis(dwi, os.path.join(output_directory, "dwi.ima"))

    # Copy bval and bvec files, with homogeneous names
    bval_copy = os.path.join(output_directory, "dwi.bval")
    bvec_copy = os.path.join(output_directory, "dwi.bvec")
    shutil.copyfile(bval, bval_copy)
    shutil.copyfile(bvec, bvec_copy)

    # Convert and rename B0 map(s)
    b0_magnitude = nifti_to_gis(b0_magnitude,
                                os.path.join(output_directory, "b0_magnitude.ima"))

    if b0_phase:
        b0_phase = nifti_to_gis(b0_phase,
                                os.path.join(output_directory, "b0_phase.ima"))
    else:
        b0_phase = None

    return (dwi, bval_copy, bvec_copy, b0_magnitude,
            b0_phase)


def dwi_data_import_and_qspace_sampling(output_directory,
                                        dwi,
                                        bval,
                                        bvec,
                                        manufacturer,
                                        b0_magnitude,
                                        b0_phase   = None,
                                        invertX    = True,
                                        invertY    = False,
                                        invertZ    = False,
                                        subject_id = None):
    """
    Wrapper to Connectomist's "DWI & Q-space" tab.

    Parameters
    ----------
    output_directory: Str, path to Connectomist's output directory.
    dwi:              Str, path to Nifti diffusion-weighted data.
    bvec:             Str, path to .bval file associated to the Nifti.
    bval:             Str, path to .bvec file associated to the Nifti.
    manufacturer:     Str, name of the manufacturer (e.g. "Siemens", "GE",
                      "Philips" or "Bruker").
    invertX:          Bool, if True invert x-axis of diffusion model.
    invertY           Same as invertX for y-axis.
    invertZ:          Same as invertX for z-axis.

    <unit>
        <output name="raw_dwi_directory" type="Directory" description="Path to
            Connectomist output directory."/>

        <input name="output_directory" type="Directory" />
        <input name="dwi"              type="File"      />
        <input name="bval"             type="File"      />
        <input name="bvec"             type="File"      />
        <input name="manufacturer" type="Str" description="Name of the MRI
            manufacturer (e.g. 'Siemens', 'GE', 'Philips' or 'Bruker')." />
        <input name="b0_magnitude"     type="File"      />
        <input name="b0_phase"         type="File"      />
        <input name="invertX" type="Bool" description="If True invert x-axis
            of diffusion model."/>
        <input name="invertY" type="Bool" description="Same as invertX but for
            y-axis."/>
        <input name="invertZ" type="Bool" description="Same as invertX but for
            z-axis."/>
        <input name="subject_id" type="Str" description="Subject's identifier." />
    </unit>
    """

    dwi, bval, bvec, b0_magnitude, b0_phase = \
        gather_and_format_input_files(output_directory,
                                      dwi,
                                      bval,
                                      bvec,
                                      b0_magnitude,
                                      b0_phase)

    algorithm = "DWI-Data-Import-And-QSpace-Sampling"

    # Dict with all parameters for connectomist
    parameters_dict = {
        # Parameters are ordered as they appear in connectomist's GUI

        # ---------------------------------------------------------------------
        # Field: "Diffusion weighted-images"
        "fileNameDwi":        dwi,  # "DW data"
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
        "qSpaceChoice5OrientationFileNames": bvec,

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

    parameters_dict["manufacturer"] = MANUFACTURERS[manufacturer]

    # Read .bval file to infer nb of T2, nb of shells...
    try:
        bvalues = np.loadtxt(bval)
        if set(bvalues) == {0}:  # If only 0s raise Exception
            raise Exception
    except:
        raise BadFileError(bval)

    nb_T2     = np.sum(bvalues == 0)  # nb of volumes where bvalue=0
    bvals_set = set(bvalues) - {0}    # set of non-zero bvalues
    nb_shells = len(bvals_set)

    parameters_dict["numberOfT2"] = nb_T2

    if nb_shells == 1:
        # Spherical single-shell custom
        parameters_dict["qSpaceSamplingType"] = 4
    else:
        raise ConnectomistError("Multiple shell models not handled. "
                                "Path to .bval file: %s" % bval)

    # Check validity of .bvec file.
    # If bvec file does not exist or filled with 0s, raise Error
    if (not os.path.isfile(bvec)) or np.loadtxt(bvec).max() == 0:
        raise BadFileError(bvec)

    parameter_file = create_parameter_file(algorithm, parameters_dict,
                                           output_directory)
    run_connectomist(algorithm, parameter_file, output_directory)

    # Capsul needs the output name to be different from input arguments
    raw_dwi_directory = output_directory

    return raw_dwi_directory
