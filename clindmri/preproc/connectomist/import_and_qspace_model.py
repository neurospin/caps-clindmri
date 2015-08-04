#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
import numpy as np

from .manufacturers import MANUFACTURERS
from .exceptions    import (ConnectomistError, BadManufacturerNameError,
                            BadFileError)
from .utils         import create_parameter_file, run_connectomist


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
        raise ConnectomistError("Multiple shell models not handled.")

    # Check validity of .bvec file.
    # If bvec file does not exist or filled with 0s, raise Error
    if (not os.path.isfile(path_bvec)) or np.loadtxt(path_bvec).max() == 0:
        raise BadFileError(path_bvec)

    parameter_file = create_parameter_file(algorithm_name,
                                           parameters_value,
                                           output_directory)
    run_connectomist(algorithm_name, parameter_file)

    return output_directory
