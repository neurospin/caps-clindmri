#!/usr/bin/env python
# -*- coding: utf-8 -*-

##########################################################################
# NSAp - Copyright (C) CEA, 2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################


from .utils import create_parameter_file, run_connectomist


def dwi_outlier_detection(output_directory,
                          raw_dwi_directory,
                          rough_mask_directory):
    """
    Wrapper to Connectomist's "Outliers" tab.

    Parameters
    ----------
    output_directory:     Str, path to Connectomist output work directory.
    raw_dwi_directory:    Str, path to Connectomist Raw DWI folder.
    rough_mask_directory: Str, path to Connectomist Rough Mask folder.

    <unit>
        <output name="outliers_directory" type="Directory" />

        <input name="output_directory" type="Directory"    />
        <input name="raw_dwi_directory" type="Directory" description="Path to
            Connectomist Raw DWI folder."                  />
        <input name="rough_mask_directory" type="Directory" description="Path
            to Connectomist Rough Mask folder."            />
    </unit>
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

    # Capsul needs the output name to be different from input arguments
    outliers_directory = output_directory
    return outliers_directory
