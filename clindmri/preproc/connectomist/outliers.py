#!/usr/bin/env python
# -*- coding: utf-8 -*-

##########################################################################
# NSAp - Copyright (C) CEA, 2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################


from .utils import create_parameter_file, run_connectomist


def dwi_outlier_detection(outdir,
                          raw_dwi_dir,
                          rough_mask_dir):
    """
    Wrapper to Connectomist's "Outliers" tab.

    Parameters
    ----------
    outdir:     Str, path to Connectomist output work directory.
    raw_dwi_dir:    Str, path to Connectomist Raw DWI folder.
    rough_mask_dir: Str, path to Connectomist Rough Mask folder.

    <unit>
        <output name="outliers_dir"  type="Directory" />

        <input name="outdir"         type="Directory" />
        <input name="raw_dwi_dir"    type="Directory" />
        <input name="rough_mask_dir" type="Directory" />
    </unit>
    """

    algorithm = "DWI-Outlier-Detection"

    parameters_dict = {'rawDwiDirectory':        raw_dwi_dir,
                       'roughMaskDirectory':  rough_mask_dir,
                       'outputWorkDirectory':         outdir,
                       '_subjectName':                    '',
                       'discardedOrientationList':        '',
                       'outlierFactor':                  3.0}

    parameter_file = create_parameter_file(algorithm, parameters_dict,
                                           outdir)
    run_connectomist(algorithm, parameter_file, outdir)

    # Capsul needs the output name to be different from input arguments
    outliers_dir = outdir

    return outliers_dir
