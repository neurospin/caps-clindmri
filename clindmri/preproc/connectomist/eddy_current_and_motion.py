#!/usr/bin/env python
# -*- coding: utf-8 -*-

##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

from .utils import create_parameter_file, run_connectomist


def dwi_eddy_current_and_motion_correction(output_directory,
                                           raw_dwi_directory,
                                           rough_mask_directory,
                                           susceptibility_directory):
    """
    Wrapper to Connectomist's "Eddy current & motion" tab.

    Parameters
    ----------
    output_directory:        Str, path to Connectomist output work directory.
    raw_dwi_directory:       Str, path to Connectomist Raw DWI directory.
    rough_mask_directory:    Str, path to Connectomist Rough Mask directory.
    susceptbility_directory: Str, path to Connectomist Susceptibility directory.

    <unit>
        <output name="eddy_motion_directory" type="Directory" />

        <input name="output_directory" type="Directory" description="Path to
            Connectomist output work directory."/>
        <input name="raw_dwi_directory" type="Directory" description="Path to
            Connectomist Raw DWI directory."/>
        <input name="rough_mask_directory" type="Directory" description="Path
            to Connectomist Rough Mask directory."/>
        <input name="susceptibility_directory" type="Directory" description="
            Path to Connectomist Susceptibility directory."/>
    </unit>
    """

    algorithm = "DWI-Eddy-Current-And-Motion-Correction"

    parameters_dict = {
        # ---------------------------------------------------------------------
        # Used parameters
        "rawDwiDirectory":              raw_dwi_directory,
        "roughMaskDirectory":        rough_mask_directory,
        "correctedDwiDirectory": susceptibility_directory,
        "outputWorkDirectory":           output_directory,
        "eddyCurrentCorrection":                        2,
        "motionCorrection":                             1,
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

    parameter_file = create_parameter_file(algorithm, parameters_dict,
                                           output_directory)
    run_connectomist(algorithm, parameter_file, output_directory)

    # Capsul needs the output name to be different from input arguments
    eddy_motion_directory = output_directory

    return eddy_motion_directory
