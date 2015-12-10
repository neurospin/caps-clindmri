#!/usr/bin/env python
# -*- coding: utf-8 -*-

##########################################################################
# NSAp - Copyright (C) CEA, 2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################


from .utils import create_parameter_file, run_connectomist


def dwi_rough_mask_extraction(output_directory, raw_dwi_directory):
    """
    Wrapper to Connectomist's "Rough mask" tab.

    Parameters
    ----------
    output_directory:  Str, path to Connectomist output work directory.
    raw_dwi_directory: Str, path to Connectomist Raw DWI folder.

    <unit>
        <output name="rough_mask_directory" type="Directory" description="Path
            to Connectomist output work directory"/>

        <input name="output_directory" type="Directory" description="Path to
            Connectomist output work directory"/>
        <input name="raw_dwi_directory" type="Directory" description="Path to
            Connectomist Raw DWI folder."/>
    </unit>
    """

    algorithm = "DWI-Rough-Mask-Extraction"

    parameters_dict = {
        # ---------------------------------------------------------------------
        # Used parameters
        'outputWorkDirectory': output_directory,
        'rawDwiDirectory':    raw_dwi_directory,
        # ---------------------------------------------------------------------
        # Parameters not used/handled by the code
        '_subjectName': '',
        'anatomy': '',
        'dwToT1RegistrationParameter': {
            'applySmoothing':                                1,
            'floatingLowerThreshold':                      0.0,
            'initialParametersRotationX':                    0,
            'initialParametersRotationY':                    0,
            'initialParametersRotationZ':                    0,
            'initialParametersScalingX':                   1.0,
            'initialParametersScalingY':                   1.0,
            'initialParametersScalingZ':                   1.0,
            'initialParametersShearingXY':                 0.0,
            'initialParametersShearingXZ':                 0.0,
            'initialParametersShearingYZ':                 0.0,
            'initialParametersTranslationX':                 0,
            'initialParametersTranslationY':                 0,
            'initialParametersTranslationZ':                 0,
            'initializeCoefficientsUsingCenterOfGravity': True,
            'levelCount':                                   32,
            'maximumIterationCount':                      1000,
            'maximumTestGradient':                      1000.0,
            'maximumTolerance':                           0.01,
            'optimizerName':                                 0,
            'optimizerParametersRotationX':                  5,
            'optimizerParametersRotationY':                  5,
            'optimizerParametersRotationZ':                  5,
            'optimizerParametersScalingX':                0.05,
            'optimizerParametersScalingY':                0.05,
            'optimizerParametersScalingZ':                0.05,
            'optimizerParametersShearingXY':              0.05,
            'optimizerParametersShearingXZ':              0.05,
            'optimizerParametersShearingYZ':              0.05,
            'optimizerParametersTranslationX':              30,
            'optimizerParametersTranslationY':              30,
            'optimizerParametersTranslationZ':              30,
            'referenceLowerThreshold':                     0.0,
            'resamplingOrder':                               1,
            'similarityMeasureName':                         1,
            'stepSize':                                    0.1,
            'stoppingCriterionError':                     0.01,
            'subSamplingMaximumSizes':                    '64',
            'transform3DType':                               0
        },
        'maskClosingRadius':        0.0,
        'maskDilationRadius':       4.0,
        'morphologistBrainMask':     '',
        'noiseThresholdPercentage': 2.0,
        'strategyRoughMaskFromT1':    0,
        'strategyRoughMaskFromT2':    1
    }

    parameter_file = create_parameter_file(algorithm, parameters_dict,
                                           output_directory)
    run_connectomist(algorithm, parameter_file, output_directory)

    # Capsul needs the output name to be different from input arguments
    rough_mask_directory = output_directory

    return rough_mask_directory
