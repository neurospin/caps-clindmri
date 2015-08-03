#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

from .manufacturers import MANUFACTURERS
from .exceptions    import BadManufacturerNameError, MissingParametersError
from .utils         import create_parameter_file, run_connectomist


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
                          if args_map[manufacturer][p] is None]

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
