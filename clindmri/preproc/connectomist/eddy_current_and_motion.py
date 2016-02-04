#!/usr/bin/env python
# -*- coding: utf-8 -*-

##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

import os
import numpy as np

from .exceptions import BadFileError
from .utils      import ptk_gis_to_nifti, ptk_concatenate_volumes

from .utils import create_parameter_file, run_connectomist


def dwi_eddy_current_and_motion_correction(outdir,
                                           raw_dwi_dir,
                                           rough_mask_dir,
                                           corrected_dir):
    """
    Wrapper to Connectomist's "Eddy current & motion" tab.

    Parameters
    ----------
    outdir:         Str, path to Connectomist output work directory.
    raw_dwi_dir:    Str, path to Connectomist Raw DWI directory.
    rough_mask_dir: Str, path to Connectomist Rough Mask directory.
    corrected_dir:  Str, path to Connectomist Susceptibility or Outlier directory.
                    Depending whether you make Eddy Current correction before
                    or after susceptibility correction.

    <unit>
        <output name="eddy_motion_dir" type="Directory" />

        <input name="outdir" type="Directory" description="Path to
            Connectomist output work directory."/>
        <input name="raw_dwi_dir" type="Directory" description="Path to
            Connectomist Raw DWI directory."/>
        <input name="rough_mask_dir" type="Directory" description="Path
            to Connectomist Rough Mask directory."/>
        <input name="corrected_dir" type="Directory" description="
            Connectomist Susceptibility or Outlier directory."/>
    </unit>
    """

    algorithm = "DWI-Eddy-Current-And-Motion-Correction"

    parameters_dict = {
        # ---------------------------------------------------------------------
        # Used parameters
        "rawDwiDirectory":          raw_dwi_dir,
        "roughMaskDirectory":    rough_mask_dir,
        "correctedDwiDirectory":  corrected_dir,
        "outputWorkDirectory":           outdir,
        "eddyCurrentCorrection":              2,
        "motionCorrection":                   1,
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
                                           outdir)
    run_connectomist(algorithm, parameter_file, outdir)

    # Capsul needs the output name to be different from input arguments
    eddy_motion_dir = outdir

    return eddy_motion_dir


def export_eddy_motion_results_to_nifti(eddy_motion_dir,
                                        outdir   = None,
                                        filename = "dwi"):
    """
    After Connectomist has done Eddy current and motion correction, convert
    the result to Nifti with bval/bvec files (bvec with corrrected directions).

    Parameters
    ----------
    eddy_motion_dir: Str, path to the Connectomist "Eddy current and motion"
                     directory.
    outdir:          Str, path to directory where to output:
                         <outdir>/<filename>.nii.gz
                         <outdir>/<filename>.bval
                         <outdir>/<filename>.bvec
                     By default <outdir> is <eddy_motion_dir>.
    filename:        Str, to change output filenames, by default "dwi".
    """
    # Step 0 - Set outdir path and check directory existence
    if outdir is None:
        outdir = eddy_motion_dir
    else:
        if not os.path.isdir(outdir):  # If outdir does not exist, create it
            os.makedirs(outdir)

    # Step 1 - Concatenate preprocessed T2 and preprocessed DW volumes

    # Set input and output paths (Gis files) without extension (.ima)
    t2    = os.path.join(eddy_motion_dir, "t2_wo_eddy_current_and_motion.ima")
    dw    = os.path.join(eddy_motion_dir, "dw_wo_eddy_current_and_motion.ima")
    t2_dw = os.path.join(eddy_motion_dir, "t2_dw_wo_eddy_current_and_motion.ima")

    # Check existence of input files
    for path in t2, dw:
        if not os.path.isfile(path):
            raise BadFileError(path)

    # Apply concatenation: result is a Gis file
    ptk_concatenate_volumes([t2, dw], t2_dw)

    # Step 2 - Convert to Nifti
    dwi = ptk_gis_to_nifti(t2_dw, os.path.join(outdir, "%s.nii.gz" % filename))

    # Step 3 - Create .bval and .bvec (with corrected directions)

    # The new directions of gradients (modified by the Eddy current and motion
    # correction) are found in the .ima.minf (Gis format) file associated to
    # diffusion weighted data.
    try:
        # The .ima.minf is a text file that defines a python dict
        path_minf = dw + ".minf"
        exec_dict = dict()  # To store variables created by execfile() call.
        execfile(path_minf, exec_dict)

        # Verify dict structure by checking 2 keys.
        if "bvalues" not in exec_dict["attributes"]:
            raise Exception
        if "diffusion_gradient_orientations" not in exec_dict["attributes"]:
            raise Exception
    except:
        raise BadFileError(path_minf)

    # Get bvalues and create .bval file
    bvalues = np.array(exec_dict["attributes"]["bvalues"], dtype=np.int)

    # Add 0 in bvalues for b=0 associated to T2
    bvalues = np.concatenate(([0], bvalues))

    # Create "dwi.bval"
    bval = os.path.join(outdir, "%s.bval" % filename)
    np.savetxt(bval, bvalues, newline=" ", fmt="%d")

    # Get gradient directions to create .bvec file
    directions = np.array(exec_dict["attributes"]["diffusion_gradient_orientations"])

    # Normalize vectors
    magnitudes = np.linalg.norm(directions, axis=1)
    for i, vector in enumerate(directions):
        directions[i, :] = directions[i, :] / magnitudes[i]

    # Add null vector for direction corresponding to bvalue=0 volume
    directions = np.concatenate(([[0, 0, 0]], directions))

    # Create "dwi.bvec"
    bvec = os.path.join(outdir, "%s.bvec" % filename)
    np.savetxt(bvec, directions.T, fmt="%.10f")

    return dwi, bval, bvec

