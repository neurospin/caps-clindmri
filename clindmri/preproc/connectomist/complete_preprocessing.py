#!/usr/bin/env python
# -*- coding: utf-8 -*-

##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

import os

# Wrappers of Connectomist's tabs
from .import_and_qspace_model import dwi_data_import_and_qspace_sampling
from .mask                    import dwi_rough_mask_extraction
from .outliers                import dwi_outlier_detection
from .susceptibility          import dwi_susceptibility_artifact_correction
from .eddy_current_and_motion import dwi_eddy_current_and_motion_correction
from .export_results          import export_results


def complete_preprocessing(output_directory,
                           subject_id,
                           dwi,
                           bval,
                           bvec,
                           manufacturer,
                           delta_TE,
                           partial_fourier_factor,
                           parallel_acceleration_factor,
                           b0_magnitude,
                           b0_phase        = None,
                           invertX         = True,
                           invertY         = False,
                           invertZ         = False,
                           negative_sign   = False,
                           echo_spacing    = None,
                           EPI_factor      = None,
                           b0_field        = 3.0,
                           water_fat_shift = 4.68,
                           delete_steps    = False):
    """
    Function that runs all preprocessing tabs from Connectomist.

    Parameters
    ----------
    output_directory: Str, path to folder where all the preprocessing will be done.
    subject_id:       Str, subject identifier.
    dwi               Str, path to input Nifti DW data.
    bval:             Str, path to Nifti's associated .bval file.
    bvec:             Str, path to Nifti's associated .bval file.
    manufacturer:     Str, manufacturer name (e.g. "Siemens", "GE"...).
    delta_TE:         Float, difference in seconds between the 2 echoes in B0
                      magnitude map acquisition.
    partial_fourier_factor:
                      Float (]0;1]), percentage of k-space plane acquired.
    parallel_acceleration_factor:
                      Int, nb of parallel acquisition in k-space plane.
    b0_magnitude:     Str, path to B0 magnitude map, also contains phase for GE.
    b0_phase:         Str, not for GE, path to B0 phase map.
    invertX:          Bool, if True invert x-axis in diffusion model.
    invertY:          Bool, same as invertX but for y-axis.
    invertZ:          Bool, same as invertX but for z-axis.
    negative_sign:    Bool, if True invert direction of unwarping in
                      susceptibility-distortion correction.
    echo_spacing:     Float, not for Philips, acquisition time in ms between
                      2 centers of 2 consecutively acquired lines in k-space.
    EPI_factor:       Int, nb of echoes after one RF pulse, i.e. echo train length.
    b0_field:         Float, Philips only, B0 field intensity, by default 3.0.
    water_fat_shift:  Float, Philips only, default 4.68 pixels.
    delete_steps:     Bool, if True remove all intermediate files and
                      directories at the end of preprocessing, to keep only
                      4 files:
                      preprocessed Nifti + bval + bvec + outliers.py

    Returns
    -------
    preproc_dwi, preproc_bval, preproc_bvec: Paths to output files.

    <unit>
        <output name="preproc_dwi"                 type="File"      />
        <output name="preproc_bval"                type="File"      />
        <output name="preproc_bvec"                type="File"      />

        <input name="output_directory"             type="Directory" />
        <input name="subject_id"                   type="Str"       />
        <input name="dwi"                          type="File"      />
        <input name="bval"                         type="File"      />
        <input name="bvec"                         type="File"      />
        <input name="manufacturer"                 type="Str"       />
        <input name="delta_TE"                     type="Float"     />
        <input name="partial_fourier_factor"       type="Float"     />
        <input name="parallel_acceleration_factor" type="Float"     />
        <input name="b0_magnitude"                 type="File"      />
        <input name="b0_phase"                     type="File"      />
        <input name="invertX"                      type="Bool"      />
        <input name="invertY"                      type="Bool"      />
        <input name="invertZ"                      type="Bool"      />
        <input name="negative_sign"                type="Bool"      />
        <input name="echo_spacing"                 type="Float"     />
        <input name="EPI_factor"                   type="Int"       />
        <input name="b0_field"                     type="Float"     />
        <input name="water_fat_shift"              type="Float"     />
        <input name="delete_steps"                 type="Bool"      />
    </unit>
    """
    # Step 0 - Create the preprocessing output directory if not existing
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # Step 1 - Import files to Connectomist and choose q-space model
    raw_dwi_directory = os.path.join(output_directory, "01-Import_and_qspace_model")
    dwi_data_import_and_qspace_sampling(raw_dwi_directory,
                                        dwi,
                                        bval,
                                        bvec,
                                        manufacturer,
                                        b0_magnitude,
                                        b0_phase,
                                        invertX,
                                        invertY,
                                        invertZ,
                                        subject_id)

    # Step 2 - Create a brain mask
    rough_mask_directory = os.path.join(output_directory, "02-Rough_mask")
    dwi_rough_mask_extraction(rough_mask_directory,
                              raw_dwi_directory)

    # Step 3 - Detect and correct outlying diffusion slices
    outliers_directory = os.path.join(output_directory, "03-Outliers")
    dwi_outlier_detection(outliers_directory,
                          raw_dwi_directory,
                          rough_mask_directory)

    # Step 4 - Susceptibility correction
    susceptibility_directory = os.path.join(output_directory, "04-Suceptibility")
    dwi_susceptibility_artifact_correction(susceptibility_directory,
                                           raw_dwi_directory,
                                           rough_mask_directory,
                                           outliers_directory,
                                           manufacturer,
                                           delta_TE,
                                           partial_fourier_factor,
                                           parallel_acceleration_factor,
                                           b0_magnitude,
                                           b0_phase,
                                           negative_sign,
                                           echo_spacing,
                                           EPI_factor,
                                           b0_field,
                                           water_fat_shift)

    # Step 5 - Eddy current and motion correction
    eddy_motion_directory = os.path.join(output_directory, "05-Eddy_current_and_motion")
    dwi_eddy_current_and_motion_correction(eddy_motion_directory,
                                           raw_dwi_directory,
                                           rough_mask_directory,
                                           susceptibility_directory)

    # Step 6 - Export result as a Nifti with a .bval and a .bvec
    preproc_dwi, preproc_bval, preproc_bvec = \
        export_results(output_directory,
                       raw_dwi_directory,
                       rough_mask_directory,
                       outliers_directory,
                       susceptibility_directory,
                       eddy_motion_directory,
                       delete_steps)

    return preproc_dwi, preproc_bval, preproc_bvec
