#!/usr/bin/env python
# -*- coding: utf-8 -*-

##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

import os
import shutil

# Wrappers of Connectomist's tabs
from .import_and_qspace_model import dwi_data_import_and_qspace_sampling
from .mask                    import dwi_rough_mask_extraction
from .outliers                import dwi_outlier_detection
from .susceptibility          import dwi_susceptibility_artifact_correction
from .eddy_current_and_motion import (dwi_eddy_current_and_motion_correction,
                                      export_eddy_motion_results_to_nifti)


def complete_preprocessing(outdir,
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
    outdir:          Str, path to folder where all the preprocessing will be done.
    subject_id:      Str, subject identifier.
    dwi              Str, path to input Nifti DW data.
    bval:            Str, path to Nifti's associated .bval file.
    bvec:            Str, path to Nifti's associated .bval file.
    manufacturer:    Str, manufacturer name (e.g. "Siemens", "GE"...).
    delta_TE:        Float, difference in seconds between the 2 echoes in B0
                     magnitude map acquisition.
    partial_fourier_factor:
                     Float (]0;1]), percentage of k-space plane acquired.
    parallel_acceleration_factor:
                     Int, nb of parallel acquisition in k-space plane.
    b0_magnitude:    Str, path to B0 magnitude map, also contains phase for GE.
    b0_phase:        Str, not for GE, path to B0 phase map.
    invertX:         Bool, if True invert x-axis in diffusion model.
    invertY:         Bool, same as invertX but for y-axis.
    invertZ:         Bool, same as invertX but for z-axis.
    negative_sign:   Bool, if True invert direction of unwarping in
                     usceptibility-distortion correction.
    echo_spacing:    Float, not for Philips, acquisition time in ms between
                     2 centers of 2 consecutively acquired lines in k-space.
    EPI_factor:      Int, nb of echoes after one RF pulse, i.e. echo train length.
    b0_field:        Float, Philips only, B0 field intensity, by default 3.0.
    water_fat_shift: Float, Philips only, default 4.68 pixels.
    delete_steps:    Bool, if True remove all intermediate files and
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

        <input name="outdir"                       type="Directory" />
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
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Step 1 - Import files to Connectomist and choose q-space model
    raw_dwi_dir = os.path.join(outdir, "01-Import_and_qspace_model")
    dwi_data_import_and_qspace_sampling(raw_dwi_dir,
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
    rough_mask_dir = os.path.join(outdir, "02-Rough_mask")
    dwi_rough_mask_extraction(rough_mask_dir, raw_dwi_dir)

    # Step 3 - Detect and correct outlying diffusion slices
    outliers_dir = os.path.join(outdir, "03-Outliers")
    dwi_outlier_detection(outliers_dir, raw_dwi_dir, rough_mask_dir)

    # Step 4 - Susceptibility correction
    susceptibility_dir = os.path.join(outdir, "04-Suceptibility")
    dwi_susceptibility_artifact_correction(susceptibility_dir,
                                           raw_dwi_dir,
                                           rough_mask_dir,
                                           outliers_dir,
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
    eddy_motion_dir = os.path.join(outdir, "05-Eddy_current_and_motion")
    dwi_eddy_current_and_motion_correction(eddy_motion_dir,
                                           raw_dwi_dir,
                                           rough_mask_dir,
                                           susceptibility_dir)

    # Step 6 - Export result as a Nifti with a .bval and a .bvec
    preproc_dwi, preproc_bval, preproc_bvec = \
        export_eddy_motion_results_to_nifti(eddy_motion_dir,
                                            outdir   = outdir,
                                            filename = "dwi")

    # Step 7 - Export outliers.py
    path_outliers_py = os.path.join(outliers_dir, "outliers.py")
    shutil.copy(path_outliers_py, outdir)

    # Step 8 - Delete intermediate files and directories if requested,
    if delete_steps:
        intermediate_dirs = [raw_dwi_dir, rough_mask_dir, outliers_dir,
                             susceptibility_dir, eddy_motion_dir]
        for directory in intermediate_dirs:
            shutil.rmtree(directory)

    return preproc_dwi, preproc_bval, preproc_bvec
