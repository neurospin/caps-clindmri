#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os
import shutil
import numpy as np
import nibabel

from .exceptions import BadFileError
from .utils      import nifti_to_gis, gis_to_nifti, concatenate_volumes

# Wrappers of Connectomist's tabs
from .import_and_qspace_model import dwi_data_import_and_qspace_sampling
from .mask                    import dwi_rough_mask_extraction
from .outliers                import dwi_outlier_detection
from .susceptibility          import dwi_susceptibility_artifact_correction
from .eddy_current_and_motion import dwi_eddy_current_and_motion_correction


# TO BE COMPLETED: CAPSUL
def gather_and_format_input_files(output_directory,
                                  path_nifti,
                                  path_bval,
                                  path_bvec,
                                  path_b0_magnitude,
                                  path_b0_phase = None):
    """
    Gather all files needed to start the preprocessing in the right format
    (Gis for images and B0 maps).

    Parameters
    ----------
    output_directory:  Str, path to directory where to gather all files.
    path_nifti:        Str, path to input Nifti file to be preprocessed.
    path_bval:         Str, path to .bval file associated to Nifti.
    path_bvec:         Str, path to .bvec file associated to Nifti.
    path_b0_magnitude: Str, path to B0 magnitude map, may also contain phase.
    path_b0_phase:     Str, not required if phase is already contained in
                       path_b0_magnitude.

    Returns
    -------
    path_gis, path_bval, path_bvec, path_b0_maps: New paths.
    """

    # Check that files exist
    files = [path_nifti, path_bval, path_bvec, path_b0_magnitude]
    if path_b0_phase:
        files.append(path_b0_phase)

    for path in files:
        if not os.path.isfile(path):
            raise BadFileError(path)

    # Create the directory if not existing
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # If there is only one b0 map file and if this file contains 2 volumes,
    # split it in 2 files: magnitude and phase
    if path_b0_magnitude and not path_b0_phase:
        b0_maps = nibabel.load(path_b0_magnitude)
        if b0_maps.shape[-1] == 2:
            b0_maps   = nibabel.load(path_b0_magnitude)
            voxels    = b0_maps.get_data()
            header    = b0_maps.get_header()
            affine    = b0_maps.get_affine()
            magnitude = nibabel.Nifti1Image(voxels[:,:,:,0], affine, header)
            phase     = nibabel.Nifti1Image(voxels[:,:,:,1], affine, header)
            path_b0_magnitude = os.path.join(output_directory, "b0_magnitude.nii.gz")
            path_b0_phase     = os.path.join(output_directory, "b0_phase.nii.gz")
            magnitude.to_filename(path_b0_magnitude)
            phase    .to_filename(path_b0_phase)

    # Convert Nifti to Gis
    filename_wo_ext = os.path.basename(path_nifti).rsplit(".gz", 1)[0].rsplit(".nii",1)[0]
    path_gis        = os.path.join(output_directory, filename_wo_ext)
    nifti_to_gis(path_nifti, path_gis)

    # Copy bval and bvec files, with homogeneous names
    path_bval_copy = path_gis + ".bval"
    path_bvec_copy = path_gis + ".bvec"
    shutil.copyfile(path_bval, path_bval_copy)
    shutil.copyfile(path_bvec, path_bvec_copy)

    # Convert and rename B0 map(s)
    path_b0_magnitude_gis = os.path.join(output_directory, "b0_magnitude.gis")
    nifti_to_gis(path_b0_magnitude, path_b0_magnitude_gis)

    if path_b0_phase:
        path_b0_phase_gis = os.path.join(output_directory, "b0_phase.gis")
        nifti_to_gis(path_b0_phase, path_b0_phase_gis)
    else:
        path_b0_phase_gis = None

    return (path_gis, path_bval_copy, path_bvec_copy, path_b0_magnitude_gis,
            path_b0_phase_gis)


def export_and_format_corrected_files(eddy_motion_directory,
                                      path_nifti,
                                      path_bval=None,
                                      path_bvec=None):
    """
    After Eddy current and motion correction, export the result as Nifti with
    metafiles (.bval and .bvec with corrrected directions).

    Parameters
    ----------
    eddy_motion_directory: Str, path to Connectomist "Eddy Current & Motion"
                           directory.
    path_nifti: Str, path to output Nifti file.
    path_bval:  Str (optional), path to output .bval file.
                By default same as path_nifti but with .bval extension.
    path_bvec:  Same a path_bval for .bvec file.

    Returns
    -------
    path_nifti, path_bval, path_bvec: (Str, Str, Str) paths to output files.

    <process>
        <return name="path_nifti" type="File" />
        <return name="path_bval"  type="File" />
        <return name="path_bvec"  type="File" />
        <input name="eddy_motion_directory" type="Directory"
               desc="Path to Connectomist "Eddy Current & Motion" directory."/>
        <input name="path_nifti" type="File" desc="Path to output Nifti file."/>
        <input name="path_bval" optional="True" type="File"
               desc="Path to the output .bval file. By default it is the same
                     as path_nifti with .bval extension."/>
        <input name="path_bvec" optional="True" type="File"
               desc="path to the output .bval file. By default it is the same
                     as path_nifti with .bvec extension."/>
    </process>
    """

    # Concatenate preprocessed T2 and preprocessed DW volumes
    path_t2        = os.path.join(eddy_motion_directory, "t2_wo_eddy_current_and_motion")
    path_dw        = os.path.join(eddy_motion_directory, "dw_wo_eddy_current_and_motion")
    path_t2_dw_gis = os.path.join(eddy_motion_directory, "t2_dw_wo_eddy_current_and_motion")
    concatenate_volumes([path_t2, path_dw], path_t2_dw_gis)

    # Convert to Nifti
    gis_to_nifti(path_t2_dw_gis, path_nifti)

    # Load a dict from the metadata file (<nifti path>.minf).
    try:
        path_minf = path_dw + ".ima.minf"
        exec_dict = dict()  # To store variables created by execfile call.
        execfile(path_minf, exec_dict)

        # Verify dict structure
        if "bvalues" not in exec_dict["attributes"]:
            raise Exception

        if "diffusion_gradient_orientations" not in exec_dict["attributes"]:
            raise Exception
    except:
        raise BadFileError(path_minf)

    # Get bvalues and create .bval file
    bvalues = np.array(exec_dict["attributes"]["bvalues"], dtype=np.int)
    bvalues = np.concatenate(([0], bvalues))  # don't forget b=0 associated to T2
    if not path_bval:
        path_bval = path_nifti.rsplit(".gz", 1)[0].rsplit(".nii", 1)[0] + ".bval"
    np.savetxt(path_bval, bvalues, newline=" ", fmt="%d")

    # Get gradient directions and create .bvec file
    directions = np.array(exec_dict["attributes"]["diffusion_gradient_orientations"])
    if not path_bvec:
        path_bvec = path_nifti.rsplit(".gz", 1)[0].rsplit(".nii", 1)[0] + ".bvec"
    np.savetxt(path_bvec, directions.T, fmt="%.10f")

    return path_nifti, path_bval, path_bvec


###############################################################################
# Main function
###############################################################################

# TO BE COMPLETED: capsul
def complete_preprocessing(path_nifti,
                           path_bval,
                           path_bvec,
                           manufacturer,
                           delta_TE,
                           partial_fourier_factor,
                           parallel_acceleration_factor,
                           path_b0_magnitude,
                           path_b0_phase             = None,
                           invertX                   = True,
                           invertY                   = False,
                           invertZ                   = False,
                           negative_sign             = False,
                           echo_spacing              = None,
                           EPI_factor                = None,
                           b0_field                  = 3.0,
                           water_fat_shift_per_pixel = 4.68,
                           output_directory          = None,
                           delete_steps              = False):
    """
    Function that runs all preprocessing tabs from Connectomist.

    Parameters
    ----------
    path_nifti         Str, path to input Nifti DWI data.
    path_bval:         Str, path to Nifti's associated .bval file.
    path_bvec:         Str, path to Nifti's associated .bval file.
    manufacturer:      Str, manufacturer name (e.g. "Siemens", "GE"...).
    delta_TE:          Float, difference in seconds between the 2 echoes in B0
                       magnitude map acquisition.
    partial_fourier_factor: Float (]0;1]), percentage of k-space plane acquired.
    parallel_acceleration_factor: Int, nb of parallel acquisition in k-space plane.
    path_b0_magnitude: Str, path to B0 magnitude map, also contains phase for GE.
    path_b0_phase:     Str, not for GE, path to B0 phase map.
    invertX:           Bool, if True invert x-axis in diffusion model.
    invertY:           Bool, same as invertX but for y-axis.
    invertZ:           Bool, same as invertX but for z-axis.
    negative_sign:     Bool, if True invert direction of unwarping in
                       susceptibility-distortion correction.
    echo_spacing:      Float, not for Philips, acquisition time in ms between
                       2 centers of 2 consecutively acquired lines in k-space.
    EPI_factor:        Int, nb of echoes after one excitation (90 degrees),
                       i.e. echo train length.
    b0_field:          Float, Philips only, B0 field intensity, by default 3.0.
    water_fat_shift_per_pixel: Float, Philips only, default 4.68Hz
    output_directory:  Str, path to folder where all the preprocessing will be
                       done.
    delete_steps:      Bool, if True remove all intermediate files and
                       directories at the end of preprocessing, to keep only
                       the "07-Preprocessed" directory with 4 files inside:
                       preprocessed Nifti + bval + bvec + outliers.py


    Returns
    -------
    path_nifti, path_bval, path_bvec: Paths to output files.

    <process>
    </process>
    """
    # Step 0 - Create the preprocessing output directory if not existing
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # Step 1 - Bring all files in Gis format in a directory: "01-Input_Data"
    input_directory = os.path.join(output_directory, "01-Input_Data")
    path_gis, path_bval, path_bvec, path_b0_magnitude, path_b0_phase = \
        gather_and_format_input_files(input_directory,
                                      path_nifti,
                                      path_bval,
                                      path_bvec,
                                      path_b0_magnitude,
                                      path_b0_phase)

    # Step 2 - Import files to Connectomist and choose diffusion model
    raw_dwi_directory = os.path.join(output_directory, "02-Raw_Dwi")
    dwi_data_import_and_qspace_sampling(path_gis,
                                        path_bval,
                                        path_bvec,
                                        raw_dwi_directory,
                                        manufacturer,
                                        invertX,
                                        invertY,
                                        invertZ)

    # Step 3 - Create a brain mask
    rough_mask_directory = os.path.join(output_directory, "03-Rough_Mask")
    dwi_rough_mask_extraction(raw_dwi_directory, rough_mask_directory)

    # Step 4 - Detect and filter outlying diffusion slices
    outliers_directory = os.path.join(output_directory, "04-Outliers")
    dwi_outlier_detection(raw_dwi_directory,
                          rough_mask_directory,
                          outliers_directory)

    # Step 5 - Susceptibility correction
    susceptibility_directory = os.path.join(output_directory, "05-Susceptiblity")
    dwi_susceptibility_artifact_correction(raw_dwi_directory,
                                           rough_mask_directory,
                                           outliers_directory,
                                           susceptibility_directory,
                                           manufacturer,
                                           delta_TE,
                                           partial_fourier_factor,
                                           parallel_acceleration_factor,
                                           path_b0_magnitude,
                                           path_b0_phase,
                                           negative_sign,
                                           echo_spacing,
                                           EPI_factor,
                                           b0_field,
                                           water_fat_shift_per_pixel)

    # Step 6 - Eddy current and motion correction
    eddy_motion_directory = os.path.join(output_directory, "06-Eddy_Current_And_Motion")
    dwi_eddy_current_and_motion_correction(raw_dwi_directory,
                                           rough_mask_directory,
                                           outliers_directory,
                                           eddy_motion_directory)

    # Step 7 - Export result as a Nifti with a .bval and a .bvec
    preprocessed_directory = os.path.join(output_directory, "07-Preprocessed")
    if not os.path.isdir(preprocessed_directory):
        os.mkdir(preprocessed_directory)
    nifti_name = os.path.basename(path_gis).rsplit(".ima", 1)[0] + ".nii"
    path_nifti = os.path.join(preprocessed_directory, nifti_name)
    path_nifti, path_bval, path_bvec = \
        export_and_format_corrected_files(eddy_motion_directory, path_nifti)

    # Delete intermediate files and directories if requested
    if delete_steps:
        path_outliers_py = os.path.join(outliers_directory, "outliers.py")
        shutil.move(path_outliers_py, preprocessed_directory)

        intermediate_directories = [input_directory, raw_dwi_directory,
                                    rough_mask_directory, outliers_directory,
                                    susceptibility_directory,
                                    eddy_motion_directory]
        for directory in intermediate_directories:
            shutil.rmtree(directory)

        # Remove useless metadata file
        os.remove(path_nifti + ".minf")

    return path_nifti, path_bval, path_bvec
