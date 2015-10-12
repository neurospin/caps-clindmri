
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

import os
import shutil
import numpy as np

from .exceptions import BadFileError
from .utils      import gis_to_nifti, concatenate_volumes

def export_results(export_directory,
                   raw_dwi_directory,
                   rough_mask_directory,
                   outliers_directory,
                   susceptibility_directory,
                   eddy_motion_directory,
                   delete_steps):
    """
    After all preprocessing steps from Connectomist have be done, export the
    result as Nifti with metafiles (.bval and .bvec with corrrected directions).
    outliers.py is also exported.
    Optionnaly the directories of preprocessing steps can be deleted, to save
    disk space.

    Parameters
    ----------
    export_directory:        Str, path to directory where to export the files:
                             "dwi.nii.gz", "dwi.bval", "dwi.bvec" and "outliers.py".
    raw_dwi_directory:       Str, path to Connectomist "Import_and_qspace_model"
                             directory.
    rough_mask_directory:    Str, path to Connectomist "Rough mask" directory.
    outliers_directory:      Str, path to Connetomist "Outliers" directory.
    susceptibility_directory Str, path to Connectomist "Susceptibility" directory.
    eddy_motion_directory:   Str, path to Connectomist "Eddy current and motion"
                             directory.
    delete_steps             Bool, if True delete the Connectomist directories.
                             By default False.

    Returns
    -------
    path_dwi, path_bval, path_bvec: (Str, Str, Str) paths to output files.

    <unit>
        <output name="path_dwi"                type="File"      />
        <output name="path_bval"               type="File"      />
        <output name="path_bvec"               type="File"      />

        <input name="export_directory"         type="Directory" />
        <input name="raw_dwi_directory"        type="Directory" />
        <input name="rough_mask_directory"     type="Directory" />
        <input name="outliers_directory"       type="Directory" />
        <input name="susceptibility_directory" type="Directory" />
        <input name="eddy_motion_directory"    type="Directory" />
        <input name="delete_steps"             type="Bool"      />
    </unit>
    """
    # Step 0 - If export directory does not exist, create it
    if not os.path.isdir(export_directory):
        os.makedirs(export_directory)

    # Step 1 - Concatenate preprocessed T2 and preprocessed DW volumes

    # Set input and output paths (Gis files) without extension (.ima)
    path_t2    = os.path.join(eddy_motion_directory, "t2_wo_eddy_current_and_motion.ima")
    path_dw    = os.path.join(eddy_motion_directory, "dw_wo_eddy_current_and_motion.ima")
    path_t2_dw = os.path.join(eddy_motion_directory, "t2_dw_wo_eddy_current_and_motion.ima")

    # Check existence of input files
    for path in path_t2, path_dw:
        if not os.path.isfile(path):
            raise BadFileError(path)

    # Apply concatenation: result is a Gis file
    concatenate_volumes([path_t2, path_dw], path_t2_dw)

    # Step 2 - Convert to Nifti
    path_dwi = gis_to_nifti(path_t2_dw, os.path.join(export_directory, "dwi.nii.gz"))

    # Step 3 - Create .bval and .bvec (with corrected directions)

    # The new directions of gradients (modified by the Eddy current and motion
    # correction) are found in the .ima.minf (Gis format) file associated to
    # diffusion weighted data.
    try:
        # The .ima.minf is a text file that defines a python dict
        path_minf = path_dw + ".minf"
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
    path_bval = os.path.join(export_directory, "dwi.bval")
    np.savetxt(path_bval, bvalues, newline=" ", fmt="%d")

    # Get gradient directions and create .bvec file
    directions = np.array(exec_dict["attributes"]["diffusion_gradient_orientations"])

    # Normalize vectors
    norms      = np.linalg.norm(directions, axis=0)
    directions = directions / norms

    # Add null vector for direction corresponding to b=0
    directions = np.concatenate(([[0, 0, 0]], directions))

    # Create "dwi.bvec"
    path_bvec = os.path.join(export_directory, "dwi.bvec")
    np.savetxt(path_bvec, directions.T, fmt="%.10f")
    
    # Export outliers.py
    path_outliers_py = os.path.join(outliers_directory, "outliers.py")
    shutil.copy(path_outliers_py, export_directory)
    
    # Delete intermediate files and directories if requested,
    if delete_steps:
        intermediate_directories = [raw_dwi_directory,
                                    rough_mask_directory,
                                    outliers_directory,
                                    susceptibility_directory,
                                    eddy_motion_directory]
        for directory in intermediate_directories:
            shutil.rmtree(directory)

    return path_dwi, path_bval, path_bvec
