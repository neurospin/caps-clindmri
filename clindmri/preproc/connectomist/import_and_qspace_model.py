##########################################################################
# NSAp - Copyright (C) CEA, 2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

# System import
import os
import shutil
import numpy as np
import nibabel

# Clindmri import
from clindmri.extensions.connectomist.manufacturers import MANUFACTURERS
from clindmri.extensions.connectomist.exceptions import ConnectomistError
from clindmri.extensions.connectomist.exceptions import (
    ConnectomistBadManufacturerNameError)
from clindmri.extensions.connectomist.exceptions import (
    ConnectomistBadFileError)
from clindmri.extensions.connectomist import ConnectomistWrapper
from .utils import ptk_nifti_to_gis


def gather_and_format_input_files(outdir,
                                  dwi,
                                  bval,
                                  bvec,
                                  b0_magnitude=None,
                                  b0_phase=None):
    """
    Gather all files needed to start the preprocessing in the right format
    (Gis format for images and B0 maps).

    Parameters
    ----------
    outdir: str
        path to directory where to gather all files.
    dwi: str
        path to input Nifti file to be preprocessed.
    bval: str
        path to .bval file associated to Nifti.
    bvec: str
        path to .bvec file associated to Nifti.
    b0_magnitude: str (optional, default None)
        path to B0 magnitude map, may also contain phase.
        Required only if you want to make fieldmap-based
        correction of susceptibility distortions.
    b0_phase: str (optional, default None)
        not required if phase is already contained in
        b0_magnitude or if you don't want to make fieldmap-based
        correction of susceptibility distortions.

    Returns
    -------
    raw_dwi_dir:
        Connectomist's import data directory with file in Gis format.
    """
    # Check that files exist
    files = [dwi, bval, bvec]
    if b0_magnitude is not None:
        files.append(b0_magnitude)
    if b0_phase is not None:
        files.append(b0_phase)
    for path in files:
        if not os.path.isfile(path):
            raise ConnectomistBadFileError(path)

    # Create the directory if not existing
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # If there is only one b0 map file and if this file contains 2 volumes,
    # split it in 2 files: magnitude and phase, assuming the first one is
    # magnitude
    if b0_magnitude and not b0_phase:
        b0_maps = nibabel.load(b0_magnitude)
        if b0_maps.shape[-1] == 2:
            voxels = b0_maps.get_data()
            header = b0_maps.get_header()
            affine = b0_maps.get_affine()
            magnitude = nibabel.Nifti1Image(voxels[..., 0], affine, header)
            phase = nibabel.Nifti1Image(voxels[..., 1], affine, header)
            b0_magnitude = os.path.join(outdir, "b0_magnitude.nii.gz")
            b0_phase = os.path.join(outdir, "b0_phase.nii.gz")
            magnitude.to_filename(b0_magnitude)
            phase.to_filename(b0_phase)

    # Convert Nifti to Gis
    dwi = ptk_nifti_to_gis(dwi, os.path.join(outdir, "dwi.ima"))

    # Copy bval and bvec files, with homogeneous names
    bval_copy = os.path.join(outdir, "dwi.bval")
    bvec_copy = os.path.join(outdir, "dwi.bvec")
    shutil.copyfile(bval, bval_copy)
    shutil.copyfile(bvec, bvec_copy)

    # Convert and rename B0 map(s) if they are given
    if b0_magnitude is not None:
        b0_magnitude = ptk_nifti_to_gis(
            b0_magnitude, os.path.join(outdir, "b0_magnitude.ima"))

    if b0_phase is not None:
        b0_phase = ptk_nifti_to_gis(
            b0_phase, os.path.join(outdir, "b0_phase.ima"))

    return dwi, bval_copy, bvec_copy, b0_magnitude, b0_phase


def dwi_data_import_and_qspace_sampling(
        outdir,
        dwi,
        bval,
        bvec,
        manufacturer,
        invertX=True,
        invertY=False,
        invertZ=False,
        subject_id=None,
        b0_magnitude=None,
        b0_phase=None,
        nb_tries=10,
        path_connectomist=(
            "/i2bm/local/Ubuntu-14.04-x86_64/ptk/bin/connectomist")):
    """
    Wrapper to Connectomist's "DWI & Q-space" tab.

    Parameters
    ----------
    outdir: str
        path to Connectomist's output directory.
    dwi: str
        path to Nifti diffusion-weighted data.
    bvec: str
        path to .bval file associated to the Nifti.
    bval: str
        path to .bvec file associated to the Nifti.
    manufacturer: str
        name of the manufacturer (e.g. "Siemens", "GE", "Philips" or "Bruker").
    invertX: bool (optional, default True)
        if True invert x-axis of diffusion model.
    invertY, invertZ: bool (optional, default False)
        if True invert y or z-axis of diffusion model.
    b0_magnitude: str (optional, default None)
        path to the magnitude fieldmap (if fieldmap-based correction of
        susceptibility distortions is to be used).
    b0_phase: str (optional, default None)
        path to phase fieldmap (if fieldmap-based correction of susceptibility
        distortions is to be used).
    nb_tries: int (optional, default 10)
        nb of times to try an algorithm if it fails.
        It often crashes when running in parallel. The reason
        why it crashes is unknown.
    path_connectomist: str (optional)
        path to the Connectomist executable.

    Returns
    -------
    outdir: str
        path to Connectomist's output directory.

    <unit>
        <output name="raw_dwi_dir" type="Directory" description="Path to
            Connectomist output directory."/>

        <input name="outdir"       type="Directory" />
        <input name="dwi"          type="File"      />
        <input name="bval"         type="File"      />
        <input name="bvec"         type="File"      />
        <input name="manufacturer" type="Str"  description="'Siemens', 'GE', 'Philips' or 'Bruker'" />
        <input name="invertX"      type="Bool" description="If True invert x-axis of diffusion model."/>
        <input name="invertY"      type="Bool" description="Same as invertX but for y-axis."/>
        <input name="invertZ"      type="Bool" description="Same as invertX but for z-axis."/>
        <input name="subject_id"   type="Str"       />
        <input name="b0_magnitude" type="File"      />
        <input name="b0_phase"     type="File"      />
    </unit>
    """
    # Create the Connectomist's import data directory
    (dwi, bval, bvec, b0_magnitude,
     b0_phase) = gather_and_format_input_files(outdir,
                                               dwi,
                                               bval,
                                               bvec,
                                               b0_magnitude,
                                               b0_phase)

    # Dict with all parameters for connectomist
    algorithm = "DWI-Data-Import-And-QSpace-Sampling"
    parameters_dict = {
        # Parameters are ordered as they appear in connectomist's GUI

        # ---------------------------------------------------------------------
        # Field: "Diffusion weighted-images"
        "fileNameDwi":        dwi,  # "DW data"
        "sliceAxis":            2,  # "Slice axis", default "Z-axis"
        "phaseAxis":            1,  # "Phase axis", default "Y-axis"
        "manufacturer":      None,

        # Subfield: "Advanced parameters"
        "flipAlongX":           0,  # "Flip data along x"
        "flipAlongY":           0,
        "flipAlongZ":           0,
        "numberOfDiscarded":    0,  # "#discarded images at beginning"
        "numberOfT2":        None,  # "#T2"
        "numberOfRepetitions":  1,  # "#repetitions"
        # ---------------------------------------------------------------------
        # Field: "Rotation of field of view", default is identity matrix
        "qSpaceTransform_xx": 1.0,
        "qSpaceTransform_xy": 0.0,
        "qSpaceTransform_xz": 0.0,
        "qSpaceTransform_yx": 0.0,
        "qSpaceTransform_yy": 1.0,
        "qSpaceTransform_yz": 0.0,
        "qSpaceTransform_zx": 0.0,
        "qSpaceTransform_zy": 0.0,
        "qSpaceTransform_zz": 1.0,
        # ---------------------------------------------------------------------
        # Field: "Q-space sampling"
        "qSpaceSamplingType":     4,  # default "spherical single-shell custom"
        "qSpaceChoice5BValue": 1300,
        "qSpaceChoice5OrientationFileNames": bvec,

        # Apparently Connectomist uses 2 as True, and 0 as False.
        "invertXAxis": 2 if invertX else 0,
        "invertYAxis": 2 if invertY else 0,
        "invertZAxis": 2 if invertZ else 0,

        # In this field but not used/handled parameters
        "qSpaceChoice1MaximumBValue":       1000,  # case Cartesian
        "qSpaceChoice2BValue":              1000,
        "qSpaceChoice3BValue":              1000,
        "qSpaceChoice4BValue":              1000,
        "qSpaceChoice6BValues":               "",
        "qSpaceChoice7BValues":               "",
        "qSpaceChoice8BValues":               "",
        "qSpaceChoice9BValues":               "",
        "qSpaceChoice10BValues":              "",
        "qSpaceChoice11BValues":              "",
        "qSpaceChoice12BValues":              "",
        "qSpaceChoice13BValues":              "",
        "qSpaceChoice1NumberOfSteps":         11,
        "qSpaceChoice2NumberOfOrientations":   6,
        "qSpaceChoice3NumberOfOrientations":   6,
        "qSpaceChoice4NumberOfOrientations":   6,
        "qSpaceChoice6NumberOfOrientations":   6,
        "qSpaceChoice7NumberOfOrientations":   6,
        "qSpaceChoice8NumberOfOrientations":   6,
        "qSpaceChoice9OrientationFileNames":  "",
        "qSpaceChoice10NumberOfOrientations": "",
        "qSpaceChoice11NumberOfOrientations": "",
        "qSpaceChoice12NumberOfOrientations": "",
        "qSpaceChoice13OrientationFileNames": "",
        # ---------------------------------------------------------------------
        # Field: "Diffusion time (in ms)"
        "diffusionTime": 1.0,
        # ---------------------------------------------------------------------
        # Field: "Work directory"
        "outputWorkDirectory": outdir,
        # ---------------------------------------------------------------------
        # unknown parameter
        "_subjectName": subject_id if subject_id else "",
    }

    # Map the manufacturer name with Connectomist convention
    if manufacturer not in MANUFACTURERS:
        raise ConnectomistBadManufacturerNameError(manufacturer)
    parameters_dict["manufacturer"] = MANUFACTURERS[manufacturer]

    # Read .bval file to infer nb of T2, nb of shells...
    # Check the at least one bval is not null
    try:
        bvalues = np.loadtxt(bval)
        if set(bvalues) == {0}:  # If only 0s raise Exception
            raise Exception
    except:
        raise ConnectomistBadFileError(bval)
    nb_t2 = np.sum(bvalues == 0)  # nb of volumes where bvalue=0
    bvals_set = set(bvalues) - {0}    # set of non-zero bvalues
    nb_shells = len(bvals_set)

    # Update Connectomist step description
    parameters_dict["numberOfT2"] = nb_t2
    if nb_shells == 1:
        # Spherical single-shell custom
        parameters_dict["qSpaceSamplingType"] = 4
    else:
        raise ConnectomistError("Multiple shell models not handled. "
                                "Path to .bval file: %s" % bval)

    # Check validity of .bvec file.
    # If bvec file does not exist or filled with 0s, raise Error
    if (not os.path.isfile(bvec)) or np.loadtxt(bvec).max() == 0:
        raise ConnectomistBadFileError(bvec)

    # Call with Connectomist
    connprocess = ConnectomistWrapper(path_connectomist)
    parameter_file = ConnectomistWrapper.create_parameter_file(
        algorithm, parameters_dict, outdir)
    connprocess(algorithm, parameter_file, outdir, nb_tries=nb_tries)

    # When there are multiple t2 (nodif) volumes, Connectomist merges them
    # rewrite bvec, bval file accordingly (remove extra T2 values)
    if nb_t2 > 1:
        dw_indexes = np.where(bvalues != 0)[0]
        new_bvalues = np.concatenate(([0], bvalues[dw_indexes]))
        np.savetxt(bval, new_bvalues)

        bvecs = np.loadtxt(bvec)
        new_bvecs = np.concatenate(([[0], [0], [0]], bvecs[:, dw_indexes]),
                                   axis=1)
        np.savetxt(bvec, new_bvecs)

    return outdir
