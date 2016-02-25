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
import numpy as np
import nibabel

from clindmri.preproc.utils                    import select_first_b0
from clindmri.registration.utils               import extract_image
from clindmri.segmentation.fsl                 import bet2
from clindmri.tractography.fsl                 import bedpostx, probtrackx2
from clindmri.extensions.freesurfer.wrappers   import FSWrapper
from clindmri.extensions.freesurfer.exceptions import FreeSurferRuntimeError
from clindmri.extensions.configuration         import environment
from clindmri.plot.slicer                      import plot_image, animate_image

import matplotlib.pyplot as plt

from scipy.ndimage.morphology import binary_dilation

from clindmri.connectivity.freesurfer_color_lut import (APARC_ASEG_RAW_LUT,
                                                        APARC_A2009S_ASEG_RAW_LUT,
                                                        ASEG_RAW_LUT)

# To avoid using strings everywhere, associate a class attribute to each type
# of tracto mask
class TractoMaskTypes:

    nodif_brain            = "nodif_brain"
    wm                     = "wm"
    wm_dilated_1vox_6conn  = "wm_dilated_1vox_6conn"
    wm_dilated_1vox_14conn = "wm_dilated_1vox_14conn"

    choices = {nodif_brain, wm, wm_dilated_1vox_6conn, wm_dilated_1vox_14conn}


def register_diffusion_to_anatomy(outdir,
                                  nodif_brain,
                                  subject_id,
                                  fs_subjects_dir = None,
                                  subdir          = "diff_to_anat",
                                  fsl_init        = "/etc/fsl/5.0/fsl.sh"):
    """
    Register the diffusion to the anatomy (T1) using Freesurfer bbregister
    (boundary-based registration).
    <unit>
        <input name="outdir"          type="Directory"           />
        <input name="nodif_brain"     type="File"                />
        <input name="subject_id"      type="Str"                 />
        <input name="fs_subjects_dir" type="Directory"           />
        <input name="subdir"          type="Str"                 />
        <input name="fsl_init"        type="File"                />

        <output name="dif2anat_dat"   type="File" description="
            The .dat file created by tkregister/bbregister cmd." />
    </unit>
    """

    if fs_subjects_dir is None:
        if "SUBJECTS_DIR" in os.environ:
            fs_subjects_dir = os.environ["SUBJECTS_DIR"]
        else:
            raise ValueError("Missing <SUBJECTS_DIR>: set the $SUBJECTS_DIR "
                             "environment variable for Freesurfer or pass it "
                             "as an argument.")

    # If a subdirectory name has been passed, adapt outdir
    if subdir:
        outdir = os.path.join(outdir, subdir)

    # Create outdir if it does not exist
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    dif2anat_dat = os.path.join(outdir, "dif2anat.dat")
    cmd = ["bbregister", "--s", subject_id, "--mov", nodif_brain,
           "--reg", dif2anat_dat, "--dti", "--init-fsl"]
    fsprocess = FSWrapper(cmd)
    fsprocess.environment["SUBJECTS_DIR"] = fs_subjects_dir

    # bbregister requires FSL, so we have to add the FSL environment
    fsl_env = environment(fsl_init)
    for k, i in fsl_env.items():
        if k not in fsprocess.environment:
            fsprocess.environment[k] = i
        else:
            # A variable that exists in both FS and FSL environments
            if k == "PATH":
                fsprocess.environment["PATH"] += ":" + fsl_env["PATH"]
            else:
                pass  # ignore this variable

    fsprocess()  # Run
    if fsprocess.exitcode != 0:
        raise FreeSurferRuntimeError(cmd[0], " ".join(cmd[1:]))

    return dif2anat_dat


def qc_dif2anat_registration(outdir,
                             nodif_brain,
                             dif2anat_dat,
                             subject_id,
                             fs_subjects_dir = None,
                             subdir          = "qc"):
    """
    Function meant to help quality check (qc) the registration between
    the diffusion and the anatomy (T1 from Freesurfer recon-all).
    It creates snap shots to help with the visualization:
        - T1 brain registered in diffusion + contour of nodif brain volume.

    The snap shot is saved in <outdir>/<subdir>/"t1_to_diff.pdf". By default
    <subdir> is "qc". To write in outdir directly, set subdir to anything that
    evaluates to False (empty string or None).
    Directories are automatically created if they don't exist.

    <unit>
        <input name="outdir"          type="Directory" />
        <input name="nodif_brain"     type="File"      />
        <input name="dif2anat_dat"    type="File"      />
        <input name="subject_id"      type="Str"       />
        <input name="fs_subjects_dir" type="Directory" />
        <input name="subdir"          type="Str"       />

        <output name="qc_dir"         type="Directory" />
    </unit>
    """

    if fs_subjects_dir is None:
        if "SUBJECTS_DIR" in os.environ:
            fs_subjects_dir = os.environ["SUBJECTS_DIR"]
        else:
            raise ValueError("Missing <SUBJECTS_DIR>: set the $SUBJECTS_DIR "
                             "environment variable for Freesurfer or pass it "
                             "as an argument.")

    # If requested use a subdirectory in outdir
    if subdir:
        outdir = os.path.join(outdir, subdir)

    # If outdir does not exist, create it
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Paths
    t1        = os.path.join(fs_subjects_dir, subject_id, "mri/brainmask.mgz")
    t1_to_dif = os.path.join(outdir, "t1_to_dif.nii.gz")

    # Project T1 in dif
    cmd = ["mri_vol2vol", "--mov", nodif_brain, "--targ", t1,
           "--inv", "--interp", "nearest", "--o", t1_to_dif,
           "--reg", dif2anat_dat, "--no-save-reg"]
    fsprocess = FSWrapper(cmd)
    fsprocess()  # Run
    if fsprocess.exitcode != 0:
        raise FreeSurferRuntimeError(cmd[0], " ".join(cmd[1:]))

    nb_slices_in_z = nibabel.load(nodif_brain).get_shape()[2]

    # Gif of the T1 in diffusion space with diffusion edges
    animate_image(t1_to_dif, edge_file=nodif_brain, outdir=outdir,
                  name="t1_to_dif", cut_coords=nb_slices_in_z/2, clean=True)

    # First png: T1 registered in diffusion with nodif edges
    t1_with_nodif_edges_png = os.path.join(outdir, "t1_with_nodif_edges.png")
    plot_image(t1_to_dif,
               edge_file  = nodif_brain,
               snap_file  = t1_with_nodif_edges_png,
               name       = "T1 in diffusion + edges of nodif",
               cut_coords = nb_slices_in_z/2)

    # Second png: nodif with edges of T1 registered in diffusion
    nodif_with_t1_edges_png = os.path.join(outdir, "nodif_with_t1_edges.png")
    plot_image(nodif_brain,
               edge_file  = t1_to_dif,
               snap_file  = nodif_with_t1_edges_png,
               name       = "nodif + edges of registered T1",
               cut_coords = nb_slices_in_z/2)

#    # Third png: nodif with WM edges (from T1 segmentation)
#    nodif_with_wm_edges_png = os.path.join(outdir, "nodif_with_wm_edges.png")
#    plot_image(nodif_brain,
#               edge_file  = t1_to_dif,
#               snap_file  = nodif_with_t1_edges_png,
#               name       = "nodif + edges of registered T1",
#               cut_coords = nb_slices_in_z/2)

    # Return something for Capsul
    qc_dir = outdir

    return qc_dir


def dilate_mask_by_one_voxel(input_nifti, output_nifti=None):
    """
    Dilate Nifti image by one voxel, using a 6-neighborhood structuring element
    (voxels in diagonal are not included in the dilation).
    """

    if not input_nifti.endswith(".nii.gz"):
        raise ValueError("Input has to be .nii.gz file, passed: %s" % input_nifti)

    image = nibabel.load(input_nifti)
    dtype = image.get_data_dtype()

    if output_nifti is None:
        output_nifti = input_nifti.split(".nii.gz")[0] + "_dilated_1vox_6conn.nii.gz"

    structuring_element = np.array([[[0, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 0]],

                                    [[0, 1, 0],
                                     [1, 1, 1],
                                     [0, 1, 0]],

                                    [[0, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 0]]])
    data         = image.get_data()
    data_dilated = binary_dilation(data, structuring_element).astype(dtype)

    # Create and save Nifti
    image_dilated = nibabel.Nifti1Image(data_dilated, image.get_affine(),
                                        header=image.get_header())
    image_dilated.to_filename(output_nifti)

    return output_nifti


def create_masks_for_tractography(outdir,
                                  nodif_brain,
                                  dif2anat_dat,
                                  subject_id,
                                  cortical_atlas   = "Desikan",
                                  add_subcortical  = True,
                                  tracto_mask_type = "nodif_brain",
                                  fs_subjects_dir  = None):
    """
    Create the volume masks required for the probabilist tractography:
        - tractography mask (nodif_brain or white matter mask)
        - seed masks, to be used as seeds in probtrackx2, the seeds can include
          the subcortical regions (from Freesurfer 'aseg.mgz')

    Parameters
    ----------
    outdir:           Directory where to output.
    nodif_brain:      Brain-only volume extracted from the preprocessed DWI.
    dif2anat_dat:     The registration .dat file, diffusion to anatomy transformation.
    cortical_atlas:   Select which freesurfer cortical parcellation to use, either
                      "Desikan" (default) or "Destrieux".
    add_subcortical:  If True, create masks for the subcortical regions from "aseg.mgz".
    tracto_mask_type: The type of tractography mask to create, allowed types:
                      "wm", "wm_dilated_1vox_6conn" (default), "wm_dilated_1vox_14conn"
                      or "nodif_brain" (whole brain).
                      Two of the proposed white matter masks are dilated because a
                      non-dilated white matter mask does not overlap with the "gray"
                      subcortical regions, therefore the tracts will never get there.
                      Moreover the right and left white matter regions are much less
                      connected without dilation, therefore the connectogram shows
                      few interhemisphere connections with a simple white matter mask.
    fs_subjects_dir:  If the Freesurfer $SUBJECTS_DIR environment variable is
                      not set, or to bypass it, pass the path.

    <unit>

        <input name="outdir"           type="Directory"           />
        <input name="nodif_brain"      type="File"                />
        <input name="dif2anat_dat"     type="File"                />
        <input name="subject_id"       type="Str"                 />
        <input name="cortical_atlas"   type="Str"                 />
        <input name="add_subcortical"  type="Bool"                />
        <input name="tracto_mask_type" type="Str"                 />
        <input name="fs_subjects_dir"  type="Directory"           />

        <output name="seed_masks"      type="List" content="File" />
        <output name="tracto_mask"     type="File"                />
        <output name="stop_mask"       type="File"                />
        <output name="avoid_mask"      type="File"                />

    </unit>
    """

    # -------------------------------------------------------------------------
    # Check arguments
    if tracto_mask_type not in TractoMaskTypes.choices:
        raise ValueError("Bad argument 'tracto_mask_type': {}, should be in {}"
                         .format(tracto_mask_type, TractoMaskTypes.choices))

    if fs_subjects_dir is None:
        if "SUBJECTS_DIR" in os.environ:
            fs_subjects_dir = os.environ["SUBJECTS_DIR"]
        else:
            raise ValueError("Missing <SUBJECTS_DIR>: set the $SUBJECTS_DIR "
                             "environment variable for Freesurfer or pass it "
                             "as an argument.")

    # -------------------------------------------------------------------------
    # Set the paths according to the requested options

    subject_dir = os.path.join(fs_subjects_dir, subject_id)

    # Create a subdirectory for cortical masks in outdir, if not existing
    cortical_masks_dir = os.path.join(outdir, "cortical_masks_" + cortical_atlas)
    if not os.path.isdir(cortical_masks_dir):
        os.mkdir(cortical_masks_dir)

    # Create a subdirectory for other masks in outdir, if not existing
    # It is used for subcortical masks, 'stop' mask , 'avoid' mask...
    other_masks_dir = os.path.join(outdir, "other_masks")
    if not os.path.isdir(other_masks_dir):
        os.mkdir(other_masks_dir)

    # Set the right paths according to the chosen cortical atlas
    if cortical_atlas == "Desikan":
        cortical_seg     = os.path.join(subject_dir, "mri/aparc+aseg.mgz")
        cortical_seg2dif = os.path.join(cortical_masks_dir, "aparc+aseg2dif.nii.gz")
        raw_cortical_LUT = APARC_ASEG_RAW_LUT
    elif cortical_atlas == "Destrieux":
        cortical_seg     = os.path.join(subject_dir, "mri/aparc.a2009s+aseg.mgz")
        cortical_seg2dif = os.path.join(cortical_masks_dir, "aparc.a2009s+aseg2dif.mgz")
        raw_cortical_LUT = APARC_A2009S_ASEG_RAW_LUT
    else:
        raise ValueError("Bad 'cortical_atlas' name: {}".format(cortical_atlas))

    # -------------------------------------------------------------------------
    # Project cortical segmentation in diffusion
    cmd = ["mri_vol2vol", "--mov", nodif_brain, "--targ", cortical_seg,
           "--inv", "--interp", "nearest", "--o", cortical_seg2dif,
           "--reg", dif2anat_dat, "--no-save-reg"]
    fsprocess = FSWrapper(cmd)
    fsprocess()  # Run
    if fsprocess.exitcode != 0:
        raise FreeSurferRuntimeError(cmd[0], " ".join(cmd[1:]))

    # ------------------------------------------------------------------------
    # Project subcortical segmentation in diffusion

    # Set paths for subcortical segmentation, and the projection in diffusion
    subcortical_seg     = os.path.join(subject_dir, "mri/aseg.mgz")
    subcortical_seg2dif = os.path.join(other_masks_dir, "aseg2dif.nii.gz")

    cmd = ["mri_vol2vol", "--mov", nodif_brain, "--targ", subcortical_seg,
           "--inv", "--interp", "nearest", "--o", subcortical_seg2dif,
           "--reg", dif2anat_dat, "--no-save-reg"]
    fsprocess = FSWrapper(cmd)
    fsprocess()  # Run
    if fsprocess.exitcode != 0:
        raise FreeSurferRuntimeError(cmd[0], " ".join(cmd[1:]))

    # -------------------------------------------------------------------------
    # Create the tracto, according to the requested tracto mask type
    if tracto_mask_type == TractoMaskTypes.wm:
        tracto_mask = os.path.join(other_masks_dir, "wm_mask.nii.gz")
        cmd = ["mri_binarize", "--i", subcortical_seg2dif, "--wm",
               "--o", tracto_mask]

        fsprocess = FSWrapper(cmd)
        fsprocess()  # Run
        if fsprocess.exitcode != 0:
            raise FreeSurferRuntimeError(cmd[0], " ".join(cmd[1:]))

    elif tracto_mask_type == TractoMaskTypes.wm_dilated_1vox_6conn:
        wm_mask = os.path.join(other_masks_dir, "wm_mask.nii.gz")
        cmd = ["mri_binarize", "--i", subcortical_seg2dif, "--wm", "--o", wm_mask]

        fsprocess = FSWrapper(cmd)
        fsprocess()  # Run
        if fsprocess.exitcode != 0:
            raise FreeSurferRuntimeError(cmd[0], " ".join(cmd[1:]))

        tracto_mask = dilate_mask_by_one_voxel(wm_mask)

    elif tracto_mask_type == TractoMaskTypes.wm_dilated_1vox_14conn:
        tracto_mask = os.path.join(other_masks_dir, "wm_mask_dilated_1vox_14conn.nii.gz")
        cmd = ["mri_binarize", "--i", subcortical_seg2dif, "--wm",
               "--dilate", str(1), "--o", tracto_mask]

        fsprocess = FSWrapper(cmd)
        fsprocess()  # Run
        if fsprocess.exitcode != 0:
            raise FreeSurferRuntimeError(cmd[0], " ".join(cmd[1:]))

    elif tracto_mask_type == TractoMaskTypes.nodif_brain:
        tracto_mask = nodif_brain

    else:
        raise ValueError("'tracto_mask_type': {}, should be in {}"
                         .format(tracto_mask_type, TractoMaskTypes.choices))

    # -------------------------------------------------------------------------

    # List to accumulate all paths to the created masks
    # Creating multiple lists helps ordering the masks
    left_cortical_masks, right_cortical_masks = [], []
    left_subcortical_masks, right_subcortical_masks = [], []

    # -------------------------------------------------------------------------
    # Parse the cortical Look Up Table and create cortical mask volumes
    for line in raw_cortical_LUT:
        lh_ctx_int_label, lh_ctx_str_label, _, _, _, _ = line.strip().split()

        lh_str_label = lh_ctx_str_label[len("ctx-"):]  # Remove prefix ("ctx-" or "ctx_")
        lh_str_label = lh_str_label.replace("lh",   "left",  1)
        rh_str_label = lh_str_label.replace("left", "right", 1)

        # Infer right cortex int label from left label
        lh_ctx_int_label = int(lh_ctx_int_label)
        rh_ctx_int_label = lh_ctx_int_label + 1000

        # Infer white matter labels from cortical labels
        lh_wm_int_label = lh_ctx_int_label + 2000
        rh_wm_int_label = rh_ctx_int_label + 2000

        # Create left mask
        path_lh_mask = os.path.join(cortical_masks_dir, lh_str_label + ".nii.gz")
        cmd = ["mri_binarize", "--i", cortical_seg2dif,
               "--match", str(lh_ctx_int_label), str(lh_wm_int_label),
               "--o", path_lh_mask]
        fsprocess = FSWrapper(cmd)
        fsprocess()  # Run
        if fsprocess.exitcode != 0:
            raise FreeSurferRuntimeError(cmd[0], " ".join(cmd[1:]))
        left_cortical_masks.append(path_lh_mask)

        # Create right mask
        path_rh_mask = os.path.join(cortical_masks_dir, rh_str_label + ".nii.gz")
        cmd = ["mri_binarize", "--i", cortical_seg2dif,
               "--match", str(rh_ctx_int_label), str(rh_wm_int_label),
               "--o", path_rh_mask]
        fsprocess = FSWrapper(cmd)
        fsprocess()  # Run
        if fsprocess.exitcode != 0:
            raise FreeSurferRuntimeError(cmd[0], " ".join(cmd[1:]))
        right_cortical_masks.append(path_rh_mask)

    # -------------------------------------------------------------------------
    # Create "avoid" mask: mask of the ventricles

    avoid_mask = os.path.join(other_masks_dir, "ventricles.nii.gz")
    cmd = ["mri_binarize", "--i", subcortical_seg2dif, "--ventricles",
           "--o", avoid_mask]
    fsprocess = FSWrapper(cmd)
    fsprocess()  # Run
    if fsprocess.exitcode != 0:
        raise FreeSurferRuntimeError(cmd[0], " ".join(cmd[1:]))

    # -------------------------------------------------------------------------
    # Create tracto stop mask: cortex + ganglia
    stop_mask = os.path.join(other_masks_dir, "cortex_and_ganglia.nii.gz")
    stop_mask_labels = ["3", "42"]  # Initialization with cortex labels
    for line in ASEG_RAW_LUT:  # Add ganglio labels to the list
        label = line.strip().split()[0]
        stop_mask_labels.append(label)
    cmd = ["mri_binarize", "--i", subcortical_seg2dif,
           "--match", " ".join(stop_mask_labels),
           "--o", stop_mask]
    fsprocess = FSWrapper(cmd)
    fsprocess()  # Run
    if fsprocess.exitcode != 0:
        raise FreeSurferRuntimeError(cmd[0], " ".join(cmd[1:]))

    # -------------------------------------------------------------------------
    # If user has requested subcortical masks: create them
    if add_subcortical:

        # Parse the cortical Look Up Table and create mask volumes
        for line in ASEG_RAW_LUT:
            int_label, str_label, _, _, _, _ = line.strip().split()

            # Create left mask
            path_mask = os.path.join(other_masks_dir, str_label + ".nii.gz")
            cmd = ["mri_binarize", "--i", subcortical_seg2dif,
                   "--match", int_label, "--o", path_mask]
            fsprocess = FSWrapper(cmd)
            fsprocess()  # Run
            if fsprocess.exitcode != 0:
                raise FreeSurferRuntimeError(cmd[0], " ".join(cmd[1:]))
            if str_label.startswith("Left"):
                left_subcortical_masks.append(path_mask)
            else:
                right_subcortical_masks.append(path_mask)

    # -------------------------------------------------------------------------

    seed_masks  = (left_cortical_masks + left_subcortical_masks +
                   right_cortical_masks + right_subcortical_masks)

    return seed_masks, tracto_mask, stop_mask, avoid_mask


# To complete with more snap shots, in particular aseg2dif.nii.gz with colors
# + ventricles (avoid mask) + stop mask (cortex + ganglia)
def qc_tractography_masks(outdir,
                          nodif_brain,
                          tracto_mask,
                          seed_masks,
                          subject_id,
                          cortical_atlas  = "Desikan",
                          fs_subjects_dir = None,
                          subdir          = "qc"):
    """
    Function meant to help quality check (qc) the masks created by the
    create_masks_for_tractography function.
    It creates snap shots to visualize the quality of the registration
    of the tractography mask (white matter) and seed masks in the diffusion space.

    The snap shots are saved in <outdir>/<subdir>/<snap shot>. By default
    <subdir> is "qc". To write in outdir directly, set subdir to anything that
    evaluates to False (empty string or None).
    Directories are automatically created if they don't exist.

    <unit>
        <input name="outdir"          type="Directory"           />
        <input name="nodif_brain"     type="File"                />
        <input name="tracto_mask"     type="File"                />
        <input name="seed_masks"      type="List" content="File" />
        <input name="subject_id"      type="Str"                 />
        <input name="cortical_atlas"  type="Str"                 />
        <input name="fs_subjects_dir" type="Directory"           />
        <input name="subdir"          type="Str"                 />

        <output name="qc_dir"         type="Directory"           />
    </unit>
    """

    if fs_subjects_dir is None:
        if "SUBJECTS_DIR" in os.environ:
            fs_subjects_dir = os.environ["SUBJECTS_DIR"]
        else:
            raise ValueError("Missing <SUBJECTS_DIR>: set the $SUBJECTS_DIR "
                             "environment variable for Freesurfer or pass it "
                             "as an argument.")

    # If requested use a subdirectory in outdir
    if subdir:
        outdir = os.path.join(outdir, subdir)

    # If outdir does not exist, create it
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    nb_slices_in_z = nibabel.load(nodif_brain).get_shape()[2]

    # Snap shot for registration and white matter mask checking
    wm_mask_to_dif_png = os.path.join(outdir, "tracto_mask.png")
    plot_image(nodif_brain,
               overlay_file = tracto_mask,
               snap_file    = wm_mask_to_dif_png,
               name         = "White matter mask in diffusion",
               cut_coords   = nb_slices_in_z/2)

    # Gif of white matter mask in diffusion
    animate_image(nodif_brain, overlay_file=tracto_mask, outdir=outdir,
                  name="tracto_mask", cut_coords=nb_slices_in_z/2,
                  clean=True)

    # Return something for Capsul
    qc_dir = outdir
    return qc_dir


def extract_nodif_volume(outdir, dwi, bval):
    """

    Parameters
    ----------
    outdir: Str, path to directory where to write "b0.nii.gz"
    dwi:    Str, path to DW data in which at least one volume was acquired
            with bvalue=0.
    bval:   Str, path to .bval file associated to the DW data.

    Return
    ------
    nodif_volume: Str, path to a/the volume for which bvalue is 0.
    <unit>
        <input name="outdir" type="Directory" />
        <input name="dwi"    type="File"      />
        <input name="bval"   type="File"      />

        <output name="nodif_volume" type="File" />
    </unit>
    """
    # Get index of the first volume acquired with bvalue=0
    b0_index = select_first_b0(bval)

    # Extract volume to a temporary Nifti
    nodif_volume = os.path.join(outdir, "nodif.nii.gz")
    extract_image(dwi, index=b0_index, out_file=nodif_volume)

    return nodif_volume


def bet2_nodif_brain(outdir, dwi, bval, subdir="bet2_nodif_brain", qc=True):
    """
    Extract brain from b0 volume, i.e. in DW data for a volume where bvalue=0.

    Parameters
    ----------
    outdir: Str, path to directory where to output.
    dwi:    Str, path to DW data in which at least one volume was acquired
            with bvalue=0.
    bval:   Str, path to .bval file associated to the DW data.
    subdir: Str, if you want the result files to be written in a subdirectory,
            specify the name, by default "bet2_nodif_brain".

    Return
    ------
    nodif_brain:      Str, path to the brain only volume.
    nodif_brain_mask: Str, path to the brain-only binary mask.

    <unit>
        <output name="nodif_brain"      type="File"      />
        <output name="nodif_brain_mask" type="File"      />

        <input name="outdir"            type="Directory" />
        <input name="dwi"               type="File"      />
        <input name="bval"              type="File"      />
        <input name="subdir"            type="Str"       />
        <input name="qc"                type="Bool"      />
    </unit>
    """
    if subdir:
        outdir = os.path.join(outdir, subdir)

    # Create outdir if it does not exist
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Get a/the volume with bvalue=0.
    nodif_volume = extract_nodif_volume(outdir, dwi, bval)

    # Set output path with desired prefix name
    output_prefix = os.path.join(outdir, "nodif_brain")

    # Run FSL bet2
    bet2(nodif_volume, output_prefix, f=0.25, m=True)

    # Output paths
    nodif_brain      = output_prefix + ".nii.gz"
    nodif_brain_mask = output_prefix + "_mask.nii.gz"

    # Check existence of resulting files
    if not (os.path.isfile(nodif_brain) and os.path.isfile(nodif_brain_mask)):
        raise Exception("bet2: output file(s) missing: %s or/and %s ."
                        % (nodif_brain, nodif_brain_mask))

    # If Quality Check, generate a PNG snapshot
    if qc:
        # Snap shot of brain-only contour on T2 image
        brain_contour_png = os.path.join(outdir, "nodif_brain_mask.png")
        plot_image(nodif_volume, contour_file=nodif_brain_mask,
                   snap_file=brain_contour_png, name="nodif_brain_mask")

    return nodif_brain, nodif_brain_mask


def create_bedpostx_indir(bedpostx_indir, dwi, bval, bvec, nodif_brain_mask):
    """
    Bedpostx requires an input directory with specific filenames inside to be run.
    The function moves the required files to a directory with the right names.
    """
    # Create directory if it does not exist
    if not os.path.isdir(bedpostx_indir):
        os.makedirs(bedpostx_indir)

    # Extension is either ".nii" or ".nii.gz"
    data_ext = ".nii" if dwi.endswith(".nii") else ".nii.gz"
    mask_ext = ".nii" if nodif_brain_mask.endswith(".nii") else ".nii.gz"

    shutil.copy2(dwi,  os.path.join(bedpostx_indir, "data" + data_ext))
    shutil.copy2(bval, os.path.join(bedpostx_indir, "bvals"))
    shutil.copy2(bvec, os.path.join(bedpostx_indir, "bvecs"))
    shutil.copy2(nodif_brain_mask,
                 os.path.join(bedpostx_indir, "nodif_brain_mask" + mask_ext))

    return bedpostx_indir


def bedpostx_diffusion_model(outdir, dwi, bval, bvec, nodif_brain_mask):
    """
    <unit>
        <output name="bedpostx_dir"    type="Directory" />

        <input name="outdir"           type="Directory" />
        <input name="dwi"              type="File"      />
        <input name="bval"             type="File"      />
        <input name="bvec"             type="File"      />
        <input name="nodif_brain_mask" type="File"      />
    </unit>
    """
    bedpostx_indir = os.path.join(outdir, "bedpostx")

    create_bedpostx_indir(bedpostx_indir, dwi, bval, bvec, nodif_brain_mask)
    bedpostx_dir = bedpostx(bedpostx_indir)[0]

    return bedpostx_dir


def probtrackx2_connectogram(outdir,
                             bedpostx_dir,
                             seed_masks,
                             tracto_mask,
                             stop_mask  = None,
                             avoid_mask = None,
                             subdir     = "probtrackx2",
                             nsamples   = 5000,
                             nsteps     = 2000,
                             #cthr=0.2,
                             #loopcheck=None,
                             #onewaycondition=None,
                             #usef=None,
                             #simple=None,
                             #seedref=None,
                             #steplength=0.5,
                             #fibthresh=0.01,
                             #distthresh=0.0,
                             #sampvox=0.0,
                             #network=None
                             ):
    """
    <unit>
        <output name="proba_matrix"   type="File"                />
        <output name="network_matrix" type="File"                />

        <input name="outdir"          type="Directory"           />
        <input name="bedpostx_dir"    type="Directory"           />
        <input name="seed_masks"      type="List" content="File" />
        <input name="tracto_mask"     type="File"                />
        <input name="stop_mask"       type="File"                />
        <input name="avoid_mask"      type="File"                />
        <input name="subdir"          type="Str"                 />
        <input name="nsamples"        type="Int"                 />
        <input name="nsteps"          type="Int"                 />
    </unit>
    """

    if subdir:
        outdir = os.path.join(outdir, subdir)

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    file_listing_masks = os.path.join(outdir, "masks.txt")
    np.savetxt(file_listing_masks, np.asarray(seed_masks), delimiter=" ", fmt="%s")

    basepath_samples = os.path.join(bedpostx_dir, "merged")

    proba_files, network_matrix = probtrackx2(network         = True,
                                              seed            = file_listing_masks,
                                              loopcheck       = True,
                                              onewaycondition = True,
                                              samples         = basepath_samples,
                                              mask            = tracto_mask,
                                              dir             = outdir,
                                              nsamples        = nsamples,
                                              nsteps          = nsteps)
    proba_matrix = proba_files[0]

    return proba_matrix, network_matrix


def plot_connectogram(outdir, network_matrix, seed_masks, transform=None):
    """
    Inspired from the plot_matrix() function, rewritten to make plots for the
    'small' connectogram, with the cortex region labels..

    Parameters
    ----------
    outdir:         Str, path to directory where to write 'connectogram.pdf'.
    network_matrix: Str, path.
    seed_masks:     List of str, paths to the seeds (masks).
    transform:      Callable (optional, default None), a callable function
                    to apply on the matrix.

    Returns
    -------
    snap_file: Str, path to the output snap shot: 'connectogram.png'.
    """

    # Check existence of the input files
    required_files = [network_matrix] + seed_masks
    for path in required_files:
        if not os.path.isfile(path):
            raise ValueError("File does not exist: %s" % path)

    # Load track counts
    matrix = np.loadtxt(network_matrix)

    # Apply transformation if requested
    if transform is not None:
        matrix = transform(matrix)

    # Extract names of seeds
    seed_names = [os.path.basename(l)[:-len(".nii.gz")] for l in seed_masks]

    ##################
    # Create the figure with matplotlib
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(matrix, cmap=plt.cm.Reds)

    ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticks(np.arange(matrix.shape[0])+0.5, minor=False)
    ax.set_yticks(np.arange(matrix.shape[1])+0.5, minor=False)
    ax.set_xticklabels(seed_names, minor=False, size=4)
    ax.set_yticklabels(seed_names, minor=False, size=4)

    ax.set_xticklabels(seed_names, minor=False)
    ax.set_yticklabels(seed_names, minor=False)

    plt.xticks(rotation=90)

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    ax.set_aspect("equal")
    plt.gcf().subplots_adjust(top=0.8)

    colorbar = fig.colorbar(heatmap)
    colorbar.set_label("Log(# of tracks)", rotation=270, labelpad=20)

    fig.tight_layout()
    #################

    # Path to the output PNG file
    snap_file = os.path.join(outdir, "connectogram.png")
    fig.savefig(snap_file, dpi=200)

    return snap_file


def qc_connectogram(outdir,
                    tracto_mask,
                    proba_matrix,
                    network_matrix,
                    seed_masks,
                    subdir="qc"):
    """
    Function meant to help quality check (qc) the tractography and the
    connectogram computed by the probtrackx2_probabilist_tractography function.
    It creates snap shots to visualize the connectogram and the fiber density
    in the diffusion space.

    The snap shots are saved in <outdir>/<subdir>/<snap shot>. By default
    <subdir> is "qc". To write in outdir directly, set subdir to anything that
    evaluates to False (empty string or None).
    Directories are automatically created if they don't exist.

    <unit>
        <input name="outdir"         type="Directory"           />
        <input name="tracto_mask"    type="File"                />
        <input name="proba_matrix"   type="File"                />
        <input name="network_matrix" type="File"                />
        <input name="seed_masks"     type="List" content="File" />
        <input name="subdir"         type="Str"                 />

        <output name="qc_dir"        type="Directory"           />
    </unit>
    """
    # If requested use a subdirectory in outdir
    if subdir:
        outdir = os.path.join(outdir, subdir)

    # If outdir does not exist, create it
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    nb_slices_in_z = nibabel.load(tracto_mask).get_shape()[2]

    fiber_density_png = os.path.join(outdir, "fiber_density_map.png")
    plot_image(tracto_mask,
               overlay_file = proba_matrix,
               snap_file    = fiber_density_png,
               name         = "fiber density map",
               overlay_cmap = "cold_hot",
               cut_coords   = nb_slices_in_z/2)

    # Second connectogram snap shot as PNG with cortex region labels
    plot_connectogram(outdir, network_matrix, seed_masks, transform=np.log1p)

    # Return something for Capsul
    qc_dir = outdir
    return qc_dir
