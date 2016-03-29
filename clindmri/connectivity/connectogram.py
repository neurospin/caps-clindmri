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
import tempfile
import re
import json

from clindmri.preproc.utils                    import select_first_b0
from clindmri.registration.utils               import extract_image
from clindmri.segmentation.fsl                 import bet2
from clindmri.tractography.fsl                 import bedpostx, probtrackx2
from clindmri.extensions.freesurfer.wrappers   import FSWrapper
from clindmri.extensions.freesurfer.exceptions import FreeSurferRuntimeError
from clindmri.extensions.fsl.wrappers          import FSLWrapper
from clindmri.extensions.fsl.exceptions        import FSLRuntimeError
from clindmri.extensions.configuration         import environment
from clindmri.plot.slicer                      import plot_image, animate_image

import matplotlib.pyplot as plt

from scipy.ndimage.morphology import binary_dilation


###############################################################################
# MODULE GLOBAL VARIABLES

# Get Freesurfer Look Up Table path
if "FREESURFER_HOME" in os.environ:
    PATH_LUT = os.path.join(os.environ["FREESURFER_HOME"], "FreeSurferColorLUT.txt")
else:
    raise Exception("Environment variable 'FREESURFER_HOME' is not set.")

# Load Freesurfer Look Up Table and create a dict mapping: <ROI name> -> label
try:
    LABEL_OF_ROI = dict(np.loadtxt(PATH_LUT, dtype=str, usecols=[1, 0]))
except:
    raise Exception("Failed to load Freesurfer Look Up Table: {}".format(PATH_LUT))


# Desikan atlas left regions without corpus callosum ordered like Lausanne 2008 atlas
DESIKAN_LH_ROIS = [
    'ctx-lh-lateralorbitofrontal',
    'ctx-lh-parsorbitalis',
    'ctx-lh-frontalpole',
    'ctx-lh-medialorbitofrontal',
    'ctx-lh-parstriangularis',
    'ctx-lh-parsopercularis',
    'ctx-lh-rostralmiddlefrontal',
    'ctx-lh-superiorfrontal',
    'ctx-lh-caudalmiddlefrontal',
    'ctx-lh-precentral',
    'ctx-lh-paracentral',
    'ctx-lh-rostralanteriorcingulate',
    'ctx-lh-caudalanteriorcingulate',
    'ctx-lh-posteriorcingulate',
    'ctx-lh-isthmuscingulate',
    'ctx-lh-postcentral',
    'ctx-lh-supramarginal',
    'ctx-lh-superiorparietal',
    'ctx-lh-inferiorparietal',
    'ctx-lh-precuneus',
    'ctx-lh-cuneus',
    'ctx-lh-pericalcarine',
    'ctx-lh-lateraloccipital',
    'ctx-lh-lingual',
    'ctx-lh-fusiform',
    'ctx-lh-parahippocampal',
    'ctx-lh-entorhinal',
    'ctx-lh-temporalpole',
    'ctx-lh-inferiortemporal',
    'ctx-lh-middletemporal',
    'ctx-lh-bankssts',
    'ctx-lh-superiortemporal',
    'ctx-lh-transversetemporal',
    'ctx-lh-insula'
]

# Desikan atlas left regions without corpus callosum ordered like Lausanne 2008 atlas
DESIKAN_RH_ROIS = [
    'ctx-rh-lateralorbitofrontal',
    'ctx-rh-parsorbitalis',
    'ctx-rh-frontalpole',
    'ctx-rh-medialorbitofrontal',
    'ctx-rh-parstriangularis',
    'ctx-rh-parsopercularis',
    'ctx-rh-rostralmiddlefrontal',
    'ctx-rh-superiorfrontal',
    'ctx-rh-caudalmiddlefrontal',
    'ctx-rh-precentral',
    'ctx-rh-paracentral',
    'ctx-rh-rostralanteriorcingulate',
    'ctx-rh-caudalanteriorcingulate',
    'ctx-rh-posteriorcingulate',
    'ctx-rh-isthmuscingulate',
    'ctx-rh-postcentral',
    'ctx-rh-supramarginal',
    'ctx-rh-superiorparietal',
    'ctx-rh-inferiorparietal',
    'ctx-rh-precuneus',
    'ctx-rh-cuneus',
    'ctx-rh-pericalcarine',
    'ctx-rh-lateraloccipital',
    'ctx-rh-lingual',
    'ctx-rh-fusiform',
    'ctx-rh-parahippocampal',
    'ctx-rh-entorhinal',
    'ctx-rh-temporalpole',
    'ctx-rh-inferiortemporal',
    'ctx-rh-middletemporal',
    'ctx-rh-bankssts',
    'ctx-rh-superiortemporal',
    'ctx-rh-transversetemporal',
    'ctx-rh-insula'
]

# Ordered left subcortical regions of Lausanne 2008 scale 33 atlas
LH_SUBCORTICAL_ROIS = [
    'Left-Thalamus-Proper',
    'Left-Caudate',
    'Left-Putamen',
    'Left-Pallidum',
    'Left-Accumbens-area',
    'Left-Hippocampus',
    'Left-Amygdala'
]

# Ordered right subcortical regions of Lausanne 2008 scale 33 atlas
RH_SUBCORTICAL_ROIS = [
    'Right-Thalamus-Proper',
    'Right-Caudate',
    'Right-Putamen',
    'Right-Pallidum',
    'Right-Accumbens-area',
    'Right-Hippocampus',
    'Right-Amygdala',
]

# Non-hemispheric subcortical region of Lausanne 2008 scale 33 atlas
AXIAL_SUBCORTICAL_ROIS = ['Brain-Stem']


# ATLASES

# Lausanne2008 scale33 atlas is the combination of the Desikan cortical atlas,
# without the corpus callosum, with selected subcortical regions.
# The list is ordered
LAUSANNE2008_SCALE33_ROIS = (DESIKAN_LH_ROIS + LH_SUBCORTICAL_ROIS +
                             DESIKAN_RH_ROIS + RH_SUBCORTICAL_ROIS +
                             AXIAL_SUBCORTICAL_ROIS)

# Destrieux cortical atlas list of labels in the Freesurfer LUT
DESTRIEUX_CTX_LH_MIN_LABEL, DESTRIEUX_CTX_LH_MAX_LABEL = 11101, 11176
DESTRIEUX_LH_LABELS = [str(x)      for x in range(DESTRIEUX_CTX_LH_MIN_LABEL,
                                                  DESTRIEUX_CTX_LH_MAX_LABEL)]
DESTRIEUX_RH_LABELS = [str(x+1000) for x in range(DESTRIEUX_CTX_LH_MIN_LABEL,
                                                  DESTRIEUX_CTX_LH_MAX_LABEL)]

DESTRIEUX_LH_ROIS = [x for x in LABEL_OF_ROI if LABEL_OF_ROI[x] in set(DESTRIEUX_LH_LABELS)]
DESTRIEUX_RH_ROIS = [x for x in LABEL_OF_ROI if LABEL_OF_ROI[x] in set(DESTRIEUX_RH_LABELS)]

DESTRIEUX_WITH_SUBCORTICAL_ROIS = (DESTRIEUX_LH_ROIS + LH_SUBCORTICAL_ROIS +
                                   DESTRIEUX_RH_ROIS + RH_SUBCORTICAL_ROIS +
                                   AXIAL_SUBCORTICAL_ROIS)

# Set of available atlas
CORTICAL_ATLASES = frozenset(["Desikan", "Destrieux"])

# Set of available types of stop masks when using --target3 tractography
STOP_MASK_TYPES = {"target_rois", "inverse_wm"}

# Set of available types of tracto masks when using --network tractography
TRACTO_MASK_TYPES = frozenset(["nodif_brain", "wm", "wm_dilated_1vox_6conn",
                               "wm_dilated_1vox_14conn"])

###############################################################################
# Utility functions

def get_or_check_freesurfer_subjects_dir(subjects_dir=None):
    """
    Factorize the code to check the Freesurfer 'subjects_dir' or if not passed,
    to read the $SUBJECTS_DIR environment variable.
    If not passed and not in environment, raise an Exception.
    """
    if subjects_dir is not None:
        if not os.path.isdir(subjects_dir):
            raise ValueError("Non existing directory: {}".format(subjects_dir))
    elif "SUBJECTS_DIR" in os.environ:
        subjects_dir = os.environ["SUBJECTS_DIR"]
    else:
        raise ValueError("Missing 'subjects_dir': set the $SUBJECTS_DIR "
                         "environment variable for Freesurfer or pass it "
                         "as an argument.")

    return subjects_dir


def run_freesurfer_cmd(cmd, subjects_dir=None, add_fsl_env=False,
                       fsl_init="/etc/fsl/5.0/fsl.sh"):
    """
    To avoid repeating the code to run Freesurfer and check exitcode everywhere.
    Step:
        - add $SUBJECTS_DIR to the env if requested
        - add FSL environment if requested (some Freesurfer commands require FSL)
        - run the Freesurfer cmd
        - check exit code

    Parameters
    ----------
    cmd: list of str
        the command to run (subprocess like).
    subjects_dir: str, default None.
        To set the $SUBJECTS_DIR environment variable.
    add_fsl_env:  bool, default False
        To activate the FSL environment, required for commands like bbregister.
    fsl_init: str
        Path to the Bash script setting the FSL environment, if needed.
    """
    fsprocess = FSWrapper(cmd)

    if subjects_dir is not None:
        fsprocess.environment["SUBJECTS_DIR"] = subjects_dir

    if add_fsl_env:
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

    return fsprocess


def run_fsl_cmd(cmd):
    fslprocess = FSLWrapper()
    fslprocess()
    if fslprocess.exitcode != 0:
        raise FSLRuntimeError(fslprocess.cmd[0], " ".join(fslprocess.cmd[1:]),
                              fslprocess.stderr)

###############################################################################

def register_diffusion_to_anatomy(outdir,
                                  nodif_brain,
                                  subject_id,
                                  subjects_dir = None,
                                  subdir       = "diff_to_anat",
                                  fsl_init     = "/etc/fsl/5.0/fsl.sh"):
    """
    Register the diffusion to the anatomy (T1) using Freesurfer bbregister
    (boundary-based registration).

    The resulting .dat file is saved in <outdir>/<subdir>/dif2anat.dat
    To write in 'outdir' directly, set 'subdir' to anything that evaluates
    to False (empty string or None).

    'subjects_dir' has to be passed if not set as an environnement variable ($SUBJECTS_DIR).

    <unit>
        <input name="outdir"          type="Directory"           />
        <input name="nodif_brain"     type="File"                />
        <input name="subject_id"      type="Str"                 />
        <input name="subjects_dir"    type="Directory"           />
        <input name="subdir"          type="Str"                 />
        <input name="fsl_init"        type="File"                />

        <output name="dif2anat_dat"   type="File" description="
            The .dat file created by tkregister/bbregister." />
    </unit>
    """

    # Freesurfer 'subjects_dir' has to be passed or set as environment variable
    subjects_dir = get_or_check_freesurfer_subjects_dir(subjects_dir)

    # If a subdirectory name has been passed, adapt outdir
    if subdir:
        outdir = os.path.join(outdir, subdir)

    # Create outdir if it does not exist
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    dif2anat_dat = os.path.join(outdir, "dif2anat.dat")
    cmd = ["bbregister", "--s", subject_id, "--mov", nodif_brain,
           "--reg", dif2anat_dat, "--dti", "--init-fsl"]
    run_freesurfer_cmd(cmd, subjects_dir=subjects_dir, add_fsl_env=True)

    return dif2anat_dat


def qc_dif2anat_registration(outdir,
                             nodif_brain,
                             dif2anat_dat,
                             subject_id,
                             subjects_dir = None,
                             subdir       = "qc"):
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
        <input name="outdir"       type="Directory" />
        <input name="nodif_brain"  type="File"      />
        <input name="dif2anat_dat" type="File"      />
        <input name="subject_id"   type="Str"       />
        <input name="subjects_dir" type="Directory" />
        <input name="subdir"       type="Str"       />

        <output name="qc_dir"      type="Directory" />
    </unit>
    """

    # Freesurfer 'subjects_dir' has to be passed or set as environment variable
    subjects_dir = get_or_check_freesurfer_subjects_dir(subjects_dir)

    # If requested use a subdirectory in outdir
    if subdir:
        outdir = os.path.join(outdir, subdir)

    # If outdir does not exist, create it
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Paths
    t1        = os.path.join(subjects_dir, subject_id, "mri/brainmask.mgz")
    t1_to_dif = os.path.join(outdir, "t1_to_dif.nii.gz")

    # Project T1 in dif
    cmd = ["mri_vol2vol", "--mov", nodif_brain, "--targ", t1,
           "--inv", "--interp", "nearest", "--o", t1_to_dif,
           "--reg", dif2anat_dat, "--no-save-reg"]
    run_freesurfer_cmd(cmd)

    nb_slices_in_z = nibabel.load(nodif_brain).get_shape()[2]

    # GIF: registered T1 (in diffusion) + nodif_brain edges
    animate_image(t1_to_dif, edge_file=nodif_brain, outdir=outdir,
                  name="t1_with_nodif_edges.gif", cut_coords=nb_slices_in_z-2, clean=True)

    # GIF: nodif_brain + T1 edges
    animate_image(nodif_brain, edge_file=t1_to_dif, outdir=outdir,
                  name="nodif_with_t1_edges.gif", cut_coords=nb_slices_in_z-2, clean=True)

    # PNG: registered T1 (in diffusion) + nodif edges
    t1_with_nodif_edges_png = os.path.join(outdir, "t1_with_nodif_edges.png")
    plot_image(t1_to_dif,
               edge_file  = nodif_brain,
               snap_file  = t1_with_nodif_edges_png,
               name       = "T1 in diffusion + edges of nodif",
               cut_coords = nb_slices_in_z-2)

    # PNG: nodif brain + registered T1 edges
    nodif_with_t1_edges_png = os.path.join(outdir, "nodif_with_t1_edges.png")
    plot_image(nodif_brain,
               edge_file  = t1_to_dif,
               snap_file  = nodif_with_t1_edges_png,
               name       = "nodif + edges of registered T1",
               cut_coords = nb_slices_in_z-2)

    # Return something for Capsul
    qc_dir = outdir

    return qc_dir


def dilate_mask_by_one_voxel(input_nifti, connexity=6, output_nifti=None):
    """
    Dilate Nifti image by one voxel, using a 6 or 14-neighborhood structuring
    element (with 6-connexity voxels in diagonal are not included in the dilation).
    """

    if connexity not in {6, 14}:
        raise ValueError("Bad argument 'connexity': has to be 6 or 14.")

    if not input_nifti.endswith(".nii.gz"):
        raise ValueError("Input has to be .nii.gz file, passed: %s" % input_nifti)

    image = nibabel.load(input_nifti)
    dtype = image.get_data_dtype()

    if output_nifti is None:
        suffix       = "_dilated_1vox_%iconn.nii.gz" % connexity
        output_nifti = input_nifti.split(".nii.gz")[0] + suffix

    if connexity == 6:
        structuring_element = np.array([[[0, 0, 0],
                                         [0, 1, 0],
                                         [0, 0, 0]],

                                        [[0, 1, 0],
                                         [1, 1, 1],
                                         [0, 1, 0]],

                                        [[0, 0, 0],
                                         [0, 1, 0],
                                         [0, 0, 0]]])
    else:  # connexity == 14
        structuring_element = np.ones((3,3,3))

    data         = image.get_data()
    data_dilated = binary_dilation(data, structuring_element).astype(dtype)

    # Create and save Nifti
    image_dilated = nibabel.Nifti1Image(data_dilated, image.get_affine(),
                                        header=image.get_header())
    image_dilated.to_filename(output_nifti)

    return output_nifti


def project_aparc_and_aseg_to_diffusion(outdir,
                                        dif2anat_dat,
                                        nodif_brain,
                                        subject_id,
                                        subjects_dir   = None,
                                        cortical_atlas = "Desikan",
                                        outext         = ".nii.gz"):
    """
    Apply the transform specified by dif2anat_dat to the cortical and
    subcortical Freesurfer segmentations (aparc and aseg).

    Parameters
    ----------
    outdir: str
        Directory where to output the 2 projections.
    dif2anat_dat: str
        Path to .dat file generated by the registration process (bbregister, tkregister...).
        Created by register_diffusion_to_anatomy().
    nodif_brain: str
        Path to the brain-only volume in diffusion space used when the registration was done.
    subject_id: str
        Id of the subject in the Freesurfer subjects_dir.
    subjects_dir: str, optional default None
        If the Freesurfer $SUBJECTS_DIR environment variable is not set, or to
        bypass it, pass the path.
    cortical_atlas: str, optional default 'Desikan'
        Name of the Freesurfer cortical atlas, either 'Desikan' or 'Destrieux'.
    outext: str, optional default ".nii.gz"
        By setting the 'output extension' you can specify the format.
    """

    # Freesurfer 'subjects_dir' has to be passed or set as environment variable
    subjects_dir = get_or_check_freesurfer_subjects_dir(subjects_dir)

    subj_dir = os.path.join(subjects_dir, subject_id)

    if not cortical_atlas in CORTICAL_ATLASES:
        raise ValueError("Bad 'cortical_atlas': {}, should be in {}"
                         .format(cortical_atlas, CORTICAL_ATLASES))

    # Set paths according to given atlas
    if cortical_atlas == "Desikan":
        aparc_mgz = os.path.join(subj_dir, "mri/aparc+aseg.mgz")
        aparc2dif = os.path.join(outdir, "aparc+aseg2dif%s" % outext)
    else:  # Destrieux atlas
        aparc_mgz = os.path.join(subj_dir, "mri/aparc.a2009s+aseg.mgz")
        aparc2dif = os.path.join(outdir, "aparc.a2009s+aseg2dif%s" % outext)

    # Project cortical segmentation in diffusion
    cmd = ["mri_vol2vol", "--mov", nodif_brain, "--targ", aparc_mgz,
           "--inv", "--interp", "nearest", "--o", aparc2dif,
           "--reg", dif2anat_dat, "--no-save-reg"]
    run_freesurfer_cmd(cmd)

    # Project subcortical segmentation in diffusion
    aseg_mgz = os.path.join(subj_dir, "mri/aseg.mgz")
    aseg2dif = os.path.join(outdir, "aseg2dif%s" % outext)
    cmd = ["mri_vol2vol", "--mov", nodif_brain, "--targ", aseg_mgz,
           "--inv", "--interp", "nearest", "--o", aseg2dif,
           "--reg", dif2anat_dat, "--no-save-reg"]
    run_freesurfer_cmd(cmd)

    return aparc2dif, aseg2dif


def create_white_matter_mask(outdir, path_aseg, outext=".nii.gz"):
    """
    Create a mask of the white matter (both hemispheres), from the
    Freesurfer subcortical segmentation (aseg).
    If you want the mask to be in a specific space, project the Freesurfer
    segmentation before calling this function.
    """
    wm_mask = os.path.join(outdir, "wm_mask%s" % outext)
    cmd = ["mri_binarize", "--i", path_aseg, "--wm", "--o", wm_mask]
    run_freesurfer_cmd(cmd)

    return wm_mask


def create_ventricles_mask(outdir, path_aseg, outext=".nii.gz"):
    """
    Create a mask of the ventricles (both hemispheres), from the Freesurfer
    subcortical segmentation (aseg).
    It can be used as 'avoid' mask in tractography.
    If you want the mask to be in a specific space, project the Freesurfer
    segmentation before calling this function.
    """
    ventricles_mask = os.path.join(outdir, "ventricles.nii.gz")
    cmd = ["mri_binarize", "--i", path_aseg, "--ventricles", "--o", ventricles_mask]
    run_freesurfer_cmd(cmd)

    return ventricles_mask


def create_target_masks(outdir, target_rois, path_aparc, outext=".nii.gz"):
    """
    Create a mask for each target region.
    If you want the masks to be in a specific space, project the segmentation
    (aparc) before calling this function.

    Parameters
    ----------
    target_rois: list of str
        The names of regions for which a mask is needed. The names should be
        the ones used in the FreesurferColorLUT.
    path_aparc: str
        The path to the 'aparc+aseg' segmentation of Freesurfer.

    Returns
    -------
    roi_masks: list of str
        The paths to the created masks.
    """
    # Create a mask for each target ROI
    roi_masks = []
    for roi in target_rois:
        mask_path = os.path.join(outdir, "%s.nii.gz" % roi)
        cmd = ["mri_binarize", "--i", path_aparc, "--match", LABEL_OF_ROI[roi],
               "--o", mask_path]
        run_freesurfer_cmd(cmd)
        roi_masks.append(mask_path)

    return roi_masks


def create_masks_for_tracto_seeding_endpoints(outdir,
                                              nodif_brain,
                                              nodif_brain_mask,
                                              dif2anat_dat,
                                              subject_id,
                                              cortical_atlas   = "Desikan",
                                              tracto_mask_type = "nodif_brain",
                                              subjects_dir     = None,
                                              subdir           = "masks"):
    """
    Create the volume masks required for the probabilist tractography
    when using the --omatrix1 option in probtrackx2 (seeding in ROIs):
        - tractography mask (nodif_brain or white matter mask (dilated or not))
        - ROI masks, depends on the cortical_atlas, can include the
          subcortical regions (from Freesurfer 'aseg.mgz')
        - stop mask: gray matter (cortex + ganglia)
        - avoid_mask: mask of ventricles.

    Note that is you want subcortical ROIs and want to use white matter as
    tractography mask, you have to select a dilated one otherwise there will be
    no overlap with some ROIs and thus the samples will never end up there.

    Parameters
    ----------
    outdir: str
        Path to the directory where to output.
    nodif_brain: str
        Path to the brain-only volume extracted from the preprocessed DWI.
    nodif_brain_mask: str
        Path to the brain-only binary mask.
    dif2anat_dat: str
        Path to the registration .dat file, diffusion to anatomy transformation.
    subject_id: str
        Id of the subject in the Freesurfer subjects_dir.
    cortical_atlas: str, default "Desikan"
        Name of the freesurfer cortical parcellation to use, either "Desikan"
        (default) or "Destrieux". The corpus callosum is not included.
    tracto_mask_type: str, default "nodif_brain
        "The type of tractography mask to create, allowed types:
        "wm", "wm_dilated_1vox_6conn" (default), "wm_dilated_1vox_14conn"
        or "nodif_brain" (whole brain).
        Two of the proposed white matter masks are dilated because a non-dilated
        white matter mask does not overlap with the "gray" subcortical regions,
        therefore the tracts will never get there. Moreover the right and left
        white matter regions are much less connected without dilation, therefore
        the connectogram shows few interhemisphere connections with a simple
        white matter mask.
    subjects_dir: str, defaut None
        If the Freesurfer $SUBJECTS_DIR environment variable is not set, or to
        bypass it, pass the path.
    subdir: str, default "masks"
        If 'subdir' is set the masks are saved in <outdir>/<subdir>/.
        Otherwise they are saved in <outdir>/.

    <unit>

        <input name="outdir"           type="Directory" />
        <input name="nodif_brain"      type="File"      />
        <input name="nodif_brain_mask" type="File"      />
        <input name="dif2anat_dat"     type="File"      />
        <input name="subject_id"       type="Str"       />
        <input name="cortical_atlas"   type="Str"       />
        <input name="tracto_mask_type" type="Str"       />
        <input name="subjects_dir"     type="Directory" />
        <input name="subdir"           type="Str"       />

        <output name="roi_masks_txt"   type="File"      />
        <output name="tracto_mask"     type="File"      />
        <output name="wm_mask"         type="File"      />
        <output name="stop_mask"       type="File"      />
        <output name="avoid_mask"      type="File"      />

    </unit>
    """

    # Check arguments
    if not cortical_atlas in CORTICAL_ATLASES:
        raise ValueError("Bad 'cortical_atlas': {}, should be in {}"
                         .format(cortical_atlas, CORTICAL_ATLASES))

    if tracto_mask_type not in TRACTO_MASK_TYPES:
        raise ValueError("Bad argument 'tracto_mask_type': {}, should be in {}"
                         .format(tracto_mask_type, TRACTO_MASK_TYPES))

    # Freesurfer 'subjects_dir' has to be passed or set as environment variable
    subjects_dir = get_or_check_freesurfer_subjects_dir(subjects_dir)

    # If requested use a subdirectory in outdir
    if subdir:
        outdir = os.path.join(outdir, subdir)

    # If outdir does not exist, create it
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Project cortical and subcortical segmentation in diffusion
    aparc2dif, aseg2dif = \
        project_aparc_and_aseg_to_diffusion(outdir         = outdir,
                                            dif2anat_dat   = dif2anat_dat,
                                            nodif_brain    = nodif_brain,
                                            subject_id     = subject_id,
                                            subjects_dir   = subjects_dir,
                                            cortical_atlas = cortical_atlas)

    # Create the target ROI masks
    if cortical_atlas == "Desikan":
        target_ROIs = LAUSANNE2008_SCALE33_ROIS
    else:
        target_ROIs = DESTRIEUX_WITH_SUBCORTICAL_ROIS
    roi_masks = create_target_masks(outdir, target_ROIs, aparc2dif)

    # Create white matter mask
    wm_mask = create_white_matter_mask(outdir, aseg2dif)

    # Create the tracto, according to the requested tracto mask type
    if tracto_mask_type == "wm":
        tracto_mask = wm_mask
    elif tracto_mask_type == "wm_dilated_1vox_6conn":
        tracto_mask = dilate_mask_by_one_voxel(wm_mask, connexity=6)
    elif tracto_mask_type == "wm_dilated_1vox_14conn":
        tracto_mask = dilate_mask_by_one_voxel(wm_mask, connexity=14)
    else:  # tracto_mask_type == TractoMaskTypes.nodif_brain
        tracto_mask = nodif_brain_mask

    # Create "avoid" mask: mask of the ventricles
    avoid_mask = create_ventricles_mask(outdir, aseg2dif)

    # Create tracto stop mask:
    stop_mask = None

    # Write the list in a txt file (probtrackx2 takes a txt list as input)
    roi_masks_txt = os.path.join(outdir, "roi_masks.txt")
    np.savetxt(roi_masks_txt, roi_masks, fmt="%s")

    return roi_masks_txt, tracto_mask, wm_mask, stop_mask, avoid_mask


def create_masks_for_tracto_seeding_wm(outdir,
                                       nodif_brain,
                                       nodif_brain_mask,
                                       dif2anat_dat,
                                       subject_id,
                                       cortical_atlas = "Desikan",
                                       stop_mask_type = "target_rois",
                                       subjects_dir   = None,
                                       subdir         = "masks"):
    """
    Create the volume masks required for the probabilist tractography
    when using the --omatrix3 option in probtrackx2 (seeding in white matter):
        - tractography mask is nodif_brain_mask
        - ROI masks, depends on the cortical_atlas, can include the
          subcortical regions (from Freesurfer 'aseg.mgz')
        - seed mask: eroded white matter.
        - stop mask: inverse of dilated white matter

    Parameters
    ----------
    outdir: str
        Path to directory where to output.
    nodif_brain: str
        Path to the brain-only volume extracted from the preprocessed DWI.
    nodif_brain_mask: str
        Path to the brain-only binary mask.
    dif2anat_dat: str
        Path to the registration .dat file, diffusion to anatomy transformation.
    subject_id: str
        Id of the subject in the Freesurfer subjects_dir.
    cortical_atlas: str, default "Desikan"
        Name of the freesurfer cortical parcellation to use, either "Desikan"
        (default) or "Destrieux". The corpus callosum is not included.
    stop_mask_type: str, default "target_rois"
        What type of stop mask to create:
        - "target_rois": stop a sample as soon as it reaches a target region
        - "inverse_wm":  stop a sample as soon as it leaves the white matter
    subjects_dir: str, defaut None
        If the Freesurfer $SUBJECTS_DIR environment variable is not set, or to
        bypass it, pass the path.
    subdir: str, default "masks"
        If 'subdir' is set the masks are saved in <outdir>/<subdir>/.
        Otherwise they are saved in <outdir>/.

    <unit>

        <input name="outdir"           type="Directory" />
        <input name="nodif_brain"      type="File"      />
        <input name="nodif_brain_mask" type="File"      />
        <input name="dif2anat_dat"     type="File"      />
        <input name="subject_id"       type="Str"       />
        <input name="cortical_atlas"   type="Str"       />
        <input name="stop_mask_type"   type="Str"       />
        <input name="subjects_dir"     type="Directory" />
        <input name="subdir"           type="Str"       />

        <output name="roi_masks_txt"   type="File"      />
        <output name="tracto_mask"     type="File"      />
        <output name="wm_mask"         type="File"      />
        <output name="stop_mask"       type="File"      />
        <output name="avoid_mask"      type="File"      />

    </unit>
    """
    # Check arguments
    if not cortical_atlas in CORTICAL_ATLASES:
        raise ValueError("Bad 'cortical_atlas': {}, should be in {}"
                         .format(cortical_atlas, CORTICAL_ATLASES))

    if not stop_mask_type in STOP_MASK_TYPES:
        raise ValueError("Bad 'stop_mask_type': {}, should be in {}"
                         .format(stop_mask_type, STOP_MASK_TYPES))

    # Freesurfer 'subjects_dir' has to be passed or set as environment variable
    subjects_dir = get_or_check_freesurfer_subjects_dir(subjects_dir)

    # If requested use a subdirectory in outdir
    if subdir:
        outdir = os.path.join(outdir, subdir)

    # If outdir does not exist, create it
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Project cortical and subcortical segmentation in diffusion
    aparc2dif, aseg2dif = \
        project_aparc_and_aseg_to_diffusion(outdir         = outdir,
                                            dif2anat_dat   = dif2anat_dat,
                                            nodif_brain    = nodif_brain,
                                            subject_id     = subject_id,
                                            subjects_dir   = subjects_dir,
                                            cortical_atlas = cortical_atlas)

    # Create the target ROI masks
    if cortical_atlas == "Desikan":
        target_ROIs = LAUSANNE2008_SCALE33_ROIS
    else:
        target_ROIs = DESTRIEUX_WITH_SUBCORTICAL_ROIS
    roi_masks = create_target_masks(outdir, target_ROIs, aparc2dif)

    # Create seed mask: white matter mask
    wm_mask = create_white_matter_mask(outdir, aseg2dif)

    # Create "avoid" mask: mask of the ventricles
    avoid_mask = create_ventricles_mask(outdir, aseg2dif)

    # Create the tractography stop mask
    stop_mask = os.path.join(outdir, "%s_mask.nii.gz" % stop_mask_type)
    if stop_mask_type == "target_rois":
        # The stop mask is the combination of all target_ROIs, so that the
        # samples stop propagating as soon as they reach a target region.
        labels = [LABEL_OF_ROI[x] for x in target_ROIs]
        cmd    = ["mri_binarize", "--i", aparc2dif, "--o", stop_mask,
                   "--match"] + labels
    else:  #  stop_mask_type == "inverse_wm":
        cmd = ["mri_binarize", "--i", aseg2dif, "--o", stop_mask, "--wm", "--inv"]
    run_freesurfer_cmd(cmd)

    # Write the list in a txt file (probtrackx2 takes a txt list as input)
    roi_masks_txt = os.path.join(outdir, "roi_masks.txt")
    np.savetxt(roi_masks_txt, roi_masks, fmt="%s")

    tracto_mask = nodif_brain_mask

    return roi_masks_txt, tracto_mask, wm_mask, stop_mask, avoid_mask


def compute_surface_area_of_roi_masks(roi_masks_txt):
    """
    We want the surface area of all ROI masks because we use them to normalize
    the connectogram.

    Steps for each ROI mask:
        - tessellate with Freesurfer mri_tessellate cmd to get a surface
        - call mris_info to get the total_area of the mask.

    Parameters
    ----------
    roi_masks_txt: str
        Path to .txt file listing all the paths to the masks.
    filename: str, default 'area_of_roi'
        The name without extension of the output JSON.

    Returns
    -------
    area_of_roi: dict
        Maps <ROI name> -> <surface area of mask>
    """
    # Create the surfaces in a temp dir
    tempdir = tempfile.mkdtemp()

    # Load the paths
    roi_masks = np.loadtxt(roi_masks_txt, dtype=str)

    # Dict mapping <ROI> -> <surface area>
    area_of_roi = dict()

    # Pattern to use to extract the total_area from the mris_info stdout
    pattern = re.compile("^total_area[ ]+(\d+.\d+)", flags=re.M|re.I)

    for roi_mask in roi_masks:
        roi_name     = os.path.basename(roi_mask).split(".nii")[0]
        path_surface = os.path.join(tempdir, roi_name)
        cmd = ["mri_tessellate", roi_mask, "1", path_surface]
        run_freesurfer_cmd(cmd)

        cmd_info = ["mris_info", path_surface]
        stdout   = run_freesurfer_cmd(cmd_info).stdout

        # Look for the pattern
        matches = pattern.findall(stdout)
        try:
            area = float(matches[0])
        except:
            raise Exception("Failed to compute the surface area of mask: %s " % roi_mask)

        area_of_roi[roi_name] = area

    # Remove temp dir
    shutil.rmtree(tempdir)

    return area_of_roi


# TODO: complete with more snap shots
def qc_tracto_masks(outdir,
                    nodif_brain,
                    tracto_mask,
                    wm_mask,
                    stop_mask,
                    roi_masks_txt,
                    subdir = "qc"):
    """
    Function meant to help quality check (qc) the masks created by the
    create_masks_for_omatrix#_tractography() functions.
    It creates snap shots to visualize the quality of the registration
    of the tractography masks in the diffusion space.

    The snap shots are saved in <outdir>/<subdir>/<snap shot>. By default
    <subdir> is "qc". To write in outdir directly, set subdir to anything that
    evaluates to False (empty string or None).
    Directories are automatically created if they don't exist.

    <unit>
        <input name="outdir"        type="Directory" />
        <input name="nodif_brain"   type="File"      />
        <input name="tracto_mask"   type="File"      />
        <input name="wm_mask"       type="File"      />
        <input name="stop_mask"     type="File"      />
        <input name="roi_masks_txt" type="File"      />
        <input name="subdir"        type="Str"       />

        <output name="qc_dir"       type="Directory" />
    </unit>
    """
    # If requested use a subdirectory in outdir
    if subdir:
        outdir = os.path.join(outdir, subdir)

    # If outdir does not exist, create it
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    nb_slices_in_z = nibabel.load(nodif_brain).get_shape()[2]

    # Snap shot: nodif_brain with stop mask
    snapshot1_png = os.path.join(outdir, "nodif_brain_with_stop_mask.png")
    plot_image(nodif_brain,
               overlay_file = stop_mask,
               snap_file    = snapshot1_png,
               name         = "nodif_brain_with_stop_mask",
               cut_coords   = nb_slices_in_z-2)

    # Snap shot: nodif_brain with wm mask
    snapshot2_png = os.path.join(outdir, "nodif_brain_with_wm_mask.png")
    plot_image(nodif_brain,
               overlay_file = wm_mask,
               snap_file    = snapshot2_png,
               name         = "nodif_brain_with_wm_mask",
               cut_coords   = nb_slices_in_z-2)

    # Return something for Capsul
    qc_dir = outdir
    return qc_dir


def extract_nodif_volume(outdir, dwi, bval):
    """

    Parameters
    ----------
    outdir: str
        Path to directory where to write "nodif.nii.gz"
    dwi: str
        Path to DW data in which at least one volume was acquired with bvalue=0.
    bval: str
        Path to .bval file associated to the DW data.

    Return
    ------
    nodif_volume: str
        Path to a/the volume for which bvalue is 0.
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
    outdir: str
        Path to the directory where to output.
    dwi: str
        Path to DW data in which at least one volume was acquired with bvalue=0.
    bval: str
        Path to .bval file associated to the DW data.
    subdir: str, default "bet2_nodif_brain"
        If you want the result files to be written in a subdirectory, specify
        the name.

    Return
    ------
    nodif_brain: str
        Path to the brain only volume.
    nodif_brain_mask: str
        Path to the brain-only binary mask.

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
    Bedpostx requires an input directory with specific filenames inside to run.
    The function copies the required files to a directory with the right names.
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


def omatrix3_to_roi_network(probtrackx2_dir, outdir=None):
    """
    When using the --omatrix3 option in probtrackx2, the result is a
    matrix: VOXELxVOXEL but we want ROIxROI. This function is meant to make
    the conversion.

    Parameters
    ----------
    probtrackx2_dir: str
        Path to dir where to find the files created by probtrackx2 when using
        --omatrix3 option, i.e "fdt_matrix3.dot" and "coords_for_fdt_matrix3"
    outdir: str, default None
        Path to dir where to output file (matrix as text file).
        By default in the same directory as in the input.

    Returns
    -------
    matrix_txt: str
        Path to the output file.

    """

    # Check input and output dirs
    if not os.path.isdir(probtrackx2_dir):
        raise ValueError("Input directory does not exist: {}".format(probtrackx2_dir))

    if outdir is None:
        outdir = probtrackx2_dir
    else:
        if not os.path.isdir(outdir):
            raise ValueError("Outdir does not exist: {}".format(probtrackx2_dir))

    path_matrix3 = os.path.join(probtrackx2_dir, "fdt_matrix3.dot")
    path_coords  = os.path.join(probtrackx2_dir, "coords_for_fdt_matrix3")

    matrix3 = np.loadtxt(path_matrix3, dtype=int)
    coords  = np.loadtxt(path_coords, dtype=int)
    nb_ROIs = coords[:, 3].max() + 1  # index start at 0, add 1

    # Map <voxel id> to <ROI id>
    ROI_of_voxel = dict()
    nb_voxels = coords.shape[0]
    for i in range(nb_voxels):
        ROI_of_voxel[coords[i, 4]] = coords[i, 3]

    # ROIxROI fiber counts
    matrix = np.zeros((nb_ROIs, nb_ROIs), dtype=int)

    # For each pair of voxels add the fiber count to the associated ROIs
    nb_voxel_pairs = matrix3.shape[0]
    for i in range(nb_voxel_pairs):
        ROI_1 = ROI_of_voxel[matrix3[i,0]]
        ROI_2 = ROI_of_voxel[matrix3[i, 1]]
        fiber_count = matrix3[i, 2]

        # Ignore internal ROI connections and symmetrize the matrix
        if ROI_1 != ROI_2:
            matrix[ROI_1, ROI_2] += fiber_count
            matrix[ROI_2, ROI_1] += fiber_count

    # Write output
    matrix_txt = os.path.join(outdir, "matrix3")
    np.savetxt(matrix_txt, matrix, fmt="%i")

    return matrix_txt


def normalize_connectogram_by_target_surf_area(outdir,
                                               roi_masks_txt,
                                               matrix_txt,
                                               labels_txt,
                                               suffix="_norm"):
    """
    """
    area_of_roi = compute_surface_area_of_roi_masks(roi_masks_txt)
    matrix      = np.loadtxt(matrix_txt)
    labels      = np.loadtxt(labels_txt, dtype=str)

    assert matrix.shape[0] == matrix.shape[1] == labels.shape[0]

    # Normalize each value by the area of the target ROI
    for i, roi in enumerate(labels):
        matrix[:, i] = matrix[:, i] / area_of_roi[roi]

    # Save output matrix
    filename = "%s%s" % (os.path.basename(matrix_txt), suffix)
    matrix_norm_txt = os.path.join(outdir, filename)
    np.savetxt(matrix_norm_txt, matrix, fmt="%s")

    # Save area_of_roi as a JSON file
    out_json = os.path.join(outdir, "area_of_roi.json")
    with open(out_json, "w") as f:
        json.dump(area_of_roi, f)

    return matrix_norm_txt, out_json


def probtrackx2_connectogram_seeding_wm(outdir,
                                        bedpostx_dir,
                                        roi_masks_txt,
                                        tracto_mask,
                                        wm_mask,
                                        stop_mask  = None,
                                        avoid_mask = None,
                                        subdir     = "probtrackx2",
                                        nsamples   = 5000,
                                        nsteps     = 2000,
                                        cthr       = None,
                                        loopcheck  = True,
                                        steplength = 0.5,
                                        fibthresh  = None,
                                        distthresh = None,
                                        sampvox    = None
                                        ):
    """
    Generate a connectogram using probtrackx2 with the following approach:
        - seed the white matter with <nsamples> in each voxel
        - for each sample, bidirectionally propagate a fiber until both ends
          of the fiber end in stop mask.

    Parameters
    ----------
    outdir: str
        Path to the directory where to output.
    bedpostx_dir: str
        Path to the bedpostx output dir.
    roi_masks_txt: str
        Path to the txt file listing the paths to the ROI masks.
    tracto_mask: str
        Path to the tractography mask (nodif_brain_mask.nii.gz).
    wm_mask: str
        Path to the white matter mask.
    stop_mask: str
        Path to the stop mask (inverse of white matter).
    subdir: str, default "probtrackx2"
        If set, the outputs are written in a subdirectory of outdir.
    <other args>:
        'probtrackx2 --help'

    <unit>
        <output name="fiber_density"   type="File"      />
        <output name="matrix_txt"      type="File"      />
        <output name="matrix_norm_txt" type="File"      />
        <output name="labels_txt"      type="File"      />

        <input name="outdir"           type="Directory" />
        <input name="bedpostx_dir"     type="Directory" />
        <input name="roi_masks_txt"    type="File"      />
        <input name="tracto_mask"      type="File"      />
        <input name="wm_mask"          type="File"      />
        <input name="stop_mask"        type="File"      />
        <input name="avoid_mask"       type="File"      />
        <input name="subdir"           type="Str"       />
        <input name="nsamples"         type="Int"       />
        <input name="nsteps"           type="Int"       />
        <input name="cthr"             type="Float"     />
        <input name="loopcheck"        type="Bool"      />
        <input name="steplength"       type="Float"     />
        <input name="fibthresh"        type="Float"     />
        <input name="distthresh"       type="Float"     />
        <input name="sampvox"          type="Float"     />
    </unit>
    """

    if subdir:
        outdir = os.path.join(outdir, subdir)

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    probtrackx2(dir        = outdir,
                samples    = os.path.join(bedpostx_dir, "merged"),
                mask       = tracto_mask,
                seed       = wm_mask,
                omatrix3   = True,
                target3    = roi_masks_txt,
                stop       = stop_mask,
                avoid      = avoid_mask,
                nsamples   = nsamples,
                nsteps     = nsteps,
                cthr       = cthr,
                loopcheck  = loopcheck,
                steplength = steplength,
                fibthresh  = fibthresh,
                distthresh = distthresh,
                sampvox    = sampvox
                )

    fiber_density = os.path.join(outdir, "fdt_paths.nii.gz")

    # Convert omatrix3 output (VOXELxVOXEL fiber count) to ROIxROI fiber count
    matrix_txt = omatrix3_to_roi_network(outdir)

    # Create a list of labels from the list of mask paths
    roi_masks  = np.loadtxt(roi_masks_txt, dtype=str)
    labels     = [os.path.basename(x).split(".nii")[0] for x in roi_masks]
    labels_txt = os.path.join(outdir, "labels.txt")
    np.savetxt(labels_txt, labels, fmt="%s")

    # Create a normalized connectogram in the same directory
    matrix_norm_txt, _ = normalize_connectogram_by_target_surf_area(outdir,
                                                                    roi_masks_txt,
                                                                    matrix_txt,
                                                                    labels_txt)

    return fiber_density, matrix_txt, matrix_norm_txt, labels_txt


def probtrackx2_connectogram_seeding_endpoints(outdir,
                                               bedpostx_dir,
                                               roi_masks_txt,
                                               tracto_mask,
                                               network         = True,
                                               stop_mask       = None,
                                               avoid_mask      = None,
                                               subdir          = "probtrackx2",
                                               nsamples        = 5000,
                                               nsteps          = 2000,
                                               cthr            = None,
                                               loopcheck       = True,
                                               onewaycondition = None,
                                               steplength      = 0.5,
                                               fibthresh       = None,
                                               distthresh      = None,
                                               sampvox         = None
                                               ):
    """
    Generate a connectogram using probtrackx2 with the following approach:
        - seed the ROIs (the endpoints) with <nsamples> in each voxel
        - for each sample, propagate until stop
        - the --network option automatically generates a ROIxROI matrix
          called 'fdt_network_matrix'
    <unit>
        <output name="fiber_density"   type="File"      />
        <output name="matrix_txt"      type="File"      />
        <output name="matrix_norm_txt" type="File"      />
        <output name="labels_txt"      type="File"      />

        <input name="outdir"           type="Directory" />
        <input name="bedpostx_dir"     type="Directory" />
        <input name="roi_masks_txt"    type="File"      />
        <input name="tracto_mask"      type="File"      />
        <input name="network"          type="Bool"      />
        <input name="stop_mask"        type="File"      />
        <input name="avoid_mask"       type="File"      />
        <input name="subdir"           type="Str"       />
        <input name="nsamples"         type="Int"       />
        <input name="nsteps"           type="Int"       />
        <input name="cthr"             type="Float"     />
        <input name="loopcheck"        type="Bool"      />
        <input name="onewaycondition"  type="Bool"      />
        <input name="steplength"       type="Float"     />
        <input name="fibthresh"        type="Float"     />
        <input name="distthresh"       type="Float"     />
        <input name="sampvox"          type="Float"     />
    </unit>
    """

    if subdir:
        outdir = os.path.join(outdir, subdir)

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    basepath_samples = os.path.join(bedpostx_dir, "merged")

    list_of_files, matrix_txt = probtrackx2(network         = network,
                                            seed            = roi_masks_txt,
                                            stop            = stop_mask,
                                            avoid           = avoid_mask,
                                            loopcheck       = loopcheck,
                                            onewaycondition = onewaycondition,
                                            samples         = basepath_samples,
                                            mask            = tracto_mask,
                                            dir             = outdir,
                                            nsamples        = nsamples,
                                            nsteps          = nsteps,
                                            cthr            = cthr,
                                            steplength      = steplength,
                                            fibthresh       = fibthresh,
                                            distthresh      = distthresh,
                                            sampvox         = sampvox)
    fiber_density = list_of_files[0]

    # Create a list of labels from the list of mask paths
    roi_masks  = np.loadtxt(roi_masks_txt, dtype=str)
    labels     = [os.path.basename(x).split(".nii")[0] for x in roi_masks]
    labels_txt = os.path.join(outdir, "labels.txt")
    np.savetxt(labels_txt, labels, fmt="%s")

    # Create a normalized connectogram
    matrix_norm_txt, _ = normalize_connectogram_by_target_surf_area(outdir,
                                                                    roi_masks_txt,
                                                                    matrix_txt,
                                                                    labels_txt)

    return fiber_density, matrix_txt, matrix_norm_txt, labels_txt


def plot_connectogram(outdir, matrix, labels=None, transform=None,
                      prefix="connectogram", dpi=200):
    """
    Inspired from the plot_matrix() function, rewritten to make plots for the
    'small' connectogram, with the region labels.

    Parameters
    ----------
    outdir: str
        Path to the dir where to save the snapshot.
    matrix: numpy array or str
        Connectivity matrix or path to txt file storing the connectivity matrix.
    labels: list of str or str, default None
        Labels or path to txt file listing the labels.
        By default no labels. Should be ordered like the rows/cols of the matrix.
    transform: callable, default None
        A Callable function to apply on the matrix (e.g. numpy.log1p).
        By default no transformation is applied.
    prefix: str, default "connectogram"
        Filename without extension of output pdf.
    dpi: int, default 200
        "Dot Per Inch", set higher for better resolution.

    Returns
    -------
    snap_file: str
        Path to the output snap shot.
    """

    # matrix should be either a txt file storing the matrix or a numpy array
    if isinstance(matrix, basestring):
        matrix = np.loadtxt(matrix)

    # check matrix dimensions
    assert matrix.shape[0] == matrix.shape[1], "Square matrix only"

    # Apply transformation if requested
    if transform is not None:
        matrix = transform(matrix)

    ##################
    # Create the figure with matplotlib
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(matrix, cmap=plt.cm.Reds)

    ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticks(np.arange(matrix.shape[0])+0.5, minor=False)
    ax.set_yticks(np.arange(matrix.shape[1])+0.5, minor=False)

    # Add the labels if passed
    if labels is not None:
        if isinstance(labels, basestring):
            labels = np.loadtxt(labels, dtype=str)

        assert len(labels) == matrix.shape[0], "Wrong number of labels."

        ax.set_xticklabels(labels, minor=False, size=4)
        ax.set_yticklabels(labels, minor=False, size=4)

        ax.set_xticklabels(labels, minor=False)
        ax.set_yticklabels(labels, minor=False)

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
    snap_file = os.path.join(outdir, prefix)
    fig.savefig(snap_file, dpi=dpi)

    # Release memory
    fig.clear()
    plt.close()

    return snap_file


def qc_connectogram(outdir,
                    tracto_mask,
                    fiber_density,
                    matrix_txt,
                    matrix_norm_txt,
                    labels_txt,
                    subdir="qc"):
    """
    Function meant to help quality check (qc) the tractography and the
    connectogram computed by the probtrackx2_probabilist_tractography function.
    It creates snap shots to visualize the connectogram and the fiber density.

    The snap shots are saved in <outdir>/<subdir>/<snap shot>. By default
    <subdir> is "qc". To write in outdir directly, set subdir to anything that
    evaluates to False (empty string or None).
    Directories are automatically created if they don't exist.

    'fiber_density' and 'matrix' are outputs of the
    probtrackx2_omatrix#_connectogram() functions.

    <unit>
        <input name="outdir"          type="Directory" />
        <input name="tracto_mask"     type="File"      />
        <input name="fiber_density"   type="File"      />
        <input name="matrix_txt"      type="File"      />
        <input name="matrix_norm_txt" type="File"      />
        <input name="labels_txt"      type="File"      />
        <input name="subdir"          type="Str"       />

        <output name="qc_dir"         type="Directory" />
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
               overlay_file = fiber_density,
               snap_file    = fiber_density_png,
               name         = "fiber density map",
               overlay_cmap = "cold_hot",
               cut_coords   = nb_slices_in_z-2)

    # connectogram snapshot
    plot_connectogram(outdir, matrix_txt, labels=labels_txt, transform=np.log1p)

    # Normalized by surface target area connectogram snapshot
    plot_connectogram(outdir, matrix_norm_txt, labels=labels_txt,
                      transform=np.log1p, prefix="connectogram_normalized")

    # Return something for Capsul
    qc_dir = outdir
    return qc_dir
