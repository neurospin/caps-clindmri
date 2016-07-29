# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################


import os
import shutil
import tempfile
import re
import json

import numpy as np
import nibabel
import matplotlib.pyplot as plt

from clindmri.extensions.configuration import environment
from clindmri.extensions.freesurfer.wrappers import FSWrapper
from clindmri.extensions.freesurfer.exceptions import FreeSurferRuntimeError
from clindmri.estimation.fsl import dtifit
from clindmri.tractography.fsl import probtrackx2
from clindmri.plot.slicer import plot_image


###############################################################################
# MODULE GLOBAL VARIABLES

# Get Freesurfer Look Up Table path
if "FREESURFER_HOME" in os.environ:
    PATH_FREESURFER_HOME = os.environ["FREESURFER_HOME"]
    PATH_LUT = os.path.join(PATH_FREESURFER_HOME, "FreeSurferColorLUT.txt")
else:
    raise Exception("Environment variable 'FREESURFER_HOME' is not set.")

# Load Freesurfer Look Up Table and create a dict mapping: <ROI name> -> label
try:
    LABEL_OF_ROI = dict(np.loadtxt(PATH_LUT, dtype=str, usecols=[1, 0]))
except:
    raise Exception("Failed to load Freesurfer Look Up Table: %s" % PATH_LUT)


# Desikan atlas left regions without corpus callosum ordered like Lausanne
# 2008 atlas
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

# Desikan atlas left regions without corpus callosum ordered like Lausanne
# 2008 atlas
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
DESTRIEUX_LH_LABELS = [str(x) for x in range(DESTRIEUX_CTX_LH_MIN_LABEL,
                                             DESTRIEUX_CTX_LH_MAX_LABEL)]
DESTRIEUX_RH_LABELS = [str(x + 1000) for x in range(
    DESTRIEUX_CTX_LH_MIN_LABEL, DESTRIEUX_CTX_LH_MAX_LABEL)]

DESTRIEUX_LH_ROIS = [x for x in LABEL_OF_ROI
                     if LABEL_OF_ROI[x] in set(DESTRIEUX_LH_LABELS)]
DESTRIEUX_RH_ROIS = [x for x in LABEL_OF_ROI
                     if LABEL_OF_ROI[x] in set(DESTRIEUX_RH_LABELS)]

DESTRIEUX_WITH_SUBCORTICAL_ROIS = (DESTRIEUX_LH_ROIS + LH_SUBCORTICAL_ROIS +
                                   DESTRIEUX_RH_ROIS + RH_SUBCORTICAL_ROIS +
                                   AXIAL_SUBCORTICAL_ROIS)

# Set of available atlas
CORTICAL_ATLASES = frozenset(["Desikan", "Destrieux"])

# Set of available types of stop masks when using --target3 tractography
STOP_MASK_TYPES = {"target_rois", "inverse_wm"}

###############################################################################
# Utility functions


def get_or_check_freesurfer_subjects_dir(subjects_dir=None):
    """
    If 'subjects_dir' is passed, check whether the directory exists, otherwise
    look for the $SUBJECTS_DIR environment variable. If 'subjects_dir' is not
    passed and $SUBJECTS_DIR not in the environment, raise an Exception.
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
    To avoid repeating the code to run Freesurfer and check exitcode
    everywhere.
    Step:
        - add $SUBJECTS_DIR to the environment if requested
        - add FSL's environment if requested (some Freesurfer commands require
          FSL)
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


###############################################################################


def register_diffusion_to_anatomy(outdir,
                                  nodif_brain,
                                  subject_id,
                                  subjects_dir=None,
                                  subdir="diff_to_anat",
                                  fsl_init="/etc/fsl/5.0/fsl.sh"):
    """
    Register the diffusion to the anatomy (T1) using Freesurfer bbregister
    (boundary-based registration).

    The resulting .dat file is saved in <outdir>/<subdir>/dif2anat.dat
    To write in 'outdir' directly, set 'subdir' to anything that evaluates
    to False (empty string or None).

    'subjects_dir' has to be passed if not set as an environnement variable
    ($SUBJECTS_DIR).
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
    cmd = ["bbregister",
           "--s",   subject_id,
           "--mov", nodif_brain,
           "--reg", dif2anat_dat,
           "--dti",
           "--init-fsl"]
    run_freesurfer_cmd(cmd, subjects_dir=subjects_dir, add_fsl_env=True)

    return dif2anat_dat


def qc_dif2anat_registration(outdir,
                             nodif_brain,
                             dif2anat_dat,
                             subject_id,
                             subjects_dir=None,
                             subdir="qc"):
    """
    Function meant to help quality check (qc) the registration between
    the diffusion and the anatomy (T1 from Freesurfer recon-all).
    It creates snap shots:
        - T1 brain registered in diffusion + contour of nodif brain volume.

    The snap shot is saved in <outdir>/<subdir>/"t1_to_diff.png". By default
    <subdir> is "qc". To write in outdir directly, set subdir to anything that
    evaluates to False (empty string or None).
    Directories are automatically created if they don't exist.
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
    t1 = os.path.join(subjects_dir, subject_id, "mri/brainmask.mgz")
    t1_to_dif = os.path.join(outdir, "t1_to_dif.nii.gz")

    # Project T1 in dif
    cmd = ["mri_vol2vol",
           "--mov",    nodif_brain,
           "--targ",   t1,
           "--inv",
           "--interp", "nearest",
           "--o",      t1_to_dif,
           "--reg",    dif2anat_dat,
           "--no-save-reg"]
    run_freesurfer_cmd(cmd)

    nb_slices_in_z = nibabel.load(nodif_brain).get_shape()[2]

    # PDF snapshots: registered T1 (in diffusion) + nodif edges
    pdf_t1_with_nodif_edges = os.path.join(outdir, "t1_with_nodif_edges.pdf")
    plot_image(t1_to_dif,
               edge_file=nodif_brain,
               snap_file=pdf_t1_with_nodif_edges,
               name="T1 in diffusion + edges of nodif",
               cut_coords=nb_slices_in_z - 2)

    # PDF snapshots: nodif brain + registered T1 edges
    pdf_nodif_with_t1_edges = os.path.join(outdir, "nodif_with_t1_edges.pdf")
    plot_image(nodif_brain,
               edge_file=t1_to_dif,
               snap_file=pdf_nodif_with_t1_edges,
               name="nodif + edges of registered T1",
               cut_coords=nb_slices_in_z - 2)

    return outdir


def project_aparc_and_aseg_to_diffusion(outdir,
                                        dif2anat_dat,
                                        nodif_brain,
                                        subject_id,
                                        subjects_dir=None,
                                        cortical_atlas="Desikan",
                                        outext=".nii.gz"):
    """
    Apply the transform specified by dif2anat_dat to the cortical and
    subcortical Freesurfer segmentations (aparc and aseg).

    Parameters
    ----------
    outdir: str
        Directory where to output the 2 projections.
    dif2anat_dat: str
        Path to .dat file generated by the registration process (bbregister,
        tkregister...). Created by register_diffusion_to_anatomy().
    nodif_brain: str
        Path to the brain-only volume in diffusion space used when the
        registration was done.
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

    if cortical_atlas not in CORTICAL_ATLASES:
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
    cmd = ["mri_vol2vol",
           "--mov",    nodif_brain,
           "--targ",   aparc_mgz,
           "--inv",
           "--interp", "nearest",
           "--o",      aparc2dif,
           "--reg",    dif2anat_dat,
           "--no-save-reg"]
    run_freesurfer_cmd(cmd)

    # Project subcortical segmentation in diffusion
    aseg_mgz = os.path.join(subj_dir, "mri/aseg.mgz")
    aseg2dif = os.path.join(outdir, "aseg2dif%s" % outext)
    cmd = ["mri_vol2vol",
           "--mov",    nodif_brain,
           "--targ",   aseg_mgz,
           "--inv",
           "--interp", "nearest",
           "--o",      aseg2dif,
           "--reg",    dif2anat_dat,
           "--no-save-reg"]
    run_freesurfer_cmd(cmd)

    return aparc2dif, aseg2dif


def create_masks_for_tracto_seeding_wm(outdir,
                                       dwi,
                                       bval,
                                       bvec,
                                       nodif_brain,
                                       nodif_brain_mask,
                                       dif2anat_dat,
                                       subject_id,
                                       cortical_atlas="Desikan",
                                       stop_mask_type="target_rois",
                                       subjects_dir=None,
                                       subdir="masks"):
    """
    Create the volume masks required for the probabilist tractography when
    using the --omatrix3 option in probtrackx2 (seeding in white matter):
        - tractography mask is nodif_brain_mask
        - ROI masks, depends on the cortical_atlas, includes a set of selected
          subcortical regions (from Freesurfer 'aseg.mgz')
        - seed mask: white matter.
        - stop mask: inverse of white matter or target ROIs (stop as soon as a
                     target region is reached)

    Parameters
    ----------
    outdir: str
        Path to directory where to output.
    dwi: str
        Path to the diffusion-weighted images (Nifti required).
    bval: str
        Path to the bvalue list.
    bvec: str
        Path to the list of diffusion-sensitized directions.
    nodif_brain: str
        Path to the brain-only volume extracted from the preprocessed DWI.
    nodif_brain_mask: str
        Path to the brain-only binary mask.
    dif2anat_dat: str
        Path to the registration .dat file, diffusion to anatomy
        transformation.
    subject_id: str
        Id of the subject in the Freesurfer subjects_dir.
    cortical_atlas: str, optional
        Name of the freesurfer cortical parcellation to use, either "Desikan"
        (default) or "Destrieux". The corpus callosum is not included.
    stop_mask_type: str, optional
        What type of stop mask to create:
        - "target_rois": stop a sample as soon as it reaches a target region
        - "inverse_wm":  stop a sample as soon as it leaves the white matter
    subjects_dir: None or str, optional
        If the Freesurfer $SUBJECTS_DIR environment variable is not set, or to
        bypass it, pass the path.
    subdir: str, optional
        If 'subdir' is set the masks are saved in <outdir>/<subdir>/.
        Otherwise they are saved in <outdir>/.
    """
    # Check arguments
    if cortical_atlas not in CORTICAL_ATLASES:
        raise ValueError("Bad 'cortical_atlas': {}, should be in {}"
                         .format(cortical_atlas, CORTICAL_ATLASES))

    if stop_mask_type not in STOP_MASK_TYPES:
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

    # -------------------------------------------------------------------------
    # Project cortical and subcortical segmentation in diffusion
    aparc2dif, aseg2dif = \
        project_aparc_and_aseg_to_diffusion(outdir=outdir,
                                            dif2anat_dat=dif2anat_dat,
                                            nodif_brain=nodif_brain,
                                            subject_id=subject_id,
                                            subjects_dir=subjects_dir,
                                            cortical_atlas=cortical_atlas)

    # -------------------------------------------------------------------------
    # Create a mask for each target ROI
    if cortical_atlas == "Desikan":
        target_ROIs = LAUSANNE2008_SCALE33_ROIS
    else:
        target_ROIs = DESTRIEUX_WITH_SUBCORTICAL_ROIS

    roi_masks = []
    for roi in target_ROIs:
        mask_path = os.path.join(outdir, "%s.nii.gz" % roi)
        cmd = ["mri_binarize",
               "--i",     aparc2dif,
               "--match", LABEL_OF_ROI[roi],
               "--o",     mask_path]
        run_freesurfer_cmd(cmd)
        roi_masks.append(mask_path)

    # -------------------------------------------------------------------------
    # Create seed mask: mask of white matter voxels where FA > 0.2

    # 1 - Create white matter mask
    wm_mask = os.path.join(outdir, "wm_mask.nii.gz")
    cmd = ["mri_binarize", "--i", aseg2dif, "--wm", "--o", wm_mask]
    run_freesurfer_cmd(cmd)

    # 2 - Compute FA in white matter
    dtifit_dir = os.path.join(outdir, "dtifit_wm")
    wm_fa = dtifit(k=dwi, m=wm_mask, b=bval, r=bvec, o=dtifit_dir)[7]

    # 3 - Binarize to create a mask: the seed mask
    seed_mask = os.path.join(outdir, "wm_fa_thresh_0p2.nii.gz")
    cmd = ["mri_binarize",
           "--i",   wm_fa,
           "--min", "0.2",
           "--o",   seed_mask]
    run_freesurfer_cmd(cmd)

    # Create "avoid" mask: mask of the ventricles
    avoid_mask = os.path.join(outdir, "ventricles.nii.gz")
    cmd = ["mri_binarize",
           "--i", aseg2dif,
           "--ventricles",
           "--o", avoid_mask]
    run_freesurfer_cmd(cmd)

    # -------------------------------------------------------------------------
    # Create the tractography stop mask
    stop_mask = os.path.join(outdir, "%s_mask.nii.gz" % stop_mask_type)
    if stop_mask_type == "target_rois":
        # The stop mask is the combination of all target_ROIs, so that the
        # samples stop propagating as soon as they reach a target region.
        labels = [LABEL_OF_ROI[x] for x in target_ROIs]
        cmd = ["mri_binarize",
               "--i", aparc2dif,
               "--o", stop_mask,
               "--match"] + labels
    else:  # stop_mask_type == "inverse_wm":
        cmd = ["mri_binarize",
               "--i", aseg2dif,
               "--o", stop_mask,
               "--wm",
               "--inv"]
    run_freesurfer_cmd(cmd)

    # -------------------------------------------------------------------------
    # Write the paths to the ROI masks in a txt file
    # (probtrackx2 takes a txt list as input)
    txt_roi_masks = os.path.join(outdir, "roi_masks.txt")
    np.savetxt(txt_roi_masks, roi_masks, fmt="%s")

    tracto_mask = nodif_brain_mask

    return txt_roi_masks, tracto_mask, seed_mask, stop_mask, avoid_mask


def compute_surface_area_of_roi_masks(txt_roi_masks):
    """
    We want the surface area of all ROI masks because we use them to normalize
    the connectogram.

    Steps for each ROI mask:
        - tessellate with Freesurfer mri_tessellate cmd to get a surface
        - call mris_info to get the total_area of the mask.

    Parameters
    ----------
    txt_roi_masks: str
        Path to .txt file listing all the paths to the masks.

    Returns
    -------
    area_of_roi: dict
        Maps <ROI name> -> <surface area of mask>
    """
    # Create the surfaces in a temp dir
    tempdir = tempfile.mkdtemp()

    # Load the paths
    roi_masks = np.loadtxt(txt_roi_masks, dtype=str)

    # Dict mapping <ROI> -> <surface area>
    area_of_roi = dict()

    # Pattern to use to extract the total_area from the mris_info stdout
    pattern = re.compile("^total_area[ ]+(\d+.\d+)", flags=re.M | re.I)

    for roi_mask in roi_masks:
        roi_name = os.path.basename(roi_mask).split(".nii")[0]
        path_surface = os.path.join(tempdir, roi_name)
        cmd = ["mri_tessellate", roi_mask, "1", path_surface]
        run_freesurfer_cmd(cmd)

        cmd_info = ["mris_info", path_surface]
        stdout = run_freesurfer_cmd(cmd_info).stdout

        # Look for the pattern
        matches = pattern.findall(stdout)
        try:
            area = float(matches[0])
        except:
            raise Exception(
                "Failed to compute the surface area of mask: %s " % roi_mask)

        area_of_roi[roi_name] = area

    # Remove temp dir
    shutil.rmtree(tempdir)

    return area_of_roi


def qc_tracto_masks(outdir,
                    tracto_mask,
                    seed_mask,
                    stop_mask,
                    txt_roi_masks,
                    subdir="qc"):
    """
    Function meant to help quality check (qc) the masks created by the
    create_masks_for_tracto_seeding_*() functions.
    It creates snap shots to visualize the quality of the registration
    of the tractography masks in the diffusion space.

    The snap shots are saved in <outdir>/<subdir>/<snap shot>. By default
    <subdir> is "qc". To write in outdir directly, set subdir to anything that
    evaluates to False (empty string or None).
    Directories are automatically created if they don't exist.
    """
    # If requested use a subdirectory in outdir
    if subdir:
        outdir = os.path.join(outdir, subdir)

    # If outdir does not exist, create it
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    nb_slices_in_z = nibabel.load(tracto_mask).get_shape()[2]

    # PDF snapshot
    basename = "tracto_mask_and_stop_mask"
    pdf_snapshot_1 = os.path.join(outdir, "%s.pdf" % basename)
    plot_image(tracto_mask,
               overlay_file=stop_mask,
               snap_file=pdf_snapshot_1,
               name=basename,
               cut_coords=nb_slices_in_z - 2)

    # PDF snapshot
    basename = "tracto_mask_and_seed_mask"
    pdf_snapshot_2 = os.path.join(outdir, "%s.pdf" % basename)
    plot_image(tracto_mask,
               overlay_file=seed_mask,
               snap_file=pdf_snapshot_2,
               name=basename,
               cut_coords=nb_slices_in_z - 2)

    return outdir


def omatrix3_to_roi_network(probtrackx2_dir, outdir=None):
    """
    When using the --omatrix3 option in probtrackx2, the result is a
    matrix: VOXELxVOXEL but we want ROIxROI. This function does the conversion.

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
    txt_matrix: str
        Path to the output txt file.
    """

    # Check input and output dirs
    if not os.path.isdir(probtrackx2_dir):
        raise ValueError(
            "Input directory does not exist: {}".format(probtrackx2_dir))

    if outdir is None:
        outdir = probtrackx2_dir
    else:
        if not os.path.isdir(outdir):
            raise ValueError(
                "Outdir does not exist: {}".format(probtrackx2_dir))

    path_matrix3 = os.path.join(probtrackx2_dir, "fdt_matrix3.dot")
    path_coords = os.path.join(probtrackx2_dir, "coords_for_fdt_matrix3")

    matrix3 = np.loadtxt(path_matrix3, dtype=int)
    coords = np.loadtxt(path_coords, dtype=int)
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
        ROI_1 = ROI_of_voxel[matrix3[i, 0]]
        ROI_2 = ROI_of_voxel[matrix3[i, 1]]
        fiber_count = matrix3[i, 2]

        # Ignore intra ROI connections and symmetrize the matrix
        if ROI_1 != ROI_2:
            matrix[ROI_1, ROI_2] += fiber_count
            matrix[ROI_2, ROI_1] += fiber_count

    # Write output
    txt_matrix = os.path.join(outdir, "matrix3")
    np.savetxt(txt_matrix, matrix, fmt="%i")

    return txt_matrix


def normalize_connectogram_by_target_surf_area(outdir,
                                               txt_roi_masks,
                                               txt_matrix,
                                               txt_labels,
                                               suffix="_norm"):
    """
    If we call M the input connectivity and N the output normalized matrix
    Nij = Mij / (surface area of region i + surface area of region j)

    Parameters
    ----------
    txt_roi_masks: str
        Path to the txt file listing the paths to the ROI masks.
    txt_matrix: str
        Path to the txt file storing the connectivity matrix.
    txt_labels: str
        Path to the txt file listing the labels (=ROI names), ordered like the
        rows/columns of the matrix.
    """
    area_of_roi = compute_surface_area_of_roi_masks(txt_roi_masks)
    matrix = np.loadtxt(txt_matrix)
    labels = np.loadtxt(txt_labels, dtype=str)

    assert matrix.shape[0] == matrix.shape[1] == labels.shape[0]

    # Normalize each value by the area of the target ROI
    for i, roi_i in enumerate(labels):
        for j, roi_j in enumerate(labels):
            sum_areas = area_of_roi[roi_i] + area_of_roi[roi_j]
            matrix[i, j] = matrix[i, j] / sum_areas

    # Save output matrix
    filename = "%s%s" % (os.path.basename(txt_matrix), suffix)
    txt_matrix_normalized = os.path.join(outdir, filename)
    np.savetxt(txt_matrix_normalized, matrix, fmt="%s")

    # Save area_of_roi as a JSON file
    out_json = os.path.join(outdir, "area_of_roi.json")
    with open(out_json, "w") as f:
        json.dump(area_of_roi, f)

    return txt_matrix_normalized, out_json


def probtrackx2_connectogram_seeding_wm(outdir,
                                        bedpostx_dir,
                                        txt_roi_masks,
                                        tracto_mask,
                                        seed_mask,
                                        stop_mask,
                                        avoid_mask,
                                        nsamples,
                                        nsteps,
                                        steplength,
                                        cthr=None,
                                        loopcheck=True,
                                        fibthresh=None,
                                        distthresh=None,
                                        sampvox=None,
                                        waypoints=None,
                                        subdir="probtrackx2"):
    """
    Generate a connectogram using probtrackx2 with the following approach:
        - seed the white matter with <nsamples> in each voxel
        - for each sample, bidirectionally propagate a fiber until both ends
          of the fiber reach the stop mask.

    Parameters
    ----------
    outdir: str
        Path to the directory where to output.
    bedpostx_dir: str
        Path to the bedpostx output dir.
    txt_roi_masks: str
        Path to the txt file listing the paths to the ROI masks.
    tracto_mask: str
        Path to the tractography mask (nodif_brain_mask.nii.gz).
    seed_mask: str
        Path to the mask where to start the tractography. It corresponds to
        the voxels of white matter where the FA > 0.2.
    stop_mask: str
        Path to the stop mask (inverse of white matter).
    subdir: str, default "probtrackx2"
        If set, the outputs are written in a subdirectory of outdir.
    <other args>:
        'probtrackx2 --help'
    """

    if subdir:
        outdir = os.path.join(outdir, subdir)

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    probtrackx2(dir=outdir,
                samples=os.path.join(bedpostx_dir, "merged"),
                mask=tracto_mask,
                seed=seed_mask,
                omatrix3=True,
                target3=txt_roi_masks,
                stop=stop_mask,
                avoid=avoid_mask,
                nsamples=nsamples,
                nsteps=nsteps,
                cthr=cthr,
                loopcheck=loopcheck,
                steplength=steplength,
                fibthresh=fibthresh,
                distthresh=distthresh,
                sampvox=sampvox,
                waypoints=waypoints)

    fiber_density = os.path.join(outdir, "fdt_paths.nii.gz")

    # Convert omatrix3 output (VOXELxVOXEL fiber count) to ROIxROI fiber count
    txt_matrix = omatrix3_to_roi_network(outdir)

    # Create a list of labels from the list of mask paths
    roi_masks = np.loadtxt(txt_roi_masks, dtype=str)
    labels = [os.path.basename(x).split(".nii")[0] for x in roi_masks]
    txt_labels = os.path.join(outdir, "labels.txt")
    np.savetxt(txt_labels, labels, fmt="%s")

    # Create a normalized connectogram in the same directory
    txt_matrix_normalized, _ = normalize_connectogram_by_target_surf_area(
        outdir, txt_roi_masks, txt_matrix, txt_labels)

    return fiber_density, txt_matrix, txt_matrix_normalized, txt_labels


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
        Connectivity matrix or path to txt file storing the connectivity
        matrix.
    labels: list of str or str, default None
        Labels or path to txt file listing the labels.
        By default no labels. Should be ordered like the rows/cols of the
        matrix.
    transform: callable, default None
        A Callable function to apply on the matrix (e.g. numpy.log1p).
        By default no transformation is applied.
    prefix: str, default "connectogram"
        Filename without extension of output pdf.
    dpi: int, default 200
        "Dot Per Inch", set higher for better resolution.

    Returns
    -------
    connectogram_snapshot: str
        Path to the output snapshot.
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

    ax.set_xticks(np.arange(matrix.shape[0]) + 0.5, minor=False)
    ax.set_yticks(np.arange(matrix.shape[1]) + 0.5, minor=False)

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

    # Path to the output PDF file
    connectogram_snapshot = os.path.join(outdir, prefix)
    fig.savefig(connectogram_snapshot, dpi=dpi)

    # Release memory
    fig.clear()
    plt.close()

    return connectogram_snapshot


def qc_connectogram(outdir,
                    tracto_mask,
                    fiber_density,
                    txt_matrix,
                    txt_matrix_normalized,
                    txt_labels,
                    subdir="qc"):
    """
    Create 3 snapshots:
        - fiber density
        - raw connetivity matrix
        - normalized connectivity matrix

    The snap shots are saved in <outdir>/<subdir>/<snap shot>. By default
    <subdir> is "qc". To write in outdir directly, set subdir to anything that
    evaluates to False (empty string or None).
    Directories are automatically created if they don't exist.

    'fiber_density' and 'txt_matrix' are outputs of probtrackx2.
    """
    # If requested use a subdirectory in outdir
    if subdir:
        outdir = os.path.join(outdir, subdir)

    # If outdir does not exist, create it
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    nb_slices_in_z = nibabel.load(tracto_mask).get_shape()[2]

    fiber_density_snapshot = os.path.join(outdir, "fiber_density_map.pdf")
    plot_image(tracto_mask,
               overlay_file=fiber_density,
               snap_file=fiber_density_snapshot,
               name="fiber density map",
               overlay_cmap="cold_hot",
               cut_coords=nb_slices_in_z - 2)

    # connectogram snapshot
    connectogram_snapshot = plot_connectogram(outdir,
                                              txt_matrix,
                                              labels=txt_labels,
                                              transform=np.log1p)

    # Normalized by surface target area connectogram snapshot
    connectogram_norm_snapshot = \
        plot_connectogram(outdir,
                          txt_matrix_normalized,
                          labels=txt_labels,
                          transform=np.log1p,
                          prefix="connectogram_normalized")

    return connectogram_snapshot, connectogram_norm_snapshot
