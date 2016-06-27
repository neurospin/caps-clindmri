# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

import os
import subprocess
import time
import shutil

import nibabel

from clindmri.segmentation.fsl import bet2
from clindmri.plot.slicer import plot_image

from clindmri.connectivity.connectogram import (
    run_freesurfer_cmd, get_or_check_freesurfer_subjects_dir)

# Wrappers of Connectomist's tabs
from connectomist.import_and_qspace_model import (
    dwi_data_import_and_qspace_sampling)
from connectomist.mask import dwi_rough_mask_extraction
from connectomist.outliers import dwi_outlier_detection
from connectomist.eddy_current_and_motion import (
    dwi_eddy_current_and_motion_correction,
    export_eddy_motion_results_to_nifti)


def check_brainsuite_installation():
    """
    Check that Brainsuite commands, bfc and bdp.sh, are executable.
    If not raise an Exception.
    """
    # devnull: to kill output, not print in the console
    devnull = open(os.devnull, "w")

    try:
        subprocess.check_call(["bfc"], stdout=devnull)
    except:
        raise Exception("Brainsuite is not installed or bfc is not in $PATH.")

    try:
        subprocess.check_call(["bdp.sh"], stdout=devnull)
    except:
        raise Exception(
            "Brainsuite is not installed or bdp.sh is not in $PATH.")


def brainsuite_susceptibility_correction(outdir,
                                         dwi,
                                         bval,
                                         bvec,
                                         subject_id,
                                         phase_enc_dir,
                                         subjects_dir=None,
                                         qc_dir=None,
                                         bdp_nthread=4):
    """
    Assuming the beginning of the preprocessing was done with Connectomist
    up to Eddy current and motion correction, we now want to make susceptbility
    distortion correction using Brainsuite.

    Parameters
    ----------
    outdir: str
        Path to directory where to output.
    dwi: str
        Path to the preprocessed diffusion-weighted Nifti.
    bval: str
        Path to the bval file associated to the Nifti.
    bvec: str
        Path to the bvec file associated to the Nifti.
    subject_id: str
        Subject identifier used in Freesurfer.
    phase_enc_dir: str
        In plane phase encoding direction, "y", "y-", "x" or "x-".
    subjects_dir: str, default None
        If the Freesurfer $SUBJECTS_DIR environment variable is not set, or to
        bypass it, pass the path.
    qc_dir: str, default None
        Path to directory where to output snapshots for QC. By default in
        outdir.
    bdp_nthread: int, default 4
        Number of threads for bdp (see bdp.sh --thread flag)
    """

    # Freesurfer 'subjects_dir' has to be passed or set as environment variable
    subjects_dir = get_or_check_freesurfer_subjects_dir(subjects_dir)

    # Check in plane phase encoding direction
    # The sign of the phase encoding direction has no impact in this strategy
    # e.g. for "+y" or "-y" use "y"
    possible_phase_enc_dirs = {"y", "x"}
    if phase_enc_dir not in possible_phase_enc_dirs:
        raise ValueError("Bad argument 'phase_enc_dir': {}, should be in {}."
                         .format(phase_enc_dir, possible_phase_enc_dirs))

    # Set and check path to T1 brain-only volume from Freesurfer (mgz format)
    t1_brain_mgz = os.path.join(subjects_dir, subject_id, "mri/brain.mgz")
    if not os.path.isfile(t1_brain_mgz):
        raise Exception("Missing file: {}".format(t1_brain_mgz))

    # Convert Freesurfer T1 to Nifti with reorientation to RAS (to be in the
    # same orientation as the diffusion data)
    t1_brain_RAS_nii = os.path.join(outdir, "t1_brain.nii.gz")
    cmd = ["mri_convert", t1_brain_mgz, t1_brain_RAS_nii,
           "--out_orientation", "RAS"]
    run_freesurfer_cmd(cmd, subjects_dir=subjects_dir)

    # Run bfc (bias correction: required by BrainSuite)
    t1_bfc = os.path.join(outdir, "t1_brain.bfc.nii.gz")
    cmd = ["bfc", "-i", t1_brain_RAS_nii, "-o", t1_bfc]
    subprocess.check_call(cmd)

    # Extract brain from the nodif volume with FSL bet2
    nodif_brain = os.path.join(outdir, "nodif_brain.nii.gz")
    bet2(dwi, nodif_brain, f=0.25, m=False)

    # Run bdp.sh: registration + diffusion model
    cmd = ["bdp.sh", t1_bfc, "--nii", dwi, "--bval", bval, "--bvec", bvec,
           "--dwi-mask", nodif_brain, "--dir=%s" % phase_enc_dir,
           "--threads=%i" % bdp_nthread]
    subprocess.check_call(cmd)

    # Path to files of interest, created by BrainSuite bdp.sh
    dwi_wo_susceptibility = os.path.join(
        outdir, "t1_brain.dwi.RAS.correct.nii.gz")

    ###############
    # Quality check: create snapshots to visually assert the registration
    # quality

    if qc_dir is None:
        qc_dir = outdir

    # The snapshots won't contain all slices, half of them
    nb_slices_in_z = nibabel.load(nodif_brain).get_shape()[2]

    # Path to registered T1
    t1_to_dif = os.path.join(outdir, "t1_brain.D_coord.nii.gz")

    # First png: T1 registered in diffusion with nodif edges
    t1_with_nodif_edges_png = os.path.join(qc_dir, "t1_with_nodif_edges.png")
    plot_image(t1_to_dif,
               edge_file=nodif_brain,
               snap_file=t1_with_nodif_edges_png,
               name="T1 in diffusion + edges of nodif",
               cut_coords=nb_slices_in_z - 2)

    # Second png: nodif with edges of T1 registered in diffusion
    nodif_with_t1_edges_png = os.path.join(qc_dir, "nodif_with_t1_edges.png")
    plot_image(nodif_brain,
               edge_file=t1_to_dif,
               snap_file=nodif_with_t1_edges_png,
               name="nodif + edges of registered T1",
               cut_coords=nb_slices_in_z - 2)

    return dwi_wo_susceptibility, bval, bvec


def complete_preproc_wo_fieldmap(outdir,
                                 dwi,
                                 bval,
                                 bvec,
                                 phase_enc_dir,
                                 subject_id,
                                 subjects_dir=None,
                                 invertX=True,
                                 invertY=False,
                                 invertZ=False,
                                 delete_steps=False):
    """
    Function that runs all preprocessing steps using Connectomist but
    with BrainSuite for the correction of susceptibility distortions.

    Parameters
    ----------
    outdir: str
        Path to directory where all the preprocessing will be done.
    dwi: str
        path to input diffusion-weighted Nifti data.
    bval: str
        Path to the Nifti's associated .bval file.
    bvec: str
        Ppath to Nifti's associated .bval file.
    phase_enc_dir: str
        In plane phase encoding direction, "y", "y-", "x" or "x-".
    subject_id: str
        Subject identifier used in Freesurfer.
    subjects_dir: Str, default None
        If the Freesurfer $SUBJECTS_DIR environment variable is not set, or to
        bypass it, pass the path.
    invertX: bool, default True
        If True invert x-axis in diffusion model.
    invertY: bool, default False
        Same as invertX but for y-axis.
    invertZ: bool, default False
        Same as invertX but for z-axis.
    delete_steps: bool
        If True remove all intermediate files and directories at the end of
        preprocessing, to keep only selected files:
        - preprocessed Nifti + bval + bvec + outliers.py + nodif_brain*.nii.gz

    Returns
    -------
    outdir: str,
        Path to directory with the preprocessed files.

    <unit>
        <output name="preproc_dwi"  type="File"      />
        <output name="preproc_bval" type="File"      />
        <output name="preproc_bvec" type="File"      />

        <input name="outdir"        type="Directory" />
        <input name="dwi"           type="File"      />
        <input name="bval"          type="File"      />
        <input name="bvec"          type="File"      />
        <input name="phase_enc_dir" type="Str"       />
        <input name="subject_id"    type="Str"       />
        <input name="subjects_dir"  type="Str"       />
        <input name="invertX"       type="Bool"      />
        <input name="invertY"       type="Bool"      />
        <input name="invertZ"       type="Bool"      />
        <input name="delete_steps"  type="Bool"      />
    </unit>
    """

    # Step 0 - Initialization

    # Freesurfer 'subjects_dir' has to be passed or set as environment variable
    subjects_dir = get_or_check_freesurfer_subjects_dir(subjects_dir)

    # Raise an Exception if BrainsuiteCheck is not installed
    check_brainsuite_installation()

    # Create the preprocessing output directory if not existing
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Step 1 - Import files to Connectomist and choose q-space model
    raw_dwi_dir = os.path.join(outdir, "01-Import_and_qspace_model")
    # The manufacturer option is unused (use brainsuite for susceptibility
    # correction) but required
    dwi_data_import_and_qspace_sampling(raw_dwi_dir,
                                        dwi=dwi,
                                        bval=bval,
                                        bvec=bvec,
                                        manufacturer="Siemens",
                                        invertX=invertX,
                                        invertY=invertY,
                                        invertZ=invertZ,
                                        subject_id=subject_id)

    # Step 2 - Create a brain mask
    rough_mask_dir = os.path.join(outdir, "02-Rough_mask")
    dwi_rough_mask_extraction(rough_mask_dir, raw_dwi_dir)

    # Step 3 - Detect and correct outlying diffusion slices
    outliers_dir = os.path.join(outdir, "03-Outliers")
    dwi_outlier_detection(outliers_dir, raw_dwi_dir, rough_mask_dir)

    # Export outliers.py
    path_outliers_py = os.path.join(outliers_dir, "outliers.py")
    shutil.copy(path_outliers_py, outdir)

    # Step 4 - Eddy current and motion correction
    eddy_motion_dir = os.path.join(outdir, "04-Eddy_current_and_motion")
    dwi_eddy_current_and_motion_correction(eddy_motion_dir,
                                           raw_dwi_dir,
                                           rough_mask_dir,
                                           outliers_dir)

    # Step 5 - Convert Connectomist result to Nifti with bval/bvec
    dwi, bval, bvec = export_eddy_motion_results_to_nifti(eddy_motion_dir,
                                                          filename="dwi_ecc")

    # Step 6 - Susceptibility correction using BrainSuite
    brainsuite_dir = os.path.join(outdir, "05-Suceptibility_BrainSuite")
    if not os.path.isdir(brainsuite_dir):
        os.mkdir(brainsuite_dir)

    dwi_wo_susceptibility, bval, bvec = \
        brainsuite_susceptibility_correction(outdir=brainsuite_dir,
                                             dwi=dwi,
                                             bval=bval,
                                             bvec=bvec,
                                             phase_enc_dir=phase_enc_dir,
                                             subject_id=subject_id,
                                             subjects_dir=subjects_dir,
                                             qc_dir=outdir)

    # Step 7 - move corrected diffusion to outdir
    dwi_preproc = os.path.join(outdir, "dwi.nii.gz")
    bval_preproc = os.path.join(outdir, "dwi.bval")
    bvec_preproc = os.path.join(outdir, "dwi.bvec")
    shutil.copyfile(dwi_wo_susceptibility, dwi_preproc)
    shutil.copyfile(bval, bval_preproc)
    shutil.copyfile(bvec, bvec_preproc)

    # Step 8 - Create a T2 brain mask of preprocessed DWI data
    bet2_prefix = os.path.join(outdir, "nodif_brain")
    nodif_brain = os.path.join(bet2_prefix + ".nii.gz")
    nodif_brain_mask = os.path.join(bet2_prefix + "_mask.nii.gz")
    bet2(dwi, bet2_prefix, f=0.25, m=True)

    # Step 9 - clean intermediate directories if requested
    if delete_steps:
        intermediate_directories = [raw_dwi_dir,
                                    rough_mask_dir,
                                    outliers_dir,
                                    eddy_motion_dir,
                                    brainsuite_dir]
        for directory in intermediate_directories:
            shutil.rmtree(directory)

    return (outdir, dwi_preproc, bval_preproc, bvec_preproc, nodif_brain,
            nodif_brain_mask)


###############################################################################
# Code meant to run multiple call to complete_preproc_wo_fieldmap() in parallel

import logging
import traceback

from multiprocessing import Manager, Process

# Messages for communication between processes (parallel processing)
FLAG_STOP_PROCESS = "STOP_WORK"
FLAG_PROCESS_FINISHED = "PROCESS_HAS_FINISHED"


def parallel_worker(work_queue, result_queue):
    """ Function to make complete_preprocessing work in parallel processing.
    """
    while True:
        new_work = work_queue.get()
        if new_work == FLAG_STOP_PROCESS:
            result_queue.put(FLAG_PROCESS_FINISHED)
            break
        kwargs = new_work
        try:
            outputs = complete_preproc_wo_fieldmap(**kwargs)
            result_queue.put(outputs)
        except Exception as e:
            e.message += "\nPreprocessing failed for {}".format(new_work)
            e.message += "\n" + traceback.format_exc()
            result_queue.put(e.message)


def parallel_complete_preproc_wo_fieldmap(nb_procs,
                                          list_of_kwargs,
                                          log_dir=None):
    """
    Parameters
    ----------
    nb_procs: int
        Number of processes to run in parallel.
    list_of_kwargs: list of dicts
        For each dict (kwargs) will be run:
        complete_preprocessing_wo_fieldmap(**kwargs).
    log_dir: str
        Path to directory where to write the log file.
        If None, log file is written in current dir.
    """

    ###########################################################################
    # SETTING UP LOGGING SYSTEM
    ###########################################################################
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    # Check if file handler already exist, otherwise create it.
    has_file_handler = False
    for handler in logger.handlers:
        if type(handler) is logging.FileHandler:
            has_file_handler = True
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    if not has_file_handler:
        log_filename = "preproc_wo_fieldmap_%s" % time.strftime(
            "%Y-%m-%d_%H:%M:%S")

        if log_dir is None:
            path_log = log_filename
        else:
            path_log = os.path.join(log_dir, log_filename)

        file_handler = logging.FileHandler(path_log, mode="a")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        logger.info("Path to log file: %s" % path_log)

    ###########################################################################

    # Data structures for parallel processing
    manager = Manager()  # multiprocessing.Manager()
    work_queue, result_queue = manager.Queue(), manager.Queue()

    # Add work to the work_queue
    for kwargs in list_of_kwargs:
        work_queue.put(kwargs)

    # Add poison pills to stop the remote workers
    # When a process gets this job, it will stop
    for n in range(nb_procs):
        work_queue.put(FLAG_STOP_PROCESS)

    # Define processes
    workers = []
    for i in range(nb_procs):
        worker = Process(target=parallel_worker,
                         args=(work_queue, result_queue))
        worker.daemon = True
        workers.append(worker)
        worker.start()

    # Process results and log everything
    # try/except: to stop processes if user uses ctrl+c
    nb_finished_processes = 0
    try:
        while True:
            output = result_queue.get()
            if output == FLAG_PROCESS_FINISHED:
                nb_finished_processes += 1
                logger.info("Finished processes: %d/%d" % (
                    nb_finished_processes, nb_procs))
                if nb_finished_processes == nb_procs:
                    break
            elif isinstance(output, tuple):
                logger.info(
                    "Successful preprocessing, outdir: {}".format(output[0]))
            else:
                logger.warning("{}".format(output))

    except KeyboardInterrupt:  # To stop if user uses ctrl+c
        logger.info("KeyboardInterrupt: stopping processes.")
        for worker in workers:
            worker.terminate()
            worker.join()
