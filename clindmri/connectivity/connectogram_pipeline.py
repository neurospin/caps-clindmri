# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

import os
import time

from clindmri.connectivity.connectogram import (
    register_diffusion_to_anatomy,
    qc_dif2anat_registration,
    create_masks_for_tracto_seeding_wm,
    qc_tracto_masks,
    probtrackx2_connectogram_seeding_wm,
    qc_connectogram)


def connectogram_seeding_wm_pipeline(outdir,
                                     nodif_brain,
                                     nodif_brain_mask,
                                     bedpostx_dir,
                                     subject_id,
                                     cortical_atlas="Desikan",
                                     stop_mask_type="target_rois",
                                     subjects_dir=None,
                                     nsamples=5000,
                                     nsteps=2000,
                                     cthr=None,
                                     loopcheck=True,
                                     steplength=0.5,
                                     fibthresh=None,
                                     distthresh=None,
                                     sampvox=None,
                                     fsl_init="/etc/fsl/5.0/fsl.sh"):
    """
    Assuming you have a preprocessed DWI and that you have already run
    Bedpostx for this DWI, this function computes the rest of the steps to
    get a connectogram, including calls to the QC functions.
    The connectogram is computed using the --omatrix3 option in probtrackx2,
    seeding in white matter.

    Parameters
    ----------
    outdir: str
        Directory where to output.
    nodif_brain: str
        Path to the preprocessed brain-only DWI volume.
    nodif_brain_mask: str
        Path to the brain binary mask.
    bedpostx_dir: str
        Directory where Bedpostx has outputted.
    subject_id: str
        Subject id used with Freesurfer 'recon-all' command.
    cortical_atlas: {["Desikan"], "Destrieux"}
        The Freesurfer cortical atlas to use.
    stop_mask_type: {["target_rois"], "inverse_wm"}
        Strategy to use for stopping samples in probtrackx2:
        - "target_rois": stop samples as soon as they reach a target region
        - "inverse_wm": stop samples as soon as they leave the white matter
    subjects_dir: str or None, default None
        Path to the Freesurfer subjects directory. Required if the Freesurfer
        environment variable (i.e. $SUBJECTS_DIR) is not set.
    nsamples, nsteps, cthr, steplength: int, optional
        Probtrackx2 options.
    fibthresh, distthresh, sampvox: int, optional
        Probtrackx2 options.
    loopcheck: bool, optional
        Probtrackx2 option.
    fsl_init: str, optional.
        Path to the Bash script setting the FSL environment.
    """

    # STEP 1 - Compute the affine transformation between T1 and nodif_brain
    dif2anat_dat = register_diffusion_to_anatomy(outdir=outdir,
                                                 nodif_brain=nodif_brain,
                                                 subject_id=subject_id,
                                                 subjects_dir=subjects_dir,
                                                 subdir="diff_to_anat",
                                                 fsl_init=fsl_init)

    # STEP 2 - QC the affine registration
    qc_dif2anat_registration(outdir=outdir,
                             nodif_brain=nodif_brain,
                             dif2anat_dat=dif2anat_dat,
                             subject_id=subject_id,
                             subjects_dir=subjects_dir,
                             subdir="qc")

    # STEP 3 - Create the needed files to run probtrackx2: the masks
    txt_roi_masks, tracto_mask, wm_mask, stop_mask, avoid_mask = \
        create_masks_for_tracto_seeding_wm(outdir=outdir,
                                           nodif_brain=nodif_brain,
                                           nodif_brain_mask=nodif_brain_mask,
                                           dif2anat_dat=dif2anat_dat,
                                           subject_id=subject_id,
                                           cortical_atlas=cortical_atlas,
                                           stop_mask_type=stop_mask_type,
                                           subjects_dir=subjects_dir,
                                           subdir="masks")

    # STEP 4 - QC the tractography masks
    qc_tracto_masks(outdir,
                    nodif_brain,
                    tracto_mask,
                    wm_mask,
                    stop_mask,
                    txt_roi_masks,
                    subdir="qc")

    # STEP 5 - Run probtrackx2
    fiber_density, txt_matrix, txt_matrix_normalized, txt_labels = \
        probtrackx2_connectogram_seeding_wm(outdir=outdir,
                                            bedpostx_dir=bedpostx_dir,
                                            txt_roi_masks=txt_roi_masks,
                                            tracto_mask=tracto_mask,
                                            wm_mask=wm_mask,
                                            stop_mask=stop_mask,
                                            avoid_mask=avoid_mask,
                                            nsamples=nsamples,
                                            nsteps=nsteps,
                                            cthr=cthr,
                                            loopcheck=loopcheck,
                                            steplength=steplength,
                                            fibthresh=fibthresh,
                                            distthresh=distthresh,
                                            sampvox=sampvox,
                                            subdir="probtrackx2")

    # STEP 6 - QC the connectograms (default and normalized)
    qc_connectogram(outdir=outdir,
                    tracto_mask=tracto_mask,
                    fiber_density=fiber_density,
                    txt_matrix=txt_matrix,
                    txt_matrix_normalized=txt_matrix_normalized,
                    txt_labels=txt_labels,
                    subdir="qc")


###############################################################################
# To run multiple call to connectogram_seeding_wm_pipeline() in parallel

import logging
import traceback

from multiprocessing import Manager, Process

# Messages for communication between processes (parallel processing)
FLAG_STOP_PROCESS = "STOP_WORK"
FLAG_PROCESS_FINISHED = "PROCESS_HAS_FINISHED"


def parallel_worker(work_queue, result_queue):
    """ Each process executes this function.
    """
    while True:
        new_work = work_queue.get()
        if new_work == FLAG_STOP_PROCESS:
            result_queue.put(FLAG_PROCESS_FINISHED)
            break
        kwargs = new_work
        try:
            connectogram_seeding_wm_pipeline(**kwargs)
            result_queue.put("Successful run: %s" % kwargs["outdir"])
        except Exception as e:
            e.message += "\nProcessing failed for {}".format(new_work)
            e.message += "\n" + traceback.format_exc()
            result_queue.put(e.message)


def parallel_connectogram_seeding_wm_pipeline(nb_procs,
                                              list_of_kwargs,
                                              log_dir=None):
    """
    Parameters
    ----------
    nb_procs: int
        Number of processes to run in parallel.
    list_of_kwargs: list of dicts
        For each dict (kwargs) will be run:
        connectogram_seeding_wm_pipeline(**kwargs).
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
        log_filename = "connectogram_seeding_wm_pipeline_%s" % time.strftime(
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
