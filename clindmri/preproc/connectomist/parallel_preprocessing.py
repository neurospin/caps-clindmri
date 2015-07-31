#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" Module that allows multiple complete_preprocessing() calls to be run in
    parallel.
"""

from multiprocessing import Process, Manager
import logging
import tempfile
import traceback

from .exceptions             import ConnectomistError
from .complete_preprocessing import complete_preprocessing

###############################################################################
# MODULE VARIABLES
###############################################################################

# Messages for communication between processes
FLAG_STOP_PROCESS     = "STOP_WORK"
FLAG_PROCESS_FINISHED = "PROCESS_HAS_FINISHED"

###############################################################################


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
            path_nifti, path_bval, path_bvec = complete_preprocessing(**kwargs)
            result_queue.put((path_nifti, path_bval, path_bvec))
        except ConnectomistError as e:
            result_queue.put(e.message)
        except Exception as e:
            e.message += "\nUnknown error happened for %s" % kwargs["path_nifti"]
            e.message += "\n" + traceback.format_exc()
            result_queue.put(e.message)


# TO BE COMPLETED: capsul + test
def parallel_preprocessing(nb_processes, list_kwargs, log_path=None):
    """
    Function to make complete_preprocessing() run in parallel.

    Parameters
    ----------
    list_kwargs: List of dicts, each dict stores the arguments for one call to
            complete_preprocessing(), i.e. one job.
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
        if not log_path:
            log_path = tempfile.mkstemp(dir="/volatile/logs", suffix=".log",
                                        prefix="preprocessing_")[1]
        file_handler = logging.FileHandler(log_path, mode="a")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        logger.info("Path to log file: %s" % log_path)

    ###########################################################################

    # Data structures for parallel processing
    manager = Manager()  # multiprocessing.Manager()
    work_queue, result_queue = manager.Queue(), manager.Queue()

    # Add jobs in work_queue
    for kwargs in list_kwargs:
        work_queue.put(kwargs)

    # Add poison pills to stop the remote workers
    # When a process gets this job, it will stop
    for n in range(nb_processes):
        work_queue.put(FLAG_STOP_PROCESS)

    # Define processes
    workers = []
    for i in range(nb_processes):
        worker = Process(target=parallel_worker, args=(work_queue, result_queue))
        worker.daemon = True
        workers.append(worker)
        worker.start()

    # Process results and log everything
    nb_finished_processes = 0
    try:
        while True:
            new_result = result_queue.get()
            if new_result == FLAG_PROCESS_FINISHED:
                nb_finished_processes += 1
                logger.info("Finished processes: %d/%d" % (nb_finished_processes,
                                                           nb_processes))
                if nb_finished_processes == nb_processes:
                    break
            elif type(new_result) is str:
                logger.warning(new_result)
            elif not (type(new_result) is tuple and len(new_result) == 3):
                logger.warning("Something went wrong; result: {}".format(new_result))
            else:
                logger.info("Successfull preprocessing, resulting files:"
                            "\n%s\n%s\n%s" % new_result)
    except KeyboardInterrupt:  # To stop if user uses ctrl+c
        logger.info("KeyboardInterrupt: stopping processes.")
        for worker in workers:
            worker.terminate()
            worker.join()
