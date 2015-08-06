#! /usr/bin/env python
##########################################################################
# CAPS - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import multiprocessing
import itertools


def multi_processing(function_multi, nb_workers, *parameters):
    """ Multiprocessing for diffusion estimation processes

    Parameters
    ----------
    function: @callable
        a function we want to optimize
    nb_workers: int
        the number of worker processes to use.
        If processes is None then all the cpus are used.
        If initializer is not None then *nb_workers* workers are used.
    parameters: list
        list with function parameters
    """
    results = []
    pool = multiprocessing.Pool(processes=nb_workers)
    r = pool.map_async(function_multi, itertools.izip(*parameters),
                       callback=results.append)
    r.wait()
    return results[0]
