#! /usr/bin/env python
##########################################################################
# CAPS - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
from __future__ import division
import logging
import numpy as numpy

# Caps import
from caps.diffusion_estimation.resources import get_sphere

# Define factorial function
factorial = lambda n: reduce(lambda a, b: a * (b + 1), range(n), 1)


def construct_matrix_of_integrals(g, order, delta):

    logging.debug("Constructing matrix of Basis-Monomial Integrals...")

    # Unit vectors [321,3]
    g2 = numpy.loadtxt(get_sphere("symmetric321"))

    # Construct matrix of integrals
    N = g.shape[0]
    Bg = numpy.zeros((N, factorial(2 + order) / (2 * factorial(order))),
                     dtype=numpy.single)
    for k in range(N):
        c = 0
        for i in range(order + 1):
            for j in range(order + 1 - i):
                Bg[k, c] = 0
                for integr in range(321):
                    Bg[k, c] = (Bg[k, c] + numpy.exp(- delta * numpy.dot(
                         g2[integr, :], g[k, :]) ** 2) *
                         (g2[integr, 0] ** i) * (g2[integr, 1] ** j) *
                         (g2[integr, 2] ** (order - i - j)))
                Bg[k, c] = Bg[k, c] / 321
                c = c + 1
    return Bg