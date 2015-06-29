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

# Define factorial function
factorial = lambda n: reduce(lambda a, b: a * (b + 1), range(n), 1)


def construct_matrix_of_monomials(g, order):
    """ Construct all possible monomials in 3 variables of a given order.

    Parameters
    ----------
    g: array [N,3] (mandatory)
        gradient directions ( without the reference gradient (0,0,0) )
    order: odd int (mandatory)
        reconstruction order

    Returns
    -------
    A: array [N,(2+order)!/2(order)!]
        list of all possible monomials
    """

    N = g.shape[0]
    e = factorial(2 + order) / (2 * factorial(order))
    A = numpy.zeros((N, e), dtype=numpy.single)

    logging.debug("Constructing matrix of monomials...")
    for k in range(N):
        c = 0
        for i in range(order + 1):
            for j in range(order - i + 1):
                A[k, c] = ((g[k, 0] ** i) * (g[k, 1] ** j) *
                           (g[k, 2] ** (order - i - j)))
                c += 1
    return A


def generate_fourth_order_basis(g):
    """ From a given set of gradient g, generates the fourth order polynomial
    basis.
    """
    N = g.shape[0]
    A = numpy.zeros([N - 1, 15], dtype=numpy.single)

    A[:, 0] = g[1:, 2] ** 4
    A[:, 1] = g[1:, 2] ** 3 * g[1:, 1]
    A[:, 2] = g[1:, 2] ** 2 * g[1:, 1] ** 2
    A[:, 3] = g[1:, 1] ** 3 * g[1:, 2]
    A[:, 4] = g[1:, 1] ** 4
    A[:, 5] = g[1:, 2] ** 3 * g[1:, 0]
    A[:, 6] = g[1:, 2] ** 2 * g[1:, 0] * g[1:, 1]
    A[:, 7] = g[1:, 1] ** 2 * g[1:, 0] * g[1:, 2]
    A[:, 8] = g[1:, 1] ** 3 * g[1:, 0]
    A[:, 9] = g[1:, 0] ** 2 * g[1:, 2] ** 2
    A[:, 10] = g[1:, 0] ** 2 * g[1:, 1] * g[1:, 2]
    A[:, 11] = g[1:, 0] ** 2 * g[1:, 1] ** 2
    A[:, 12] = g[1:, 0] ** 3 * g[1:, 2]
    A[:, 13] = g[1:, 0] ** 3 * g[1:, 1]
    A[:, 14] = g[1:, 0] ** 4

    return A

if __name__ == '__main__':

    g = numpy.array([[0, 0, 0], [0.1639, 0.5115, 0.8435]])
    print g.shape

    print construct_matrix_of_monomials(g[1:], 4)
    print generate_fourth_order_basis(g)


