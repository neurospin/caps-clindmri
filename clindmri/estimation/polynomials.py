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
import numpy

# Caps import
from caps.diffusion_estimation.resources import get_sphere

# Define factorial function
factorial = lambda n: reduce(lambda a, b: a * (b + 1), range(n), 1)


def construct_set_of_polynomials(order):
    """ Construct the coefficients of homogenous polynomials in 3
    variables of a given order, which correspond to squares of lower
    (half) order polynomials.

    Parameters
    ----------
    order: odd int (mandatory)
        reconstruction order

    Returns
    -------
    C: array [M,(2+order)!/2(order)!]
        the computed list of polynomial coefficients
    """
    logging.debug("Constructing matrix of 321 polynomials...")

    # first load the unit vectors [321,3]
    g = numpy.loadtxt(get_sphere("symmetric321"))
    M = 321

    # Get the multiplicity of each independant element
    multiplicity = numpy.zeros((order + 1, order + 1, order + 1),
                               dtype=numpy.single)
    for i in range(order + 1):
        for j in range(order - i + 1):
            multiplicity[i, j, order - i - j] = get_coefficient_multiplicity(
                order, i, j, order - i - j)

    e = factorial(2 + order) / (2 * factorial(order))
    C = numpy.zeros((M, e), dtype=numpy.single)

    for k in range(M):
        c = 0
        for i in range(order + 1):
            for j in range(order - i + 1):
                C[k, c] = (multiplicity[i, j, order - i - j] *
                           (g[k, 0] ** i) * (g[k, 1] ** j) *
                           (g[k, 2] ** (order - i - j)))
                c += 1

    return C


def get_coefficient_multiplicity(order, i, j, k):
    """ Computes how many identical copies of a tensor coefficient appear
    in a fully symmetric high order tensor in 3 variables of spesific order.

    Parameters
    ----------
    order: odd int (mandatory)
        the order of the symmetric tensor in 3 variables (x,y,z)
    i,j,k: ints (mandatory)
        a triplet that represent the coefficient which is weighted by
        the monomial x^i*y^j*z^k, where i+j+k=order

    Returns
    -------
    counter	: int
        number of identical copies
    """
    size = 3 ** order
    counter = 0

    for z in range(size):
        ret = numpy.zeros((3, ), dtype=numpy.single)
        c = population_basis(z, order, ret)
        if ((c[0] == i) and (c[1] == j) and (c[2] == k)):
            counter += 1

    return counter


def population_basis(i, order, c):
    """ Generate basis coefficient
    """
    if (order == 0):
        ret = c
    else:
        j = i % 3
        c[j] += 1
        ret = population_basis((i - j) / 3, order - 1, c)

    return ret


def print_polynomial_basis(order):
    """ Print the polynomial basis

    Parameters
    ----------
    order: odd int (mandatory)
        the order of the symmetric tensor in 3 variables (x,y,z)
    """
    for i in range(order + 1):
        for j in range(order - i + 1):
            print ("Coefficient (x,y,z)=({0},{1},{2}) - Multiplicity: "
                  "{3}".format(i, j, order - i - j,
                  get_coefficient_multiplicity(order, i, j, order - i - j)))


if __name__ == '__main__':

    order = 2
    print_polynomial_basis(order)
    print construct_set_of_polynomials(order)
