#! /usr/bin/env python
##########################################################################
# NSAP - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import numpy
import nibabel

# Clindmri import
from clindmri.tractography import loadtxt


def connectivity_matrix(track_file, label_file, symmetric=True):
    """ Counts the tracks that start and end at each label pair.

    Parameters
    ----------
    track_file: str (mandatory)
        a text file containing tracks.
    label_file: str (mandatory)
        a file containing labels that represent a segmentation of the cortex
        surface.
    symmetric: bool (optional, default True)
        symmetric means we don't distinguish between start and end points. If
        symmetric is True, 'matrix[i, j] == matrix[j, i]'.

    Returns
    -------
    matrix: array
        the number of connection between each pair of regions defined in the
        'label_file'.
    """
    # Load the dataset
    label_array = nibabel.load(label_file).get_data()
    tracks = loadtxt(track_file)

    # Check the validity of the label array
    kind = label_array.dtype.kind
    label_positive = (
        (kind == "u") or ((kind == "i") and (label_array.min() >= 0)))
    if not (label_positive and label_array.ndim == 3):
        raise ValueError("Label array must be a 3d integer array with "
                         "non-negative label values.")

    # To compute the connectivity matrix we consider only the first and last
    # point of each track
    endpoints = [line[0::len(line)-1] for line in tracks]
    endpoints = numpy.asarray(endpoints).astype(int)
    pointsx, pointsy, pointsz = endpoints.T

    # Get labels associted to track end points
    endlabels = label_array[pointsx, pointsy, pointsz]
    if symmetric:
        endlabels.sort(axis=0)
    matrix = ndbincount(endlabels)
    if symmetric:
        matrix = numpy.maximum(matrix, matrix.T)

    # Remove the connectivity associated to the background
    matrix = matrix[1:, 1:]

    return matrix


def ndbincount(x, weights=None, shape=None):
    """ Count number of occurrences of each value in array of non-negative ints
    for nd-indicies.

    Parameters
    ----------
    x: array (N, M)
        N input 1-dim array of size M
    weights : array (M,)
        weights associated with 1-dim array
    shape: 2-uplet (optional, default None)
        the shape of the output
    """
    if shape is None:
        shape = x.max(axis=1) + 1
    x = numpy.ravel_multi_index(x, shape)
    out = numpy.bincount(x, weights)
    out.resize(shape)

    return out
