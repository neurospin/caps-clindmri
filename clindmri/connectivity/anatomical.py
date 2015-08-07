#! /usr/bin/env python
##########################################################################
# NSAP - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import numpy
import nibabel

# Clindmri import
from clindmri.tractography import Tractogram
from clindmri.extensions.fsl import flirt2aff


def diffusion_connectivity_matrix(track_file, label_file, outdir,
                                  symmetric=True):
    """ Counts the tracks that start and end at each label pair in the 
    diffusion space.

    Parameters
    ----------
    track_file: str (mandatory)
        a text file containing tracks.
    label_file: str (mandatory)
        a file containing labels that represent a segmentation of the cortex
        surface.
    outdir: str (mandatory)
        the output directory.
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
    label_image = nibabel.load(label_file)
    label_array = label_image.get_data()
    label_shape = label_image.get_shape()
    tractogram = Tractogram(track_file)

    # Check the validity of the label array
    kind = label_array.dtype.kind
    label_positive = (
        (kind == "u") or ((kind == "i") and (label_array.min() >= 0)))
    if not (label_positive and label_array.ndim == 3):
        raise ValueError("Label array must be a 3d integer array with "
                         "non-negative label values.")

    # To compute the connectivity matrix we consider only the first and last
    # point of each track
    endpoints = tractogram.endpoints().astype(int)
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

    # Compute the fiber density map
    density_map = tractogram.density(shape=label_shape)

    # Save the resulting connectivity matrix and density map
    proba_file = os.path.join(outdir, "det_paths.nii.gz")
    density_image = nibabel.Nifti1Image(density_map, label_image.get_affine())
    nibabel.save(density_image, proba_file)
    network_file = os.path.join(outdir, "det_network_matrix")
    numpy.savetxt(network_file, matrix)

    return proba_file, network_file


def anatomical_connectivity_matrix(track_file, label_file, t1_file,
                                   diffusion_file, trf_file,  outdir,
                                   symmetric=True):
    """ Counts the tracks that start and end at each label pair in the 
    anatomical space.

    Parameters
    ----------
    track_file: str (mandatory)
        a text file containing tracks.
    label_file: str (mandatory)
        a file containing labels that represent a segmentation of the cortex
        surface.
    t1_file: str (mandatory)
        a file containing the t1 image used in FreeSurfer for the segmentation.
    diffusion_file: str (optional, default None)
        a file containing the diffusion b0 3d image.
    trf_file: str (mandatory)
        a file with the FSL flirt transformation from the diffusion to the
        t1 spaces.
    outdir: str (mandatory)
        the output directory.
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
    label_image = nibabel.load(label_file)
    label_array = label_image.get_data()
    label_shape = label_image.get_shape()
    tractogram = Tractogram(track_file)
    affine = flirt2aff(trf_file, diffusion_file, t1_file)

    # Check the validity of the label array
    kind = label_array.dtype.kind
    label_positive = (
        (kind == "u") or ((kind == "i") and (label_array.min() >= 0)))
    if not (label_positive and label_array.ndim == 3):
        raise ValueError("Label array must be a 3d integer array with "
                         "non-negative label values.")

    # To compute the connectivity matrix we consider only the first and last
    # point of each track
    tractogram.apply_affine(affine)
    endpoints = tractogram.endpoints().astype(int)
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

    # Compute the fiber density map
    density_map = tractogram.density(shape=label_shape)

    # Save the resulting connectivity matrix and density map
    proba_file = os.path.join(outdir, "det_paths.nii.gz")
    density_image = nibabel.Nifti1Image(density_map, label_image.get_affine())
    nibabel.save(density_image, proba_file)
    network_file = os.path.join(outdir, "det_network_matrix")
    numpy.savetxt(network_file, matrix)

    return proba_file, network_file


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
