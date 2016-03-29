##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import numpy


def decompose_second_order_tensor(dt33):
    """ Returns eigenvalues and eigenvectors of a given diffusion tensor

    Negative eigenvalues are set to an arbitrary small value.

    Parameters
    ----------
    dt33: numpy array (..., 3, 3)
        a second order symmetric tensor.

    Returns
    -------
    eigvals: array [..., 3]
        the eigenvalues from the tensor eigen decomposition sorted from
        the largest to smallest.
    eigvecs: array [..., 3, 3]
        associated eigenvectors
    """
    # Flattened the tensor array
    dt33_flat = dt33.reshape((-1, dt33.shape[-2], dt33.shape[-1]))

    # Eigen decomposition
    eigenvals_flat, eigenvecs_flat = numpy.linalg.eigh(dt33_flat)

    # Need to sort the eigenvalues and associated eigenvectors:
    # note that each direction is a column
    order = eigenvals_flat.argsort()[..., ::-1]
    indices = numpy.indices(eigenvals_flat.shape)
    eigenvals_flat = eigenvals_flat[indices[0], order]
    indices = numpy.indices(eigenvecs_flat.shape)
    order = order.repeat(3, axis=0).reshape(eigenvecs_flat.shape)
    eigenvecs_flat = eigenvecs_flat[indices[0], indices[1], order]

    # Eigenvalues are positive definite
    eigenvals_flat = eigenvals_flat.clip(min=1e-8)

    # Create arrays
    eigenvals = eigenvals_flat.reshape(dt33.shape[:-1])
    eigenvecs = eigenvecs_flat.reshape(dt33.shape)

    return eigenvals, eigenvecs


def dti6to33(dt6):
    """ Full second order symmetric tensor from the six independent components.

    Parameters
    ----------
    dt6: array (..., 6)
        a second order symmetric tensor represented only with the six unique
        components.

    Returns
    -------
    dt33: array (..., 3, 3)
        a second order symmetric tensor.
    """
    dt33 = numpy.zeros(dt6.shape[:-1] + (3, 3), dtype=dt6.dtype)
    dt33[..., 0, 0] = dt6[..., 5]  # dxx
    dt33[..., 0, 1] = dt33[..., 1, 0] = dt6[..., 4] / 2  # dxy
    dt33[..., 0, 2] = dt33[..., 2, 0] = dt6[..., 3] / 2  # dxz
    dt33[..., 1, 1] = dt6[..., 2]  # dyy
    dt33[..., 1, 2] = dt33[..., 2, 1] = dt6[..., 1] / 2  # dyz
    dt33[..., 2, 2] = dt6[..., 0]  # dzz

    return dt33


def dti33to6(dt33):
    """ Six independent components from the full second order symmetric tensor.

    Parameters
    ----------
    dt33: numpy array (...,3,3)
        a second order symmetric tensor.

    Returns
    -------
    dt6: numpy array (...,6)
        a second order symmetric tensor represented only with the six unique
        components.
    """
    dt6 = numpy.zeros(dt33.shape[:-2] + (6, ), dtype=dt33.dtype)
    dt6[..., 5] = dt33[..., 0, 0]
    dt6[..., 4] = dt33[..., 0, 1] + dt33[..., 1, 0]
    dt6[..., 3] = dt33[..., 0, 2] + dt33[..., 2, 0]
    dt6[..., 2] = dt33[..., 1, 1]
    dt6[..., 1] = dt33[..., 1, 2] + dt33[..., 2, 1]
    dt6[..., 0] = dt33[..., 2, 2]

    return dt6


def eigen_compose(eigenvals, eigenvecs):
    """ Recover a second order tensor image [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz]
    """
    tensor = numpy.zeros(eigenvals.shape[:-1] + (6,))
    tensor[..., 5] = numpy.sum(
        eigenvecs[..., 0, i] * eigenvals[..., i] * eigenvecs[..., 0, i]
        for i in range(3))
    tensor[..., 2] = numpy.sum(
        eigenvecs[..., 1, i] * eigenvals[..., i] * eigenvecs[..., 1, i]
        for i in range(3))
    tensor[..., 0] = numpy.sum(
        eigenvecs[..., 2, i] * eigenvals[..., i] * eigenvecs[..., 2, i]
        for i in range(3))
    tensor[..., 4] = numpy.sum(
        eigenvecs[..., 0, i] * eigenvals[..., i] * eigenvecs[..., 1, i]
        for i in range(3))
    tensor[..., 3] = numpy.sum(
        eigenvecs[..., 0, i] * eigenvals[..., i] * eigenvecs[..., 2, i]
        for i in range(3))
    tensor[..., 1] = numpy.sum(
        eigenvecs[..., 1, i] * eigenvals[..., i] * eigenvecs[..., 2, i]
        for i in range(3))

    return tensor


if __name__ == "__main__":

    import nibabel
    import os

    output_directory = "/volatile/nsap/diffusion/estimation/"

    dt6_im = nibabel.load(os.path.join(output_directory, "dwi_gdti.nii.gz"))
    dt6 = dt6_im.get_data()
    affine = dt6_im.get_affine()

    dt33 = dti6to33(dt6)

    eigvals, eigvecs = decompose_second_order_tensor(dt33)
    recompose_dt6 = eigen_compose(eigvals, eigvecs)

    for index in range(dt6.shape[2]):
        print "---------"
        print "DT6 REF: ", dt6[0, 0, index]
        print "VALS: ", eigvals[0, 0, index]
        print "VECS: ", eigvecs[0, 0, index]
        print "DT6: ", recompose_dt6[0, 0, index]
