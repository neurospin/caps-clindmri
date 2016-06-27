##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import numpy as numpy
import nibabel
import warnings

# Clindmri import
from clindmri.estimation.utils import decompose_second_order_tensor
from clindmri.estimation.utils import dti6to33


def gdti_scalars(gdtifile, output_directory):
    """ Compute scalar maps from generalized diffusion tensor model data.

    <process>
        <return name="scalars" type="List_File" desc="The A dictionary with
            scalar map name as keys and sclar map files derived from the tensor
            model as values."/>
        <input name="gdtifile" type="File" desc="A file containing the
            generalized tensor model independent coefficients."/>
        <input name="output_directory" type="Directory" desc="The destination
            folder."/>
    </process>
    """
    # Load the tensor model
    gdti = nibabel.load(gdtifile)
    gdti_params = gdti.get_data()
    affine = gdti.get_affine()
    e = gdti_params.shape[-1]

    # Switch defined on the model order
    scalars = {}
    if e == 6:
        fa, md, cl, cp, cs = compute_second_order_scalar_parameters(
            gdti_params)
        scalars["fa"] = fa
        scalars["md"] = md
        scalars["cl"] = cl
        scalars["cp"] = cp
        scalars["cs"] = cs
    elif e == 15:
        ga, gmd = compute_fourth_order_scalar_parameters(
            gdti_params, k1=5000., k2=250.)
        scalars["ga"] = ga
        scalars["gmd"] = gmd
    else:
        warnings.warn("No scalar map has been computed for image "
                      "'{0}'.".format(gdtifile))

    # Save the computed scalar maps
    fname = os.path.basename(gdtifile).split(".")[0]
    for map_name, map_array in scalars.items():
        map_image = nibabel.Nifti1Image(map_array, affine=affine)
        mapfile = os.path.join(output_directory, "{0}_{1}.nii.gz".format(
            fname, map_name))
        map_image.to_filename(mapfile)
        scalars[map_name] = mapfile

    return scalars.values()


def compute_second_order_scalar_parameters(gdti6):
    """ Compute the Fractional Anisotropy [1]_ (FA),
    Mean Diffusivity (MD) and the Westion Shapes coefficients [1]_
    (cl, cp, cs) of a second order tensor.

    .. hidden-technical-block::
        :label: [+show/hide second order tensor scalars]
        :starthidden: True

        .. include:: source/_static/technical_documentation/dt2_scalars.txt

    **References**

    .. [1] C. Westin, S. Maier, H. Mamata, A. Nabavi, F. Jolesz and
           R. Kikinis : Processing and visualization of diffusion tensor
           MRI. Medical Image Analysis, 6(2):93-108, 2002.

    Parameters
    ----------
    dt33: numpy array (..., 3, 3)
        a second order symmetric tensor.

    Returns
    -------
    fa: array [...]
        fractional anisotropy map.
    md: array [...]
        mean diffusivity map.
    cl: array [...]
        linearity coefficient map.
    cp: array [...]
        planarity coefficient map.
    cs: array [...]
        sphericity coefficient map.
    """
    # Full tensor representation
    gdti33 = dti6to33(gdti6)

    # Eigenvalue decomposition
    eigenvalues, eigenvectors = decompose_second_order_tensor(gdti33)
    ev1 = eigenvalues[..., 0]
    ev2 = eigenvalues[..., 1]
    ev3 = eigenvalues[..., 2]

    # Compute denominator: avoid division by zero
    all_zero = (eigenvalues == 0).all(axis=-1)
    ev2 = ev1*ev1 + ev2*ev2 + ev3*ev3 + all_zero

    # Compute mean diffusivity map
    md = numpy.mean(eigenvalues, axis=3)

    # Compute fractional anisotropy map
    fa = numpy.sqrt(0.5 * ((ev1 - ev2) ** 2 + (ev2 - ev3) ** 2 +
                           (ev3 - ev1) ** 2) / ev2)

    # Compute lineraity coefficient map
    cl = (ev1 - ev2) / numpy.sqrt(ev2)

    # Compute planarity coefficient map
    cp = 2 * (ev2 - ev3) / numpy.sqrt(ev2)

    # Compute sphericity coefficient map
    cs = 3 * ev3 / numpy.sqrt(ev2)

    return fa, md, cl, cp, cs


def compute_fourth_order_scalar_parameters(gdti15, k1=5000., k2=250.):
    """ Compute the Generalized Anisotropy (GA) and Meand Diffusivity (MD)
    of fourth order tensors.

    Parameters
    ----------
    gti15: numpy array (..., 15)
        a fourth order symmetric tensor independent coefficients.
    k1, k2: float
        the state of the art GA regularization parameters.

    Returns
    -------
    ga: array [...]
        generalized fractional anisotropy.
    gmd: array [...]
        generalized mean diffusivity.

    References
        Ozarslan and Barmpoutnis
    """
    # Flattened the input data
    gdti15_flat = gdti15.reshape((-1, gdti15.shape[-1]))

    # Allocate output scalar maps
    gas = numpy.zeros((len(gdti15_flat), 1))
    mds = numpy.zeros((len(gdti15_flat), 1))

    # Go through all tensors
    for tensor, ga, md in zip(gdti15_flat, gas, mds):

        # Compute the mean diffusivity
        md[:] = 0.2 * (tensor[14] + tensor[4] + tensor[0] + 1 / 3 *
                       (tensor[11] + tensor[2] + tensor[9]))

        # Variance of the normalized diffusivities as a measure of anisotropy
        if md > 0:
            dist = hellinger(tensor, numpy.zeros((15, ), dtype=numpy.single))
            var = 1 / 9 * (dist / md ** 2 - 1)
            epsi = 1 + 1 / (1 + k1 * var)
            ga[:] = 1 - 1 / (1 + (k2 * var) ** epsi)

    # Get the associated scalar arrays
    ga = gas.reshape(gdti15.shape[:-1])
    gmd = mds.reshape(gdti15.shape[:-1])

    return ga, gmd


def hellinger(D1, D2):
    """ Compute a distance measure between 4th order diffusion tensor
    D1 and D2 by computing the normalized L2 distance between the
    corresponding diffusivity functuions d1(g) and d2(g).

    Parameters
    ----------
    D1: array [15,]
        4th order diffusion tensor independent coefficients.
    D2: array [15,]
        4th order diffusion tensor independent coefficients.

    Returns
    -------
    dist: float
        Hellinger distance.
    """
    diff = D1 - D2
    dist = 0.0

    a = 1 / 9
    b = 1 / 105
    c = 1 / 63
    d = 1 / 315

    dist += d * (diff[14] + diff[4] + diff[0] + diff[11] + diff[2] +
                 diff[9]) ** 2
    dist += (b - d) * (diff[14] + diff[4] + diff[0]) ** 2
    dist += (c - d) * ((diff[14] + diff[11]) ** 2 + (diff[14] + diff[9]) ** 2)
    dist += (c - d) * ((diff[4] + diff[11]) ** 2 + (diff[4] + diff[2]) ** 2)
    dist += (c - d) * ((diff[0] + diff[2]) ** 2 + (diff[0] + diff[9]) ** 2)
    dist += (a - b - 2 * (c - d)) * (diff[14] ** 2 + diff[4] ** 2 +
                                     diff[0] ** 2)
    dist += (b - d - 2 * (c - d)) * (diff[11] ** 2 + diff[2] ** 2 +
                                     diff[9] ** 2)
    dist += d * (diff[10] + diff[3] + diff[1]) ** 2
    dist += d * (diff[7] + diff[12] + diff[5]) ** 2
    dist += d * (diff[6] + diff[13] + diff[8]) ** 2
    dist += (b - d) * ((diff[13] + diff[8]) ** 2 + (diff[12] + diff[5]) ** 2 +
                       (diff[3] + diff[1]) ** 2)
    dist += (c - b) * (diff[13] ** 2 + diff[12] ** 2 + diff[8] ** 2 +
                       diff[3] ** 2 + diff[5] ** 2 + diff[1] ** 2)

    return dist


if __name__ == "__main__":

    output_directory = "/volatile/nsap/diffusion/estimation/"

    gdtifile = os.path.join(output_directory, "dwi_gdti.nii.gz")
    gdti_scalars(gdtifile, output_directory)

    # fourth order tensor: isotropic diffusion
    tensor = numpy.zeros((1, 1, 1, 15, ), dtype=numpy.single)
    d = 0.001
    tensor[..., :3] = d
    tensor[..., 3:6] = 2 * d
    print tensor
    ga, md = compute_fourth_order_scalar_parameters(tensor)
    print ga
    print md

    # second order tensor: isotropic diffusion
    tensor = numpy.zeros((1, 1, 1, 6, ), dtype=numpy.single)
    d = 0.001
    tensor[..., :3] = d
    print tensor
    fa, md, cl, cp, cs = compute_second_order_scalar_parameters(tensor)
    print fa
    print md
    print cl
    print cp
    print cs
