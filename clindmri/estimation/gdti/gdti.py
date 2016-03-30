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
import numpy
import os
import scipy.optimize
import nibabel
import itertools


# Caps import
from monomials import construct_matrix_of_monomials
from polynomials import construct_set_of_polynomials
from integrals import construct_matrix_of_integrals
from caps.diffusion_estimation.resources import get_sphere
from diffusion_multiprocessing import multi_processing


def generalized_tensor_estimation(dfile, bvalfile, bvecfile, order,
                                  maskfile, number_of_workers,
                                  odf, output_directory):
    """ Fits a diffusion tensor given diffusion-weighted signals and
    gradient info using a quartic decomposition and non negative least square
    estimation strategy [1]_.

    **References**

    .. [1] A. Barmpoutis et B. Vermuri : A unified framework
           for estimating diffusion tensors of any order with
           symmetric positive-definite constraints.
           In IEEE Internaional Symposium on Biomedical
           Imaging, Rotterdam, ISBI, p1385-1388, 2010.

    <process>
        <return name="gdtifile" type="File" desc="The result file containing
            the tensor model independent coefficients."/>
        <input name="dfile" type="File" desc="The diffusion volume file."/>
        <input name="bvalfile" type="File" desc="The diffusion b-values
            file."/>
        <input name="bvecfile" type="File" desc="The diffusion b-vectors
            file."/>
        <input name="order" type="Int" desc="The diffusion tensor order (must
            be even)."/>
        <input name="maskfile" type="File" desc="A mask applied on the
            diffusion volume."/>
        <input name="number_of_workers" type="Int" desc="The number of CPUs
            used during the execution."/>
        <input name="odf" type="Bool" desc="If True estimate the odf."/>
        <input name="output_directory" type="Directory" desc="The destination
            folder."/>
    </process>
    """
    # Intern parameters
    # > Replace diffusion signal smaller than this threshold in order to
    #   avoid taking log(0) during the ols estimation
    min_signal = 1e-3
    # > Odf kernel
    delta = 100.

    # Convert maskfile parameter
    if not isinstance(maskfile, basestring):
        maskfile = None

    # Check that the model order is even
    if order % 2 == 1:
        raise Exception("'{0}' is not even. Expect an even model "
                        "order.".format(order))

    # Load the diffusion sequence
    (data_flat, mask_flat, ref_flat, bvals, bvecs,
     shape, affine) = load_diffusion_dataset(
        dfile, bvalfile, bvecfile, maskfile)

    # Estimate the generalized tensor model
    if odf:
        dti_parameters_flat = quartic_tensor_odf_fit(
            data_flat, mask_flat, ref_flat, bvals, bvecs, order, delta,
            min_signal, number_of_workers)
    else:
        dti_parameters_flat = quartic_tensor_fit(
            data_flat, mask_flat, ref_flat, bvals, bvecs, order,
            min_signal, number_of_workers)

    # Save the estimated tensor coefficients
    e = dti_parameters_flat.shape[-1]
    fname = os.path.basename(dfile).split(".")[0]
    dti_parameters = dti_parameters_flat.reshape(shape[:-1] + (e, ))
    dti_image = nibabel.Nifti1Image(dti_parameters, affine=affine)
    gdtifile = os.path.join(output_directory, "{0}_gdti.nii.gz".format(fname))
    dti_image.to_filename(gdtifile)

    return gdtifile


def load_diffusion_dataset(dfile, bvalfile, bvecfile, maskfile=None):
    """ Load the diffusion dataset.

    Parameters
    ----------
    dfile: string (mandatory)
        The diffusion volume file.
    bvalfile: string (mandatory)
        The diffusion b-values file.
    bvecfile: string (mandatory)
        The diffusion b-vectors file.
    maskfile:string (optional, default None)
        A mask applied on the diffusion volume.

    Returns
    -------
    data_flat: array
        flattened diffusion data array.
    mask_flat: array
        flattened mask data array.
    ref_flat: array
        flattened b0 reference array.
    bvals: array
        diffusion b-values.
    bvecs: array
        diffusion b-vectors.
    shape: array
        the diffusion volume shape.
    affine: array
        the diffusion volume affine transformation.
    """
    # Load diffusion sequence
    dwi = nibabel.load(dfile)
    bvals = numpy.loadtxt(bvalfile)
    bvecs = numpy.loadtxt(bvecfile)
    affine = dwi.get_affine()
    data = dwi.get_data()
    shape = data.shape

    # Check bvec parameter
    if bvecs.shape[1] != 3 and bvecs.shape[0] == 3:
        bvecs = bvecs.T

    # Detect b0 reference volume
    b0_indices = numpy.where(bvals == 0)[0]
    bvals = numpy.delete(bvals, b0_indices, axis=0)
    bvecs = numpy.delete(bvecs, b0_indices, axis=0)
    ref = numpy.mean(data[..., b0_indices], axis=3)
    data = numpy.delete(data, b0_indices, axis=3)

    # Create/load the mask
    if maskfile is None:
        mask = numpy.ones(shape[:-1], dtype=bool)
    else:
        mask = numpy.array(nibabel.load(maskfile).get_data(),
                           dtype=bool, copy=False)

    # Flatten the diffusion data
    mask_flat = mask.reshape((-1,))
    data_flat = data.reshape((-1, data.shape[-1]))
    ref_flat = ref.reshape((-1,))

    return data_flat, mask_flat, ref_flat, bvals, bvecs, shape, affine


def quartic_tensor_fit(data_flat, mask_flat, reference_flat, bvals, bvecs,
                       order, min_signal=1., number_of_workers=1):
    """ Fits a diffusion tensor given diffusion-weighted signals and
    gradient inforomation using a quartic decomposition and non negative least
    square estimation strategy.

    This procedure guarentees the positive-definite or at least positive
    semi-definite nature of tensors.

    Parameters
    ----------
    data_flat: array (mandatory)
        flattened diffusion data array.
    mask_flat: array (mandatory)
        flattened mask data array.
    reference_flat: array (mandatory)
        flattened reference b0 data array.
    bvals: array (mandatory)
        the diffusion b-values.
    bvecs: array (mandatory)
        the diffusion b-vectors.
    order: int (mandatory)
        the diffusion tensor order.
    min_signal: float (optional, default 1.)
        replace diffusion signal smaller than this threshold.
    number_of_workers: int (optional, default 1)
        the number of CPUs used during the execution.

    Returns
    -------
    dti_parameters: array [N, e]
        the tensor independant coefficients:
        multiplicity is included in tensor coefficients.

   Reference
       Barmpoutnis
    """
    # Construct b diag
    b_diag = numpy.diag(bvals)

    # construct basis
    G = construct_matrix_of_monomials(bvecs, order)
    C = construct_set_of_polynomials(order).T
    P = numpy.dot(G, C)
    P = numpy.dot(- b_diag, P)
    e = G.shape[1]

    # Allocate
    dti_parameters = numpy.empty((len(data_flat), e))

    # NNLS
    nb_cpus = None
    if number_of_workers > 0:
        nb_cpus = number_of_workers
    results = multi_processing(
        nnls_multi, nb_cpus,
        data_flat, reference_flat,
        itertools.repeat(P), itertools.repeat(C),
        itertools.repeat(min_signal), mask_flat,
        itertools.repeat(True))

    # Format result
    for cnt, item in enumerate(results):
        dti_parameters[cnt] = item

    return dti_parameters


def quartic_tensor_odf_fit(data_flat, mask_flat, reference_flat, bvals, bvecs,
                           order, delta=100., min_signal=1.,
                           number_of_workers=1):
    """ Compute a Cartesian tensor Odf from a given DW dataset.

    Parameters
    ----------
    data_flat: array (mandatory)
        flattened diffusion data array.
    mask_flat: array (mandatory)
        flattened mask data array.
    reference_flat: array (mandatory)
        flattened reference b0 data array.
    bvals: array (mandatory)
        the diffusion b-values.
    bvecs: array (mandatory)
        the diffusion b-vectors.
    order: int (mandatory)
        the diffusion tensor order.
    delta: float (optional, default 100.)
        the ODF kernel.
    min_signal: float (optional, default 1.)
        replace diffusion signal smaller than this threshold.
    number_of_workers: int (optional, default 1)
        the number of CPUs used during the execution.

    Returns
    -------
    dti_parameters: array [N, e]
        the odf tensor independant coefficients:
        multiplicity is included in tensor coefficients.

   Reference
       Barmpoutnis
    """
    # Contruct basis
    C = construct_set_of_polynomials(order).T
    BG = construct_matrix_of_integrals(bvecs, order, delta)
    B = numpy.dot(BG, C)
    e = C.shape[0]

    # Allocate
    dti_parameters = numpy.empty((len(data_flat), e))

    # NNLS
    nb_cpus = None
    if number_of_workers > 0:
        nb_cpus = number_of_workers
    results = multi_processing(
        nnls_multi, nb_cpus,
        data_flat, reference_flat,
        itertools.repeat(B), itertools.repeat(C),
        itertools.repeat(min_signal), mask_flat,
        itertools.repeat(False))

    # Format result
    for cnt, item in enumerate(results):
        dti_parameters[cnt] = item

    return dti_parameters


def nnls_multi(parameters):
    """ Dummy function to split input parameters.
    """
    return nnls_iter(*parameters)


def nnls_iter(dw_signal, ref_signal, P, C, min_signal, in_mask, take_log):
    """ Non negative least square tensor fit.

    Parameters
    ----------
    """
    if in_mask:
        # Throw out small signals
        dw_signal = numpy.maximum(dw_signal, min_signal)
        ref_signal = numpy.maximum(ref_signal, min_signal)
        # Signal
        y = dw_signal / ref_signal
        if take_log:
            y = numpy.log(y)
        # NNLS fit
        x = scipy.optimize.nnls(P, y)[0]
        # Tensor
        return numpy.dot(C, x)
    else:
        return numpy.zeros((C.shape[0], ))


def tensor2odf(tensors, b, order):
    """ Compute a Cartesian Tensor ODF from a given Higher-order diffusion
    tensor
    Parameters
    ----------
    tensors: list of array [e,] (mandatory)
        a list with tensor independent coefficients
    b: float (mandatory)
        the diffusion b-value

    Returns
    -------
    odfs list of array [e,]
        a list of tensor independant coefficients

    Reference
        Barmpoutnis
    """
    # Unit vectors [321,3]
    g = numpy.loadtxt(get_sphere("symmetric321"))

    # Construct basis
    G = construct_matrix_of_monomials(g, order)
    C = construct_set_of_polynomials(order).T
    BG = construct_matrix_of_integrals(g, order, 100)
    B = numpy.dot(BG, C)

    # loop over elements
    odfs = []
    for tensor in tensors:
        x = scipy.optimize.nnls(B, numpy.exp(-b * numpy.dot(G, tensor)))[0]
        odfs.append(numpy.dot(C, x))

    return odfs


def print_tensor(tensor, order):
    """ Print the tensor coefficients.

    Parameters
    ----------
    tensor: array [e,] (mandatory)
        the tensor independent coefficients.
    order: odd int (mandatory)
        the tensor order.
    """
    c = 0
    for i in range(order + 1):
        for j in range(order + 1 - i):
            print("D (x,y,z)=({0},{1},{2}) - Coefficient: "
                  "{3}".format(i, j, order - i - j, tensor[c]))
            c = c + 1

if __name__ == "__main__":

    from clindmri.plot import pvtk
    from clindmri.estimation.utils import dti6to33

    output_directory = "/volatile/nsap/diffusion/estimation/"

    order = 2
    S0 = 1
    factor = 1.
    g = numpy.array([
        [0, 0, 0],
        [0, 0, 0],
        [0.1639, 0.5115, 0.8435],
        [0.1176, -0.5388, 0.8342],
        [0.5554, 0.8278, -0.0797],
        [-0.4804, 0.8719, 0.0948],
        [0.9251, -0.0442, 0.3772],
        [0.7512, -0.0273, -0.6596],
        [0.1655, -0.0161, 0.9861],
        [0.6129, -0.3427, 0.7120],
        [0.6401, 0.2747, 0.7175],
        [-0.3724, -0.3007, 0.8780],
        [-0.3451, 0.3167, 0.8835],
        [0.4228, 0.7872, 0.4489],
        [0.0441, 0.9990, 0.0089],
        [-0.1860, 0.8131, 0.5515],
        [0.8702, 0.4606, 0.1748],
        [-0.7239, 0.5285, 0.4434],
        [-0.2574, -0.8032, 0.5372],
        [0.3515, -0.8292, 0.4346],
        [-0.7680, -0.4705, 0.4346],
        [0.8261, -0.5384, 0.1660],
        [0.9852, -0.0420, -0.1660]])
    b = numpy.asarray([0, 0] + [1500] * (g.shape[0] - 2))
    b.shape += (1, )

    sdwi = []
    dwi = numpy.zeros((3, 3, 6, g.shape[0]), dtype=numpy.single)

    # fiber1 = (pi/2,100*pi/180) - fiber2 = (pi/2,20*pi/180)
    sdwi.append(numpy.array([
        0.39481, 0.43774, 0.12879, 0.31532, 0.31744,
        0.36900, 0.59490, 0.35280, 0.36880, 0.44046, 0.48088,
        0.17118, 0.22700, 0.34665, 0.26000, 0.25414, 0.21642,
        0.34456, 0.26625, 0.20723, 0.30364]) * factor)

    # fiber1 = (pi/2,0) - fiber2 = (pi/2,pi/2)
    sdwi.append(numpy.array([0.43670, 0.43596, 0.18206, 0.20415, 0.33573,
                            0.37296, 0.59236, 0.33340, 0.34970, 0.45040,
                            0.45559, 0.24635, 0.32270, 0.33154, 0.21289,
                            0.21758, 0.30981, 0.26820, 0.23028, 0.18931,
                            0.32519]) * factor)

    # fiber1 = (pi/2,0) - fiber2 = (pi/4,pi/2)
    sdwi.append(numpy.array([
        0.30972, 0.56842, 0.27685, 0.24966, 0.29334,
        0.22803, 0.37018, 0.36464, 0.18036, 0.40585,
        0.26687, 0.22845, 0.38253, 0.30230, 0.21773,
        0.15808, 0.53451, 0.46730, 0.36886, 0.30290,
        0.30921]) * factor)

    # fiber1 = (pi/2,pi/8) - fiber2 = (pi/2,3*pi/8)
    sdwi.append(numpy.array([0.36248, 0.47541, 0.06924, 0.40454, 0.28182,
                            0.33958, 0.59539, 0.48527, 0.23793, 0.36287,
                            0.56245, 0.11984, 0.21594, 0.37846, 0.08456,
                            0.44135, 0.18359, 0.41061, 0.10898, 0.41320,
                            0.26074]) * factor)

    # one fiber (1,0,0)
    sdwi.append(numpy.array([0.57138, 0.59235, 0.26610, 0.32803, 0.06001,
                            0.13251, 0.57054, 0.22062, 0.20057, 0.42144,
                            0.44422, 0.37747, 0.61147, 0.55894, 0.08094,
                            0.14717, 0.51133, 0.43896, 0.12416, 0.09872,
                            0.03862]) * factor)

    # one fiber (0,1,0)
    sdwi.append(numpy.array([0.30201, 0.27957, 0.09802, 0.08027, 0.61146,
                            0.61340, 0.61418, 0.44618, 0.49883, 0.47935,
                            0.46697, 0.11523, 0.03394, 0.10413, 0.34484,
                            0.28799, 0.10829, 0.09744, 0.33641, 0.27990,
                            0.61176]) * factor)

    # construct dwi
    dwi[..., 0] = S0 * factor
    dwi[..., 1] = S0 * factor
    for index in range(len(sdwi)):
        dwi[:, :, index, 2:] = sdwi[index]
    dwi_image = nibabel.Nifti1Image(dwi, numpy.eye(4))
    dfile = os.path.join(output_directory, "dwi.nii.gz")
    bvecfile = os.path.join(output_directory, "dwi.bvec")
    bvalfile = os.path.join(output_directory, "dwi.bval")
    dwi_image.to_filename(dfile)
    numpy.savetxt(bvecfile, g)
    numpy.savetxt(bvalfile, b.T)

    gdtifile = generalized_tensor_estimation(
        dfile, bvalfile, bvecfile, order, maskfile=None, number_of_workers=1,
        odf=False, output_directory=output_directory)
    dti_params = nibabel.load(gdtifile).get_data()
    print "Tensor:"
    for index, params in enumerate(dti_params[0, 0]):
        print "nnls:", index
        print_tensor(params, order)
    for tensor in dti6to33(dti_params)[0, 0]:
        vals, vecs = numpy.linalg.eigh(tensor)
        print vals
        print vecs

    # Plot
    ren = pvtk.ren()
    for index, params in enumerate(dti_params[0, 0]):
        actor = pvtk.tensor(params, order, position=(index, 0, 0))
        pvtk.add(ren, actor)
    pvtk.show(ren)
