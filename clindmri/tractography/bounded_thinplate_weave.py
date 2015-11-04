#! /usr/bin/env python
##########################################################################
# NSAP - Copyright (C) CEA, 2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import numpy
from scipy.weave import ext_tools


def build_bounded_thinplate_module():
    """ Creates an extension module named 'bounded_thinplate' and
    compiles it to a shared library (.so or .pyd) that can be loaded into
    Python.
    """
    # Creates an ext_module instance that is ready to have ext_function
    # instances added to it.
    mod = ext_tools.ext_module("bounded_thinplate")

    # Define ext_function
    # Calling 'pymc' conventions for new covariance:
    #   1- C: A Fortran-contiguous (column  major) array of dimension
    #      (x.shape[0], y.shape[0]). This should be overwritten in-place.
    #   2- x, y: Arrays of input locations. These will be regularized: they 
    #      will be two-dimensional, with the first index iterating over points
    #      and the second over coordinates.
    #   3- cmin=0, cmax=-1: Optional arguments. If non-default values are 
    #      provided, only the slice C[:,cmin:cmax] of C should be overwritten.
    #   4- symm=False: An optional argument indicating whether x and y are the
    #      same array. If True, C will be square and only the upper triangle of 
    #      C should be overwritten.
    #   5- R: an extra argument needed for the fiber GP computation.
    thinplate3d_code=r"""
    //thinplate3d(C, R, cmin, cmax, symm)
    if (cmax == -1) {
        cmax = Nc[1];
    }

    if (symm) {
        for(int j=cmin; j<cmax; j++) {
            C2(j, j) = 1.;
            for(int i=0; i<j; i++) {
                C2(i, j) = thinPlate3D(C2(i, j), R);
            }
        }
    }
    else {
        for(int j=cmin; j<cmax; j++) {
            for(int i=0; i<Nc[0]; i++) {
                C2(i, j) = thinPlate3D(C2(i, j), R);
            }
        }
    }
    """
    innerproduct_thinplate3d_code = r"""
    //innerProduct_thinplate3d(C, Q, R, symm)
    if (R < Q) {
      const double aux = Q;
      Q = R;
      R = aux;
    }

    double Rs[11];
    double Qs[11];
    Rs[0] = 1;
    Qs[0] = 1;
    for(int i=0; i<10; i++ ) {
       Rs[i + 1] = Rs[i] * R;
       Qs[i + 1] = Qs[i] * Q;
    }

    if (symm) {
        for (int j=0; j<Nc[1]; j++) {
            for (int i=0; i<=j; i++){
               C2(i, j) = covIntegralThinPlateR3Normalized(C2(i, j), Q, R,
                                                           Qs, Rs);
               C2(j, i) = C2(i, j);
            }
      }
    }
    else { 
        for(int j=0; j<Nc[1]; j++) {
            for(int i=0; i<Nc[0]; i++) {
                C2(i, j) = covIntegralThinPlateR3Normalized(C2(i, j), Q, R,
                                                            Qs, Rs);
            }
        }
    }
    """

    # Define support_code
    thinplate_support_code = r"""
    double thinPlate3D(double d, double R) {
        double t;
        if (d == 0.) {
            t = 1.0e0;
        }
        else {
            if (d < R) {
                const double d2 = d * d;
                const double d3 = fabs(d) * d2;
                const double R2 = R * R;
                const double R3 = fabs(R) * R2;
                t = (2. * d3 - 3. * R * d2 + R3) / R3;
            }
            else {
                t = 0.0e0; 
            }
        }
        return t;
    }
    """
    innerproduct_thinplate3d_normalized_support_code = r"""
    double covIntTPR3NwLTRminusQ(
        double w, double Q, double R, const double* Qs, const double* Rs);
    double covIntTP3NRminusQLTw(
        double w, double Q, double R, const double* Qs, const double* Rs);
    double covIntTP3NwGTR(
        double w, double Q, double R, const double* Qs, const double* Rs);

    double covIntegralThinPlateR3Normalized(double w, double Q, double R,
                                            const double* Qs, const double* Rs) {
        if (w == 0) {
            return acos(-1.0e0) * double(Qs[3]) * double(84 * Rs[3] - 81 * R *
                   Qs[2] + 35 * Qs[3]) / double(Rs[3]) / 0.315e3;
        }
        else if (w <= (R - Q)) {
            return covIntTPR3NwLTRminusQ(w, Q, R, Qs, Rs);
        }
        else if (((R - Q) < w) && (w <= R)) {
            return covIntTP3NRminusQLTw(w, Q, R , Qs, Rs);
        }
        else if (w > R) {
            return covIntTP3NwGTR(w, Q, R, Qs, Rs);
        }
        return 0;
    }

    double covIntTPR3NwLTRminusQ(double w, double Q, double R,
                                 const double* Qs, const double* Rs) {
        if (w <= Q) {
            return acos(-1.0e0) * double(-405 * Qs[8] * R - 1260 * Qs[6] * R *
                   pow(w, 2) + 420 * Qs[6] * Rs[3] + 15 * Q * pow(w, 8) +
                   175 * Qs[9] + 900 * Qs[7] * pow(w, 2) + 378 * Qs[5] *
                   pow(w, 4) - 60 * Qs[3] * pow(w, 6) - 4 * pow(w, 9)) /
                   double(Qs[3]) / double(Rs[3]) / 0.1575e4;
        }
        else { 
            return acos(-1.0e0) * double(Qs[3]) * double(-135 * Qs[2] * w *
                   R - 420 * pow(w, 3) * R + 140 * w * Rs[3] + 180 * Qs[2] *
                   pow(w, 2) + 280 * pow(w, 4) + 8 * Qs[4]) / double(w) /
                   double(Rs[3]) / 0.525e3;
        }
    }

    double covIntTP3NRminusQLTw(double w, double Q, double R,
                                const double* Qs, const double* Rs) {
        if ( w <= Q ) { 
            return  -0.40e1 / 0.525e3 * M_PI * (pow(w, 10) / 0.2e1 + (-0.5e1 /
                    0.8e1 * R - 0.5e1 / 0.8e1 * Q) * pow(w, 9) - 0.135e3 /
                    0.64e2 * pow(w, 8) * R * Q + (0.5e1 / 0.2e1 * Rs[3] +
                    0.5e1 / 0.2e1 * Qs[3]) * pow(w, 7) + 0.105e3 / 0.16e2 *
                    Q * R * (Rs[2] + Qs[2]) * pow(w, 6) + (-0.63e2 / 0.4e1 *
                    Qs[5] - 0.63e2 / 0.4e1 * Rs[5]) * pow(w, 5) + (-0.945e3 /
                    0.32e2 * Q * Rs[5] - 0.175e3 / 0.16e2 * Qs[3] * Rs[3] -
                    0.945e3 / 0.32e2 * Qs[5] * R + 0.35e2 * Rs[6] + 0.35e2 *
                    Qs[6]) * pow(w, 4) + (0.105e3 / 0.2e1 * Rs[6] * Q -
                    0.75e2 / 0.2e1 * Rs[7] + 0.105e3 / 0.2e1 * Qs[6] * R -
                    0.75e2 / 0.2e1 * Qs[7]) * pow(w, 3) + 0.45e2 / 0.2e1 *
                    pow((Q - R), 4) * (Qs[4] + 0.17e2 / 0.8e1 * Qs[3] * R +
                    0.5e1 / 0.2e1 * Qs[2] *Rs[2] + 0.17e2 / 0.8e1 * Q *
                    Rs[3] + Rs[4]) * pow(w, 2) + (0.135e3 / 0.8e1 * Rs[8] *
                    Q - 0.175e3 / 0.24e2 * Qs[9] - 0.35e2 / 0.2e1* Qs[6] *
                    Rs[3] - 0.175e3 / 0.24e2 * Rs[9] - 0.35e2 / 0.2e1 *Rs[6] *
                    Qs[3] + 0.135e3 / 0.8e1 * Qs[8] * R) * w +
                    pow((Q - R), 6) * (Qs[4] + 0.209e3 / 0.64e2 * Qs[3] *
                    R + 0.147e3 / 0.32e2 * Qs[2] * Rs[2] + 0.209e3 / 0.64e2 *
                    Q * Rs[3] + Rs[4])) / Qs[3] / w / Rs[3];
        }
        else{
            return  M_PI * (0.192e3 * Qs[10] - 0.192e3 * Rs[10] -
                    0.32e2 * pow(w, 10) + 0.8100e4 * pow(w, 2) * Qs[7] * R -
                    0.10080e5 * pow(w, 3) * Qs[6] * R + 0.8100e4 * Q *
                    pow(w, 2) * Rs[7] - 0.3780e4 * pow(w, 2) * Qs[5] * Rs[3] +
                    0.3360e4 * Rs[6] * Qs[3] * w - 0.3240e4 * w * Qs[8] *
                    R - 0.900e3 * Qs[7] * Rs[3] + 0.480e3 * pow(w, 7) *
                    Qs[3] + 0.405e3 * pow(w, 8) * R * Q - 0.3780e4 * Qs[3] *
                    pow(w, 2) * Rs[5] - 0.1260e4 * pow(w, 6) * Qs[3] * R +
                    0.5670e4 * Q * pow(w, 4) * Rs[5] + 0.5670e4 * pow(w, 4) *
                    Qs[5] * R - 0.3240e4 * Rs[8] * Q * w + 0.3360e4 * w *
                    Qs[6] * Rs[3] + 0.2100e4 * pow(w, 4) * Qs[3] * Rs[3] -
                    0.1260e4 * pow(w, 6) * Q * Rs[3] - 0.10080e5 *
                    pow(w, 3) * Rs[6] * Q + 0.525e3 * Qs[9] * R - 0.480e3 *
                    pow(w, 7) * Rs[3] + 0.120e3 * pow(w, 9) * R - 0.120e3 *
                    pow(w, 9) * Q - 0.1400e4 * w * Qs[9] + 0.4320e4 *
                    pow(w, 2) * Qs[8] - 0.7200e4 * pow(w, 3) * Qs[7] +
                    0.6720e4 * pow(w, 4) * Qs[6] + 0.525e3 * Q * Rs[9] +
                    0.1134e4 * Qs[5] * Rs[5] - 0.900e3 * Qs[3] * Rs[7] +
                    0.7200e4 * pow(w, 3) * Rs[7] + 0.1400e4 * w * Rs[9] -
                    0.6720e4 * pow(w, 4) * Rs[6] + 0.3024e4 * pow(w, 5) *
                    Rs[5] - 0.4320e4 * pow(w, 2) * Rs[8] - 0.3024e4 *
                    pow(w, 5) * Qs[5]) / Qs[3] / w / Rs[3] / 0.25200e5;
        }

    }

    double covIntTP3NwGTR(double w, double Q, double R,
                          const double* Qs, const double* Rs) {
        if (w <= R + Q) {
            return 0.4e1 / 0.525e3 * acos(-1.0e0) * pow((R + Q - w), 6) *
                   (pow(w, 4) / 0.6e1 + (0.3e1 / 0.8e1 * Q + 0.3e1 / 0.8e1 *
                   R) * pow(w, 3) + (0.103e3 / 0.64e2 * Q * R - Qs[2] /
                   0.4e1 - Rs[2] / 0.4e1) * pow(w, 2) - 0.31e2 / 0.24e2 *
                   (Rs[2] - 0.247e3 / 0.124e3 * Q * R + Qs[2]) * (R + Q) *
                   w + Rs[4] - 0.209e3 / 0.64e2 * Qs[3] * R + 0.147e3 /
                   0.32e2 * Qs[2] * Rs[2] + Qs[4] - 0.209e3 / 0.64e2 * Q *
                   Rs[3]) / Qs[3] / w / Rs[3];
        }
        else {
            return 0;
        }

    }
    """

    # Effectively a type declaration for some parameters in the
    # following functions.
    c = numpy.eye(5, dtype=float)
    R = 1.
    Q = 1.
    cmin = 1
    cmax = 1
    symm = 1

    # Define functions args
    thinplate_arguments = ["c", "R", "cmin", "cmax", "symm"]
    innerproduct_arguments = ["c", "R", "Q", "symm"]

    # We need an extension function that calls this C function to do this.
    # This is possible by including the above code snippet as 'support code'
    # and then calling it from the extension function. 
    func = ext_tools.ext_function("thinplate3d_mat", thinplate3d_code,
                                  thinplate_arguments)
    mod.add_function(func)
    func = ext_tools.ext_function("innerproduct_thinplate3d",
        innerproduct_thinplate3d_code, innerproduct_arguments)
    mod.add_function(func) 
    mod.customize.add_support_code(thinplate_support_code)
    mod.customize.add_support_code(
        innerproduct_thinplate3d_normalized_support_code)

    # Suppress bogus #warning "Using deprecated NumPy API" warnings from weave
    # (unfortunately also suppresses any other explicit #warning)
    mod.customize.add_extra_compile_arg("-Wno-cpp")

    # Called to build the extension module. By default, the module is created
    # in the current working directory.
    mod.compile()
