##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import numpy
import math
import os


def select_first_b0(bvalfile):
    """ Select the fist b0 image in a diffusion sequence.

    <process>
        <return name="b0_index" type="Int" desc="First index of non weighted
            image."/>
        <input name="bvalfile" type="File" desc="The diffusion b-values
            file."/>
    </process>
    """
    bvals = numpy.loadtxt(bvalfile).tolist()
    try:
        b0_index = bvals.index(0)
    except:
        raise Exception("No non-weighted image can be detected in file "
                        "'{0}'.".format(bvalfile))
    return b0_index


def rotate_bvecs(bvecfile, trffile, output_directory):
    """ Reorient the diffusion b-vectors from FSL .par file.

    In FSL order of fields is three Euler angles (x,y,z in radians) then
    three translation parameters (x,y,z in mm).

    <process>
        <return name="reobvecfile" type="File" desc="The reoriented
            b-vectors."/>
        <input name="bvecfile" type="File" desc="The diffusion b-vectors
            file."/>
        <input name="trffile" type="File" desc="The transformation parameters
            used to align each volume of the diffusion sequence."/>
        <input name="output_directory" type="Directory" desc="The destination
            folder."/>
    </process>
    """
    # Generate output file name
    fname = os.path.basename(bvecfile).split(".")[0]
    reobvecfile = os.path.join(output_directory,
                               "{0}_rotated.bvec".format(fname))

    # Load b-vectors and transformation parameters
    bvecs = numpy.loadtxt(bvecfile).T
    trf = numpy.loadtxt(trffile)
    if (bvecs.shape[0] != trf.shape[0] or bvecs.shape[1] != 3 or
            trf.shape[1] != 6):
        raise Exception(
            "The b-vector file '{0}' and transformation file '{1}' "
            "are miss-formated.".format(bvecfile, trffile))
    nbdirs = bvecs.shape[0]

    # Reorient the directions
    rbvecs = numpy.zeros(bvecs.shape, dtype=bvecs.dtype)
    for index in range(nbdirs):
        rotation = euler_matrix(*trf[index][:3])[:3, :3]
        rbvecs[index] = numpy.dot(rotation, bvecs[index].T)

    # Save the reoriented b-vectors
    numpy.savetxt(reobvecfile, rbvecs.T, fmt="%0.15f")

    return reobvecfile


def euler_matrix(ax, ay, az):
    """ Return homogeneous rotation matrix from Euler angles.

    Parameters
    ----------
    ax, ay, az: float (mandatory)
        Euler's roll, pitch and yaw angles.

    Returns
    -------
    M: array
        a 4 by 4 array containing the rotation matrix in homogeneous
        coordinates.
    """
    sx, sy, sz = math.sin(ax), math.sin(ay), math.sin(az)
    cx, cy, cz = math.cos(ax), math.cos(ay), math.cos(az)
    ccxz, csxz = cx * cz, cx * sz
    scxz, ssxz = sx * cz, sx * sz

    M = numpy.identity(4)
    M[0, 0] = cy * cz
    M[0, 1] = sy * scxz - csxz
    M[0, 2] = sy * ccxz + ssxz
    M[1, 0] = cy * sz
    M[1, 1] = sy * ssxz + ccxz
    M[1, 2] = sy * csxz - scxz
    M[2, 0] = - sy
    M[2, 1] = cy * sx
    M[2, 2] = cy * cx

    return M
