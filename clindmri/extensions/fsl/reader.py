##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import numpy
import nibabel


def flirt2aff(mat_file, in_file, ref_file):
    """ Map from 'in_file' image voxels to 'ref_file' voxels given `omat` FSL
    affine transformation.

    Parameters
    ------------
    mat_file: str (mandatory)
        filename of output '-omat' transformation file from FSL flirt.
    in_file: str (mandatory)
        filename of the image passed to flirt as the '-in' image.
    ref_file: str (mandatory)
        filename of the image passed to flirt as the '-ref' image.

    Returns
    -------
    omat: array (4, 4)
        array containing the transform from voxel coordinates in image
        for 'in_file' to voxel coordinates in image for 'ref_file'.
    """
    # Load dataset
    flirt_affine = numpy.loadtxt(mat_file)
    in_img = nibabel.load(in_file)
    ref_img = nibabel.load(ref_file)
    in_hdr = in_img.get_header()
    ref_hdr = ref_img.get_header()

    # Define a function to flip x
    def _x_flipper(n):
        flipr = numpy.diag([-1, 1, 1, 1])
        flipr[0, 3] = n - 1
        return flipr

    # Check image orientation
    inspace = numpy.diag(in_hdr.get_zooms() + (1, ))
    refspace = numpy.diag(ref_hdr.get_zooms() + (1, ))
    if numpy.linalg.det(in_img.get_affine()) >= 0:
        inspace = numpy.dot(inspace, _x_flipper(in_hdr.get_data_shape()[0]))
    if numpy.linalg.det(ref_img.get_affine()) >= 0:
        refspace = numpy.dot(refspace, _x_flipper(ref_hdr.get_data_shape()[0]))

    # Get the voxel to voxel mapping
    omat = numpy.dot(numpy.linalg.inv(refspace),
                     numpy.dot(flirt_affine, inspace))

    return omat
