##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os

# Clindmri import
from clindmri.extensions.fsl import FSLWrapper
from clindmri.extensions.fsl.exceptions import FSLRuntimeError


def dtifit(k, r, b, m, o, w=None):
    """ Wraps command dtifit.

    Usage:
    dtifit -k <filename>
    dtifit --verbose

    Compulsory arguments (You MUST set one or more of):
        -k,--data    dti data file
        -o,--out    Output basename
        -m,--mask    Bet binary mask file
        -r,--bvecs    b vectors file
        -b,--bvals    b values file

    Optional arguments (You may optionally specify one or more of):
        -V,--verbose    switch on diagnostic messages
        -h,--help    display this message
        --cni    Input confound regressors
        --sse    Output sum of squared errors
        -w,--wls    Fit the tensor with weighted least squares
        --littlebit    Only process small area of brain
        --save_tensor    Save the elements of the tensor
        -z,--zmin    min z
        -Z,--zmax    max z
        -y,--ymin    min y
        -Y,--ymax    max y
        -x,--xmin    min x
        -X,--xmax    max x

    Returns
    -------
    v1_file: str
        1st eigenvector
    v2_file: str
        2nd eigenvector
    v3_file: str
        3rd eigenvector
    l1_file: str
        1st eigenvalue
    l2_file: str
        2nd eigenvalue
    l3_file: str
        3rd eigenvalue
    md_file: str
        mean diffusivity
    fa_file: str
        fractional anisotropy
    s0_file: str
        raw T2 signal with no diffusion weighting
    tensor_file: str
        the full second order coefficients
    m0_file: str
        the mode of the anisotropy
    """
    # Check that the output directory exists
    if not os.path.isdir(o):
        os.makedirs(o)
    o = os.path.join(o, "dtifit")

    # Call bedpostx
    fslprocess = FSLWrapper("dtifit --save_tensor")
    fslprocess()
    if fslprocess.exitcode != 0:
        raise FSLRuntimeError(fslprocess.cmd[0], " ".join(fslprocess.cmd[1:]))

    # Build outputs
    image_ext = fslprocess.output_ext[fslprocess.environment["FSLOUTPUTTYPE"]]
    v1_file = o + "_V1" + image_ext
    v2_file = o + "_V2" + image_ext
    v3_file = o + "_V3" + image_ext
    l1_file = o + "_L1" + image_ext
    l2_file = o + "_L2" + image_ext
    l3_file = o + "_L3" + image_ext
    md_file = o + "_MD" + image_ext
    fa_file = o + "_FA" + image_ext
    s0_file = o + "_S0" + image_ext
    tensor_file = o + "_tensor" + image_ext
    m0_file = o + "_M0" + image_ext

    return (v1_file, v2_file, v3_file, l1_file, l2_file, l3_file, md_file,
            fa_file, s0_file, tensor_file, m0_file)
