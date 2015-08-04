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

# Clindmri import
from clindmri.extensions.fsl import FSLWrapper
from clindmri.extensions.fsl.exceptions import FSLRuntimeError


def flirt(in_file, ref_file, omat=None, out=None, init=None, cost="corratio",
          usesqform=None, displayinit=None, anglerep="euler", bins=256,
          interp="trilinear", dof=12, applyxfm=None, verbose=0):
    """ Wraps command flirt.
    """
    # Set default parameters
    dirname = os.path.dirname(in_file)
    basename = os.path.basename(in_file).split(".")[0]
    if out is None:
        out = os.path.join(dirname, "flirt_out_{0}.nii.gz".format(basename))
    if omat is None:
        omat = os.path.join(dirname, "flirt_omat_{0}".format(basename))

    # Call flirt
    fslprocess = FSLWrapper("flirt")
    fslprocess()
    if fslprocess.exitcode != 0:
        raise FSLRuntimeError(fslprocess.cmd[0], " ".join(fslprocess.cmd[1:]))

    return out, omat
