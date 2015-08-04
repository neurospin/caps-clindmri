#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
Package to wrap the Connectomist software and simplify scripting calls.

Each preprocessing step (i.e. one Connectomist tab) can be run through the use
of a dedicated function of the package.

All the preprocessing steps can be done at once using the 
complete_preprocessing() function.

The package also allows parallelization of multiple complete preprocessing 
using the parallel_preprocessing() function.
"""

from .manufacturers import MANUFACTURERS
from .exceptions    import (ConnectomistError, ConnectomistRuntimeError,
                            BadManufacturerNameError, MissingParametersError,
                            BadFileError)
from .complete_preprocessing import complete_preprocessing
from .parallel_preprocessing import parallel_preprocessing
