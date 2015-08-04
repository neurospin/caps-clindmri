#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################


class FSLError(Exception):
    """ Base exception type for the package.
    """
    def __init__(self, message):
        super(FSLError, self).__init__(message)


class FSLRuntimeError(FSLError):
    """ Error thrown when call to the FSL software failed.
    """
    def __init__(self, algorithm_name, parameter_file):
        message = (
            "FSL call for '{0}' failed, with parameters: '{1}'.".format(
                algorithm_name, parameter_file))
        super(FSLRuntimeError, self).__init__(message)

