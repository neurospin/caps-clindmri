#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################


class FreeSurferError(Exception):
    """ Base exception type for the package.
    """
    def __init__(self, message):
        super(FreeSurferError, self).__init__(message)


class FreeSurferRuntimeError(FreeSurferError):
    """ Error thrown when call to the FreeSurfer software failed.
    """
    def __init__(self, algorithm_name, parameter_file, error=None):
        message = (
            "FreeSurfer call for '{0}' failed, with parameters: '{1}' "
            "failed: {2}.".format(algorithm_name, parameter_file, error))
        super(FreeSurferRuntimeError, self).__init__(message)


class FreeSurferConfigurationError(FreeSurferError):
    """ Error thrown when call to the FreeSurfer software failed.
    """
    def __init__(self, command_name):
        message = "FreeSurfer command '{0}' not found.".format(command_name)
        super(FreeSurferConfigurationError, self).__init__(message)

