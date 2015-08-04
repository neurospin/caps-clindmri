#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" Package's exceptions. """

from .manufacturers import MANUFACTURERS


class ConnectomistError(Exception):
    """ Base exception type for the package.
    """
    def __init__(self, message):
        super(ConnectomistError, self).__init__(message)


class ConnectomistRuntimeError(ConnectomistError):
    """ Error thrown when call to the Connectomist software failed.
    """
    def __init__(self, algorithm_name, parameter_file):
        message = "Connectomist call for %s failed, with parameters: %s." \
                  % (algorithm_name, parameter_file)
        super(ConnectomistRuntimeError, self).__init__(message)


class BadManufacturerNameError(ConnectomistError):
    """ Error thrown when an incorrect manufacturer name is detected.
    """
    def __init__(self, manufacturer):
        message = "Incorrect manufacturer name: %s, should be in {}." \
                  .format(set(MANUFACTURERS))
        super(BadManufacturerNameError, self).__init__(message)


class MissingParametersError(ConnectomistError):
    """ Error thrown when needed parameters were not all given.
    """
    def __init__(self, algorithm_name, missing_parameters):
        message = "Missing parameters for {}: {}.".format(algorithm_name,
                                                          missing_parameters)
        super(MissingParametersError, self).__init__(message)


class BadFileError(ConnectomistError):
    """ Error thrown when a file is missing or corrupted.
    """
    def __init__(self, file_path):
        message = "Missing or corrupted file: %s" % file_path
        super(BadFileError, self).__init__(message)
