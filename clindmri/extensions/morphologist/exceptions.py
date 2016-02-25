##########################################################################
# NSAp - Copyright (C) CEA, 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################


class MorphologistError(Exception):
    """ Base exception type for the package.
    """
    def __init__(self, message):
        super(MorphologistError, self).__init__(message)


class MorphologistRuntimeError(MorphologistError):
    """ Error thrown when call to the Morphologist software failed.
    """
    def __init__(self, algorithm_name, parameters, error=None):
        message = (
            "Morphologist call for '{0}' failed.Error:: "
            "{2}.".format(algorithm_name, parameters, error))
        super(MorphologistRuntimeError, self).__init__(message)

