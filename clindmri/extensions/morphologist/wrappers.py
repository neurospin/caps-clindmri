##########################################################################
# NSAp - Copyright (C) CEA, 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import subprocess


class MorphologistWrapper(object):
    """ Parent class for the wrapping of Morphologist functions.
    """
    def __init__(self, cmd):
        """ Initialize the MorphologistWrapper class

        Parameters
        ----------
        cmd: list of str (mandatory)
            the Morphologist command to execute.
        """
        self.cmd = cmd
        self.environment = os.environ

    def __call__(self):
        """ Run the Morphologist command.
        """
        # Execute the command
        process = subprocess.Popen(
            self.cmd,
            env=self.environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        self.stdout, self.stderr = process.communicate()
        self.exitcode = process.returncode
