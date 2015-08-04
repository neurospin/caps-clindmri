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
import subprocess
import json

# Clindmri import
from clindmri.extensions.configuration import environment


class FSWrapper(object):
    """ Parent class for the wrapping of FreeSurfer functions. 
    """   
    def __init__(self, cmd, shfile="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
        """ Initialize the FSWrapper class by setting properly the
        environment.

        Parameters
        ----------
        cmd: list of str (mandatory)
            the FreeSurfer command to execute.
        shfile: str (optional, default NeuroSpin path)
            the path to the FreeSurfer 'SetUpFreeSurfer.sh' configuration file.
        """
        self.cmd = cmd
        self.shfile = shfile
        self.environment = self._environment()

    def __call__(self):
        """ Run the FreeSurfer command.
        """
        process = subprocess.Popen(
            self.cmd,
            env=self.environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        self.stdout, self.stderr = process.communicate()
        self.exitcode = process.returncode

    def _environment(self) :
        """ Return a dictionary of the environment needed by FreeSurfer
        binaries.
        """
        # Check if FreeSurfer has already been configures
        env = os.environ.get("FREESURFER_CONFIGURED", None)

        # Configure FreeSurfer
        if env is None:

            # FreeSurfer home directory    
            fs_home = os.environ.get("FREESURFER_HOME", None)
            env = {}
            if fs_home is not None:
                env["FREESURFER_HOME"] = fs_home

            # Parse configuration file
            env = environment(self.shfile, env)

            # Save the result
            os.environ["FREESURFER_CONFIGURED"] = json.dumps(env)

        # Load configuration
        else:
            env = json.loads(env)  

        return env
