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
import inspect

# Clindmri import
from clindmri.extensions.configuration import environment


class FSLWrapper(object):
    """ Parent class for the wrapping of FSL functions. 
    """   
    def __init__(self, name, shfile="/etc/fsl/4.1/fsl.sh"):
        """ Initialize the FSLWrapper class by setting properly the
        environment.

        Parameters
        ----------
        name: str (mandatory)
            the name of the FSL binary to be called.
        shfile: str (optional, default NeuroSpin path)
            the path to the FSL 'fsl.sh' configuration file.
        """
        self.name = name
        self.cmd = [name]
        self.shfile = shfile
        self.environment = self._environment()

    def __call__(self):
        """ Run the FSL command.
        """
        # Update the command to execute
        self._update_command()

        # Execute the command
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
        env = os.environ.get("FSL_CONFIGURED", None)

        # Configure FreeSurfer
        if env is None:

            # Parse configuration file
            env = environment(self.shfile)

            # Save the result
            os.environ["FSL_CONFIGURED"] = json.dumps(env)

        # Load configuration
        else:
            env = json.loads(env)  

        return env

    def _update_command(self):
        """ Update the command that will be executed.
        """
        # Get the caller frame parameters
        caller_frame = inspect.stack()[2][0]
        args, _, _, values = inspect.getargvalues(caller_frame)

        # Update the command
        for parameter_name in args:
            if values[parameter_name] is not None:

                # File parameter name extension
                if parameter_name.endswith("_file"):
                    self.cmd.append("-{0}".format(
                        parameter_name.replace("_file", "")))
                else:
                    self.cmd.append("-{0}".format(parameter_name))

                # Boolean value
                if isinstance(values[parameter_name], bool):
                    if not values[parameter_name]:
                        self.cmd.pop()
                else:
                    self.cmd.append("{0}".format(values[parameter_name]))
