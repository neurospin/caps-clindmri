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
from .exceptions import FSLConfigurationError
from .exceptions import FSLDependencyError


class FSLWrapper(object):
    """ Parent class for the wrapping of FSL functions. 
    """  
    output_ext = {
        "NIFTI_PAIR" : ".hdr",
        "NIFTI" : ".nii",
        "NIFTI_GZ" : ".nii.gz",
        "NIFTI_PAIR_GZ" : ".hdr.gz",
    }
 
    def __init__(self, name, shfile="/etc/fsl/5.0/fsl.sh", optional=None,
                 cpus=None):
        """ Initialize the FSLWrapper class by setting properly the
        environment.
        
        Parameters
        ----------
        name: str (mandatory)
            the name of the FSL binary to be called.
        shfile: str (optional, default NeuroSpin path)
            the path to the FSL 'fsl.sh' configuration file.
        optional: list (optional, default None)
            the name of the optional parameters. If 'ALL' consider that all
            the parameter are optional.
        cpus: str (optional, default None)
            if different that None use condor wiht the specified number of
            jobs.
        """
        self.name = name
        self.cmd = name.split()
        self.shfile = shfile
        self.optional = optional or []
        self.environment = self._environment()

        # Condor specific setting
        if cpus is not None:
            self.environment["FSLPARALLEL"] = cpus
            self.environment["USER"] = os.getlogin()
            process = subprocess.Popen(
                        ["which", "condor_qsub"],
                        env=self.environment,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
            self.stdout, self.stderr = process.communicate()
            self.exitcode = process.returncode
            if self.exitcode != 0:
                raise FSLDependencyError("condor_qsub", "Condor")        

    def __call__(self):
        """ Run the FSL command.

        Note that the command is built from the parent frame and that the
        'shfile' parameter is reserved.
        """
        # Update the command to execute
        self._update_command()

        # Check FSL has been configured so the command can be found
        process = subprocess.Popen(
                    ["which", self.cmd[0]],
                    env=self.environment,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
        self.stdout, self.stderr = process.communicate()
        self.exitcode = process.returncode
        if self.exitcode != 0:
            raise FSLConfigurationError(self.cmd[0])

        # Execute the command
        process = subprocess.Popen(
            self.cmd,
            env=self.environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        self.stdout, self.stderr = process.communicate()
        self.exitcode = process.returncode

    def _environment(self) :
        """ Return a dictionary of the environment needed by FSL binaries.
        """
        # Check if FSL has already been configures
        env = os.environ.get("FSL_CONFIGURED", None)
       
        # Configure FSL
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

        # 'ALL' optional case
        if self.optional == 'ALL':
            self.optional = args

        # Update the command
        for parameter_name in args:
        
            # Skip 'shfile' and 'cpus' parameters
            if parameter_name in ["shfile", "cpus"]:
                continue

            # Get parameter value
            parameter_value = values[parameter_name]

            # Clean parameter name
            cmd_parameter_name = parameter_name
            if parameter_name.endswith("_file"):
                cmd_parameter_name = parameter_name.replace("_file", "")

            if parameter_value is not None:

                # Mandatory parameter
                if parameter_name in ["input", "output"]:
                    self.cmd.append("{0}".format(parameter_value))

                # Boolean parameter
                elif isinstance(parameter_value, bool) and parameter_value:
                    if parameter_name in self.optional:
                        self.cmd.append("--{0}".format(cmd_parameter_name))
                    else:
                        self.cmd.append("-{0}".format(cmd_parameter_name))

                # Add command parameter
                elif not isinstance(parameter_value, bool):
                    if parameter_name in self.optional:
                        self.cmd.append("--{0}={1}".format(
                            cmd_parameter_name, parameter_value))
                    else:
                        self.cmd.append("-{0}".format(cmd_parameter_name))
                        self.cmd.append("{0}".format(parameter_value))
