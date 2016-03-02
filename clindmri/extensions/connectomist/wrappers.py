###############################################################
# NSAp - Copyright (C) CEA, 2015 - 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

# System import
import os
import time
import pprint
import subprocess

# Clindmri import
from .exceptions import ConnectomistError
from .exceptions import ConnectomistConfigurationError
from .exceptions import ConnectomistRuntimeError


class ConnectomistWrapper(object):
    """ Parent class for the wrapping of Connectomist functions.
    """
    # Map algorithm name to a list of files that should be created.
    # It is meant to check that the call to Connectomist worked, because the
    # exit code is 0 even when it fails.
    files_to_check = {
        "DWI-Data-Import-And-QSpace-Sampling": [
            "t2.ima", "dw.ima", "acquisition_parameters.py"],
        "DWI-Rough-Mask-Extraction": ["mask.ima"],
        "DWI-Outlier-Detection": ["t2_wo_outlier.ima", "dw_wo_outlier.ima"],
        "DWI-Susceptibility-Artifact-Correction": [
            "t2_wo_susceptibility.ima", "dw_wo_susceptibility.ima"],
        "DWI-Eddy-Current-And-Motion-Correction": [
            "t2_wo_eddy_current_and_motion.ima",
            "dw_wo_eddy_current_and_motion.ima"],
        "DWI-To-Anatomy-Matching": [
            "talairach_to_t1.trm", "dw_to_t1.trm", "t1_to_dw.trm",
            "talairach_to_t1.trm"],
        "DWI-Local-Modeling": [],
        "DWI-Tractography-Mask": ["tractography_mask.ima"],
        "DWI-Tractography": []
    }

    def __init__(self, path_connectomist=(
            "/i2bm/local/Ubuntu-14.04-x86_64/ptk/bin/connectomist")):
        """ Initialize the ConnectomistWrapper class by setting properly the
        environment.

        Parameters
        ----------
        path_connectomist: str (optional)
            path to the Connectomist executable.

        Raises
        ------
        ConnectomistConfigurationError: If Connectomist is not configured.
        """
        # Class parameters
        self.path_connectomist = path_connectomist
        self.environment = os.environ

        # Check Connectomist has been configured so the command can be found
        cmd = "%s --help" % (self.path_connectomist)
        process = subprocess.Popen(
            cmd, shell=True,
            env=self.environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        self.stdout, self.stderr = process.communicate()
        self.exitcode = process.returncode
        if self.exitcode != 0:
            raise ConnectomistConfigurationError(self.path_connectomist)

    def __call__(self, algorithm, parameter_file, outdir, nb_tries=10):
        """ Run the Connectomist 'algorithm' (tab in ui).

        Parameters
        ----------
        algorithm: str
            name of Connectomist's tab in ui.
        paramter_file: str
            path to the parameter file for the tab in ui: executable python
            file to set the connectomist tab input parameters.
        outdir: str
            path to directory where the algorithm outputs.
        nb_tries: int (optional, default 10)
            nb of times to try an algorithm if it fails.
            It often crashes when running in parallel. The reason
            why it crashes is unknown.

        Raises
        ------
        ConnectomistError: If Connectomist call failed.
        """
        # List filenames to be checked for this call
        to_check = [os.path.join(outdir, f)
                    for f in self.files_to_check.get(algorithm, [])]

        # Command to be run
        cmd = "%s -p %s -f %s" % (self.path_connectomist, algorithm,
                                  parameter_file)

        # Run the command, multiple times if it fails.
        nb_tried = 0
        while nb_tried < nb_tries:
            nb_tried += 1

            try:
                subprocess.check_call(cmd, shell=True)
            except:
                pass

            success = all(map(os.path.isfile, to_check))
            if success:
                self.exitcode = 0
                return
            else:
                time.sleep(3)  # wait 3 sec before retrying

        # If the command failed nb_tries times, raise an exception
        self.exitcode = 1
        raise ConnectomistError("Connectomist call failed with cmd:\n%s" % cmd)

    @classmethod
    def create_parameter_file(cls, algorithm, parameters_dict, outdir):
        """  Writes the '.py' file that Connectomist uses when working in
        command line.

        Parameters
        ----------
        algorithm: str
            name of Connectomist's tab.
        parameters_dict: dict
            parameter values for the tab.
        outdir: str
            path to directory where to write the parameter file.
            If not existing the directory is created.

        Returns
        -------
        parameter_file: str
            path to the created parameter file for the tab in ui: executable
            python file to set the connectomist tab input parameters.
        """
        # If not existing create outdir
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        # Write the parameter file
        parameter_file = os.path.join(outdir, "%s.py" % algorithm)
        with open(parameter_file, "w") as f:
            f.write("algorithmName = '%s'\n" % algorithm)
            # Pretty text to write, without the first "{"
            pretty_dict = pprint.pformat(parameters_dict)[1:]
            f.write("parameterValues = {\n " + pretty_dict)

        return parameter_file


class PtkWrapper(object):
    """ Parent class for the wrapping of Connectomist Ptk functions.
    """
    def __init__(self, cmd):
        """ Initialize the PtkWrapper class

        Parameters
        ----------
        cmd: list of str (mandatory)
            the Morphologist command to execute.
        """
        # Class parameter
        self.cmd = cmd
        self.environment = os.environ

        # Check Connectomist Ptk has been configured so the command can b
        # found
        process = subprocess.Popen(
            ["which", self.cmd[0]],
            env=self.environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        self.stdout, self.stderr = process.communicate()
        self.exitcode = process.returncode
        if self.exitcode != 0:
            raise ConnectomistConfigurationError(self.cmd[0])

    def __call__(self, files_to_check, nb_tries=10):
        """ Run the Connectomist Ptk command.

        Parameters
        ----------
        files_to_check: list of str
            list of expected output files. It is meant to check that the call
            to Connectomist worked, because the exit code is 0 even when it
            fails.
        nb_tries: int (optional, default 10)
            nb of times to try an algorithm if it fails.
            It often crashes when running in parallel. The reason
            why it crashes is unknown.
        """
        # Execute the command
        nb_tried = 0
        while nb_tried < nb_tries:

            nb_tried += 1
            process = subprocess.Popen(
                self.cmd,
                env=self.environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            self.stdout, self.stderr = process.communicate()
            self.exitcode = process.returncode

            success = all(map(os.path.isfile, files_to_check))
            if success:
                self.exitcode = 0
                return
            else:
                time.sleep(3)  # wait 3 sec before retrying

        # If the command failed nb_tries times, raise an exception
        self.exitcode = 1
        raise ConnectomistRuntimeError(self.cmd[0], self.cmd[1:], self.stderr)
