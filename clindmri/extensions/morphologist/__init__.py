#! /usr/bin/env python
##########################################################################
# NSAP - Copyright (C) CEA, 2013-2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import sys
import os
import re
import subprocess


def set_bvenv(bvenv_config):
    """ Execute a new program, replacing the current process in the brainvisa
    environment. On Unix, the new executable is loaded into the current
    process, and will have the same process id as the caller.

    Parameters
    ----------
    bvenv_config: str (mandatory)
        the path to the brainvisa 'bv_env' binary.
    """
    # Execute bv_env
    command = [os.path.abspath(bvenv_config)]
    process = subprocess.Popen(command, env={},
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(
            "Could not exucute 'bv_env' {0}. Maybe you should check "
            "the path".format(stderr))

    # Parse the output : each line should be of the form
    # 'VARIABLE_NAME=value'
    bvenv = {}
    for line in stdout.split(os.linesep):
        if line.startswith("export"):
            line = line.replace("export ", "")
            line = line.replace("'", "")
        match = re.match(r"^(\w+)=(\S*)$", line)
        if match:
            name, value = match.groups()
            if name != "PWD":
                bvenv[name] = value

    # Update current environment
    if not bvenv["LD_LIBRARY_PATH"] in (os.environ.get("LD_LIBRARY_PATH", [])):
        for envname, envval in bvenv.items():
            if envname in os.environ:
                os.environ[envname] = envval + ":" + os.environ[envname] 
            else:
                os.environ[envname] = envval

        try:
            os.execve(sys.executable, [sys.executable, ] + sys.argv,
                      os.environ)
        except:
            raise Exception("The automatic configuration of brainvisa"
                            "fails. Please report the bug."
                            "You can bypass this error using a manual"
                            "setting: bv_env python ... or"
                            "bv_env ipython.")

