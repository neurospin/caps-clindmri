#! /usr/bin/env python
##########################################################################
# NSAP - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import inspect
import time
import json

# Clindmri import
import decorations


def function_decorator(func, module_name):
    """ Create a decorator that display the function signature and the
    function execution time.
    """
    def wrapper(*args, **kwargs):
        # Get the function parameters
        arg_spec = inspect.getargspec(func)
        defaults = [repr(item) for item in arg_spec.defaults or []]
        optional = dict(zip(reversed(arg_spec.args or []), reversed(defaults)))
        for name, value in kwargs.items():
            if name in optional:
                optional[name] = repr(value)
        mandatory = []
        for index in range(len(arg_spec.args) - len(optional)):
            print "--", index
            print "--", arg_spec.args[index]
            try:
                if isinstance(value, list):
                    value = repr(numpy.asarray(value))[6:-1]
                else:
                    value = repr(args[index])
                #representation = repr(args[index])
                #if len(representation) > 200:
                #    representation = (
                #        representation[:100] + "..." + representation[-100:])
                #print representation
                mandatory.append((arg_spec.args[index], value))
            except:
                mandatory.append((arg_spec.args[index], None))

        # Create the function signature
        params = ["{0}={1}".format(name, value) for name, value in mandatory]
        params.extend([
            "{0}={1}".format(name, value) for name, value in optional.items()])
        signature = "{0}({1})".format(func.__name__, ", ".join(params))

        # Display a start call message
        print("{0}\n[{1}] Calling {2}...\n{3}".format(
            80 * "_", module_name, func.__module__ + "." + func.__name__,
            signature))

        # Call
        start_time = time.time()
        returncode = func(*args, **kwargs)
        duration = time.time() - start_time

        # Display an end message
        msg = "{0:.1f}s, {1:.1f}min".format(duration, duration / 60.)
        print(max(0, (80 - len(msg))) * '_' + msg)

        return returncode

    return wrapper


decorations.register(function_decorator, ["clindmri"])

