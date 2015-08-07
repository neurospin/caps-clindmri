#! /usr/bin/env python
##########################################################################
# NSAP - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
#
# From: http://code.activestate.com/recipes/577742
##########################################################################

# System import
import sys
import importlib


_hackers = []


def register(obj):
    """ Function to register a new hacker.
    """
    _hackers.append(obj)


class Loader(object):
    """ A class that import a module like normal and then passed to a hacker
    object that gets to do whatever it wants to the module. Then the return
    value from the hack call is put into sys.modules.
    """
    def __init__(self):
        self.module = None
    
    def find_module(self, name, path):
        sys.meta_path.remove(self)
        try:
            self.module = importlib.import_module(name)
        finally:
            sys.meta_path.insert(0, self)
        return self
    
    def load_module(self, name):
        if not self.module:
            raise ImportError("Unable to load module.")
        module = self.module
        for hacker in _hackers:
            module = hacker.hack(module)
        sys.modules[name] = module
        return module


sys.meta_path.insert(0, Loader())
