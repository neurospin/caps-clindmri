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
import inspect
import types

# Clindmri import
import modulehacker


_modules = {}


def register(decorator, modules):
    """ Function to register a decorator for a list of module names.

    Parameters
    ----------
    decorator: callable (mandatory)
        a decorator function.
    modules: list of str (mandatory)
        a list of module names whose functions will be decorated.
    """
    for module in modules:
        if module not in _modules:
            _modules[module] = []
        _modules[module].append(decorator)


class Decorations(object):
    """ A class that docorate a module functions based on the factory, ie.
    the '_modules' mapping.
    """
    def hack(self, module):
        module_name = module.__name__.split(".")[0]
        for decorator in _modules.get(module_name, ()):
            self.decorate(module, decorator)
        return module

    def decorate(self, module, decorator):
        module_name = module.__name__.split(".")[0]
        for module_attr, module_object in module.__dict__.items():
            if isinstance(module_object, types.FunctionType):
                setattr(module, module_attr, decorator(module_object,
                                                       module_name))


modulehacker.register(Decorations())




