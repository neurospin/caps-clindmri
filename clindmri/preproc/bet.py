#! /usr/bin/env python
##########################################################################
# CAPS - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
from __future__ import print_function
import os
import numpy

# Trait import
try:
    from traits.trait_base import _Undefined
    from traits.api import (List, Bool, File, Int, String, Float, Directory)
except ImportError:
    from enthought.traits.trait_base import _Undefined
    from enthought.traits.api import (List, Bool, File, Int, String,
                                      Float, Directory)

# Capsul import
from capsul.process import Process
from capsul.pipeline import Pipeline


##############################################################
#         Get Reference Image Process Definition
##############################################################

class ExtractReferenceB0Image(Process):
    """ Extract the first reference (b=0) image from a diffusion timeserie
    """

    def __init__(self):
        """ Initialize ExtractReferenceB0Image class
        """
        # Inheritance
        super(ExtractReferenceB0Image, self).__init__()

        # Inputs
        self.add_trait("dw_image", File(
            _Undefined(),
            optional=False,
            output=False,
            exists=True,
            desc="an existing diffusion weighted image"))
        self.add_trait("bvals", File(
            _Undefined(),
            optional=False,
            output=False,
            exists=True,
            desc="the the diffusion b-values"))
        self.add_trait("specified_index_of_ref_image", Int(
            _Undefined(),
            optional=True,
            output=False,
            desc="index of the reference b=0 volume"))

        # Outputs
        self.add_trait("index_of_ref_image", Int(
            _Undefined(),
            output=True,
            desc="index of the reference b=0 volume"))

    def _run_process(self):
        """ ExtractReferenceB0Image execution code
        """
        # Find the reference volume
        if not self.specified_index_of_ref_image:
            self.index_of_ref_image = self.find_dwi_reference()
        else:
            self.index_of_ref_image = self.specified_index_of_ref_image

    def find_dwi_reference(self):
        """ Method to get the first reference b=0 volume in a diffusion
        timeserie
        """
        # load diffusion metadata
        bvals = numpy.loadtxt(self.bvals).tolist()

        # check if the reference exists
        if 0 in bvals:
            arg_bvals = numpy.argsort(bvals)
            reference_volume = arg_bvals[0]
        else:
            raise Exception("No reference image found in diffusion"
                            "sequence {0}.".format(self.dw_image))

        return reference_volume

    run = property(_run_process)


##############################################################
#         Brain Extraction Pipeline Definition
##############################################################

class dBET(Pipeline):
    """ Diffusion FSL Brain Extraction Tool (BET).
    Deletes non-brain tissue from a b0 image.

    For complete details, see the `BET Documentation.
    <http://fsl.fmrib.ox.ac.uk/fsl/fsl-4.1.9/bet2/index.html>`_

    This pipeline include four processing steps:

    * Split time serie
    * :ref:`Extract the reference b=0 image (if necessary) <caps.diffusion_preproc.bet.ExtractReferenceB0Image>`
    * :ref:`Pick the reference volume <caps.utils.misc.Select>`
    * :ref:`Run BET on the picked image <caps.utils.bet.BET>`
    """

    def pipeline_definition(self):
        """ Diffusion BET pipeline definition
        """

        # Create processes
        self.add_process("extractor",
                         "caps.diffusion_preproc.bet.ExtractReferenceB0Image")
        self.add_process("spliter",
                         "nipype.interfaces.fsl.Split",
                         make_optional=["terminal_output", "dimension"])
        self.add_process("picker",
                         "caps.utils.misc.Select")
        self.add_process("bet",
                         "caps.utils.bet.BET")

        # Export Inputs (automatic export of BET parameters)
        self.export_parameter("extractor", "dw_image")
        self.export_parameter("extractor", "bvals")
        self.export_parameter("extractor", "specified_index_of_ref_image")

        # Link inputs
        self.add_link("dw_image->spliter.in_file")

        # Link b0_extractor
        self.add_link("extractor.index_of_ref_image->picker.index")

        # Link spliter
        self.add_link("spliter._out_files->picker.inlist")

        # Link picker
        self.add_link("picker.outelement->bet.input_file")

        # Export outputs
        self.export_parameter("extractor", "index_of_ref_image")
        self.export_parameter("spliter", "_out_files",
                              pipeline_parameter="splited_images")
        self.export_parameter("picker", "outelement",
                              pipeline_parameter="ref_image")

        # dBET algorithm parameters
        self.dw_image = _Undefined()
        self.bvals = _Undefined()
        self.nodes["spliter"].process.dimension = "t"


##############################################################
#                            Pilot
##############################################################

def pilot(working_dir="/volatile/nsap/caps", **kwargs):
    """
    ===============================
    Diffusion Brain Extraction Tool
    ===============================
    .. topic:: Objective

        We propose to extract the brain mask from a diffusion sequence.

    Import
    ------

    First we load the function that enables us to access the toy datasets
    """
    from caps.toy_datasets import get_sample_data

    """
    From capsul we then load the class to configure the study we want to
    perform
    """
    from capsul.study_config import StudyConfig

    """
    Here two utility tools are loaded. The first one enables the creation
    of ordered dictionary and the second ensure that a directory exist.
    Note that the directory will be created if necessary.
    """
    from capsul.utils.sorted_dictionary import SortedDictionary
    from nsap.lib.base import ensure_is_dir

    """
    Load the toy dataset
    --------------------

    We want to perform BET on a diffusion sequence.
    To do so, we use the *get_sample_data* function to load this
    dataset.

    .. seealso::

        For a complete description of the *get_sample_data* function, see the
        :ref:`Toy Datasets documentation <toy_datasets_guide>`
    """
    toy_dataset = get_sample_data("dwi")

    """
    The *toy_dataset* is an Enum structure with some specific
    elements of interest *dwi*, *bvals* that contain the nifti diffusion
    image and the b-values respectively.
    """
    print(toy_dataset.dwi, toy_dataset.bvals)

    """
    Will return:

    .. code-block:: python

        /home/ag239446/git/nsap-src/nsap/data/DTI30s010.nii
        /home/ag239446/git/nsap-src/nsap/data/DTI30s010.bval

    We can see that the image has been found in a local directory

    Processing definition
    ---------------------

    Now we need to define the processing step that will perform BET on
    diffusion sequence.
    """
    bet_pipeline = dBET()

    """
    It is possible to access the ipeline input specification.
    """
    print(bet_pipeline.get_input_spec())

    """
    Will return the input parameters the user can set:

    .. code-block:: python

        INPUT SPECIFICATIONS

        dw_image: ['File']
        bvals: ['File']
        specified_index_of_ref_image: ['Int']
        terminal_output: ['Enum']
        generate_binary_mask: ['Bool']
        use_4d_input: ['Bool']
        generate_mesh: ['Bool']
        generate_skull: ['Bool']
        bet_threshold: ['Float']

    We can now tune the pipeline parameters.
    We first set the input dwi file:
    """
    bet_pipeline.dw_image = toy_dataset.dwi

    """
    And set the b-values file
    """
    bet_pipeline.bvals = toy_dataset.bvals

    """
    Study Configuration
    -------------------

    The pipeline is now set up and ready to be executed.
    For a complete description of a study execution, see the
    :ref:`Study Configuration description <study_configuration_guide>`
    """
    bet_working_dir = os.path.join(working_dir, "diffusion_bet")
    ensure_is_dir(bet_working_dir)
    default_config = SortedDictionary(
        ("output_directory", bet_working_dir),
        ("fsl_config", "/etc/fsl/4.1/fsl.sh"),
        ("use_fsl", True),
        ("use_smart_caching", True),
        ("generate_logging", True)
    )
    study = StudyConfig(default_config)
    study.run(bet_pipeline)

    """
    Results
    -------

    Finally, we print the pipeline outputs
    """
    print("\nOUTPUTS\n")
    for trait_name, trait_value in bet_pipeline.get_outputs().iteritems():
        print("{0}: {1}".format(trait_name, trait_value))

    """
    .. note::
        Since only the brain mask has been requested, all the other outputs
        are set to None.
        Only the *bet_out_file*, *splited_images*, *bet_mask_file*,
        *ref_image*, *index_of_ref_image* outputs are of interest for
        this study.
    """


if __name__ == '__main__':
    pilot()