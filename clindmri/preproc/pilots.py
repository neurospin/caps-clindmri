#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# CAPSUL import
try:
    from capsul.utils.pilot import pilotfunction
except:
    def pilotfunction(func):
        return func


@pilotfunction
def pilot_fsl_preproc():
    """
    FSL preprocessings
    ==================
    """
    # System import
    import os
    import sys
    import datetime
    import PySide.QtGui as QtGui

    # CAPSUL import
    from capsul.qt_gui.widgets import PipelineDevelopperView
    from capsul.study_config.study_config import StudyConfig
    from capsul.process.loader import get_process_instance

    """
    Study configuration
    -------------------

    We first define the working directory and guarantee this folder exists on
    the file system:
    """
    working_dir = "/volatile/nsap/clindmri/fslpreproc"
    if not os.path.isdir(working_dir):
        os.makedirs(working_dir)

    """
    And then define the study configuration (here we activate the smart
    caching module that will be able to remember which process has already been
    processed):
    """
    study_config = StudyConfig(
        modules=["SmartCachingConfig", "FSLConfig", "MatlabConfig",
                 "SPMConfig", "NipypeConfig"],
        use_smart_caching=True,
        fsl_config="/etc/fsl/4.1/fsl.sh",
        use_fsl=True,        
        output_directory=working_dir)

    # Create pipeline
    start_time = datetime.datetime.now()
    print "Start Pipeline Creation", start_time
    pipeline = get_process_instance("clindmri.preproc.fsl_preproc.xml")
    print "Done in {0} seconds.".format(datetime.datetime.now() - start_time)

    # View pipeline
    if 0:
        app = QtGui.QApplication(sys.argv)
        view1 = PipelineDevelopperView(pipeline)
        view1.show()
        app.exec_()
        del view1

    # Set pipeline input parameters
    pipeline.dfile = "/volatile/imagen/dmritest/000000022453/DTI/000000022453s011a1001.nii.gz"
    pipeline.bvalfile = "/volatile/imagen/dmritest/000000022453/DTI/000000022453s011a1001.bval"
    pipeline.bvecfile = "/volatile/imagen/dmritest/000000022453/DTI/000000022453s011a1001.bvec"
    print "Done in {0} seconds.".format(datetime.datetime.now() - start_time)

    #print pipeline.nodes["eddy"].process._nipype_interface.inputs
    print pipeline.nodes["eddy"].process._nipype_interface.cmdline

    # Execute the pipeline in the configured study
    study_config.run(pipeline, verbose=1)


if __name__ == "__main__":
    pilot_fsl_preproc()
