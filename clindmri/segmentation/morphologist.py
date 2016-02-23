##########################################################################
# NSAp - Copyright (C) CEA, 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
# 
# First activate brainvisa env:
# .bv_env
##########################################################################

# System import
import os
import time
import pickle
import shutil
import json

# Soma import
import soma.uuid

# Capsul import
from capsul.study_config.study_config import StudyConfig
from capsul.process import process_with_fom
from capsul.process import get_process_instance
from capsul.pipeline import pipeline_workflow

# Soma workflow import
from soma_workflow.client import Workflow
from soma_workflow.client import WorkflowController


def morphologist_all(t1file, sid, outdir, study="morphologist", waittime=10,
                     spmexec="/i2bm/local/spm8-standalone/run_spm8.sh",
                     spmdir="/i2bm/local/spm8-standalone"):
    """ Performs all the Morphologist steps.

    Steps:

    1- Ensure image orientation and reorient it if needed (Prepare Subject for
       Anatomical Pipeline).
    2- Computation of a brain mask (Brain Mask Segmentation).
    3- Computation of a mask for each hemisphere (Split Brain Mask).
    4- A grey/white classification of each hemisphere to perform "Voxel Based
       Morphometry" (Grey White Classification) and spherical triangulation of
       cortical hemispheres (Grey White Surface).
    5- Spherical triangulation of the external interface of the cortex of one or
       two hemispheres (Get Spherical Hemi Surface).
    6- Computation of a graph representing the cortical fold topography
       (Cortical Fold Graph).
    7- Automatic identification of the cortical sulci (Automatic Sulci
       Recognition), located in the "sulci" toolbox.

    The execution is performed with soma_workflow that has to be installed in
    the bv_env environment.

    To check the worklow submission, use the 'soma_workflow_gui' command.

    If the input 't1file' has no the expected extension, an Exception will
    be raised.
    If the $outdir/$study/$sid has already been created, an Exception will
    be raised.

    Parameters
    ----------
    t1file: str (mandatory)
        the path to a ".nii.gz" anatomical T1 weighted file.
    sid: str (mandatory)
        a subject identifier.
    outdir: str (mandatory)
        the morphologist output files will be written in $outdir/$study/$sid.
    study: str (mandatory)
        the name of the study.
    waittime: float (mandatory)
        a delay (in seconds) used to check the worflow status.
    spmexec: str (mandatory)
        the path to the standalone SPM execution file.
    spmdir: str (mandatory)
        the standalone SPM directory.

    Returns
    -------
    wffile: str
        a file containing the submitted workflow.
    wfid: int
        the submitted workflow identifier.
    wfstatus: str
        the submited worflow status afer 'waittime' seconds.
    """
    # Check roughly the input file extension
    if not t1file.endswith(".nii.gz"):
        raise Exception("'{0}' is not a COMPRESSED NIFTI file.".format(t1file))
    
    # Create a configuration for the morphologist study
    study_config = StudyConfig(
        modules=StudyConfig.default_modules + ["FomConfig", "BrainVISAConfig"])
    study_dict = {
        "name" : "morphologist_fom",
        "input_directory" : outdir,
        "output_directory" : outdir,
        "input_fom" : "morphologist-auto-nonoverlap-1.0",
        "output_fom" : "morphologist-auto-nonoverlap-1.0",
        "shared_fom" : "shared-brainvisa-1.0",
        "spm_directory" : spmdir,
        "use_soma_workflow" : True,
        "use_fom" : True,
        "spm_standalone" : True,
        "use_matlab" : False,
        "volumes_format" : "NIFTI gz",
        "meshes_format" : "GIFTI",
        "use_spm" : True,
        "spm_exec" : spmexec,
        "study_config.somaworkflow_computing_resource": "localhost",
        "somaworkflow_computing_resources_config": {
            "localhost": {
            }
        }
    }
    study_config.set_study_configuration(study_dict)

    # Create the morphologist pipeline
    pipeline = get_process_instance(
        "morphologist.capsul.morphologist.Morphologist")
    morphologist_pipeline = process_with_fom.ProcessWithFom(
        pipeline, study_config)
    morphologist_pipeline.attributes = dict(
        (trait_name, getattr(morphologist_pipeline, trait_name))
        for trait_name in morphologist_pipeline.user_traits())
    morphologist_pipeline.attributes["center"] = "morphologist"
    morphologist_pipeline.attributes["subject"] = sid
    morphologist_pipeline.create_completion()

    # Create a worflow from the morphologist pipeline
    workflow = Workflow(name="{0} {1}".format(study, sid),
                        jobs=[])
    workflow.root_group = []

    # Create morphologist expected tree
    # ToDo: use ImportT1 from axon
    subjectdir = os.path.join(outdir, study, sid)
    if os.path.isdir(subjectdir):
        raise Exception("Folder '{0}' already created.".format(subjectdir))
    os.makedirs(os.path.join(
        subjectdir, "t1mri", "default_acquisition",
        "default_analysis", "folds", "3.1", "default_session_auto"))
    os.makedirs(os.path.join(
        subjectdir, "t1mri", "default_acquisition",
        "registration"))
    os.makedirs(os.path.join(
        subjectdir, "t1mri", "default_acquisition",
        "segmentation", "mesh"))

    # Copy T1 file in the morphologist expected location
    destfile = os.path.join(subjectdir, "t1mri",
                            "default_acquisition", sid + ".nii.gz")
    shutil.copy(t1file, destfile)

    # Create source_referential morphologist expected file
    source_referential = {"uuid": str(soma.uuid.Uuid())}
    referential_file = os.path.join(
        subjectdir, "t1mri", "default_acquisition", "registration",
        "RawT1-{0}_default_acquisition.referential".format(sid))
    attributes = "attributes = {0}".format(json.dumps(source_referential))
    with open(referential_file, "w") as openfile:
        openfile.write(attributes)

    # Create the workflow
    wf = pipeline_workflow.workflow_from_pipeline(
        morphologist_pipeline.process, study_config=study_config)
    workflow.add_workflow(wf, as_group="{0}_{1}".format(study, sid))
    wffile = os.path.join(subjectdir, "{0}.wf".format(study))
    pickle.dump(workflow, open(wffile, "w"))
   
    # Execute workflow
    controller = WorkflowController()
    wfid = controller.submit_workflow(
        workflow=workflow, name="{0}_{1}".format(study, sid))

    # Return the worflow status after execution
    while True:
        time.sleep(waittime)
        wfstatus = controller.workflow_status(wfid)
        if wfstatus not in ["worklflow_not_started", "workflow_in_progress"]:
            break

    return wffile, wfid, wfstatus
