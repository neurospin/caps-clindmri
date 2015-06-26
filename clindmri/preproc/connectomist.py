#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import subprocess
import numpy


MANUFACTURERS = {
    "Philips HealthCare": 2,
}


def dwi_data_import_and_qspace_sampling(dfile, bvecfile, bvalfile,
                                        manufacturer, invert_axis=(0, 0, 0), 
                                        output_directory):
    """ Bla

    <process>
        <return name="reobvecfile" type="File" desc="The reoriented b-vectors."/>
        <input name="bvecfile" type="File" desc="The diffusion b-vectors file."/>
        <input name="trffile" type="File" desc="The transformation parameters
            used to align each volume of the diffusion sequence."/>
        <input name="output_directory" type="Directory" desc="The destination
            folder."/>
    </process>
    """
    # Set default parameters
    algorithm_name = "DWI-Data-Import-And-QSpace-Sampling"
    parameter_values = {
        "qSpaceChoice3BValue": 1000, 
        "invertXAxis": 2, 
        "qSpaceChoice2NumberOfOrientations": 6, 
        "qSpaceChoice7NumberOfOrientations": 6, 
        "qSpaceChoice10NumberOfOrientations": "", 
        "diffusionTime": 1.0, 
        "fileNameDwi": "/volatile/imagen/AIMS_Philips/Analysis/01-InputData/dwi.ima", 
        "qSpaceTransform_yz": 0.0, 
        "qSpaceTransform_yy": 1.0, 
        "qSpaceTransform_yx": 0.0, 
        "qSpaceChoice8BValues": "", 
        "qSpaceChoice12BValues": "", 
        "qSpaceChoice12NumberOfOrientations": "", 
        "qSpaceChoice9BValues": "", 
        "qSpaceSamplingType": 4, 
        "invertZAxis": 0, 
        "qSpaceChoice9OrientationFileNames": "", 
        "qSpaceChoice10BValues": "", 
        "qSpaceChoice6NumberOfOrientations": 6, 
        "qSpaceChoice4BValue": 1000, 
        "qSpaceTransform_xx": 1.0, 
        "qSpaceTransform_xy": 0.0, 
        "qSpaceTransform_xz": 0.0, 
        "qSpaceChoice4NumberOfOrientations": 6, 
        "qSpaceChoice11NumberOfOrientations": "", 
        "qSpaceChoice3NumberOfOrientations": 6, 
        "numberOfT2": 1, 
        "_subjectName": "", 
        "sliceAxis": 2, 
        "qSpaceChoice13BValues": "", 
        "outputWorkDirectory": "/volatile/imagen/AIMS_Philips/Analysis/02-DataImportAndQSpace", 
        "qSpaceChoice1MaximumBValue": 1300, 
        "qSpaceChoice2BValue": 1000, 
        "qSpaceChoice13OrientationFileNames": "", 
        "qSpaceChoice5BValue": 1300, 
        "qSpaceChoice7BValues": "", 
        "numberOfRepetitions": 1, 
        "qSpaceChoice6BValues": "", 
        "qSpaceTransform_zz": 1.0, 
        "qSpaceChoice11BValues": "", 
        "flipAlongZ": 0, 
        "flipAlongX": 0, 
        "flipAlongY": 0, 
        "manufacturer": 2, 
        "qSpaceChoice8NumberOfOrientations": 6, 
        "qSpaceChoice1NumberOfSteps": 11, 
        "invertYAxis": 0, 
        "phaseAxis": 1, 
        "qSpaceChoice5OrientationFileNames": "/volatile/imagen/AIMS_Philips/Analysis/01-InputData/dwi.bvec", 
        "qSpaceTransform_zx": 0.0, 
        "qSpaceTransform_zy": 0.0, 
        "numberOfDiscarded": 0
    }
    # ToDo: check manufacturer in [Philips HealthCare ...]
    if manufacturer not in MANUFACTURERS:
        raise
    parameter_values["manufacturer"] = MANUFACTURERS[manufacturer]
    
    # Get all b-values
    # ToDo: check we have a b0 and a bval
    bvals = numpy.loadtxt(bvalfile)
    bvalsset = set(bvals.tolist()) - {0}

    # Set sampling
    # Case: spherical single-shell custom
    if len(bvalsset) == 1:
        parameter_values["qSpaceSamplingType"] = BVALS[1]

    # Execute process
    # ToDo write DWI-Data-Import-And-QSpace-Sampling.py in output_directory
    cmd = ["connectomist", "-p", ""]
    subprocess.check_call(cmd)

    return output_directory, ....
