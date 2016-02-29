#########################################################################
# NSAp - Copyright (C) CEA, 2015 - 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

# Clindmri import
from clindmri.extensions.connectomist import ConnectomistWrapper


def dwi_outlier_detection(
        outdir,
        raw_dwi_dir,
        rough_mask_dir,
        nb_tries=10,
        path_connectomist=(
            "/i2bm/local/Ubuntu-14.04-x86_64/ptk/bin/connectomist")):
    """
    Wrapper to Connectomist's "Outliers" tab.

    Parameters
    ----------
    outdir: str
        path to Connectomist output work directory.
    raw_dwi_dir: str
        path to Connectomist Raw DWI folder.
    rough_mask_dir: str
        path to Connectomist Rough Mask folder.
    nb_tries: int (optional, default 10)
        nb of times to try an algorithm if it fails.
        It often crashes when running in parallel. The reason
        why it crashes is unknown.
    path_connectomist: str (optional)
        path to the Connectomist executable.

    Returns
    -------
    outdir: str
        path to Connectomist's output directory.

    <unit>
        <output name="outliers_dir"  type="Directory" />

        <input name="outdir"         type="Directory" />
        <input name="raw_dwi_dir"    type="Directory" />
        <input name="rough_mask_dir" type="Directory" />
    </unit>
    """
    # Dict with all parameters for connectomist
    algorithm = "DWI-Outlier-Detection"
    parameters_dict = {'rawDwiDirectory':        raw_dwi_dir,
                       'roughMaskDirectory':  rough_mask_dir,
                       'outputWorkDirectory':         outdir,
                       '_subjectName':                    '',
                       'discardedOrientationList':        '',
                       'outlierFactor':                  3.0}

    # Call with Connectomist
    connprocess = ConnectomistWrapper(path_connectomist)
    parameter_file = ConnectomistWrapper.create_parameter_file(
        algorithm, parameters_dict, outdir)
    connprocess(algorithm, parameter_file, outdir, nb_tries=nb_tries)

    return outdir
