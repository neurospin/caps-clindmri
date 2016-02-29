##########################################################################
# NSAp - Copyright (C) CEA, 2015 - 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

"""
Package to wrap the Connectomist software and simplify scripting calls.

Each preprocessing step (i.e. one Connectomist tab) can be run through the use
of a dedicated function of the package.

All the preprocessing steps can be done at once using the
complete_preprocessing() function.

The package also allows parallelization of multiple complete preprocessing
using the parallel_preprocessing() function.
"""

from clindmri.extensions.connectomist.manufacturers import MANUFACTURERS
from clindmri.extensions.connectomist.exceptions import (
    ConnectomistError, ConnectomistBadManufacturerNameError,
    ConnectomistMissingParametersError, ConnectomistBadFileError)

# Wrappers of Connectomist's tabs
from .import_and_qspace_model import dwi_data_import_and_qspace_sampling
from .mask                    import dwi_rough_mask_extraction
from .outliers                import dwi_outlier_detection
from .susceptibility          import dwi_susceptibility_artifact_correction
from .eddy_current_and_motion import (dwi_eddy_current_and_motion_correction,
                                      export_eddy_motion_results_to_nifti)
from .complete_preprocessing import complete_preprocessing
from .utils import (ptk_gis_to_nifti,
                    ptk_nifti_to_gis,
                    ptk_concatenate_volumes,
                    ptk_split_t2_and_diffusion)
