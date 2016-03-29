##########################################################################
# NSAp - Copyright (C) CEA, 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
Mocking Popen directly - need to construct a Mock to return, and adjust its
communicate() return_value.
The benefit of this approach is in not needing to do the strip/split on your
fake return string.
"""

# System import
import unittest
import mock
from mock import patch

# Clindmri import
from clindmri.segmentation.fsl import bet2


class FSLSegmentation(unittest.TestCase):
    """ Test the 'segmentation.fsl' module.
    """
    mock_fslenv = {
        "FSLOUTPUTTYPE": "NIFTI_GZ"
    }

    def setUp(self):
        """ Run before each test - the mock_popen will be available and in the
        right state in every test<something> function.
        """
        self.popen_patcher = patch(
            "clindmri.extensions.fsl.wrappers.subprocess.Popen")
        self.mock_popen = self.popen_patcher.start()
        self.mock_process = mock.Mock()
        attrs = {
            "communicate.return_value": ("mock_OK", "mock_NONE"),
            "returncode": 0
        }
        self.mock_process.configure_mock(**attrs)
        self.mock_popen.return_value = self.mock_process

    def tearDown(self):
        """ Run after each test.
        """
        self.popen_patcher.stop()

    @mock.patch("clindmri.extensions.fsl.wrappers.FSLWrapper._environment")
    def test_bet2(self, mock_setenv):
        """ Test the FSL 'bet2' command wrapping.
        """
        expected_files = (
            "/my/path/mock_output.nii.gz",
            "/my/path/mock_output_mask.nii.gz",
            "/my/path/mock_output_mesh.vtk",
            None,
            "/my/path/mock_output_inskull_mask.nii.gz",
            "/my/path/mock_output_inskull_mesh.nii.gz",
            "/my/path/mock_output_outskull_mask.nii.gz",
            "/my/path/mock_output_outskull_mesh.nii.gz",
            "/my/path/mock_output_outskin_mask.nii.gz",
            "/my/path/mock_output_outskin_mesh.nii.gz",
            "/my/path/mock_output_skull_mask.nii.gz")
        mock_setenv.return_value = self.mock_fslenv
        mock_input = "/my/path/mock_input.nii.gz"
        mock_output = "/my/path/mock_output"
        output_files = bet2(
            mock_input, mock_output, f=0.5, m=True, shfile="mock.sh", s=True,
            e=True)
        self.assertTrue(self.mock_popen.called)
        self.assertTrue(output_files == expected_files)


if __name__ == "__main__":
    """ The unittest runner looks for any functions called test<something> in
    classes that extend FSLSegmentation so you can add as many test<something>
    functions to the classes above as you like, and they'll all get run.
    """
    unittest.main()
