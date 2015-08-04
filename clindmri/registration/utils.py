#! /usr/bin/env python
##########################################################################
# NSAP - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import nibabel


def extract_image(in_file, index, out_file=None):
    """ Extract the image at 'index' position.

    Parameters
    ----------
    in_file: str (mandatory)
        the input image.
    index: int (mandatory)
        the index of last image dimention to extract.
    out_file: str (optional, default None)
        the name of the extracted image file.

    Returns
    -------
    out_file: str
        the name of the extracted image file.
    """
    # Set default output if necessary
    dirname = os.path.dirname(in_file)
    basename = os.path.basename(in_file).split(".")[0]
    if out_file is None:
        out_file = os.path.join(
            dirname, "extract{0}_{1}.nii.gz".format(index, basename))

    # Extract the image of interest
    image = nibabel.load(in_file)
    affine = image.get_affine()
    extracted_array = image.get_data()[..., index]
    extracted_image = nibabel.Nifti1Image(extracted_array, affine)
    nibabel.save(extracted_image, out_file)
    
    return out_file
