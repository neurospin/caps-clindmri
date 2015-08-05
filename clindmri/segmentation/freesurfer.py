#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import nibabel
import numpy

# Clindmri import
from clindmri.extensions.freesurfer import read_cortex_surface_segmentation
from clindmri.registration.fsl import flirt
from clindmri.extensions.fsl import flirt2aff
from clindmri.registration.utils import extract_image


def cortex(t1_file, fsdir, outdir, dest_file=None, prefix="cortex",
           mask=True):
    """ Compute a white matter mask and gyri labelization from the FreeSurfer
    'white' surface.

    Parameters
    ----------
    t1_file: str (mandatory)
        a file containing the t1 image used in FreeSurfer for the segmentation.
    fsdir: str( mandatory)
        the subject freesurfer segmentation directory.
    outdir: str (mandatory)
        the output directory.
    dest_file: str (optional, default None)
        a file containing an image where we want to project the segmentations:
        an affine transform is used to align this image to the t1 image.
    prefix: str (optional, default 'wm')
        the output files prefix.
    mask: bool (optional, default True)
        if True generate a white matter binary mask.

    Returns
    -------
    mask_file: str
        the white matter mask image file.
    label_file: str
        the gyri label image file.
    """
    # Load the dataset
    t1_image = nibabel.load(t1_file)
    t1_affine = t1_image.get_affine()

    # If a destination file is specified register it to the t1
    if dest_file is not None:

        # Load dataset
        basename = os.path.basename(dest_file).split(".")[0]
        dest_image = nibabel.load(dest_file)
        dest_affine = dest_image.get_affine()
        dest_shape = dest_image.get_shape()

        # In case of temporal serie extract the first volume
        if len(dest_shape) > 3:
            temporal_dest_file = dest_file
            dest_file = os.path.join(outdir, "e" + basename + "-0.nii.gz")
            extract_image(temporal_dest_file, index=0, out_file=dest_file)
            dest_shape = dest_shape[:3]

        # Register destination image to t1 image
        trf_file = os.path.join(outdir, "dest_to_t1_" + basename + ".trf")
        reg_file = os.path.join(outdir, "dest_to_t1_" + basename + ".nii.gz")
        flirt(dest_file, t1_file, omat=trf_file, out=reg_file, usesqform=True,
              cost="normmi", dof=12)
        voxel_dest_to_t1 = flirt2aff(trf_file, dest_file, t1_file)
        voxel_t1_to_dest = numpy.linalg.inv(voxel_dest_to_t1)

    # Otherwise use identity transformation
    else:
        basename = os.path.basename(t1_file).split(".")[0]
        dest_affine = t1_affine
        dest_shape = t1_image.get_shape()
        voxel_t1_to_dest = numpy.identity(4)

    # Load the FreeSurfer surface in the 'dest_file' voxel coordinates or
    # 't1_file' coordinates if not specified
    t1_physical_to_voxel = numpy.linalg.inv(t1_affine)
    seg = read_cortex_surface_segmentation(fsdir, t1_physical_to_voxel,
                                           voxel_t1_to_dest)

    # Create a mask of the white matter of both hemisphere
    if mask:
        mask_array = seg["lh"].voxelize(dest_shape)
        mask_array += seg["rh"].voxelize(dest_shape)

    # Create a gyri label image of both hemisphere
    label_array_lh, shift_lh = seg["lh"].labelize(dest_shape)
    label_array_rh, shift_rh = seg["rh"].labelize(dest_shape, shift_lh)
    label_array = label_array_lh + label_array_rh

    # Save the results
    mask_file = None
    if mask:
        mask_file = os.path.join(outdir, prefix + basename + "-mask.nii.gz")
        mask_image = nibabel.Nifti1Image(mask_array, dest_affine)
        nibabel.save(mask_image, mask_file)       
    label_file = os.path.join(outdir, prefix + basename + "-labels.nii.gz")
    label_image = nibabel.Nifti1Image(label_array, dest_affine)
    nibabel.save(label_image, label_file)

    return mask_file, label_file
