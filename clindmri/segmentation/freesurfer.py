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
from clindmri.extensions.fsl.exceptions import FSLResultError
from clindmri.extensions.fsl.exceptions import FreeSurferRuntimeError
from clindmri.extensions.freesurfer.wrappers import FSWrapper
from clindmri.extensions.freesurfer import read_cortex_surface_segmentation
from clindmri.registration.fsl import flirt
from clindmri.extensions.fsl import flirt2aff
from clindmri.registration.utils import extract_image


def recon_all(anatfile, output_directory, sid,
              fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
    """ Performs all the FreeSurfer cortical reconstruction process.

    <unit>
        <input name="anatfile" type="File" desc="The input anatomical image
            to be segmented with freesurfer."/>
        <input name="output_directory" type="Directory" description="The
            freesurfer destination folder."/>
        <input name="sid" type="Str" description="The current subject
            identifier."/>
        <input name="fsconfig" type="File" description="The freesurfer
            configuration batch."/>
        <output name="subjfsdir" type="Directory" description="Path to the
            resulting freesurfer segmentation."/>
    </unit>
    """
    # Create the fs output directory if necessary
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    # Call freesurfer
    cmd = ["recon-all", "-all", "-subjid", sid, "-i", anatfile,
           "-sd", output_directory]
    recon = FSWrapper(cmd, shfile=fsconfig)
    recon()
    if recon.exitcode != 0:
        raise FreeSurferRuntimeError(recon.cmd[0], " ".join(recon.cmd[1:]))
    subjfsdir = os.path.join(output_directory, sid)

    return subjfsdir


def cortex(t1_file, fsdir, outdir, dest_file=None, prefix="cortex",
           generate_mask=True, generate_seeds=True):
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
    prefix: str (optional, default 'cortex')
        the output files prefix.
    generate_mask: bool (optional, default True)
        if True generate a white matter binary mask.
    generate_seeds: boll (optional, default False)
        if True create a 'seeds' directory containing all the gyri mask as
        idenpendent files.

    Returns
    -------
    mask_file: str
        the white matter mask image file.
    label_file: str
        the gyri label image file.
    seeds: list of str
        a list with the seed volumes.
    """
    # Create the output directory if necessary
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Load the dataset
    t1_image = nibabel.load(t1_file)
    t1_affine = t1_image.get_affine()

    # If a destination file is specified register it to the t1
    if dest_file is not None:

        # Load dataset
        dest_image = nibabel.load(dest_file)
        dest_affine = dest_image.get_affine()
        dest_shape = dest_image.get_shape()

        # In case of temporal serie extract the first volume
        if len(dest_shape) > 3:
            temporal_dest_file = dest_file
            dest_file = os.path.join(outdir, prefix + "_volume-0.nii.gz")
            extract_image(temporal_dest_file, index=0, out_file=dest_file)
            dest_shape = dest_shape[:3]

        # Register destination image to t1 image
        trf_file = os.path.join(outdir, prefix + "_dest_to_t1.trf")
        reg_file = os.path.join(outdir, prefix + "_dest_to_t1.nii.gz")
        flirt(dest_file, t1_file, omat=trf_file, out=reg_file, usesqform=False,
              cost="normmi", dof=6)
        voxel_dest_to_t1 = flirt2aff(trf_file, dest_file, t1_file)
        voxel_t1_to_dest = numpy.linalg.inv(voxel_dest_to_t1)

    # Otherwise use identity transformation
    else:
        trf_file = None
        reg_file = None
        dest_affine = t1_affine
        dest_shape = t1_image.get_shape()
        voxel_t1_to_dest = numpy.identity(4)

    # Load the FreeSurfer surface in the 'dest_file' voxel coordinates or
    # 't1_file' coordinates if not specified
    t1_physical_to_voxel = numpy.linalg.inv(t1_affine)
    seg = read_cortex_surface_segmentation(fsdir, t1_physical_to_voxel,
                                           voxel_t1_to_dest)

    # Create a mask of the white matter of both hemisphere
    if generate_mask:
        mask_array = seg["lh"].voxelize(dest_shape)
        mask_array += seg["rh"].voxelize(dest_shape)

    # Create a gyri label image of both hemisphere
    label_array = {}
    try:
        label_array["lh"], shift_lh = seg["lh"].labelize(dest_shape)
        label_array["rh"], shift_rh = seg["rh"].labelize(dest_shape, shift_lh)
    except:
        if reg_file is not None:
            raise FSLResultError("flirt")
        raise

    # Create the seeds
    seeds = []
    if generate_seeds:
        seedsdir = os.path.join(outdir, "gyri")
        if not os.path.isdir(seedsdir):
            os.mkdir(seedsdir)
        for hemi in ["lh", "rh"]:
            surf = seg[hemi]
            hemi_label_array = label_array[hemi]
            seed_array = numpy.zeros(hemi_label_array.shape,
                                     dtype=hemi_label_array.dtype)
            for index, item in surf.metadata.items():
                if index != 0:
                    if hemi == "rh":
                        index += shift_lh
                    seed_array[numpy.where(hemi_label_array == index)] = 1
                    seed_file = os.path.join(
                        seedsdir, "{0}-{1}.nii.gz".format(hemi, item["region"]))
                    seed_image = nibabel.Nifti1Image(seed_array, dest_affine)
                    nibabel.save(seed_image, seed_file)
                    seed_array[...] = 0
                    seeds.append(seed_file)

    # Save the mask and label images
    mask_file = None
    if generate_mask:
        mask_file = os.path.join(outdir, prefix + "_mask.nii.gz")
        mask_image = nibabel.Nifti1Image(mask_array, dest_affine)
        nibabel.save(mask_image, mask_file)
    label_array = label_array["lh"] + label_array["rh"]     
    label_file = os.path.join(outdir, prefix + "_gyri_labels.nii.gz")
    label_image = nibabel.Nifti1Image(label_array, dest_affine)
    nibabel.save(label_image, label_file)

    return mask_file, label_file, seeds, reg_file, trf_file

