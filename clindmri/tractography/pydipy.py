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
import numpy

# Clindmri import
from clindmri.tractography import savetxt

# Dipy import
from dipy.core.gradients import gradient_table
from dipy.io.gradients import read_bvals_bvecs
from dipy.tracking.eudx import EuDX
from dipy.reconst import peaks, shm
from dipy.tracking import utils


def deterministic(diffusion_file, bvecs_file, bvals_file, outdir,
                  mask_file=None, order=4, nb_seeds_per_voxel=1, step=0.5,
                  fmt="%.4f"):
    """ Compute a deterministic tractography using an ODF model.

    Parameters
    ----------
    diffusion_file: str (mandatory)
        a file containing the preprocessed diffusion data.
    bvecs_file: str (mandatory)
        a file containing the diffusion gradient directions.
    bvals_file: str (mandatory)
        a file containing the diffusion b-values.
    outdir: str (mandatory)
        the output directory.
    mask_file: str (optional, default None)
        an image used to mask the diffusion data during the tractography. If
        not set, all the image voxels are considered.
    order: int (optional, default 4)
        the order of the ODF model.
    nb_seeds_per_voxel: int (optional, default 1)
        the number of seeds per voxel used during the propagation.
    step: float (optional, default 0.5)
        the integration step in voxel fraction used during the propagation.
    fmt: str (optional, default '%.4f')
        the saved track elements format.

    Returns
    -------
    track_file: str
        a determinist model of the white matter organization.
    """
    # Read diffusion sequence
    bvals, bvecs = read_bvals_bvecs(bvals_file, bvecs_file)
    gtab = gradient_table(bvals, bvecs)
    diffusion_array = nibabel.load(diffusion_file).get_data()
    if mask_file is not None:
        mask_array = nibabel.load(mask_file).get_data()
    else:
        mask_array = numpy.ones(diffusion_array.shape[:3], dtype=numpy.uint8)

    # Estimate ODF model
    csamodel = shm.CsaOdfModel(gtab, order)
    csapeaks = peaks.peaks_from_model(
        model=csamodel, data=diffusion_array, sphere=peaks.default_sphere,
        relative_peak_threshold=.8, min_separation_angle=45,
        mask=mask_array)

    # Compute deterministic tractography in voxel space so affine is equal
    # to identity
    seeds = utils.seeds_from_mask(mask_array, density=nb_seeds_per_voxel)
    streamline_generator = EuDX(
        csapeaks.peak_values, csapeaks.peak_indices,
        odf_vertices=peaks.default_sphere.vertices, a_low=.05, step_sz=step,
        seeds=seeds)
    affine = streamline_generator.affine
    streamlines = list(streamline_generator)

    # Save the tracks
    track_file = os.path.join(outdir, "fibers.txt")
    savetxt(track_file, streamlines, fmt=fmt)

    return track_file

    
