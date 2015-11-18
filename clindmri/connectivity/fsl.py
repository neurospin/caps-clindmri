#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013-2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import nibabel

# Clindmri imports
from clindmri.tractography.fsl import probtrackx2
from clindmri.registration.fsl import flirt
from clindmri.segmentation.freesurfer import mri_vol2surf
from clindmri.extensions.freesurfer.reader import TriSurface


def get_profile(ico_order, nodif_file, nodifmask_file, seed_file,
                bedpostx_samples, outdir, t1_file, trf_file, dat_file, fsdir,
                sid, fsconfig):
    """ Probabilistic profile

    Computes the tractography using FSL probtrackx2 and projects the result
    on the cortical surface using FS mri_vol2surf.

    Parameters
    ----------
    ico_order: int (mandatory)
        icosahedron order in [0, 7] that will be used to generate the cortical
        surface texture at a specific tessalation (the corresponding cortical
        surface can be resampled using the
        'clindmri.segmentation.freesurfer.resample_cortical_surface' function).
    nodif_file: str (mandatory)
        file for probtrackx2 containing the no diffusion volume and associated
        space information.
    nodifmask_file: str (mandatory)
        file for probtrackx2 containing the tractography mask (ie., a mask of
        the white matter).
    seed_file: str (mandatory)
        text file for probtrackx2 containing seed coordinates.
    bedpostx_samples: str (mandatory)
        path prefix for bedpostX model samples files injected in probtrackx2
        (eg., fsl.bedpostX/merged).
    outdir: str (mandatory)
        output directory.
    t1_file : str (mandatory)
        T1 image file used to align the produced probabilitic tractography map
        to the T1 space.
    trf_file : str (mandatory)
        diffusion to t1 space affine transformation matrix file.
    dat_file: str (mandatory)
        structural to FreeSurfer space affine transformation matrix '.dat'
        file as computed by 'tkregister2'.
    fsdir: str (mandatory)
        FreeSurfer subjects directory 'SUBJECTS_DIR'.
    sid: str (mandatory)
        FreeSurfer subject identifier.
    fsconfig: str (mandatory)
        the FreeSurfer '.sh' config file.

    Returns
    -------
    proba_file: str
        the seed probabilistic tractography volume.
    textures: dict
        a dictionary containing the probabilist texture for each hemisphere.
    """
    # Generates the diffusion probability map
    proba_files, _ = probtrackx2(
        simple=True, seedref=nodif_file, out="fdt_paths", seed=seed_file,
        loopcheck=True, onewaycondition=True, samples=bedpostx_samples,
        mask=nodifmask_file, dir=outdir)

    # Check that only one 'fdt_paths' has been generated
    if len(proba_files) != 1:
        raise Exception("One probabilistic tractography file expected at this "
                        "point: {0}".format(proba_files))
    proba_file = proba_files[0]
    proba_fname = os.path.basename(proba_file).replace(".nii.gz", "")

    # Apply 'trf_file' affine transformation matrix using FSL flirt function:
    # probability map (diffusion space) -> probability map (T1 space).
    flirt_t1_file = os.path.join(outdir, proba_fname + "_t1_flirt.nii.gz")
    flirt(proba_file, t1_file, out=flirt_t1_file, applyxfm=True, init=trf_file)

    # Project the volumic probability map (T1 space) generated with FSL flirt
    # on the cortical surface (Freesurfer space) (both hemispheres) using
    # Freesurfer's mri_vol2surf and applying the 'dat_file' transformation.
    textures = {}
    for hemi in ["lh", "rh"]:
        prob_texture_file = os.path.join(
            outdir, "{0}.{1}_vol2surf.mgz".format(hemi, proba_fname))
        mri_vol2surf(hemi, flirt_t1_file, prob_texture_file, ico_order,
                     dat_file, fsdir, sid, surface_name="white",
                     fsconfig=fsconfig)
        textures[hemi] = prob_texture_file

    return proba_file, textures


def qc_profile(nodif_file, proba_file, proba_texture,  ico_order,
               fsdir, sid, outdir, fsconfig, actor_ang=(0., 0., 0.)):
    """ Connectivity profile QC.

    Generates views of:
    - the superposition of the nodif image with tractography result volume.
    - the connected points on the cortical surface
    Resample cortical meshes if needed.
    Results output are available as gif and png.

    Parameters
    ----------
    nodif_file: str (mandatory)
        file for probtrackx2 containing the no diffusion volume and associated
        space information.
    proba_file: str (mandatory)
        the protrackx2 output seeding probabilistic path volume.
    proba_texture: dict (mandatory)
        the FreeSurfer mri_vol2surf '.mgz' 'lh' and 'rh' textrue that contains
        the cortiacal connection strength.
    ico_order: int (mandatory)
        icosahedron order in [0, 7] that will be used to generate the cortical
        surface texture at a specific tessalation (the corresponding cortical
        surface can be resampled using the
        'clindmri.segmentation.freesurfer.resample_cortical_surface' function).
    fsdir: str (mandatory)
        FreeSurfer subjects directory 'SUBJECTS_DIR'.
    sid: str (mandatory)
        FreeSurfer subject identifier.
    outdir: str (mandatory)
        The QC output directory.
    fsconfig: str (mandatory)
        the FreeSurfer '.sh' config file.
    actor_ang: 3-uplet (optinal, default (0, 0, 0))
        the actor x, y, z position (in degrees).

    Returns
    -------
    snaps: list of str
        two gifs images, one showing the connection profile as a texture on
        the cortical surface, the other a volumic representation of the
        deterministic tractography.
    """
    import clindmri.plot.pvtk as pvtk
    from clindmri.plot.slicer import animate_image

    # Construct/check the subject directory
    subjectdir = os.path.join(fsdir, sid)
    if not os.path.isdir(subjectdir):
        raise ValueError(
            "'{0}' is not a FreeSurfer subject directory.".format(subjectdir))

    # Check that the output QC directory exists
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Superpose the nodif and probabilistic tractography volumes
    proba_shape = nibabel.load(proba_file).shape
    snaps = []
    snaps.append(
        animate_image(nodif_file, overlay_file=proba_file, clean=True,
                      overlay_cmap="Spectral", cut_coords=proba_shape[2],
                      outdir=outdir))

    # Define a renderer
    ren = pvtk.ren()

    # For each hemisphere
    for hemi in ["lh", "rh"]:

        # Get the the white mesh on the desired icosphere
        meshfile = os.path.join(
            subjectdir, "convert", "{0}.white.{1}.native".format(
                hemi, ico_order))
        if not os.path.isfile(meshfile):
            raise ValueError(
                "'{0}' is not a valid white mesh. Generate it through the "
                "'clindmri.scripts.freesurfer_conversion' script.".format(
                    meshfile))

        # Check texture has the expected extension, size
        texture_file = proba_texture[hemi]
        if not texture_file.endswith(".mgz"):
            raise ValueError("'{0}' is not a '.mgz' file. Format not "
                             "supported.".format(texture_file))
        profile_array = nibabel.load(texture_file).get_data()
        profile_dim = profile_array.ndim
        profile_shape = profile_array.shape
        if profile_dim != 3:
            raise ValueError(
                "Expected profile texture array of dimension 3 not "
                "'{0}'".format(profile_dim))
        if (profile_shape[1] != 1) or (profile_shape[2] != 1):
            raise ValueError(
                "Expected profile texture array of shape (*, 1, 1) not "
                "'{0}'.".format(profile_shape))

        # Flatten the profile texture array
        texture = profile_array.ravel()

        # Load the white mesh
        surface = TriSurface.load(meshfile)

        # Define a textured surface actor
        actor = pvtk.surface(surface.vertices, surface.triangles, texture)
        actor.RotateX(actor_ang[0])
        actor.RotateY(actor_ang[1])
        actor.RotateZ(actor_ang[2])
        pvtk.add(ren, actor)

    # Create a animaton with the generated surface
    qcname = "profile_as_texture"
    snaps.extend(
        pvtk.record(ren, outdir, qcname, n_frames=36, az_ang=10, animate=True,
                    delay=10))

    return snaps
