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

# Clindmri imports
from clindmri.tractography.fsl import probtrackx2
from clindmri.registration.fsl import flirt
from clindmri.segmentation.freesurfer import mri_vol2surf


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
        The FreeSurfer '.sh' config file.

    Returns
    -------
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

    return textures


def qc_profile(hemi, mri_vol2surf_out, nodif_file, overlay_file, level,
               fs_subjects_dir, fs_subject, output_dir, fsconfig):
    """ Tractography QC function.
        Generates views of :
        - the superposition of the nodif image with tractography result volume.
        - the connected points on the cortical surface
        Resample cortical meshes if needed.
        Results output are available as gif and png.

        Parameters
        ----------
        hemi: str (mandatory)
            hemisphere (lh or rh).
        mri_vol2surf_out: str (mandatory)
            The FS mri_vol2surf mgz volume, or connection array.
        nodif_file: str (mandatory)
            The nodif image file.
        overlay_file: str (mandatory)
            The protrackx2 output volume.
        level: int (mandatory)
            order of icosahedron when FS option trgsubject=ico.
            specifies the order of the icosahedral tesselation. the number of
            vertices is given by the formula 10*2^n+2. In general, it is best
            to use the largest size available.
        fs_subjects_dir: str (mandatory)
            Freesurfer subjects directory SUBJECTS_DIR.
        fs_subject: str (mandatory)
            Freesurfer subject name.
        output_dir: str (mandatory)
            The QC output directory.
        fsconfig: str (mandatory)
            The Freesurfer .sh config file.

        """
    subject_dir = os.path.join(fs_subjects_dir, fs_subject)

    # Define QC images directories
    png_outdir = os.path.join(output_dir, 'PNGs')
    gif_outdir = os.path.join(output_dir, 'GIFs')
    if not os.path.isdir(png_outdir):
        os.mkdir(png_outdir)
    if not os.path.isdir(gif_outdir):
        os.mkdir(gif_outdir)

    # Generate gif of the fdt_paths overlay on the nodif image if needed
    assert overlay_file.endswith('.nii.gz')
    overlay_name = os.path.basename(overlay_file)[:-len('.nii.gz')]
    prefix = 'nodif_' + overlay_name
    output_gif_file = os.path.join(gif_outdir, prefix + '.gif')
    if not os.path.isfile(output_gif_file):

        # Get the fdt_paths image shape
        overlay_shape = nibabel.load(overlay_file).shape

        # Superpose the nodif and fdt_paths volumes with nilearn
        nodif_and_overlay = plotting.plot_anat(nodif_file,
                                               cut_coords=overlay_shape[2],
                                               display_mode='z')
        # Choose a colormap
        cmap = cm.get_cmap('Spectral')

        # Add fdt_paths overlay (from probtrackx2) on the nodif image
        nodif_and_overlay.add_overlay(overlay_file, cmap=cmap)

        # Save all the z cuts of the nodif and overlay superposition
        all_z_cuts_png = os.path.join(png_outdir,
                                      prefix + '.png')
        nodif_and_overlay.savefig(all_z_cuts_png)

        # Get the large image dimensions
        width, height = get_image_dimensions(all_z_cuts_png)

        # Split the large image generated by nilearn into multiple layers to
        # build the gif file
        z_cuts = split_image(all_z_cuts_png, png_outdir, prefix,
                             width/overlay_shape[2], height)

        # Combine the png layers to obtain the gif file
        images_to_gif(
            z_cuts,
            output_gif_file,
            delay=10
        )

    # Resample the white mesh on the icosphere if needed
    meshfile = os.path.join(subject_dir, "convert", "{0}.white.{1}"
                            .format(hemi, level))
    annotfile = os.path.join(subject_dir, "convert",
                             "{0}.aparc.annot.{1}".format(hemi, level))
    if (not os.path.isfile(meshfile)) or (not os.path.isfile(annotfile)):
        resample_cortical_surface(fs_subjects_dir, output_dir, orders=[level],
                                  surface_name="white", fsconfig=fsconfig)

    # Check that the input volume has the expected properties
    assert(mri_vol2surf_out.endswith('.mgz'))

    connection_array = nibabel.load(mri_vol2surf_out).get_data()

    connection_array_dim = connection_array.ndim
    connection_array_shape = connection_array.shape
    if connection_array_dim != 3:
        raise ValueError("Expected connection array dim : 3. Found : {0}"
                         .format(connection_array_dim))
    if (connection_array_shape[1] != 1) or (connection_array_shape[2] != 1):
        raise ValueError("Expected connection array shape : (*, 1, 1). "
                         "Found : {0}".format(connection_array_shape))

    # Flatten the connection array into 1D
    labels = connection_array.ravel()

    # Load the vtk mesh
    surface = TriSurface.load(meshfile, annotfile=annotfile)

    # Define a vtk renderer + actor and add the labels on the vtk surface
    ren = pvtk.ren()
    actor = pvtk.surface(surface.vertices, surface.triangles, labels)
    pvtk.add(ren, actor)

    mri_vol2surf_basename = os.path.splitext(
        os.path.basename(mri_vol2surf_out))[0]

    # Generate the png layers to be combined into gif
    rotatez_ang = 0
    snaps = []
    for _ in range(2):
        if rotatez_ang > 0:
            # Rotate the vtk actor along the Z axis
            ren.GetActors().GetLastActor().RotateZ(rotatez_ang)
        snaps += pvtk.record(ren,
                             png_outdir,
                             mri_vol2surf_basename + '_{0}'.format(rotatez_ang),
                             n_frames=36,
                             az_ang=10)
        rotatez_ang += 90
    snaps = sorted(snaps)

    # Combine the png layers to obtain the gif file
    images_to_gif(
        snaps,
        os.path.join(gif_outdir, mri_vol2surf_basename + '.gif'),
        delay=10
    )

def get_image_dimensions(input_image):
    """ Get image dimensions using ImageMagick identify.

        Parameters
        ----------
        input_image: str (mandatory)
            Input image.

        Returns
        -------
        width, height: (int, int)
            The image dimensions.
    """
    cmd = ["identify", "-ping", "-format", "'%wx%h'", input_image]
    # Execute the ImageMagick command
    magick = ImageMagickWrapper(cmd)
    dim_str = magick()
    width, height = (int(dimension) for dimension
                     in re.search('(\d+)x(\d+)', dim_str).groups())
    return width, height


def split_image(input_image, output_dir, output_prefix, output_width,
                output_height):
    """ Split one image into multiple images using ImageMagick convert.

        Parameters
        ----------
        input_image: str (mandatory)
            Path to the image.
        output_dir: str (mandatory)
            The output directory.
        output_prefix: str (mandatory)
            The output images prefix.
        output_width: int (mandatory)
            The width in pixels of the output images.
        output_height: int (mandatory)
            The height in pixels of the output images.

        Returns
        -------
        output_img_list: list of str
            List of the output images.
    """
    input_extension = os.path.splitext(os.path.basename(input_image))[-1]
    output_pattern = os.path.join(output_dir, "{0}-%03d{1}"
                                  .format(output_prefix, input_extension))
    cmd = ["convert", input_image,
           "-crop", "{0}x{1}".format(output_width, output_height),
           output_pattern]
    # Execute the ImageMagick command
    magick = ImageMagickWrapper(cmd)
    magick()
    return sorted(glob.glob(output_pattern.replace("%03d", "*")))
