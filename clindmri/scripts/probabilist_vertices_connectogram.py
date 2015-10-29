#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" A script for cortical connectivity using FreeSurfer cortical surface
vertices and FSL probabilistic tractography.

We will describe a method for performing probabilistic tractography in a single
subject using automatically generated cortical surface vertices from
Freesurfer resulting in a probabilistic vertices connectivity matrix.

Prerequisites
=============

This analysis requires the following MR data:

* A T1 image: this should cover the entire brain at high-resolution
  (~1mm3 voxels).
* The calculated FA image.
* The nodif image.
* The seedref file for probtrackx2 tractography
* The samples file for probtrackx2 tractography

Software
========

This analysis is using Freesurfer & FSL, all of them beeing freely
available.

List options
========

>> python probabilist_vertices_connectogram.py --help

Example
========

>> python probabilist_vertices_connectogram.py
          --sd freesurfer_results/
          --srcsubject subject1
          --t1 probaconnect/t1-1mm-1-001.nii.gz
          --fa probaconnect/dtifit/dtifit_FA.nii.gz
          --nodif probaconnect/nodif.nii.gz
          --seedref probaconnect/bedpostx.bedpostX/nodif_brain_mask
          --samples probaconnect/bedpostx.bedpostX/merged
          --outputdir probaconnect/CONNECTIVITY_RESULTS
"""

# system imports
from __future__ import print_function
import sys
sys.path.append('/home/rc244162/git/caps-clindmri')
from nilearn import plotting
from matplotlib import cm
import numpy
import nibabel
import glob
import re
import subprocess
import os
import argparse
# clindmri FS imports
from clindmri.plot import pvtk
from clindmri.extensions.freesurfer import read_cortex_surface_segmentation
from clindmri.extensions.freesurfer.reader import TriSurface
from clindmri.extensions.freesurfer.exceptions import FreeSurferRuntimeError
from clindmri.extensions.freesurfer.wrappers import FSWrapper
from clindmri.segmentation.freesurfer import (conformed_to_native_space,
                                              resample_cortical_surface)
# clindmri FSL imports
from clindmri.tractography.fsl import probtrackx2
from clindmri.extensions.fsl import flirt2aff
from clindmri.registration.fsl import flirt


def get_parser():
    """ Create the arguments parser.
    """
    parser = argparse.ArgumentParser(description='NxN connectogram generator.')

    # Required arguments
    requirednamed = parser.add_argument_group('required named arguments')
    requirednamed.add_argument("--sd",
                               dest="fs_subjects_dir",
                               help=("Freesurfer subjects directory "
                                     "SUBJECTS_DIR"),
                               required=True)
    requirednamed.add_argument("--srcsubject",
                               dest="fs_subject",
                               help="Freesurfer subject name",
                               required=True)
    requirednamed.add_argument("--t1",
                               dest="t1_file",
                               help="T1 image file",
                               metavar="FILE",
                               required=True)
    requirednamed.add_argument("--fa",
                               dest="fa_file",
                               help="FA image file",
                               metavar="FILE",
                               required=True)
    requirednamed.add_argument("--nodif",
                               dest="nodif_file",
                               help="nodif image file",
                               metavar="FILE",
                               required=True)
    requirednamed.add_argument("--seedref",
                               dest="seedref",
                               help="seedref file for probtrackx2",
                               metavar="FILE",
                               required=True)
    requirednamed.add_argument("--samples",
                               dest="samples",
                               help="samples file for probtrackx2",
                               metavar="FILE",
                               required=True)
    requirednamed.add_argument("--outputdir",
                               dest="output_dir",
                               help="output directory",
                               required=True)

    # Optional arguments
    parser.add_argument("--trf",
                        dest="trf_file",
                        default=None,
                        help=("TRF file (diffusion space -> structural space "
                              "affine matrix)"),
                        metavar="FILE")
    parser.add_argument("--srcreg",
                        dest="register_dat_file",
                        default=None,
                        help=(".dat registration file as computed by "
                              "tkregister2 (structural space -> FS space "
                              "affine matrix)"),
                        metavar="FILE")
    parser.add_argument("--icoorder",
                        dest="level",
                        default=7,
                        type=int,
                        help=("specifies the order of the icosahedral "
                              "tesselation. the number of vertices is given by "
                              "the formula 10*2^n+2. In general, it is best to "
                              "use the largest size available."))
    parser.add_argument("--fsconfig",
                        dest="fsconfig",
                        default="/i2bm/local/freesurfer/SetUpFreeSurfer.sh",
                        help="Freesurfer .sh config file",
                        metavar="FILE")

    return parser


class ImageMagickWrapper(object):
    """ Parent class for the wrapping of ImageMagick commands.
    """
    def __init__(self, cmd):
        """ Initialize the ImageMagickWrapper class.

        Parameters
        ----------
        cmd: list of str (mandatory)
            the ImageMagick command to execute.
        """
        self.cmd = cmd
        self.available_cmd = ["animate", "compare", "composite", "conjure",
                              "convert", "display", "identify", "import",
                              "mogrify", "montage", "stream"]
        self.check_cmd()
        self.check_installation()

    def __call__(self):
        """ Run the ImageMagick command.

        Returns
        -------
        output: str
            The ImageMagick process function.
        """
        try:
            output = subprocess.check_output(self.cmd)
        except subprocess.CalledProcessError:
            print('***')
            print("Command {0} failed with parameters : {1}\n"
                  .format(self.cmd[0], " ".join(self.cmd[1:])))
            print('***')
            raise
        else:
            return output

    def check_installation(self):
        """ Check if ImageMagick is installed.
        """
        try:
            subprocess.check_output(["convert", "--version"])
        except OSError:
            print("ImageMagick was not found, please install it first.")
            raise

    def check_cmd(self):
        """ Check if it's an ImageMagick command.
        """
        program = self.cmd[0]
        if program not in self.available_cmd:
            raise ValueError("{0} is not a known ImageMagick command."
                             .format(program))


def images_to_gif(input_img_list, output_file, delay=100):
    """ Convert input images to a unique gif file using ImageMagick.

        Parameters
        ----------
        input_img_list: list of str (mandatory)
            List of the input images to combine.
        output_file: str (mandatory)
            The output gif file path.
        delay: int (optional, default 100)
            This option is useful for regulating the animation of image
            sequences ticks/ticks-per-second seconds must expire before the
            display of the next image. The default is no delay between each
            showing of the image sequence. The default ticks-per-second is 100.

    """
    gif_extension = ".gif"
    if not output_file.endswith(gif_extension):
        output_file += gif_extension
    cmd = ["convert"]
    cmd += input_img_list
    cmd += ["+repage", "-set", "delay", "{0}".format(delay), output_file]
    # Execute the ImageMagick command
    magick = ImageMagickWrapper(cmd)
    magick()


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


def mri_vol2surf(hemi, input_volume, output_volume, register_dat_file,
                 fs_subjects_dir, fs_subject, level, fsconfig):
    """ Freesurfer mri_vol2surf wrapper.

        mri_vol2surf - assigns values from a volume to each surface vertex

        Parameters
        ----------
        hemi: str (mandatory)
            hemisphere (lh or rh).
        input_volume: str (mandatory)
            input volume path.
        output_volume: str (mandatory)
            output path.
        register_dat_file: str (mandatory)
            source registration file.
        fs_subjects_dir: str (mandatory)
            Freesurfer subjects directory SUBJECTS_DIR.
        fs_subject: str (mandatory)
            The FS subject name.
        level: int (mandatory)
            order of icosahedron when FS option trgsubject=ico.
            specifies the order of the icosahedral tesselation. the number of
            vertices is given by the formula 10*2^n+2. In general, it is best
            to use the largest size available.
        fsconfig: str (mandatory)
            The Freesurfer .sh config file.

        """
    cmd = ["mri_vol2surf", "--src", input_volume, "--out", output_volume,
           "--srcreg", register_dat_file, "--hemi", hemi, "--trgsubject",
           "ico", "--icoorder", "{0}".format(level), "--surf", "white",
           "--sd", fs_subjects_dir, "--srcsubject", fs_subject, "--noreshape",
           '--out_type', 'mgz']
    # Execute the FS command
    recon = FSWrapper(cmd, shfile=fsconfig)
    recon()
    if recon.exitcode != 0:
        raise FreeSurferRuntimeError(
            recon.cmd[0], " ".join(recon.cmd[1:]), recon.stderr)


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


def get_profile(level, seedref, seed_file, samples, output_dir, t1_file,
                trf_file, register_dat_file, nodif_file, fs_subjects_dir,
                fs_subject, fsconfig):

    """ Computes the tractography using FSL probtrackx2 and projects the result
        on the cortical surface using FS mri_vol2surf.

        Parameters
        ----------
        level: int (mandatory)
            order of icosahedron when FS option trgsubject=ico.
            specifies the order of the icosahedral tesselation. the number of
            vertices is given by the formula 10*2^n+2. In general, it is best
            to use the largest size available.
        seedref: str (mandatory)
            seedref file for probtrackx2.
        seed_file: str (mandatory)
            text file containing seed coordinates.
        samples: str (mandatory)
            samples file for probtrackx2.
        output_dir: str (mandatory)
            output directory.
        t1_file : str (mandatory)
            T1 image file
        trf_file : str (mandatory)
            TRF file
        register_dat_file: str (mandatory)
            .dat registration file as computed by tkregister2
            (structural space -> FS space affine matrix)
        nodif_file: str (mandatory)
            The nodif image file.
        fs_subjects_dir: str (mandatory)
            Freesurfer subjects directory SUBJECTS_DIR.
        fs_subject: str (mandatory)
            Freesurfer subject name.
        fsconfig: str (mandatory)
            The Freesurfer .sh config file."""

    # Generates the diffusion probability map
    proba_files, _ = probtrackx2(simple=True,
                                 seedref=seedref,
                                 out="fdt_paths",
                                 seed=seed_file,
                                 loopcheck=True,
                                 onewaycondition=True,
                                 samples=samples,
                                 mask=seedref,
                                 dir=output_dir)
    assert len(proba_files) == 1

    fdt_paths_file = proba_files[0]
    proba_file_name = os.path.basename(fdt_paths_file)[:-len(".nii.gz")]
    flirt_out = os.path.join(output_dir, proba_file_name + '_flirt.nii.gz')

    # probability map (diffusion space) --> probability map (T1 space)
    # using FSL flirt function and trf file.
    flirt(proba_files[0], t1_file, out=flirt_out, applyxfm=True, init=trf_file)

    # Project the volumic probability map (T1 space) generated with FSL flirt
    # on the cortical surface (Freesurfer space) (both hemispheres) using
    # Freesurfer's mri_vol2surf with the previously generated registration file.

    for hemi in ["lh", "rh"]:
        mri_vol2surf_out = os.path.join(output_dir, '{0}.{1}_vol2surf.mgz'
                                        .format(hemi, proba_file_name))
        mri_vol2surf(hemi, flirt_out, mri_vol2surf_out, register_dat_file,
                     fs_subjects_dir, fs_subject, level, fsconfig)
        qc_profile(hemi, mri_vol2surf_out, nodif_file, fdt_paths_file, 7,
                   fs_subjects_dir, fs_subject, output_dir, fsconfig)


def get_connectogramm(fs_subjects_dir, fs_subject, t1_file, fa_file, nodif_file,
                      seedref, samples, output_dir, level, trf_file,
                      register_dat_file, fsconfig):
    """ Tractography entry point.
        If no trf file is provided, register the fa image on the t1 image
        to get the it.
        If no .dat registration file is provided, generates it.
        Computes the cortical surface segmentation.

        Parameters
        ----------
        fs_subjects_dir: str (mandatory)
            Freesurfer subjects directory SUBJECTS_DIR.
        fs_subject: str (mandatory)
            Freesurfer subject name.
        t1_file : str (mandatory)
            T1 image file
        fa_file : str (mandatory)
            FA image file
        nodif_file: str (mandatory)
            The nodif image file.
        seedref: str (mandatory)
            seedref file for probtrackx2.
        samples: str (mandatory)
            samples file for probtrackx2.
        output_dir: str (mandatory)
            output directory.
        level: int (mandatory)
            order of icosahedron when FS option trgsubject=ico.
            specifies the order of the icosahedral tesselation. the number of
            vertices is given by the formula 10*2^n+2. In general, it is best
            to use the largest size available.
        trf_file : str (mandatory)
            TRF file
        register_dat_file: str (mandatory)
            .dat registration file as computed by tkregister2
            (structural space -> FS space affine matrix)
        fsconfig: str (mandatory)
            The Freesurfer .sh config file."""

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    subject_output_dir = os.path.join(output_dir,
                                      "{0}_connectivity_results"
                                      .format(fs_subject))
    if not os.path.isdir(subject_output_dir):
        os.mkdir(subject_output_dir)

    if trf_file is None:
        # Register fa image to t1 image to get the trf file
        trf_file = os.path.join(subject_output_dir, "cortex_dest_to_t1.trf")
        reg_file = os.path.join(subject_output_dir, "cortex_dest_to_t1.nii.gz")
        flirt(fa_file, t1_file, omat=trf_file, out=reg_file, usesqform=False,
              cost="normmi", dof=6)

    subject_dir = os.path.join(fs_subjects_dir, fs_subject)

    if register_dat_file is None:
        # Generate the dat registration file (T1 space --> Freesurfer space).
        register_dat_file = os.path.join(subject_dir, "convert",
                                         "register.native.dat")
        conformed_to_native_space(fs_subjects_dir, subject_output_dir)

    # Get physical space to index space matrix
    physical_to_index = numpy.linalg.inv(nibabel.load(t1_file).get_affine())

    # Register FA on T1 to get the affine matrix
    voxel_dest_to_t1 = flirt2aff(trf_file, fa_file, t1_file)
    voxel_t1_to_dest = numpy.linalg.inv(voxel_dest_to_t1)

    # Read the cortex surface segmentation of Freesurfer using previously
    # generated affine matrices.
    seg = read_cortex_surface_segmentation(subject_dir, physical_to_index,
                                           fsconfig, voxel_t1_to_dest)

    # Launch the tractography on each point of the cortical surface
    # Go through the 2 hemispheres
    for hemi in ['lh', 'rh']:
        # Go through all the hemisphere vertices :
        for point in seg[hemi].vertices:
            # Create the output directory for this seed
            point_outdir = os.path.join(subject_output_dir,
                                        "{0}_{1}_{2}_{3}".format(hemi, *point))
            os.mkdir(point_outdir)
            # Write seed coordinates to file
            seed_file = os.path.join(point_outdir, "fdt_coordinates.txt")
            with open(seed_file, 'w') as f:
                for coordinate in point:
                    print(coordinate, file=f)
            # Launch the tractography on this seed
            get_profile(level, seedref, seed_file, samples, point_outdir,
                        t1_file, trf_file, register_dat_file, nodif_file,
                        fs_subjects_dir, fs_subject, fsconfig)


if __name__ == "__main__":

    # Parse terminal arguments
    result = get_parser().parse_args()
    arguments = dict(result._get_kwargs())

    # Required arguments
    fs_subjects_dir = arguments["fs_subjects_dir"]
    fs_subject = arguments["fs_subject"]
    t1_file = arguments["t1_file"]
    fa_file = arguments["fa_file"]
    nodif_file = arguments["nodif_file"]
    seedref = arguments["seedref"]
    samples = arguments["samples"]
    output_dir = arguments["output_dir"]

    # Optional arguments
    trf_file = arguments["trf_file"]
    register_dat_file = arguments["register_dat_file"]
    level = arguments["level"]
    fsconfig = arguments["fsconfig"]

    # Launch the main program
    get_connectogramm(fs_subjects_dir, fs_subject, t1_file, fa_file, nodif_file,
                      seedref, samples, output_dir, level, trf_file,
                      register_dat_file, fsconfig)
