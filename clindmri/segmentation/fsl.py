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

# Clindmri import
from clindmri.extensions.fsl import FSLWrapper
from clindmri.extensions.fsl.exceptions import FSLRuntimeError


def fslreorient2std(input, output, shfile="/etc/fsl/5.0/fsl.sh"):
    """ FSL tool for reorienting the image to match the
    approximate orientation of the standard template images (MNI152).
    It only applies 0, 90, 180 or 270 degree rotations.
    It is not a registration tool.
    It requires NIfTI images with valid orientation information
    in them (seen by valid labels in FSLView).  This tool
    assumes the labels are correct - if not, fix that before using this.
    If the output name is not specified the equivalent transformation
     matrix is written to the standard output

    Usage: fslreorient2std <input_image> [output_image]  
    """
    # Call fslreorient2std
    fslprocess = FSLWrapper("fslreorient2std", shfile=shfile)
    fslprocess()
    if fslprocess.exitcode != 0:
        raise FSLRuntimeError(fslprocess.cmd[0], " ".join(fslprocess.cmd[1:]),
                              fslprocess.stderr)


def bet2(input, output, o=None, m=None, s=None, n=None, f=0.5, g=0, r=None,
         c=None, t=None, e=None, shfile="/etc/fsl/5.0/fsl.sh"):
    """ Wraps command bet2.

    The basic usage is:
        bet2 <input> <output> [options]
    Main bet2 options:
        -o              generate brain surface outline overlaid onto original
                        image
        -m              generate binary brain mask
        -s              generate rough skull image (not as clean as what
                        betsurf generates)
        -n              don't generate the default brain image output
        -f <f>          fractional intensity threshold (0->1); default=0.5;
                        smaller values give larger brain outline estimates
        -g <g>          vertical gradient in fractional intensity threshold
                        (-1->1); default=0; positive values give larger brain
                        outline at bottom, smaller at top
        -r <r>          head radius (mm not voxels); initial surface sphere is
                        set to half of this
        -c < x y z>     centre-of-gravity (voxels not mm) of initial mesh
                        surface.
        -t              apply thresholding to segmented brain image and mask
        -e              generates brain surface as mesh in .vtk format.

    Returns
    -------
    output: str
        the extracted brain volume.
    mask_file: str
        the binary mask of the extracted brain volume.
    mesh_file: str
        the brain surface as a vtk mesh.
    outline_file: str
        the brain surface outline overlaid onto original image.
    inskull_mask_file, inskull_mesh_file,
    outskull_mask_file, outskull_mesh_file, 
    outskin_mask_file, outskin_mesh_file,
    skull_mask_file: str
        rough skull image.
    shfile: str (optional, default NeuroSpin path)
        the path to the FSL 'fsl.sh' configuration file.
    """
    # Call bet2
    fslprocess = FSLWrapper("bet2", shfile=shfile)
    fslprocess()
    if fslprocess.exitcode != 0:
        raise FSLRuntimeError(fslprocess.cmd[0], " ".join(fslprocess.cmd[1:]),
                              fslprocess.stderr)

    # Format outputs
    image_ext = fslprocess.output_ext[fslprocess.environment["FSLOUTPUTTYPE"]]
    mask_file = None
    if m is not None:
        mask_file = output + "_mask" + image_ext
    mesh_file = None
    if e is not None:
        mesh_file = basename + "_mesh.vtk"
    outline_file = None
    if o is not None:
        outline_file = output + "_outline" + image_ext
    inskull_mask_file = None
    inskull_mesh_file = None
    outskull_mask_file = None
    outskull_mesh_file = None
    outskin_mask_file = None
    outskin_mesh_file = None
    skull_mask_file = None    
    if s is not None:
        inskull_mask_file = output + "_inskull_mask" + image_ext
        inskull_mesh_file = output + "_inskull_mesh" + image_ext
        outskull_mask_file = output + "_outskull_mask" + image_ext
        outskull_mesh_file = output + "_outskull_mesh" + image_ext
        outskin_mask_file = output + "_outskin_mask" + image_ext
        outskin_mesh_file = output + "_outskin_mesh" + image_ext
        skull_mask_file = output + "_skull_mask" + image_ext 
    if n is not None:
        output = None 

    return (output, mask_file, mesh_file, outline_file, inskull_mask_file,
            inskull_mesh_file, outskull_mask_file, outskull_mesh_file, 
            outskin_mask_file, outskin_mesh_file, skull_mask_file)
