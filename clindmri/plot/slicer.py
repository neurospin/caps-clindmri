#! /usr/bin/env python
##########################################################################
# Nsap - Neurospin - Berkeley - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import numpy

# Nilearn import
from nilearn import plotting

# Matplotlib import
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm

# Clindmri import
from .animate import split_image
from .animate import get_image_dimensions
from .animate import images_to_gif


def plot_image(input_file, edge_file=None, overlay_file=None,
               contour_file=None, snap_file=None, name=None,
               overlay_cmap=None, cut_coords=None):
    """ Plot image with edge/overlay/contour on top (useful for checking
    registration).

    Parameters
    ----------
    input_file: str (mandatory)
        An image to display.
    edge_file: str (optional, default None)
        The target image to extract the edges from.
    overlay_file: str (optional, default None)
        Image to superimpose to the input_file image.
    contour_file: str (optional, default None)
        Superimpose the contour of the image to the input_file image.
    snap_file: str (optional, default None)
        The destination file: if not specified will be the 'input_file'
        name with a 'png' extension.
    name: str (optional, default None)
        The name of the plot.
    overlay_cmap: str (optional, default None)
        The color map to use: 'cold_hot' or 'blue_red' or None,
        or a N*4 array with values in 0-1 (r, g, b, a).
    cut_coords: 3-uplet (optional, default None)
        The MNI coordinates of the point where the cut is performed.
        If None is given, the cuts is calculated automaticaly.

    Returns
    -------
    snap_file: str
        A png snap of the image.
    """
    # Check the input images exist on the file system
    for in_file in [input_file, edge_file, overlay_file, contour_file]:
        if in_file is not None and not os.path.isfile(in_file):
            raise ValueError("'{0}' is not a valid filename.".format(in_file))

    # Check that the snap_file has been specified
    if snap_file is None:
        snap_file = input_file.split(".")[0] + ".png"

    # Create the plot
    kwargs = {}
    if isinstance(cut_coords, int):
        kwargs["display_mode"] = "z"
    display = plotting.plot_anat(input_file, title=name or "",
                                 cut_coords=cut_coords, **kwargs)
    if edge_file is not None:
        display.add_edges(edge_file)
    if overlay_file is not None:

        # Create a custom discrete colormap
        if isinstance(overlay_cmap, numpy.ndarray):
            cmap = colors.LinearSegmentedColormap.from_list(
                "my_colormap", overlay_cmap)
            #cmap, _ = colors.from_levels_and_colors(
            #    range(overlay_cmap.shape[0] + 1), overlay_cmap)
        # This is probably a matplotlib colormap
        elif overlay_cmap in dir(plotting.cm):
            cmap = getattr(plotting.cm, overlay_cmap)
        elif overlay_cmap in dir(cm):
            cmap = cm.get_cmap(overlay_cmap)
        # Default
        else:
            cmap = plotting.cm.alpha_cmap((1, 1, 0))

        display.add_overlay(overlay_file, cmap=cmap)
    if contour_file is not None:
        display.add_contours(contour_file, alpha=0.6, filled=True,
                             linestyles="solid")
    display.savefig(snap_file)
    display.close()

    return snap_file


def animate_image(input_file, cut_coords, edge_file=None, overlay_file=None,
                  contour_file=None, snap_file=None, outdir=None, name=None,
                  overlay_cmap=None, clean=True):
    """ Animate an image with edge/overlay/contour on top (useful for checking
    registration).

    Parameters
    ----------
    input_file: str (mandatory)
        An image to display.
    edge_file: str (optional, default None)
        The target image to extract the edges from.
    overlay_file: str (optional, default None)
        Superimpose the image to the input_file image.
    contour_file: str (optional, default None)
        Superimpose the contour of the image to the input_file image.
    snap_file: str (optional, default None)
        The destination file: if not specified will be the 'input_file'
        name with a 'gif' extension.
    outdir: str (optional, default None)
        The destination directory: if not specified will be the 'input_file'
        directory.
    name: str (optional, default None)
        The name of the plot.
    overlay_cmap: str (optional, default None)
        The color map to use: 'cold_hot' or 'blue_red' or None,
        or a N*4 array with values in 0-1 (r, g, b, a).
    cut_coords: int (mandatory)
        Number of slices to extract in z-direction (nb slices in the GIF).
    clean: bool (default True)
        It True delete the temporary snaps.

    Returns
    -------
    gif_image: str
        A gif snap of the image.
    """
    # Check that the snap_file has been specified
    if outdir is None:
        outdir = os.path.dirname(input_file)
        
    # Check that the snap_file has been specified
    if snap_file is None:
        snap_file = input_file.split(".")[0] + ".png"
    
    # Make sure there is the right extension
    if not snap_file.endswith(".png"):
        snap_file += ".png"

    # Create an image stacking the slices
    plot_image(input_file, edge_file=edge_file, overlay_file=overlay_file,
               contour_file=contour_file, snap_file=snap_file, name=name,
               overlay_cmap=overlay_cmap, cut_coords=cut_coords)

    # Get the large image dimensions
    width, height = get_image_dimensions(snap_file)

    # Split the generated large image into multiple layers to build the gif
    # file
    snap_files = split_image(
        snap_file, outdir, "tmp_", width / cut_coords, height)

    # Combine the png layers to obtain the gif file
    # The GIF has the same path as tne PNG but with .gif extension
    gif_image = os.path.join(outdir, snap_file[:-len(".png")] + ".gif")
    images_to_gif(snap_files, gif_image, delay=10)

    # Clean
    if clean:
        os.remove(snap_file)
        for fname in snap_files:
            os.remove(fname)

    return gif_image


def plot_matrix(input_file, snap_file=None, name=None, transform=None):
    """ Display a matrix.

    Parameters
    ----------
    input_file: str (mandatory)
        A matrix to display.
    snap_file: str (optional, default None)
        The destination file: if not specified will be the 'input_file'
        name with a 'pdf' extension.
    name: str (optional, default None)
        The name of the plot.
    transform: callable (optional, default None)
        A callable function applied on the matrix.

    Returns
    -------
    snap_file: str
        A pdf snap of the matrix.
    """
    # Check the input image exists
    if not os.path.isfile(input_file):
        raise ValueError("'{0}' is not a valid filename.".format(input_file))

    # Check that the snap_file has been specified
    if snap_file is None:
        snap_file = input_file.split(".")[0] + ".pdf"
    if not snap_file.endswith(".pdf"):
        snap_file += ".pdf"

    # Create the plot
    matrix = numpy.loadtxt(input_file)
    if transform is not None:
        matrix = transform(matrix)
    pdf = PdfPages(snap_file)
    try:
        fig = plt.figure()
        plt.title(name or "")
        ax = fig.add_subplot(111)
        ax.imshow(matrix, interpolation="nearest", cmap=plt.cm.jet)
        ax.set_aspect("auto")
        pdf.savefig(fig)
        plt.close()
        pdf.close()
    except:
        pdf.close()
        raise

    return snap_file
