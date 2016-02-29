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
import numbers
import types
import nibabel

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
               overlay_cmap=None, cut_coords=None, cutlower=0,
               cutupper=0, display_mode="z", figsize=(10, 10),
               nbcols=6):
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
        name with a 'pdf' extension.
    name: str (optional, default None)
        The name of the plot.
    overlay_cmap: str (optional, default None)
        The color map to use: 'cold_hot' or 'blue_red' or None,
        or a N*4 array with values in 0-1 (r, g, b, a).
    cut_coords: 3-uplet (optional, default None)
        The MNI coordinates of the point where the cut is performed.
        If None is given, the cuts is calculated automaticaly. If an integer
        is given slice the volume in the 'display_mode' direction.
    cutlower, cutupper: int (optional default 0)
        cut values used to reject the begening or the end of the volume.
    display_mode: str (optional, default 'z')
        select the slicing direction in 'x', 'y' or 'z'.
    figsize: 2-uplet (optional, default (10, 10))
        the slice size.
    nbcols: str (optional, default 6)
        the number of slices per line.

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
        snap_file = input_file.split(".")[0] + ".pdf"
    if not snap_file.endswith(".pdf"):
        snap_file += ".pdf"

    # Set options for the two slicing configurations.
    if isinstance(cut_coords, numbers.Number):
        nbrows = cut_coords // nbcols
        axis = "xyz".index(display_mode)
        this_shape = nibabel.load(input_file).get_data().shape[axis]
        samples = numpy.linspace(cutlower, this_shape - cutupper -1, cut_coords,
                                 dtype=int)
        waves = numpy.array_split(samples[:(nbrows * nbcols)], nbrows)
        if cut_coords % nbcols > 0:
            waves.append(samples[(nbrows * nbcols):])
            nbrows += 1
    else:
        nbrows = None
        display_mode = "ortho"

    # Save in Pdf
    pdf = PdfPages(snap_file)
    try:
        if isinstance(cut_coords, numbers.Number):
            for wave in waves:
                cutlower = wave[0]
                cutupper = wave[-1]
                cut_coords = len(wave)
                display = _nilearn_display(
                    input_file, edge_file, overlay_file, contour_file, name,
                    overlay_cmap, cut_coords, display_mode, cutlower, cutupper,
                    figsize)
                pdf.savefig(display.frame_axes.figure)
                display.close()
        else:
            display = _nilearn_display(
                input_file, edge_file, overlay_file, contour_file, name,
                overlay_cmap, cut_coords, display_mode, cutlower, cutupper,
                figsize)
            pdf.savefig(display.frame_axes.figure)
            display.close()
        pdf.close()
    except:
        pdf.close()
        raise   

    return snap_file


class BoundaryNorm(colors.BoundaryNorm):
    """ Litlle fix in matplotlib BoundaryNorm causing
    'IndexError: arrays used as indices must be of integer (or boolean) type'
    error.
    """
    def __init__(self, boundaries, ncolors, clip=False):
        super(BoundaryNorm, self).__init__(boundaries, ncolors, clip)

    def __call__(self, x, clip=None):
        if isinstance(self.vmin, numpy.ma.core.MaskedConstant):
            self.vmin = self.boundaries[0]
        if isinstance(self.vmax, numpy.ma.core.MaskedConstant):
            self.vmax = self.boundaries[-1]
        return super(BoundaryNorm, self).__call__(x, clip)


def _nilearn_display(input_file, edge_file, overlay_file, contour_file,
                     name, overlay_cmap, cut_coords, display_mode,
                     lowerbound, upperbound, figsize):
    """ Create a nice dispaly with nilearn.
    """
    # Monkey patch the 'find_cut_coords' class method in case of int
    # 'cut_coords'
    if isinstance(cut_coords, numbers.Number):
        def monkeypatch(cls):
            def decorator(func):
                setattr(cls, func.__name__, types.MethodType(func, cls))
                return func
            return decorator

        @monkeypatch(plotting.displays.SLICERS[display_mode])
        def find_cut_coords(cls, img=None, threshold=None, cut_coords=None):
            """ Instanciate the slicer and find cut coordinates.
            """
            # Fig size
            cls._default_figsize = list(figsize)

            # Translation
            direction = cls._direction
            n_cuts = cut_coords

            # Misc
            if direction not in "xyz":
                raise ValueError(
                    "'direction' must be one of 'x', 'y', or 'z'. Got '%s'" % (
                        direction))
            affine = img.get_affine()

            # Compute slice location in physical space
            cut_coords = numpy.unique(
                numpy.linspace(lowerbound, upperbound, n_cuts, dtype=int))

            return plotting.find_cuts._transform_cut_coords(
                        cut_coords, direction, affine)
    else:
        plotting.displays.SLICERS[display_mode]._default_figsize = list(figsize)

    # Create the renderer
    display = plotting.plot_anat(
        input_file, title=name or "", cut_coords=cut_coords,
        display_mode=display_mode)
    if edge_file is not None:
        display.add_edges(edge_file)
    if overlay_file is not None:
        # Create a custom discrete colormap
        norm = None
        if isinstance(overlay_cmap, numpy.ndarray):
            cmap = colors.ListedColormap(overlay_cmap.tolist())
            # cmap.set_over((1., 0., 0., 1.))
            # cmap.set_under((1., 0., 0., 1.))
            bounds = numpy.asarray(range(overlay_cmap.shape[0] + 1)) - 0.5
            norm = BoundaryNorm(bounds, cmap.N)
            # m = cm.ScalarMappable(cmap=cmap, norm=norm)
            # print(m.to_rgba(41), m.to_rgba(214))
        elif overlay_cmap in dir(plotting.cm):
            cmap = getattr(plotting.cm, overlay_cmap)
        elif overlay_cmap in dir(cm):
            cmap = cm.get_cmap(overlay_cmap)
        # Default
        else:
            cmap = plotting.cm.alpha_cmap((1, 1, 0))
        display.add_overlay(overlay_file, cmap=cmap, norm=norm)
    if contour_file is not None:
        display.add_contours(contour_file, alpha=0.6, filled=True,
                             linestyles="solid")
    return display


def xyz_mosaics(reffile, overlayfile, nbslices, name, cmap, outdir,
                cutupper=0, cutlower=0):
    """ Compute triplanar mosaic and mosaic in 'xyz' directions.

    Parameters
    ----------
    reffile: str (mandatory)
        the reference volume.
    overlayfile: str (mandatory)
        the volume that will be overlayed (must be in the 'reffile' volume
        coordinates).
    nbslices: int (mandatory)
        the number of slices to slect in the 'xyz' directions.
    name: str (mandatory)
        the mosaic name.
    cmap: array (mandatory)
        a RGBA color map.
    outdir: str (mandatory)
        the result folder where the mosaic are generated.
    cutupper, cutlower: int (optional, default 0)
        the number of slices to skip at the beging/end of the volume in the
        cut direction. If 0, empty slices are automatically rejected.

    Returns
    -------
    mosaics: list of str
        the generated mosaic files.
    """
    # A parameter to store the generated snaps
    mosaics = []

    # Compute overlay mosaics
    array = nibabel.load(overlayfile).get_data()
    if cutupper == 0:
        auto_cutupper = True
    if cutlower == 0:
        auto_cutlower = True
    for display_mode in ["z", "y", "x"]:

        # Compute the cutlower and cutupper bounds: slice with signal
        axis = "xyz".index(display_mode)
        sum_axes = tuple([elem for elem in range(3) if elem != axis])
        weights = array.sum(axis=sum_axes)
        if auto_cutupper:
            cutupper = 0
            for value in weights[::-1]:
                if value == 0:
                    cutupper += 1
                else:
                    break      
        if auto_cutlower:
            cutlower = 0
            for value in weights:
                if value == 0:
                    cutlower += 1
                else:
                    break

        # Display sliced images
        qcname = "{0}-{1}".format(display_mode, name)
        snap_file = os.path.join(outdir, qcname + ".pdf")
        plot_image(reffile, overlay_file=overlayfile, snap_file=snap_file,
                   name=qcname, cut_coords=nbslices, overlay_cmap=cmap,
                   cutlower=cutlower, cutupper=cutupper,
                   display_mode=display_mode)
        mosaics.append(snap_file)

    # Display triplanar
    qcname = "triplanar-{0}".format(name)
    snap_file = os.path.join(outdir, qcname + ".pdf")
    plot_image(reffile, overlay_file=overlayfile, snap_file=snap_file,
               name=qcname, cut_coords=(0, 0, 0), overlay_cmap=cmap,
               cutlower=cutlower, cutupper=cutupper,
               display_mode=display_mode)
    mosaics.append(snap_file)

    return mosaics


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
