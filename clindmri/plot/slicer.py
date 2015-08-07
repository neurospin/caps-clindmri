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

# Matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


def plot_image(input_file, edge_file=None, overlay_file=None, contour_file=None,
               snap_file=None, name=None, overlay_cmap=None):
    """ Plot image with edge/overlay/contour on top (useful for checking
    registration).

    Parameters
    ----------
    input_file: str (mandatory)
        An image to display.
    outline_file: str (optional, default None)
        The target image to extract the edges from.
    snap_file: str (optional, default None)
        The destination file: if not specified will be the 'input_file'
        name with a 'pdf' extension.
    name: str (optional, default None)
        The name of the plot.
    overlay_cmap: str (optional, default None)
        The color map to use: 'cold_hot' or 'blue_red' or None.

    Returns
    -------
    snap_file: str
        A pdf snap of the image.
    """
    # Check the input images exist on the file system
    for in_file in [input_file, edge_file, overlay_file, contour_file]:
        if in_file is not None and not os.path.isfile(in_file):
            raise ValueError("'{0}' is not a valid filename.".format(in_file))

    # Check that the snap_file has been specified
    if snap_file is None:
        snap_file = input_file.split(".")[0] + ".pdf"
    if not snap_file.endswith(".pdf"):
        snape_file += ".pdf"

    # Create the plot
    display = plotting.plot_anat(input_file, title=name or "")
    if edge_file is not None:
        display.add_edges(edge_file)
    if overlay_file is not None:
        cmap = plotting.cm.__dict__.get(
            overlay_cmap, plotting.cm.alpha_cmap((1, 1, 0)))
        display.add_overlay(overlay_file, cmap=cmap)
    if contour_file is not None:
        display.add_contours(contour_file, alpha=0.6, filled=True,
                             linestyles="solid")
    display.savefig(snap_file)
    display.close()

    return snap_file


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
        snape_file += ".pdf"

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

    
