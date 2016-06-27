##########################################################################
# NSAP - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import numpy


def savetxt(fname, tracks, fmt="%.4f", delimiter=","):
    """ Save tracks in text format.

    Parameters
    ----------
    fname: str (mandatory)
        the destination file name.
    tracks: list of array (mandatory)
        a list of tracks with shape (n_i, 3).
    fmt: str (optional, default '%.4f')
        the track elements format in the destination file.
    delimiter: str (optional, default ',')
        character separating columns.
    """
    with open(fname, "w") as open_file:
        for track in tracks:
            for row in track:
                line = delimiter.join(fmt % value for value in row)
                open_file.write(line + "\n")
            open_file.write("\n")


def loadtxt(fname, delimiter=","):
    """ Load tracks saved in text format.

    Parameters
    ----------
    fname: str (mandatory)
        the file name containing the tracks.
    delimiter: str (optional, default ',')
        character separating columns.

    Returns
    -------
    tracks: list of array (mandatory)
        a list of tracks with shape (n_i, 3).
    """
    tracks = []
    array = []
    with open(fname, "r") as open_file:
        for line in open_file:
            if line != "\n":
                array.append(numpy.fromstring(line.rstrip("\n"), dtype=float,
                                              sep=delimiter))
            else:
                tracks.append(numpy.asarray(array))
                array = []
    return tracks
