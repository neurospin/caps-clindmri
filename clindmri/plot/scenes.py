##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import numpy
import logging

# Define the logger
logger = logging.getLogger(__name__)

# Caps import
import pvtk


def bundle_representative_track_scene(tracks, representative_track_indx):
    """ Scene that shows the bundle and its most representative element.

    Parameters
    ----------
    tracks : sequence (N, )
       of tracks as arrays, shape (N1,3) .. (Nm,3).
    representative_track_indx: int
       index of the representative track of the bundle.

    Returns
    ----------
    actors: list of vtkActor
        the scene actors.
    """
    bundle_actor = pvtk.line(tracks, 1)
    representative_track_actor = pvtk.line(
        tracks[representative_track_indx], 0, linewidth=2)
    return [bundle_actor, representative_track_actor]


def field_directions(field):
    """ Scene the shows the directions of a vector field.

    Parameters
    ----------
    field: array (X, Y, N, 3)
        the vector field to plot where N is the number of peaks.

    Returns
    ----------
    actors: list of vtkActor
        the scene actors.
    """
    actors = []
    for x in range(field.shape[0]):
        for y in range(field.shape[1]):
            line = numpy.zeros((2, 3), dtype=numpy.single)
            for vector in field[x, y]:
                line[1] = vector
                actors.append(pvtk.line(line, 0, linewidth=2))
                actors[-1].SetPosition((x, y, 0))
    return actors
