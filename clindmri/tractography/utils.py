#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import numpy
import json
import copy

# Clindmri import
from clindmri.tractography import loadtxt


class Tractogram(object):
    """ A class representing white matter fibers, ie. tractogram.

    Parameters
    ----------
    track_file: str (mandatory)
        a text file containing tracks.
    label_file: str (optional, default None)
        a file containing labels that represent bundles.
    """
    def __init__(self, track_file, label_file=None):
        """ Create a new tractogram.
        """
        self.tracks = loadtxt(track_file)
        if label_file is None:
            self.labels = [0, ] * self.shape()
        else:
            self.labels = {}
            with open(label_file) as open_file:
                clusters = json.load(open_file)
                for clusterid, clusteritem in clusters.items():
                    for trkindex in clusteritem["indices"]:
                        self.labels[trkindex] = clusterid    
        self.length = [None, ] * self.shape()

    def shape(self):
        """ Tractogram shape.

        Returns
        -------
        out: int
            the number of track in the tractogram.
        """
        return len(self.tracks)

    def track_length(self, index=None):
        """ Tracks lengths.

        Parameters
        ----------
        index: int (optional, default None)
            compute the length of the index th track if specified,
            all of them otherwise.

        Returns
        -------
        length: list of float
            the requested track length.
        """
        if index is not None and index < self.shape():
            if self.length[index] is None:
                self.length[index] = length(self.tracks[index])
            return [self.length[index]]

        if index is None:
            for index in range(self.shape()):
                if self.length[index] is None:
                    self.length[index] = length(self.tracks[index])
            return self.length

        return []

    def downsample(self, index=None, nb_points=3):
        """ Downsample tracks.

        Parameters
        ----------
        index: int (optional, default None)
            downsample thet index th track if specified,
            all of them otherwise.
        nb_points: int (optional, default 3)
           the number of sampling points along the curve.

        Returns
        -------
        dowsample_tracks: list of array
            the downsampled tracks.
        """
        if index is not None and index < self.shape():
            return dowsample(self.tracks[index], nb_points)

        if index is None:
            dowsample_tracks = []
            for index in range(self.shape()):
                dowsample_tracks.append(
                    dowsample(self.tracks[index], nb_points))
            return dowsample_tracks

        return []

    def density(self, shape):
        """ Compute a track density map.

        Parameters
        ----------
        shape: 3-uplet
            the shape of the output density map.

        Returns
        -------
        density: array
            the track density map.
        """
        density = numpy.zeros(shape, dtype=numpy.single)
        for track in self.tracks:
            indices = copy.deepcopy(track).astype(int)
            density[indices.T.tolist()] += 1
        return density

    def endpoints(self):
        """ Return the tracks end points.
    
        Returns
        -------
        endpoints: array (n, 2, 3)
            the track end points.
        """
        return numpy.asarray([line[0::len(line)-1] for line in self.tracks])

    def apply_affine(self, affine):
        """ Apply an affine transformation on the track points.

        Parameters
        ----------
        affine: array (4, 4)
            the affine transformation to applied in voxel coordinates.

        Returns
        -------
        affine_tracks: list of array
            the tractogram tracks in the new coordiantes.
        """
        affine_tracks = []
        for track in self.tracks:
            ones = numpy.ones((track.shape[0], 1), dtype=track.dtype)
            track_homogeneous = numpy.concatenate((track, ones), axis=1)
            new_track = numpy.dot(affine, track_homogeneous.T)
            affine_tracks.append(new_track.T[:, :3])
        self.tracks = affine_tracks
        return affine_tracks

    def apply_affine_on_endpoints(self, affine):
        """ Apply an affine transformation on the track end points.

        Parameters
        ----------
        affine: array (4, 4)
            the affine transformation to applied.

        Returns
        -------
        affine_endpoints: array (n, 2, 3)
            the tractogram tracks end points in the new coordiantes.
        """
        affine_endpoints = []
        for track in self.endpoints():
            ones = numpy.ones((track.shape[0], 1), dtype=track.dtype)
            track_homogeneous = numpy.concatenate((track, ones), axis=1)
            new_track = numpy.dot(affine, track_homogeneous.T)
            affine_endpoints.append(new_track.T[:, :3])
        return numpy.asarray(affine_endpoints)


def length(xyz, along=False):
    """ Euclidean length of track line.

    Parameters
    ----------
    xyz : array-like shape (N,3)
       array representing x,y,z of N points in a track.
    along : bool, optional
       If True, return array giving cumulative length along track,
       otherwise (default) return scalar giving total length.

    Returns
    -------
    L : scalar or array shape (N-1,)
       scalar in case of 'along' == False, giving total length, array if
       'along' == True, giving cumulative lengths.
    """
    xyz = numpy.asarray(xyz)
    if xyz.shape[0] < 2:
        if along:
            return numpy.array([0])
        return 0
    dists = numpy.sqrt((numpy.diff(xyz, axis=0) ** 2).sum(axis=1))
    if along:
        return numpy.cumsum(dists)
    return numpy.sum(dists)


def downsample(xyz, n_pols=3):
    """ Downsample for a specific number of points along the curve.

    Uses the length of the curve. It works in as similar fashion to
    midpoint and arbitrarypoint.
    
    Parameters
    ----------
    xyz : array-like shape (N,3)
       array representing x,y,z of N points in a track
    n_pol : int
       integer representing number of points (poles) we need along the curve.

    Returns
    -------
    xyz2 : array shape (M,3)
       array representing x,z,z of M points that where extrapolated. M
       should be equal to n_pols
    """
    # Ensure we are working with numpy array
    xyz = numpy.asarray(xyz)

    # Special cases when no track are passed or when the track is represented by
    # one point
    n_pts = xyz.shape[0]
    if n_pts == 0:
        raise ValueError("xyz array cannot be empty.")
    if n_pts == 1:
        return xyz.copy().squeeze()

    # Compute the cumulative fiber euclidean length
    cumlen = numpy.zeros(n_pts)
    cumlen[1:] = length(xyz, along=True)

    # Compute the output fiber resampling step
    if n_pols <= 2:
        raise ValueError("The given number of points 'n_pols={0}' needs to be "
                         "higher than 2.".format(n_pols))
    step = cumlen[-1] / (n_pols - 1)
    if cumlen[-1] < step:
        raise ValueError("The given number of points 'n_pols={0}' is "
                         "incorrect.".format(n_pols))

    # The resampling procedure
    xyz2 = [_extrap(xyz, cumlen, distance) 
            for distance in numpy.arange(0, cumlen[-1], step)]

    return numpy.vstack((numpy.array(xyz2), xyz[-1]))


def _extrap(xyz, cumlen, distance):
    """ Helper function for extrapolation.
    """
    # Find where is the new point: interval  
    ind = numpy.where((cumlen - distance) > 0)[0][0]
    len0 = cumlen[ind - 1]        
    len1 = cumlen[ind]

    # Linear interpolation
    d0 = (distance - len0) / (len1 - len0)
    return d0 * xyz[ind] + (1 - d0) * xyz[ind - 1]
