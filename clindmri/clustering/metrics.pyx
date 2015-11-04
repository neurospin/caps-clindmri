#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
# 
# From dipy: http://nipy.org/dipy/
##########################################################################

# System import
cimport cython
import numpy as np
cimport numpy as cnp
cimport openmp
from cython cimport parallel
from libc.stdlib cimport abort, malloc, free, realloc

# Float 32 dtype for casting
cdef cnp.dtype f32_dt = np.dtype(np.float32)
cdef cnp.float32_t inf = np.inf

# Some cdefs
cdef extern from "stdlib.h" nogil:
    ctypedef unsigned long size_t
cdef extern from "math.h" nogil:
    float sqrt(float x)

# Metric type maping
metric_map = {
    "average": 0,
    "min": 1,
    "max": 2
}

# Define a track structure
ctypedef struct Track:
    float *track
    unsigned long nb_of_points


###############################################################################
# Python callable functions
###############################################################################

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def clustering_metric(tracks,
                      object[float, ndim=1] metric=None,
                      distance_metric="min",
                      float zhang_thr=0.5,
                      int nb_of_threads=1):
    """ Efficient parallel clustering metric computation.

    Use the method described in S. Zhang: 'Identifying White-Matter Fiber
    Bundles in DTI Data Using an Automated Proximity-Based Fiber Clustering
    Method.'
    This method is based on average minimume distance between two tracks.

    A condensed distance matrix is returned.
    A condensed distance matrix is a flat array containing the upper
    triangular of the symetric distance matrix.

    Parameters
    ----------
    tracks : sequence
        of tracks as arrays, shape (N1,3) .. (Nm,3).
    distance_metric: str
       metric to calculate: {'avg', 'min', 'max'}.
    zhang_thr: float
        a threshold distance below which point to point distance is not
        considered.
    nb_of_threads: int (default 1)
        the number of threads to use.

    Returns
    -------
    metric : array
        a (M * (M - 1) / 2) array containing the condensed track to track
        distances.
    """
    cdef:
        unsigned long number_of_tracks=0, condensed_matrix_size=0
        int metric_type
        unsigned long i, j

    # Some function parameters
    metric_type = metric_map[distance_metric]
    number_of_tracks = len(tracks)
    condensed_matrix_size = int(
        number_of_tracks * (number_of_tracks - 1.) / 2.)

    # Allocate the output if necessary
    if metric is None:
        metric = np.zeros((condensed_matrix_size, ), dtype=np.single)
    if metric.shape[0] != condensed_matrix_size:
        raise TypeError("Metric array must have shape: {0}".format(
            condensed_matrix_size))

    # Format input tracks to be compliant with c loops: make a copy
    # Find the longest track number of samples: used to build memory layout.
    cdef:
        cnp.ndarray[float, ndim=2] track
        Track *ctracks
    ctracks = <Track *>malloc(number_of_tracks * sizeof(Track))
    for i from 0 <= i < number_of_tracks:

        # Track copy
        ctracks[i].nb_of_points = tracks[i].shape[0]
        ctracks[i].track = <float *>malloc(ctracks[i].nb_of_points * 3 *
                                           sizeof(float))
        track = np.ascontiguousarray(tracks[i], dtype=np.float32)
        ptr = <float *>track.data
        for j from 0 <= j < (ctracks[i].nb_of_points * 3):
            ctracks[i].track[j] = ptr[j]

        # Max search
        #if ctracks[i].nb_of_points > max_number_of_points:
        #    max_number_of_points = ctracks[i].nb_of_points

    # Parallel distance computation: loop over tracks
    cdef:
        unsigned long x, y, inner, outer, metric_indx
        float dist=0
    outer = number_of_tracks - 1
    inner = number_of_tracks
    with nogil, cython.boundscheck(False), cython.wraparound(False):

        # Outer loop
        for x in xrange(outer):

            # Inner parallel loop
            for y in parallel.prange(x + 1, inner, schedule="static",
                                     num_threads=nb_of_threads):

                # Preallocate the min buffer array for track distance
                # calculations
                min_buffer = <float *>malloc(
                    (ctracks[x].nb_of_points + ctracks[y].nb_of_points) *
                    sizeof(float))

                # Average minimum distance between two tracks                   
                dist = czhang(ctracks[x].nb_of_points, ctracks[x].track,
                              ctracks[y].nb_of_points, ctracks[y].track,
                              min_buffer, metric_type, zhang_thr)

                # Store the distance between the two tracks in the output
                # metric array
                metric_indx = (
                    x * (number_of_tracks - 1) - x * (x - 1) / 2 + y - x - 1)
                metric[metric_indx] = dist

                # Free memory
                free(min_buffer)

    # Free memory
    with nogil:
        for i from 0 <= i < number_of_tracks:
            free(ctracks[i].track)
        free(ctracks)

    return metric

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def clustering_metric_resample(object[float, ndim=3] tracks not None,
                               object[float, ndim=1] metric=None,
                               int nb_of_threads=1):
    """ Efficient parallel clustering metric computation.

    All the tracks are assumed to have the same sampling scheme in order to
    speed up the point association during the track distance computation.

    A condensed distance matrix is returned.
    A condensed distance matrix is a flat array containing the upper
    triangular of the symetric distance matrix.

    Parameters
    ----------
    tracks : array (M, N, 3)
        M track, each track is expected to be of the form
        (N, 3) where N=number_of_samples.
    nb_of_threads: int (default 1)
        the number of threads to use.

    Returns
    -------
    metric : array
        a (M * (M + 1) / 2) - M array containing the condensed track to track
        distances.
    """
    cdef:
        unsigned long number_of_tracks=0, number_of_samples=0
        unsigned long condensed_matrix_size=0

    # Some function parameters
    number_of_tracks = tracks.shape[0]
    number_of_samples = tracks.shape[1]
    condensed_matrix_size = int(
        number_of_tracks * (number_of_tracks + 1.) / 2. - number_of_tracks)

    # Allocate the output if necessary
    if metric is None:
        metric = np.zeros((condensed_matrix_size, ), dtype=np.single)
    if metric.shape[0] != condensed_matrix_size:
        raise TypeError("Metric array must have shape: {0}".format(condensed_matrix_size))

    # Parallel distance computation: loop over tracks
    cdef:
        float tmp_dist=0, tmp_flip_dist=0
        float diff=0, flip_diff=0
        float dist=0, flip_dist=0
        unsigned long x, y, inner, outer, i, j, flip_i, metric_indx
    outer = number_of_tracks - 1
    inner = number_of_tracks
    with nogil, cython.boundscheck(False), cython.wraparound(False):

        # Outer loop
        for x in xrange(outer):

            # Inner parallel loop
            for y in parallel.prange(x + 1, inner, schedule="static",
                                     num_threads=nb_of_threads):

                # Direct and flip average distance between two tracks
                dist = 0
                flip_dist = 0
                for i from 0 <= i < number_of_samples:
                    tmp_dist = 0
                    tmp_flip_dist = 0
                    for j from 0 <= j < 3:
                        flip_i = number_of_samples - 1 - i
                        diff = tracks[x, i, j] - tracks[y, i, j]
                        flip_diff = (tracks[x, i, j] - 
                                     tracks[y, flip_i, j])
                        diff = diff * diff
                        flip_diff = flip_diff * flip_diff 
                        tmp_dist = tmp_dist + diff
                        tmp_flip_dist = tmp_flip_dist + flip_diff
                    dist += sqrt(tmp_dist)
                    flip_dist += sqrt(tmp_flip_dist)

                dist = dist / <float>number_of_samples
                flip_dist = flip_dist / <float>number_of_samples

                # 'dist' should contains the smallest distance
                if flip_dist < dist:
                    dist = flip_dist

                # Store the distance between the two tracks in the output
                # metric array
                metric_indx = (
                    x * (number_of_tracks - 1) - x * (x - 1) / 2 + y - x - 1)
                metric[metric_indx] = dist

    return metric


def mam_distances(xyz1, xyz2, zhang_thr=0., metric="all"):
    """ Min/Max/Mean Average Minimume Distance between tracks xyz1 and xyz2.
    
    Based on the metrics in Zhang, et al 2008: 'Identifying White-Matter Fiber
    Bundles in DTI Data Using an Automated Proximity-Based Fiber Clustering
    Method.' which in turn are based on those of Corouge et al. 2004
    
    Parameters
    ----------
    xyz1: array, shape (N1,3), dtype float32
       array representing the N2 samples (x, y, z) of the track.
    xyz2: array, shape (N2,3), dtype float32
       array representing the N1 samples (x, y, z) of the track.
	zhang_thr: float
       minimum threshold so that distances below it are not considered.
    metric: str
       metric to calculate: {'avg', 'min', 'max'} return a scalar; 'all'
       returns a tuple.
       
    Returns
    -------
    avg_mcd: float
       average mean closest distance.
    min_mcd: float
       minimum mean closest distance.
    max_mcd: float
       maximum mean closest distance.
                   
    Notes
    -----
    Algorithmic description:
    
    Lets say we have curves A and B.
    
    For every point in A calculate the minimum distance from every point
    in B stored in minAB
    
    For every point in B calculate the minimum distance from every point
    in A stored in minBA
    
    Find average of minAB stored as avg_minAB
    Find average of minBA stored as avg_minBA
    
    If metric is 'avg' then return (avg_minAB + avg_minBA) / 2.0
    If metric is 'min' then return min(avg_minAB, avg_minBA)
    If metric is 'max' then return max(avg_minAB, avg_minBA)
    """
    cdef:
        cnp.ndarray[cnp.float32_t, ndim=2] track1 
        cnp.ndarray[cnp.float32_t, ndim=2] track2
        size_t t1_len, t2_len

    # Format each track (ctype + float) and store each track number of points
    track1 = np.ascontiguousarray(xyz1, dtype=f32_dt)
    t1_len = track1.shape[0]
    track2 = np.ascontiguousarray(xyz2, dtype=f32_dt)
    t2_len = track2.shape[0]

    # Preallocate buffer array for track distance calculations
    cdef:
        cnp.float32_t *min_t2t1, *min_t1t2
        cnp.ndarray [cnp.float32_t, ndim=1] distances_buffer
    distances_buffer = np.zeros((t1_len + t2_len,), dtype=np.float32)
    min_t2t1 = <cnp.float32_t *> distances_buffer.data
    min_t1t2 = min_t2t1 + t2_len

    # Calculate the minimum distance from every points
    min_distances(t1_len, <cnp.float32_t *>track1.data,
                  t2_len, <cnp.float32_t *>track2.data,
                  min_t2t1,
                  min_t1t2)

    # Compute minimum distances averages
    cdef:
        size_t t1_pi, t2_pi, t1_zhang_len=0, t2_zhang_len=0
        cnp.float32_t mean_t2t1 = 0, mean_t1t2 = 0
    for t1_pi from 0 <= t1_pi < t1_len:
        if min_t1t2[t1_pi] >= zhang_thr:
            mean_t1t2 += min_t1t2[t1_pi]
            t1_zhang_len += 1
    mean_t1t2 = mean_t1t2 / t1_zhang_len
    for t2_pi from 0 <= t2_pi < t2_len:
        if min_t2t1[t2_pi] >= zhang_thr:
            mean_t2t1 += min_t2t1[t2_pi]
            t2_zhang_len += 1
    mean_t2t1 = mean_t2t1 / t2_zhang_len

    # Return the appropriate metric
    if metric == "all":
        return (
            (mean_t2t1 + mean_t1t2) / 2.0,
            np.min((mean_t2t1, mean_t1t2)),
            np.max((mean_t2t1, mean_t1t2)))
    elif metric == "avg":
        return (mean_t2t1 + mean_t1t2) / 2.0
    elif metric == "min":            
        return np.min((mean_t2t1, mean_t1t2))
    elif metric == "max":
        return np.max((mean_t2t1, mean_t1t2))
    else:
        raise ValueError("Wrong metric type '{0}'.".format(metric))


###############################################################################
# Intern functions
###############################################################################

@cython.cdivision(True)
cdef inline void min_distances(size_t t1_len,
                               cnp.float32_t *track1_ptr,
                               size_t t2_len,
                               cnp.float32_t *track2_ptr,
                               cnp.float32_t *min_t2t1,
                               cnp.float32_t *min_t1t2) nogil:
    """ Calculate the minimum distance from every points.

    Note ``nogil`` - no python calls allowed in this function.
    """
    cdef:
        cnp.float32_t *t1_pt, *t2_pt, d0, d1, d2
        cnp.float32_t delta2
        int t1_pi, t2_pi

    # Initilaize the distances to inf
    for t2_pi from 0 <= t2_pi < t2_len:
        min_t2t1[t2_pi] = inf
    for t1_pi from 0 <= t1_pi < t1_len:
        min_t1t2[t1_pi] = inf

    # Create a pointer to current point in track 1
    t1_pt = track1_ptr

    # Calculate min squared Euclidean distance.
    for t1_pi from 0<= t1_pi < t1_len:

        # Create a pointer to current point in track 2
        t2_pt = track2_ptr
        for t2_pi from 0 <= t2_pi < t2_len:

            # Squared Euclidean distance between two 3d points
            d0 = t1_pt[0] - t2_pt[0]
            d1 = t1_pt[1] - t2_pt[1]
            d2 = t1_pt[2] - t2_pt[2]
            delta2 = d0 * d0 + d1 * d1 + d2 * d2

            # Update buffers if necessary
            if delta2 < min_t2t1[t2_pi]:
                min_t2t1[t2_pi] = delta2
            if delta2 < min_t1t2[t1_pi]:
                min_t1t2[t1_pi] = delta2

            # Move pointer to next point in track 2
            t2_pt += 3

        # Move pointer to next point in track 1
        t1_pt += 3

    # Sqrt to get Euclidean distance from squared distance
    for t1_pi from 0 <= t1_pi < t1_len:
        min_t1t2[t1_pi] = sqrt(min_t1t2[t1_pi])
    for t2_pi from 0 <= t2_pi < t2_len:
        min_t2t1[t2_pi] = sqrt(min_t2t1[t2_pi])    


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef float czhang(int track1_size,
                  float *track1,
                  int track2_size,
                  float *track2,
                  float *min_distances_pointer,
                  int metric_type,
                  float zhang_thr) nogil:
    """ Calculate the zhang track distance.

    Note ``nogil`` - no python calls allowed in this function.
    """
    # Compute the fiber to fiber distances: set some pointers to store the
    # distances from t1 to t2 and from t2 to t1
    cdef:
        cnp.float32_t *min_t2t1, *min_t1t2

    min_t2t1 = min_distances_pointer
    min_t1t2 = min_distances_pointer + track2_size
    min_distances(track1_size, track1,
                  track2_size, track2,
                  min_t2t1, min_t1t2)
    # Compute minimum distances averages
    cdef:
        unsigned int t1_indx, t2_indx, t1_zhang_len=0, t2_zhang_len=0
        float mean_t2t1=0, mean_t1t2=0
    for t1_indx from 0 <= t1_indx < track1_size:
        if min_t1t2[t1_indx] >= zhang_thr:
            mean_t1t2 += min_t1t2[t1_indx]
            t1_zhang_len += 1
    mean_t1t2 = mean_t1t2 / t1_zhang_len
    for t2_indx from 0 <= t2_indx < track2_size:
        if min_t2t1[t2_indx] >= zhang_thr:
            mean_t2t1 += min_t2t1[t2_indx]
            t2_zhang_len += 1
    mean_t2t1 = mean_t2t1 / t2_zhang_len

    # Return the appropriate metric
    cdef:
        float distance=0
    if metric_type == 0:                
        distance = (mean_t2t1 + mean_t1t2) / 2.0
    elif metric_type == 1:        
        if mean_t2t1 < mean_t1t2:
            distance = mean_t2t1
        else:
            distance = mean_t1t2
    elif metric_type == 2:                
        if mean_t2t1 > mean_t1t2:
            distance = mean_t2t1
        else:
            distance = mean_t1t2                    
    return distance




