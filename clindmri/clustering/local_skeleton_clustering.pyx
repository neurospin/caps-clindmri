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
import clindmri.tractography.utils as utils
import sys
from cython cimport parallel
from libc.stdlib cimport abort, malloc, free, realloc

# Float 32 dtype for casting
cdef cnp.dtype f32_dt = np.dtype(np.float32)
cdef cnp.float32_t inf = np.inf
cdef cnp.float32_t biggest_float = np.finfo("f4").max

# Some cdefs
cdef extern from "math.h" nogil:
    float sqrt(float x)

# Define a cluster structure
ctypedef struct Cluster:
    unsigned long *indices
    float *hidden
    long N


###############################################################################
# Python callable functions
###############################################################################

def local_skeleton_clustering(tracks,
                              float cutoff=10.,
                              int nb_samples=10,
                              bar_length=40,
                              int nb_of_threads=1):
    """ Local skeleton clustering.

    Used in the HBM2010 abstract 'Fast Dimensionality Reduction for Brain
    Tractography Clustering' by E.Garyfallidis et.al.

    Parameters
    ----------
    tracks: sequence
        of tracks as arrays, shape (N1,3) .. (Nm,3)
    cutoff: float (default 10)
        average euclidean distance threshold.
    nb_samples: float (default 10)
        the number of points in the downsample tracks used in order to speed 
        up the compututation of the fiber to fiber metric.
    bar_length: int (default 40)
        the length of the progress bar that will be ploted.
    nb_of_threads: int (default 1)
        the number of threads to use.

    Returns
    -------
    organized_cluster: dict of dict
        the first key corresponds to the cluster label, then each cluster is
        represented by a dictionnay with three keys: 'indices' - a reference to
        the fibers included in the bundle, 'N' - the number of fibers in the
        bundle, and 'hidden' - a private parameter that contains the mean
        bundle fiber (update during the clustering).
    """
    cdef :
        cnp.ndarray[cnp.float32_t, ndim=2] track
        #cnp.ndarray[cnp.float32_t, ndim=2] skeleton
        cnp.ndarray[cnp.float32_t, ndim=1] alld
        cnp.ndarray[cnp.int_t, ndim=1] flip
        unsigned long number_of_tracks, number_of_clusters, number_of_points
        unsigned long cluster_it, fiber_it, min_index, i, j
        float *ptr, *skeleton
        float average_distances[2], min_distance

    # Downsample all the tracks
    downsample_tracks = [utils.downsample(t, nb_samples) for t in tracks]
    number_of_tracks = len(downsample_tracks)
    number_of_points = nb_samples * 3

    # Initialized the clustering
    if number_of_tracks == 0:
        raise ValueError("Expect some tracts as input of the clustering.")
    track = np.ascontiguousarray(downsample_tracks[0], dtype=f32_dt)
    cluster = <Cluster *>realloc(NULL, sizeof(Cluster))
    cluster[0].indices = <unsigned long *>realloc(NULL, sizeof(unsigned long))
    cluster[0].hidden= <float *>realloc(NULL, number_of_points * sizeof(float))
    cluster[0].indices[0] = 0
    cluster[0].N = 1
    ptr = <float *>track.data
    for i from 0 <= i < number_of_points:
        cluster[0].hidden[i] = ptr[i]
    number_of_clusters = 1
    
    # Associate in one wave each track to one cluster
    for fiber_it from 1 <= fiber_it < number_of_tracks by 1:

        # Generate a progress bar
        ratio = (fiber_it + 1) / <float>number_of_tracks
        progress = int(ratio * 100.)
        block = int(round(bar_length * ratio))
        text = "\rLocal skeleton clustering in progress: [{0}] {1}%".format(
            "=" * block + " " * (bar_length - block), progress)
        sys.stdout.write(text)
        sys.stdout.flush()

        # Get/format the current track   
        track = np.ascontiguousarray(downsample_tracks[fiber_it], dtype=f32_dt)
        ptr = <float *>track.data

        # Allocate objects to store distance between the current track and
        # the clusters             
        # > store all distance        
        alld = np.zeros(number_of_clusters, dtype=np.float32)
        # > store which fibers are fliped
        flip = np.zeros(number_of_clusters, dtype=np.int)  
        
        # Compute the track/clusters distances
        with nogil, cython.boundscheck(False), cython.wraparound(False):

            # Parallel loop
            for cluster_it in parallel.prange(number_of_clusters,
                                              schedule="static",
                                              num_threads=nb_of_threads):

                # Get/format the current cluster representative virtual fiber
                skeleton = <float *>realloc(
                    NULL, number_of_points * sizeof(float))
                for i from 0 <= i < number_of_points:
                    skeleton[i] = cluster[cluster_it].hidden[i]    

                 # Direct and flip average distance between two tracks
                average_direct_flip_distance(
                    ptr, skeleton, nb_samples, <float *>average_distances)

                # 'dist' should contains the smallest distance
                if average_distances[1] < average_distances[0]:                
                    flip[cluster_it] = 1
                    alld[cluster_it] = average_distances[1]
                else:             
                    alld[cluster_it] = average_distances[0]

                # Free memory
                free(skeleton)

            # Get the cluster that is the most representative of the current
            # fiber
            min_distance = biggest_float
            for i from 0 <= i < number_of_clusters:
                if alld[i] < min_distance:
                    min_distance = alld[i]
                    min_index = i

            # If the distance is small enough, attached the current fiber to
            # this cluster
            if min_distance < cutoff:

                # Upate the cluster skeleton
                # > correct if flipping is needed            
                if flip[min_index] == 1:
                    for i from 0 <= i < nb_samples:
                        for j from 0 <= j < 3:
                            cluster[min_index].hidden[i * 3 + j] = (
                                (cluster[min_index].N * 
                                 cluster[min_index].hidden[i * 3 + j] +
                                 ptr[(nb_samples - 1 - i) * 3 + j]) / 
                                 (cluster[min_index].N + 1))
                else:
                    for i from 0 <= i < number_of_points:
                        cluster[min_index].hidden[i] = (
                            (cluster[min_index].N *
                             cluster[min_index].hidden[i] + 
                             ptr[i]) / (cluster[min_index].N + 1))

                # Update the cluster elements
                cluster[min_index].N += 1
                cluster[min_index].indices = <unsigned long *>realloc(
                    cluster[min_index].indices,
                    cluster[min_index].N * sizeof(unsigned long))
                cluster[min_index].indices[cluster[min_index].N - 1] = fiber_it

            # Otherwise, create a new cluster that contains the current fiber 
            else:
                number_of_clusters += 1
                cluster = <Cluster *>realloc(
                    cluster, number_of_clusters * sizeof(Cluster))
                cluster[number_of_clusters - 1].indices= (
                    <unsigned long *>realloc(NULL, sizeof(unsigned long)))
                cluster[number_of_clusters - 1].hidden = (
                    <float *>realloc(NULL, number_of_points * sizeof(float)))
                cluster[number_of_clusters - 1].indices[0] = fiber_it
                cluster[number_of_clusters - 1].N = 1
                for i from 0 <= i < number_of_points:
                    cluster[number_of_clusters - 1].hidden[i] = ptr[i]

    # Copy results to a dictionary
    organized_cluster = {}
    for i from 0 <= i < number_of_clusters:
        organized_cluster[i] = {}
        #organized_cluster[i]["hidden"] = np.zeros(
        #    number_of_points, dtype=np.float32)
        #for j from 0 <= j < number_of_points:
        #    organized_cluster[i]["hidden"][j] = cluster[i].hidden[j]
        #organized_cluster[i]["hidden"].shape = (nb_samples, 3)
        organized_cluster[i]["N"] = cluster[i].N
        organized_cluster[i]["indices"] = []
        for j from 0 <= j < cluster[i].N:
            organized_cluster[i]["indices"].append(cluster[i].indices[j])

    # Free memory
    with nogil:
        for i from 0 <= i < number_of_clusters:
            free(cluster[i].indices)
            free(cluster[i].hidden)
        free(cluster)

    return organized_cluster


###############################################################################
# Intern functions
###############################################################################

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void average_direct_flip_distance(float *t1,
                                       float *t2,
                                       unsigned long nb_samples,
                                       float *out) nogil:
    """ Calculate the euclidean distance between two tracks of size nb_samples.
    Both direct and flip are given as output
    
    
    Parameters
    ----------------
    t1: float[]
        arrays representing the first track
    t2: float[]
        arrays representing the second track

     t1           t2
     0*     a     *0
      \           |
       \          |
       1*   b     *1
        |         \ 
        |          \
       2*   c       *2
    
    Returns
    -----------
    out: float[2]
        array containing the euclidean distance and the fliped euclidean
        distance.
    """
    cdef:
        float dist=0, flip_dist=0, tmp_dist, tmp_flip_dist, diff, flip_diff
        unsigned long i, j

    # Direct and flip average distance between two tracks
    for i from 0 <= i < nb_samples:
        tmp_dist = 0
        tmp_flip_dist = 0
        for j from 0 <= j < 3:
            diff = t1[i * 3 + j] - t2[i * 3 + j]
            flip_diff = t1[i * 3 + j] - t2[(nb_samples - 1 - i) * 3 + j]
            diff = diff * diff
            flip_diff = flip_diff * flip_diff 
            tmp_dist = tmp_dist + diff
            tmp_flip_dist = tmp_flip_dist + flip_diff
        dist += sqrt(tmp_dist)
        flip_dist += sqrt(tmp_flip_dist)

    out[0] = dist / <float>nb_samples
    out[1] = flip_dist / <float>nb_samples
