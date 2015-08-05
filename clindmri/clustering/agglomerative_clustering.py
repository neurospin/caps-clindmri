#! /usr/bin/env python
##########################################################################
# CAPS - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
# 
# From dipy: http://nipy.org/dipy/
##########################################################################

# System import
from __future__ import division
import numpy
import multiprocessing
import itertools
import time
import sys
import logging
from scipy.cluster import hierarchy

# Define the logger
logger = logging.getLogger(__name__)

# Dmri import
from clindmri.tractography.utils import downsample
from clindmri.clustering.metrics import mam_distances
from clindmri.clustering.metrics import clustering_metric_resample
from clindmri.clustering.metrics import clustering_metric


class Cluster(object):
    """ Structure that represent a fiber cluster.

    Attributes
    ----------
    tracks_index: list of int
        the index of the tracks involved in the cluster.
    linkage: str
        define the agglomerative strategy: 
        single linkage agglomerative clustering 's',
        complete linkage agglomerative clustering 'c',
        or average linkage agglomerative clustering 't'.
    z_thr: float (default 0.5)
        threshold used in the fiber to fiber metric that will enables us to
        differentiate fibers that only differe in a relatively small portion
        of there tracks.
    centroid: numpy array
        the mean track of the cluster.
    _tracks: list of numpy array
        usually a downsample version of the original track used to compute the
        metric.
    """
    def __init__(self, tracks_index, linkage, z_thr, tracks):
        """ Initialize a new cluster.

        Parameters
        ----------
        tracks_index: list of int
            the index of the tracks involved in the cluster.
        linkage: str
            define the agglomerative strategy: 
            single linkage agglomerative clustering 's',
            complete linkage agglomerative clustering 'c',
            or average linkage agglomerative clustering 't'.
        z_thr: float (default 0.5)
            threshold used in the fiber to fiber metric that will enables us to
            differentiate fibers that only differe in a relatively small portion
            of there tracks.
        tracks: list of numpy array
            used to initilize the fiber clustering metric
        """
        # Empty cluster don't make sense!
        if len(tracks_index) <= 0:
            raise ValueError("Illegal empty cluster initilization.")

        # Define class parameters
        self.tracks_index = tracks_index
        self._tracks = tracks
        self.linkage = linkage
        self.z_thr = z_thr
        self.centroid = None

        # Figure out what the centroid of this Cluster should be only for the
        # average linkage
        if linkage == "t":
            self.centroid = self.calculate_centroid()  
 
    def calculate_centroid(self):
        """ Calculates the centroid fiber.
        """
        s, si = most_similar_track_mam(self._tracks, metric="avg")
        return s

    def get_distance(self, cluster):
        """ Calculates the appropriate linkage distance between this and another
        cluster.
        """
        if self.linkage == "s":
            return self.get_single_lnk_distance(cluster)
        elif self.linkage == "c":
            return self.get_complete_lnk_distance(cluster)
        else:
            return self.get_centroid_distance(cluster)       

    def get_single_lnk_distance(self, cluster):
        """ Calculates the single-linkage distance between this and another
        cluster.
        """
        min_distance = float("inf")
        for Q in self._tracks:
            for R in cluster._tracks:
                distance = mam_distances(Q, R, self.z_thr, metric="max")
                if distance < min_distance: min_distance = distance
        return min_distance
 
    def get_complete_lnk_distance(self, cluster):
        """ Calculates the complete-linkage distance between this and another
        cluster.
        """
        max_distance = float("-inf")
        for Q in self._tracks:
            for R in cluster._tracks:
                distance = mam_distances(Q, R, self.z_thr, metric="max")
                if distance > max_distance: max_distance = distance
        return max_distance

    def get_centroid_distance(self, cluster):
        """ Calculates the centroid-linkage distance between this and another
        cluster.
        """
        return mam_distances(self._tracks[self.centroid],
                             cluster._tracks[cluster.centroid],
                             self.z_thr, metric="max")

    def fuse(self, cluster):
        """ Calculates the fusion of this and another cluster.
        """
        tracks_index = self.tracks_index
        intern_tracks = self._tracks
        tracks_index.extend(cluster.tracks_index)
        intern_tracks.extend(cluster._tracks)
        return Cluster(tracks_index, self.linkage, self.z_thr, intern_tracks)


def make_cluster_distance_matrix(clusters, nb_cpu=None):
    """ Compute the distance matrix capturing the distances between all
    clusters.

    Parameters
    ----------
    clusters: list of Cluster
        all the clusters we want to process.
    nb_cpu: int (default None - cpu_count())
        the number of worker processes to use

    Returns
    -------
    metric: numpy array
        (N, N) array containing the cluster to cluster distances.    
    """
    # Create a full square distance matrix
    # ToDo: only create triangular matrix to save memory since it is symetric
    nb_of_clusters = len(clusters)
    metric = numpy.ones((nb_of_clusters, nb_of_clusters), dtype=numpy.single)
    metric *= float("inf")

    # Compute a list of paired index to compute the fiber metric on.
    # Deal only with unique elements
    cluster_pairs = []
    for i in range(nb_of_clusters):
        for j in range(i + 1, nb_of_clusters):
            cluster_pairs.append((i, j))

    # Compute the fiber metric: create a group of CPUs to run on
    pool = multiprocessing.Pool(processes=nb_cpu) 
    async_result = pool.map_async(
        make_distance,
        [(x, clusters) for x in cluster_pairs])
    result_list = async_result.get()
    for pair, fdist in result_list:
        metric[pair[0], pair[1]] = metric[pair[1], pair[0]] = fdist

    return metric


def make_distance(params):
    """ Computes distance between two clusters
    """
    # Unpack parameters
    pair, clusters = params
    c1 = clusters[pair[0]]
    c2 = clusters[pair[1]]

    # Get the appropriate distance
    return (pair, c1.get_distance(c2))


def agglomerative_clustering(tracks, linkage="single", cutoff=10.,
                             nb_samples=10, nb_cpu=1, preserve_input=False):
    """ Agglomerative clustering.

    Create fiber bundles using an agglomerative clustering strategy that
    recursively merges the pair of clusters that minimally increases
    a given linkage distance.

    This method has been proposed by S. Zhang in 'Identifying White-Matter Fiber
    Bundles in DTI Data Using an Automated Proximity-Based Fiber Clustering
    Method'.

    Different linkage:

     * linkage='single'

       .. math::
          d(u,v) = \min(dist(u[i],v[j]))

       for all points :math:`i` in cluster :math:`u` and
       :math:`j` in cluster :math:`v`. This is also known as the
       Nearest Point Algorithm.

     * linkage='complete'

       .. math::
          d(u, v) = \max(dist(u[i],v[j]))

       for all points :math:`i` in cluster u and :math:`j` in
       cluster :math:`v`. This is also known by the Farthest Point
       Algorithm or Voor Hees Algorithm.

     * linkage='average'

       .. math::
          d(u,v) = \sum_{ij} \frac{d(u[i], v[j])}
                                  {(|u|*|v|)}

       for all points :math:`i` and :math:`j` where :math:`|u|`
       and :math:`|v|` are the cardinalities of clusters :math:`u`
       and :math:`v`, respectively. This is also called the UPGMA
       algorithm.

    Parameters
    ----------
    tracks: list of numpy array
        the list of fibers we want ot cluster.
    linkage: str {'single', 'complete', 'average'} (default 'single')
        define which linkage criterion to use. It define the method used
        for calculating the distance between the newly formed cluster:

    cutoff:    float 
        the stopping criterion in the clustering.
    nb_samples: float (default 10)
        the number of points in the downsample tracks used in order to speed 
        up the compututation of the fiber to fiber metric.
    nb_cpu: int (default 1)
        the number of threads to use.
    preserve_input: bool (default False)
        specifies whether the method makes a working copy of the input tracks 
        or writes resampled data into the existing stricture.

    Returns
    -------
    organized_clusters: dict
        keys indentify clusters and contains numbers pointing to the original
        observation positions.
    """
    # Want to keep a trace of the execution time, start timer
    now = time.time()
 
    # Currently, we only allow single, complete, or average linkage
    if linkage not in ["single", "complete", "average"]: 
        raise ValueError("Invalide linkage '{0}'.".format(linkage))

    # Udersample each track
    if preserve_input:
        undersample_tracks = []
    for indx in range(len(tracks)):
        if preserve_input:
            undersample_tracks.append(downsample(tracks[indx], nb_samples))
        else:
            tracks[indx] = downsample(tracks[indx], nb_samples)
    if preserve_input:
        undersample_tracks = numpy.asarray(undersample_tracks).astype(numpy.single)
    else:
        undersample_tracks = numpy.asarray(tracks).astype(numpy.single)

    # Compute a distance matrix from the input tracks
    fdist = clustering_metric_resample(undersample_tracks, nb_of_threads=nb_cpu)
    #fdist = clustering_metric(undersample_tracks, nb_of_threads=nb_cpu)
    #fdist = make_cluster_distance_matrix(clusters, nb_cpu)
    logger.debug("Done in {0} s.".format(time.time() - now))
    logger.debug("Metric:\n{0}".format(fdist))

    # Start the agglomerative clustering
    linkage_matrix = hierarchy.linkage(fdist, method=linkage)
    logger.debug("Done in {0} s.".format(time.time() - now))
    logger.debug("Linkage Matrix:\n{0}".format(linkage_matrix))

    # Foms clusters given a threshold to apply when forming the clusters:
    # clusters[i] is the flat cluster number to which original observation
    # i belongs
    clusters = hierarchy.fcluster(linkage_matrix, cutoff, criterion="distance")
    logger.debug("Done in {0} s.".format(time.time() - now))
    logger.debug("Clusters:\n{0}".format(clusters))

    # Organize the cluster result in a dictionary
    organized_clusters = dict((x, numpy.where(clusters == x)[0])
                              for x in set(clusters))
    logger.debug("Done in {0} s.".format(time.time() - now))
    logger.debug("Organized Clusters:\n{0}".format(organized_clusters))

    return organized_clusters

