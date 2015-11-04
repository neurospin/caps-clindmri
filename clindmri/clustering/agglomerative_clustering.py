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
import time
import sys
import json
from scipy.cluster import hierarchy

# Clindmri import
from clindmri.clustering.metrics import mam_distances
from clindmri.clustering.metrics import clustering_metric


def fdist(fibers, cutoff, distfile):
    """ Compute the fiber to fiber condensed distance matrix using a geometric
    distance.

    Parameters
    ----------
    fibers: n size list of (Ni, 3) array
        the fibers to compare.
    cutoff: float 
        minimum threshold so that distances below it are not considered.
    distfile: str
        a file path where the computed condensed distance matrix will be
        saved.

    Returns
    -------
    dist: array (n(n-1)/2,)
        a condensed fiber to fiber distance matrix.
    """
    # Compute the fiber to fiber distances based on the estimated GPs
    dist = clustering_metric(fibers, distance_metric="min", zhang_thr=cutoff,
                             nb_of_threads=20)

    # Save the distance array
    with open(distfile, "w") as open_file:
        numpy.save(open_file, dist)

    return dist


def agglomerative_clustering(dist, clusterfile, linkage="single", cutoff=10.):
    """ Agglomerative clustering.

    Create fiber bundles using an agglomerative clustering strategy that
    recursively merges the pair of clusters that minimally increases
    a given linkage distance.

    This method has been proposed by S. Zhang in 'Identifying White-Matter
    Fiber Bundles in DTI Data Using an Automated Proximity-Based Fiber
    Clustering Method'.

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
    dist: array (n(n-1)/2,)
        a condensed fiber to fiber distance matrix.
    clusterfile: str
        a path to a file where the clustering result will be saved in Json
        format.
    linkage: str {'single', 'complete', 'average'} (default 'single')
        define which linkage criterion to use. It define the method used
        for calculating the distance between the newly formed cluster:
    cutoff:    float 
        the stopping criterion in the clustering.

    Returns
    -------
    organized_clusters: dict
        keys indentify clusters and contains numbers pointing to the original
        observed fiber positions.
    """
    # Want to keep a trace of the execution time, start timer
    now = time.time()
 
    # Currently, we only allow single, complete, or average linkage
    if linkage not in ["single", "complete", "average"]: 
        raise ValueError("Invalide linkage '{0}'.".format(linkage))

    # Start the agglomerative clustering
    linkage_matrix = hierarchy.linkage(dist, method=linkage)

    # Forms clusters given a threshold to apply when forming the clusters:
    # clusters[i] is the flat cluster number to which original observation
    # i belongs
    clusters = hierarchy.fcluster(linkage_matrix, cutoff, criterion="distance")
    print clusters

    # Organize the cluster result in a dictionary
    print numpy.where(clusters == 2)[0]
    print set(clusters)
    organized_clusters = dict(
        (str(x), (int(x), numpy.where(clusters == x)[0].tolist()))
        for x in set(clusters))

    # Save the result
    with open(clusterfile, "w") as open_file:
        json.dump(organized_clusters, open_file, indent=4)

    return organized_clusters

