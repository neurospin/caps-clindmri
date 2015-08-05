#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
########################################################################## 

import os
import json
import numpy as np
import clindmri.plot.pvtk as pvtk
from clindmri.clustering import local_skeleton_clustering
from clindmri.clustering import agglomerative_clustering
from clindmri.tractography import Tractogram


np.set_printoptions(precision=4)
key = 3


tracks = [np.array([[0,0,0],[1,0,0,],[2,0,0]], dtype=np.single),            
          np.array([[3,0,0],[3.5,1,0],[4,2,0]], dtype=np.single),
          np.array([[3.2,0,0],[3.7,1,0],[4.4,2,0]], dtype=np.single),
          np.array([[3.4,0,0],[3.9,1,0],[4.6,2,0]], dtype=np.single),
          np.array([[0,0.2,0],[1,0.2,0],[2,0.2,0]], dtype=np.single),
          np.array([[2,0.2,0],[1,0.2,0],[0,0.2,0]], dtype=np.single),
          np.array([[0,0,0],[0,1,0],[0,2,0]], dtype=np.single)]


if key == 0:
    clusters = agglomerative_clustering(
        tracks, "single", cutoff=0.5, nb_samples=10, nb_cpu=1,
        preserve_input=False)
    ren = pvtk.ren()
    nb_of_clusters = len(clusters)
    for cluster_id, tracks_index in clusters.iteritems():
        fibers = []
        for index in tracks_index:
            fibers.append(tracks[index])
        actor = pvtk.line(fibers, float(cluster_id) / nb_of_clusters)
        pvtk.add(ren, actor)
    pvtk.show(ren)

elif key==1:
    tracks *= 1000
    print len(tracks)
    print clustering_metric_resample(np.asarray(tracks), nb_of_threads=24)
    print clustering_metric(tracks, nb_of_threads=24, zhang_thr=0.,
                            distance_metric="max")
    
elif key==2:
    tracks *= 30
    clusters = local_skeleton_clustering(tracks, cutoff=0.7, nb_samples=10,
                                         nb_of_threads=1)
    print clusters

    # Plot
    ren = pvtk.ren()
    nb_of_clusters = float(len(clusters))
    for cl_id, cl in clusters.iteritems():
        fibers = []
        for fiber_ind in cl["indices"]:
            fibers.append(tracks[fiber_ind])
        actor = pvtk.line(fibers, float(cl_id) / nb_of_clusters)
        pvtk.add(ren, actor)
    pvtk.show(ren)

elif key==3:
    outdir = "/volatile/imagen/dmritest/000000022453/processed"
    tractogram = Tractogram(os.path.join(outdir, "000000022453s011a1001-reduce.trk"))
    cluster = os.path.join(outdir, "000000022453s011a1001-cluster.json")
    if 1:
        clusters = local_skeleton_clustering(tractogram.tracks, cutoff=8,
                                             nb_samples=10, nb_of_threads=12)
    else:
        clusters = agglomerative_clustering(
            tractogram.tracks, "single", cutoff=4, nb_samples=10, nb_cpu=12,
            preserve_input=False)
    json.dump(clusters, open(cluster, "w"), indent=4)
    print "nb clusters = ", len(clusters)
    ren = pvtk.ren()
    nb_of_clusters = float(len(clusters))
    for cl_id, cl in clusters.iteritems():
        fibers = []
        for fiber_ind in cl["indices"]:
            fibers.append(tractogram.tracks[fiber_ind])
        actor = pvtk.line(fibers, np.random.rand())
        pvtk.add(ren, actor)
    pvtk.show(ren)
    pvtk.clear(ren)
    
