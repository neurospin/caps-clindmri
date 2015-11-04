#! /usr/bin/env python
##########################################################################
# NSAP - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import numpy
import os

os.environ["OMP_NUM_THREADS"] = "1"

# Pymc (Markov Chain Monte Carlo sampling toolkit) import gaussian process
# module
from pymc import gp

# Clindmri import
from .bounded_thinplate import innerproduct_thinplate3d


# Define 'pymc' new covariance using 'covariance_function_bundles' method:
# Init method takes three arguments:
#     1- cov_fun_name: The name of the covariance function.
#     2- cov_fun_module: The name of the module in which the covariance
#        function can be found. This must be somewhere on your PYTHONPATH.
#     3- extra_cov_params: A dictionary whose keys are the extra parameters 
#        the covariance function takes, and whose values are brief sentences
#        explaining the roles of those parameters. This dictionary will be used
#        to generate the docstring of the bundle.
cov_thinplate3d = gp.covariance_function_bundle(
    "thinplate3d_mat", "clindmri.tractography.bounded_thinplate",
    {"R": ("the radius of a shpere in R3 inside which the soultion is "
           "calculated.")})

# Lambda function to compute the fiber sampling deltas: if a determinist
# tractography with constant integration step has been used to produce the
# fiber, all the deltas are the same and are equal to the integration step.
fiber_deltas = lambda f: numpy.sqrt(((f[1:] - f[:-1]) ** 2).sum(1))


class FiberGP(object):
    """ Single fiber as a gaussian process (GP).

    The GP of a fiber ~ the blurred indicator function (BIF) of a fiber is 
    characterized by the BIF for a smooth trajectory as a GP with covariance
    Cs (geometric) and by the anisotrpoic blurring of a smooth function as a
    GP with covariance Cd (diffusion).

    Reference:
    D. Wassermann, L. Bloy, E. Kanterakis, R. Verma, and R. Deriche.
    Unsupervised white matter fiber clustering and tract probability map
    generation: Applications of a gaussian process framework for white matter
    fibers.
    NeuroImage, 51(1):228-241, 2010.
    """
    def __init__(self, fiber, rfactor=1, spacing=(1., 1., 1.)):
        """ Initialize the FiberGP class.

        Parameters
        ----------
        fiber: array (N, 3)
            a fiber representation in voxel coordiantes.
        rfactor: float
            the radius factor of a shpere in R3 inside which the soultion is
            calculated. It will be multiplied by the largest fiber sampling
            delta.
        spacing: 3-uplet
            the x, y, z spacing used to evaluate the fiber GP.
        """
        # Check that we have a 3d fiber object
        if fiber.shape[1] != 3 or fiber.ndim != 2:
            raise ValueError(
                "Unsupported fiber representation with shape '{0}', must be of "
                "the form '(N, 3)' where N is the number of fiber "
                "samples.".format(fiber.shape))

        # Set the geometric covariance function
        self._cs_pymc = cov_thinplate3d.euclidean
        self._cd_pymc = None

        # Store the fiber representation in a class parameter and compute the
        # largest fiber sampling delta
        self._fiber = fiber
        delta = fiber_deltas(fiber).max()

        # Compute the radius of the shpere in R3 inside which the soultion
        # will be calculated.
        self._r = float(rfactor * delta)

        # Store the GP sampling deltas in the x, y and z directions
        self._spacing = spacing

        # Define a bounding vox (bbox) around the fiber by stacking the lower
        # and upper fiber samples and declare an adjusted one that is
        # calculated from the search radius
        self.bbox = numpy.vstack((fiber.min(axis=0), fiber.max(axis=0)))
        self.adjusted_bbox = self.bbox.copy()
        self.adjusted_bbox[0] -= self._r
        self.adjusted_bbox[1] += self._r

        # Define the observed variance on the GP distribution GP(mean, cov).
        self._observed_variance = delta

        # Compute the smoothness GP gps and the diffusion-associated
        # blurring GP gpd
        self.gps = None
        self.gpd = None
        self.covs_matrix = None
        self.covd_matrix = None
        self.create_gps()  

    def create_gps(self):
        """ Create the smoothness GP gps and the diffusion-associated
        blurring GP gpd.

        ToDo: compute gpd
        """
        # Creating a mean function
        self._means = gp.Mean(constant_mean, val=0)

        # Creating a covairance function
        # The covariance function is multiplied by amp**2, and this effectively
        # multiplies realizations by amp. In other words, a larger amp
        # parameter means that realizations will deviate further from their
        # mean
        self._covs = gp.Covariance(self._cs_pymc, amp=1, R=self._r)
        self.covs_matrix = self._covs(self._fiber, self._fiber)

        # Some parameters needed to compute the inner product
        self._invcovs_matrix = numpy.linalg.inv(self.covs_matrix)
        self.alphas = numpy.asmatrix(
            numpy.dot(numpy.ones(len(self._fiber)), self._invcovs_matrix )).T

        # Normally-distributed observations on gaussian process distribution
        self.gps = gp.observe(
            self._means, self._covs, obs_mesh=self._fiber,
            obs_vals=numpy.ones(len(self._fiber)),
            obs_V=numpy.zeros(len(self._fiber)) + self._observed_variance)

    @classmethod
    def intersect_bboxes(self, bbox1, bbox2):
        """ Intersect two bounding boxes.

        Logical statement about how to determine if two bounding boxes
        overlap:

            1- the radius from the center of bounding box 1 to bounding
               box 2 in the direction d must be less than or equal to the sum
               of the radius of 1 in the d direction and the radius of 2 in
               the d direction, d in {x, y, z}.
    
        Parmaeters
        ----------
        bbox1, bbox2: array (2, 3)
            the bboxes to intersect.

        Retruns
        -------
        bbox: array (2, 3) or None
            the intersection between the two input bboxes. If None is returned,
            there is no intersection between the two bbox.
        """
        # Check if the two bboxes intersect themselves
        is_intersected = True
        for index in range(3):
            center_to_center = numpy.abs(
                bbox1[0, index] + bbox1[1, index] - bbox2[0, index] -
                bbox2[1, index])
            sum_radius = (
                bbox1[1, index] - bbox1[0, index] + bbox2[1, index] -
                bbox2[0, index])
            is_intersected = (
                is_intersected and (center_to_center <= sum_radius))

        # Compute the intersection
        if is_intersected:
            return numpy.vstack((numpy.vstack((bbox1[0],bbox2[0])).max(axis=0),
                                 numpy.vstack((bbox1[1],bbox2[1])).min(axis=0)))
        else:
            return None

    @classmethod
    def innerproduct_thinplate(self, fibergp1, fibergp2):
        """ Compute the inner product between two fiber GPs.

        Compute a similarity measure described by this inner product that is
        larger when the fibers are more similar and 0 when there is no overlap.

        Parameters
        ----------
        fibergp1, fibergp2: FiberGP
            two fiber GPs.

        Returns
        -------
        innerprod: float
            the similarity between two fibers in term of overlaps done by
            computing the inner product between two GPs.
        """
        # Get the two fibers meshes
        fiber1 = fibergp1._fiber
        fiber2 = fibergp2._fiber

        # Get the two resolution radius
        R = fibergp1._r
        Q = fibergp2._r

        # Compute the fiber point to point distance matrix
        difference = (fiber2[numpy.newaxis, ..., :] -
                      fiber1[..., numpy.newaxis, :])
        distmatrix = numpy.sqrt((difference ** 2).sum(axis=-1))      

        # Compute analytically the integral of the inner product
        if fiber1 is fiber2:
            innerproduct_thinplate3d(distmatrix, R, Q, symm=True)
        else:
            innerproduct_thinplate3d(distmatrix, R, Q, symm=False)

        # Compute the inner product by applying the fiber GPs covariances
        innerprod = float(numpy.dot(numpy.dot(fibergp1.alphas.T, distmatrix),
                          fibergp2.alphas))

        return innerprod

    @classmethod
    def innerproduct_thinplate_normalized(self, fibergp1, fibergp2):
        """ Compute the normalized inner product between two fiber GPs.

        Compute a similarity measure bounded by 1 (fibergp1 and fibergp2 are
        the same) and 0 (no intersection).

        Compute a similarity measure described by this inner product bounded
        by 1 (fibergp1 and fibergp2 are the same) and 0 (no intersection).

        Parameters
        ----------
        fibergp1, fibergp2: FiberGP
            two fiber GPs.

        Returns
        -------
        norm_innerprod: float
            the normalized similarity between two fibers: 1 = the same;
            0 = no overlap.
        """
        # Compute the GPs norms
        norm_fibergp1 = numpy.sqrt(
            FiberGP.innerproduct_thinplate(fibergp1, fibergp1))
        norm_fibergp2 = numpy.sqrt(
            FiberGP.innerproduct_thinplate(fibergp2, fibergp2))

        # Compute the GPs inner product
        innerprod = FiberGP.innerproduct_thinplate(fibergp1, fibergp2)

        # Compute the normalized inner product
        norm_innerprod = innerprod / (norm_fibergp1 * norm_fibergp2)

        return norm_innerprod

    @classmethod
    def fdist(self, fibers, rfactor, spacing, distfile):
        """ Compute the fiber to fiber condensed distance matrix using GPs.

        Parameters
        ----------
        fibers: n size list of (Ni, 3) array
            the fibers to compare.
        rfactor: float
            the radius factor of a shpere in R3 inside which the GPs soultions
            are calculated.
        spacing: 3-uplet
            the x, y, z image spacing used to evaluate the fiber GPs.
        distfile: str
            a file path where the computed condensed distance matrix will be
            saved.

        Returns
        -------
        dist: array (n(n-1)/2,)
            a condensed fiber to fiber distance matrix.
        """
        from clindmri.clustering.metrics import mam_distances

        # Method parameters
        nb_fibers = len(fibers)
        condensed_matrix_size = int(nb_fibers * (nb_fibers - 1) / 2)
        dist = numpy.zeros((condensed_matrix_size, ), dtype=numpy.single)

        # Compute all the GPs
        gps = [FiberGP(track, rfactor) for track in fibers]

        # Compute the fiber to fiber distances based on the estimated GPs
        for i in range(nb_fibers - 1):
            for j in range(i + 1, nb_fibers):
                index = i * (nb_fibers - 1.5 - 0.5 * i) + j - 1
                dist[index] = FiberGP.intersection_normalized(
                    gps[i], gps[j], spacing)

        # Save the distance array
        with open(distfile, "w") as open_file:
            numpy.save(open_file, dist)

        return dist

    @classmethod
    def intersection_normalized(self, fibergp1, fibergp2, shape):
        """ Compute the intersection between two fiber GPs.

        Compute an intersection measure bounded by 1 (fibergp1 and fibergp2 are
        the same) and 0 (no intersection).

        Parameters
        ----------
        fibergp1, fibergp2: FiberGP
            two fiber GPs.
        shape: 3-uplet
            the image shape.

        Returns
        -------
        norm_intersection: float
            the normalized intersection measure between two fibers:
            1 = the same; 0 = no overlap.
        """
        # Compute the GPs means
        mean1 = fibergp1.get_mean_field(shape=shape)
        mean2 = fibergp2.get_mean_field(shape=shape)

        # Compute the GPs intersection
        denominator = numpy.sqrt(numpy.sum(mean1**2) * numpy.sum(mean2**2))
        if denominator == 0:
            return 0.
        norm_intersection = numpy.sum(mean1 * mean2) / denominator

        return norm_intersection

    def get_variance_field(self, bbox=None, spacing=None, shape=None):
        """ Compute the GP estimated variance field on a sampling grid.

        Parameters
        ----------
        bbox: array (2, 3) (optional, default None)
            the bbox used to compute sampling grid. If None the class adjusted
            bbox is used.
        spacing: tuple (optional, default None)
            the sampling grid spacing. If None the class spacing is used.
        shape: tuple (3, ) (optional, default None)
            the image shape (must contains integer only). If specified, force
            the spacing to (1, 1, 1) and use an optimized evaluation that
            return a map of the parent shape.

        Returns
        -------
        cov_array: array (X, Y, Z)
            the evaluated covariances.
        """
        # If not specified set default parameters
        cov_array = None
        if spacing is None:
            spacing = self._spacing
        if bbox is None:
            bbox = self.adjusted_bbox
        if shape is not None:
            cov_array = numpy.ones(shape, dtype=numpy.single)
            spacing = (1, 1, 1)
            parent_bbox = numpy.array([
                [0, 0, 0], numpy.asarray(shape).astype(int) - 1])
            eval_bbox = numpy.array([numpy.floor(bbox[0]),
                                     numpy.ceil(bbox[1])])
            bbox = FiberGP.intersect_bboxes(eval_bbox, parent_bbox)
            if bbox is None:
                return cov_array

        # Compute the GP evaluation grid
        grid = numpy.mgrid[
            bbox[0, 0]: bbox[1, 0]: spacing[0],
            bbox[0, 1]: bbox[1, 1]: spacing[1],
            bbox[0, 2]: bbox[1, 2]: spacing[2]].T

        # Evaluate the covariance
        if len(grid) > 0:
            cov_values = self._covs(grid)
        else:
            return cov_array

        # If requested place the computed mean in the image space
        if shape is not None:
            grid = grid.astype(int)
            cov_array[
                grid[..., 0], grid[..., 1], grid[..., 2]] = cov_values
        else:
            cov_array = cov_values

        return cov_array

    def get_mean_field(self, bbox=None, spacing=None, shape=None):
        """ Compute the GP estimated mean field on a sampling grid.

        Parameters
        ----------
        bbox: array (2, 3) (optional, default None)
            the bbox used to compute sampling grid. If None the class adjusted
            bbox is used.
        spacing: tuple (optional, default None)
            the sampling grid spacing. If None the class spacing is used.
        shape: tuple (3, ) (optional, default None)
            the image shape (must contains integer only). If specified, force
            the spacing to (1, 1, 1) and use an optimized evaluation that
            return a map of the parent shape.

        Returns
        -------
        mean_array: array (X, Y, Z)
            the evaluated means.
        """
        # If not specified set default parameters
        mean_array = None
        if spacing is None:
            spacing = self._spacing
        if bbox is None:
            bbox = self.adjusted_bbox
        if shape is not None:
            mean_array = numpy.zeros(shape, dtype=numpy.single)
            spacing = (1, 1, 1)
            parent_bbox = numpy.array([
                [0, 0, 0], numpy.asarray(shape).astype(int) - 1])
            eval_bbox = numpy.array([numpy.floor(bbox[0]),
                                     numpy.ceil(bbox[1])])
            bbox = FiberGP.intersect_bboxes(eval_bbox, parent_bbox)
            if bbox is None:
                return mean_array

        # Compute the GP evaluation grid
        grid = numpy.mgrid[
            bbox[0, 0]: bbox[1, 0]: spacing[0],
            bbox[0, 1]: bbox[1, 1]: spacing[1],
            bbox[0, 2]: bbox[1, 2]: spacing[2]].T

        # Evaluate the mean
        if len(grid) > 0:
            mean_values = self._means(grid)
        else:
            return mean_array

        # If requested place the computed mean in the image space
        if shape is not None:
            grid = grid.astype(int)
            mean_array[
                grid[..., 0], grid[..., 1], grid[..., 2]] = mean_values
        else:
            mean_array = mean_values

        return mean_array


def constant_mean(x, val=0):
    """ Function that generates a constant mean that is used as a prior guess
    for GPs.
    """
    return numpy.zeros(x.shape[:-1], dtype=float) + val

