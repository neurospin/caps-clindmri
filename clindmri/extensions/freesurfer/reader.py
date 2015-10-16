#! /usr/bin/env python
##########################################################################
# NSAP - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import copy
import re
import numpy
import json
from nibabel import freesurfer

# Clindmri import
from .wrappers import FSWrapper
from .exceptions import FreeSurferRuntimeError


def read_cortex_surface_segmentation(fsdir, physical_to_index, fsconfig,
                                     affine=None):
    """ Read the cortex gyri surface segmentatation of freesurfer.

    Give access to the right and left hemisphere segmentations that can be
    projected on the cortical and inflated cortical surfaces.
    The vertex are expressed in the voxel coordinates.

    Parameters
    ----------
    fsdir: str( mandatory)
        the subject freesurfer segmentation directory.
    physical_to_index: array (mandatory)
        the transformation to project a physical point in an array.
    fsconfig: str (mandatory)
        the freesurfer configuration file.
    affine: array (optional, default None)
        an affine transformation in voxel coordinates that will be applied on
        the output vertex of the cortex surface.

    Returns
    -------
    segmentation: dict
        contain the two hemisphere 'lh' and 'rh' triangular surfaces and
        inflated surfaces represented in a TriSurface structure.
    """
    # Construct the path to the surface segmentation results and associated
    # labels
    meshdir = os.path.join(fsdir, "surf")
    labeldir = os.path.join(fsdir, "label")
    segfile = os.path.join(fsdir, "mri")

    # Get deformation between the ras and ras-tkregister spaces
    asegfile = os.path.join(segfile, "aseg.mgz")
    translation = tkregister_translation(asegfile, fsconfig)

    # Construct the deformation to apply on the cortex mesh
    if affine is None:
        affine = numpy.identity(4)
    deformation = numpy.dot(affine, numpy.dot(physical_to_index, translation))

    # Create an dictionary to contain all the surfaces and labels
    segmentation = {}

    # Select the hemisphere
    for hemi in ["lh", "rh"]:

        # Get annotation id at each vertex (if a vertex does not belong
        # to any label and orig_ids=False, its id will be set to -1) and
        # the names of the labels
        annotfile = os.path.join(labeldir, "{0}.aparc.annot".format(hemi))
        labels, ctab, regions = freesurfer.read_annot(
            annotfile, orig_ids=False)
        meta = dict((index, {"region": item[0], "color": item[1][:4].tolist()})
                    for index, item in enumerate(zip(regions, ctab)))

        # Select the surface type
        hemisegmentation = {}
        for surf in ["white", "inflated"]:

            # Load the mesh: a 2-uplet with vertex (x, y, z) coordinates and
            # mesh triangles
            meshfile = os.path.join(meshdir, "{0}.{1}".format(hemi, surf))
            mesh = freesurfer.read_geometry(meshfile)
            hemisegmentation[surf] = {
                "vertices": apply_affine_on_mesh(mesh[0], deformation),
                "triangles": mesh[1]
            }

        # Save the segmentation result
        segmentation[hemi] = TriSurface(
            vertices=hemisegmentation["white"]["vertices"],
            inflated_vertices=hemisegmentation["inflated"]["vertices"],
            triangles=hemisegmentation["white"]["triangles"],
            labels=labels,
            metadata=meta)

    return segmentation


def apply_affine_on_mesh(vertex, affine):
    """ Apply an affine transformation on each vetex of the mesh.

    Parameters
    ----------
    vertex: array (N, 3)
        N vertex.
    affine: array (4, 4)
        an affine transformation to applied.

    Results
    -------
    warp_vertex: array (N, 3)
        N interpolated vertex.
    """
    N, _ = vertex.shape
    ones = numpy.ones((N, 1), dtype=vertex.dtype)
    homogenous_vertex = numpy.concatenate((vertex, ones), axis=1)
    warp_vertex = numpy.dot(affine, homogenous_vertex.T).T[..., :3]
    return warp_vertex


def tkregister_translation(mgzfile, fsconfig):
    """ Get the tkregister translation.
  
    FreeSurfer use a special origin for the Right-Anterior-Superior
    (anatomical coordinates) space. To get the standard, freesurfer scanner
    space in RAS coordinates we can use the 'mri_info --vox2ras aseg.mgz' or
    'mri_info --vox2ras-trk aseg.mgz' commands respectively

    Parameters
    ----------
    mgzfile: str (mandatory)
        a FreeSurfer '.mgz' file.
    fsconfig: str (mandatory)
        the freesurfer configuration file.

    Returns
    -------
    translation: array
        the translation matrix betwwen the ras and ras-tkregister spaces.
    """
    # Get the affine matrices corresponding to the the ras or ras-tkregister
    # spaces
    affines = {}
    for tkregister in [True, False]:

        # Execute the FreeSurfer command
        command = ["mri_info", "--vox2ras", mgzfile]
        if tkregister:
            command[1] = "--vox2ras-tkr"
        fsprocess = FSWrapper(command, shfile=fsconfig)
        fsprocess()
        if fsprocess.exitcode != 0:
            raise FreeSurferRuntimeError(command[0], " ".join(command[1:]))

        # Get the affine matrix displayed in the stdout
        affine = fsprocess.stdout.splitlines()
        affine = ",".join([line.strip() for line in affine])
        affine = re.sub(r"  *", ",", affine)
        affine = numpy.fromstring(affine, dtype=float, sep=",").reshape(4, 4)
        affines[tkregister] = affine

    # Compute the translation
    translation = numpy.eye(4)
    translation += (affines[False] - affines[True])

    return translation


class TriSurface(object):
    """ A class representing a triangulated 3D surface.

    The surface contains 'ntrian' triangles, each having 3 vertices with
    3 coordinates. The surface can be initialized from one of the following.

    Parameters
    ----------
    vertices: array (nvert, 3)
        the mesh vertices.
    triangles: array (ntrian, 3)
        the indices of the vertices of each triangle.
    labels: array (nvert)
        a label associated to each vertex.
    metadata: dict
        a mapping between each label and associated 'color' and 'region' name.
    inflated_vertices: array (nvert, 3)
        the mesh inflated vertices.
    """

    def __init__(self, vertices, triangles, labels=None, metadata=None,
                 inflated_vertices=None):
        """ Create a new surface.
        """
        self.vertices = vertices
        self.triangles = triangles
        if labels is None:
            self.labels = numpy.asarray([0, ] * vertices.shape[0])
        else:
            self.labels = labels
        self.metadata = metadata
        self.inflated_vertices = inflated_vertices

    def save(self, outdir, outname):
        """ Export a mesh in freesurfer format.

        Parameters
        ----------
        outdir: str (mandatory)
            the location where the mesh will be written.
        outname: str (mandatory)
            the name of the file(s) that will be written.
        """
        meshfile = os.path.join(outdir, outname)
        freesurfer.write_geometry(meshfile, self.vertices, self.triangles)
        if self.inflated_vertices is not None:
            meshfile = os.path.join(outdir, outname + ".inflated")
            freesurfer.write_geometry(meshfile, self.inflated_vertices,
                                      self.triangles)

    @classmethod
    def load(self, meshfile, inflatedmeshpath=None, annotfile=None):
        """ Load a FreeSurfer surface.

        Parameters
        ----------
        meshfile: str (mandatory)
            the location of the file containing the FreeSurfer mesh to be
            loaded.
        inflatedmeshpath: str (mandatory)
            the location of the file containing the FreeSurfer inflated mesh
            to be loaded.

        Returns
        -------
        surf: TriSurface
            a triangular surface instance
        """
        vertices, triangles = freesurfer.read_geometry(meshfile)
        if inflatedmeshpath is not None:
            inflated_vertices, _triangles = freesurfer.read_geometry(
                inflatedmeshpath)
            if not numpy.allclose(triangles, _triangles):
                raise ValueError("'{0}' and '{1}' do not represent the same "
                                 "surface.".format(meshfile, inflatedmeshpath))
        else:
            inflated_vertices = None
        if annotfile is not None:
            labels, ctab, regions = freesurfer.read_annot(
                annotfile, orig_ids=False)
            meta = dict((index, {"region": item[0], "color": item[1][:4].tolist()})
                        for index, item in enumerate(zip(regions, ctab)))
        else:
            labels = None
            meta = None

        return TriSurface(vertices=vertices, triangles=triangles, labels=labels,
                          metadata=meta, inflated_vertices=inflated_vertices)

    def save_vtk(self, outfile, inflated=False):
        """ Export a mesh in .vtk format

        Parameters
        ----------
        outfile: str (mandatory)
            the location where the mesh will be written.
        inflated: bool (optional, default False)
            if True write the inflated volume.
        """
        import vtk

        print self.metadata

        # Check that the inflated
        if inflated and self.inflated_vertices is None:
            raise ValueError("Can't save inflated volume '{0}' since it has "
                             "not been specified.".format(outfile))

        # Create the desired polydata
        polydata = self._polydata(inflated=inflated)

        # Write the polydata
        writer = vtk.vtkPolyDataWriter()
        writer.SetDataModeToAscii()
        writer.SetFileName(outfile)
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(polydata)
        else:
            writer.SetInputData(polydata)
        writer.Write()

    def nedges(self):
        """ Using Euler's formula for triangle mesh return an approximation of
        the number of edges of the TriSurface.
        """
        return 3 * self.vertices.shape[0]

    def shape(self):
        """ TriSurface shape.

        Returns
        -------
        out: 3-uplet
            the number of points, edges, faces of the TriSurface.
        """
        return self.vertices.shape[0], self.nedges(), self.triangles.shape[0]

    def labelize(self, shape, shift=0):
        """ Compute a label image of the TriSurface.

        Parameters
        ----------
        shape: 3-uplet (mandatory)
            the image shape.
        shift: int (optional, default 0)
            shift the labels of this number.

        Returns
        -------
        label_array: array
            an array with the surface labels.
        nb_of_labels: int
            the number of valid labels.
        """
        label_array = numpy.zeros(shape, dtype=numpy.int16)
        indices = numpy.round(self.vertices)
        nb_of_labels = 0
        for label in set(self.labels):
            if label != -1:
                nb_of_labels += 1
                label_indices = indices[numpy.where(self.labels == label)]
                label_array[label_indices.T.tolist()] = label + shift
        return label_array, nb_of_labels
        

    def voxelize(self, shape, tol=0):
        """ Compute the enclosed points of the TriSurface.
        This code uses vtk.

        Parameters
        ----------
        shape: 3-uplet
            the image shape.

        Returns
        -------
        inside_array: array
            a mask array with the enclosed voxels.
        """
        import vtk
        from vtk.util.numpy_support import vtk_to_numpy

        # Construct the mesh grid from shape
        nx, ny, nz = shape
        gridx, gridy, gridz = numpy.meshgrid(numpy.linspace(0, nx - 1, nx),
                                             numpy.linspace(0, ny - 1, ny),
                                             numpy.linspace(0, nz - 1, nz))
       
        # Create polydata
        vtk_points = vtk.vtkPoints()
        for point in zip(gridx.flatten(), gridy.flatten(), gridz.flatten()):
            vtk_points.InsertNextPoint(point)
        points_polydata = vtk.vtkPolyData()
        points_polydata.SetPoints(vtk_points)
        surf_polydata = self._polydata()
        
        # Compute enclosed points
        enclosed_pts = vtk.vtkSelectEnclosedPoints()
        enclosed_pts.SetInput(points_polydata)
        enclosed_pts.SetTolerance(tol)
        enclosed_pts.SetSurface(surf_polydata)
        enclosed_pts.SetCheckSurface(1)
        enclosed_pts.Update()
        inside_points = enclosed_pts.GetOutput().GetPointData().GetArray(
            "SelectedPoints")
        enclosed_pts.ReleaseDataFlagOn()
        enclosed_pts.Complete()

        # Convert result as a numpy array
        inside_array = vtk_to_numpy(inside_points).reshape(ny, nx, nz)
        inside_array = numpy.swapaxes(inside_array, 1, 0)

        return inside_array

    def _polydata(self, inflated=False):
        """ Compute a vtk polydata of the TriSurface.
        This code uses vtk.

        Parameters
        ----------
        inflated: bool (optional, default False)
            If True use the inflated vertices.            

        Returns
        -------
        polydata: vtkPolyData
            the TriSurface vtk polydata.
        """
        import vtk

        # Select the vertices to use
        labels = copy.deepcopy(self.labels)
        if inflated:
            vertices = self.inflated_vertices
        else:
            vertices = self.vertices

        # First setup points, triangles and colors.
        vtk_points = vtk.vtkPoints()
        vtk_triangles = vtk.vtkCellArray()
        vtk_colors = vtk.vtkUnsignedCharArray()
        vtk_colors.SetNumberOfComponents(1)
        nb_of_labels = len(set(self.labels))
        labels[numpy.where(labels < 0)] = 0
        for index in range(len(vertices)):
            vtk_points.InsertNextPoint(vertices[index])
            vtk_colors.InsertNextTuple1(labels[index])
        for triangle in self.triangles:
            vtk_triangle = vtk.vtkTriangle()
            vtk_triangle.GetPointIds().SetId(0, triangle[0])
            vtk_triangle.GetPointIds().SetId(1, triangle[1])
            vtk_triangle.GetPointIds().SetId(2, triangle[2])
            vtk_triangles.InsertNextCell(vtk_triangle)

        # Create (geometry and topology) the associated polydata
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(vtk_points)
        polydata.GetPointData().SetScalars(vtk_colors)
        polydata.SetPolys(vtk_triangles)

        return polydata

            
