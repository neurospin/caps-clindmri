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
import types
import logging

# Caps import
from clindmri.estimation.gdti.monomials import construct_matrix_of_monomials
from .colors import *

# VTK import
try:
    import vtk
except ImportError:
    raise ImportError("VTK is not installed.")


def ren():
    """ Create a renderer

    Returns
    --------
    ren: vtkRenderer() object

    Examples
    --------
    >>> import plot_vtk
    >>> ren = plot_vtk.ren()
    >>> plot_vtk.add(ren, actor)
    >>> plot_vtk.show(ren)
    """
    return vtk.vtkRenderer()


def add(ren, actor):
    """ Add a specific actor
    """
    if isinstance(actor, vtk.vtkVolume):
        ren.AddVolume(actor)
    else:
        ren.AddActor(actor)


def rm(ren, actor):
    """ Remove a specific actor
    """
    ren.RemoveActor(actor)


def clear(ren):
    """ Remove all actors from the renderer
    """
    ren.RemoveAllViewProps()


def show(ren, title="pvtk", size=(300, 300)):
    """ Show window

    Parameters
    ----------
    ren : vtkRenderer() object
        as returned from function ren()
    title : string
        a string for the window title bar
    size : (int, int)
        (width,height) of the window
    """
    ren.ResetCameraClippingRange()

    window = vtk.vtkRenderWindow()
    window.AddRenderer(ren)
    window.SetWindowName(title)
    window.SetSize(size)

    style = vtk.vtkInteractorStyleTrackballCamera()
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(window)
    iren.SetInteractorStyle(style)
    iren.Initialize()

    window.Render()
    iren.Start()


def tensor(coeff, order, position=(0, 0, 0),
           radius=0.5, thetares=20, phires=20, opacity=1, tessel=0):
    """ Generate a generic tensor actor.
    """
    # Create a sphere that we will deform
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(radius)
    sphere.SetLatLongTessellation(tessel)
    sphere.SetThetaResolution(thetares)
    sphere.SetPhiResolution(phires)

    # Get the polydata
    poly = sphere.GetOutput()
    poly.Update()

    # Get the mesh
    numPts = poly.GetNumberOfPoints()
    mesh = numpy.zeros((numPts, 3), dtype=numpy.single)
    for i in range(numPts):
        mesh[i, :] = (poly.GetPoint(i)[0], poly.GetPoint(i)[1],
                      poly.GetPoint(i)[2])

    # Deform mesh
    design_matrix = construct_matrix_of_monomials(mesh, order)
    signal = numpy.dot(design_matrix, coeff)
    #signal = np.maximum(signal, 0.0)
    signal /= signal.max()
    signal *= 0.5

    scalars = vtk.vtkFloatArray()
    pts = vtk.vtkPoints()
    pts.SetNumberOfPoints(numPts)
    for i in range(numPts):
        pts.SetPoint(i, signal[i] * mesh[i, 0], signal[i] * mesh[i, 1],
                     signal[i] * mesh[i, 2])
        scalars.InsertTuple1(i, signal[i])

    poly.SetPoints(pts)
    poly.GetPointData().SetScalars(scalars)
    poly.Update()

    lut = vtk.vtkLookupTable()
    lut.SetHueRange(0.667, 0.0)
    lut.Build()

    spherem = vtk.vtkPolyDataMapper()
    spherem.SetInput(poly)
    spherem.SetLookupTable(lut)
    spherem.ScalarVisibilityOn()
    spherem.SetColorModeToMapScalars()
    spherem.SetScalarRange(0.0, 0.5)

    actor = vtk.vtkActor()
    actor.SetMapper(spherem)
    actor.SetPosition(position)
    actor.GetProperty().SetOpacity(opacity)

    return actor


def line(lines, scalar, lut=None, opacity=1, linewidth=1):
    """ Create a line actor for one or more lines.    
    
    Parameters
    ----------
    lines : list
        a list of array representing a line as 3d points (N, 3)
    scalar : a float   
        0 <= scalar <= 1 to associate a color to the bloc of lines.         
    opacity : float (default = 1)
        the transparency of the bloc of lines: 0 <= transparency <= 1.
    linewidth : float (default = 1)
        the line thickness.              
    
    Returns
    ----------
    actor: vtkActor
        the bloc of lines actor.    
    """
    # Consteruct le lookup table if necessary
    if lut is None:
		lut = vtk.vtkLookupTable()
		lut.SetHueRange(0.667, 0.0)
		lut.Build()
    
    # If one line is passed as a numpy array, create virtually a list around
    # this structure
    if not isinstance(lines, types.ListType):
        lines = [lines]  
      
    # Create vtk structures
    vtk_points = vtk.vtkPoints()
    vtk_line = vtk.vtkCellArray()
    vtk_scalars = vtk.vtkFloatArray()
  
    # Go through all lines for the rendering
    point_id = 0
    for line in lines:

        # Get the line size
        nb_of_points, line_dim = line.shape

        # Associate one scalar to each point of the line for color rendering 
        scalars = [scalar] * nb_of_points

        # Fill the vtk structure for the curretn line
        for point_position in range(nb_of_points - 1):

            # Get the segment [p0, p1]
            p0 = line[point_position] 
            p1 = line[point_position + 1]

            # Set line points
            vtk_points.InsertNextPoint(p0)               
            vtk_points.InsertNextPoint(p1)

            # Set color property
            vtk_scalars.SetNumberOfComponents(1)
            vtk_scalars.InsertNextTuple1(scalars[point_position])
            vtk_scalars.InsertNextTuple1(scalars[point_position])

            # Set line segment         
            vtk_line.InsertNextCell(2)
            vtk_line.InsertCellPoint(point_id)
            vtk_line.InsertCellPoint(point_id + 1)
                    
            point_id += 2            

    # Create the line polydata
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetLines(vtk_line)
    polydata.GetPointData().SetScalars(vtk_scalars)
  
    # Create the line mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInput(polydata)
    mapper.SetLookupTable(lut)
    mapper.SetColorModeToMapScalars() 
    mapper.SetScalarRange(0.0, 1.0)  
    mapper.SetScalarModeToUsePointData()
       
    # Create the line actor
    actor=vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetOpacity(opacity)
    actor.GetProperty().SetLineWidth(linewidth)
    
    return actor


def dots(points, color=(1,0,0), psize=1, opacity=1):
    """ Create one or more 3d dot points.

    Returns
    -------
    actor: vtkActor
        one actor handling all the points.
    """
    if points.ndim==2:
        points_no=points.shape[0]
    else:
        points_no=1

    polyVertexPoints = vtk.vtkPoints()
    polyVertexPoints.SetNumberOfPoints(points_no)
    aPolyVertex = vtk.vtkPolyVertex()
    aPolyVertex.GetPointIds().SetNumberOfIds(points_no)

    cnt=0
    if points.ndim>1:
        for point in points:
            polyVertexPoints.InsertPoint(cnt, point[0], point[1], point[2])
            aPolyVertex.GetPointIds().SetId(cnt, cnt)
            cnt+=1
    else:
        polyVertexPoints.InsertPoint(cnt, points[0], points[1], points[2])
        aPolyVertex.GetPointIds().SetId(cnt, cnt)
        cnt+=1

    aPolyVertexGrid = vtk.vtkUnstructuredGrid()
    aPolyVertexGrid.Allocate(1, 1)
    aPolyVertexGrid.InsertNextCell(aPolyVertex.GetCellType(), aPolyVertex.GetPointIds())

    aPolyVertexGrid.SetPoints(polyVertexPoints)
    aPolyVertexMapper = vtk.vtkDataSetMapper()
    aPolyVertexMapper.SetInput(aPolyVertexGrid)
    aPolyVertexActor = vtk.vtkActor()
    aPolyVertexActor.SetMapper(aPolyVertexMapper)

    aPolyVertexActor.GetProperty().SetColor(color)
    aPolyVertexActor.GetProperty().SetOpacity(opacity)
    aPolyVertexActor.GetProperty().SetPointSize(psize)
    return aPolyVertexActor


def surface(points, triangles, labels, ctab, opacity=1):
    """ Create a colored triangular surface.

    Parameters
    ----------
    points: array (n_vertices, 3)
        the surface vertices.
    triangles: array
        nfaces x 3 array defining mesh triangles.
    labels: array (n_vertices)
        Annotation id at each vertex.
        If a vertex does not belong to any label its id must be negative.
    ctab: ndarray (n_labels, 5)
        RGBA + label id color table array.
    opacity: float
        the actor global opacity.

    Returns
    -------
    actor: vtkActor
        one actor handling the surface.
    """
    # First setup points, triangles and colors
    vtk_points = vtk.vtkPoints()
    vtk_triangles = vtk.vtkCellArray()
    vtk_colors = vtk.vtkUnsignedCharArray()
    vtk_colors.SetNumberOfComponents(1)
    nb_of_labels = len(set(labels))
    labels[numpy.where(labels < 0)] = 0
    for index in range(len(points)):
        vtk_points.InsertNextPoint(points[index])
        vtk_colors.InsertNextTuple1(labels[index])
    for cnt, triangle in enumerate(triangles):
        vtk_triangle = vtk.vtkTriangle()
        vtk_triangle.GetPointIds().SetId(0, triangle[0])
        vtk_triangle.GetPointIds().SetId(1, triangle[1])
        vtk_triangle.GetPointIds().SetId(2, triangle[2])
        vtk_triangles.InsertNextCell(vtk_triangle)

    # Make a lookup table using vtkColorSeries
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfColors(nb_of_labels + 1)
    lut.Build()
    for cnt, lut_element in enumerate(ctab):
        lut.SetTableValue(cnt, lut_element[0] / 255., lut_element[1] / 255.,
                          lut_element[2] / 255., lut_element[3] / 255.)
    lut.SetNanColor(1, 0, 0, 1)

    # Create (geometry and topology) the associated polydata
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.GetPointData().SetScalars(vtk_colors)
    polydata.SetPolys(vtk_triangles)

    # Create the mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInput(polydata)
    mapper.SetLookupTable(lut)
    mapper.SetColorModeToMapScalars()
    mapper.SetScalarRange(0, nb_of_labels)
    mapper.SetScalarModeToUsePointData()
    
    # Create the actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetOpacity(opacity)

    return actor
    
