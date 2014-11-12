# $Id: SceneObject.py 216 2012-04-06 03:50:25Z kato $

import vtk
import numpy as n

class SceneObject(object):

	def __init__(self, xyz):
		"""
		xyz contains a surface point set that defines the patch of interest.
		xyz itself is a 3-D array of shape (3, :, :, :) but one of the spatial
		index range is one because it's a surface.
		"""

		self.ActorSurface, self.ActorOutline = MakeActors(xyz)

def MakeActors(xyz):

	npts = xyz.size / xyz.shape[0]

	floats = vtk.vtkFloatArray()
	floats.SetNumberOfComponents(3)
	floats.SetNumberOfTuples(npts)

	data = xyz.reshape((3, npts), order = "F").transpose()
	for i, p in enumerate(data):
		floats.SetTuple3(i, p[0], p[1], p[2])

	points = vtk.vtkPoints()
	points.SetData(floats)

	grid = vtk.vtkStructuredGrid()
	grid.SetDimensions(xyz.shape[1], xyz.shape[2], xyz.shape[3])
	grid.SetPoints(points)

	polydataSurface = vtk.vtkStructuredGridGeometryFilter()
	polydataSurface.SetInput(grid)

	mapperSurface = vtk.vtkPolyDataMapper()
	mapperSurface.SetInput(polydataSurface.GetOutput())

	actorSurface = vtk.vtkActor()
	actorSurface.SetMapper(mapperSurface)
	actorSurface.GetProperty().SetOpacity(0.3)

	polydataOutline = vtk.vtkStructuredGridOutlineFilter()
	polydataOutline.SetInput(grid)
	mapperOutline = vtk.vtkPolyDataMapper()
	mapperOutline.SetInput(polydataOutline.GetOutput())

	actorOutline = vtk.vtkActor()
	actorOutline.SetMapper(mapperOutline)
	actorOutline.GetProperty().SetColor(0.0, 0.0, 0.0)

	return actorSurface, actorOutline

