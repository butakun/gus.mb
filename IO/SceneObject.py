#    gus.mb, an open source flow solver.
#    Copyright (C) 2016 Hiromasa Kato <hiromasa at gmail.com>
#
#    This file is part of gus.mb.
#
#    gus.mb is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    gus.mb is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
	#polydataSurface.SetInputData(grid)
	polydataSurface.SetInputData(grid)

	mapperSurface = vtk.vtkPolyDataMapper()
	#mapperSurface.SetInputData(polydataSurface.GetOutput())
	mapperSurface.SetInputConnection(polydataSurface.GetOutputPort())

	actorSurface = vtk.vtkActor()
	actorSurface.SetMapper(mapperSurface)
	actorSurface.GetProperty().SetOpacity(0.3)

	polydataOutline = vtk.vtkStructuredGridOutlineFilter()
	#polydataOutline.SetInputData(grid)
	polydataOutline.SetInputData(grid)
	mapperOutline = vtk.vtkPolyDataMapper()
	#mapperOutline.SetInputData(polydataOutline.GetOutput())
	mapperOutline.SetInputConnection(polydataOutline.GetOutputPort())

	actorOutline = vtk.vtkActor()
	actorOutline.SetMapper(mapperOutline)
	actorOutline.GetProperty().SetColor(0.0, 0.0, 0.0)

	return actorSurface, actorOutline

