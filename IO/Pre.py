#!/usr/bin/env python
# $Id: Pre.py 320 2014-08-25 02:28:50Z kato $

import wx
import vtk
from vtk.wx.wxVTKRenderWindow import *
from wxVTKRenderWindowInteractor import *
import numpy as n
import CGNSFile
import Roster
from SceneObject import SceneObject
from PatchSelectionPanel import *

class MainFrame(wx.Frame):

	def __init__(self, parent):
		wx.Frame.__init__(self, parent, -1, "Albastru Pre", size = (1024, 768))

		splitter = wx.SplitterWindow(self, wx.ID_ANY)

		# Scene
		self.vtkPanel = wxVTKRenderWindowInteractor(splitter, -1)
		ren = vtk.vtkRenderer()
		ren.SetBackground(0.9, 0.9, 1.0);
		self.vtkPanel.GetRenderWindow().AddRenderer(ren)
		Roster.Renderer = ren
		Roster.RenderWindow = self.vtkPanel.GetRenderWindow()

		# Patch list
		self.patchPanel = PatchSelectionPanel(splitter)

		splitter.SplitVertically(self.patchPanel, self.vtkPanel, 200)

class App(wx.App):

	MenuData = [
		["File", -1, [
			["&New", -1, None],
			["&Open", -1, None],
			["&Save", -1, None],
			["Save &as", -1, None],
			["&Quit", -1, None],
			]
		],
		["Edit", -1, [
			["Dummy", -1, None],
			]
		],
	]

	def OnInit(self):

		frame = MainFrame(None)

		menuBar = wx.MenuBar()
		self._InstallMenu(menuBar, self.MenuData)
		frame.SetMenuBar(menuBar)

		tb = frame.CreateToolBar((wx.TB_HORIZONTAL | wx.NO_BORDER))
		new_bmp = wx.ArtProvider.GetBitmap(wx.ART_NEW, wx.ART_TOOLBAR)
		open_bmp = wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN, wx.ART_TOOLBAR)
		save_bmp = wx.ArtProvider.GetBitmap(wx.ART_FILE_SAVE, wx.ART_TOOLBAR)
		tb.AddLabelTool(wx.NewId(), "New", new_bmp)
		tb.AddLabelTool(wx.NewId(), "Open", open_bmp)
		tb.AddLabelTool(wx.NewId(), "Save", save_bmp)
		tb.Realize()

		self.SetTopWindow(frame)
		frame.Show()

		return True

	def _InstallMenu(self, parent, data):

		for name, actionID, action in data:
			if isinstance(action, list):
				menu = wx.Menu()
				self._InstallMenu(menu, action)
				parent.Append(menu, name)
			else:
				parent.Append(actionID, name)

def CanonizeRange(r):

	rr = list(r)
	for i, j in [[0, 3], [1, 4], [2, 5]]:
		if rr[i] > rr[j]:
			rr[i], rr[j] = rr[j], rr[i]
	return rr

def PickEventCallback(obj, event):
	actor = obj.GetActor()
	if actor:
		picked = Roster.FindPatchByActor(actor)
		patchName = picked[1]
		print patchName
		Roster.UnselectAllPatches()
		Roster.SelectPatch(patchName)

def RenamePatch(patch):

	i = 2
	while True:
		candidateName = "%s_%d" % (patch["Name"], i)
		if Roster.FindPatchByName(candidateName) == None:
			print "Duplicate patch %s was renamed to %s" % (patch["Name"], candidateName)
			patch["Name"] = candidateName
			break
		i += 1

def CoarsenPatchMesh(xyz):

	shape = xyz.shape
	shapeCoarse = list(shape)
	stepI = max(1, shape[1] / 10)
	stepJ = max(1, shape[2] / 10)
	stepK = max(1, shape[3] / 10)
	xyzCoarse = xyz[:, ::stepI, ::stepJ, ::stepK]
	print "Coarsened:", xyz.shape, " -> ", xyzCoarse.shape
	return xyzCoarse

def Main(cgnsMeshFileName):

	app = App()
	#app.MainLoop()
	#return

	style = vtk.vtkInteractorStyleTrackballCamera()
	Roster.RenderWindow.GetInteractor().SetInteractorStyle(style)

	"""
	picker = vtk.vtkPropPicker()
	picker.AddObserver("EndPickEvent", PickEventCallback)
	Roster.RenderWindow.GetInteractor().SetPicker(picker)
	"""

	axesActor = vtk.vtkAxesActor()
	markerWidget = vtk.vtkOrientationMarkerWidget()
	markerWidget.SetOrientationMarker(axesActor)
	markerWidget.SetInteractor(Roster.RenderWindow.GetInteractor())

	if cgnsMeshFileName != None:
		cgnsf = CGNSFile.CGNSFile(cgnsMeshFileName)
		zones = cgnsf.ReadZones()
		Roster.Zones = zones

		#ren = vtk.vtkRenderer()
		ren = Roster.Renderer

		B = 1
		for zone in zones:

			xyz = cgnsf.ReadZoneCoord(B, zone["Zone"])

			for boco in zone["Bocos"]:
				if Roster.FindPatchByName(boco["Name"]) != None:
					RenamePatch(boco)
				r = CanonizeRange(boco["Range"])
				xyz2 = xyz[:, r[0] - 1:r[3], r[1] - 1:r[4], r[2] - 1:r[5]]
				xyzCoarse = CoarsenPatchMesh(xyz2)
				so = SceneObject(xyzCoarse)
				actors = so.ActorSurface, so.ActorOutline
				actors[0].GetProperty().SetColor(1.0, 1.0, 1.0)
				for actor in actors:
					ren.AddActor(actor)
				#Roster.Patches.append([zone, boco, actors])
				Roster.RegisterPatch(zone, boco, actors)

			if True:
				for c1to1 in zone["1to1s"]:
					if not c1to1.IsPeriodic():
						continue
					if Roster.FindPatchByName(c1to1["Name"]) != None:
						RenamePatch(c1to1)
						"""
						for i in range(2,100):
							candidateName = "%s_%d" % (c1to1["Name"], i)
							if Roster.FindPatchByName(candidateName) == None:
								print "Duplicate patch %s was renamed to %s" % (c1to1["Name"], candidateName)
								c1to1["Name"] = candidateName
								break
						"""
					r = CanonizeRange(c1to1["Range"])
					xyz2 = xyz[:, r[0] - 1:r[3], r[1] - 1:r[4], r[2] - 1:r[5]]
					xyzCoarse = CoarsenPatchMesh(xyz2)
					so = SceneObject(xyzCoarse)
					actors = so.ActorSurface, so.ActorOutline
					actors[0].GetProperty().SetColor(0.0, 0.0, 1.0)
					for actor in actors:
						ren.AddActor(actor)
					Roster.RegisterPatch(zone, c1to1, actors)

			del xyz

	markerWidget.SetEnabled(True)
	markerWidget.InteractiveOn()

	Roster.RosterModified()

	app.MainLoop()

if __name__ == "__main__":
	import sys
	if len(sys.argv) > 1:
		Main(sys.argv[1])
	else:
		Main(None)

