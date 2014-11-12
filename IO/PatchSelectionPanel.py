# $Id: PatchSelectionPanel.py 288 2013-07-17 00:19:30Z kato $

import wx
import wx.lib.scrolledpanel as scrolled
import Roster

class PatchSelectionPanel(scrolled.ScrolledPanel):

	def __init__(self, parent):
		scrolled.ScrolledPanel.__init__(self, parent, -1)

		#self.ListBox = wx.ListBox(self, -1)
		self.Tree = wx.TreeCtrl(self, -1, style = wx.TR_DEFAULT_STYLE | wx.TR_HIDE_ROOT)
		self.PopulateTree()

		sizer = wx.BoxSizer(wx.VERTICAL)
		#sizer.Add(self.ListBox, 1, wx.EXPAND, 0)
		sizer.Add(self.Tree, 1, wx.EXPAND, 0)

		self.SetSizer(sizer)
		self.Layout()

		Roster.Observers.append(self)

		#self.Bind(wx.EVT_LISTBOX, self.EvtListBox, self.ListBox)

	def RosterModified(self):

		patchNames = map(lambda p: [p["Name"], p], Roster.GetBocoPatches())
		patchNames.sort(lambda a, b: cmp(a[0], b[0]))
		self.Tree.DeleteChildren(self.PatchRoot)
		for name, patch in patchNames:
			item = self.Tree.AppendItem(self.PatchRoot, name)
			self.Tree.SetPyData(item, ["Patch", name])

		c1to1Names = map(lambda p: [p["Name"], p], Roster.Get1to1Patches())
		c1to1Names.sort(lambda a, b: cmp(a[0], b[0]))
		self.Tree.DeleteChildren(self.C1to1Root)
		for name, patch in c1to1Names:
			print "C1to1: ", name
			item = self.Tree.AppendItem(self.C1to1Root, name)
			self.Tree.SetPyData(item, ["Patch", name])

		self.Tree.DeleteChildren(self.BlockRoot)
		for zone in Roster.Zones:
			item = self.Tree.AppendItem(self.BlockRoot, "(%d)%s" % (zone["Zone"], zone["Name"]))
			self.Tree.SetPyData(item, ["Block", zone["Name"]])

	"""
	def EvtListBox(self, event):
		name = event.GetString()
		print name
		#patch = Roster.FindPatchByName(name)
		#print patch[0:2]
		Roster.UnselectAllPatches()
		Roster.SelectPatch(name)
		Roster.RenderWindow.Render()
	"""

	def PopulateTree(self):

		root = self.Tree.AddRoot("root")
		self.BlockRoot = self.Tree.AppendItem(root, "Blocks")
		self.PatchRoot = self.Tree.AppendItem(root, "Patches")
		self.C1to1Root = self.Tree.AppendItem(root, "1to1 Connectivities")
		self.Bind(wx.EVT_TREE_SEL_CHANGED, self.OnTreeSelChanged, self.Tree)

	def OnTreeSelChanged(self, event):

		item = event.GetItem()
		print item
		label = self.Tree.GetItemText(item)
		data = self.Tree.GetPyData(item)
		if data == None:
			return
		elif data[0] == "Patch":
			Roster.UnselectAllPatches()
			Roster.SelectPatch(data[1])
		elif data[0] == "Block":
			Roster.UnselectAllPatches()
			Roster.SelectZone(data[1])
		Roster.RenderWindow.Render()

