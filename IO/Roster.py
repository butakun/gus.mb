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
# $Id: Roster.py 234 2012-05-25 15:48:32Z kato $

import CGNSFile

Renderer = None
RenderWindow = None

Zones = None

# patch = patchInfo with "ZoneInfo" and "Actors" keys added.
Patches = []

# {zoneName:patches} where patches are as in Patches above.
ZoneNameMap = {}

Observers = []

def SelectPatch(name):

	patch = FindPatchByName(name)
	actors = patch["Actors"]
	actors[0].GetProperty().SetColor(1.0, 0.0, 0.0)
	actors[0].GetProperty().SetOpacity(0.8)

def SelectZone(zoneName):

	patches = ZoneNameMap[zoneName]
	for patch in patches:
		SelectPatch(patch["Name"])

def UnselectAllPatches():

	for patch in Patches:
		actors = patch["Actors"]
		actors[0].GetProperty().SetColor(1.0, 1.0, 1.0)
		actors[0].GetProperty().SetOpacity(0.3)

def FindPatchByName(name):

	matched = filter(lambda p: p["Name"] == name, Patches)
	if len(matched) == 1:
		return matched[0]
	elif len(matched) > 1:
		print "Warning: there are multiple patches named ", name
		return matched[0]
	else:
		return None

def FindPatchByActor(actor):

	matched = filter(lambda p: actor in p["Actors"], Patches)
	assert(len(matched) == 1)
	return matched[0]

def GetBocoPatches():

	return filter(lambda p: isinstance(p, CGNSFile.PatchBoco), Patches)

def Get1to1Patches():

	return filter(lambda p: isinstance(p, CGNSFile.Patch1to1), Patches)

def RegisterPatch(zoneInfo, patchInfo, actors):

	patchInfo["ZoneInfo"] = zoneInfo
	patchInfo["Actors"] = actors
	Patches.append(patchInfo)

	if not zoneInfo["Name"] in ZoneNameMap:
		zonePatches = []
		ZoneNameMap[zoneInfo["Name"]] = zonePatches
	else:
		zonePatches = ZoneNameMap[zoneInfo["Name"]]
	zonePatches.append(patchInfo)

def RosterModified():

	for observer in Observers:
		observer.RosterModified()

