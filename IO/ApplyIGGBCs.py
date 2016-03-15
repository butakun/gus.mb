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
# $Id: ApplyIGGBCs.py 320 2014-08-25 02:28:50Z kato $

import argparse
import CGNSFile, CGNS, PointRange
import IGGBCS

def Main(cgnsFileName, bcsFileName, namingScheme):

	cgnsFile = CGNSFile.CGNSFile(cgnsFileName, "w")

	zones = cgnsFile.ReadZones()
	allBocos = zones.AllBocos()

	for zone in zones:
		print "Zone %d, %d bc patches" % (zone["Zone"], len(zone["Bocos"]))
		for boco in zone["Bocos"]:
			print "  %s" % boco.Path

	print "Reading IGG bcs file %s" % bcsFileName
	iggbcs = IGGBCS.IGGBCS(bcsFileName)

	hubNextID = 1
	shroudNextID = 1
	bladeNextID = 1
	inletNextID = 1
	outletNextID = 1
	interfaceNextID = 1
	unknownNextID = 1

	print "Converting patch names"
	bocosRenamed = []
	for iggblock in iggbcs.Blocks:
		blockName = iggblock["Name"]
		zone = zones.FindZoneByName(blockName)
		assert(zone != None)
		print "* found block %s at Z = %d" % (blockName, zone["Zone"])
		for iggPatch in iggblock["Patches"]:

			iggname = iggPatch["Name"]
			iggtype = iggPatch["Type"]
			if namingScheme == "ag5":
				if "hub" in iggname:
					iggPatch["CGNSPatchName"] = "Hub_%03d" % hubNextID
					hubNextID += 1
				elif "shroud" in iggname:
					iggPatch["CGNSPatchName"] = "Shroud_%03d" % shroudNextID
					shroudNextID += 1
				elif "skin" in iggname or "blunt" in iggname:
					iggPatch["CGNSPatchName"] = "Blade_%03d" % bladeNextID
					bladeNextID += 1
				elif iggtype == "INL":
					iggPatch["CGNSPatchName"] = "Inlet_%03d" % inletNextID
					inletNextID += 1
				elif iggtype == "OUT":
					iggPatch["CGNSPatchName"] = "Outlet_%03d" % outletNextID
					outletNextID += 1
				elif iggtype == "ROT":
					iggPatch["CGNSPatchName"] = "Interface_%03d" % interfaceNextID
					interfaceNextID += 1
				else:
					iggPatch["CGNSPatchName"] = "Unknown_%03d" % unknownNextID
					unknownNextID += 1
			elif namingScheme == "igg":
				iggPatch["CGNSPatchName"] = iggname
			else:
				raise ValueError
			print "* * IGG bc \"%s\" will be renamed to %s" % (iggname, iggPatch["CGNSPatchName"])

			iggPointRange = PointRange.PointRange2D(iggPatch["Start"], iggPatch["End"])
			print "* * trying to find patch %s %s" % (iggPatch["Name"], str(iggPointRange))
			bocoFound = None
			for boco in zone["Bocos"]:
				if boco.PointRange == iggPointRange:
					bocoFound = boco
					break
			if bocoFound != None:
				print "* * * found cgns patch %s %s" % (boco["Name"], boco.PointRange)
				CGNS.GoPath(cgnsFile.fn, boco.ParentPath())
				CGNS.DeleteNode(boco.LeafName())
				bocosRenamed.append(boco)
			B = 1
			CGNS.BocoWrite(cgnsFile.fn, B, zone["Zone"], iggPatch["CGNSPatchName"], iggPatch["Type"], iggPointRange.Sextuplets())

	# AG5 might have created duplicate boco names (multiple bocos using the same name).
	# We detect them and rename them so the boco names are all unique globally.
	count = 1
	bocosToRename = []
	for boco in allBocos:
		print boco["Name"]
		#if not boco in bocosRenamed:
		renamed = filter(lambda b: b["Zone"]["Zone"] == boco["Zone"]["Zone"] and b["Name"] == boco["Name"], bocosRenamed)
		if len(renamed) == 0:
			CGNS.GoPath(cgnsFile.fn, boco.ParentPath())
			CGNS.DeleteNode(boco.LeafName())
			newName = "Unknown_%03d" % count
			bocosToRename.append([boco, newName])
			count += 1
	for boco, newName in bocosToRename:
		print boco["Zone"]["Zone"], boco["Type"], boco["Range"]
		CGNS.BocoWrite(cgnsFile.fn, B, boco["Zone"]["Zone"], newName, boco["Type"], list(boco["Range"]))
		print "Renamed %s to %s" % (boco["Name"], newName)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "Apply IGG bcs info into a CGNS file")
	parser.add_argument("cgnsFile")
	parser.add_argument("bcsFile")
	parser.add_argument("--naming", default = "ag5", choices = ["ag5", "igg"])

	args = parser.parse_args()

	Main(args.cgnsFile, args.bcsFile, args.naming)

