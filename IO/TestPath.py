# $Id: TestPath.py 104 2011-06-03 04:21:30Z kato $

import CGNSFile
import CGNS

def PointRangeSextuplets(key, zone):

	imax, jmax, kmax = zone["Size"][:3]
	if key == "IMIN":
		ptRange = [1, 1, 1, 1, jmax, kmax]
	elif key == "IMAX":
		ptRange = [imax, 1, 1, imax, jmax, kmax]
	elif key == "JMIN":
		ptRange = [1, 1, 1, imax, 1, kmax]
	elif key == "JMAX":
		ptRange = [1, jmax, 1, imax, jmax, kmax]
	elif key == "KMIN":
		ptRange = [1, 1, 1, imax, jmax, 1]
	elif key == "KMAX":
		ptRange = [1, 1, kmax, imax, jmax, kmax]
	return ptRange

def Main(filename):

	f = CGNSFile.CGNSFile(filename, "w")

	baseName = "fluid"

	zones = f.ReadZones()

	# Delete all the BC patches.
	for zone in zones:
		bocoNames = [ boco["Name"] for boco in zone["Bocos"] ]
		print bocoNames
		CGNS.GoPath(f.fn, "/%s/%s/ZoneBC/" % (baseName, zone["Name"]))
		for bocoName in bocoNames:
			CGNS.DeleteNode(bocoName)

	hubPatches = [ [1, "JMIN"], [2, "JMIN"], [3, "JMIN"], [4, "JMIN"], [5, "JMIN"], [6, "KMIN"], [7, "KMIN"] ]
	shroudPatches = [ [1, "JMAX"], [2, "JMAX"], [3, "JMAX"], [4, "JMAX"], [5, "JMAX"], [6, "KMAX"], [7, "KMAX"] ]
	bladePatches = [ [2, "KMIN"], [3, "KMAX"], [5, "IMAX"] ]
	inletPatches = [ [6, "IMAX"] ]
	outletPatches = [ [7, "IMIN"] ]

	B = 1

	for stem, patches in [["Hub", hubPatches], ["Shroud", shroudPatches], ["Blade", bladePatches], ["Inlet", inletPatches], ["Outlet", outletPatches] ]:
		for i, ZKey in enumerate(patches):
			Z, key = ZKey
			zone = zones[Z - 1]
			assert(zone["Zone"] == Z)
			ptRange = PointRangeSextuplets(key, zone)
			CGNS.BocoWrite(f.fn, B, Z, "%s_%d" % (stem, i), "Dummy", ptRange)

	zones = f.ReadZones()

	for zone in zones:
		Z = zone["Zone"]
		"""
		if Z == 2:
			patches = zone.UnassignedPatches()
		continue
		"""
		patches = zone.UnassignedPatches()

		for patch in patches:
			print "Zone ", Z, ", orphans = ", patch.Min, patch.Max
		for i, patch in enumerate(patches):
			name = "%d_%d" % (Z, i + 1)
			ptRange = patch.Sextuplets()
			CGNS.BocoWrite(f.fn, B, Z, name, "Unknown", ptRange)

if __name__ == "__main__":
	Main("./test.cgns")

