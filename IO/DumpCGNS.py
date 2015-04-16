#!/usr/bin/env python

import CGNS
import CGNSFile
import numpy as n

def Main(filename):

	f = CGNSFile.CGNSFile(filename)

	families = f.ReadFamilies()
	print "*** Families ***"
	for key, value in families.items():
		print "Family: ", key
		print "  ", value

	zones = f.ReadZones()
	print "*** Zones ***"
	for zone in zones:
		print "Zone: ", zone["Zone"]
		print "Name: ", zone["Name"]
		print "Size: ", zone["Size"]
		print "Family: ", zone["FamilyName"]
		print "Bocos:"
		for boco in zone["Bocos"]:
			print "  Name: \"%s\"" % boco["Name"]
			print "  Type: ", boco["Type"]
			print "  Range: ", boco["Range"]
			if boco["Type"] == "FamilySpecified":
				print "  FamilyName: ", boco["FamilyName"]
			print
		print "1to1s:"
		for c1to1 in zone["1to1s"]:
			print "  Name: ", c1to1["Name"]
			print "  DonorZoneName: ", c1to1["DonorZoneName"]
			print "  Range (self): ", c1to1["Range"]
			print "  Range (donor): ", c1to1["DonorRange"]
			print "  Transform (self to donor): ", c1to1["Transform"]
			print "  RotationCenter: ", c1to1["RotationCenter"]
			print "  RotationAngle: ", c1to1["RotationAngle"]
			print "  Translation: ", c1to1["Translation"]
			print
		print "CNMBs:"
		for conn in zone["CNMBs"]:
			print "  Name: ", conn["Name"]
			print "  DonorZoneName: ", conn["DonorZoneName"]
			print "  Range (self): ", conn["Range"]
			#print "  donorData: ", conn["donorData"]
			print "  DonorPatch: ", conn["DonorPatch"]
			print "  ConnectionId: ", conn["ConnectionId"]
		print

if __name__ == "__main__":
	import sys
	Main(sys.argv[1])

