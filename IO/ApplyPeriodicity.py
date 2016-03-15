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
#!/usr/bin/env python

import CGNS
import CGNSFile
import numpy as n

class Patch(object):

	def __init__(self, name, xyz, range_, boco):
		self.Name = name
		self.XYZ = xyz
		self.Begin1 = n.array(range_[0])
		self.End1 = n.array(range_[1])
		self.Boco = boco

	def Get2D(self):
		xyz = self.XYZ
		s = xyz.shape
		if s[1] == 1:
			return xyz[:, 0, :, :]
		elif s[2] == 1:
			return xyz[:, :, 0, :]
		elif s[3] == 1:
			return xyz[:, :, :, 0]
		raise Error

	def PatchCorners(self):
		pp = self.Get2D()
		pp = n.array([pp[:, 0, 0], pp[:, -1, 0], pp[:, 0, -1], pp[:, -1, -1]])
		return pp

	def CornerIndices(self):
		xyz = self.XYZ
		s = xyz.shape
		if s[1] == 1:
			return [0, 0, 0], [0, -1, 0], [0, 0, -1], [0, -1, -1]
		elif s[2] == 1:
			return [0, 0, 0], [-1, 0, 0], [0, 0, -1], [-1, 0, -1]
		elif s[3] == 1:
			return [0, 0, 0], [-1, 0, 0], [0, -1, 0], [-1, -1, 0]
		raise Error

	def CornerIndices2(self):
		xyz = self.XYZ
		s = xyz.shape
		if s[1] == 1:
			return n.array([[0, 0, 0], [0, s[2] - 1, 0], [0, 0, s[3] - 1], [0, s[2] - 1, s[3] - 1]]) + self.Begin1
		elif s[2] == 1:
			return n.array([[0, 0, 0], [s[1] - 1, 0, 0], [0, 0, s[3] - 1], [s[1] - 1, 0, s[3] - 1]]) + self.Begin1
		elif s[3] == 1:
			return n.array([[0, 0, 0], [s[1] - 1, 0, 0], [0, s[2] - 1, 0], [s[1] - 1, s[2] - 1, 0]]) + self.Begin1
		raise Error

	def PatchDirection(self):
		""" 0 and 1 are the axes spanning the patch surface, 2 is the axis normal to the patch surface """
		xyz = self.XYZ
		s = xyz.shape
		if s[1] == 1:
			return n.array([2, 0, 1])
		elif s[2] == 1:
			return n.array([0, 2, 1])
		elif s[3] == 1:
			return n.array([0, 1, 2])
		raise Error

	def BoundingBox(self):
		xyz = self.XYZ
		xmin = n.min(xyz[0, :, :, :])
		xmax = n.max(xyz[0, :, :, :])
		ymin = n.min(xyz[1, :, :, :])
		ymax = n.max(xyz[1, :, :, :])
		zmin = n.min(xyz[2, :, :, :])
		zmax = n.max(xyz[2, :, :, :])
		return n.array([[xmin, ymin, zmin], [xmax, ymax, zmax]])

	def AngleTo(self, patch2, axis):
		bb = self.BoundingBox()
		bb2 = patch2.BoundingBox()
		theta1 = AngleBetween(bb[0], bb2[0], axis)
		theta2 = AngleBetween(bb[1], bb2[1], axis)
		return theta1, theta2

def GatherBocoPatches(f, zones):

	patches = []

	B = 1
	for zone in zones:
		xyzZone = f.ReadZoneCoord(B, zone["Zone"])
		for boco in zone["Bocos"]:
			bocoType = boco["Type"]
			if bocoType == "Outflow" or bocoType == "Inflow" or bocoType == "Wall":
				continue
			r = boco["Range"]
			xyzPatch = xyzZone[:, r[0]-1:r[3], r[1]-1:r[4], r[2]-1:r[5]]
			patch = Patch(boco["Name"], xyzPatch, [[r[0] - 1, r[1] - 1, r[2] - 1], [r[3] - 1, r[4] - 1, r[5] - 1]], boco)
			patches.append(patch)

	return patches

def RotateCorners(pp, axis, theta):

	# FIXME
	assert(n.all(axis == n.array([0.0, 0.0, 1.0])))

	c, s = n.cos(theta), n.sin(theta)
	m = n.array([[c, s], [-s, c]])
	pp2 = n.array(pp)
	for i in range(4):
		pp2[i, :2] = n.dot(m, pp[i, :2])
	return pp2

def CornerDistances(pp1, pp2):

	dd = n.zeros(4)
	ii = []
	for i in range(4):
		dmin, jmin = 1.0e30, -1
		for j in range(4):
			d = pp1[i, :] - pp2[j, :]
			d = n.dot(d, d)
			if d < dmin:
				dmin, jmin = d, j
		dd[i] = dmin
		ii.append(jmin)
	return dd, ii

def TransformFullToShort(full):

	short = n.array([1, 2, 3])
	for i in range(3):
		d = n.nonzero(full[:, i])[0][0] + 1
		if full[d - 1, i] < 0:
			d = -d
		short[i] = d
	return short

def IsPeriodic(patch1, patch2, axis, theta, distMax):

	p1s = patch1.PatchCorners()
	p2s = patch2.PatchCorners()
	p2s1 = RotateCorners(p2s, axis, theta)
	p2s2 = RotateCorners(p2s, axis, -theta)
	dd1, ii1 = CornerDistances(p1s, p2s1)
	dd2, ii2 = CornerDistances(p1s, p2s2)
	dd1max, dd2max = n.max(dd1), n.max(dd2)
	if dd1max < dd2max:
		dd, ii, ddmax, dtheta = dd1, ii1, dd1max, -theta
	else:
		dd, ii, ddmax, dtheta = dd2, ii2, dd2max, theta

	if ddmax <= distMax:
		return True, dd, ii, dtheta
	else:
		return False, dd, ii, dtheta

def MatchPatches(periodicity):

	patch1 = periodicity["Patch1"]
	patch2 = periodicity["Patch2"]
	ii = periodicity["Indices1to2"]
	icorners1 = n.array(patch1.CornerIndices2())
	icorners2 = n.array(patch2.CornerIndices2())
	dcorners1 = n.abs(n.array(patch1.CornerIndices()))
	dcorners2 = n.abs(n.array(patch2.CornerIndices()))

	Begin1 = icorners1[0]
	Begin2 = icorners2[ii[0]]

	# Patch surface
	i11 = dcorners1[1] - dcorners1[0]
	i12 = dcorners2[ii[1]] - dcorners2[ii[0]]
	i21 = dcorners1[2] - dcorners1[0]
	i22 = dcorners2[ii[2]] - dcorners2[ii[0]]
	assert(len(n.nonzero(i11)) == 1)
	transform = n.zeros((3, 3), int)
	for ii1, ii2 in [[i11, i12], [i21, i22]]:
		i1 = n.nonzero(ii1)[0][0]
		i2 = n.nonzero(ii2)[0][0]
		sgn = ii1[i1] * ii2[i2]
		transform[i2, i1] = sgn

	# Patch normal
	i1 = n.where(patch1.PatchDirection() == 2)[0][0]
	i2 = n.where(patch2.PatchDirection() == 2)[0][0]
	ix1 = patch1.Begin1[i1]
	ix2 = patch2.Begin1[i2]
	if (ix1 == 0 and ix2 == 0) or (ix1 != 0 and ix2 != 0):
		transform[i2, i1] = -1
	else:
		transform[i2, i1] = 1

	print "Matching Patch %s and Patch %s" % (patch1.Name, patch2.Name)
	print "complete transform = \n", transform

	transformShort = TransformFullToShort(transform)
	print "shorthand transform = ", transformShort

	one = n.array([1, 1, 1])
	begin2 = n.dot(transform, patch1.Begin1 - patch1.Begin1) + patch2.Begin1
	end2 = n.dot(transform, patch1.End1 - patch1.Begin1) + patch2.Begin1
	print patch1.Begin1 + one, " -> ", begin2 + one
	print patch1.End1 + one, " -> ", end2 + one
	print patch1.Boco["Range"], " -> ", patch2.Boco["Range"]

	periodicity["Transform"] = transformShort

def FindPeriodicity(f, zones, axis, theta, distMax):

	patches = GatherBocoPatches(f, zones)

	npatches = len(patches)

	periodicities = []

	for i in range(npatches):
		for j in range(npatches): #range(i + 1, npatches):
			if i == j:
				continue
			periodic, dd, ii, dtheta = IsPeriodic(patches[i], patches[j], axis, theta, distMax)
			if periodic:
				print patches[i].Name, patches[j].Name, dd, ii, n.degrees(dtheta)
				periodicities.append({"Patch1":patches[i], "Patch2":patches[j], "Indices1to2":ii, "Axis":axis, "Angle1to2":dtheta})

	for periodicity in periodicities:
		print
		MatchPatches(periodicity)

	return periodicities

def FindMinimalCellSize(xyz, dsmin):

	dxyz = xyz[:, 1:, 1:, 1:] - xyz[:, :-1, :-1, :-1]
	dxyz = n.abs(dxyz)
	dsmin_ = n.min(dxyz)
	return min(dsmin_, dsmin)

def FindBocoByName(zones, name):

	for zone in zones:
		bocos = zone["Bocos"]
		boco = filter(lambda bc: bc["Name"] == name, bocos)
		if boco:
			boco = boco[0]
			return zone, boco
	raise Error

def Main(filename):

	f = CGNSFile.CGNSFile(filename, "w")
	zones = f.ReadZones()

	B = 1

	dsmin = 1.0e30
	for zone in zones:
		print "Zone: ", zone["Zone"]
		print "Name: ", zone["Name"]
		print "Size: ", zone["Size"]
		print "Bocos:"
		for boco in zone["Bocos"]:
			print "  Name: ", boco["Name"]
			print "  Type: ", boco["Type"]
			print "  Range: ", boco["Range"]
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
		print

		xyz = f.ReadZoneCoord(B, zone["Zone"])
		dsmin = FindMinimalCellSize(xyz, dsmin)
		del xyz

	print "dsmin = ", dsmin

	print
	print "Applying the periodic interfaces"
	print
	NBLADES = 33
	one = n.array([1, 1, 1])
	periodicities = FindPeriodicity(f, zones, n.array([0.0, 0.0, 1.0]), 2.0 * n.pi / float(NBLADES), dsmin * 0.99)
	for periodicity in periodicities:
		print
		patch1, patch2 = periodicity["Patch1"], periodicity["Patch2"]
		zone1, boco1 = FindBocoByName(zones, patch1.Name)
		zone2, boco2 = FindBocoByName(zones, patch2.Name)
		print zone1["Name"], boco1["Name"], " <-> ",  zone2["Name"], boco2["Name"]

		connectName = patch1.Name + "_1to1"
		begin1 = patch1.Begin1 + one
		end1 = patch1.End1 + one
		selfRange = []
		selfRange.extend(begin1)
		selfRange.extend(end1)
		begin2 = patch2.Begin1 + one
		end2 = patch2.End1 + one
		donorRange = []
		donorRange.extend(begin2)
		donorRange.extend(end2)
		transform = list(periodicity["Transform"])
		print connectName, zone1["Name"], zone2["Name"], selfRange, donorRange, transform

		I = CGNS.C1to1Write(f.fn, 1, zone1["Zone"], connectName, zone2["Name"], selfRange, donorRange, transform)
		print I
		# FIXME
		center = [0.0, 0.0, 0.0]
		axis = list(periodicity["Axis"] * periodicity["Angle1to2"])
		translation = [0.0, 0.0, 0.0]
		CGNS.C1to1PeriodicWrite(f.fn, 1, zone1["Zone"], I, center, axis, translation)

if __name__ == "__main__":
	import sys
	Main(sys.argv[1])

