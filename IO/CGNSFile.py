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
# $Id: CGNSFile.py 314 2013-12-24 15:47:30Z kato $

import numpy as np
import math as m
import fnmatch
import CGNS
from PointRange import *

class Zones(list):

	def AllBocos(self):

		""" deprecated """
		assert(false)
		patches = []
		for zone in self:
			patches.extend(zone["Bocos"])
		return patches

	def FindBocoByName(self, name):

		""" deprecated """
		assert(false)
		for zone in self:
			match = filter(lambda b: b["Name"] == name, zone["Bocos"])
			assert(len(match) <= 1)
			if match:
				return zone, match[0]
		return None

	def FindZoneByID(self, Z):

		match = filter(lambda z: z["Zone"] == Z, self)
		assert(len(match) <= 1)
		if match:
			return match[0]
		return None

	def FindZoneByName(self, name):

		match = filter(lambda z: z["Name"] == name, self)
		assert(len(match) <= 1)
		if match:
			return match[0]
		return None

class Node(dict):

	def __init__(self, path):
		self.Path = path

	def LeafName(self):
		names = self.Path.split('/')
		return names[-1]

	def ParentPath(self):
		names = self.Path.split('/')
		return reduce(lambda a, b: a + '/' + b, names[:-1])

class Base(Node):

	def __init__(self, baseName):
		Node.__init__(self, "/%s" % baseName)

class Zone(Node):

	def __init__(self, base, zone, name, size):
		Node.__init__(self, "%s/%s" % (base.Path, name))
		self["Zone"] = zone
		self["Name"] = name
		self["Size"] = size
		self["Bocos"] = None
		self["1to1s"] = None

	def UnassignedPatches(self):

		imax, jmax, kmax = np.array(self["Size"][0:3]) - 1
		print imax, jmax, kmax

		pr1 = PointRange2D([0, 0, 0], [imax, jmax, 0])
		pr2 = PointRange2D([0, 0, kmax], [imax, jmax, kmax])
		pr3 = PointRange2D([0, 0, 0], [imax, 0, kmax])
		pr4 = PointRange2D([0, jmax, 0], [imax, jmax, kmax])
		pr5 = PointRange2D([0, 0, 0], [0, jmax, kmax])
		pr6 = PointRange2D([imax, 0, 0], [imax, jmax, kmax])
		prs = [pr1, pr2, pr3, pr4, pr5, pr6]

		for patches in [self["Bocos"], self["1to1s"]]:
			for patch in patches:
				r = np.array(patch["Range"]) - 1
				lo = np.minimum(r[0:3], r[3:6])
				hi = np.maximum(r[0:3], r[3:6])
				pr1to1 = PointRange2D(lo, hi)
				for pr in prs:
					pr.Subtract(pr1to1)

		patches = []
		for pr in prs:
			subpatches = pr.SplitToRects()
			patches.extend(subpatches)
		return patches

	def FindBocoByName(self, name):

		match = filter(lambda b: b["Name"] == name, self["Bocos"])
		assert(len(match) <= 1)
		if match:
			return match[0]

class Patch(Node):

	def __init__(self, pathRoot, dic):
		self.update(dic)
		lo, hi = self.MinMax()
		self.PointRange = PointRange2D(lo, hi)

		Node.__init__(self, "%s/%s" % (pathRoot, self["Name"]))

	def MinMax(self):

		r = np.array(self["Range"]) - 1
		lo = np.minimum(r[0:3], r[3:6])
		hi = np.maximum(r[0:3], r[3:6])
		return lo, hi

	def XYZ(self, xyz):

		lo, hi = self.MinMax()
		return xyz[:, lo[0]:hi[0]+1, lo[1]:hi[1]+1, lo[2]:hi[2]+1]

	def Index3DFrom2D(self, ij):

		ijk = np.zeros(3, int)
		ijk[self.PointRange.Idx] = ij
		ijk[self.PointRange.Idx0] = self.PointRange.Min[self.PointRange.Idx0]
		return ijk

	def CanonicalRange(self):

		r = self["Range"]
		start = np.min([r[:3], r[3:]], 0)
		end   = np.max([r[:3], r[3:]], 0)
		return [start[0], start[1], start[2], end[0], end[1], end[2]]

class PatchBoco(Patch):

	def __init__(self, pathRoot, dic):
		Patch.__init__(self, pathRoot, dic)

class Patch1to1(Patch):

	def __init__(self, pathRoot, dic):
		Patch.__init__(self, pathRoot, dic)

	def IsPeriodic(self):
		rotAngle = np.array(self["RotationAngle"])
		trans = np.array(self["Translation"])
		return np.abs(rotAngle).sum() > 0.0 or np.abs(trans).sum() > 0.0

class PatchNonMatching(Patch):

	def __init__(self, pathRoot, dic):
		Patch.__init__(self, pathRoot, dic)

class CGNSFile(object):

	def __init__(self, filename, mode = "r"):

		self.fn = CGNS.Open(filename, mode)

		self.BaseName, cellDim, physDim = CGNS.BaseRead(self.fn, 1)
		mass, length, time, temp, angle = CGNS.UnitsRead()
		self.Units = {"Mass":mass, "Length":length, "Time":time, "Angle":angle}

	def ReadZones(self, B = 1):

		nzones = CGNS.NZones(self.fn, B)

		base = Base(self.BaseName)

		zones = Zones()
		for Z in range(1, nzones + 1):

			name, size = CGNS.ZoneRead(self.fn, B, Z)
			#zone = {"Zone":Z, "Name":name, "Size":size}
			zone = Zone(base, Z, name, size)
			CGNS.GoPath(self.fn, zone.Path)
			famName = CGNS.FamNameRead()
			zone["FamilyName"] = famName

			nbocos = CGNS.NBocos(self.fn, B, Z)

			bocos = []
			for iboco in range(1, nbocos + 1):
				bocoName, bocoType, r = CGNS.BocoInfo(self.fn, B, Z, iboco)
				patchBoco = PatchBoco("%s/ZoneBC" % zone.Path, {"Zone":zone, "Name":bocoName, "Type":bocoType, "Range":r})
				if bocoType == "FamilySpecified":
					CGNS.GoPath(self.fn, "%s" % patchBoco.Path)
					famName = CGNS.FamNameRead()
					patchBoco["FamilyName"] = famName
				bocos.append(patchBoco)
			zone["Bocos"] = bocos

			n1to1s = CGNS.N1to1(self.fn, B, Z)

			c1to1s = []
			for I in range(1, n1to1s + 1):
				connectName, donorName, selfRange, donorRange, transform = CGNS.C1to1Read(self.fn, B, Z, I)
				rotCenter, rotAngle, trans = CGNS.C1to1PeriodicRead(self.fn, B, Z, I)
				c1to1 = Patch1to1("%s/__NOTIMPLMENTED__" % zone.Path, {"Zone":zone, "Name":connectName, "DonorZoneName":donorName, "Range":selfRange, "DonorRange":donorRange, "Transform":transform, "RotationCenter":rotCenter, "RotationAngle":rotAngle, "Translation":trans})
				c1to1s.append(c1to1)
			zone["1to1s"] = c1to1s

			nconns = CGNS.NConns(self.fn, B, Z)
			conns = [] # general connectivity (NMB, FNMB, etc.)
			for I in range(1, nconns + 1):
				connectName, gridLocation, connectType, ptsetType, donorName, donorZoneType, donorPtsetType, donorDataType, pnts, donorData = CGNS.ConnRead(self.fn, B, Z, I)
				assert(ptsetType == "PointRange")
				# Numeca-specific data
				# first check if this is an NMB connectivity. (it could also be an inter-stage mixing plane, which we ignore)
				CGNS.GoPath(self.fn, "%s/ZoneGridConnectivity/%s" % (zone.Path, connectName))
				ndescriptors = CGNS.NDescriptors()
				if ndescriptors == 0:
					continue
				descriptors = {}
				for D in range(1, ndescriptors + 1):
					name, text = CGNS.DescriptorRead(D)
					descriptors[name] = text
				if not descriptors.has_key("ConnectionId"):
					continue
				arrays = {}
				CGNS.GoPath(self.fn, "%s/ZoneGridConnectivity/%s/DonorPatch" % (zone.Path, connectName))
				narrays = CGNS.NArrays()
				for A in range(1, narrays + 1):
					arrayName, dataType, dataDim, dimVec = CGNS.ArrayInfo(A)
					arrayData = CGNS.ArrayReadInteger(A)
					arrays[arrayName] = arrayData
				cnmbDict = {"Zone":zone, "Name":connectName, "DonorZoneName":donorName, "Range":pnts, "donorData":donorData, "DonorPatch":arrays}
				cnmbDict.update(descriptors)
				conn =  PatchNonMatching("%s/ZoneGridConnectivity" % zone.Path, cnmbDict)
				conns.append(conn)
			zone["CNMBs"] = conns

			zones.append(zone)

		self.Zones = zones

		return zones

	def ReadZoneCoord(self, B, Z, xyzAxis = 0):

		xyz = CGNS.CoordRead2(self.fn, B, Z, xyzAxis)
		return xyz

	def ReadPatchCoord(self, B, Z, patch):

		xyz = self.ReadZoneCoord(B, Z)
		r = patch.CanonicalRange()
		xyz2 = np.array(xyz[:, r[0] - 1:r[3], r[1] - 1:r[4], r[2] - 1:r[5]])
		del xyz
		return xyz2

	def ReadSolution(self, B, Z):

		S = 1
		nFields = CGNS.NFields(self.fn, B, Z, S)
		Qs = {}
		for F in range(1, nFields + 1):
			dataType, fieldName = CGNS.FieldInfo(self.fn, B, Z, S, F)
			Q = CGNS.FieldRead2(self.fn, B, Z, S, fieldName)
			Qs[fieldName] = Q

		return Qs

	def ReadRigidBodyMotion(self, B, Z):

		nrigids = CGNS.NRigidMotions(self.fn, B, Z)
		if nrigids == 0:
			return None

		assert(nrigids == 1)
		R = 1
		rmName, rmType, origin, omega = CGNS.RigidMotionRead(self.fn, B, Z, R)

		return origin, omega

	def Close(self):

		CGNS.Close(self.fn)

	def ReadFamilies(self, B = 1):

		families = {}

		nFamilies = CGNS.NFamilies(self.fn, B);
		for Fam in range(1, nFamilies + 1):
			familyName, nFamBC, nGeo = CGNS.FamilyRead(self.fn, B, Fam)
			assert(not families.has_key(familyName))
			assert(nFamBC <= 1)
			family = {}
			# Read family bc
			if nFamBC == 1:
				famBCName, bocoType = CGNS.FamilyBCRead(self.fn, B, Fam, 1)
				family["BC"] = {"Name":famBCName, "Type":bocoType}
			# Read other family properties
			CGNS.GoPath(self.fn, "/Base#1/%s" % familyName)
			nUserData = CGNS.NUserData()
			if nUserData > 0:
				for i in range(1, nUserData + 1):
					userDataName = CGNS.UserDataRead(i)
					family[userDataName] = None
			families[familyName] = family

		self.Families = families

		return families

	def AllBocos(self):

		allBocos = []
		for zone in self.Zones:
			allBocos.extend(zone["Bocos"])
		return allBocos

	def MatchBocosByName(self, pattern):

		matchedBocos = filter(lambda boco: fnmatch.fnmatch(boco["Name"], pattern), self.AllBocos())
		return matchedBocos

	def MatchBocosByFamilyName(self, familyName):

		matchedBocos = filter(lambda boco: fnmatch.fnmatch(boco["FamilyName"], familyName), self.AllBocos())
		return matchedBocos

	def MatchZonesByName(self, pattern):

		matchedZones = filter(lambda zone: fnmatch.fnmatch(zone["Name"], pattern), self.Zones)
		return matchedZones

	def MatchZonesByFamilyName(self, familyName):

		matchedZones = filter(lambda zone: fnmatch.fnmatch(zone["FamilyName"], familyName), self.Zones)
		return matchedZones

	def FindZoneByName(self, name):

		return self.Zones.FindZoneByName(name)

	def FindZoneAndBocoByName(self, zoneName, bocoName):

		zone = self.FindZoneByName(zoneName)
		boco = zone.FindBocoByName(bocoName)
		return zone, boco

def ComputeCellCenters(XYZ):

	XYZC = np.zeros(np.array(XYZ.shape) - np.array((1, 1, 1, 0)))
	print XYZ.shape
	print XYZC.shape

	XYZC[:, :, :, :] = (XYZ[:-1, :-1, :-1, :] + XYZ[1:, :-1, :-1, :] + XYZ[:-1, 1:, :-1, :] + XYZ[1:, 1:, :-1, :]
		+ XYZ[:-1, :-1, 1:, :] + XYZ[1:, :-1, 1:, :] + XYZ[:-1, 1:, 1:, :] + XYZ[1:, 1:, 1:, :]) / 8.0
	return XYZC

def ComputeROmegas(XYZ, rot, Qs):

	print "ComputeROmegas"
	XYZC = ComputeCellCenters(XYZ)
	if rot:
		origin, omega = rot
		omega = np.array(omega)
		no = omega / m.sqrt(np.dot(omega, omega))
	else:
		origin = np.zeros(3)
		omega = np.zeros(3)
		no = np.zeros(3)

	p = XYZC
	pa = np.dot(p, no)
	pa = pa[:, :, :, np.newaxis] * no
	pr = XYZC - pa
	rsq = (pr * pr).sum(axis = -1)
	r = np.sqrt(rsq)
	Qs["R"] = r

	# ROmega = mag(Ve), so we don't need to store this quantity.
	#Qs["ROmega"] = r * np.sqrt(np.dot(omega, omega))

	# entrainment velocity vector
	ve = np.cross(omega, p)
	Qs["Ve"] = ve

	del XYZC


