# $Id: CGNSFile.py 314 2013-12-24 15:47:30Z kato $

import numpy as np
import math as m
import CGNS
from PointRange import *

class Zones(list):

	def AllBocos(self):

		patches = []
		for zone in self:
			patches.extend(zone["Bocos"])
		return patches

	def FindBocoByName(self, name):

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

class CGNSFile(object):

	def __init__(self, filename, mode = "r"):

		self.fn = CGNS.Open(filename, mode)

		self.BaseName, cellDim, physDim = CGNS.BaseRead(self.fn, 1)

	def ReadZones(self, B = 1):

		nzones = CGNS.NZones(self.fn, B)

		base = Base(self.BaseName)

		zones = Zones()
		for Z in range(1, nzones + 1):

			name, size = CGNS.ZoneRead(self.fn, B, Z)
			#zone = {"Zone":Z, "Name":name, "Size":size}
			zone = Zone(base, Z, name, size)

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
			zones.append(zone)

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

		return families

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


