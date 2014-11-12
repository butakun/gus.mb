#!/usr/bin/env python

import CGNSFile
import CGNS
import numpy as n

class BoundingBox(object):

	def __init__(self, pmin, pmax):
		self.Min = pmin
		self.Max = pmax

	def __add__(self, bb2):
		pmin = n.minimum(self.Min, bb2.Min)
		pmax = n.maximum(self.Max, bb2.Max)
		return BoundingBox(pmin, pmax)

	def __str__(self):
		return str(self.Min) + " - " + str(self.Max)

class Transform(object):

	def __init__(self, T, Begin1, End1, Begin2):
		self.T = T
		self.TT = self.__TT(self.T)
		self.Begin1 = Begin1
		self.End1 = End1
		self.Begin2 = Begin2
		self.End2 = self.Index2(End1)

	def Index2(self, Index1):
		return n.dot(self.TT, (Index1 - self.Begin1)) + self.Begin2

	def __TT(self, T):
		TT = n.zeros((3, 3), int)
		for i in range(3):
			ja = abs(T[i])
			j = ja - 1
			TT[j, i] = n.sign(T[i])
		return TT

class Patch(object):

	def __init__(self, zone, boco):
		self.XYZ = zone["XYZ"]
		self.Boco = boco
		self.Name = boco["Name"]
		self.Zone = zone["Zone"]
		r = self.Boco["Range"]
		i0 = n.array([r[0]-1, r[1]-1, r[2]-1])
		i1 = n.array([r[3], r[4], r[5]])
		self.Range = i1 - i0
		self.Begin1 = i0
		self.End1 = i1 - n.array([1, 1, 1])
		if self.Range[0] == 1:
			self.Delta1 = n.array([0, 1, 0])
			self.Delta2 = n.array([0, 0, 1])
			if i0[0] == 0:
				self.Normal = n.array([1, 0, 0])
			else:
				self.Normal = n.array([-1, 0, 0])
		elif self.Range[1] == 1:
			self.Delta1 = n.array([1, 0, 0])
			self.Delta2 = n.array([0, 0, 1])
			if i0[1] == 0:
				self.Normal = n.array([0, 1, 0])
			else:
				self.Normal = n.array([0, -1, 0])
		elif self.Range[2] == 1:
			self.Delta1 = n.array([1, 0, 0])
			self.Delta2 = n.array([0, 1, 0])
			if i0[2] == 0:
				self.Normal = n.array([0, 0, 1])
			else:
				self.Normal = n.array([0, 0, -1])
		self.Rotation = None
		self.Translation = None

	def __ApplyTransform(self, xyz):

		if self.Rotation and self.Translation:
			raise ValueError

		if self.Translation:
			return self.__ApplyTranslation(xyz)
		elif self.Rotation:
			return self.__ApplyRotation(xyz)
		else:
			return n.array(xyz)

	def __ApplyTranslation(self, xyz):
		if not self.Translation:
			return n.array(xyz)

		xyz2 = n.array(xyz)
		t = n.array(self.Translation)
		xyz2[0, :, :, :] += t[0]
		xyz2[1, :, :, :] += t[1]
		xyz2[2, :, :, :] += t[2]

		return xyz2

	def __ApplyRotation(self, xyz):
		if not self.Rotation:
			return n.array(xyz)

		axis, degrees = self.Rotation
		theta = n.radians(degrees)
		c = n.cos(theta)
		s = n.sin(theta)

		xyz2 = n.array(xyz)
		if axis == "X":
			raise ValueError
		elif axis == "Y":
			raise ValueError
		elif axis == "Z":
			xyz2[0, :, :, :] = c * xyz[0, :, :, :] - s * xyz[1, :, :, :]
			xyz2[1, :, :, :] = s * xyz[0, :, :, :] + c * xyz[1, :, :, :]
		return xyz2

	def BoundingBox(self):
		xyz = self.PatchCoord()
		xmin = n.min(xyz[0, :, :, :])
		xmax = n.max(xyz[0, :, :, :])
		ymin = n.min(xyz[1, :, :, :])
		ymax = n.max(xyz[1, :, :, :])
		zmin = n.min(xyz[2, :, :, :])
		zmax = n.max(xyz[2, :, :, :])
		return BoundingBox(n.array([xmin, ymin, zmin]), n.array([xmax, ymax, zmax]))

	def PatchCoord(self):
		r = self.Boco["Range"]
		xyz = self.XYZ[:, r[0]-1:r[3], r[1]-1:r[4], r[2]-1:r[5]]
		return self.__ApplyTransform(xyz)

	def PatchIndexToGlobalIndex(self, ip):
		return n.array(ip) + self.Begin1

def FindMinimalCellSize(xyz, dsmin):

	dsI = xyz[:, 1:, :, :] - xyz[:, :-1, :, :]
	dsI = dsI * dsI
	dsI = n.add.reduce(dsI, 0)
	dsI = n.sqrt(n.min(dsI))

	dsJ = xyz[:, :, 1:, :] - xyz[:, :, :-1, :]
	dsJ = dsJ * dsJ
	dsJ = n.add.reduce(dsJ, 0)
	dsJ = n.sqrt(n.min(dsJ))

	dsK = xyz[:, :, :, 1:] - xyz[:, :, :, :-1]
	dsK = dsK * dsK
	dsK = n.add.reduce(dsK, 0)
	dsK = n.sqrt(n.min(dsK))

	dsmin_= dsmin
	dsmin_ = min(dsmin_, dsI)
	dsmin_ = min(dsmin_, dsJ)
	dsmin_ = min(dsmin_, dsK)

	return dsmin_

def PointContainedInPatch(p, patch, tol):

	tolsq = tol * tol
	xyzD = patch.PatchCoord()
	d = xyzD.reshape(3, -1) - p.reshape(3, 1)
	d = d.reshape(xyzD.shape)
	d = d * d
	d = n.add.reduce(d, 0)
	dsqmin = n.min(d)
	contained = dsqmin < tolsq

	minloc = None
	for i in range(d.shape[0]):
		for j in range(d.shape[1]):
			for k in range(d.shape[2]):
				if d[i, j, k] == dsqmin:
					minloc = (i, j, k)
					break

	if contained:
		return minloc
	else:
		return None

def CheckTransform(delta11, delta12, delta21, delta22, patch1, patch2):
	"""
	Find a linear transform between corner indices of Patch1 (corner1)
	and those of Patch2 (corner2)
	"""

	T = [0, 0, 0]
	for d1, d2 in [[delta11, delta21], [delta12, delta22]]:
		dd1 = n.abs(d1)
		if dd1[0] == 1:
			i = 1
		elif dd1[1] == 1:
			i = 2
		elif dd1[2] == 1:
			i = 3
		else:
			raise ValueError
		dd2 = n.abs(d2)
		if dd2[0] == 1:
			j = 1
		elif dd2[1] == 1:
			j = 2
		elif dd2[2] == 1:
			j = 3
		else:
			raise ValueError

		sign1 = n.sign(d1[i - 1])
		sign2 = n.sign(d2[j - 1])
		sign = sign1 * sign2
		T[i - 1] = sign * j

	ij = [0, 0]
	for l, normal in [[0, patch1.Normal], [1, patch2.Normal]]:
		nn = n.abs(normal)
		if nn[0] == 1:
			ij[l] = n.sign(normal[0]) * 1
		elif nn[1] == 1:
			ij[l] = n.sign(normal[1]) * 2
		elif nn[2] == 1:
			ij[l] = n.sign(normal[2]) * 3
	i, j = ij
	T[abs(i) - 1] = -1 * n.sign(i) * n.sign(j) * abs(j)

	return T

def FindDonorPatch(patch, donorPatches, tol):
	"""
	Find the donor patch from the candidate patches (donorPatches)
	"""

	xyz = patch.PatchCoord()
	delta11 = n.array(patch.Delta1)
	delta12 = n.array(patch.Delta2)
	diag = delta11 + delta12
	p11 = xyz[:, 0, 0, 0]
	p12 = xyz[:, diag[0], diag[1], diag[2]]
	p21 = xyz[:, -1, -1, -1]
	p22 = xyz[:, -2 * diag[0], -2 * diag[1], -2 * diag[2]]

	donors = []
	for donorPatch in donorPatches:
		flag11 = PointContainedInPatch(p11, donorPatch, tol)
		flag12 = PointContainedInPatch(p12, donorPatch, tol)
		flag21 = PointContainedInPatch(p21, donorPatch, tol)
		flag22 = PointContainedInPatch(p22, donorPatch, tol)
		if flag11 and flag12 and flag21 and flag22:
			print "Patch %s overlaps Donor Patch %s" % (patch.Name, donorPatch.Name)
			#print "  ", flag11, flag12, flag21, flag22

			pp1 = xyz[:, delta11[0], delta11[1], delta11[2]]
			pp2 = xyz[:, delta12[0], delta12[1], delta12[2]]
			i0 = n.array(flag11)
			i1 = n.array(PointContainedInPatch(pp1, donorPatch, tol))
			i2 = n.array(PointContainedInPatch(pp2, donorPatch, tol))
			delta21 = i1 - i0
			delta22 = i2 - i0
			#print "i0, i1, i2", i0, i1, i2

			T = CheckTransform(delta11, delta12, delta21, delta22, patch, donorPatch)
			print "  Transform = ", T

			transform = Transform(T, patch.Begin1, patch.End1, donorPatch.PatchIndexToGlobalIndex(i0))

			#donors.append([donorPatch, T, i0])
			donors.append([donorPatch, transform])

	assert(len(donors) <= 1)
	if donors:
		return donors[0]
	else:
		return None, None

def WritePLOT3D(f, xyz):

	idim, jdim, kdim = xyz.shape[1:]
	print >> f, idim, jdim, kdim
	for l in range(3):
		for k in range(kdim):
			for j in range(jdim):
				for i in range(idim):
					print >> f, xyz[l, i, j, k]

def IdentifyConnectivities(zones, patchNames1, patchNames2, rotation, translation, tol):

	patches1 = []
	for patchName in patchNames1:
		zone, boco = zones.FindBocoByName(patchName)
		patch = Patch(zone, boco)
		patches1.append(patch)

	patches2 = []
	for patchName in patchNames2:
		zone, boco = zones.FindBocoByName(patchName)
		patch = Patch(zone, boco)
		patch.Rotation = rotation
		patch.Translation = translation
		patches2.append(patch)

	pairs = []
	for patch1 in patches1:
		donor, transform = FindDonorPatch(patch1, patches2, tol)
		pairs.append([patch1, donor, transform])

	for patch2 in patches2:
		donor, transform = FindDonorPatch(patch2, patches1, tol)
		pairs.append([patch2, donor, transform])

	patches = {}
	for patch, donor, transform in pairs:
		if donor:
			if not patches.has_key(donor.Name):
				patches[donor.Name] = [patch.Name]
			else:
				patches[donor.Name].append(patch.Name)

	print "Identified the following periodic connectivities:"
	for patch, donor, transform in pairs:
		if donor:
			Begin1 = transform.Begin1
			End1 = transform.End1
			Begin2 = transform.Begin2
			End2 = transform.End2
			print "  %s (%d) <-> %s (%d)" % (patch.Name, patch.Zone, donor.Name, donor.Zone), " T = ", transform.T, Begin1, End1, Begin2, End2
		else:
			print "  ", patch.Name, " <-> None found"

	return pairs

def WriteConnectivities(f, pairs, zones):

	for patch, donor, transform in pairs:
		begin1 = transform.Begin1 + n.array([1, 1, 1])
		end1 = transform.End1 + n.array([1, 1, 1])
		begin2 = transform.Begin2 + n.array([1, 1, 1])
		end2 = transform.End2 + n.array([1, 1, 1])
		selfRange = list(begin1)
		selfRange.extend(end1)
		donorRange = list(begin2)
		donorRange.extend(end2)
		donorZoneName = zones.FindZoneByID(donor.Zone)["Name"]
		print "writing ranges = ", selfRange, donorRange, donorZoneName
		CGNS.C1to1Write(f.fn, 1, patch.Zone, patch.Name, donorZoneName, selfRange, donorRange, transform.T)

def Main(cgnsMeshFileName, action):

	f = CGNSFile.CGNSFile(cgnsMeshFileName, "w")
	zones = f.ReadZones()

	dsmin = 1.0e30
	print dsmin
	B = 1
	for zone in zones:
		print "Reading Zone %d, Name \"%s\"" % (zone["Zone"], zone["Name"])
		xyz = f.ReadZoneCoord(B, zone["Zone"])

		dsmin = FindMinimalCellSize(xyz, dsmin)
		print "dsmin = %e" % dsmin

		for boco in zone["Bocos"]:
			print "Boco \"%s\", Type %s, Range %s" % (boco["Name"], boco["Type"], str(boco["Range"]))

		zone["XYZ"] = xyz

	tol = dsmin * 0.49
	print "tol = ", tol

	"""
	patchNames1 = ["Unspecified", "Unspecified4", "Unspecified5", "Unspecified8", "Unspecified13"]
	patchNames2 = ["Unspecified2", "Unspecified3", "Unspecified9", "Unspecified14"]
	rotation = ["Z", -360.0 / 36.0]
	translation = None
	"""

	patchNames1 = ["Periodic1", "Periodic12", "Periodic13", "Periodic14", "Periodic15", "Periodic16"]
	patchNames2 = ["Periodic2", "Periodic22", "Periodic23", "Periodic24", "Periodic25", "Periodic26"]
	rotation = None
	translation = [0.0, -175.0, 0.0]

	pairs = IdentifyConnectivities(zones, patchNames1, patchNames2, rotation, translation, tol)

	if action == "write":
		WriteConnectivities(f, pairs, zones)

if __name__ == "__main__":
	import sys
	Main(sys.argv[1], sys.argv[2])

