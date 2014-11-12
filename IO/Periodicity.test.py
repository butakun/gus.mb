# $Id: Periodicity.py -1   $

import CGNSFile, CGNS
import PeriodicityUtils
import numpy as np
import math

class Rotate(object):

	def __init__(self, axis, angle):
		"""
		axis: "X", "Y", or "Z"
		angle: in radians
		"""
		assert(axis == "X")
		self.Axis = axis
		self.Angle = angle

	def __call__(self, XYZ):

		c = math.cos(self.Angle)
		s = math.sin(self.Angle)

		if self.Axis == "X":
			# Rotation around X-axis.
			# [ 1  0  0 ]
			# [ 0  c -s ]
			# [ 0  s  c ]
			R = np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
		elif self.Axis == "Y":
			raise ValueError
		elif self.Axis == "Z":
			raise ValueError
		else:
			raise ValueError

		XYZr = np.zeros(XYZ.shape)
		for i in range(XYZ.shape[1]):
			for j in range(XYZ.shape[2]):
				XYZr[:, i, j] = np.dot(R, XYZ[:, i, j])
		return XYZr

	def Inverse(self):

		return Rotate(self.Axis, -self.Angle)

class Translate(object):

	def __init__(self, offset):
		self.Offset = np.array(offset)

	def __call__(self, XYZ):
		XYZt = XYZ + self.Offset[:, np.newaxis, np.newaxis]
		return XYZt

	def Inverse(self):

		return Translate(-self.Offset)

"""
def Rotate(origin, axis, angle, XYZ):

	c = math.cos(angle)
	s = math.sin(angle)

	# Rotation around X-axis.
	# [ 1  0  0 ]
	# [ 0  c -s ]
	# [ 0  s  c ]
	R = np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])

	XYZr = np.zeros(XYZ.shape)
	for i in range(XYZ.shape[1]):
		for j in range(XYZ.shape[2]):
			XYZr[:, i, j] = np.dot(R, XYZ[:, i, j])
	return XYZr

def ExtractPatches(f, patchDescs):

	XYZs = []
	for zone, desc in patchDescs:
		XYZ = f.ReadZoneCoord(zone)
		if desc == "IMIN":
			XYZs.append(XYZ[:, 0, :, :])
		elif desc == "IMAX":
			XYZs.append(XYZ[:, -1, :, :])
		elif desc == "JMIN":
			XYZs.append(XYZ[:, :, 0, :])
		elif desc == "JMAX":
			XYZs.append(XYZ[:, :, -1, :])
		elif desc == "KMIN":
			XYZs.append(XYZ[:, :, :, 0])
		elif desc == "KMAX":
			XYZs.append(XYZ[:, :, :, -1])
		else:
			raise ValueError
	return XYZs

def Main(filename):

	f = CGNSFile.CGNSFile(filename)

	nblades = 47
	PS = [[1, "KMIN"], [3, "KMIN"], [4, "IMAX"], [6, "JMIN"], [7, "JMIN"]]
	SS = [[1, "KMAX"], [2, "KMAX"], [6, "JMAX"], [7, "JMAX"]]
	PeriodicPS = ExtractPatches(f, PS)
	PeriodicSS = ExtractPatches(f, SS)
	AnglePS2SS = 2.0 * math.pi / float(nblades)

	pairs = [[PS, PeriodicPS, SS, PeriodicSS, AnglePS2SS], [SS, PeriodicSS, PS, PeriodicPS, -AnglePS2SS]]
	for selfDescs, selfPatches, donorDescs, donorPatches, angleSelfToDonor in pairs:
		for i, XYZ in enumerate(selfPatches):
			print selfDescs[i]
			XYZr = Rotate(None, None, angleSelfToDonor, XYZ)
			matches = PeriodicityUtils.MatchCells(XYZr, donorPatches)
			print "1-to-1 connectivity to "
			for irange, jrange, donorPatch, donorIRange, donorJRange in matches:
				print "  ", irange, jrange, " => ", donorDescs[donorPatch], donorIRange, donorJRange
"""

class PatchTriplets(object):

	def __init__(self, Z, lo, hi):
		self.Z = Z
		self.Lo = lo
		self.Hi = hi
		self.Name = None

	def __eq__(self, p):
		return self.Z == p.Z and np.all(self.Lo == p.Lo) and np.all(self.Hi == p.Hi)

	def TransformWithDonor(self, donor):
		"""
		Index2 = T * (Index1 - Begin1) + Begin2
		"""
		pass

	def RangeTriplet(self):
		"""
		Bias-One, suitable for CGNS function calls.
		"""
		return np.array([self.Lo[0], self.Lo[1], self.Lo[2], self.Hi[0], self.Hi[1], self.Hi[2]], int) + 1

def TransformTriplet(T):
	"""
	T = 3x3 index transform matrix, this function returns a shortened 3-component vector
	"""
	rows, cols = np.where(np.abs(T) == 1)
	t = T[rows, cols] * (cols + 1)
	return t

def CellRangeToNodeRange(irange, jrange):

	irangeNode = list(irange)
	jrangeNode = list(jrange)
	for range_ in [irangeNode, jrangeNode]:
		if range_[1] > range_[0]:
			range_[1] += 1
		elif range_[1] < range_[0]:
			range_[0] += 1
		else:
			# Does not suppor one-cell-width patch yet
			raise ValueError
	return irangeNode, jrangeNode

def ApplyPeriodicity(filename, PS, SS, transformPS2SS):

	transformSS2PS = transformPS2SS.Inverse()

	f = CGNSFile.CGNSFile(filename)

	zones = f.ReadZones()

	PSxyz = []
	SSxyz = []
	for patchXYZ, names in [[PSxyz, PS], [SSxyz, SS]]:
		for name in names:
			zone, patch = zones.FindBocoByName(name)
			xyzZone = f.ReadZoneCoord(zone["Zone"])
			xyz = patch.XYZ(xyzZone)
			patchXYZ.append([zone, patch, xyz])

	f.Close()

	pairs = []
	for selfPatches, donorPatches, transform in [[PSxyz, SSxyz, transformPS2SS], [SSxyz, PSxyz, transformSS2PS]]:
		for zoneSelf, patchSelf, xyzSelf in selfPatches:
			XYZ = xyzSelf.squeeze()
			XYZr = transform(XYZ)
			XYZsDonor = [ xyzD.squeeze() for ZD, patchD, xyzD in donorPatches ]
			matches = PeriodicityUtils.MatchCells(XYZr, XYZsDonor)
			print "Self: ", zoneSelf["Zone"]
			assert(len(matches) > 0)
			for irange, jrange, donorIndex, donorIRange, donorJRange, dI, dJ in matches:
				irangeNode, jrangeNode = CellRangeToNodeRange(irange, jrange)
				ZS = zoneSelf["Zone"]
				StartS = patchSelf.Index3DFrom2D([irangeNode[0], jrangeNode[0]])
				EndS   = patchSelf.Index3DFrom2D([irangeNode[1], jrangeNode[1]])

				donorIRangeNode, donorJRangeNode = CellRangeToNodeRange(donorIRange, donorJRange)
				zoneDonor, patchDonor, xyzDonor = donorPatches[donorIndex]
				ZD = zoneDonor["Zone"]
				StartD = patchDonor.Index3DFrom2D([donorIRangeNode[0], donorJRangeNode[0]])
				EndD = patchDonor.Index3DFrom2D([donorIRangeNode[1], donorJRangeNode[1]])

				T = np.zeros((3, 3), int)
				T[patchDonor.PointRange.Idx, patchSelf.PointRange.Idx[0]] = dI
				T[patchDonor.PointRange.Idx, patchSelf.PointRange.Idx[1]] = dJ
				KSelf  = StartS[patchSelf.PointRange.Idx0]
				KDonor = StartD[patchDonor.PointRange.Idx0]
				if (KSelf == 0 and KDonor == 0) or (KSelf > 0 and KDonor > 0):
					dK = -1
				else:
					dK = 1
				T[patchDonor.PointRange.Idx0, patchSelf.PointRange.Idx0] = dK
				StartDCheck = StartD + np.dot(T, StartS - StartS)
				EndDCheck = StartD + np.dot(T, EndS - StartS)

				print "  ", ZS, StartS, EndS, " => ", ZD, StartD, EndD, ", dI = ", dI, ", dJ = ", dJ
				print "     ", T[0, :], T[1, :], T[2, :], "(", TransformTriplet(T), ") : ", StartDCheck, EndDCheck
				pairs.append([PatchTriplets(ZS, StartS, EndS), PatchTriplets(ZD, StartD, EndD), TransformTriplet(T), transform])

	i = 0
	patches = {}
	for patch1, patch2, T, transform in pairs:
		if not patches.has_key(patch1):
			patches[patch1] = "Patch%02d" % i
			i += 1

	for patch in patches.keys():
		print patches[patch], patch.Z, patch.Lo, patch.Hi

	f = CGNSFile.CGNSFile(filename, "w")
	zones = f.ReadZones()

	# Delete patches
	patchesToDelete = list(PS)
	patchesToDelete.extend(list(SS))
	for patchName in patchesToDelete:
		zone, patch = zones.FindBocoByName(patchName)
		path = "/%s/%s/ZoneBC/" % (f.BaseName, zone["Name"])
		CGNS.GoPath(f.fn, path)
		CGNS.DeleteNode(patchName)

	# Add new 1-to-1 connectivities
	for patch1, patch2, T, transform in pairs:
		Z1 = patch1.Z
		Z2 = patch2.Z
		Range1 = list(patch1.RangeTriplet())
		Range2 = list(patch2.RangeTriplet())
		ZoneName2 = zones.FindZoneByID(Z2)["Name"]
		PatchName1 = patches[patch1]
		I = CGNS.C1to1Write(f.fn, 1, Z1, PatchName1, ZoneName2, Range1, Range2, list(T))
		if isinstance(transform, Rotate):
			center = [0.0, 0.0, 0.0]
			axis = [math.degrees(transform.Angle), 0.0, 0.0]
			translation = [0.0, 0.0, 0.0]
		elif isinstance(transform, Translate):
			center = [0.0, 0.0, 0.0]
			axis = [0.0, 0.0, 0.0]
			translation = list(transform.Offset)
		else:
			raise ValueError
		CGNS.C1to1PeriodicWrite(f.fn, 1, Z1, I, center, axis, translation)

def Main2(filename):

	PS = ["1_1", "3_1", "4_1", "6_1", "7_1"]
	SS = ["1_2", "2_1", "6_2", "7_2"]

	nblades = 47
	AnglePS2SS = 2.0 * math.pi / float(nblades)
	transformPS2SS = Rotate(axis = "X", angle = AnglePS2SS)

	ApplyPeriodicity(filename, PS, SS, transformPS2SS)

def Main3(filename):

	config = {}
	execfile(filename, globals(), config)

	Patches1 = config["Patches1"]
	Patches2 = config["Patches2"]
	Transform1to2 = config["Transform1to2"]
	CGNSFileName = config["CGNSFileName"]
	ApplyPeriodicity(CGNSFileName, Patches1, Patches2, Transform1to2)

if __name__ == "__main__":
	import sys
	#Main(sys.argv[1])
	#Main2(sys.argv[1])
	Main3(sys.argv[1])

