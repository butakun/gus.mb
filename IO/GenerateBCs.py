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
# $Id: GenerateBCs.py 316 2014-02-04 01:10:52Z kato $

import CGNSFile, LoadBalance
import numpy as np
import fnmatch

class C1to1(object):
	def __init__(self, zone, rang, donorZone, donorRange):
		self.Zone = zone
		self.Range = rang
		self.DonorZone = donorZone
		self.DonorRange = donorRange
		self.UniqueID = self.__ComputeUniqueID()

	def __ComputeUniqueID(self):
		ID = [0] * 14
		if (self.Zone < self.DonorZone):
			ID[0] = self.Zone
			ID[1:7] = CanonizeRange(self.Range)
			ID[7] = self.DonorZone
			ID[8:14] = CanonizeRange(self.DonorRange)
		elif (self.Zone > self.DonorZone):
			ID[0] = self.DonorZone
			ID[1:7] = CanonizeRange(self.DonorRange)
			ID[7] = self.Zone
			ID[8:14] = CanonizeRange(self.Range)
		elif (self.Zone == self.DonorZone):
			rZ = CanonizeRange(self.Range)
			rD = CanonizeRange(self.DonorRange)
			if (rZ[0] < rD[0]):
				r1 = rZ
				r2 = rD
			elif (rZ[0] > rD[0]):
				r1 = rD
				r2 = rZ
			elif (rZ[1] < rD[1]):
				r1 = rZ
				r2 = rD
			elif (rZ[1] > rD[1]):
				r1 = rD
				r2 = rZ
			elif (rZ[2] < rD[2]):
				r1 = rZ
				r2 = rD
			elif (rZ[2] > rD[2]):
				r1 = rD
				r2 = rZ
			else:
				raise ValueError # this is impossible
			ID[0] = self.Zone
			ID[1:7] = r1
			ID[7] = self.Zone
			ID[8:14] = r2
		return np.array(ID)

	def __hash__(self):
		return hash(np.sum(self.UniqueID))

	def __eq__(self, c):
		return np.all(self.UniqueID == c.UniqueID)

class CNMB(C1to1):
	def __init__(self, zone, rang, donorZone, donorRange):
		C1to1.__init__(self, zone, rang, donorZone, donorRange)

class Family(object):
	def __init__(self, name):
		self.FamilyName = name

def CanonizeRange(r):
	start = np.min([r[:3], r[3:]], 0)
	end   = np.max([r[:3], r[3:]], 0)
	return [start[0], start[1], start[2], end[0], end[1], end[2]]

def RangeString(r):

	rr = map(lambda i: i - 1, r) # FIXME: stupid. isn't it better to do "np.array(r) - np.one(len(r))"?
	return "%d %d %d %d %d %d" % (rr[0], rr[1], rr[2], rr[3], rr[4], rr[5])

def TripletString(t):
	return "%f %f %f" % (t[0], t[1], t[2])

def PatchNormalAndPos(xyz, axis):

	pk = ["X", "Y", "Z"].index(axis)
	pij = [0, 1, 2]
	del pij[pk]
	pi, pj = pij

	# Reshape the 3-d array of coordinates into a planar (2-d) array
	s = list(xyz.shape)
	k = s.index(1)
	del s[k]
	xyz2d = xyz.reshape(s)

	p1 = xyz2d[:, 0, 0]
	p2 = xyz2d[:, -1, 0]
	p3 = xyz2d[:, 0, -1]
	p4 = xyz2d[:, -1, -1]
	r1 = p1[pi] * p1[pi] + p1[pj] * p1[pj]
	r2 = p2[pi] * p2[pi] + p2[pj] * p2[pj]
	r3 = p3[pi] * p3[pi] + p3[pj] * p3[pj]
	r4 = p4[pi] * p4[pi] + p4[pj] * p4[pj]

	r = min(r1, r2, r3, r4)

#***DEPRECATED***
#def KOmega(tke, mut0, RhoRef):
#
#	omega = RhoRef * tke / mut0
##	return RhoRef * tke, RhoRef * omega

def TurbulenceSpecString(spec):

	if spec[0] == "NoTurbulence":
		return spec[0]
	elif spec[0] == "IntensityViscosityRatio":
		I, ViscRatio, Rho0, U0, T0 = spec[1:]
		return "%s %e %e %e %e %e" % (spec[0], I, ViscRatio, Rho0, U0, T0)
	else:
		raise ValueError

def ExpandWildcards(pattern, candidates):

	if isinstance(pattern, list):
		allMatches = []
		for pat in pattern:
			matches = filter(lambda a: fnmatch.fnmatch(a, pat), candidates)
			if len(matches) == 0:
				print "no match found for pattern %s" % pat
				raise ValueError
			allMatches.extend(matches)
		return allMatches
	else:
		return filter(lambda a: fnmatch.fnmatch(a, pattern), candidates)

def ExpandPatchNames(names, cgnsfile):

	matchedBocos = []
	for name in names:
		if isinstance(name, str):
			matched = cgnsfile.MatchBocosByName(name)
			if len(matched) > 0:
				matchedBocos.extend(matched)
		elif isinstance(name, Family):
			matched = cgnsfile.MatchBocosByFamilyName(name.FamilyName)
			if len(matched) > 0:
				matchedBocos.extend(matched)
		else:
			assert(false)
	return matchedBocos

def ExpandZones(names, cgnsfile):

	matchedZones = []
	for name in names:
		if isinstance(name, str):
			matched = cgnsfile.MatchZonesByName(name)
			if len(matched) > 0:
				matchedZones.extend(matched)
		elif isinstance(name, Family):
			matched = cgnsfile.MatchZonesByFamilyName(name.FamilyName)
			if len(matched) > 0:
				matchedZones.extend(matched)
		else:
			assert(false)
	return matchedZones

def Main(filename, commsize):

	config = {}
	execfile(filename, globals(), config)

	ref = config["Reference"] # Gamma, RhoRef, TRef, RGAS

	mesh = CGNSFile.CGNSFile(config["MeshFileName"])
	zones = mesh.ReadZones()
	families = mesh.ReadFamilies()

	# zone attributes
	zoneNames = [ zone["Name"] for zone in zones ]
	B = 1

	ZoneAttrs = config["ZoneAttributes"]
	for zoneAttr in ZoneAttrs:
		#zoneNamesAll = ExpandWildcards(zoneAttr["Zones"], zoneNames)
		zonesAll = ExpandZones(zoneAttr["Zones"], mesh)
		attrs = filter(lambda pair: pair[0] != "Zones", zoneAttr.items())
		for zone in zonesAll:
			#zone["RigidBodyMotion"] = zoneAttr["RigidBodyMotion"]
			for key, value in attrs:
				assert(key not in zone)
				zone[key] = value

	# fill in default attributes (so that the user doesn't have to explicitly set them)
	defaultAttrs = {"RigidBodyMotion":{"Type":"NONE", "Data":[0.0, 0.0, 0.0], "InitialOffset":[0.0, 0.0, 0.0]}, "Periodicity":{"Type":"NONE", "Data":[0.0, 0.0, 00]}}
	for zone in zones:
		for key, value in defaultAttrs.items():
			if key not in zone:
				zone[key] = value

	# extract boco patch names
	allBocos = mesh.AllBocos()
	bocoNames = [ patch["Name"] for patch in allBocos ]

	BCs = config["BCs"]
	for bc in BCs:
		patchNames = bc["Patches"]
		#patchNamesAll = ExpandWildcards(patchNames, bocoNames)
		#assert(len(patchNamesAll) > 0)
		bocos = ExpandPatchNames(patchNames, mesh)
		assert(len(bocos) > 0)
		bcType = bc["Type"]
		if bc.has_key("Data"):
			bcData = bc["Data"]
		else:
			bcData = None
		#for patchName in patchNamesAll:
		#	zone, boco = zones.FindBocoByName(patchName)
		for boco in bocos:
			zone = boco["Zone"]
			bocoInfo = [boco, bcType, bcData]
			if zone.has_key("BCs"):
				zone["BCs"].append(bocoInfo)
			else:
				zone["BCs"] = [bocoInfo]

	# assign unique ids to all bocos (FIXME: should assign to connectivities as well?)
	for i, boco in enumerate(allBocos):
		boco["UniqueID"] = i + 1;

	if config.has_key("Interfaces"):
		Interfaces = config["Interfaces"]
	else:
		Interfaces = []
	for zone in zones:
		zone["Interfaces"] = {}
	for intf in Interfaces:
		name, patches1, patches2 = intf
		for patchesSelf, patchesDonor in [[patches1, patches2], [patches2, patches1]]:
			#patchNamesSelf  = ExpandWildcards(patchesSelf, bocoNames)
			#patchNamesDonor = ExpandWildcards(patchesDonor, bocoNames)
			bocosSelf  = ExpandPatchNames(patchesSelf, mesh)
			bocosDonor = ExpandPatchNames(patchesDonor, mesh)
			for boco in bocosSelf:
				zone = boco["Zone"]
				if zone["Interfaces"].has_key(name):
					zone["Interfaces"][name][0].append(boco)
				else:
					zone["Interfaces"][name] = [[boco], bocosDonor]

	Solver = config["Solver"]
	TimeAdv = config["TimeAdvancement"]
	CFL = config["CFL"]
	Iterations = config["Iterations"]

	Execution = config["Execution"]
	if Execution["Type"] == "Parallel":
		if Execution["Distribution"] == "Simple":
			for iz in range(len(zones)):
				zone = zones[iz]
				zone["Rank"] = iz
		elif Execution["Distribution"] == "Automatic":
			if len(zones) < commsize:
				print "Number of zones < world communicator size"
				raise ValueError
			bins, bin_sizes = LoadBalance.Balance(zones, commsize)
			for rank, abin in enumerate(bins):
				for z in abin:
					zone = zones.FindZoneByID(z)
					zone["Rank"] = rank
		elif Execution["Distribution"] == "Custom":
			for irank, Zs in enumerate(Execution["ZoneGroups"]):
				for z in Zs:
					zone = zones.FindZoneByID(z)
					zone["Rank"] = irank
		else:
			raise ValueError
	else:
		for zone in zones:
			zone["Rank"] = 0

	initFlow = config["Initialization"]["Flow"]
	if initFlow["Type"] == "ISENTROPICEXPANSION":
		direction = initFlow["Direction"]
		pSpecs = initFlow["StaticPressureSpecs"]
		ppSpecs = [] # [pos, p]
		for zoneName, bocoName, pres in pSpecs:
			zone, boco = mesh.FindZoneAndBocoByName(zoneName, bocoName)
			xyz = mesh.ReadPatchCoord(B, zone["Zone"], boco)
			xyzIndex = ["X", "Y", "Z"].index(direction)
			pos = xyz[xyzIndex, 0, 0, 0] * config["MeshScale"]
			ppSpecs.append([pos, pres])
		initFlow["PositionAndStaticPressure"] = ppSpecs
	elif initFlow["Type"] == "CENTRIFUGAL":
		raise ValueError # FIXME: NOT IMPLEMENTED YET. NOT EVEN SURE IF I SHOULD EVEN BOTHER
		axis = initFlow["Axis"]
		specs = initFlow["Specs"]
		arpSpecs = [] # [axial, radial, pressure]
		for patchName, pres in specs:
			zone, patch = zones.FindBocoByName(patchName)
			xyz = mesh.ReadPatchCoord(B, zone["Zone"], patch)
			PatchNormalAndPos(xyz)
			raise Error
			p = xyz[:, 0, 0, 0]
			if axis == "X":
				a = p[0]
				r = m.sqrt(p[1] * p[1] + p[2] * p[2])
			elif axis == "Y":
				a = p[1]
				r = m.sqrt(p[0] * p[0] + p[2] * p[2])
			elif axis == "Z":
				a = p[2]
				r = m.sqrt(p[0] * p[0] + p[1] * p[1])
			arpSpecs.append([a, r, pres])
		initFlow["AxialRadialPressure"] = arpSpecs

	print "# Mesh File"
	print config["MeshFileName"]
	print "# Scale"
	print config["MeshScale"]
	print "# Reference"
	print ref[0], ref[1], ref[2], ref[3]
	print "# Time advancement"
	if TimeAdv["Type"] == "Steady":
		print 0
		print "N/A"
		Iterations["Subiterations"] = 1
	elif TimeAdv["Type"] == "Unsteady":
		print 1
		print TimeAdv["DeltaTime"], TimeAdv["Time0"]
	else:
		raise ValueError
	print "# Solver"
	print 0, CFL["CFL0"], CFL["CFLMAX"], CFL["SER"], CFL["StartSERAt"], int(CFL["LocalTimeStepping"])
	print "# Convective Flux"
	print config["Flux"]["Type"], config["Flux"]["Order"], config["Flux"]["Limiter"]
	print "# # of iterations, start turbulence at, subiterations"
	print Iterations["MaxIterations"], Iterations["TurbulenceStartsAt"], Iterations["Subiterations"]
	print "# # of blocks"
	print len(zones)
	print "# Blocks (BlockID, Rank, Zone, ZoneIndices)"
	for zone in zones:
		B = zone["Zone"]
		Z = B
		size = zone["Size"]
		meshRangeStr = RangeString([1, 1, 1, size[0], size[1], size[2]])
		print "# Block %d" % B
		print "%d\t%d\t%d\t%s" % (B, zone["Rank"], Z, meshRangeStr) # FIXME: need to write zone's rigid body motion info
		print "# Rigid body motion [Type(NONE|ROTATIONAL|TRANSLATIONAL], VecX, VecY, VecZ, initial offset...]"
		rbm = zone["RigidBodyMotion"]
		offset = rbm["InitialOffset"]
		print "%s %f %f %f %f %f %f" % (rbm["Type"], rbm["Data"][0], rbm["Data"][1], rbm["Data"][2], offset[0], offset[1], offset[2])
		print "# Periodicity [Type(NONE|ROTATIONAL|TRANSLATIONAL], VecX, VecY, VecZ]"
		per = zone["Periodicity"]
		print "%s %f %f %f" % (per["Type"], per["Data"][0], per["Data"][1], per["Data"][2])
	print "# Patch Unique IDs"
	nbocos = 0
	for zone in zones:
		nbocos += len(zone["Bocos"])
	print nbocos
	for zone in zones:
		Z = zone["Zone"]
		for boco in zone["Bocos"]:
			r = RangeString(boco.CanonicalRange())
			print "%d %s %d \"%s\""  % (Z, r, boco["UniqueID"], boco.Path)

	print "# Patch Families"
	print len(families)
	for i, family in enumerate(families.items()):
		familyID = i + 1
		familyName, familyData = family
		familyData["ID"] = familyID
		familyBocoIDs = []
		for boco in allBocos:
			if boco.has_key("FamilyName") and boco["FamilyName"] == familyName:
				familyBocoIDs.append(boco["UniqueID"])
		if len(familyBocoIDs) == 0:
			print familyID, " # ", familyName, " (no bocos found)"
		else:
			print familyID, reduce(lambda s1, s2: s1 + " " + s2, map(str, familyBocoIDs)), " # ", familyName

	InterfaceDB = {}
	nextTag = 1

	for zone in zones:
		print "# Zone"
		print zone["Zone"], "\t", zone["Name"]
		"""
		print "# Rigid motion [Type(NONE|ROTATIONAL|TRANSLATIONAL], VecX, VecY, VecZ]"
		rbm = zone["RigidBodyMotion"]
		rbmStr = "%s %f %f %f" % (rbm["Type"], rbm["Data"][0], rbm["Data"][1], rbm["Data"][2])
		print rbmStr
		print "# Periodicity [0:not periodic, 1:rotational, 2:translational]"
		print 0, 0.0, 0.0, 0.0
		"""
		bcs = zone["BCs"]
		print "# BCs (name, family, patch ID, range, type, spec)"
		print len(bcs)
		for bc in bcs:
			boco, bocoType, bocoData = bc
			familyID = -1
			if boco.has_key("FamilyName"):
				familyID = families[boco["FamilyName"]]["ID"]
			dataStr = "-1"
			if bocoType == "InletTotal":
				dataStr = "%f %f %f %f %f" % tuple(bocoData[0:5])
				turbSpec = bocoData[5]
				dataStr += " %s" % TurbulenceSpecString(turbSpec)
			elif bocoType == "Fix":
				dataStr = "%f %f %f %f %f" % tuple(bocoData[0:5])
				turbSpec = bocoData[5]
				dataStr += " %s" % TurbulenceSpecString(turbSpec)
			elif bocoType == "ViscousWall":
				dataStr = "%s %f" % (bocoData[0], bocoData[1])
			elif bocoType == "OutletStaticPressure":
				dataStr = "%f %d" % (bocoData[0], bocoData[1])
			else:
				if bocoData:
					dataStr = ""
					for v in bocoData:
						dataStr += " %f" % v
			print "\"%s\"\t%d\t%d\t%s\t%s\t%s  " % (boco.Path, familyID, boco["UniqueID"], RangeString(boco["Range"]), bocoType, dataStr)

		c1to1s = zone["1to1s"]
		print "# 1-to-1 Connectivities (name, tag, donorZoneID, selfRange, donorRange, transform, periodicity)"
		print len(c1to1s)
		for c1to1 in c1to1s:
			name = c1to1["Name"]
			donorZoneName = c1to1["DonorZoneName"]
			donorZoneID = zones.FindZoneByName(donorZoneName)["Zone"]
			sr = c1to1["Range"]
			dr = c1to1["DonorRange"]
			tr = c1to1["Transform"]
			interf = C1to1(zone["Zone"], sr, donorZoneID, dr)
			if not InterfaceDB.has_key(interf):
				InterfaceDB[interf] = nextTag;
				nextTag += 1
			tag = InterfaceDB[interf]
			selfRange = RangeString(sr)
			donorRange = RangeString(dr)
			transform = "%d %d %d" % (tr[0], tr[1], tr[2])
			periodicity = "%s %s %s" % (TripletString(c1to1["RotationCenter"]), TripletString(c1to1["RotationAngle"]), TripletString(c1to1["Translation"]))
			print "%s %d %d %s %s %s %s" % (name, tag, donorZoneID, selfRange, donorRange, transform, periodicity)

		cnmbs = zone["CNMBs"]
		print "# Non-matching Connectivities (name, tag, donorZoneID, selfRange, donorRange, transform)"
		print len(cnmbs)
		for cnmb in cnmbs:
			name = cnmb["Name"]
			donorZoneName = cnmb["DonorZoneName"]
			donorZoneID = zones.FindZoneByName(donorZoneName)["Zone"]
			sr = cnmb["Range"]
			dr = cnmb["DonorPatch"]["DonorPointRange"]
			tr = cnmb["DonorPatch"]["TransformMatrix"]
			interf = CNMB(zone["Zone"], sr, donorZoneID, dr)
			if not InterfaceDB.has_key(interf):
				InterfaceDB[interf] = nextTag;
				nextTag += 1
			tag = InterfaceDB[interf]
			selfRange = RangeString(sr)
			donorRange = RangeString(dr)
			transform = "%d %d %d" % (tr[0], tr[1], tr[2])
			print "%s %d %d %s %s %s" % (name, tag, donorZoneID, selfRange, donorRange, transform)

		print "# Interfaces"
		print len(zone["Interfaces"])
		for ifname, bocos_self_donor in zone["Interfaces"].items():
			print "\"%s\"" % ifname
			#for names in [patches, donors]:
			for bocos in bocos_self_donor:
				print len(bocos)
				for boco in bocos:
					z = boco["Zone"]
					print "\"%s\" %d %s" % (name, z["Zone"], RangeString(boco["Range"]))

		print "# Flowfield Initialization"
		print initFlow["Type"]
		if initFlow["Type"] == "UNIFORM":
			data = initFlow["Data"]
			print 1
			print "%f %f %f %f %f" % (data[0], data[1], data[2], data[3], data[4])
		elif initFlow["Type"] == "READFROMFILE":
			print 1
			print initFlow["FileName"]
		elif initFlow["Type"] == "ISENTROPICEXPANSION":
			ppSpecs = initFlow["PositionAndStaticPressure"]
			print len(ppSpecs) + 1
			print "%s %f %f # Direction, TotalPressure, TotalTemperature" % (initFlow["Direction"], initFlow["TotalPressure"], initFlow["TotalTemperature"])
			for pos, pres in ppSpecs:
				print "%f %f" % (pos, pres)
		elif initFlow["Type"] == "CENTRIFUGAL":
			arpSpecs = initFlow["AxialRadialPressure"]
			print len(arpSpecs) + 1
			print "%s %f %f # Axis, TotalPressure, TotalTemperature" % (initFlow["Axis"], initFlow["TotalPressure"], initFlow["TotalTemperature"])
			for axial, radial, pressure in arpSpecs:
				print "%f %f %f" % (axial, radial, pressure)
		else:
			raise ValueError

		initTurb = config["Initialization"]["Turbulence"]
		print "# Turbulence Initialization"
		print initTurb["Type"]
		if initTurb["Type"] == "UNIFORM":
			#data = initTurb["Data"]
			#print "%f %f" % (data[0], data[1])  # UT0[0], UT0[1]
			turbSpec = initTurb["Data"]
			print TurbulenceSpecString(turbSpec)
		elif initTurb["Type"] == "READFROMFILE":
			print initTurb["FileName"]

if __name__ == "__main__":
	import sys
	Main(sys.argv[1], int(sys.argv[2]))

