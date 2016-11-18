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

import CGNSFile
import numpy as np
import math as m
import os, optparse, sys, re

def ComputeCellCenters(XYZ):

	XYZC = np.zeros(np.array(XYZ.shape) - np.array((1, 1, 1, 0)))
	print XYZ.shape
	print XYZC.shape

	XYZC[:, :, :, :] = (XYZ[:-1, :-1, :-1, :] + XYZ[1:, :-1, :-1, :] + XYZ[:-1, 1:, :-1, :] + XYZ[1:, 1:, :-1, :]
		+ XYZ[:-1, :-1, 1:, :] + XYZ[1:, :-1, 1:, :] + XYZ[:-1, 1:, 1:, :] + XYZ[1:, 1:, 1:, :]) / 8.0
	return XYZC

def ComputeROmegas(XYZC, rot, Qs):

	print "ComputeROmegas"
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

	print "ComputeROmegas done"

def ProcessMomentum(Qs):

	px = Qs["MomentumX"]
	py = Qs["MomentumY"]
	pz = Qs["MomentumZ"]

	s1 = px.shape
	p = np.zeros((s1[0], s1[1], s1[2], 3))
	p[:, :, :, 0] = px
	p[:, :, :, 1] = py
	p[:, :, :, 2] = pz

	Qs["Momentum"] = p
	Qs.pop("MomentumX")
	Qs.pop("MomentumY")
	Qs.pop("MomentumZ")

def Flatten(XYZ, axis):
	""" Flatten annular mesh into linear cascade mesh """

	assert(axis == "Z")

	R = np.sqrt(XYZ[:, :, :, 0] * XYZ[:, :, :, 0] + XYZ[:, :, :, 1] * XYZ[:, :, :, 1])
	theta = np.arctan2(XYZ[:, :, :, 1], XYZ[:, :, :, 0])
	theta = 0.5 * np.pi - theta

	XYZ[:, :, :, 0] = R * theta
	XYZ[:, :, :, 1] = R

def WriteVTS(o, XYZ, Qs):

	pshape = XYZ.shape[:-1]
	cshape = np.array(pshape) - np.array((1, 1, 1))

	print >>o, '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'
	print >>o, '<StructuredGrid WholeExtent="0 %d 0 %d 0 %d">' % (pshape[0] - 1, pshape[1] - 1, pshape[2] - 1)
	print >>o, '<Piece Extent="0 %d 0 %d 0 %d">' % (pshape[0] - 1, pshape[1] - 1, pshape[2] - 1)
	print >>o, '<Points>'
	print >>o, '<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">'
	for k in range(pshape[2]):
		for j in range(pshape[1]):
			for i in range(pshape[0]):
				print >>o, "%12.4e %12.4e %12.4e" % (XYZ[i, j, k, 0], XYZ[i, j, k, 1], XYZ[i, j, k, 2])
	print >>o, '</DataArray>'
	print >>o, '</Points>'
	print >>o, '<CellData>'
	for fieldName, Q in Qs.items():
		q = Q
		if len(Q.shape) < 4:
			s = Q.shape
			q = Q.reshape(s[0], s[1], s[2], 1)
		numComponents = q.shape[3]
		print >>o, '<DataArray type="Float32" Name="%s" NumberOfComponents="%d" format="ascii">' % (fieldName, numComponents)
		count = 0
		for k in range(cshape[2]):
			for j in range(cshape[1]):
				for i in range(cshape[0]):
					for n in range(numComponents):
						o.write("%12.4e" % q[i, j, k, n])
						count += 1
						if count % 6 == 0:
							count = 0
							o.write('\n')
		print >>o, '</DataArray>'
	print >>o, '</CellData>'
	print >>o, '</Piece>'
	print >>o, '</StructuredGrid>'
	print >>o, '</VTKFile>'

def WriteVTM(stem, XYZQs, makedir, flatten):

	path_prefix = ""
	if makedir:
		if os.path.exists(stem):
			if not os.path.isdir(stem):
				print "A file named %s exists and it's not directory" % stem
				raise IOError
		else:
			os.mkdir(stem)
		path_prefix = stem + os.sep

	datasets = []
	for Z, XYZ, Qs in XYZQs:
		#if flatten: # FIXME: should not be here
		#	Flatten(XYZ, flatten)
		vtsFileName = "%s%s.%03d.vts" % (path_prefix, stem, Z)
		fvts = open(vtsFileName, "w")
		WriteVTS(fvts, XYZ, Qs)
		fvts.close()
		datasets.append([Z, vtsFileName])

	datasets.sort(lambda a, b: cmp(a[0], b[0]))

	fvtm = open("%s.vtm" % stem, "w")
	print >>fvtm, '<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">'
	print >>fvtm, '<vtkMultiBlockDataSet>'
	for i, z_filename in enumerate(datasets):
		Z, filename = z_filename
		print >>fvtm, '<DataSet index="%d" file="%s" />' % (i, filename)
	print >>fvtm, '</vtkMultiBlockDataSet>'
	print >>fvtm, '</VTKFile>'
	fvtm.close()

def WriteVTMUnstructured(stem, XYZQPieces, zones, makedir):

	path_prefix = ""
	if makedir:
		if os.path.exists(stem):
			if not os.path.isdir(stem):
				print "A file named %s exists and it's not directory" % stem
				raise IOError
		else:
			os.mkdir(stem)
		path_prefix = stem + os.sep

	datasets = []
	for i, XYZQPiece in enumerate(XYZQPieces):
		print "Piece: ", [ a[0] for a in XYZQPiece ]
		Points2ZIJK, ZIJK2Points, ZCells = StructuredToUnstructured(XYZQPiece, zones)
		vtuFileName = "%s%s.%03d.vtu" % (path_prefix, stem, i + 1)
		WriteVTKUnstructured(open(vtuFileName, "w"), XYZQPiece, Points2ZIJK, ZIJK2Points, ZCells)
		datasets.append([i + 1, vtuFileName])

	datasets.sort(lambda a, b: cmp(a[0], b[0]))

	fvtm = open("%s.vtm" % stem, "w")
	print >>fvtm, '<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">'
	print >>fvtm, '<vtkMultiBlockDataSet>'
	for i, z_filename in enumerate(datasets):
		Z, filename = z_filename
		print >>fvtm, '<DataSet index="%d" file="%s" />' % (i, filename)
	print >>fvtm, '</vtkMultiBlockDataSet>'
	print >>fvtm, '</VTKFile>'
	fvtm.close()

def StructuredToPieces(XYZQs, zones):
	""" collects zones connected by 1-to-1 connectivity """

	pieces = [] # sets of Zones

	for zone in zones:
		Z = zone["Zone"]
		if len(pieces) == 0:
			pieces.append(set([Z]))
		for c1to1 in zone["1to1s"]:
			dzn = c1to1["DonorZoneName"]
			dZ = zones.FindZoneByName(dzn)["Zone"]
			if Z == dZ:
				continue
			piece = None
			matches = []
			for p in pieces:
				if Z in p or dZ in p:
					matches.append(p)
			if len(matches) == 1:
				piece = matches[0]
			elif len(matches) > 1:
				piece = matches[0]
				for p in matches[1:]:
					piece |= p
				for p in matches[1:]:
					pieces.remove(p)
			if piece == None:
				piece = set([Z, dZ])
				pieces.append(piece)
			else:
				piece.add(dZ)
	print "Pieces = ", pieces

	XYZQPieces = []
	for piece in pieces:
		XYZQPiece = filter(lambda XYZQ: XYZQ[0] in piece, XYZQs)
		XYZQPieces.append(XYZQPiece)

	return XYZQPieces

def StructuredToUnstructured(XYZQs, zones):

	Points2ZIJK = []
	ZIJK2Points = {} # indexed by Z
	ZCells = {} # indexed by Z

	for Z, XYZ, Qs in XYZQs:
		zoneShape = XYZ.shape[:-1]
		IJK2Points = np.zeros(zoneShape, int)
		IJK2Points[:, :, :] = -1
		ZIJK2Points[Z] = IJK2Points

	vid = 0
	for Z, IJK2Points in ZIJK2Points.items():
		zone = zones.FindZoneByID(Z)
		shape = IJK2Points.shape
		for k in range(shape[2]):
			for j in range(shape[1]):
				for i in range(shape[0]):
					if IJK2Points[i, j, k] < 0:
						IJK2Points[i, j, k] = vid
						Points2ZIJK.append((Z, i, j, k))
						vid += 1

		for c1to1 in zone["1to1s"]:
			# We skip periodic connectivities
			rotAngle = np.array(c1to1["RotationAngle"])
			if np.any(rotAngle != 0.0):
				continue

			sr = np.array(c1to1["Range"]) - np.ones(6, int)
			dr = np.array(c1to1["DonorRange"]) - np.ones(6, int)
			dzn = c1to1["DonorZoneName"]
			dZ = zones.FindZoneByName(dzn)["Zone"]
			transform = c1to1["Transform"]
			a, b, c = transform
			T = np.zeros((3, 3), int)
			T[abs(a) - 1, 0] = a / abs(a)
			T[abs(b) - 1, 1] = b / abs(b)
			T[abs(c) - 1, 2] = c / abs(c)

			donorZIJK2Points = ZIJK2Points[dZ]
			sri = np.sign(sr[3:] - sr[:3])
			sri[sri == 0] = 1
			dri = np.sign(dr[3:] - dr[:3])
			dri[dri == 0] = 1
			srS = sr[:3]
			srE = sr[3:] + sri
			drS = dr[:3]
			drE = dr[3:] + dri
			assert(np.all(sri > 0))
			for k in range(srS[2], srE[2]):
				for j in range(srS[1], srE[1]):
					for i in range(srS[0], srE[0]):
						di, dj, dk = np.array(drS) + np.dot(T, np.array([i, j, k]) - np.array(srS))
						donorZIJK2Points[di, dj, dk] = IJK2Points[i, j, k]

	for Z, XYZ, Qs in XYZQs:
		zoneShape = XYZ.shape[:-1]
		shape = (zoneShape[0] - 1, zoneShape[1] - 1, zoneShape[2] - 1)
		IJK2Points = ZIJK2Points[Z]
		zoneCells = np.zeros((shape[0], shape[1], shape[2], 8), int)
		zoneCells[:, :, :, 0] = IJK2Points[:-1, :-1, :-1]
		zoneCells[:, :, :, 1] = IJK2Points[1:, :-1, :-1]
		zoneCells[:, :, :, 2] = IJK2Points[1:, 1:, :-1]
		zoneCells[:, :, :, 3] = IJK2Points[:-1, 1:, :-1]
		zoneCells[:, :, :, 4] = IJK2Points[:-1, :-1, 1:]
		zoneCells[:, :, :, 5] = IJK2Points[1:, :-1, 1:]
		zoneCells[:, :, :, 6] = IJK2Points[1:, 1:, 1:]
		zoneCells[:, :, :, 7] = IJK2Points[:-1, 1:, 1:]
		ZCells[Z] = zoneCells.reshape(shape[0] * shape[1] * shape[2], 8)

	return Points2ZIJK, ZIJK2Points, ZCells

def WriteVTKUnstructured(o, XYZQs, Points2ZIJK, ZIJK2Points, ZCells):

	Z2XYZQs = {}
	for Z, XYZ, Qs in XYZQs:
		Z2XYZQs[Z] = [XYZ, Qs]

	nCells = 0
	for Z, Cells in ZCells.items():
		nCells += len(Cells)
	print "# of cells = ", nCells

	print >>o, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
	print >>o, '<UnstructuredGrid>'
	print >>o, '<Piece NumberOfPoints="%d" NumberOfCells="%d">' % (len(Points2ZIJK), nCells)
	print >>o, '<Points>'
	print >>o, '<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">'
	for Z, i, j, k in Points2ZIJK:
		XYZ = Z2XYZQs[Z][0]
		print >>o, "%12.4e %12.4e %12.4e" % (XYZ[i, j, k, 0], XYZ[i, j, k, 1], XYZ[i, j, k, 2])
	print >>o, '</DataArray>'
	print >>o, '</Points>'

	print >>o, '<Cells>'
	print >>o, '<DataArray type="Int32" Name="connectivity" format="ascii">'
	for Z, Cells in ZCells.items():
		for vids in Cells:
			print >>o, '%d %d %d %d %d %d %d %d' % (vids[0], vids[1], vids[2], vids[3], vids[4], vids[5], vids[6], vids[7])
	print >>o, '</DataArray>'
	print >>o, '<DataArray type="Int32" Name="offsets" format="ascii">'
	for i in range(nCells):
		print >>o, 8 * (i + 1)
	print >>o, '</DataArray>'
	print >>o, '<DataArray type="UInt8" Name="types" format="ascii">'
	for i in range(nCells):
		print >>o, 12
	print >>o, '</DataArray>'
	print >>o, '</Cells>'

	print >>o, '<CellData>'
	fieldNames = XYZQs[0][2].keys()
	for fieldName in fieldNames:
		Q0 = XYZQs[0][2][fieldName]
		if len(Q0.shape) == 3:
			numComponents = 1
		elif len(Q0.shape) == 4:
			numComponents = Q0.shape[3]
		else:
			raise ValueError
		print >>o, '<DataArray type="Float32" Name="%s" NumberOfComponents="%d" format="ascii">' % (fieldName, numComponents)
		for Z, XYZ, Qs in XYZQs:
			Q = Qs[fieldName]
			if len(Q.shape) == 4:
				q = Q
			elif len(Q.shape) == 3:
				s = Q.shape
				q = Q.reshape(s[0], s[1], s[2], 1)
			print "fieldName, Q, q = ", fieldName, Q.shape, q.shape
			pshape = XYZ.shape[:-1]
			cshape = np.array(pshape) - np.array((1, 1, 1))
			count = 0
			for i in range(cshape[0]):
				for j in range(cshape[1]):
					for k in range(cshape[2]):
						for n in range(numComponents):
							o.write("%12.4e" % q[i, j, k, n])
							count += 1
							if count % 6 == 0:
								count = 0
								o.write('\n')
		print >>o, '</DataArray>'
	print >>o, '</CellData>'

	print >>o, '</Piece>'
	print >>o, '</UnstructuredGrid>'
	print >>o, '</VTKFile>'

def Main(meshFileName, solnFileName, stem, makedir, flatten, unstructured):

	fm = CGNSFile.CGNSFile(meshFileName)
	fs = CGNSFile.CGNSFile(solnFileName)

	zones = fm.ReadZones()

	B = 1
	XYZQs = []
	for zone in zones:
		Z = zone["Zone"]
		XYZ = fm.ReadZoneCoord(B, Z, xyzAxis = -1)
		XYZC = ComputeCellCenters(XYZ)
		Qs = fs.ReadSolution(B, Z)
		print "Zone %d:" % Z, Qs.keys()
		ProcessMomentum(Qs)

		rot = fs.ReadRigidBodyMotion(B, Z)
		print rot
		print "Q keys", Qs.keys()
		ComputeROmegas(XYZC, rot, Qs)
		print "Q keys", Qs.keys()

		XYZQs.append([Z, XYZ, Qs])

	if unstructured:
		XYZQPieces = StructuredToPieces(XYZQs, zones)
		WriteVTMUnstructured(stem, XYZQPieces, zones, makedir)
	else:
		WriteVTM(stem, XYZQs, makedir, flatten)

if __name__ == "__main__":

	parser = optparse.OptionParser()
	parser.add_option("-d", action = "store_true", dest = "makedir", default = True, help = "store vts files in a child directory")
	parser.add_option("-s", action = "store", type = "string", dest = "stem", help = "output file name (without extension)")
	parser.add_option("-i", action = "store", type = "string", dest = "infofile", help = "info file name (.in file used by the solver)")
	parser.add_option("-l", action = "store", type = "string", dest = "flatten", help = "flatten mesh into linear cascade")
	parser.add_option("-u", action = "store_true", dest = "unstructured", default = False, help = "convert mesh into unstructured mesh")
	(options, args) = parser.parse_args()

	meshFileName = args[0]
	solnFileName = args[1]
	if not options.stem:
		options.stem = os.path.splitext(os.path.basename(solnFileName))[0]

	if options.infofile:
		execfile(options.infofile, globals(), config)

	if options.flatten:
		if options.flatten != "Z":
			raise ValueError

	Main(meshFileName, solnFileName, options.stem, options.makedir, flatten = options.flatten, unstructured = options.unstructured)

