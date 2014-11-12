# $Id: ExportToVTK.py 308 2013-10-02 10:08:55Z kato $

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
		if flatten:
			Flatten(XYZ, flatten)
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

def Main(meshFileName, solnFileName, stem, makedir, flatten):

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

	WriteVTM(stem, XYZQs, makedir, flatten)

if __name__ == "__main__":

	parser = optparse.OptionParser()
	parser.add_option("-d", action = "store_true", dest = "makedir", default = True, help = "store vts files in a child directory")
	parser.add_option("-s", action = "store", type = "string", dest = "stem", help = "output file name (without extension)")
	parser.add_option("-i", action = "store", type = "string", dest = "infofile", help = "info file name (.in file used by the solver)")
	parser.add_option("-l", action = "store", type = "string", dest = "flatten", help = "flatten mesh into linear cascade")
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

	Main(meshFileName, solnFileName, options.stem, options.makedir, flatten = options.flatten)

