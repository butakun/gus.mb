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
# $Id: PeriodicityUtils.py 111 2011-07-06 08:23:44Z kato $

import numpy as np

def SmallestCellSize(XYZ):

	DI = XYZ[:, 1:, :] - XYZ[:, :-1, :]
	DJ = XYZ[:, :, 1:] - XYZ[:, :, :-1]

	DI = DI * DI
	DI = DI.sum(0)

	DJ = DJ * DJ
	DJ = DJ.sum(0)

	dmin = min(DI.min(), DJ.min())
	dmin = np.sqrt(dmin)
	return dmin

def ComputeTolerance(XYZs):

	dmin = SmallestCellSize(XYZs[0])
	for XYZ in XYZs[1:]:
		dmin = min(dmin, SmallestCellSize(XYZ))
	return 0.5 * dmin

def MatchPointAndPatch(p, donorXYZ):
	"""
	Find the index of the donor point closest to the point, p.
	"""

	shapeDonor = donorXYZ.shape
	d = donorXYZ - p[:, np.newaxis, np.newaxis]
	d = d * d
	d = d.sum(0)
	dmin = d.min()
	dmin = np.sqrt(dmin)
	dminloc1d = d.argmin()
	dminloc = (dminloc1d / shapeDonor[2], dminloc1d % shapeDonor[2])
	return dmin, dminloc

def ComputeCellCenters(XYZ):

	return 0.25 * (XYZ[:, :-1, :-1] + XYZ[:, 1:, :-1] + XYZ[:, 1:, 1:] + XYZ[:, :-1, 1:])

def ExtractMaskRanges(connectivity):

	patchIndices = np.unique(connectivity[0, :, :])

	ranges = []
	for ipatch in patchIndices:
		if ipatch < 0:
			continue
		mask = connectivity[0, :, :] == ipatch
		rowMask = map(lambda row: np.any(row), mask)
		colMask = map(lambda row: np.any(row), mask.T)
		imin = rowMask.index(True)
		try:
			imax = rowMask.index(False, imin) - 1
		except ValueError:
			imax = len(rowMask) - 1
		jmin = colMask.index(True)
		try:
			jmax = colMask.index(False, jmin) - 1
		except ValueError:
			jmax = len(colMask) - 1
		assert(np.all(mask[min(imax + 1, len(rowMask)):, min(jmax + 1, len(colMask)):] == False))

		IRange = [imin, imax]
		JRange = [jmin, jmax]
		donorBegin = connectivity[1:, imin, jmin]
		donorEnd = connectivity[1:, imax, jmax]
		donorIRange = [donorBegin[0], donorEnd[0]]
		donorJRange = [donorBegin[1], donorEnd[1]]

		donorIEnd = connectivity[1:, imax, jmin]
		donorJEnd = connectivity[1:, imin, jmax]
		dI = (donorIEnd - donorBegin) / (imax - imin)
		dJ = (donorJEnd - donorBegin) / (jmax - jmin)

		ranges.append([IRange, JRange, ipatch, donorIRange, donorJRange, dI, dJ])

	return ranges

def ConvertCellRangeToVertexRange(IRange, JRange, donorIRange, donorJRange):

	dI = IRange[1] - IRange[0]
	dJ = JRange[1] - JRange[0]
	assert(abs(dI) > 0 and abs(dJ) > 0)

def MatchCells(XYZ, donorXYZs, tol = None):

	if tol == None:
		tol = ComputeTolerance(donorXYZs)

	print "Tolerance = ", tol

	XYZc = ComputeCellCenters(XYZ)
	donorXYZcs = map(ComputeCellCenters, donorXYZs)

	shapeXYZc = XYZc.shape
	connectivityc = np.zeros((3, shapeXYZc[1], shapeXYZc[2]), np.int)
	connectivityc[:, :, :] = -1
	for i in range(shapeXYZc[1]):
		for j in range(shapeXYZc[2]):
			p = XYZc[:, i, j]
			found = False
			duplicate = False
			for idonor, donorXYZc in enumerate(donorXYZcs):
				dmin, dminloc = MatchPointAndPatch(p, donorXYZc)
				if dmin < tol:
					if found:
						duplicate = True
					found = True
					connectivityc[0, i, j] = idonor
					connectivityc[1:, i, j] = dminloc
			if duplicate:
				connectivityc[:, i, j] = -1

	ranges = ExtractMaskRanges(connectivityc)
	return ranges

def MatchGridPoints(XYZ, donorXYZs, tol = None):

	if tol == None:
		tol = ComputeTolerance(donorXYZs)

	print "Tolerance = ", tol

	shapeXYZ = XYZ.shape
	connectivity = np.zeros((3, shapeXYZ[1], shapeXYZ[2]), np.int)
	connectivity[:, :, :] = -1
	print "shape = ", shapeXYZ
	for i in range(shapeXYZ[1]):
		for j in range(shapeXYZ[2]):
			p = XYZ[:, i, j]
			found = False
			duplicate = False
			for idonor, donorXYZ in enumerate(donorXYZs):
				dmin, dminloc = MatchPointAndPatch(p, donorXYZ)
				if dmin < tol:
					if found:
						duplicate = True
					found = True
					connectivity[0, i, j] = idonor
					connectivity[1:, i, j] = dminloc
			assert(found)
			if duplicate:
				connectivity[:, i, j] = -1

			for idonor, donorXYZ in enumerate(donorXYZs):
				dmin, dminloc = MatchPointAndPatch(p, donorXYZ)
				if dmin < tol:
					if found:
						raise Error

