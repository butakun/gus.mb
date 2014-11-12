# $Id: PointRange.py 210 2012-03-29 22:25:59Z kato $

import numpy as np

class PointRange2D(object):

	def __init__(self, lo, hi):
		""" (0, 0, 0) -> (imax, jmax, kmax) (includes imax, jmax, kmax-th vertices) """

		self.Min = np.array(lo)
		self.Max = np.array(hi)

		d = (self.Max - self.Min) == 0
		if d[0] == True:
			self.Idx = np.array([1, 2])
			self.Idx0 = 0
			self.CellMin = self.Min
			self.CellMax = self.Max - np.array([0, 1, 1])
		elif d[1] == True:
			self.Idx = np.array([0, 2])
			self.Idx0 = 1
			self.CellMin = self.Min
			self.CellMax = self.Max - np.array([1, 0, 1])
		elif d[2] == True:
			self.Idx = np.array([0, 1])
			self.Idx0 = 2
			self.CellMin = self.Min
			self.CellMax = self.Max - np.array([1, 1, 0])
		else:
			raise ValueError

		self.Mask = np.ones(self.CellShape(), dtype = bool)

	def Overlaps(self, pr):

		if np.any(self.Idx != pr.Idx) or self.Min[self.Idx0] != pr.Min[self.Idx0]:
			return False

		idx = self.Idx
		return np.all(pr.CellMax[idx] >= self.CellMin[idx]) and np.all(pr.CellMin[idx] <= self.CellMax[idx])

	def Shape(self):

		s = self.Max - self.Min + np.array([1, 1, 1])
		return s[self.Idx]

	def CellShape(self):

		return self.Shape() - np.array([1, 1])

	def Subtract(self, pr):

		if not self.Overlaps(pr):
			return

		lo = self.CellMin[self.Idx]
		hi = self.CellMax[self.Idx]

		lo2 = np.maximum(self.CellMin, pr.CellMin)
		hi2 = np.minimum(self.CellMax, pr.CellMax)
		lo2 = lo2[self.Idx]
		hi2 = hi2[self.Idx]

		self.Mask[lo2[0] - lo[0]:hi2[0] + 1 - lo[0], lo2[1] - lo[1]:hi2[1] + 1 - lo[1]] = False

	def SplitToRects(self):

		rects = []
		rectsActive = []
		numRows = self.Mask.shape[0]
		for irow in range(numRows):
			row = list(self.Mask[irow, :])
			# Find True-strips within this row
			strips = []
			i = 0
			while i < len(row):
				try:
					start = row[i:].index(True)
					try:
						end = row[i + start:].index(False)
					except:
						end = len(row[i + start:])
					strips.append([i + start, i + start + end])
					i = i + start + end
				except:
					i = i + len(row[i:])

			# Compare the strips with the active rects
			rectsToKeep = []
			for strip in strips:
				found = None
				for rect in rectsActive:
					strip0, row0 = rect
					if strip0[0] == strip[0] and strip0[1] == strip[1]:
						found = rect
						break
				if found == None:
					# new rect
					rectsToKeep.append([strip, irow])
				else:
					rectsToKeep.append(found)
					rectsActive.remove(found)
			# Those left in rectsActive are to be archived as they are no longer active.
			for rect in rectsActive:
				strip, rowStart = rect
				rects.append([strip, rowStart, irow])
			rectsActive = rectsToKeep
		for rect in rectsActive:
			strip, rowStart = rect
			rects.append([strip, rowStart, numRows])

		prs = []
		for rect in rects:
			strip, rowStart, rowEnd = rect
			lo2D = [rowStart, strip[0]]
			hi2D = [rowEnd, strip[1]]
			lo, hi = self.Min, self.Max
			lo[self.Idx] = lo2D
			hi[self.Idx] = hi2D
			prs.append(PointRange2D(lo, hi))
		return prs

	def Sextuplets(self):

		lo = self.Min + 1
		hi = self.Max + 1
		return [lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]]

	def __str__(self):

		return "PointRange: (%d, %d, %d) - (%d, %d, %d)" % (self.Min[0], self.Min[1], self.Min[2], self.Max[0], self.Max[1], self.Max[2])

	def __eq__(self, other):

		return np.all(self.Min == other.Min) and np.all(self.Max == other.Max)

def Test():

	pr = PointRange2D([0, 0, 0], [4, 0, 8])
	pr2 = PointRange2D([2, 0, 1], [5, 0, 8])

	pr.Subtract(pr2)

	pr.SplitToRects()

if __name__ == "__main__":
	Test()

