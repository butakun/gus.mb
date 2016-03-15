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
# $Id: IGGBCS.py 212 2012-03-29 22:27:11Z kato $

import numpy as np

class IGGBCS(object):

	def __init__(self, filename):

		f = open(filename)

		tokens = f.readline().split()
		assert(tokens[0] == "NUMBER_OF_BLOCKS")
		numBlocks = int(tokens[1])
		print numBlocks

		blocks = []
		for i in range(numBlocks):
			iggBlockName, cgnsBlockName, blockType = f.readline().split()
			idim, jdim, kdim = map(int, f.readline().split())
			f.readline() # periodicity
			patches = []
			# faces are ordered as [imin, jmin, kmin, kmax, jmax, imax]
			faceOrder = [
				[[1, 0, 0], [1, 2]],	# imin
				[[0, 1, 0], [0, 2]],	# jmin
				[[0, 0, 1], [0, 1]],	# kmin
				[[0, 0, kdim], [0, 1]],	# kmax
				[[0, jdim, 0], [0, 2]],	# jmax
				[[idim, 0, 0], [1, 2]],	# imax
				]
			for facei, principal in faceOrder:
				numPatchesInThisFace = int(f.readline())
				for k in range(numPatchesInThisFace):
					tokens = f.readline().split()
					patchName, bcType = tokens[0:2]
					if bcType == "CON" or bcType == "NMB" or bcType == "PER" or bcType == "PERNM":
						continue
					start = np.array(facei)
					end = np.array(facei)
					start[principal] = map(int, tokens[2:4])
					end[principal] = map(int, tokens[4:6])
					start -= np.array([1, 1, 1])
					end -= np.array([1, 1, 1])
					patches.append({"Name":patchName, "Type":bcType, "Start":start, "End":end})
			block = {"Name":cgnsBlockName, "Type":blockType, "Patches":patches}
			blocks.append(block)

		self.Blocks = blocks

def Test(filename):

	iggbcs = IGGBCS(filename)

	print iggbcs.Blocks

if __name__ == "__main__":
	import sys
	Test(sys.argv[1])

