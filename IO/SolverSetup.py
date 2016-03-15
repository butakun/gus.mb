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
# $Id: SolverSetup.py 294 2013-08-23 14:11:39Z kato $

#import shutil
import os
from GenerateBCs import Family

def Main(infile):

	OUTFILENAME = "out.cgns"
	TMPFILENAME = "tmp.cgns"

	config = {}
	execfile(infile, globals(), config)

	#shutil.copy(config["MeshFileName"], OUTFILENAME)
	#shutil.copy(config["MeshFileName"], TMPFILENAME)
	#os.system("cp %s %s" % (config["MeshFileName"], OUTFILENAME))
	#os.system("cp %s %s" % (config["MeshFileName"], TMPFILENAME))

if __name__ == "__main__":
	import sys
	Main(sys.argv[1])

