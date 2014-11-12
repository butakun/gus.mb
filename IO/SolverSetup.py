# $Id: SolverSetup.py 294 2013-08-23 14:11:39Z kato $

#import shutil
import os

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

