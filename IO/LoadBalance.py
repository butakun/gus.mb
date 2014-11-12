# $Id: LoadBalance.py 314 2013-12-24 15:47:30Z kato $

import CGNSFile
import numpy as np

def Balance(zones, num_bins):

	dimensions = [ (zone["Zone"], zone["Size"][3:6]) for zone in zones ]
	sizes = [ (z[0], z[1][0] * z[1][1] * z[1][2]) for z in dimensions ]
	sizes.sort(lambda a, b: -cmp(a[1], b[1]))

	#print sizes

	bin_sizes = np.zeros(num_bins)
	bins = []
	for i in range(num_bins):
		bins.append([])

	while len(sizes) > 0:
		z = sizes.pop(0)
		ibin = bin_sizes.argmin()
		#print "len(sizes) = ", len(sizes), " appending ", z[0], " to bin ", ibin
		bins[ibin].append(z[0])
		bin_sizes[ibin] += z[1]

	return bins, bin_sizes

def Main(filename, np):

	cgnsfile = CGNSFile.CGNSFile(filename)
	zones = cgnsfile.ReadZones()
	dist = Balance(zones, np)
	print dist

if __name__ == "__main__":
	import sys
	Main(sys.argv[1], int(sys.argv[2]))

