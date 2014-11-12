import glob, re

def MakeVTM(o, stem, blocks, timestep):

	print >>o, '<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">'
	print >>o, '<vtkMultiBlockDataSet>'

	for i, b in enumerate(blocks):
		print >>o, '<DataSet index="%d" file="%s.%s.%s.vts"/>' % (i, stem, b, timestep)

	print >>o, '</vtkMultiBlockDataSet>'
	print >>o, '</VTKFile>'

def Main(stem):

	files = glob.glob("%s.*.*.vts" % stem)

	compiled = re.compile("%s.([0-9]+).([0-9]+)" % stem)

	blocks, timesteps = [], []
	for f in files:
		m = compiled.match(f)
		blocks.append(m.group(1))
		timesteps.append(m.group(2))

	blocks = list(set(blocks))
	blocks.sort()
	timesteps = list(set(timesteps))
	timesteps.sort()

	print blocks
	print timesteps

	for timestep in timesteps:
		f = open("%s.%s.vtm" % (stem, timestep), "w")
		MakeVTM(f, stem, blocks, timestep)

if __name__ == "__main__":
	Main("out")

