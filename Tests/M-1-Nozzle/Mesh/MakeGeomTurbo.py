import numpy as np
import StringIO

def ReadAirfoil(f):

	curves = []
	for i in range(2):
		curve = []
		blank = 0
		while True:
			line = f.readline().strip()
			line = line.split()
			if len(line) == 0:
				blank += 1
			if blank == 1:
				continue
			elif blank >= 2:
				break
			x, y = 0.001 * float(line[0]), 0.001 * float(line[1])
			curve.append([x, y])
		curves.append(np.array(curve))

	return curves

def MakeXYZSection(sections):

	sio = StringIO.StringIO()

	print >>sio, "SECTIONAL"
	print >>sio, len(sections)
	for isec, section in enumerate(sections):
		print >>sio, "# %d" % (isec + 1)
		print >>sio, "XYZ"
		print >>sio, len(section)
		z = 0.066
		for x, y, z in section:
			print >>sio, x, y, z
	return sio.getvalue()

def MakeZRCurve(zr):

	sio = StringIO.StringIO()
	print >>sio, len(zr)
	for z, r in zr:
		print >>sio, z, r
	return sio.getvalue()

def GeomTurboRow(row):

	buf = open("template.row.geomTurbo").read()

	buf = buf.replace("$ROW_NAME", row["Name"])
	buf = buf.replace("$PERIODICITY", "%d" % row["Periodicity"])

	bufXYZ = MakeXYZSection(row["Suction"]).strip()
	buf = buf.replace("$SUCTION_SECTIONS", bufXYZ)

	bufXYZ = MakeXYZSection(row["Pressure"]).strip()
	buf = buf.replace("$PRESSURE_SECTIONS", bufXYZ)

	return buf

def GeomTurbo(channel, rows):

	buf = open("template.geomTurbo").read()

	buf = buf.replace("$HUB_CURVE", MakeZRCurve(channel["Hub"]).strip())
	buf = buf.replace("$SHROUD_CURVE", MakeZRCurve(channel["Shroud"]).strip())

	row_specs = ""
	for row in rows:
		row_specs += GeomTurboRow(row)

	buf = buf.replace("$ROW_SPECS", row_specs.strip())

	return buf

def Make3DBlade(curves, rr):
	""" rr is a list of R values """

	curves3D = []
	for curve in curves:
		curvesXYZ = np.zeros([len(rr), curve.shape[0], 3])
		for k, r in enumerate(rr):
			curvesXYZ[k, :, 0:2] = curve[:, :]
			curvesXYZ[k, :, 0] = curve[:, 1]
			curvesXYZ[k, :, 2] = curve[:, 0]
			curvesXYZ[k, :, 1] = r
		curves3D.append(curvesXYZ)
	return curves3D

def Main():

	MidSpan = 0.152 * 0.5
	ScalingCFD = 44.0 / 47.0
	Cax_Nozzle = 0.012 * ScalingCFD
	Gap = 0.003775

	curvesN = ReadAirfoil(open("nozzle.dat"))
	curvesN3D = Make3DBlade(curvesN, [0.074, 0.0765])
	nozzle = {"Name":"Nozzle", "Suction":curvesN3D[0], "Pressure":curvesN3D[1], "Periodicity":44}

	channelCurves = ReadAirfoil(open("channel.dat"))
	channel = {"Hub":channelCurves[0], "Shroud":channelCurves[1]}

	buf = GeomTurbo(channel, [nozzle])
	print buf

if __name__ == "__main__":
	Main()

