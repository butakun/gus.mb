# $Id: CalculatePoints.py 112 2011-07-20 11:29:06Z kato $

import math
import numpy as np

def Camber(Cax, ThetaTE):

	R = Cax / math.sin(ThetaTE)

	PLE = np.array([0.0, 0.0])
	PTE = np.array([Cax, R * (1.0 - math.cos(ThetaTE))])

	return PLE, PTE, R

def Arc(PLE, PTE, ThetaLE, ThetaTE):

	C = PTE - PLE
	C = math.sqrt(np.dot(C, C))

	dTheta = ThetaTE - ThetaLE

	R = 0.5 * C / math.sin(0.5 * dTheta)

	Pc = np.array([-R * math.sin(ThetaLE), R * math.cos(ThetaLE)])

	return R, Pc

def Main():

	Cax = 1.0
	ThetaLE = 0.0
	ThetaTE = math.radians(30.0)
	deltaLE = math.radians(8.0)
	deltaTE = math.radians(5.0)

	PLE, PTE, R = Camber(Cax, ThetaTE)

	RSS, PcSS = Arc(PLE, PTE, ThetaLE - 0.5 * deltaLE, ThetaTE + 0.5 * deltaTE)
	RPS, PcPS = Arc(PLE, PTE, ThetaLE + 0.5 * deltaLE, ThetaTE - 0.5 * deltaTE)

	print PLE[0], PLE[1], 0.0
	print PTE[0], PTE[1], 0.0
	print RSS, PcSS[0], PcSS[1], 0.0
	print RPS, PcPS[0], PcPS[1], 0.0

if __name__ == "__main__":
	Main()

