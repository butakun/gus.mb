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
def Main():

	RMid = 152.0 * 0.5

	GAP_N1_R1 = 3.775
	X_N1_In = -10.0
	X_N1_R1 = 13.887
	X_R1_TE = 12.0 + GAP_N1_R1 + 7.5
	X_R1_Out = 40.0

	SpanNozzle = 7.4
	SpanRotor1Out = 8.8

	SlabHeight = 0.08
	ExpR1 = SpanRotor1Out / SpanNozzle

	XR_Hub = [
		[X_N1_In, RMid - 0.5 * SlabHeight],
		[X_N1_R1, RMid - 0.5 * SlabHeight],
		[X_R1_TE, RMid - 0.5 * ExpR1 * SlabHeight],
		[X_R1_Out, RMid - 0.5 * ExpR1 * SlabHeight]
		]

	XR_Shroud = [
		[X_N1_In, RMid + 0.5 * SlabHeight],
		[X_N1_R1, RMid + 0.5 * SlabHeight],
		[X_R1_TE, RMid + 0.5 * ExpR1 * SlabHeight],
		[X_R1_Out, RMid + 0.5 * ExpR1 * SlabHeight]
		]

	f = open("channel.dat", "w")
	for xrs in [XR_Hub, [], XR_Shroud]:
		if len(xrs) == 0:
			print >>f
			print >>f
			continue
		for x, r in xrs:
			print >>f, x, r
	f.close()

if __name__ == "__main__":
	Main()

