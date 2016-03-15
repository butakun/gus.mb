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
# $Id$

import CGNSFile

class Range(list):

	def __init__(self, range_):
		list.__init__(self, range_)

def Intersection(r1, r2):

	intersects = r1[0] <= r2[3] and r2[0] <= r1[3] and r1[1] <= r2[4] and r2[1] <= r1[4] and r1[2] <= r2[5] and r2[2] <= r1[5]
	if not intersects:
		return None

class Zone(object):

	def __init__(self, name, range_):
		self.Name = name
		self.Range = range_
		self.Patches = []

	def FindIntersectingPatches(self, range_):
		pass

class Block(object):

	def __init__(self, parent, zone, range_):
		self.Zone = zone
		self.Range = range_ # range in the zone, not local to this block.
		self.Parent = parent
		self.Children = []

class Patch(object):

	def __init__(self, zone1, zone2, begin1, end1, transform, begin2):
		self.Zone1 = zone1
		self.Zone2 = zone2
		self.Begin1 = begin1
		self.End1 = end1
		self.Transform = transform
		self.Begin2 = begin2

