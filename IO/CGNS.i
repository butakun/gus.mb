%module CGNS
%{
#include "CGNSFuncs.h"
#include <iostream>
%}

%include "carrays.i"
%array_class(double, doubleArray)
%array_class(int, intArray)
%array_class(char, charArray)

%include "std_string.i"

// UnitsRead
%typemap(in, numinputs=0) (std::string& mass, std::string& length, std::string& time, std::string& temperature, std::string& angle)
	(std::string _mass, std::string _length, std::string _time, std::string _temperature, std::string _angle) {
	$1 = &_mass;
	$2 = &_length;
	$3 = &_time;
	$4 = &_temperature;
	$5 = &_angle;
}
%typemap(argout) (std::string& mass, std::string& length, std::string& time, std::string& temperature, std::string& angle) {
	PyObject *o1, *o2, *o3, *o4, *o5;
	o1 = PyString_FromString($1->c_str());
	o2 = PyString_FromString($2->c_str());
	o3 = PyString_FromString($3->c_str());
	o4 = PyString_FromString($4->c_str());
	o5 = PyString_FromString($5->c_str());
	$result = PyTuple_New(5);
	PyTuple_SetItem($result, 0, o1);
	PyTuple_SetItem($result, 1, o2);
	PyTuple_SetItem($result, 2, o3);
	PyTuple_SetItem($result, 3, o4);
	PyTuple_SetItem($result, 4, o5);
}

// BaseRead
%typemap(in, numinputs=0) (char* baseName, int* cell_dim, int* phys_dim) (char tempName[33], int temp1, int temp2) {
	$1 = tempName;
	$2 = &temp1;
	$3 = &temp2;
}
%typemap(argout) (char* baseName, int* cell_dim, int* phys_dim) {
	PyObject *o1, *o2, *o3;
	o1 = PyString_FromString($1);
	o2 = PyInt_FromLong(*$2);
	o3 = PyInt_FromLong(*$3);
	$result = PyTuple_New(3);
	PyTuple_SetItem($result, 0, o1);
	PyTuple_SetItem($result, 1, o2);
	PyTuple_SetItem($result, 2, o3);
}

// ZoneRead
%typemap(in, numinputs=0) (char* zoneNameOut, CGNS::cgsize_t* sizeOut) (char tempName[33], CGNS::cgsize_t tempSize[9]) {
	$1 = tempName;
	$2 = tempSize;
}
%typemap(argout) (char* zoneNameOut, CGNS::cgsize_t* sizeOut) {
	PyObject *o1, *o2;
	o1 = PyString_FromString($1);
	o2 = PyTuple_New(9);
	for (int i = 0; i < 9; ++i)
	{
		PyObject* o = PyInt_FromLong($2[i]);
		PyTuple_SetItem(o2, i, o);
	}
	$result = PyTuple_New(2);
	PyTuple_SetItem($result, 0, o1);
	PyTuple_SetItem($result, 1, o2);
}

// ZoneWrite
%typemap(in) int* sizeIn (int tmp[6]) {
	tmp[0] = PyLong_AsLong(PyTuple_GetItem($input, 0));
	tmp[1] = PyLong_AsLong(PyTuple_GetItem($input, 1));
	tmp[2] = PyLong_AsLong(PyTuple_GetItem($input, 2));
	$1 = tmp;
}

// SolInfo
%typemap(in, numinputs=0) (std::string& solname, std::string& location) (std::string tempSolName, std::string tempLocation) {
	$1 = &tempSolName;
	$2 = &tempLocation;
}
%typemap(argout) (std::string& solname, std::string& location) {
	PyObject *o1, *o2;
	o1 = PyString_FromString($1->c_str());
	o2 = PyString_FromString($2->c_str());
	$result = PyTuple_New(2);
	PyTuple_SetItem($result, 0, o1);
	PyTuple_SetItem($result, 1, o2);
}

// FieldInfo
%typemap(in, numinputs=0) (std::string& dataType, std::string& fieldName) (std::string tempDataType, std::string tempFieldName) {
	$1 = &tempDataType;
	$2 = &tempFieldName;
}
%typemap(argout) (std::string& dataType, std::string& fieldName) {
	PyObject *o1, *o2;
	o1 = PyString_FromString($1->c_str());
	o2 = PyString_FromString($2->c_str());
	$result = PyTuple_New(2);
	PyTuple_SetItem($result, 0, o1);
	PyTuple_SetItem($result, 1, o2);
}

// BocoInfo
%typemap(in, numinputs=0) (std::string& bocoName, std::string& bocoType, std::vector<int>& range)
	(std::string tempBocoName, std::string tempBocoType, std::vector<int> tempRange) {
	$1 = &tempBocoName;
	$2 = &tempBocoType;
	$3 = &tempRange;
}
%typemap(argout) (std::string& bocoName, std::string& bocoType, std::vector<int>& range) {
	PyObject *o1, *o2, *o3;
	o1 = PyString_FromString($1->c_str());
	o2 = PyString_FromString($2->c_str());
	size_t rangeSize = $3->size();
	o3 = PyTuple_New(rangeSize);
	for (int i = 0; i < rangeSize; ++i)
	{
		PyObject* o = PyInt_FromLong($3->at(i));
		PyTuple_SetItem(o3, i, o);
	}
	$result = PyTuple_New(3);
	PyTuple_SetItem($result, 0, o1);
	PyTuple_SetItem($result, 1, o2);
	PyTuple_SetItem($result, 2, o3);
}

// C1to1Read
%typemap(in, numinputs=0) (std::string& connectName, std::string& donorName, CGNS::cgsize_t* range, CGNS::cgsize_t* donorRange, int* transform)
	(std::string tempConnectName, std::string tempDonorName, CGNS::cgsize_t tempRange[6], CGNS::cgsize_t tempDonorRange[6], int tempTransform[3] ) {
	$1 = &tempConnectName;
	$2 = &tempDonorName;
	$3 = tempRange;
	$4 = tempDonorRange;
	$5 = tempTransform;
}
%typemap(argout) (std::string& connectName, std::string& donorName, CGNS::cgsize_t* range, CGNS::cgsize_t* donorRange, int* transform) {
	PyObject *o1, *o2, *o3, *o4, *o5;
	o1 = PyString_FromString($1->c_str());
	o2 = PyString_FromString($2->c_str());
	o3 = PyTuple_New(6);
	o4 = PyTuple_New(6);
	o5 = PyTuple_New(3);
	for (int i = 0; i < 6; ++i)
	{
		PyObject* o = PyInt_FromLong($3[i]);
		PyTuple_SetItem(o3, i, o);
		o = PyInt_FromLong($4[i]);
		PyTuple_SetItem(o4, i, o);
	}
	for (int i = 0; i < 3; ++i)
	{	
		PyObject* o = PyInt_FromLong($5[i]);
		PyTuple_SetItem(o5, i, o);
	}

	$result = PyTuple_New(5);
	PyTuple_SetItem($result, 0, o1);
	PyTuple_SetItem($result, 1, o2);
	PyTuple_SetItem($result, 2, o3);
	PyTuple_SetItem($result, 3, o4);
	PyTuple_SetItem($result, 4, o5);
}

// C1to1PeriodicRead
%typemap(in, numinputs=0) (float* rotationCenter, float* rotationAngle, float* translation) (float rotCenter[3], float rotAngle[3], float trans[3]) {
	$1 = rotCenter;
	$2 = rotAngle;
	$3 = trans;
}
%typemap(argout) (float* rotationCenter, float* rotationAngle, float* translation) {
	PyObject *o1, *o2, *o3;
	o1 = PyTuple_New(3);
	o2 = PyTuple_New(3);
	o3 = PyTuple_New(3);
	for (int i = 0; i < 3; ++i)
	{
		PyObject *o;
		o = PyFloat_FromDouble($1[i]);
		PyTuple_SetItem(o1, i, o);
		o = PyFloat_FromDouble($2[i]);
		PyTuple_SetItem(o2, i, o);
		o = PyFloat_FromDouble($3[i]);
		PyTuple_SetItem(o3, i, o);
	}
	$result = PyTuple_New(3);
	PyTuple_SetItem($result, 0, o1);
	PyTuple_SetItem($result, 1, o2);
	PyTuple_SetItem($result, 2, o3);
}

// C1to1Write
%typemap(in) int* range (int tmpRange[6]) {
	for (int i = 0; i < 6; ++i)
	{
		tmpRange[i] = PyInt_AsLong(PyList_GetItem($input, i));
	}
	$1 = tmpRange;
}
%typemap(in) int* donorRange = int* range;
%typemap(in) int* transform (int tmpTransform[3]) {
	for (int i = 0; i < 3; ++i)
	{
		tmpTransform[i] = PyInt_AsLong(PyList_GetItem($input, i));
	}
	$1 = tmpTransform;
}
%typemap(in, numinputs=0) int* I (int tmpI) {
	$1 = &tmpI;
}
%typemap(argout) int* I {
	$result = PyInt_FromLong(*$1);
}

// C1to1PeriodicWrite
%typemap(in) float* {
	int len = PyList_Size($input);
	float* floats = new float[len];
	for (int i = 0; i < len; ++i)
	{
		PyObject* o = PyList_GetItem($input, i);
		float f = float(PyFloat_AsDouble(o));
		floats[i] = f;
	}
	$1 = floats;
}
%typemap(freearg) float* {
	delete[] $1;
}

// FamilyNameRead
%typemap(in, numinputs=0) (std::string& nodeName, std::string& familyName)
	(std::string nodeName_, std::string familyName_) {
	$1 = &nodeName_;
	$2 = &familyName_;
}
%typemap(argout) (std::string& nodeName, std::string& familyName) {
	PyObject *o1, *o2, *o3;
	o1 = PyString_FromString($1->c_str());
	o2 = PyString_FromString($2->c_str());
	$result = PyTuple_New(2);
	PyTuple_SetItem($result, 0, o1);
	PyTuple_SetItem($result, 1, o2);
}

// FamilyRead
%typemap(in, numinputs=0) (std::string& familyName, int* nFamBC, int* nGeo)
	(std::string familyName_, int nFamBC_, int nGeo_) {
	$1 = &familyName_;
	$2 = &nFamBC_;
	$3 = &nGeo_;
}
%typemap(argout) (std::string& familyName, int* nFamBC, int* nGeo) {
	PyObject *o1, *o2, *o3;
	o1 = PyString_FromString($1->c_str());
	o2 = PyInt_FromLong(*$2);
	o3 = PyInt_FromLong(*$3);
	$result = PyTuple_New(3);
	PyTuple_SetItem($result, 0, o1);
	PyTuple_SetItem($result, 1, o2);
	PyTuple_SetItem($result, 2, o3);
}

// FamilyBCRead
%typemap(in, numinputs=0) (std::string& famBCName, std::string& bocoType)
	(std::string famBCName_, std::string bocoType_)
{
	$1 = &famBCName_;
	$2 = &bocoType_;
}
%typemap(argout) (std::string& famBCName, std::string& bocoType) {
	PyObject *o1, *o2;
	o1 = PyString_FromString($1->c_str());
	o2 = PyString_FromString($2->c_str());
	$result = PyTuple_New(2);
	PyTuple_SetItem($result, 0, o1);
	PyTuple_SetItem($result, 1, o2);
}

// FamNameRead
%typemap(in, numinputs=0) (std::string& famName)
	(std::string famName_)
{
	$1 = &famName_;
}
%typemap(argout) (std::string& famName) {
	$result = PyString_FromString($1->c_str());
}

// RigidMotionRead
%typemap(in, numinputs=0) (std::string& rmName, std::string& rmType, double* origin, double* angularVel)
	(std::string rmName_, std::string rmType_, double origin_[3], double angularVel_[3]) {
	$1 = &rmName_;
	$2 = &rmType_;
	$3 = origin_;
	$4 = angularVel_;
}
%typemap(argout) (std::string& rmName, std::string& rmType, double* origin, double* angularVel) {
	PyObject *o1, *o2, *o3, *o4;
	o1 = PyString_FromString($1->c_str());
	o2 = PyString_FromString($2->c_str());
	o3 = PyTuple_New(3);
	o4 = PyTuple_New(3);
	for (int i = 0; i < 3; ++i)
	{
		PyTuple_SetItem(o3, i, PyFloat_FromDouble($3[i]));
		PyTuple_SetItem(o4, i, PyFloat_FromDouble($4[i]));
	}

	$result = PyTuple_New(4);
	PyTuple_SetItem($result, 0, o1);
	PyTuple_SetItem($result, 1, o2);
	PyTuple_SetItem($result, 2, o3);
	PyTuple_SetItem($result, 3, o4);
}

// Where
%typemap(in, numinputs = 0) (std::vector<std::string>& labels, std::vector<int>& indices) (std::vector<std::string> labels_, std::vector<int> indices_) {
	$1 = &labels_;
	$2 = &indices_;
}
%typemap(argout) (std::vector<std::string>& labels, std::vector<int>& indices) {
	PyObject* o1;
	o1 = PyTuple_New($1->size());
	for (int i = 0; i < $1->size(); ++i)
	{
		PyObject *o = PyString_FromString($1->at(i).c_str());
		PyTuple_SetItem(o1, i, o);
	}
	$result = o1;
}

// UserDataRead
%typemap(in, numinputs = 0) std::string& userDefinedDataName (std::string name_) {
	$1 = &name_;
}
%typemap(argout) std::string& userDefinedDataName {
	$result = PyString_FromString($1->c_str());
}

%include "CGNSFuncs.h"

%pythoncode %{

import numpy

def ConvertToString(array, length):
	buf = ""
	for i in range(length):
		if array[i] == '\0':
			break
		buf += array[i]
	return buf

def ConvertToSequence(array, length):
	seq = []
	for i in range(length):
		seq.append(array[i])
	return seq

def CoordRead2(f, B, Z, xyzAxis = 0):
	name, size = ZoneRead(f, B, Z)
	imax, jmax, kmax = size[0:3]
	xx_ = doubleArray(imax * jmax * kmax)
	yy_ = doubleArray(imax * jmax * kmax)
	zz_ = doubleArray(imax * jmax * kmax)
	CoordRead(f, B, Z, imax, jmax, kmax, xx_, yy_, zz_)

	# not efficient but I copy the data to the numpy array.
	if xyzAxis == 0:
		xyz = numpy.zeros((3, imax, jmax, kmax))
		for l, p in enumerate([xx_, yy_, zz_]):
			for i in range(imax):
				for j in range(jmax):
					for k in range(kmax):
						ii = (imax * jmax) * k + imax * j + i
						xyz[l, i, j, k] = p[ii]
	elif xyzAxis == -1:
		xyz = numpy.zeros((imax, jmax, kmax, 3))
		for l, p in enumerate([xx_, yy_, zz_]):
			for i in range(imax):
				for j in range(jmax):
					for k in range(kmax):
						ii = (imax * jmax) * k + imax * j + i
						xyz[i, j, k, l] = p[ii]

	del xx_
	del yy_
	del zz_

	return xyz

def FieldRead2(f, B, Z, S, fieldName):
	name, size = ZoneRead(f, B, Z)
	imax, jmax, kmax = size[0:3]
	imax -= 1
	jmax -= 1
	kmax -= 1
	Q_ = doubleArray(imax * jmax * kmax)
	FieldRead(f, B, Z, S, fieldName, "RealDouble", imax, jmax, kmax, Q_)

	Q = numpy.zeros((imax, jmax, kmax))
	for i in range(imax):
		for j in range(jmax):
			for k in range(kmax):
				ii = (imax * jmax) * k + imax * j + i
				Q[i, j, k] = Q_[ii]

	del Q_

	return Q

%}

