// $Id: CGNSFuncs.h 321 2014-08-26 09:52:44Z kato $
#ifndef INCLUDE_CGNS_FUNCS_H__
#define INCLUDE_CGNS_FUNCS_H__

#include <string>
#include <vector>
namespace CGNS
{
#include "cgnslib.h"
}

namespace CGNSFuncs
{

int Open(const char* filename, const char* mode);
int Close(int fn);

void UnitsRead(std::string& mass, std::string& length, std::string& time, std::string& temperature, std::string& angle);

int NBases(int fn);
int NZones(int fn, int B);
int NBocos(int fn, int B, int Z);
int N1to1(int fn, int B, int Z);
int N1to1Global(int fn, int B);
int NConns(int fn, int B, int Z);

void BaseRead(int fn, int B, char* baseName, int* cell_dim, int* phys_dim);
int BaseWrite(int fn, char* baseName, int cell_dim, int phys_dim);

void ZoneRead(int fn, int B, int Z, char* zoneNameOut, CGNS::cgsize_t* sizeOut);
int ZoneWrite(int fn, int B, char* zoneName, CGNS::cgsize_t* sizeIn);

void CoordRead(int fn, int B, int Z, int imax, int jmax, int kmax, double* xx, double* yy, double* zz);
void CoordWrite(int fn, int B, int Z, double* xx, double* yy, double* zz);

int NSols(int fn, int B, int Z);
void SolInfo(int fn, int B, int Z, int S, std::string& solname, std::string& location);

int NFields(int fn, int B, int Z, int S);
void FieldInfo(int fn, int B, int Z, int S, int F, std::string& dataType, std::string& fieldName);
void FieldRead(int fn, int B, int Z, int S, const char* fieldName, const char* dataType, int imax, int jmax, int kmax, void* array);

void BocoInfo(int fn, int B, int Z, int BC, std::string& bocoName, std::string& bocoType, std::vector<int>& range);
int BocoWrite(int fn, int B, int Z, char* bocoName, char* bocoType, CGNS::cgsize_t* range);

void C1to1Read(int fn, int B, int Z, int I, std::string& connectName, std::string& donorName, CGNS::cgsize_t* range, CGNS::cgsize_t* donorRange, int* transform);
void C1to1Write(int fn, int B, int Z, const char* connectName, const char* donorName, CGNS::cgsize_t* range, CGNS::cgsize_t* donorRange, int* transform, int* I);
//void C1to1ReadGlobal(int fn, int B, cgns_names& connectNames, cgns_names& zoneNames, cgns_names& donorNames, int* range, int* donorRange, int* transform);

void C1to1PeriodicRead(int fn, int B, int Z, int I, float* rotationCenter, float* rotationAngle, float* translation);
void C1to1PeriodicWrite(int fn, int B, int Z, int I, float* rotationCenterIn, float* rotationAngleIn, float* translationIn);

void ConnRead(int fn, int B, int Z, int I, std::string& connName, std::string& gridLocation, std::string& connectType, std::string& ptsetType, std::string& donorName, std::string& donorZoneType, std::string& donorPtsetType, std::string& donorDataType, std::vector<int>& pnts, std::vector<int>& donorData);

int NFamilies(int fn, int B);
int NFamilyNames(int fn, int B, int Fam);
void FamilyNameRead(int fn, int B, int Fam, int N, std::string& nodeName, std::string& familyName);
void FamilyRead(int fn, int B, int Fam, std::string& familyName, int* nFamBC, int* nGeo);

void FamilyBCRead(int fn, int B, int Fam, int BC, std::string& famBCName, std::string& bocoType);

void FamNameRead(std::string& famName);

int NRigidMotions(int fn, int B, int Z);
void RigidMotionRead(int fn, int B, int Z, int R, std::string& rmName, std::string& rmType, double* origin, double* angularVel);

void GoPath(int fn, char* path);
void Where(std::vector<std::string>& labels, std::vector<int>& indices);
void DeleteNode(char* nodeName);

int NDescriptors();
void DescriptorRead(int D, std::string& name, std::string& text);

int NArrays();
void ArrayInfo(int A, std::string& arrayName, std::string& dataType, int* dataDimension, CGNS::cgsize_t* dimensionVector);
void ArrayReadInteger(int A, std::vector<int>& integerArrayData);
void ArrayReadFloat(int A, std::vector<double>& floatArrayData);

int NUserData();
void UserDataRead(int index, std::string& userDefinedDataName);

}

#endif // INCLUDE_CGNS_FUNCS_H__

