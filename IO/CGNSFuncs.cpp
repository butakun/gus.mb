// $Id: CGNSFuncs.cpp 321 2014-08-26 09:52:44Z kato $

#include "CGNSFuncs.h"
#include <iostream>
#include <cassert>

namespace CGNS
{
#include "cgnslib.h"
}

const size_t CGNS_NAME_LENGTH = 33;

typedef struct {
    CGNS::CG_BCType_t type;
    const char* name;
} boco_type_to_string_t;

boco_type_to_string_t boco_type_to_string_map[] = {
    { (CGNS::CG_BCType_t)CG_UserDefined, "UserDefined" },
    { CGNS::CG_BCAxisymmetricWedge, "AxisymmetricWedge" },
    { CGNS::CG_BCDegenerateLine, "DegenerateLine" },
    { CGNS::CG_BCExtrapolate, "Extrapolate" },
    { CGNS::CG_BCDegeneratePoint, "DegeneratePoint" },
    { CGNS::CG_BCDirichlet, "Dirichlet" },
    { CGNS::CG_BCFarfield, "Farfield" },
    { CGNS::CG_BCNeumann, "Neumann" },
    { CGNS::CG_BCGeneral, "General" },
    { CGNS::CG_BCInflow, "Inflow" },
    { CGNS::CG_BCOutflow, "Outflow" },
    { CGNS::CG_BCInflowSubsonic, "InflowSubsonic" },
    { CGNS::CG_BCOutflowSubsonic, "OutflowSubsonic" },
    { CGNS::CG_BCInflowSupersonic, "InflowSupersonic" },
    { CGNS::CG_BCOutflowSupersonic, "OutflowSupersonic" },
    { CGNS::CG_BCSymmetryPlane, "SymmetryPlane" },
    { CGNS::CG_BCTunnelInflow, "TunnelInflow" },
    { CGNS::CG_BCSymmetryPolar, "SymmetryPolar" },
    { CGNS::CG_BCTunnelOutflow, "TunnelOutflow" },
    { CGNS::CG_BCWallViscous, "WallViscous" },
    { CGNS::CG_BCWall, "Wall" },
    { CGNS::CG_BCWallViscousHeatFlux, "WallViscousHeatFlux" },
    { CGNS::CG_BCWallInviscid, "WallInviscid" },
    { CGNS::CG_BCWallViscousIsothermal, "WallViscousIsothermal" },
    { CGNS::CG_FamilySpecified, "FamilySpecified" },
    { (CGNS::CG_BCType_t)0, "Null" }
};

typedef struct {
    CGNS::CG_DataType_t type;
    const char* name;
} data_type_to_string_t;

data_type_to_string_t data_type_to_string_map[] = {
    { CGNS::CG_Integer, "Integer" },
    { CGNS::CG_RealSingle, "RealSingle" },
    { CGNS::CG_RealDouble, "RealDouble" },
    { (CGNS::CG_DataType_t)0, "Null" }
};

const char* boco_type_to_string(CGNS::CG_BCType_t type)
{
    boco_type_to_string_t* p = boco_type_to_string_map;
    for (; p->type != 0; ++p)
        if (p->type == type)
            return p->name;
    return "Unknown";
}

CGNS::CG_BCType_t string_to_boco_type(std::string name)
{
    boco_type_to_string_t* p = boco_type_to_string_map;
    for (; p->type != 0; ++p)
        if (p->name == name)
            return p->type;
    std::cerr << "string_to_boco_type: unknown boco type " << name << ", using CG_UserDefined instead." << std::endl;
    return (CGNS::CG_BCType_t)CG_UserDefined;
}

const char* data_type_to_string(CGNS::CG_DataType_t type)
{
    data_type_to_string_t* p = data_type_to_string_map;
    for (; p->type != 0; ++p)
        if (p->type == type)
            return p->name;
    return "Unknown";
}

const char* grid_location_to_string(CGNS::CG_GridLocation_t location)
{
    switch (location)
    {
    case CGNS::CG_Vertex:
        return "Vertex";
        break;
    case CGNS::CG_CellCenter:
        return "CellCenter";
        break;
    case CGNS::CG_IFaceCenter:
        return "IFaceCenter";
        break;
    case CGNS::CG_JFaceCenter:
        return "JFaceCenter";
        break;
    case CGNS::CG_KFaceCenter:
        return "KFaceCenter";
        break;
    default:
        return "Unknown";
        break;
    }
}

const char* grid_connectivity_type_to_string(CGNS::CG_GridConnectivityType_t type)
{
    switch (type)
    {
    case CGNS::CG_Overset:
        return "Overset";
        break;
    case CGNS::CG_Abutting:
        return "Abutting";
        break;
    case CGNS::CG_Abutting1to1:
        return "Abutting1to1";
        break;
    default:
        return "Unknown";
        break;
    }
}

const char* point_set_type_to_string(CGNS::CG_PointSetType_t type)
{
    switch (type)
    {
    case CGNS::CG_PointRange:
        return "PointRange";
        break;
    case CGNS::CG_PointList:
        return "PointList";
        break;
    default:
        return "Unknown";
        break;
    }
}

const char* zone_type_to_string(CGNS::CG_ZoneType_t type)
{
    switch (type)
    {
    case CGNS::CG_Structured:
        return "Structured";
        break;
    case CGNS::CG_Unstructured:
        return "Unstructured";
        break;
    default:
        return "Unknown";
        break;
    }
}

void cgns_assert(int err)
{
    if (err == 0)
        return;
    std::cerr << CGNS::cg_get_error() << std::endl;
    assert(false);
}

int
CGNSFuncs::Open(const char* filename, const char* mode)
{
    int err, fn, mode_ = CG_MODE_READ;
    std::string modestr(mode);
    if (modestr == "r")
        mode_ = CG_MODE_READ;
    else if (modestr == "w")
        mode_ = CG_MODE_MODIFY;
    else
        assert(false);
    err = CGNS::cg_open(filename, mode_, &fn);
    cgns_assert(err);
    return fn;
}

int
CGNSFuncs::Close(int fn)
{
    int err;
    err = CGNS::cg_close(fn);
    cgns_assert(err);
    return err;
}

void
CGNSFuncs::UnitsRead(std::string& mass, std::string& length, std::string& time, std::string& temperature, std::string& angle)
{
    int err;

    CGNS::CG_MassUnits_t mass_;
    CGNS::CG_LengthUnits_t length_;
    CGNS::CG_TimeUnits_t time_;
    CGNS::CG_TemperatureUnits_t temperature_;
    CGNS::CG_AngleUnits_t angle_;
    err = CGNS::cg_units_read(&mass_, &length_, &time_, &temperature_, &angle_);
    cgns_assert(err);

    switch (mass_)
    {
    case CGNS::CG_Kilogram:
        mass = "Kilogram"; break;
    case CGNS::CG_Gram:
        mass = "Gram"; break;
    case CGNS::CG_Slug:
        mass = "Slug"; break;
    case CGNS::CG_PoundMass:
        mass = "PoundMass"; break;
    case CG_Null:
        mass = "Null"; break;
    case CG_UserDefined:
        mass = "UserDefined"; break;
    default:
        mass = "Unknown"; break;
    }

    switch (length_)
    {
    case CGNS::CG_Meter:
        length = "Meter"; break;
    case CGNS::CG_Centimeter:
        length = "Centimeter"; break;
    case CGNS::CG_Millimeter:
        length = "Millimeter"; break;
    case CGNS::CG_Foot:
        length = "Foot"; break;
    case CGNS::CG_Inch:
        length = "Inch"; break;
    case CG_Null:
        length = "Null"; break;
    case CG_UserDefined:
        length = "UserDefined"; break;
    default:
        length = "Unknown"; break;
    }

    switch (time_)
    {
    case CGNS::CG_Second:
        time = "Second"; break;
    case CG_Null:
        time = "Null"; break;
    case CG_UserDefined:
        time = "UserDefined"; break;
    default:
        time = "Unknown"; break;
    }

    switch (temperature_)
    {
    case CGNS::CG_Kelvin:
        temperature = "Kelvin"; break;
    case CGNS::CG_Celsius:
        temperature = "Celsius"; break;
    case CGNS::CG_Rankine:
        temperature = "Rankine"; break;
    case CGNS::CG_Fahrenheit:
        temperature = "Fahrenheit"; break;
    case CG_Null:
        temperature = "Null"; break;
    case CG_UserDefined:
        temperature = "UserDefined"; break;
    default:
        temperature = "Unknown"; break;
    }

    switch (angle_)
    {
    case CGNS::CG_Degree:
        angle = "Degree"; break;
    case CGNS::CG_Radian:
        angle = "Radian"; break;
    case CG_Null:
        angle = "Null"; break;
    case CG_UserDefined:
        angle = "UserDefined"; break;
    default:
        angle = "Unknown"; break;
    }
}

int
CGNSFuncs::NBases(int fn)
{
    int err, nbases;
    err = CGNS::cg_nbases(fn, &nbases);
    cgns_assert(err);
    return nbases;
}

int
CGNSFuncs::NZones(int fn, int B)
{
    int err, nzones;
    err = CGNS::cg_nzones(fn, B, &nzones);
    cgns_assert(err);
    return nzones;
}

int
CGNSFuncs::NBocos(int fn, int B, int Z)
{
    int err, nbocos;
    err = CGNS::cg_nbocos(fn, B, Z, &nbocos);
    cgns_assert(err);
    return nbocos;
}

int
CGNSFuncs::N1to1(int fn, int B, int Z)
{
    int err, n1to1s;
    err = CGNS::cg_n1to1(fn, B, Z, &n1to1s);
    cgns_assert(err);
    return n1to1s;
}

int
CGNSFuncs::N1to1Global(int fn, int B)
{
    int err, n1to1s;
    err = CGNS::cg_n1to1_global(fn, B, &n1to1s);
    cgns_assert(err);
    return n1to1s;
}

int
CGNSFuncs::NConns(int fn, int B, int Z)
{
    int err, nconns;
    err = CGNS::cg_nconns(fn, B, Z, &nconns);
    cgns_assert(err);
    return nconns;
}

void
CGNSFuncs::BaseRead(int fn, int B, char* baseName, int* cell_dim, int* phys_dim)
{
    int err;
    err = CGNS::cg_base_read(fn, B, baseName, cell_dim, phys_dim);
    cgns_assert(err);
}

int
CGNSFuncs::BaseWrite(int fn, char* baseName, int cell_dim, int phys_dim)
{
    int err;
    int B;
    err = CGNS::cg_base_write(fn, baseName, cell_dim, phys_dim, &B);
    cgns_assert(err);
    return B;
}

void
CGNSFuncs::ZoneRead(int fn, int B, int Z, char* zoneName, CGNS::cgsize_t* size)
{
    int err;
    err = CGNS::cg_zone_read(fn, B, Z, zoneName, size);
    cgns_assert(err);
}

int
CGNSFuncs::ZoneWrite(int fn, int B, char* zoneName, CGNS::cgsize_t* size)
{
    int err;
    int Z;
    err = CGNS::cg_zone_write(fn, B, zoneName, size, CGNS::CG_Structured, &Z);
    cgns_assert(err);
    return Z;
}

void
CGNSFuncs::CoordRead(int fn, int B, int Z, int imax, int jmax, int kmax, double* xx, double* yy, double* zz)
{
    int err;
    CGNS::cgsize_t rangemin[3] = { 1, 1, 1 };
    CGNS::cgsize_t rangemax[3] = { imax, jmax, kmax };
    err = CGNS::cg_coord_read(fn, B, Z, "CoordinateX", CGNS::CG_RealDouble, rangemin, rangemax, xx);
    cgns_assert(err);
    err = CGNS::cg_coord_read(fn, B, Z, "CoordinateY", CGNS::CG_RealDouble, rangemin, rangemax, yy);
    cgns_assert(err);
    err = CGNS::cg_coord_read(fn, B, Z, "CoordinateZ", CGNS::CG_RealDouble, rangemin, rangemax, zz);
    cgns_assert(err);
}

void
CGNSFuncs::CoordWrite(int fn, int B, int Z, double* xx, double* yy, double* zz)
{
    int err;
    int C;
    err = CGNS::cg_coord_write(fn, B, Z, CGNS::CG_RealDouble, "CoordinateX", xx, &C);
    cgns_assert(err);
    err = CGNS::cg_coord_write(fn, B, Z, CGNS::CG_RealDouble, "CoordinateY", yy, &C);
    cgns_assert(err);
    err = CGNS::cg_coord_write(fn, B, Z, CGNS::CG_RealDouble, "CoordinateZ", zz, &C);
    cgns_assert(err);
}

int
CGNSFuncs::NSols(int fn, int B, int Z)
{
    int err;
    int nsols;
    err = CGNS::cg_nsols(fn, B, Z, &nsols);
    cgns_assert(err);

    return nsols;
}

void
CGNSFuncs::SolInfo(int fn, int B, int Z, int S, std::string& solname, std::string& location)
{
    int err;
    char solname_[CGNS_NAME_LENGTH];
    CGNS::CG_GridLocation_t location_;
    err = CGNS::cg_sol_info(fn, B, Z, S, solname_, &location_);
    cgns_assert(err);

    solname = solname_;

#if 0
    switch (location_)
    {
    case CGNS::CG_Vertex:
        location = "Vertex";
        break;
    case CGNS::CG_CellCenter:
        location = "CellCenter";
        break;
    case CGNS::CG_IFaceCenter:
        location = "IFaceCenter";
        break;
    case CGNS::CG_JFaceCenter:
        location = "JFaceCenter";
        break;
    case CGNS::CG_KFaceCenter:
        location = "KFaceCenter";
        break;
    default:
        location = "Unknown";
        break;
    }
#else
    location = grid_location_to_string(location_);
#endif
}

int
CGNSFuncs::NFields(int fn, int B, int Z, int S)
{
    int err;
    int nfields;
    err = CGNS::cg_nfields(fn, B, Z, S, &nfields);
    cgns_assert(err);

    return nfields;
}

void
CGNSFuncs::FieldInfo(int fn, int B, int Z, int S, int F, std::string& dataType, std::string& fieldName)
{
    int err;
    CGNS::CG_DataType_t dataType_;
    char fieldName_[CGNS_NAME_LENGTH];
    err = CGNS::cg_field_info(fn, B, Z, S, F, &dataType_, fieldName_);
    cgns_assert(err);

#if 0
    switch (dataType_)
    {
    case CGNS::CG_Integer:
        dataType = "Integer";
        break;
#if 0
    case CGNS::LongInteger:
        dataType = "LongInteger";
        break;
#endif
    case CGNS::CG_RealSingle:
        dataType = "RealSingle";
        break;
    case CGNS::CG_RealDouble:
        dataType = "RealDouble";
        break;
    default:
        dataType = "Unknown";
        break;
    }
#else
    dataType = data_type_to_string(dataType_);
#endif

    fieldName = fieldName_;
}

void
CGNSFuncs::FieldRead(int fn, int B, int Z, int S, const char* fieldName, const char* dataType, int imax, int jmax, int kmax, void* array)
{
    CGNS::cgsize_t rangemin[3] = { 1, 1, 1 };
    CGNS::cgsize_t rangemax[3] = { imax, jmax, kmax };

    std::string tmp(dataType);
    CGNS::CG_DataType_t dataType_ = CGNS::CG_RealSingle;
    if (tmp == "Integer")
        dataType_ = CGNS::CG_Integer;
    else if (tmp == "RealSingle")
        dataType_ = CGNS::CG_RealSingle;
    else if (tmp == "RealDouble")
        dataType_ = CGNS::CG_RealDouble;
    else
        assert(false);

    int err;
    err = CGNS::cg_field_read(fn, B, Z, S, const_cast<char*>(fieldName), dataType_, rangemin, rangemax, array);
    cgns_assert(err);
}

void
CGNSFuncs::BocoInfo(int fn, int B, int Z, int BC, std::string& bocoName, std::string& bocoType, std::vector<int>& range)
{
    int err;
    char bocoName_[CGNS_NAME_LENGTH];
    CGNS::CG_BCType_t bocoType_;
    CGNS::CG_PointSetType_t ptsetType_;
    CGNS::cgsize_t npnts;
    int NormalIndex[3];
    CGNS::cgsize_t NormalListSize;
    CGNS::CG_DataType_t NormalDataType;
    int ndataset;
    err = CGNS::cg_boco_info(fn, B, Z, BC, bocoName_, &bocoType_, &ptsetType_, &npnts, NormalIndex, &NormalListSize, &NormalDataType, &ndataset);
    cgns_assert(err);

    bocoName = bocoName_;
    bocoType = boco_type_to_string(bocoType_);

    CGNS::cgsize_t* pnts = new CGNS::cgsize_t[npnts * 3];
    err = CGNS::cg_boco_read(fn, B, Z, BC, pnts, NULL);
    cgns_assert(err);
    range.clear();
    for (int i = 0; i < npnts * 3; ++i)
    {
        range.push_back(int(pnts[i]));
    }

    delete[] pnts;
}

int
CGNSFuncs::BocoWrite(int fn, int B, int Z, char* bocoName, char* bocoType, CGNS::cgsize_t* range)
{
    // FIXME: ptsetType are ignored.

    int err;
    CGNS::CG_BCType_t bocoType_ = string_to_boco_type(bocoType);
    CGNS::CG_PointSetType_t ptsetType_ = CGNS::CG_PointRange;
    CGNS::cgsize_t npts = 2;
    int BC;
    err = CGNS::cg_boco_write(fn, B, Z, bocoName, bocoType_, ptsetType_, npts, range, &BC);
    cgns_assert(err);

    return BC;
}

void
CGNSFuncs::C1to1Read(int fn, int B, int Z, int I, std::string& connectName, std::string& donorName, CGNS::cgsize_t* range, CGNS::cgsize_t* donorRange, int* transform)
{
    int err;
    char connectName_[CGNS_NAME_LENGTH];
    char donorName_[CGNS_NAME_LENGTH];
    err = CGNS::cg_1to1_read(fn, B, Z, I, connectName_, donorName_, range, donorRange, transform);
    cgns_assert(err);
    connectName = connectName_;
    donorName = donorName_;
}

void
CGNSFuncs::C1to1PeriodicRead(int fn, int B, int Z, int I, float* rotationCenter, float* rotationAngle, float* translation)
{
    int err;
    err = CGNS::cg_1to1_periodic_read(fn, B, Z, I, rotationCenter, rotationAngle, translation);
    if (err == CG_NODE_NOT_FOUND)
    {
        rotationCenter[0] = rotationCenter[1] = rotationCenter[2] = 0.0;
        rotationAngle[0] = rotationAngle[1] = rotationAngle[2] = 0.0;
        translation[0] = translation[1] = translation[2] = 0.0;
    }
    else
    {
        cgns_assert(err);
    }
}

void
CGNSFuncs::C1to1Write(int fn, int B, int Z, const char* connectName, const char* donorName, CGNS::cgsize_t* range, CGNS::cgsize_t* donorRange, int* transform, int* I)
{
    std::cout << "fn = " << fn << ", B = " << B << ", Z = " << Z << std::endl;
    std::cout << ", connectName = " << connectName << ", donorName = " << donorName << std::endl;
    std::cout << ", range = " << range[0] << ", " << range[1] << ", " << range[2] << ", " << range[3] << ", " << range[4]  << ", " << range[5] << std::endl;
    std::cout << ", donorRange = " << donorRange[0] << ", " << donorRange[1] << ", " << donorRange[2] << ", " << donorRange[3] << ", " << donorRange[4]  << ", " << donorRange[5] << std::endl;
    std::cout << ", transform = " << transform[0] << ", " << transform[1] << ", " << transform[2] << std::endl;

    int err;
    err = CGNS::cg_1to1_write(fn, B, Z, connectName, donorName, range, donorRange, transform, I);
    cgns_assert(err);
}

void
CGNSFuncs::C1to1PeriodicWrite(int fn, int B, int Z, int I, float* rotationCenter, float* rotationAngle, float* translation)
{
    int err;
    err = CGNS::cg_1to1_periodic_write(fn, B, Z, I, rotationCenter, rotationAngle, translation);
    cgns_assert(err);
}

void
CGNSFuncs::ConnRead(int fn, int B, int Z, int I, std::string& connName, std::string& gridLocation, std::string& connectType, std::string& ptsetType, std::string& donorName, std::string& donorZoneType, std::string& donorPtsetType, std::string& donorDataType, std::vector<int>& pnts, std::vector<int>& donorData)
{
    int err;
    char connName_[CGNS_NAME_LENGTH], donorName_[CGNS_NAME_LENGTH];
    CGNS::CG_GridLocation_t gridLocation_;
    CGNS::CG_GridConnectivityType_t connectType_;
    CGNS::CG_PointSetType_t ptsetType_, donorPtsetType_;
    CGNS::cgsize_t npnts_, nDataDonor_;
    CGNS::CG_ZoneType_t donorZoneType_;
    CGNS::CG_DataType_t donorDataType_;
    err = CGNS::cg_conn_info(fn, B, Z, I, connName_, &gridLocation_, &connectType_, &ptsetType_, &npnts_, donorName_, &donorZoneType_, &donorPtsetType_, &donorDataType_, &nDataDonor_);
    cgns_assert(err);

    CGNS::cgsize_t* pnts_ = new CGNS::cgsize_t[npnts_ * 3];
    CGNS::cgsize_t* donorData_ = new CGNS::cgsize_t[nDataDonor_ * 3];
    err = CGNS::cg_conn_read(fn, B, Z, I, pnts_, donorDataType_, donorData_);
    for (CGNS::cgsize_t i = 0; i < npnts_; ++i)
    {
        pnts.push_back(pnts_[3 * i    ]);
        pnts.push_back(pnts_[3 * i + 1]);
        pnts.push_back(pnts_[3 * i + 2]);
    }
    for (CGNS::cgsize_t i = 0; i < nDataDonor_; ++i)
    {
        donorData.push_back(donorData_[3 * i    ]);
        donorData.push_back(donorData_[3 * i + 1]);
        donorData.push_back(donorData_[3 * i + 2]);
    }
    delete[] pnts_;
    delete[] donorData_;

    connName = connName_;
    gridLocation = grid_location_to_string(gridLocation_);
    connectType = grid_connectivity_type_to_string(connectType_);
    ptsetType = point_set_type_to_string(ptsetType_);
    donorName = donorName_;
    donorZoneType = zone_type_to_string(donorZoneType_);
    donorPtsetType = point_set_type_to_string(donorPtsetType_);
    donorZoneType = data_type_to_string(donorDataType_);
}

int
CGNSFuncs::NFamilies(int fn, int B)
{
    int err, nfamilies;
    err = CGNS::cg_nfamilies(fn, B, &nfamilies);
    cgns_assert(err);

    return nfamilies;
}

int
CGNSFuncs::NFamilyNames(int fn, int B, int Fam)
{
    int err, nNames;
    err = CGNS::cg_nfamily_names(fn, B, Fam, &nNames);
    cgns_assert(err);

    return nNames;
}

void
CGNSFuncs::FamilyNameRead(int fn, int B, int Fam, int N, std::string& nodeName, std::string& familyName)
{
    int err;
    char nodeName_[CGNS_NAME_LENGTH], familyName_[CGNS_NAME_LENGTH];
    err = CGNS::cg_family_name_read(fn, B, Fam, N, nodeName_, familyName_);
    cgns_assert(err);
    nodeName = nodeName_;
    familyName = familyName_;
    std::cout << "Fam = " << Fam << ", N = " << N << ", nodeName = " << nodeName << ", familyName = " << familyName << std::endl;
}

void
CGNSFuncs::FamilyRead(int fn, int B, int Fam, std::string& familyName, int* nFamBC, int* nGeo)
{
    int err;
    char familyName_[CGNS_NAME_LENGTH];
    err = CGNS::cg_family_read(fn, B, Fam, familyName_, nFamBC, nGeo);
    cgns_assert(err);
    familyName = familyName_;
}

void
CGNSFuncs::FamilyBCRead(int fn, int B, int Fam, int BC, std::string& famBCName, std::string& bocoType)
{
    int err;
    char famBCName_[CGNS_NAME_LENGTH];
    CGNS::CG_BCType_t bocoType_;
    err = CGNS::cg_fambc_read(fn, B, Fam, BC, famBCName_, &bocoType_);
    cgns_assert(err);

    famBCName = famBCName_;
    bocoType = boco_type_to_string(bocoType_);
}

void
CGNSFuncs::FamNameRead(std::string& famName)
{
    int err;
    char famName_[CGNS_NAME_LENGTH];
    err = CGNS::cg_famname_read(famName_);
    cgns_assert(err);

    famName = famName_;
}

int
CGNSFuncs::NRigidMotions(int fn, int B, int Z)
{
    int err, nmotions;
    err = CGNS::cg_n_rigid_motions(fn, B, Z, &nmotions);
    cgns_assert(err);

    return nmotions;
}

void
CGNSFuncs::RigidMotionRead(int fn, int B, int Z, int R, std::string& rmName, std::string& rmType, double* origin, double* angularVel)
{
    int err;
    char rmName_[CGNS_NAME_LENGTH];
    CGNS::CG_RigidGridMotionType_t rmType_;
    err = CGNS::cg_rigid_motion_read(fn, B, Z, R, rmName_, &rmType_);
    cgns_assert(err);

    rmName = rmName_;
    switch (rmType_)
    {
    case CG_UserDefined:
        rmType = "UserDefined";
        break;
    case CGNS::CG_ConstantRate:
        rmType = "ConstantRate";
        break;
    case CGNS::CG_VariableRate:
        rmType = "VariableRate";
        break;
    default:
        rmType = "Unknown";
        break;
    }

    err = CGNS::cg_goto(fn, B, "Zone_t", Z, "RigidGridMotion_t", R, "end");
    cgns_assert(err);
    int narrays;
    err = CGNS::cg_narrays(&narrays);
    cgns_assert(err);
    for (int A = 1; A <= narrays; ++A)
    {
        char arrayName[CGNS_NAME_LENGTH];
        CGNS::CG_DataType_t dataType;
        int dim;
        CGNS::cgsize_t dims[12];
        err = CGNS::cg_array_info(A, arrayName, &dataType, &dim, dims);
        cgns_assert(err);
        if (std::string(arrayName) == "OriginLocation")
        {
            double origins[6];
            err = CGNS::cg_array_read(A, origins);
            cgns_assert(err);
            origin[0] = origins[0];
            origin[1] = origins[1];
            origin[2] = origins[2];
        }
        if (std::string(arrayName) == "RigidRotationRate")
        {
            double omega[3];
            err = CGNS::cg_array_read(A, omega);
            cgns_assert(err);
            angularVel[0] = omega[0];
            angularVel[1] = omega[1];
            angularVel[2] = omega[2];
        }
    }
}

void
CGNSFuncs::GoPath(int fn, char* path)
{
    int err;
    err = CGNS::cg_gopath(fn, path);
    cgns_assert(err);
}

void
CGNSFuncs::Where(std::vector<std::string>& labels, std::vector<int>& indices)
{
    int err;
    int fn, B, depth;
    char labels_[CG_MAX_GOTO_DEPTH][CGNS_NAME_LENGTH];
    int indices_[CG_MAX_GOTO_DEPTH];
    err = CGNS::cg_where(&fn, &B, &depth, (char**)labels_, indices_);
    cgns_assert(err);

    labels.clear();
    indices.clear();
    for (int i = 0; i < depth; ++i)
    {
        labels.push_back(std::string(labels_[i]));
        indices.push_back(indices_[i]);
    }
}

void
CGNSFuncs::DeleteNode(char* nodeName)
{
    int err;
    err = CGNS::cg_delete_node(nodeName);
    cgns_assert(err);
}

int
CGNSFuncs::NDescriptors()
{
    int err, nd;
    err = CGNS::cg_ndescriptors(&nd);
    cgns_assert(err);
    return nd;
}

void
CGNSFuncs::DescriptorRead(int D, std::string& name, std::string& text)
{
    int err;
    char name_[CGNS_NAME_LENGTH];
    char* text_;
    err = CGNS::cg_descriptor_read(D, name_, &text_);
    cgns_assert(err);

    name = name_;
    text = text_;

    CGNS::cg_free(text_);
}

int
CGNSFuncs::NArrays()
{
    int err, narrays;
    err = CGNS::cg_narrays(&narrays);
    cgns_assert(err);

    return narrays;
}

void
CGNSFuncs::ArrayInfo(int A, std::string& arrayName, std::string& dataType, int* dataDimension, CGNS::cgsize_t* dimensionVector)
{
    int err;
    char arrayName_[CGNS_NAME_LENGTH];
    //CGNS::cgsize_t dimensionVector_;
    CGNS::CG_DataType_t dataType_;
    err = CGNS::cg_array_info(A, arrayName_, &dataType_, dataDimension, dimensionVector);
    cgns_assert(err);

    arrayName = arrayName_;
    dataType = data_type_to_string(dataType_);
}

void
CGNSFuncs::ArrayReadInteger(int A, std::vector<int>& integerArrayData)
{
    int err;
    char arrayName[CGNS_NAME_LENGTH];
    CGNS::CG_DataType_t dataType;
    int dataDimension;
    CGNS::cgsize_t dimensionVector;
    err = CGNS::cg_array_info(A, arrayName, &dataType, &dataDimension, &dimensionVector);
    cgns_assert(err);
    assert(dataType == CGNS::CG_RealSingle);

    size_t nData = dataDimension * dimensionVector;
    int* data = new int[nData];
    err = CGNS::cg_array_read(A, (void *)data);
    cgns_assert(err);

    for (size_t i = 0; i < nData; ++i)
    {
        integerArrayData.push_back(data[i]);
    }

    delete[] data;
}

void
CGNSFuncs::ArrayReadFloat(int A, std::vector<double>& floatArrayData)
{
}

int
CGNSFuncs::NUserData()
{
    int err, nuserdata;
    err = CGNS::cg_nuser_data(&nuserdata);
    cgns_assert(err);

    return nuserdata;
}

void
CGNSFuncs::UserDataRead(int index, std::string& userDefinedDataName)
{
    int err;
    char name_[CGNS_NAME_LENGTH];
    err = CGNS::cg_user_data_read(index, name_);
    cgns_assert(err);

    userDefinedDataName = name_;
}

