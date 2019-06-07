#ifndef __ESG_DEFN__
#define __ESG_DEFN__

#include "esg.hh"
#include "iwinfo.hh"

std::string IWaveInfo::iwave_model="esg";

FIELD IWaveInfo::iwave_fields[]
={
    //                      (z,x,y)
    {"gridinfo",    0,  0,  {0,0,0}},   //  pass "buoyancy" to "gridinfo"
    {"buoyancy",    1,  0,  {0,0,0}},
    {"c11",         1,  0,  {0,0,0}},
    {"c22",         1,  0,  {0,0,0}},
    {"c33",         1,  0,  {0,0,0}},
    {"c23",         1,  0,  {0,0,0}},
    {"c13",         1,  0,  {0,0,0}},
    {"c12",         1,  0,  {0,0,0}},
    {"c44",         1,  0,  {0,0,0}},
    {"c55",         1,  0,  {0,0,0}},
    {"c66",         1,  0,  {0,0,0}},
    {"vx",          1,  0,  {0,0,0}},
    {"vy",          1,  0,  {0,1,1}},
    {"vz",          1,  0,  {1,1,0}},
    {"sxx",         1,  1,  {0,1,0}},
    {"syy",         1,  1,  {0,1,0}},
    {"szz",         1,  1,  {0,1,0}},
    {"syz",         1,  1,  {1,1,1}},
    {"sxz",         1,  1,  {1,0,0}},
    {"sxy",         1,  1,  {0,0,1}},
    {"trace",       1,  0,  {0,0,0}},
    {"",            0,  0,  {0,0,0}}
};


FD_MODELINIT IWaveInfo::minit = esg_modelinit;
FD_MODELDEST IWaveInfo::mdest = esg_modeldest;
FD_TIMEGRID IWaveInfo::timegrid = esg_timegrid;
FD_TIMESTEP IWaveInfo::timestep = esg_timestep;
FD_STENCIL IWaveInfo::createstencil = esg_create_sten;
FD_CHECK IWaveInfo::check = esg_check;
#endif
