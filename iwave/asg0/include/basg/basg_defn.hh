#ifndef __BASG_DEFN__
#define __BASG_DEFN__

#include "basg.hh"
#include "iwinfo.hh"

std::string IWaveInfo::iwave_model = "basg";
FIELD IWaveInfo::iwave_fields[]
= {
  {"bulkmod",    0,    0,  {0, 0, 0}},
  {"buoyancy",   0,    0,  {0, 0, 0}},
  {"p0",         1,    0,  {0, 0, 0}},
  {"p1",         1,    0,  {0, 0, 0}},
  {"p2",         1,    0,  {0, 0, 0}},
  {"v0",         1,    0,  {1, 0, 0}},
  {"v1",         1,    0,  {0, 1, 0}},
  {"v2",         1,    0,  {0, 0, 1}},
  {"dbulkmod",   0,    0,  {0, 0, 0}},  
  {"dbuoyancy",  0,    0,  {0, 0, 0}},
  {"dp0",        1,    0,  {0, 0, 0}},
  {"dp1",        1,    0,  {0, 0, 0}},
  {"dp2",        1,    0,  {0, 0, 0}},
  {"dv0",        1,    0,  {1, 0, 0}},
  {"dv1",        1,    0,  {0, 1, 0}},
  {"dv2",        1,    0,  {0, 0, 1}},
  {"",           0,    0,  {0, 0, 0}}
};

FD_MODELINIT IWaveInfo::minit = basg_modelinit;
FD_MODELDEST IWaveInfo::mdest = asg_modeldest;
FD_TIMESTEP IWaveInfo::timestep = basg_timestep;
FD_TIMEGRID IWaveInfo::timegrid = asg_timegrid;
FD_STENCIL IWaveInfo::createstencil = basg_create_sten;
FD_CHECK IWaveInfo::check = asg_check;
FD_LOOPDEF IWaveInfo::loopdef = asg_loop_refine;

#endif

