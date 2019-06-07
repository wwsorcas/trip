#ifndef __IWAVE_GRIDTRACE__
#define __IWAVE_GRIDTRACE__

#define IWAVE_USE_FMGR

#include "utils.h"
#include "cubic.h"

#ifdef IWAVE_USE_FMGR
#include "iwave_fopen.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif
#include "su.h"
#include "header.h"
#include "segy.h"
#ifdef __cplusplus
}
#endif

#include "gridio.h"
#include "except.hh"

void createExplReflTraces(std::string gridname,
			  std::string tracename,
			  std::string dataname,
			  std::string cwproot,
			  FILE * stream);


#endif
