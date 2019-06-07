// MPS_includes.hh
// Author: Mario J. Bencomo
// last modified 11/10/16

/**
 * \file MPS_includes.hh
 * File with handy includes.
 */

#ifndef __MPS_INCLUDES_HH_
#define __MPS_INCLUDES_HH_

#include <math.h>
#include "op.hh"
#include "productspace.hh"
#include "blockop.hh"
#include "ocdc.hh"
#include "iwop.hh"
#include "iwsim.hh"
#include "su.h"
#include "segy.h"
#include "parser.h"
#include "segyops.hh"
#include "gridops.hh"

#ifdef IWAVE_USE_MPI
  #include "mpigridpp.hh"
  #include "mpisegypp.hh"
#else
  #include "gridpp.hh"
  #include "segypp.hh"
#endif

namespace TSOpt{
  
  using std::vector;
  using std::string;
  using std::shared_ptr;
  using std::pow;

  using RVL::Space;
  using RVL::StdSpace;
  using RVL::ProductSpace;
  using RVL::StdProductSpace;
  using RVL::SpaceDCF;
  using RVL::DataContainer;
  using RVL::StdProductDataContainer;
  using RVL::ConstContainer;
  using RVL::Operator;
  using RVL::LinearOp;
  using RVL::LinOpValOp;
  using RVL::LinCompLOVOp;
  using RVL::valparse;
  using RVL::Vector;
  using RVL::Components;
  using RVL::AssignParams;
  using RVL::AssignFilename;
  using RVL::RVLException;
}

#endif
