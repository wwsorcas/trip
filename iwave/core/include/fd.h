/*
  Author: XW 01-15-2010
  wangxin.tom@gmail.com
  fd.h 
  explicit finite difference framework

  Heavily modified by WWS 2012-13
*/
/*===========================================================================*/
/**
 * \file fd.h
 * Explicit finite difference modeling framework
 *
 * Contains the \ref FD_MODEL and \ref FD_TS_PARS structs
 */

#ifndef __XWAVE_FD_H_
#define __XWAVE_FD_H_
/*---------------------------------------------------------------------------*/

#include "std_cpp_includes.hh"
#include "rdomain.h"
#include "exchange.h"
#include "model.h"
#include "gridio.h"
#include "iwinfo.hh"

/* #define VERBOSE */
/** computational domain includes the boundary points or not */
#define INCLUDE_BOUNDARY_PNTS 0
  
#define DUAL_GRID   1
#define PRIMAL_GRID 0

using RVL::parse;

  /** 
      @param[in] ia - index of RARR in RDOM
      @param[in] iv - substep number
      @param[in] ic - iwave info object
      @return - true if RARR ia is updated in substep iv
  */
bool fd_update(int ia, int iv, IWaveInfo const & ic);

/** (implemented) returns true if arg is index in rdomain of a dynamic field
    (i.e. a field updated in the simulation), else false. For
    example, in the constant density acoustic example, the current
    and previous time levels of pressure are dynamic, whereas the
    sound velocity is not. So isdyn(i)=1 if i=0 or 1, 0 else. May
    depend on spatial dimension, so same preconditions as isarr.
    
    Preconditions: fd_model_init has been called
    
    Postconditions: none

    @param[in] fd - FD_MODEL defining scheme
    @param[in] i - array index
    @return 1 if i is index of dynamic field, else 0
    
    Called in fd_modelcrea (only dynamic arrays need ghost cell
    augmentation) and fprint_weqn (iwave/src/model/fd.c)
*/
//int isdyn(FD_MODEL * fd, int i);
// new interface 11.13 - depends only on static data
int fd_isdyn(int i, IWaveInfo const & ic);
int fd_isarr(int i, IMODEL & model, IWaveInfo const & ic);

/** number of substeps defined in scheme */
int fd_numsubsteps(IWaveInfo const & ic);

/** number of dynamic arrays defined in scheme

    @param[in] fd - FD_MODEL defining scheme
    @return number of dynamic arrays

Called in setrecvexchange (iwaveinfo.c)
*/
/*
int fd_narr(FD_MODEL * fd, IWaveInfo const & ic);
*/
  /** reads spatial grid info and records it in data members of \ref
      IMODEL struct. IWAVE presumes that all grids involved in the
      simulation are derived from a common primal grid, in one of two 
      ways:


      (1) Physical grid belonging to the data structure named first in
      the FIELDS array (see iwinfo.h) - this would be the non-extended
      axes of this grid

      (2) indicates whether each non-dynamic grid is (internally)
      extended or not.

      Precondition: all extended non-dynamic fields share the same
      grid, but non-dynamic grids may also be non-extended. All
      dynamic fields and non-dynamic non-extended fields share the
      same non-extended grid, derived from the primal grid.

      Initializes the dynamic grid IMODEL.g (possibly extended),
      IMODEL.gp (non-extended, physical primal grid), also array of
      extend flags.

      Postconditions: IMODEL.g, IMODEL.gp, and IMODEL.extd are initialized.

      @param[in] par - associative array struct containing simulation
      parameters as key=value pairs 

      @param[in] stream - verbose output stream
      
      @param[in] mdl - \ref IMODEL struct, declared at driver level and
      passed as part of IWAVE struct through iwave_construct.

      @param[in] ic - IWaveInfo struct carrying field info

      @return 0 on successful completion, else nonzero error code

      Called in iwave_construct (iwave/core/lib/iwave.c).

      Implemented in iwave/core/lib/fd.cc */

int fd_readgrid(PARARRAY * par, FILE * stream, IMODEL * mdl, IWaveInfo const & ic);

/*----------------------------------------------------------------------------*/
/**
 * General FD model creator (IMODEL struct), suitable for 1st order
 * wave equations FD modeling.  Mimics sgn_modelcrea but more
 * general, does the following operations:
 * 
 *    - read grid info
 *    - make action list
 *    - create stencil
 *    - compute size of the domain on its host processor according
 *    to the parallel info Cartesian grid dimensions (cdims) and
 *    Cartesian rank (crank)
 *    - allocated memory, NOTE: STORAGE ALLOCATION occurs only once,
 *    that is for \ref IMODEL::ld_a and set the dimension
 *    information for IMODEL::ld_r[i], IMODEL::ld_s[i], but let \n
 *    IMODEL::ld_r[i]._s[j]._s = IMODEL::ld_s[i]._s[j]._s =
 *    IMODEL::ld_a[i]._s[j]._s,\n where the first '_s' is RARR type
 *    and the second '_s' is (ireal *) type - etc
 *
 * @param [in] cdims  - (IPNT) cartesian grid dimensions in MPI communicator.
 * @param [in] crank  - (IPNT) cartesian rank of this processor in MPI communicator.
 * @param [in] par    - (PARARRAY *) parameter arrary pointer.
 * @param [in] stream - (FILE *) stream to output comments (created by driver).
 * @paoram [out] model - (IMODEL *)  IMODEL pointer.
 * 
 * @return 0 on successful completion, else error code as in 
 * base/include/utils.h.
 */
int fd_modelcrea(IPNT cdims, IPNT crank, PARARRAY * par, FILE * stream, IMODEL * model, IWaveInfo const & ic);

/*--------------*/
#endif
