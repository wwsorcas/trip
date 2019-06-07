#ifndef __ESG__
#define __ESG__

#include "iwave.h"
#include "except.hh"

int esg_modelinit(PARARRAY pars,
                  FILE *stream,
                  IMODEL & model);

void esg_modeldest(void ** specs);

int esg_timegrid(PARARRAY * pars,
                 FILE * stream,
                 grid const & g,
                 float & dt,
                 float & rhs);

void esg_timestep(std::vector<RDOM *> dom,
                  bool fwd,
                  int iv,
                  void* fdpars);

int esg_create_sten(void *,
                    FILE *,
                    int,
                    IPNT[RDOM_MAX_NARR],
                    STENCIL *);

void esg_check(RDOM * dom,
               void * specs,
               FILE * stream);

typedef struct s_esg_ts_pars {
    
    ireal cp; //  May not need in esg
    ireal cs; //  May not need in esg
    int k;      //  scheme order
    int ndim;   //  dimension, copied from IMODEL.grid
    ireal dt;   //  time step - copied from IMODEL.tsinfo
    IPNT lbc;   //  flag left boundary conditions
    IPNT rbc;   //  flag right boundary conditions
    ireal ** coeffs;    // FD coeffs
    int eflag;  //  1: print out energy(t), 0: no print
    int effective_media;  //  1: using effective media pars, 0: not using them
    ireal CellSize; //  grid cell size
    ireal h;    // grid step: if using the same grid step
    
} ESG_TS_PARS;

/* default fraction of max time step */
#define CFL_DEF 0.95

#endif

