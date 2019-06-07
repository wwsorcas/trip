#ifndef __ASG__
#define __ASG__

#include "iwave.h"
#include "except.hh"

int asg_modelinit(PARARRAY pars,
		  FILE *stream,
		  IMODEL & model);

void asg_modeldest(void ** specs);

int asg_timegrid(PARARRAY * pars, 
		 FILE * stream, 
		 grid const & g, 
		 float & dt,
		 float & rhs);

void asg_timestep(std::vector<RDOM *> dom, 
		  bool fwd, 
		  int iv, 
		  void* fdpars);

int asg_create_sten(void *, 
		    FILE *, 
		    int, 
		    IPNT[RDOM_MAX_NARR], 
		    STENCIL *);

void asg_check(RDOM * dom,
	       void * specs,
	       FILE * stream);

void asg_loop_refine(int const * gmin, int const * gmax,
		     float tmax, int input,
		     FILE * stream, void * specs);

typedef struct s_asg_ts_pars {
  float dt;      /* time step - copied from IMODEL.tsinfo */
  RPNT lam;      /* courant params */
  int k;         /* scheme order */
  int ndim;      /* dimension, copied from IMODEL.grid */
  IPNT lbc;      /* flag left boundary conditions */
  IPNT rbc;      /* flag right boundary conditions */
  IPNT nls;      /* number points left pml zone */
  IPNT nrs;      /* number points right pml zone */
  int extd;      /* extended coeff array and adjoint buffer flag*/
  FILE * stream; /* convenience for parallel debugging */
  // FD coeffs
  float ** coeffs;//[RARR_MAX_NDIM];
  // test params
  float cmax;
  float cmin;
  // PML quantities
  float amp;
  float * ep[RARR_MAX_NDIM];     /* precomputed pml arrays */
  float * epp[RARR_MAX_NDIM];    /* p = (1-eta*dt^2) */
  float * ev[RARR_MAX_NDIM];     /* pp = p/(1+eta*dt^2) */
  float * evp[RARR_MAX_NDIM];
  // model extension indicators
  IPNT ihmin;  /* start indices on extended axes */
  IPNT ihmax;  /* stop indices on extended axes */
  RPNT dxs;
} ASG_TS_PARS;  

/* default fraction of max time step */
#define CFL_DEF 0.95

/* special to this app */
void asg_pmlaxis(int n0, int nl, int nr,
                 ireal amp, ireal dt, int gtype,
                 ireal ** ep, ireal ** epp);

#endif


