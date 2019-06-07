#ifndef __IWAVE_ACD_INFO
#define __IWAVE_ACD_INFO

//#include "iwinfo.hh"
#include "iwave.h"
#include "except.hh"
#include "acd_gfdm.h"
#include "acd_gfdm2.h"

#define OLD

/*
  int acd_modelinit(PARARRAY *pars,
  FILE *stream,
  grid const & g,
  ireal dt,
  std::vector<std::string> & active,
  void ** specs);
*/
int acd_modelinit(PARARRAY pars,
		  FILE *stream,
		  IMODEL & model);

void acd_modeldest(void ** specs);

int acd_timegrid(PARARRAY * pars, 
		 FILE * stream, 
		 grid const & g, 
		 ireal & dt,
		 ireal & rhs);

void acd_timestep(std::vector<RDOM *> dom, 
		  bool fwd, 
		  int iv, 
		  void* fdpars);

int acd_create_sten(void *, 
		    FILE *, 
		    int, 
		    IPNT[RDOM_MAX_NARR], 
		    STENCIL *);

void acd_check(RDOM * dom,
	       void * specs,
	       FILE * stream);

void acd_loop_refine(int const * gmin, int const * gmax,
		     float tmax, int input,
		     FILE * stream, void * specs);

typedef struct {
  ireal dt;      /* time step - copied from IMODEL.tsinfo */
  RPNT lam;      /* courant params */
  int k;         /* scheme half-order */
  int ndim;      /* dimension, copied from IMODEL.grid */
  IPNT lbc;      /* flag left boundary conditions */
  IPNT rbc;      /* flag right boundary conditions */
  // FD coefficients - old style
  ireal c0;
  RPNT c1;
  RPNT c2;
  RPNT c3; 
  RPNT c4;
  // FD coefficients - new style
  ireal ** c;
  // workspace for inner loops
  ireal * lap;
  ireal * lap1;
  ireal * lap2;
  ireal * lap3;
  // number of samples in each axis
  IPNT n;
  // axis origins in local indices
  IPNT z;
  // axis steps
  RPNT dx;
  // number of extended axes
  int next;
  // physical background flag
  int pbg;
  // test params
  ireal cmax;
  ireal cmin;
  // additional loop limit workspace
  float looprad;
  IPNT rs; // receiver min index - set to s0 initially
  IPNT re; // receiver max index - set to e0 initially
  IPNT ss; // source min index - set to s0 initially
  IPNT se; // source max index - set to e0 initially
  IPNT dyns; // min loop limits - based on tmax, cmax, dx, source and receiver boxes
  IPNT dyne; // min loop limits - based on tmax, cmax, dx, source and receiver boxes
} ACD_TS_PARS;  

/*
  Indices of arrays - const density acoustic. 
  D_UC  :  current acoustic potential
  D_UP  :  past acoustic potential
  D_CSQ :  square velocity (multiplier of Laplacian)
*/
#define D_UC  1
#define D_UP  2
#define D_CSQ 0

/* default fraction of max time step */
#define CFL_DEF 0.95

#endif


