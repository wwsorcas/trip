#ifndef __BASG__
#define __BASG__

#include "asg.hh"

int basg_modelinit(PARARRAY pars,
		  FILE *stream,
		  IMODEL & model);
void basg_timestep(std::vector<RDOM *> dom, 
		  bool fwd, 
		  int iv, 
		  void* fdpars);
int basg_create_sten(void * specs,
		     FILE * stream,
		     int ndim,
		     IPNT gtype[RDOM_MAX_NARR],
		     STENCIL * sten);


#endif


