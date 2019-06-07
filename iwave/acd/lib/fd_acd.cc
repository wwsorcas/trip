#include "acd.hh"

#define OLD
//#undef OLD

using RVL::parse;

//#define IWAVE_VERBOSE

/*--- time step functions ---------------------------------------------------*/

extern std::vector<double> fdcoeff(int,
				   double,
				   std::vector<double>);

extern "C" void acd_2d(float **, 
		       float **, 
		       float **, 
		       int *, 
		       int *, 
		       float **,
		       int,
		       int *, 
		       int *,
		       float *);
			   
extern void acd_2d_2(float **, 
		     float **, 
		     float **, 
		     int *, 
		     int *, 
		     float, 
		     float *);

extern void acd_2d_4(float **, 
		     float **, 
		     float **, 
		     int *, 
		     int *, 
		     float, 
		     float *,
		     float *,
		     int *,
		     int *);

extern void acd_2d_8(float **, 
		     float **, 
		     float **, 
		     int *, 
		     int *, 
		     float, 
		     float *,
		     float *,
		     float *,
		     float *,
		     int *,
		     int *);

extern void acd_3d_2(float ***, 
		     float ***, 
		     float ***, 
		     int *, 
		     int *, 
		     float, 
		     float *);

extern void acd_3d_4(float ***, 
		     float ***, 
		     float ***, 
		     int *, 
		     int *, 
		     float, 
		     float *,
		     float *,
		     int *,
		     int *);

extern void acd_3d_8(float ***, 
		     float ***, 
		     float ***, 
		     int *, 
		     int *, 
		     float, 
		     float *,
		     float *,
		     float *,
		     float *,
		     int *,
		     int *);

extern "C" void acde_2d(float **, 
			float **, 
			float ***, 
			int *, 
			int *, 
			int *,
			float **,
			int,
			int,
			int *,
			int *,
			float *,
			float *,
			float *);

extern "C" void acde_2d_4(float **, 
			  float **, 
			  float ***, 
			  int *, 
			  int *, 
			  int *,
			  float, 
			  float *,
			  float *,
			  int *,
			  int *,
			  float *,
			  float *,
			  float *);

extern int acd_step(RDOM*, int, void *);

/*----------------------------------------------------------------------------*/
/* no-ops for this implementation                                             */
/*----------------------------------------------------------------------------*/
int acd_build_sten_dep(FILE * stream, 
		       int ndim, 
		       int stendep[RDOM_MAX_NARR][RDOM_MAX_NARR]) {
  return 0;
}

/*----------------------------------------------------------------------------*/
/* working functions                                                          */
/*----------------------------------------------------------------------------*/

/*int acd_modelinit(PARARRAY *pars, 
  FILE *stream, 
  grid const & g,
  float dt,
  std::vector<std::string> & active,
  void ** fdpars) {
*/
int acd_modelinit(PARARRAY pars, 
		  FILE *stream,
		  IMODEL & model) {

  int err=0;           /* return value */
  ACD_TS_PARS *acdpars;   /* model pars */
  
  IPNT cdims;          /* workspace for cartesian grid dim info */
  IPNT crank;          /* workspace for cartesian grid rank */
#ifdef IWAVE_USE_MPI
  IPNT cpers;          /* workspace for periodic wrap info  - currently not used */
#endif
  
  int i;       /* counter */
  RPNT dxs;    /* grid steps */
  float lam;   /* slownesses dt/dx */
  int idim;    /* counter */
  
  /* allocate sgn model ----------------------------------------------------*/   
  acdpars = (ACD_TS_PARS*)usermalloc_(sizeof(ACD_TS_PARS));
  if ( acdpars == NULL ) { 
    err=E_ALLOC;
    fprintf(stream,"ERROR: acd_modelinit\n");
    fprintf(stream,"failed to allocate SGN_TS_PARS object\n");
    return err;
  }
  
  /* decode dimensions, parallel rank - read grid dimn on rank 0, broadcast */
  
  IASN(cdims, IPNT_1); /* default grid size */ 
  IASN(crank, IPNT_0); /* default cartisian ranks */ 
  
#ifdef IWAVE_USE_MPI
  MPI_Comm cm=retrieveComm();
  
  if ( MPI_Cart_get(cm, RARR_MAX_NDIM, cdims, cpers, crank) != MPI_SUCCESS )  {
    fprintf(stream, "ERROR. Internal: cannot get Cartesian coordinates.\n");
    return E_INTERNAL;
  }
  
  MPI_Bcast((void*)(&(model.g.dim)),1,MPI_INT,0,cm);
#endif
  
  /* set boundary flags */
  IASN(acdpars->lbc, IPNT_0); /* default left bc flag */ 
  IASN(acdpars->rbc, IPNT_0); /* default right bc flag */ 
  
  for (i=0;i<model.g.dim;i++) {
    if (crank[i]==0) acdpars->lbc[i]=1;
    if (crank[i]==cdims[i]-1) acdpars->rbc[i]=1;
  }
  
  /* decode half-order - default=1 */
  acdpars->k = RVL::valparse<int>(pars,"order",1);
#ifdef IWAVE_VERBOSE
  fprintf(stream,"NOTE: initializing ACD with half-order = %d\n",acdpars->k);
#endif
  
  /* decode physical background flag - meaningful only for extended 
     modeling, if set then uses zero offset slice of extended model 
     as physical model in stencils */
  acdpars->pbg = RVL::valparse<int>(pars,"physbg",0);
#ifdef IWAVE_VERBOSE
  fprintf(stream,"NOTE: initializing ACD with physbg flag = %d\n",acdpars->pbg);
#endif
  
  /* set model dimn par */
  acdpars->ndim = model.g.dim;

  /* initialize scaled Courant arrays - old style */
  acdpars->c0=REAL_ONE;
  RASN(acdpars->c1,RPNT_0);
  RASN(acdpars->c2,RPNT_0);
  RASN(acdpars->c3,RPNT_0);
  RASN(acdpars->c4,RPNT_0);

  /* allocate scaled Courant arrays - new style */
  acdpars->c = (float **)usermalloc_(acdpars->ndim * sizeof(float*));
  for (int i=0;i<acdpars->ndim;i++) 
    (acdpars->c)[i]=(float*)usermalloc_((acdpars->k + 1)*sizeof(float));

  /* initialize bound check params */
  acdpars->cmax = RVL::valparse<float>(pars,"cmax");
  acdpars->cmin = RVL::valparse<float>(pars,"cmin");

  /* extract grid steps from grid */
  get_d(dxs, model.g);

  /* assign scaled courant arrays */
  acdpars->c0 = 0.0;

  // new style basic fd coeffs
  double xbar = 0.0;
  std::vector<double> x(2*(acdpars->k)+1);
  for (int i=0;i<2*(acdpars->k)+1;i++) x[i]=double(i-(acdpars->k));
  std::vector<double> cbase = fdcoeff(2,xbar,x);
#ifdef IWAVE_VERBOSE
  fprintf(stream,"RL fdcoeffs\n");
  for (int i=0;i<cbase.size();i++)
    fprintf(stream,"i=%d cbase=%g\n",i,cbase[i]);
#endif

  for (idim = 0;idim < acdpars->ndim;idim ++) {

    if (dxs[idim] <= 0.0) {
      fprintf(stream, "Error: bad input: wrong grid space step, dim=%d, step=%g\n",
	      idim, dxs[idim]);
      return E_BADINPUT;
    }
    lam = model.tsind.dt / dxs[idim];
      
    /* OLD */
    /* assign scaled Courant numbers for orders 2, 4, and 8 - these are the only */
    /* choices implemented */
    if (acdpars->k==1) {
      acdpars->c1[idim]   = lam*lam;
      acdpars->c0        += lam*lam*(-2.0);
    }
    else if (acdpars->k==2) {
      acdpars->c1[idim]   = lam*lam*(4.0/3.0);
      acdpars->c2[idim]   = lam*lam*(-1.0/12.0);
      acdpars->c0        += lam*lam*(-5.0/2.0);
    }
    else if (acdpars->k==4) {
      acdpars->c1[idim]   = lam*lam*(8.0/5.0);
      acdpars->c2[idim]   = lam*lam*(-1.0/5.0);
      acdpars->c3[idim]   = lam*lam*(8.0/315.0);
      acdpars->c4[idim]   = lam*lam*(-1.0/560.0);
      acdpars->c0        += lam*lam*(-205.0/72.0);
    }
    else {
#ifdef OLD
      fprintf(stream,"ERROR: fd_acd parameter setup\n");
      fprintf(stream,"assigned scheme half-order = %d not defined\n",acdpars->k);
      fprintf(stream,"currently defined schemes: half-orders 1, 2, and 4\n");
      fflush(stream);
      return E_BADINPUT;
#else
      fprintf(stream,"NOTE: old style fd coeffs not defined\n");
#endif
    }

    // initialize new style
    for (int i=0;i<=(acdpars->k);i++) (acdpars->c)[idim][i]=lam*lam*cbase[i + (acdpars->k)];

    /* NEW */
#ifdef IWAVE_VERBOSE
    fprintf(stream, "k=%d lam[%d] = %g\n", acdpars->k, idim, lam);
    if (acdpars->k==1)
      fprintf(stream,"c1[%d]=%g\n",
	      idim,acdpars->c1[idim]);
    if (acdpars->k==2) 
      fprintf(stream,"c1[%d]=%g c2[%d]=%g\n",
	      idim,acdpars->c1[idim],idim,acdpars->c2[idim]);
    if (acdpars->k==4) 
      fprintf(stream,"c1[%d]=%g c2[%d]=%g c3[%d]=%g c4[%d]=%g\n",
	      idim,acdpars->c1[idim],idim,acdpars->c2[idim],
	      idim,acdpars->c3[idim],idim,acdpars->c4[idim]);
    fprintf(stream,"new style coefficients:\n");
    for (int i=0; i<=acdpars->k; i++) fprintf(stream,"c[%d][%d]=%g\n",idim,i,(acdpars->c)[idim][i]);
#endif
  }
#ifdef IWAVE_VERBOSE
  fprintf(stream,"c0=%g\n",acdpars->c0);
#endif
  /* reserve a copy of dt for use in source scaling */
  acdpars->dt = model.tsind.dt;
    
  /* identify active fields */
  model.active.resize(3);
  model.active[0]="csq";
  model.active[1]="uc";
  model.active[2]="up";

#ifdef IWAVE_VERBOSE
  for (size_t i=0; i<model.active.size(); i++) 
    cerr<<"active["<<i<<"] = "<<model.active[i]<<endl;
#endif
	
  /* laplacian workspace for inner loops */
  get_n(acdpars->n, model.g);
  acdpars->lap  = (float*)usermalloc_(2*acdpars->n[0] * sizeof(float));
  acdpars->lap1 = (float*)usermalloc_(2*acdpars->n[0] * sizeof(float));
  acdpars->lap2 = (float*)usermalloc_(2*acdpars->n[0] * sizeof(float));
  acdpars->lap3 = (float*)usermalloc_(2*acdpars->n[0] * sizeof(float));

  /* spatial origin indices - NOTE THIS ONLY WORKS FOR 
     EXTENDED AXES BECAUSE MODEL GRID IS PROTOTYPE FOR ALL */
  /* also count internal extended axes */
  acdpars->next=0;
  get_gs(acdpars->z,model.g);
  RASN(acdpars->dx,RPNT_1);
  for (int i=0;i<model.g.gdim;i++) { 
    acdpars->z[i] *= -1;
    acdpars->dx[i]=model.g.axes[i].d;
    if (model.g.axes[i].id > EXTINT-1) (acdpars->next)++;
#ifdef IWAVE_VERBOSE
    fprintf(stream,"n[%d]=%d ",i,acdpars->n[i]);
    fprintf(stream,"z[%d]=%d ",i,acdpars->z[i]);
    fprintf(stream,"d[%d]=%e\n ",i,acdpars->dx[i]);
#endif
  }
#ifdef IWAVE_VERBOSE
  fprintf(stream,"next=%d\n",acdpars->next);
#endif

  /* initialize dynamic loop limits workspace - based on physical grid */
  get_gs(acdpars->rs,model.g);
  get_gs(acdpars->ss,model.g);
  get_gs(acdpars->dyns,model.g);
  get_ge(acdpars->re,model.g);
  get_ge(acdpars->se,model.g);
  get_ge(acdpars->dyne,model.g);
  // default is off
  acdpars->cmin = RVL::valparse<float>(pars,"looprad",0.0f);
  
  /* sanity test for number of extended axes */
  if (acdpars->next > acdpars->ndim-1) {
    fprintf(stream,"ERROR: fd_acd parameter setup\n");
    fprintf(stream,"  too many internal extended axes:\n");
    fprintf(stream,"  have %d but only spatial dim - 1 = %d allowed\n",acdpars->next,acdpars->ndim-1);
    fflush(stream);
    return E_BADINPUT;
  }
  /* assign param object pointer */
  model.specs = (void*)acdpars;
  return 0;
}

/*----------------------------------------------------------------------------*/
void acd_modeldest(void ** fdpars) {

  ACD_TS_PARS *acdpars = (ACD_TS_PARS *)(*fdpars);   /* model pars */  
  userfree_(acdpars->lap);
  userfree_(acdpars->lap1);
  userfree_(acdpars->lap2);
  userfree_(acdpars->lap3);
  for (int i=0;i<acdpars->ndim;i++) userfree_((acdpars->c)[i]);
  userfree_(acdpars->c);
  userfree_(*fdpars);
}

/*----------------------------------------------------------------------------*/
/* this variant so simple that might as well write it from
// scratch, so no need for sten_dep_mat - note that all arrays in this
// app are primal, so no
// need for gtype. Clearly this interface should be refactored so that
// these things can be hidden.
// ndim (from fdpars) gtype and sten_dep_mat should be internal details */

int acd_create_sten(void * fdm,
		    FILE * stream,
		    //      		    IWaveInfo const & ic,
		    int ndim,
		    IPNT gtype[RDOM_MAX_NARR], 
		    STENCIL * sten) {
  ACD_TS_PARS * acdpars = (ACD_TS_PARS *)(fdm);
  STENCIL_MASK mask;/* workspace */
  int nmask;        /* number of masks - dependent pairs of dynamic arrays */
  int ipair[2][2];  /* index pairs for use in defining mask */
  int imask;        /* mask counter */
  int idim;         /* dim counter */
  int iv;           /* mask entry counter */
  int len;          /* length of mask - number of entries */
  int j;            /* counter */
  int k;            /* scheme order */
  IPNT ind;         /* workspace for mask entry */
  int err = 0;

  /* set order variable */
  k = acdpars->k;

  /* initialize index pairs */
  /* first - uc->up */
  ipair[0][0] = D_UC;
  ipair[0][1] = D_UP;
  /* second - up->uc */
  ipair[1][0] = D_UP;
  ipair[1][1] = D_UC;

  /* initialize STENCIL to null stencil */
  sten_setnull(sten);
  
  /* sanity check */
  if (k < 1) {
    fprintf(stream,"ERROR: acd_create_sten - illegal value of k = %d\n",k);
    return E_BADINPUT;
  }

  /* declare number of masks: one for shape of uc->up stencil, one for
  // shape of up->uc stencil */
  nmask=2;

  /* nontrivial STENCIL initialization */
  if ((err = sten_create(sten,nmask))) {
    fprintf(stream,"ERROR: acd_create_sten - failed to create stencil\n");
    return err;
  }

  /* length of stencil is 2k+1 in each direction, but origin is common to all
  // directions, so */
  len = 2*k*ndim+1;
  for (imask=0;imask<nmask;imask++) {
    if ((err = mask_create(&mask, ipair[imask][0], ipair[imask][1], len))) {
      fprintf(stream,"ERROR: acd_create_sten from mask_create\n");
      sten_destroy(sten);
      return err;
    }
    /* "cross stencil" - same in every dimension
    // 2d 4th order - k=0
    // idim = 0
    //   j = 0
    //     iv=0: ind[0]=-1
    //     iv=2: ind[0]= 1
    //     iv=1: ind[0]=-2
    //     iv=3  ind[0]= 2
    // etc
    // eventually iv = ndim*4, ind=IPNT_0 */

    for (idim=0;idim<ndim;idim++) {
      IASN(ind,IPNT_0);
      for (j=0;j<k;j++) {
	/* left half of mask on axis idim */
	ind[idim]=-j-1;
	iv = idim*2*k+j;
	if ((err = mask_set(&mask,iv,ind))) {
	  fprintf(stream,"ERROR: acd_create_sten from mask_set\n");	
	  sten_destroy(sten);
	  return err;
	}	
	/* right half of mask on axis idim */
	ind[idim]=j+1;
	iv = idim*2*k+j+k;
	if ((err = mask_set(&mask,iv,ind))) {
	  fprintf(stream,"ERROR: acd_create_sten from mask_set\n");	
	  sten_destroy(sten);
	  return err;
	}	
      }
    
      IASN(ind,IPNT_0);
      iv=2*k*ndim;
      if ((err = mask_set(&mask,iv,ind))) {
	fprintf(stream,"ERROR: acd_create_sten from mask_set\n");	
	sten_destroy(sten);
	return err;
      }	
    }	
    /*
      fprintf(stderr,"setting mask %d\n",imask);
      fprintf(stderr,"ip = %d\n",mask.ip);
      fprintf(stderr,"ir = %d\n",mask.ir);
      fprintf(stderr,"n  = %d\n",mask.n);
      for (j=0;j<mask.n;j++) 
      fprintf(stderr,"s[%d] = [%d,%d,%d]\n",j,(mask.s[j])[0],(mask.s[j])[1],(mask.s[j])[2]);
    */
    if ((err=sten_set(sten,imask, &mask))) {
      fprintf(stream,"ERROR: acd_create_sten from sten_set\n");	
      sten_destroy(sten);
      return err;
    }
  }
  /*
    if (err=sten_out(sten,stderr,acd_ind2str)) {
    fprintf(stream,"ERROR: acd_create_sten from sten_out\n");	
    sten_destroy(sten);
    return err;
    }
  */
  return 0;
}

/*----------------------------------------------------------------------------*/
/* implements new time grid logic: choose stable time step (max_step set), 
   optionally with reduced cfl - to be called in iwave_construct BEFORE
   any data are read, so must depend on max velo, cannot check here.
*/

//int acd_readtimegrid(PARARRAY *pars, FILE * stream, IMODEL * model) {
//int acd_readtimegrid(PARARRAY *pars, FILE * stream, grid const & g, float & dt) {
int acd_timegrid(PARARRAY *pars, 
		 FILE * stream, 
		 grid const & g, 
		 float & dt,
		 float & rhs) {

  float cmax;                   /* max velo, computed or read from params */
  float cfl = CFL_DEF;          /* default cfl fraction */
  float a;                      /* accumulator for max space step */
  int i;                        /* counter */
  
  /* branch on presence of parameter dt - if dt set in param table,
     use it */
  try {
    dt = RVL::valparse<float>(*pars,"dt");
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: sg_readtimegrid - dt=%12.4e read from param table\n", dt);	
    fprintf(stream,"NOTE: NOT CHECKED FOR STABILITY!\n");
#endif
    rhs=dt*dt;
    return 0; 
  }
  catch (RVL::RVLException & e) {
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: sg_readtimegrid - dt not read from param table\n");	
    fprintf(stream,"NOTE: stable dt computed\n");
#endif      
  }
  
  /* for either method of computing dt, need cfl = proportion
     of max time step to use, also cmax */
  cfl   = RVL::valparse<float>(*pars,"cfl",CFL_DEF);
  try {
    cmax  = RVL::valparse<float>(*pars,"cmax");
  }
  catch (RVL::RVLException & e) {
    fprintf(stream,"ERROR: acd_readtimegrid - param table does not specify cmax\n");
    return E_BADINPUT;
  }      
  
#ifdef IWAVE_VERBOSE
  if (cfl==CFL_DEF) {
    fprintf(stream,"NOTE: acd_readtimegrid\n");
    fprintf(stream,"  using default cfl fraction %g\n",cfl);
  }
#endif
  
  a = g.axes[0].d;
  
  for ( i = 1; i < g.dim; ++i ) a = iwave_min(a, g.axes[i].d);
#ifdef IWAVE_VERBOSE
  fprintf(stream,"NOTE: acd_readtimegrid: ndim=%d min dx=%e cfl=%e\n",g.dim,a,cfl);
  fflush(stream);
#endif
  
  if ( (a < REAL_EPS) || (cfl < REAL_EPS) ) {
    fprintf(stream,"ERROR: acd_readtimegrid - either min dx=%e or cfl=%e "
	    " too small\n", a, cfl);
    return E_BADINPUT;
  }
  
  dt = a*cfl/(cmax*sqrt((float)(g.dim)));
  //cerr << "dt="<< dt << endl;
#ifdef IWAVE_VERBOSE
  fprintf(stream,"NOTE: acd_readtimegrid: on return, dt=%e\n",dt);
  fflush(stream);
#endif
  rhs=dt*dt;
  return 0;
}

int acd_step(RDOM* dom, int iv, void * tspars) {

  /* pointers for 2D case */
  register float ** restrict uc2;
  register float ** restrict up2;
  register float ** restrict csq2;
  /* pointers for 3D case */
  register float *** restrict uc3;
  register float *** restrict up3;
  register float *** restrict csq3;
    
  int ndim;                       /* problem dmn */
  IPNT s, s0, cs0;                /* loop starts  */
  IPNT e, e0, ce0;                /* loop ends */
  float tmp;
  IPNT i;

  /* acd struct */
  ACD_TS_PARS * acdpars = (ACD_TS_PARS *)tspars;
  /* extract dimn info */
  ra_ndim (&(dom->_s[D_UC]),&ndim);
  ra_gse  (&(dom->_s[D_UC]),s,e);
  ra_a_gse(&(dom->_s[D_CSQ]),cs0,ce0);
  ra_a_gse(&(dom->_s[D_UC]),s0,e0);

  // pass extended axes
  for (int i=0;i<acdpars->next;i++) {
    s[i+acdpars->ndim]=cs0[i+acdpars->ndim];
    e[i+acdpars->ndim]=ce0[i+acdpars->ndim];
  }

  //fprintf(stderr,"acd_step: s[0]=%d s[1]=%d s[2]=%d\n",s[0],s[1],s[2]);
  //fprintf(stderr,"acd_step: e[0]=%d e[1]=%d e[2]=%d\n",e[0],e[1],e[2]);

  if (ndim == 2) {

    /* 2D computational arrays */
    uc2   = (dom->_s)[D_UC ]._s2;
    up2   = (dom->_s)[D_UP ]._s2;
    csq2  = (dom->_s)[D_CSQ]._s2;
    csq3  = (dom->_s)[D_CSQ]._s3;
      
    // new arbitrary order interfaces
#ifndef OLD
    // extended modeling
    if (acdpars->next == 1) {
      int pbg = acdpars->pbg;
      acde_2d(uc2, up2, csq3,
	      s, e, acdpars->z,
	      acdpars->c,
	      acdpars->k,
	      pbg,
	      acdpars->lbc, 
	      acdpars->rbc,
	      acdpars->dx,
	      acdpars->lap,
	      acdpars->lap1);
    }
    // physical modeling
    else {
      acd_2d(uc2, up2, csq2,
	     s, e, 
	     acdpars->c, 
	     acdpars->k,
	     acdpars->lbc,
	     acdpars->rbc, 
	     acdpars->lap);
    }

#else // end NEW branch

      // old 2,4,8 th order	
      /* 2nd order case */
    if (acdpars->k == 1) {
      acd_2d_2(uc2, up2, csq2,
	       s, e,
	       acdpars->c0, 
	       acdpars->c1);
    }
    /* 4th order case */
    else if (acdpars->k == 2) {
      if (acdpars->next == 1) {
	acde_2d_4(uc2, up2, csq3,
		  s, e, acdpars->z,
		  acdpars->c0,
		  acdpars->c1, acdpars->c2,
		  acdpars->lbc, acdpars->rbc,
		  acdpars->dx,
		  acdpars->lap,
		  acdpars->lap1);
      } 
      else {
	acd_2d_4(uc2, up2, csq2,
		 s, e, 
		 acdpars->c0, 
		 acdpars->c1, acdpars->c2,
		 acdpars->lbc, acdpars->rbc);
      }
    }
    /* 8th order case */
    else if (acdpars->k == 4) {
      acd_2d_8(uc2, up2, csq2,
	       s, e, 
	       acdpars->c0, 
	       acdpars->c1, acdpars->c2,
	       acdpars->c3, acdpars->c4,
	       acdpars->lbc, acdpars->rbc);
    }
    else {
      fprintf(stderr,"ERROR: acd_step\n");
      fprintf(stderr,"called with half-order != 1, 2, 4\n");
      return E_BADINPUT;
    }

#endif // end OLD branch

    int s00=s0[0]; int e00=e0[0];
    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
      for (i[0]=s00;i[0]<=e00;i[0]++) {
	tmp=((dom->_s)[D_UC]._s2)[i[1]][i[0]];
	((dom->_s)[D_UC]._s2)[i[1]][i[0]]=((dom->_s)[D_UP]._s2)[i[1]][i[0]];
	((dom->_s)[D_UP]._s2)[i[1]][i[0]]=tmp;
      }
    }
      
  }
  else if (ndim == 3) {
      
    uc3   = (dom->_s)[D_UC ]._s3;
    up3   = (dom->_s)[D_UP ]._s3;
    csq3  = (dom->_s)[D_CSQ]._s3;

    /* 2nd order case */
    if (acdpars->k == 1) {
      acd_3d_2(uc3, up3, csq3, 
	       s, e, 
	       acdpars->c0, 
	       acdpars->c1);
    }
    /* 4th order case */
    else if (acdpars->k == 2) {
      acd_3d_4(uc3, up3, csq3,
	       s, e, 
	       acdpars->c0, 
	       acdpars->c1, acdpars->c2,
	       acdpars->lbc, acdpars->rbc);
    }
    /* 8th order case */
    else if (acdpars->k == 4) {
      acd_3d_8(uc3, up3, csq3,
	       s, e, 
	       acdpars->c0, 
	       acdpars->c1, acdpars->c2,
	       acdpars->c3, acdpars->c4,
	       acdpars->lbc, acdpars->rbc);
    }
    else {
      fprintf(stderr,"ERROR: acd_step\n");
      fprintf(stderr,"called with half-order != 1, 2, 4\n");
      return E_BADINPUT;
    }
      
    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	  tmp=((dom->_s)[D_UC]._s3)[i[2]][i[1]][i[0]];
	  ((dom->_s)[D_UC]._s3)[i[2]][i[1]][i[0]]=((dom->_s)[D_UP]._s3)[i[2]][i[1]][i[0]];
	  ((dom->_s)[D_UP]._s3)[i[2]][i[1]][i[0]]=tmp;
	}
      }
    }
  }
  else {
    fprintf(stderr,"ERROR: acd_step\n");
    fprintf(stderr,"called with space dim != 2 or 3\n");
    return E_BADINPUT;
  }
    
  return 0;
}

/*----------------------------------------------------------------------------*/
void acd_check(RDOM * dom, void * specs, FILE * stream) {
  // arrays to test: only csq = rarr[0]
  size_t n = 1; // length of spatial data array
  float vmax = -FLT_MAX;
  float vmin = FLT_MAX;

  /* acd struct */
  ACD_TS_PARS * acdpars = (ACD_TS_PARS *)specs;
  
  /* extract dimn info */
  IPNT s0;
  IPNT e0;
  ra_a_gse(&(dom->_s[D_CSQ]),s0,e0);
  for (int i=0;i<acdpars->ndim;i++) n *= (e0[i]-s0[i]+1);

  size_t start = 0;
  if ((acdpars->ndim==2) && (acdpars->next==1)) start=(acdpars->z[acdpars->ndim])*n;
  if (((acdpars->ndim==2) && (acdpars->next>1)) ||
      ((acdpars->ndim>2) && (acdpars->next>0))) {
    RVLException e;
    e<<"ERROR: acd_check\n";
    e<<"  for internal extension only dim=2 currently implemented\n";
    throw e;
  }
  float * v = &((dom->_s[0]._s0)[start]);
  /* max & min */
  for (size_t i=0;i<n;i++) {
    vmax=iwave_max(v[i],vmax);
    vmin=iwave_min(v[i],vmin);
    if (vmax<0.0f || sqrt(vmax)>acdpars->cmax || 
	vmin<0.0f || sqrt(vmin)<acdpars->cmin) {
      RVLException e;
      e<<"ERROR: acd_check\n";
      e<<"  velocity out of bounds\n";
      e<<"  index i = " << i << "\n";
      e<<"  from start at "<<start<<"\n";
      e<<"  csq max = " << sqrt(vmax) << "\n";
      e<<"  csq min = " << sqrt(vmin) << "\n";
      e<<"  input csq at index "<<i<<" = "<<v[i]<<"\n";
      e<<"  out of bounds ["<<acdpars->cmin<<", "<<acdpars->cmax<<"]\n";
      throw e;
      
    }
  }
  //  fprintf(stderr,"acd_check: tested range [%d,%d] min=%e max=%e\n",start,start+n-1,vmin,vmax); 
}

void acd_loop_refine(int const * gmin, int const * gmax,
		     float tmax, int input,
		     FILE * stream, void * specs) {

  /* acd struct */
  ACD_TS_PARS * ap = (ACD_TS_PARS *)specs;

  // no-op if off
  if (ap->looprad < numeric_limits<float>::epsilon()) {
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: from acd_loop_refine\n");
    fprintf(stream,"  no redefinition of loop limits: looprad=%12.4e\n",ap->looprad);
#endif
    return;
  }
  if (input) {
    for (int i=0; i<ap->ndim; i++) {
      ap->ss[i]=iwave_max(ap->ss[i],gmin[i]);
      ap->se[i]=iwave_min(ap->se[i],gmax[i]);
    }
  }
  else {
    for (int i=0; i<ap->ndim; i++) {
      ap->rs[i]=iwave_max(ap->rs[i],gmin[i]);
      ap->re[i]=iwave_min(ap->re[i],gmax[i]);
    }
  }

  // left: smin-x + rmin-x = tmax*cmax so
  // x = (smin+rmin-tmax*cmax)/2
  // right: x-smax + x-rmax = tmax*cmax
  for (int i=0; i<ap->ndim; i++) {
    ap->dyns[i]=0.5*(ap->ss[i]+ap->rs[i]-(ap->looprad*tmax*ap->cmax/(ap->dx[i]))) - 1;
    ap->dyne[i]=0.5*(ap->se[i]+ap->re[i]+(ap->looprad*tmax*ap->cmax/(ap->dx[i]))) + 1;
  }

#ifdef IWAVE_VERBOSE
  fprintf(stream,"NOTE: from acd_loop_refine\n");
  fprintf(stream,"  possible redefinition of loop limits: looprad=%12.4e\n",ap->looprad);
  fprintf(stream,"  source start indices: ");
  for (int i=0;i<ap->ndim;i++) fprintf(stream,"%d ",ap->ss[i]);
  fprintf(stream,"\n");
  fprintf(stream,"  source end indices: ");
  for (int i=0;i<ap->ndim;i++) fprintf(stream,"%d ",ap->se[i]);
  fprintf(stream,"\n");
  fprintf(stream,"  receiver start indices: ");
  for (int i=0;i<ap->ndim;i++) fprintf(stream,"%d ",ap->rs[i]);
  fprintf(stream,"\n");
  fprintf(stream,"  receiver end indices: ");
  for (int i=0;i<ap->ndim;i++) fprintf(stream,"%d ",ap->re[i]);
  fprintf(stream,"\n");
  fprintf(stream,"  dynamic loop start indices: ");
  for (int i=0;i<ap->ndim;i++) fprintf(stream,"%d ",ap->dyns[i]);
  fprintf(stream,"\n");
  fprintf(stream,"  dynamic loop end indices: ");
  for (int i=0;i<ap->ndim;i++) fprintf(stream,"%d ",ap->dyne[i]);
  fprintf(stream,"\n");
#endif
  return;
}
