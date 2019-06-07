#include "asg.hh"

//#define VERBOSE

using RVL::valparse;
using RVL::RVLException;


enum{
    
    D_SIZE=8,
    D_BULK=0,
    D_BUOY=1,
    D_P0=2,
    D_P1=3,
    D_P2=4,
    D_V0=5,
    D_V1=6,
    D_V2=7,
    
};


// max FD half-order
#define MAXK 12

extern float * sgcoeffs(int k);

extern "C" void asg2dtap(float *** bulk, float ** buoy,
			 float ** p0, float ** p1,
			 float ** v0, float ** v1,
			 float ** ep, float ** epp,
			 float ** ev, float ** evp,
			 float ** sdiv, float ** sdivb,
			 float * gradp0, float * gradp1,
			 int * gsc_p, int * gec_p, 
			 int * gsc_v0, int * gec_v0,
			 int * gsc_v1, int * gec_v1,
			 int * lbc, int * rbc,
			 float dh, int ihmax,
			 int maxoff, float ** c,
			 FILE * stream);

extern "C" void asg2dtap_d(float *** restrict bulk, float *** restrict bulkd, float ** restrict buoy, float ** restrict buoyd, float ** restrict p0, float **  restrict p0d, float ** restrict p1, float ** restrict p1d, float ** restrict v0, float ** restrict v0d, float ** restrict v1, float ** restrict v1d, float ** restrict ep, float ** restrict epp, float ** restrict ev, float ** restrict evp, float ** restrict sdiv, float ** restrict sdivd, float ** restrict sdivb, float ** restrict sdivbd, float *gradp0, float *gradp0d, float *gradp1, float *gradp1d, int *gsc_p, int *gec_p, int *gsc_v0, int *gec_v0, int *gsc_v1, int *gec_v1, int *lbc, int *rbc, float dh, int ihmax, int maxoff, float ** restrict c, FILE * stream);

extern "C" void asg2dtap_b(float *** restrict bulk, float *** restrict bulkb, float ** 
			   restrict buoy, float ** restrict buoyb, float ** restrict p0, float **
			   restrict p0b, float ** restrict p1, float ** restrict p1b, float ** 
			   restrict v0, float ** restrict v0b, float ** restrict v1, float ** 
			   restrict v1b, float ** restrict ep, float ** restrict epp, float ** 
			   restrict ev, float ** restrict evp, float ** restrict sdiv, float ** 
			   restrict sdivb0, float ** restrict sdivb, float ** restrict sdivbb, 
			   float *gradp0, float *gradp0b, float *gradp1, float *gradp1b, int *
			   gsc_p, int *gec_p, int *gsc_v0, int *gec_v0, int *gsc_v1, int *gec_v1,
			   int *lbc, int *rbc, float dh, int ihmax, int maxoff, float ** restrict c, FILE *
		stream);

/*
extern "C" void asg_pstep3d(float *** restrict bulk,
                            float *** restrict p0,
			    float *** restrict p1,
			    float *** restrict p2,
                            float *** restrict v0,
			    float *** restrict v1,
			    float *** restrict v2,
                            float ** restrict ep,
			    float ** restrict epp,
                            float * restrict sdiv,
                            int * gsc, int * gec,
                            int * lbc, int * rbc,
                            int maxoff, float ** restrict c);

extern "C" void asg_vstep3d(float *** restrict buoy,
                            float *** restrict p0,float *** restrict p1, float *** restrict p2,
                            float *** restrict v0,float *** restrict v1, float *** restrict v2,
                            float ** restrict ev, float ** restrict evp,
                            float ** restrict gradp,
                            int * gsc_v0, int * gec_v0,
                            int * gsc_v1, int * gec_v1,
                            int * gsc_v2, int * gec_v2,
                            int * lbc, int * rbc,
                            int maxoff,float ** restrict c);

*/
/*
 int asg_modelinit(PARARRAY *pars,
 FILE *stream,
 grid const & g,
 ireal dt,
 std::vector<std::string> & active,
 void ** specs) {
 */
int asg_modelinit(PARARRAY pars,
                  FILE *stream,
                  IMODEL & model) {
  try {
        
    int err=0;           /* return value */
    ASG_TS_PARS *asgpars;   /* model pars */
        
    IPNT cdims;          /* workspace for cartesian grid dim info */
    IPNT crank;          /* workspace for cartesian grid rank */
#ifdef IWAVE_USE_MPI
    IPNT cpers;          /* workspace for periodic wrap info  - currently not used */
#endif

    // moved to private class data 2017.11.10 WWS
    //    RPNT dxs;    /* grid steps */
    ireal lam;   /* slownesses dt/dx */
        
    /* allocate private data */
    asgpars = (ASG_TS_PARS*)usermalloc_(sizeof(ASG_TS_PARS));
    if ( asgpars == NULL ) {
      err=E_ALLOC;
      fprintf(stream,"ERROR: asg_modelinit\n");
      fprintf(stream,"failed to allocate SGN_TS_PARS object\n");
      return err;
    }
        
    asgpars->stream = stream;
        
    /* The next block of code determines which boundary
     * faces of the domain assigned to the current process are
     * external, and which internal.  this code should be moved to the
     * core, along with the boundary flags lbc and rbc.  
     */

     /************ BEGIN *************/
    
    /* decode dimensions, parallel rank - read grid dimn on rank 0,
       broadcast */
        
    IASN(cdims, IPNT_1); /* default grid size */
    IASN(crank, IPNT_0); /* default cartisian ranks */
        
#ifdef IWAVE_USE_MPI
    MPI_Comm cm=retrieveComm();
        
    if ( MPI_Cart_get(cm, RARR_MAX_NDIM, cdims, cpers, crank)
	 != MPI_SUCCESS )  {
      RVLException e;
      e<<"ERROR: asg_modelinit\n";
      e<<"  cannot get Cartesian coordinates.\n";
      throw e;
    }

    /* model has been initialized on root */
    MPI_Bcast((void*)(&(model.g.dim)),1,MPI_INT,0,cm);
#endif
        
    /* set boundary flags */
    IASN(asgpars->lbc, IPNT_0); /* default left bc flag */
    IASN(asgpars->rbc, IPNT_0); /* default right bc flag */
        
    for (int i=0;i<model.g.dim;i++) {
      if (crank[i]==0) asgpars->lbc[i]=1;
      if (crank[i]==cdims[i]-1) asgpars->rbc[i]=1;
    }
    /************ END *************/
        
    /* decode half-order - with version 2.0, deprecated syntax
       "scheme_phys" etc. is dropped */
    asgpars->k = valparse<int>(pars,"order",1);
#ifdef IWAVE_VERBOSE
    fprintf(stream,"NOTE: initializing ASG with half-order = %d\n",asgpars->k);
#endif
    // MAXK set in sgcoeffs.h
    if (asgpars->k<1 || asgpars->k>MAXK) {
      RVLException e;
      e<<"ERROR: asgmodelinit\n";
      e<<"finite difference half-order requested = "<<asgpars->k<<" out of available range [1,"<<MAXK<<"]\n";
      throw e;
    }
        
    /* extract grid steps from grid */
    get_d(asgpars->dxs, model.g);
        
    /* initialize scaled Courant arrays, pml layers */
    asgpars->coeffs = (float **)usermalloc_(RARR_MAX_NDIM*sizeof(float*));
    for (int idim = 0; idim<RARR_MAX_NDIM; idim++)
      (asgpars->coeffs)[idim]=NULL;
    for (int idim = 0;idim < model.g.dim;idim ++) {
      if (asgpars->dxs[idim] <= 0.0) {
	RVLException e;
	e<<"ERROR: asg_modelinit\n";
	e<<"  bad input: wrong grid space step, dim="<<idim<<", step="<<asgpars->dxs[idim]<<"\n";
	throw e;
      }
      lam = model.tsind.dt / ((asgpars->dxs)[idim]);
      (asgpars->coeffs)[idim] = sgcoeffs(asgpars->k);
      for (int i=0;i<asgpars->k;i++) {
	(asgpars->coeffs)[idim][i]*=lam;
	//	fprintf(stderr,"c[%d][%d]=%g\n",idim,i,(asgpars->coeffs)[idim][i]);
      }
            
      ireal tmp;
            
      std::stringstream ln;
      ln<<"nl";
      ln<<idim+1;
      tmp = valparse<ireal>(pars,ln.str(),0.0f);
      model.nls[idim]=iwave_max(0,(int)((ceil)(tmp/((asgpars->dxs)[idim]))));
            
      std::stringstream rn;
      rn<<"nr";
      rn<<idim+1;
      tmp = valparse<ireal>(pars,rn.str(),0.0f);
      model.nrs[idim]=iwave_max(0,(int)((ceil)(tmp/((asgpars->dxs)[idim]))));
    }
        
    /* reserve copies of model attributes to pass to timestep -
       these should really be protected attributes of base class
       in which case they would not need to be copied!
    */
    asgpars->ndim = model.g.dim;
    asgpars->dt = model.tsind.dt;
    IASN(asgpars->nls,model.nls);
    IASN(asgpars->nrs,model.nrs);
        
    /* initialize pml amplitude */
    asgpars->amp = valparse<ireal>(pars,"pmlampl",1.0f);

    /* null initialization of pml damping arrays -
       will flag initialization in timestep function,
       called after grid initialization - can't do it
       here, as comp rarrs not yet set */
    for (int i=0;i<RARR_MAX_NDIM; i++) {
      (asgpars->ep)[i]=NULL;
      (asgpars->epp)[i]=NULL;
      (asgpars->ev)[i]=NULL;
      (asgpars->evp)[i]=NULL;
    }
        
    /* unset coefficient extension flag - set on first 
       pass to indicate coefficient definition in pml zone 
    */
    asgpars->extd = 0;

    /* initialize bound check params */
    asgpars->cmax = valparse<ireal>(pars,"cmax");
    asgpars->cmin = valparse<ireal>(pars,"cmin");

    /* test for internal model extension, identify axes,
       test for symmetry. Current rules:
       1. extended axes must be disposed symmetrically about 0
       2. extended axes must follow the time axis (index dim) in 
       the model grid, contiguously. Axis dim+1 extends axis 0; 
       axis dim+2 extends axis 1; etc. There can be at most dim-1
       extended axes.
    */
    
    IASN(asgpars->ihmin,IPNT_0);
    IASN(asgpars->ihmax,IPNT_0);    
    get_gs(asgpars->ihmin,model.g);
    get_ge(asgpars->ihmax,model.g);

#ifdef VERBOSE
    fprintf(asgpars->stream,"dim    = %d\n",model.g.dim);
    fprintf(asgpars->stream,"gdim   = %d\n",model.g.gdim);
    for (int i=0;i<RARR_MAX_NDIM; i++) {
      fprintf(asgpars->stream,"i=%d ihmin=%d ihmax=%d\n",i,asgpars->ihmin[i],asgpars->ihmax[i]);
      fprint_axis(asgpars->stream,model.g.axes[i]);
    }
#endif
    int count=0;
    // assumption of symmetry - can be dropped
    for (int i=model.g.dim; i<model.g.gdim; i++) {
      if (model.g.axes[i].id > EXTINT - 1) {
	if ((asgpars->ihmax[i] < 0) ||
	    (asgpars->ihmin[i] != -asgpars->ihmax[i])) {
	  RVLException err;
	  err<<"Error: for axis "<<i<<" identified as internal extended:\n";
	  err<<"  id = "<< model.g.axes[i].id <<"\n";
	  err<<"  either limits are not symmetric = ["
	     <<asgpars->ihmin[i]<<", "
	     <<asgpars->ihmax[i]<<"]\n";
	  err<<"  (symmetry is required!) or too many extended axes:\n";
	  err<<"  count = "<<count<<" dim = "<<model.g.dim<<" gdim = "
	   <<model.g.gdim<<"\n";
	  err<<" o = "<<model.g.axes[i].o <<" n = "<<model.g.axes[i].n << " d = "<<model.g.axes[i].d<<"\n";
	  throw err;
	}
	// redundant check that extended axis step is same as that of corresp
	// physical axis - for convenience don't do it
      }
      // if not internal, then internal loop limits equal to those
      // specifying unitialized axis in rarray
      else {
	asgpars->ihmax[i]=-1;
	asgpars->ihmin[i]=0;
      }
    }
         
    /* identify active fields */
    // cerr<<"resize\n";
    model.active.resize(D_SIZE);
    // cerr<<"assign active\n";
    model.active[D_BULK]="bulkmod";
    model.active[D_BUOY]="buoyancy";
    model.active[D_P0]="p0";
    model.active[D_V0]="v0";
    // cerr<<"dim>1\n";
    if (asgpars->ndim > 1) {
      model.active[D_P1]="p1";
      model.active[D_V1]="v1";
    }
    else {
      model.active[D_P1]="fruitcake";
      model.active[D_V1]="fruitcake";
    }
    // cerr<<"dim>2\n";
    if (asgpars->ndim > 2) {
      model.active[D_P2]="p2";
      model.active[D_V2]="v2";
    }
    else {
      model.active[D_P2]="fruitcake";
      model.active[D_V2]="fruitcake";
    }
    // cerr<<"specs\n";
    /* assign param object pointer */
    model.specs = (void*)asgpars;
        
    // cerr<<"return\n";
    return 0;
  }
  catch (RVLException & e) {
    e<<"\ncalled from asg_modelinit\n";
    throw e;
  }
}

/*----------------------------------------------------------------------------*/
void asg_modeldest(void ** fdpars) {
  /* destroy asgpars - all data members allocated on stack */
  ASG_TS_PARS * asgpars = (ASG_TS_PARS *)(*fdpars);
    
  //  for (int i=0; i<asgpars->ndim; i++) {
  for (int i=0; i<RARR_MAX_NDIM; i++) {
    if ((asgpars->ep)[i]) userfree_((asgpars->ep)[i]);
    if ((asgpars->epp)[i]) userfree_((asgpars->epp)[i]);
    if ((asgpars->ev)[i]) userfree_((asgpars->ev)[i]);
    if ((asgpars->evp)[i]) userfree_((asgpars->evp)[i]);
    if ((asgpars->coeffs)[i]) userfree_((asgpars->coeffs)[i]);
  }
  userfree_(asgpars->coeffs);
  userfree_(*fdpars);
}

/*----------------------------------------------------------------------------*/
int asg_timegrid(PARARRAY * pars,
                 FILE * stream,
                 grid const & g,
                 ireal & dt,
                 ireal & rhs) {
    
  try {
        
    ireal cmax = valparse<ireal>(*pars,"cmax");
    ireal cfl  = valparse<ireal>(*pars,"cfl");
        
    /* compute max stable step, optionally scaled by cfl from table */
    /* THIS BELONGS IN CHECK
       if ((err=sg_readgustime(dom,stream,&dtgus,g,cflgus,par))) {
       fprintf(stream,"\ncalled from sg_readtimegrid\n");
       return err;
       }
    */
        
    /* inline
       if ((err=sg_readcfltime(dom,stream,dt,g,cfl,cmax,par))) {
       fprintf(stream,"\ncalled from sg_readtimegrid\n");
       return err;
       }
    */
    ireal a = g.axes[0].d;
    for (int i = 1; i < g.dim; ++i ) a = iwave_min(a, g.axes[i].d);
        
    if ( (a < REAL_EPS) || (cfl < REAL_EPS) ) {
      fprintf(stream,"ERROR: sg_readcfltime - either min dx=%e or cfl=%e "
	      " too small\n", a, cfl);
      return E_BADINPUT;
    }
    dt = a*cfl/(cmax*sqrt((float)(g.dim)));
        
    /*  MOVE TO CHECK
	if (*dt > dtgus) {
	fprintf(stream,"ERROR: sg_readtimegrid - max_step option\n");
	fprintf(stream,"computed dt based on CFL criterion only = %e\n",*dt);
	fprintf(stream,"exceeds max stable step = %e\n",dtgus);
	return E_BADINPUT;
	}
    */
        
    /* however if dt is in the parameters, use it and throw the rest away */
    ireal dttmp = dt;
    dt = valparse<ireal>(*pars,"dt",dttmp);
        
    /* finally,... */
    rhs = dt;
        
    return 0;
  }
  catch (RVLException & e) {
    e<<"\ncallef from asg_timegrid\n";
    throw e;
  }
}

void asg_pmlaxis(int n0, int nl, int nr,
                 ireal amp, ireal dt, int gtype,
                 ireal ** ep, ireal ** epp) {
    
  if (*ep || *epp) return;
  nl=iwave_max(0,nl);
  nr=iwave_max(0,nr);
    
  //    cerr<<"pmlaxis: gtype="<<gtype<<" n0="<<n0<<" nr="<<nr<<" nl="<<nl<<"\n";
    
  *ep  = (ireal*)usermalloc_(sizeof(ireal)*n0);
  *epp = (ireal*)usermalloc_(sizeof(ireal)*n0);
  for (int i=0;i<nl;i++) {
    ireal p = (ireal(nl-i)-0.5*gtype)/ireal(nl);
    p = amp*fabs(p*p*p);
    (*ep)[i] = REAL_ONE - 0.5*dt*p;
    (*epp)[i] = REAL_ONE/(REAL_ONE + 0.5*dt*p);
  }
  for (int i=nl; i<n0-nr;i++) {
    (*ep)[i] = REAL_ONE;
    (*epp)[i]= REAL_ONE;
  }
  for (int i=n0-nr; i<n0;i++) {
    ireal p =  (ireal(i-(n0-nr))+1.0-0.5*gtype)/ireal(nr);
    p = amp*fabs(p*p*p);
    (*ep)[i] = REAL_ONE - 0.5*dt*p;
    (*epp)[i] = REAL_ONE/(REAL_ONE + 0.5*dt*p);
    //        cerr<<"i="<<i<<" p="<<p<<" ep="<<(*ep)[i]<<" epp="<<(*epp)[i]<<endl;
  }
}

void asg_timestep(std::vector<RDOM *> dom,
                  bool fwd,
                  int iv,
                  void* fdpars) {
    
  IPNT np;
  IPNT nv;

  //  cerr<<"step\n";
  
  /* fd parameter struct */
  ASG_TS_PARS * asgpars = (ASG_TS_PARS *) fdpars;
    
#ifdef VERBOSE
  fprintf(asgpars->stream,"asg_timestep: fwd=%d iv=%d\n",fwd,iv);
  fflush(asgpars->stream);
#endif
    
  /* field indices - from FIELDS struct */
  IPNT i_p;
  i_p[0] = D_P0;
  i_p[1] = D_P1;
  i_p[2] = D_P2;
  IPNT i_v;
  i_v[0] = D_V0;
  i_v[1] = D_V1;
  i_v[2] = D_V2;
    
  /* fill in pml arrays; these are no-ops after the first call */
  register ireal * restrict ep[RARR_MAX_NDIM];
  register ireal * restrict epp[RARR_MAX_NDIM];
  register ireal * restrict ev[RARR_MAX_NDIM];
  register ireal * restrict evp[RARR_MAX_NDIM];

  // local pml widths for this domain
  IPNT nrsloc;
  IPNT nlsloc;
  IASN(nlsloc,IPNT_0);
  IASN(nrsloc,IPNT_0);
  
  for (int idim=0;idim<asgpars->ndim;idim++) {
    IPNT gsa;
    IPNT gea;
    rd_a_size(dom[0], i_p[idim], np);
    rd_a_size(dom[0], i_v[idim], nv);
    rd_a_gse(dom[0], i_p[idim],gsa,gea);
    // pml params
    // first check if this is a sticky-outy part
    if (asgpars->lbc[idim]) nlsloc[idim]=asgpars->nls[idim];
    if (asgpars->rbc[idim]) nrsloc[idim]=asgpars->nrs[idim];
    // then set up axis
    asg_pmlaxis(np[idim], nlsloc[idim], nrsloc[idim],
		asgpars->amp, asgpars->dt, 0,
		&((asgpars->ep)[idim]), &((asgpars->epp)[idim]));
    // re-index
    ep[idim]=(asgpars->ep)[idim] - gsa[idim];
    epp[idim]=(asgpars->epp)[idim] - gsa[idim];
    // do likewise for v
    rd_a_gse(dom[0], i_v[idim],gsa,gea);
    asg_pmlaxis(nv[idim], nlsloc[idim], nrsloc[idim],
		asgpars->amp, asgpars->dt, 1,
		&((asgpars->ev)[idim]), &((asgpars->evp)[idim]));
    ev[idim]=(asgpars->ev)[idim] - gsa[idim];
    evp[idim]=(asgpars->evp)[idim] - gsa[idim];
  }
    
  /* sizes of computational domain - P0 same as P1, P2 */
  IPNT gsc_p;
  IPNT gec_p;
  int gsc_v[RARR_MAX_NDIM][RARR_MAX_NDIM];
  int gec_v[RARR_MAX_NDIM][RARR_MAX_NDIM];
  rd_gse(dom[0], i_p[0], gsc_p, gec_p);
  rd_gse(dom[0], i_v[0], gsc_v[0], gec_v[0]);
  rd_gse(dom[0], i_v[1], gsc_v[1], gec_v[1]);
  rd_gse(dom[0], i_v[2], gsc_v[2], gec_v[2]);
    
  int ndim;
  ra_ndim(&(dom[0]->_s[i_p[0]]),&ndim);

  // extend coefficient arrays into pml zones
  // only for dom[0] - perturbations should live only
  // in non-pml zones!
  if (!(asgpars->extd)) {
    // based on principle that coeffs are first batch of
    // rarrays - D_BULK=0, D_BUOY=1, could be generalized
    //    fprintf(stderr,"extend: nlsloc[0]=%d nrsloc[0]=%d nlsloc[1]=%d nrsloc[1]=%d\n",nlsloc[0],nrsloc[0],nlsloc[1],nrsloc[1]);
    IPNT gs;
    IPNT ge;
    int ncoeff = 2;
    for (int icoeff=0;icoeff<ncoeff;icoeff++) {
      rd_a_gse(dom[0], icoeff, gs, ge);
      float ** f2 = (dom[0]->_s)[icoeff]._s2;
      float *** f3= (dom[0]->_s)[icoeff]._s3;
      //      if (ndim == 2) {
      if ((f2) && (!f3)) {
	register ireal ** restrict f = (dom[0]->_s)[icoeff]._s2;
	// extend in dim 1
	for (int i0=gs[0]+nlsloc[0]; i0<=ge[0]-nrsloc[0]; i0++) {
	  for (int i1=gs[1]; i1<=gs[1]+nlsloc[1]; i1++)
	    f[i1][i0]=f[gs[1]+nlsloc[1]][i0];
	  for (int i1=ge[1]-nrsloc[1]; i1<=ge[1]; i1++)
	    f[i1][i0]=f[ge[1]-nrsloc[1]][i0];
	}
	// extend in dim 0
	for (int i1=gs[1]; i1<=ge[1]; i1++) {
	  for (int i0=gs[0]; i0<=gs[0]+nlsloc[0]-1; i0++)
	    f[i1][i0]=f[i1][gs[0]+nlsloc[0]];
	  for (int i0=ge[0]-nrsloc[0]+1; i0<=ge[0]; i0++)
	    f[i1][i0]=f[i1][ge[0]-nrsloc[0]];
	}
      }
      //      else if (ndim == 3) {
      else if (f3) {
	register ireal *** restrict f = (dom[0]->_s)[icoeff]._s3;
	// extend in dim 2
	for (int i0=gs[0]+nlsloc[0]; i0<=ge[0]-nrsloc[0]; i0++) {
	  for (int i1=gs[1]+nlsloc[1]; i1<=ge[1]-nrsloc[1]; i1++) {
	    for (int i2=gs[2]; i2<=gs[2]+nlsloc[2]; i2++)
	      f[i2][i1][i0]=f[gs[2]+nlsloc[2]][i1][i0];
	    for (int i2=ge[2]-nrsloc[2]; i2<=ge[2]; i2++)
	      f[i2][i1][i0]=f[ge[2]-nrsloc[2]][i1][i0];
	  }
	}
	// extend in dim 1
	for (int i2=gs[2]; i2<=ge[2]; i2++) {
	  for (int i0=gs[0]+nlsloc[0]; i0<=ge[0]-nrsloc[0]; i0++) {
	    for (int i1=gs[1]; i1<=gs[1]+nlsloc[1]-1; i1++)
	      f[i2][i1][i0]=f[i2][gs[1]+nlsloc[1]][i0];
	    for (int i1=ge[1]-nrsloc[1]+1; i1<=ge[1]; i1++)
	      f[i2][i1][i0]=f[i2][ge[1]-nrsloc[1]][i0];
	  }
	}
	// extend in dim 0
	for (int i2=gs[2]; i2<=ge[2]; i2++) {
	  for (int i1=gs[1]; i1<=ge[1]; i1++) {
	    for (int i0=gs[0]; i0<=gs[0]+nlsloc[0]-1; i0++)
	      f[i2][i1][i0]=f[i2][i1][gs[0]+nlsloc[0]];
	    for (int i0=ge[0]-nrsloc[0]+1; i0<=ge[0]; i0++)
	      f[i2][i1][i0]=f[i2][i1][ge[0]-nrsloc[0]];
	  }
	}
      }
      else {
	RVLException e;
	e<<"Error: asg_timestep(): extension defined only for dims 2 and 3";
	throw e;
      }
    }

    // extension done, set flag
    asgpars->extd = 1;
  }

  // find loop limits for extended axes if any
  // execute if first possible - any internal extended axis must
  // be symmetric and match axis lengths from model grid
  // and any other extended axis should have zero loop limits
  // 2017.11.10 THIS CHECK IS REDUNDANT WWS
  /*
  IPNT gs;
  IPNT ge;
  rd_a_gse(dom[0], D_BULK, gs, ge);
  for (int i=ndim; i<RARR_MAX_NDIM; i++) {
    if (ge[i]-gs[i]+1 != asgpars->ihmax[i]-asgpars->ihmin[i]+1) {
      RVLException e;
      e<<"Error: asg_step - extended bulkmod\n";
      e<<"  illegal limits on extended axis "<<i<<": s="
       <<gs[i]<<" e="<<ge[i]<<"\n";
      e<<"  model grid limits are "<<asgpars->ihmin[i]<<", "<<asgpars->ihmax[i]<<"\n";
      throw e;
    }
  }
  */
  
  int n  = dom.size();

  // default values for extended axis step and limit
  RPNT dha;
  RASN(dha,asgpars->dxs);
  IPNT ihmax;
  IASN(ihmax,asgpars->ihmax);
  
  if (ndim==2)
      
  // deriv = 0
  if (n==1) {
    //    if (fwd == true) {
    if (ndim == 2) {

      register ireal ** restrict bulk2 = (dom[0]->_s)[D_BULK]._s2;
      register ireal ** restrict buoy2 = (dom[0]->_s)[D_BUOY]._s2;
      register ireal *** restrict bulk3 = (dom[0]->_s)[D_BULK]._s3;
      if (!bulk3) {
	// create new 3D array if necessary, also set loop limit on
	// extended axis appropriately
	bulk3 = (float ***)usermalloc_(sizeof(float **));
	bulk3[0] = bulk2;
	dha[2]=1.0f;
	ihmax[2]=0;
      }
      register ireal ** restrict p02 = (dom[0]->_s)[i_p[0]]._s2;
      register ireal ** restrict p12 = (dom[0]->_s)[i_p[1]]._s2;
      register ireal ** restrict v02 = (dom[0]->_s)[i_v[0]]._s2;
      register ireal ** restrict v12 = (dom[0]->_s)[i_v[1]]._s2;
                
#ifdef VERBOSE
      fprintf(asgpars->stream,"asg_timestep: deriv=0 extd=%d\n",
	      asgpars->next[2]);
      fflush(asgpars->stream);
      float pt=0.0;
#endif
      /* approach with div, grad buffers */
      // aligned vector lengths
      int align_words = USERMALLOC_ALIGN_BYTES/sizeof(float);
      int p_alloc = align_words*(1+((gec_p[0]-gsc_p[0]+1)/align_words));
      // allocate sdiv
      float * sdiv_alloc = (float *)usermalloc_((gec_p[1]-gsc_p[1]+1)*p_alloc*sizeof(float));
      float ** sdiv = (float **)usermalloc_((gec_p[1]-gsc_p[1]+1)*sizeof(float *));
      for (int idim=0; idim<gec_p[1]-gsc_p[1]+1; idim++) {
	sdiv[idim] = &(sdiv_alloc[idim*p_alloc]);
      }

      // allocate sdivb
      float * sdivb_alloc = (float *)usermalloc_((gec_p[1]-gsc_p[1]+1)*p_alloc*sizeof(float));
      float ** sdivb = (float **)usermalloc_((gec_p[1]-gsc_p[1]+1)*sizeof(float *));
      for (int idim=0; idim<gec_p[1]-gsc_p[1]+1; idim++) {
	sdivb[idim] = &(sdivb_alloc[idim*p_alloc]);
      }
      
      // allocate gradp
      float * gradp0 = (float *)usermalloc_((gec_v[0][0]-gsc_v[0][0]+1)*sizeof(float));
      float * gradp1 = (float *)usermalloc_((gec_v[1][0]-gsc_v[1][0]+1)*sizeof(float));
      
      for (int i1=gsc_p[1]+asgpars->nls[1]; i1<=gec_p[1]-asgpars->nrs[1];i1++) {
	for (int i0=gsc_p[0]+asgpars->nls[0]; i0<=gec_p[0]-asgpars->nrs[0];i0++) { 
	  p12[i1][i0]=p02[i1][i0];
#ifdef VERBOSE
	  pt += p02[i1][i0]*p02[i1][i0];
#endif
	}
      }
#ifdef VERBOSE
      fprintf(asgpars->stream,"asg_timestep -> asg2d\n");
      fprintf(asgpars->stream,"input p0 l2 = %g\n",pt);
      fflush(asgpars->stream);
#endif
      
      asg2dtap(bulk3,
	       buoy2,
	       p02,p12,
	       v02,v12,
	       ep, epp,
	       ev,evp,
	       sdiv,sdivb,
	       gradp0,gradp1,
	       gsc_p,gec_p,
	       &(gsc_v[0][0]),&(gec_v[0][0]),
	       &(gsc_v[1][0]),&(gec_v[1][0]),
	       asgpars->lbc,asgpars->rbc,
	       dha[2],ihmax[2],
	       asgpars->k,asgpars->coeffs,
	       asgpars->stream);

#ifdef VERBOSE
      pt=0.0f;
      for (int i1=gsc_p[1]+asgpars->nls[1]; i1<=gec_p[1]-asgpars->nrs[1];i1++) {
	for (int i0=gsc_p[0]+asgpars->nls[0]; i0<=gec_p[0]-asgpars->nrs[0];i0++) {
	  pt += p02[i1][i0]*p02[i1][i0];
	}
      }
      fprintf(asgpars->stream,"output p0 l2 = %g\n",pt);
      fflush(asgpars->stream);
#endif

      userfree_(sdiv);
      userfree_(sdiv_alloc);
      userfree_(sdivb);
      userfree_(sdivb_alloc);
      userfree_(gradp0);
      userfree_(gradp1);
      
#ifdef VERBOSE
      fprintf(asgpars->stream,"asg_timestep <- asg2d\n");
      fflush(asgpars->stream);
#endif
    }
    if (ndim == 3) {
      RVLException e;
      e<<"Error: asg_timestep(). 3D modeling not implemented yet.\n";

      /*      
      register ireal *** restrict bulk3 = (dom[0]->_s)[D_BULK]._s3;
      register ireal *** restrict buoy3 = (dom[0]->_s)[D_BUOY]._s3;
      register ireal *** restrict p03 = (dom[0]->_s)[i_p[0]]._s3;
      register ireal *** restrict p13 = (dom[0]->_s)[i_p[1]]._s3;
      register ireal *** restrict p23 = (dom[0]->_s)[i_p[2]]._s3;
      register ireal *** restrict v03 = (dom[0]->_s)[i_v[0]]._s3;
      register ireal *** restrict v13 = (dom[0]->_s)[i_v[1]]._s3;
      register ireal *** restrict v23 = (dom[0]->_s)[i_v[2]]._s3;

      ireal * sdiv_alloc = (ireal *)usermalloc_((gec_p[0]-gsc_p[0]+1)*sizeof(ireal));
      ireal * sdiv = &(sdiv_alloc[-gsc_p[0]]);
      for (int i2=gsc_p[2]+asgpars->nls[2]; i2<=gec_p[2]-asgpars->nrs[2];i2++) {
	for (int i1=gsc_p[1]+asgpars->nls[1]; i1<=gec_p[1]-asgpars->nrs[1];i1++) {
	  for (int i0=gsc_p[0]+asgpars->nls[0]; i0<=gec_p[0]-asgpars->nrs[0];i0++) {
	    p13[i2][i1][i0]=p03[i2][i1][i0];
	    p23[i2][i1][i0]=p03[i2][i1][i0];
	  }
	}
      }
      asg_pstep3d(bulk3,
		  p03,p13,p23,
		  v03,v13,v23,
		  ep, epp,
		  sdiv,
		  gsc_p,gec_p,
		  asgpars->lbc,asgpars->rbc,
		  asgpars->k,asgpars->coeffs);
      userfree_(sdiv_alloc);

      ireal * gradp[RARR_MAX_NDIM];
      ireal * gradp_alloc[RARR_MAX_NDIM];
      gradp_alloc[0] = (ireal *)usermalloc_((gec_v[0][0]-gsc_v[0][0]+1)*sizeof(ireal));
      gradp_alloc[1] = (ireal *)usermalloc_((gec_v[1][0]-gsc_v[1][0]+1)*sizeof(ireal));
      gradp_alloc[2] = (ireal *)usermalloc_((gec_v[2][0]-gsc_v[2][0]+1)*sizeof(ireal));
      gradp[0]=&(gradp_alloc[0][-gsc_v[0][0]]);
      gradp[1]=&(gradp_alloc[1][-gsc_v[1][0]]);
      gradp[2]=&(gradp_alloc[2][-gsc_v[2][0]]);
                    
      asg_vstep3d(buoy3,
		  p03,p13,p23,
		  v03,v13,v23,
		  ev, evp,
		  gradp,
		  &(gsc_v[0][0]),&(gec_v[0][0]),
		  &(gsc_v[1][0]),&(gec_v[1][0]),
		  &(gsc_v[2][0]),&(gec_v[2][0]),
		  asgpars->lbc,asgpars->rbc,
		  asgpars->k,asgpars->coeffs);
      userfree_(gradp_alloc[0]);
      userfree_(gradp_alloc[1]);
      userfree_(gradp_alloc[2]);
      */

    }
  }  
  else if (n==2) {
    //    cerr<<"step n=2 (deriv=1)\n";
        //    if (fwd == true) {
    if (ndim == 2) {
      register ireal ** restrict bulk2 = (dom[0]->_s)[D_BULK]._s2;
      register ireal ** restrict buoy2 = (dom[0]->_s)[D_BUOY]._s2;
      register ireal *** restrict bulk3 = (dom[0]->_s)[D_BULK]._s3;
      register ireal ** restrict p02 = (dom[0]->_s)[i_p[0]]._s2;
      register ireal ** restrict p12 = (dom[0]->_s)[i_p[1]]._s2;
      register ireal ** restrict v02 = (dom[0]->_s)[i_v[0]]._s2;
      register ireal ** restrict v12 = (dom[0]->_s)[i_v[1]]._s2;
      register ireal ** restrict bulk2d = (dom[1]->_s)[D_BULK]._s2;
      register ireal ** restrict buoy2d = (dom[1]->_s)[D_BUOY]._s2;
      register ireal *** restrict bulk3d = (dom[1]->_s)[D_BULK]._s3;
      register ireal ** restrict p02d = (dom[1]->_s)[i_p[0]]._s2;
      register ireal ** restrict p12d = (dom[1]->_s)[i_p[1]]._s2;
      register ireal ** restrict v02d = (dom[1]->_s)[i_v[0]]._s2;
      register ireal ** restrict v12d = (dom[1]->_s)[i_v[1]]._s2;

      // this version requires bulk and bulkd to be commensurable
      // ie representing vectors in same space
      if ((!bulk3) && (!bulk3d)) {
	bulk3 = (float ***)usermalloc_(sizeof(float **));
	bulk3[0] = bulk2;
	bulk3d = (float ***)usermalloc_(sizeof(float **));
	bulk3d[0] = bulk2d;
	dha[2]=1.0f;
	ihmax[2]=0;
      }

      /* approach with div, grad buffers */
      // aligned vector lengths
      int align_words = USERMALLOC_ALIGN_BYTES/sizeof(float);
      int p_alloc = align_words*(1+((gec_p[0]-gsc_p[0]+1)/align_words));

      // allocate sdiv
      float * sdiv_alloc = (float *)usermalloc_((gec_p[1]-gsc_p[1]+1)*p_alloc*sizeof(float));
      float ** sdiv = (float **)usermalloc_((gec_p[1]-gsc_p[1]+1)*sizeof(float *));
      for (int idim=0; idim<gec_p[1]-gsc_p[1]+1; idim++) {
	sdiv[idim] = &(sdiv_alloc[idim*p_alloc]);
      }

      // allocate sdivd
      float * sdivd_alloc = (float *)usermalloc_((gec_p[1]-gsc_p[1]+1)*p_alloc*sizeof(float));
      float ** sdivd = (float **)usermalloc_((gec_p[1]-gsc_p[1]+1)*sizeof(float *));
      for (int idim=0; idim<gec_p[1]-gsc_p[1]+1; idim++) {
	sdivd[idim] = &(sdivd_alloc[idim*p_alloc]);
      }

      // allocate sdivb
      float * sdivb_alloc = (float *)usermalloc_((gec_p[1]-gsc_p[1]+1)*p_alloc*sizeof(float));
      float ** sdivb = (float **)usermalloc_((gec_p[1]-gsc_p[1]+1)*sizeof(float *));
      for (int idim=0; idim<gec_p[1]-gsc_p[1]+1; idim++) {
	sdivb[idim] = &(sdivb_alloc[idim*p_alloc]);
      }
      
      // allocate sdivbd
      float * sdivbd_alloc = (float *)usermalloc_((gec_p[1]-gsc_p[1]+1)*p_alloc*sizeof(float));
      float ** sdivbd = (float **)usermalloc_((gec_p[1]-gsc_p[1]+1)*sizeof(float *));
      for (int idim=0; idim<gec_p[1]-gsc_p[1]+1; idim++) {
	sdivbd[idim] = &(sdivbd_alloc[idim*p_alloc]);
      }
      
      // allocate gradp
      float * gradp0 = (float *)usermalloc_((gec_v[0][0]-gsc_v[0][0]+1)*sizeof(float));
      float * gradp1 = (float *)usermalloc_((gec_v[1][0]-gsc_v[1][0]+1)*sizeof(float));

      // allocate gradp
      float * gradp0d = (float *)usermalloc_((gec_v[0][0]-gsc_v[0][0]+1)*sizeof(float));
      float * gradp1d = (float *)usermalloc_((gec_v[1][0]-gsc_v[1][0]+1)*sizeof(float));

      for (int i1=gsc_p[1]+asgpars->nls[1]; i1<=gec_p[1]-asgpars->nrs[1]; i1++) {
	for (int i0=gsc_p[0]+asgpars->nls[0]; i0<=gec_p[0]-asgpars->nrs[0]; i0++) {
	  p12[i1][i0]=p02[i1][i0];
	}
      }

      if (fwd) {

	asg2dtap_d(bulk3, bulk3d, buoy2, buoy2d, p02, p02d, p12, p12d,
		   v02, v02d, v12, v12d, ep, epp, ev, evp,
		   sdiv, sdivd, sdivb, sdivbd, 
		   gradp0, gradp0d, gradp1, gradp1d,
		   gsc_p,gec_p,
		   &(gsc_v[0][0]),&(gec_v[0][0]),
		   &(gsc_v[1][0]),&(gec_v[1][0]),
		   asgpars->lbc,asgpars->rbc,
		   dha[2],ihmax[2],		   
		   asgpars->k,asgpars->coeffs,
		   asgpars->stream);

      }
      else {
	
	asg2dtap_b(bulk3, bulk3d, buoy2, buoy2d, p02, p02d, p12, p12d,
		   v02, v02d, v12, v12d, ep, epp, ev, evp,
		   sdiv, sdivd, sdivb, sdivbd, 
		   gradp0, gradp0d, gradp1, gradp1d,
		   gsc_p,gec_p,
		   &(gsc_v[0][0]),&(gec_v[0][0]),
		   &(gsc_v[1][0]),&(gec_v[1][0]),
		   asgpars->lbc,asgpars->rbc,
		   dha[2],ihmax[2],		   
		   asgpars->k,asgpars->coeffs,
		   asgpars->stream);
      }

      userfree_(sdiv);		    
      userfree_(sdiv_alloc);
      userfree_(sdivb);
      userfree_(sdivb_alloc);
      userfree_(gradp0);
      userfree_(gradp1);
      
      userfree_(sdivd);
      userfree_(sdivd_alloc);
      userfree_(sdivbd);
      userfree_(sdivbd_alloc);
      userfree_(gradp0d);
      userfree_(gradp1d);
      
   }

    if (ndim == 3) {
      RVLException e;
      e<<"Error: asg_timestep(). 3D derivatives not implemented yet.\n";
    }
  }

  else {
    RVLException e;
    //    e<<"Error: asg_timestep(). Only 0th, 1st derivatives are implemented.\n";
    e<<"Error: asg_timestep(). Only 0th derivatives are implemented.\n";    
    throw e;
  }

  //cerr<<"exit fd_asg::asg_timestep\n";
}

/**
 * all stencils are of the same type: for each axis idim,
 * - source v_idim, target p, offset range -k,...,k-1 on axis idim
 * - source p, target v_idim, offset range --k+1,...k on axis idim
 * so there are 2*ndim masks in total, each mask storing 2k offset vectors
 */
int asg_create_sten(void * specs,
                    FILE * stream,
                    int ndim,
                    IPNT gtype[RDOM_MAX_NARR],
                    STENCIL * sten) {
    
  // cerr<<"begin stencil\n";
    
  ASG_TS_PARS * asgpars = (ASG_TS_PARS *)(specs);
  STENCIL_MASK mask;/* workspace */
  int err = 0;
    
  /* number of masks */
  int nmask = (asgpars->ndim+2)*(asgpars->ndim)+1;
    
  /* initialize STENCIL to null stencil */
  sten_setnull(sten);
    
  // cerr<<"asg_stencil_create stencil nmask="<<nmask<<endl;
  /* nontrivial STENCIL initialization */
  if ((err = sten_create(sten,nmask))) {
    fprintf(stream,"ERROR: asg_create_sten - failed to create stencil\n");
    return err;
  }
    
  /* field indices - from FIELDS struct */
  IPNT ip;
    
  ip[0] = D_P0;
  ip[1] = D_P1;
  ip[2] = D_P2;
  /*
    ip[0] = 2;
    ip[1] = 2;
    ip[2] = 2;
  */
  IPNT iv;
  iv[0] = D_V0;
  iv[1] = D_V1;
  iv[2] = D_V2;
    
  /* every mask stores 2k offsets */
  int len = 2*asgpars->k;
    
  /* loop over axes */
    
  int imask=0; // mask counter
    
  for (int idim=0; idim<asgpars->ndim; idim++) {
        
    for (int jdim=0; jdim<asgpars->ndim; jdim++) {
      // cerr<<"asg_stencil_create imask="<<imask<<" idim="<<idim<<" jdim="<<jdim<<" mask workspace\n";
      /* source = v[idim], target = p[jdim] */
      if ((err=mask_create(&mask, iv[idim], ip[jdim], len))) {
	fprintf(stream,"ERROR: asg_create_sten from mask_create\n");
	sten_destroy(sten);
	return err;
      }
      // cerr<<"asg_stencil_create: idim="<<idim<<endl;
      for (int i=0; i<len; i++) {
	IPNT offs;
	IASN(offs,IPNT_0);
	offs[idim]=i-asgpars->k;
	//offs[idim]=i-asgpars->k+1;
	if ((err=mask_set(&mask,i,offs))) {
	  fprintf(stream,"ERROR: asg_create_sten from mask_set\n");
	  sten_destroy(sten);
	  return err;
	}
      }
      /* insert mask */
      if ((err=sten_set(sten,imask, &mask))) {
	fprintf(stream,"ERROR: asg_create_sten from sten_set\n");
	sten_destroy(sten);
	return err;
      }
            
      imask++;
            
    }
        
    /* source = p, target = v[idim] */
    // cerr<<"asg_stencil_create mask workspace\n";
    // cerr<<"asg_stencil_create imask="<<imask<<" idim="<<idim<<" mask workspace\n";
    if ((err=mask_create(&mask, ip[idim], iv[idim], len))) {
      fprintf(stream,"ERROR: asg_create_sten from mask_create\n");
      sten_destroy(sten);
      return err;
    }
    for (int i=0; i<len; i++) {
      IPNT offs;
      IASN(offs,IPNT_0);
      offs[idim]=i-asgpars->k+1;
      //offs[idim]=i-asgpars->k;
      if ((err=mask_set(&mask,i,offs))) {
	fprintf(stream,"ERROR: asg_create_sten from mask_set\n");
	sten_destroy(sten);
	return err;
      }
    }
    /* insert mask */
    if ((err=sten_set(sten,imask, &mask))) {
      fprintf(stream,"ERROR: asg_create_sten from sten_set\n");
      sten_destroy(sten);
      return err;
    }
    imask++;
  }
  /* end loop over axes */

  //  [[[add mask buoy(1) -> v0(5), v1(6), v2(7)
  for (int idim = 0; idim < asgpars->ndim; idim++){
        
    if ((err=mask_create(&mask, D_BUOY, iv[idim], 2))) {
      fprintf(stream,"ERROR: asg_create_sten from mask_create\n");
      sten_destroy(sten);
      return err;
    }
        
    IPNT offs;
    IASN(offs,IPNT_0);
    offs[idim]=0;
    if ((err=mask_set(&mask,0,offs))) {
      fprintf(stream,"ERROR: asg_create_sten from mask_set\n");
      sten_destroy(sten);
      return err;
    }
    IASN(offs,IPNT_0);
    offs[idim]=1;
    if ((err=mask_set(&mask,1,offs))) {
      fprintf(stream,"ERROR: asg_create_sten from mask_set\n");
      sten_destroy(sten);
      return err;
    }
        
        
    if ((err=sten_set(sten,imask, &mask))) {
      fprintf(stream,"ERROR: asg_create_sten from sten_set\n");
      sten_destroy(sten);
      return err;
    }
        
    imask++;
  }
  //  end]]]
    
    
  // cerr<<"end stencil\n";
  return 0;
}

void asg_check(RDOM * dom,
               void * specs,
               FILE * stream) {}

// no-op for now
void asg_loop_refine(int const * gmin, int const * gmax,
                     float tmax, int input,
                     FILE * stream, void * specs) { }

