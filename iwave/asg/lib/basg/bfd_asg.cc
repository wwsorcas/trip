#include "basg.hh"

//#define VERBOSE

using RVL::valparse;
using RVL::RVLException;


enum{
    
    D_SIZE=16,
    D_BULK=0,
    D_BUOY=1,
    D_P0=2,
    D_P1=3,
    D_P2=4,
    D_V0=5,
    D_V1=6,
    D_V2=7,
    D_DBULK=8,
    D_DBUOY=9,
    D_DP0=10,
    D_DP1=11,
    D_DP2=12,
    D_DV0=13,
    D_DV1=14,
    D_DV2=15,
    
};


// max FD half-order
#define MAXK 12

extern float * sgcoeffs(int k);

extern "C" void basg2d_d(float ** restrict bulk, float *** restrict bulkd, float ** restrict buoy, float ** restrict buoyd, float ** restrict p0, float **  restrict p0d, float ** restrict p1, float ** restrict p1d, float ** restrict v0, float ** restrict v0d, float ** restrict v1, float ** restrict v1d, float ** restrict ep, float ** restrict epp, float ** restrict ev, float ** restrict evp, float ** restrict sdiv, float ** restrict sdivd, float ** restrict sdivb, float ** restrict sdivbd, float *gradp0, float *gradp0d, float *gradp1, float *gradp1d, int *gsc_p, int *gec_p, int *gsc_v0, int *gec_v0, int *gsc_v1, int *gec_v1, int *lbc, int *rbc, int ihmax, int maxoff, float ** restrict c, FILE * stream);

extern "C" void basg2d_b(float ** restrict bulk, float *** restrict bulkb, float ** 
        restrict buoy, float ** restrict buoyb, float ** restrict p0, float **
        restrict p0b, float ** restrict p1, float ** restrict p1b, float ** 
        restrict v0, float ** restrict v0b, float ** restrict v1, float ** 
        restrict v1b, float ** restrict ep, float ** restrict epp, float ** 
        restrict ev, float ** restrict evp, float ** restrict sdiv, float ** 
        restrict sdivb0, float ** restrict sdivb, float ** restrict sdivbb, 
        float *gradp0, float *gradp0b, float *gradp1, float *gradp1b, int *
        gsc_p, int *gec_p, int *gsc_v0, int *gec_v0, int *gsc_v1, int *gec_v1,
        int *lbc, int *rbc, int ihmax, int maxoff, float ** restrict c, FILE *
		stream);

int basg_modelinit(PARARRAY pars,
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
        
    RPNT dxs;    /* grid steps */
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
    get_d(dxs, model.g);
        
    /* initialize scaled Courant arrays, pml layers */
    asgpars->coeffs = (float **)usermalloc_(RARR_MAX_NDIM*sizeof(float*));
    for (int idim = 0; idim<RARR_MAX_NDIM; idim++)
      (asgpars->coeffs)[idim]=NULL;
    for (int idim = 0;idim < model.g.dim;idim ++) {
      if (dxs[idim] <= 0.0) {
	RVLException e;
	e<<"ERROR: asg_modelinit\n";
	e<<"  bad input: wrong grid space step, dim="<<idim<<", step="<<dxs[idim]<<"\n";
	throw e;
      }
      lam = model.tsind.dt / dxs[idim];
      (asgpars->coeffs)[idim] = sgcoeffs(asgpars->k);
      for (int i=0;i<asgpars->k;i++) (asgpars->coeffs)[idim][i]*=lam;
            
      ireal tmp;
            
      std::stringstream ln;
      ln<<"nl";
      ln<<idim+1;
      tmp = valparse<ireal>(pars,ln.str(),0.0f);
      model.nls[idim]=iwave_max(0,(int)((ceil)(tmp/dxs[idim])));
            
      std::stringstream rn;
      rn<<"nr";
      rn<<idim+1;
      tmp = valparse<ireal>(pars,rn.str(),0.0f);
      model.nrs[idim]=iwave_max(0,(int)((ceil)(tmp/dxs[idim])));
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

    fprintf(asgpars->stream,"dim    = %d\n",model.g.dim);
    fprintf(asgpars->stream,"gdim   = %d\n",model.g.gdim);
    for (int i=0;i<RARR_MAX_NDIM; i++) {
      fprintf(asgpars->stream,"i=%d ihmin=%d ihmax=%d\n",i,asgpars->ihmin[i],asgpars->ihmax[i]);
      fprint_axis(asgpars->stream,model.g.axes[i]);
    }
    int count=0;
    // assumption of symmetry - can be dropped
    // NOTE THIS IS ACTUALLY VOID AS MODEL GRID DOES NOT HAVE EXTENDED AXES
    // under current design, there is no way to identify extended axes, as they do
    // not occur in the model (i.e. bulkmod, in this case) grid. new iwave needs a
    // better system to relate the entire array of grid definitions

    // so remove asgpars->ihmax, infer extended axes dynamically in the kernel driver
    
    /*
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
	  throw err;
	}
      }
      // if not internal, then internal loop limits equal to those
      // specifying unitialized axis in rarray
      else {
	asgpars->ihmax[i]=-1;
	asgpars->ihmin[i]=0;
      }
    }
    */   
    /* identify active fields */
    // cerr<<"resize\n";
    model.active.resize(D_SIZE);
    // cerr<<"assign active\n";
    model.active[D_BULK]="bulkmod";
    model.active[D_BUOY]="buoyancy";
    model.active[D_DBULK]="dbulkmod";
    model.active[D_DBUOY]="dbuoyancy";
    model.active[D_P0]="p0";
    model.active[D_V0]="v0";
    model.active[D_DP0]="dp0";
    model.active[D_DV0]="dv0";
    // cerr<<"dim>1\n";
    if (asgpars->ndim > 1) {
      model.active[D_P1]="p1";
      model.active[D_V1]="v1";
      model.active[D_DP1]="dp1";
      model.active[D_DV1]="dv1";
    }
    else {
      model.active[D_P1]="fruitcake";
      model.active[D_V1]="fruitcake";
      model.active[D_DP1]="fruitcake";
      model.active[D_DV1]="fruitcake";
    }
    // cerr<<"dim>2\n";
    if (asgpars->ndim > 2) {
      model.active[D_P2]="p2";
      model.active[D_V2]="v2";
      model.active[D_DP2]="dp2";
      model.active[D_DV2]="dv2";
    }
    else {
      model.active[D_P2]="fruitcake";
      model.active[D_V2]="fruitcake";
      model.active[D_DP2]="fruitcake";
      model.active[D_DV2]="fruitcake";
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


void basg_timestep(std::vector<RDOM *> dom,
                  bool fwd,
                  int iv,
                  void* fdpars) {
    
  IPNT np;
  IPNT nv;

  // cerr<<"step\n";
  
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
  IPNT i_dp;
  i_dp[0] = D_DP0;
  i_dp[1] = D_DP1;
  i_dp[2] = D_DP2;
  IPNT i_dv;
  i_dv[0] = D_DV0;
  i_dv[1] = D_DV1;
  i_dv[2] = D_DV2;  
    
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

  // cerr<<"done with pml arrays\n";
  
  /* sizes of computational domain - P0 same as P1, P2 */
  IPNT gsc_p;
  IPNT gec_p;
  int gsc_v[RARR_MAX_NDIM][RARR_MAX_NDIM];
  int gec_v[RARR_MAX_NDIM][RARR_MAX_NDIM];
  rd_gse(dom[0], i_p[0], gsc_p, gec_p);
  rd_gse(dom[0], i_v[0], gsc_v[0], gec_v[0]);
  rd_gse(dom[0], i_v[1], gsc_v[1], gec_v[1]);
  rd_gse(dom[0], i_v[2], gsc_v[2], gec_v[2]);

  // ndim is dimension of reference p
  int ndim;
  ra_ndim(&(dom[0]->_s[i_p[0]]),&ndim);

  // cerr<<"extend coefficient arrays into pml zones\n";
  
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

  // cerr<<"find loop limits for extended axes if any\n";
  // execute if first possible - any internal extended axis must
  // be symmetric and match axis lengths from model grid
  // and any other extended axis should have zero loop limits

  // must make assumption about extended axis, since no geom info
  // is available in RARR. Current assumptions:
  // - only 2D
  // - only one extended axis, for DBULK
  // - offset interval symmetric about h=0
  if (ndim !=2) {
    RVLException e;
    e<<"Error: asg_step\n";
    e<<"  current version only handles ndim=2\n";
    throw e;
  }
  IPNT gs;
  IPNT ge;
  rd_a_gse(dom[0], D_DBULK, gs, ge);
  if (2*((ge[2]-gs[2]+1)/2)+1  != ge[2]-gs[2]+1) {
    RVLException e;
    e<<"Error: asg_step - extended bulkmod pert\n";
    e<<"  extended axis has even number of samples so \n";
    e<<"  cannot be symmetric\n";
    throw e;
  }
  int ihmax=(ge[2]-gs[2]+1)/2;
  
  int n  = dom.size();

  // deriv = 0
  if (n==1) {
    //    if (fwd == true) {

    if (ndim == 2) {
      // cerr<<"pull out arrays\n";
      register ireal ** restrict bulk2 = (dom[0]->_s)[D_BULK]._s2;
      register ireal ** restrict buoy2 = (dom[0]->_s)[D_BUOY]._s2;
      register ireal ** restrict p02 = (dom[0]->_s)[i_p[0]]._s2;
      register ireal ** restrict p12 = (dom[0]->_s)[i_p[1]]._s2;
      register ireal ** restrict v02 = (dom[0]->_s)[i_v[0]]._s2;
      register ireal ** restrict v12 = (dom[0]->_s)[i_v[1]]._s2;
      register ireal ** restrict buoy2d = (dom[0]->_s)[D_DBUOY]._s2;
      register ireal *** restrict bulk3d = (dom[0]->_s)[D_DBULK]._s3;
      register ireal ** restrict p02d = (dom[0]->_s)[i_dp[0]]._s2;
      register ireal ** restrict p12d = (dom[0]->_s)[i_dp[1]]._s2;
      register ireal ** restrict v02d = (dom[0]->_s)[i_dv[0]]._s2;
      register ireal ** restrict v12d = (dom[0]->_s)[i_dv[1]]._s2;
      //      if (bulk3d) cerr<<"3d dbulk OK\n";
      //      else cerr<<"3d dbulk bad\n";
      
      //      cerr<<"allocate buffers\n";
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
	//	cerr<<"->basg2d_d\n";
	basg2d_d(bulk2, bulk3d, buoy2, buoy2d, p02, p02d, p12, p12d,
		 v02, v02d, v12, v12d, ep, epp, ev, evp,
		 sdiv, sdivd, sdivb, sdivbd, 
		 gradp0, gradp0d, gradp1, gradp1d,
		 gsc_p,gec_p,
		 &(gsc_v[0][0]),&(gec_v[0][0]),
		 &(gsc_v[1][0]),&(gec_v[1][0]),
		 asgpars->lbc,asgpars->rbc,
		 ihmax,		   
		 asgpars->k,asgpars->coeffs,
		 asgpars->stream);
	//	cerr<<"<-basg2d_d\n";
      }
      else {
	float pn=0.0f;
	float v0n=0.0f;
	float v1n=0.0f;
	for (int i1=gsc_p[1]; i1<=gec_p[1]; i1++) {
	  for (int i0=gsc_p[0]; i0<=gec_p[0]; i0++) {
	    pn+=p02[i1][i0]*p02[i1][i0];
	  }
	}
	for (int i1=gsc_v[0][1]; i1<=gec_v[0][1]; i1++) {
	  for (int i0=gsc_v[0][0]; i0<=gec_v[0][0]; i0++) {
	    v0n+=v02[i1][i0]*v02[i1][i0];
	  }
	}
	for (int i1=gsc_v[1][1]; i1<=gec_v[1][1]; i1++) {
	  for (int i0=gsc_v[1][0]; i0<=gec_v[1][0]; i0++) {
	    v1n+=v12[i1][i0]*v12[i1][i0];
	  }
	}
	fprintf(stderr,"after step: p02=%g v02=%g v12=%g\n",pn,v0n,v1n);
	
	basg2d_b(bulk2, bulk3d, buoy2, buoy2d, p02, p02d, p12, p12d,
		 v02, v02d, v12, v12d, ep, epp, ev, evp,
		 sdiv, sdivd, sdivb, sdivbd, 
		 gradp0, gradp0d, gradp1, gradp1d,
		 gsc_p,gec_p,
		 &(gsc_v[0][0]),&(gec_v[0][0]),
		 &(gsc_v[1][0]),&(gec_v[1][0]),
		 asgpars->lbc,asgpars->rbc,
		 ihmax,		   
		 asgpars->k,asgpars->coeffs,
		 asgpars->stream);
	float bn=0.0f;
	float pdn=0.0f;
	for (int i1=gsc_p[1]; i1<=gec_p[1]; i1++) {
	  for (int i0=gsc_p[0]; i0<=gec_p[0]; i0++) {
	    bn+=bulk3d[10][i1][i0]*bulk3d[10][i1][i0];
	    pdn+=p02d[i1][i0]*p02d[i1][i0];
	  }
	}
	fprintf(stderr,"after step: bulkd3sq=%g p02dsq=%g\n",bn,pdn);
      
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
      e<<"Error: asg_timestep(). 3D not implemented yet.\n";
    }
  }

  else {
    RVLException e;
    //    e<<"Error: asg_timestep(). Only 0th, 1st derivatives are implemented.\n";
    e<<"Error: basg_timestep(). Only 0th derivatives are implemented.\n";    
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
int basg_create_sten(void * specs,
		     FILE * stream,
		     int ndim,
		     IPNT gtype[RDOM_MAX_NARR],
		     STENCIL * sten) {
  
  // cerr<<"begin stencil\n";
    
  ASG_TS_PARS * asgpars = (ASG_TS_PARS *)(specs);
  STENCIL_MASK mask;/* workspace */
  int err = 0;
    
  /* number of masks 
     - every p depends on every v (dim^2)
     - every dp depends on every v (dim^2)
     - every dp depends on every dv (dim^2)
     - every v depends on corresp. p (dim)
     - every dv depends on corresp p (dim)
     - every dv depends on corresp dp (dim)
     - v's depend on buoy for some reason
     - dv's depend on dbuoy for some reason
   */
  int nmask = 3*(asgpars->ndim)*(asgpars->ndim)+3*(asgpars->ndim)+2*(asgpars->ndim); 
    
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

  IPNT iv;
  iv[0] = D_V0;
  iv[1] = D_V1;
  iv[2] = D_V2;

    IPNT idp;    
  idp[0] = D_DP0;
  idp[1] = D_DP1;
  idp[2] = D_DP2;

  IPNT idv;
  idv[0] = D_DV0;
  idv[1] = D_DV1;
  idv[2] = D_DV2;

  /* every mask stores 2k offsets */
  int len = 2*asgpars->k;
    
  /* loop over axes */
    
  int imask=0; // mask counter
    
  for (int idim=0; idim<asgpars->ndim; idim++) {
        
    for (int jdim=0; jdim<asgpars->ndim; jdim++) {

      // pj dep on vi
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

      /* source = v[idim], target = dp[jdim] */
      if ((err=mask_create(&mask, iv[idim], idp[jdim], len))) {
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
            
      /* source = dv[idim], target = dp[jdim] */
      if ((err=mask_create(&mask, idv[idim], idp[jdim], len))) {
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

    /* source = p, target = dv[idim] */
    // cerr<<"asg_stencil_create mask workspace\n";
    // cerr<<"asg_stencil_create imask="<<imask<<" idim="<<idim<<" mask workspace\n";
    if ((err=mask_create(&mask, ip[idim], idv[idim], len))) {
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

        /* source = dp, target = dv[idim] */
    // cerr<<"asg_stencil_create mask workspace\n";
    // cerr<<"asg_stencil_create imask="<<imask<<" idim="<<idim<<" mask workspace\n";
    if ((err=mask_create(&mask, idp[idim], idv[idim], len))) {
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

    if ((err=mask_create(&mask, D_DBUOY, idv[idim], 2))) {
      fprintf(stream,"ERROR: asg_create_sten from mask_create\n");
      sten_destroy(sten);
      return err;
    }
        
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
