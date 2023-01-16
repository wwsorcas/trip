#include "iwsim.hh"

//#define IWAVE_VERBOSE
//#define IWAVE_APPLY_VERBOSE
#define IWAVE_SIM_TIME
#ifdef IWAVE_SIM_TIME
#include <time.h>
#endif


namespace TSOpt {

  // helper function - not accessible outside this file!
  void synch(IWAVE * pstate, 
	     bool fwd,
	     int it,
	     int iv,
	     IWaveInfo const & ic,
	     FILE * stream) {
    
    // no-op unless MPI is defined
    
#ifdef IWAVE_USE_MPI
    
    int err = 0;
    // step info for the state to be updated

    int iopp;     /* opposite neighbor index */
    IPNT offs;    /* neighbor offset */
    MPI_Status status;  
    int tmpsend, tmprecv;
    MPI_Datatype tmpsend_dt, tmprecv_dt;
    int tmpsend_val, tmprecv_val;
    void *tmpsend_buf, *tmprecv_buf;
    //    double time;
    //    time = 0.0; /* To avoid "unitialized variable" warning */

    //    fprintf(stream,"IWAVE::synch: printact = %d\n",pstate->printact);

    /* we use the time step internal index for the perturbed field because
       this is the field being updated - on call should be same as index for
       unperturbed field, but don't use this */

    if ( (pstate->printact > 1) ) {
      if (fwd) {
	fprintf(stream,"\n\n**************************\n");
	fprintf(stream,    "****** SYNCH: FWD ********\n");
	fprintf(stream,    "**************************\n\n");
      }
      else {
	fprintf(stream,"\n\n**************************\n");
	fprintf(stream,    "****** SYNCH: BWD ********\n");
	fprintf(stream,    "**************************\n\n");
      }
    }
    if ( (pstate->printact > 6) ) {
      fprintf(stream,"\n------ synch: before exchange, step %d substep %d\n",it,iv);
      for (int ia=0;ia<RDOM_MAX_NARR;ia++) {
	if (fd_update(ia,iv,ic)) {
	  fprintf(stream,"------ iarr = %d\n",ia);
	  rd_print(&((pstate->model).ld_a), ia, stream);
	}
      }
      fflush(stream); 
    }

    for (int ia=0;ia<RDOM_MAX_NARR;ia++) {

      int upd=fd_update(ia,iv,ic);
      int isa=fd_isarr(ia,pstate->model,ic);
      
      if (pstate->printact>5) {
	fprintf(stream,"synch: array loop, test ia=%d iv=%d upd=%d isa=%d\n",
		ia,iv,upd,isa);
	fflush(stream);
      }
      
      //      if (fd_update(ia,iv,ic) && fd_isarr(ia,pstate->model,ic)) {
      if (upd && isa) {
	
	if ( (pstate->printact > 1) ) {
	  fprintf(stream,"\n------ synch fwd=%d array=%d -------------\n",fwd,ia);
	  fflush(stream); 
	}

	IPNT gs;
	IPNT ge;
	RARR sendsave;
	ra_setnull(&sendsave);

	for ( int i = 0; i < (pstate->model).nnei; ++i ) {
	  if ( (pstate->printact > 1) ) {
	    fprintf(stream,"array=%d neighbor=%d\n",ia,i);
	  }

	  // initialize all tmps to null values
	  tmpsend = MPI_PROC_NULL;
	  tmpsend_buf = &tmpsend_val;
	  tmpsend_dt = MPI_INT;
	  tmprecv = MPI_PROC_NULL;
	  tmprecv_buf = &tmprecv_val;
	  tmprecv_dt = MPI_INT;

	  if ( (pstate->pinfo).seinfo[ia][i].type != MPI_DATATYPE_NULL ) {  
	    
	    tmpsend = (pstate->pinfo).sranks[i];
	    tmpsend_buf = (pstate->pinfo).seinfo[ia][i].buf;
	    tmpsend_dt = (pstate->pinfo).seinfo[ia][i].type;
	      
	    if ( (pstate->printact > 1) ) {
	      fprintf(stream,"  sranks=%d\n",(pstate->pinfo).sranks[i]);
	      fflush(stream);
	    }
	    // adjoint case - save a copy of send RARR, which becomes recv buffer
	    // for adjoint

	    if (!fwd) {
	      rd_gse(&(((pstate->model).ld_s)[i]),ia,gs,ge);
	      // cleanup - why is this necessary
	      for (int k=(pstate->model).g.dim; k<RARR_MAX_NDIM; k++) { gs[k]=0; ge[k]=0; }
	      if ( (pstate->printact > 1) ) {
		fprintf(stream,"\n  bwd - copy send array with\n");
		fprintf(stream,"  gs[0]=%d gs[1]=%d gs[2]=%d\n",gs[0],gs[1],gs[2]);
		fprintf(stream,"  ge[0]=%d ge[1]=%d ge[2]=%d\n",ge[0],ge[1],ge[2]);
		fflush(stream); 
	      }
	      if ((err=ra_create(&sendsave,gs,ge))) {
		fprintf(stream,"\nError: synch from ra_create err=%d\n",err);
		RVLException e;
		e<<"\nError: IwaveSynch from ra_create err="<<err<<"\n";
		throw e;
	      }
	      
	      ra_zero(&sendsave);
	      if ((err=ra_copy(&sendsave,&(((((pstate->model).ld_s)[i])._s)[ia])))) {
		fprintf(stream,"\nError: synch from ra_copy err=%d\n",err);
		RVLException e;
		e<<"\nError: IWaveSynch from ra_copy err="<<err<<"\n";
		throw e;
	      }
	    }
	  }
	  else {
	    if ( (pstate->printact > 1) ) {
	      fprintf(stream,"\n  no send buffer\n");
	    }
	  }
	  
	  // having found send index, find recv index
	  if ((err=gen_i2pnt((pstate->model).g.dim,i,offs))) {
	    fprintf(stream,"\nError: synch from gen_i2pnt\n");
	    fprintf(stream,"  dim=%d i=%d\n",(pstate->model).g.dim, i);
	    RVLException e;
	    e<<"\nError: IWaveSynch from gen_i2pnt\n";
	    throw e;
	  }
	  for (int k=0;k<(pstate->model).g.dim;k++) offs[k]*=(-1);
	  if ((err=gen_pnt2i((pstate->model).g.dim,offs,&iopp))) {
	    fprintf(stream,"\nError: synch from gen_pnt2i\n");
	    RVLException e;
	    e<<"\nError: IWaveSynch from gen_i2pnt\n";
	    throw e;
	  }

	  if ( (pstate->pinfo).reinfo[ia][iopp].type != MPI_DATATYPE_NULL ) {
	    tmprecv = (pstate->pinfo).rranks[iopp];
	    tmprecv_buf = (pstate->pinfo).reinfo[ia][iopp].buf;
	    tmprecv_dt = (pstate->pinfo).reinfo[ia][iopp].type;
	    if ( (pstate->printact > 1) ) {
	      fprintf(stream," receive from %d rranks=%d\n",iopp,(pstate->pinfo).rranks[iopp]);
	      fflush(stream);
	    }
	  }
	  else {
	    if ( (pstate->printact > 1) ) {
	      fprintf(stream,"\n  no recv buffer\n");
	    }
	  }
	  /*	  
		  else {

		  RVLException e;
		  e<<"Error: IWaveSynch rk="<<rk<<"\n";
		  e<<"  array index "<<ia<<"\n";
		  e<<"  neighbor index "<<i<<" is active send index, but\n";
		  e<<"  opposite neighbor index "<<iopp<<" is not active recv index\n";
		  fprintf(stream,"ERROR: IWaveSynch rk=%d array=%d\n",rk,ia);
		  fprintf(stream,"  neighbor index=%d is active send index but\n",i);
		  fprintf(stream,"  opposite neighbor index=%d is not active recv index\n",iopp);
		  fprintf(stream,"THROWING EXCEPTION!!!\n");
		  fflush(stream);
		  throw e; 
		  }
		  }
	  */
	  

	  if ((pstate->printact > 1)) {
	    if ((tmpsend != MPI_PROC_NULL ) || (tmprecv != MPI_PROC_NULL)) fprintf(stream,"exchange with relative position %d: ",i);
	    if (tmpsend != MPI_PROC_NULL) {
	      fprintf(stream,"send to rk %d from virtual array\n",tmpsend);
	      if (fwd) rd_dump(&(((pstate->model).ld_s)[i]), ia, stream);
	      else ra_dump(&sendsave, stream);
	    }
	    if (tmprecv != MPI_PROC_NULL) {
	      fprintf(stream,"recv fr rk %d into virtual array\n",tmprecv);
	      rd_dump(&(((pstate->model).ld_r)[iopp]), ia, stream);
	    }
	    if ((tmpsend != MPI_PROC_NULL ) || (tmprecv != MPI_PROC_NULL)) fprintf(stream,"\n");
	    fflush(stream);
	  }
	  if ((tmpsend != MPI_PROC_NULL) && (pstate->printact > 5)) {
	    fprintf(stream,"\n------ synch: send buffer before exchange\n");
	    fprintf(stream,"  step %d substep %d field %d rel pos %d\n",
		    it,iv,ia,i);
	    if (fwd) rd_print(&(((pstate->model).ld_s)[i]), ia, stream);
	    else ra_print(&sendsave, stream);
	    fflush(stream); 
	  }
	  if ((tmprecv != MPI_PROC_NULL) && (pstate->printact > 5)) {
	    fprintf(stream,"\n------ synch: receive buffer before exchange\n");
	    fprintf(stream,"  step %d substep %d field %d rel pos %d\n",
		    it,iv,ia,i);
	    rd_print(&(((pstate->model).ld_r)[iopp]), ia, stream);
	    fflush(stream); 
	  }
	
	  /* 
	     "send" is send buffer
	     "recv" is receive buffer

	     if fwd: 
	     SEND data in tmpsend_buf to rk=tmpsend
	     RECEIVE data in tmprecv_buf from rk=tmprecv

	     else:
	     SEND data in tmprecv_buf to rk=tmprecv
	     RECEIVE data in tmpsend_buf from rk=tmpsend

	  */
	  //	  if ((tmpsend != MPI_PROC_NULL) && (tmprecv != MPI_PROC_NULL)) {
	  if ((tmpsend != MPI_PROC_NULL) || (tmprecv != MPI_PROC_NULL)) {
	    if (fwd) {
	      if (pstate->printact > 1) {
		fprintf(stream,"\n------- send/recv: recv=%d send=%d substep=%d\n",
			tmprecv, tmpsend, iv);
	      }
	      err = MPI_Sendrecv(tmpsend_buf, 1, tmpsend_dt, tmpsend, iv,
				 tmprecv_buf, 1, tmprecv_dt, tmprecv, iv,
				 (pstate->pinfo).ccomm, &status);
	      if ( err != MPI_SUCCESS ) {
		fprintf(stream, 
			"ERROR. Internal: MPI_Sendrecv error #%d, nei=%d, iv=%d, arr=%d. ABORT.\n", 
			err, i, iv, ia);
		RVLException e;
		e<<"    ERROR. Internal: MPI_Sendrecv error #"<<err<<", nei="<<i<<", iv="<<iv<<", arr="<<ia<<" ABORT\n";
		throw e;	      
	      }
	    }
	    else {
	      // first receive into send buffer, which has been copied into a tmp buffer
	      // reciprocally, send receive buffer
	      err = MPI_Sendrecv(  tmprecv_buf, 1, tmprecv_dt, tmprecv, iv,
				   tmpsend_buf, 1, tmpsend_dt, tmpsend, iv,
				   (pstate->pinfo).ccomm, &status);
	      
	      if ( err != MPI_SUCCESS ) {
		fprintf(stream, 
			"ERROR. Internal: MPI_Sendrecv error #%d, nei=%d, iv=%d, arr=%d. ABORT.\n", 
			err, i, iv, ia);
		RVLException e;
		e<<"    ERROR. Internal: MPI_Sendrecv error #"<<err<<", nei="<<i<<", iv="<<iv<<", arr="<<ia<<" ABORT\n";
		throw e;	      
	      }

	      // then add tmp buffer back
	      if ( ( (pstate->pinfo).seinfo[ia][i].type != MPI_DATATYPE_NULL ) ) {
		if (!(sendsave._s0)) {
		  fprintf(stream,"\nError: synch before axpy: sendsave not initialized\n");
		  RVLException e;
		  e<<"\nError: IWaveSynch before axpy: sendsave not initialized\n";
		  throw e;
		}
		if ((pstate->printact) > 5) {
		  fprintf(stream,"\n------ synch: delta send buffer after exchange \n");
		  fprintf(stream,"  step %d substep %d field %d rel pos %d\n",
			  it,iv,ia,i);
		  rd_dump(&(((pstate->model).ld_s)[i]), ia, stream);
		  rd_print(&(((pstate->model).ld_s)[i]), ia, stream);
		}

		fflush(stream);
		if ((err=ra_axpy(&(((((pstate->model).ld_s)[i])._s)[ia]),&sendsave,REAL_ONE))) {
		  fprintf(stream,"\nError: synch from ra_axpy err=%d\n",err);
		  ra_dump(&(((((pstate->model).ld_s)[i])._s)[ia]),stream);
		  ra_dump(&sendsave,stream);
		  ra_destroy(&sendsave);
		  RVLException e;
		  e<<"Error: IWaveSynch from ra_axpy err="<<err<<"\n";
		  throw e;
		}
		ra_destroy(&sendsave);

		// now comp subdom data is correct - synchronize
		/////// 13.04.14: apparently unnecessary, but zeroing the receive buffer IS!!
		/*	      
			      err = MPI_Sendrecv(tmpsend_buf, 1, tmpsend_dt, tmpsend, iv,
			      tmprecv_buf, 1, tmprecv_dt, tmprecv, iv,
			      (pstate->pinfo).ccomm, &status);
			      if ( err != MPI_SUCCESS ) {
			      fprintf(stream, 
			      "ERROR. Internal: MPI_Sendrecv error #%d, nei=%d, iv=%d, arr=%d. ABORT.\n", 
			      err, i, iv, ia);
			      RVLException e;
			      e<<"    ERROR. Internal: MPI_Sendrecv error #"<<err<<", nei="<<i<<", iv="<<iv<<", arr="<<ia<<" ABORT\n";
			      throw e;
			      }
		*/
		//now zero the receive buffer
		if ( ( (pstate->pinfo).reinfo[ia][iopp].type != MPI_DATATYPE_NULL ) ) { 
		  ra_zero(&(((((pstate->model).ld_r)[iopp])._s)[ia])); 
		}
	      }
	    }
	  }
	  if ((tmpsend != MPI_PROC_NULL) && (pstate->printact > 5)) {
	    fprintf(stream,"\n------ synch: send buffer after exchange\n");
	    fprintf(stream,"  step %d substep %d field %d rel pos %d\n",
		    it,iv,ia,i);
	    rd_print(&(((pstate->model).ld_s)[i]),ia, stream);
	    fflush(stream); 
	  }
	  if ((tmprecv != MPI_PROC_NULL) && (pstate->printact > 5)) {
	    fprintf(stream,"\n------ synch: receive buffer after exchange\n");
	    fprintf(stream,"  step %d substep %d field %d rel pos %d\n",
		    it,iv,ia,i);
	    rd_print(&(((pstate->model).ld_r)[iopp]), ia, stream);
	    fflush(stream); 
	  }
	}
      }
    }

    if ( (pstate->printact > 6) ) {
      fprintf(stream,"\n------ synch: after exchange, step %d substep %d\n",it,iv);
      for (int ia=0;ia<RDOM_MAX_NARR;ia++) {
	if (fd_update(ia,iv,ic)) {
	  fprintf(stream,"------ iarr = %d\n",ia);
	  rd_print(&((pstate->model).ld_a), ia, stream);
	}
      }
      fflush(stream); 
    }
#endif
  }


  IWaveSim::IWaveSim(int _order, 
		     bool _fwd, 
		     PARARRAY & pars, 
		     FILE * _stream, 
		     IWaveInfo const & _ic, 
		     int _printact,
		     bool _dryrun,
		     ostream & _drystr,
		     ostream & _announce)
    : ic(_ic), fwd(_fwd), stream(_stream), cstr(stream), 
      printact(_printact), order(_order),
      cps(NULL), narr(0), 
      dryrun(_dryrun), drystr(_drystr), 
      announce(_announce) {
    try {

#ifdef IWAVE_SIM_TIME
      clock_t t_start = clock();
#endif
      
      // cerr<<"iwavesim constr\n";
      
      // cerr<<"step 1: create list of i/o tasks\n";
      IOTask(t,order,fwd,ic);
#ifdef IWAVE_VERBOSE
      fprintf(stream,"IWaveSim constructor: fwd=%d order=%d\n",fwd,order);
      {
	std::stringstream ss;
	IOTaskWriter(t,ss);
	fprintf(stream,"%s\n",ss.str().c_str());
      }
      fflush(stream);
#endif
      // cerr<<"step 2: build state\n";
      w = new IWaveTree(pars, stream, ic, order);

      // cerr<<"step 2a: build grid.\n";
      // start with spatial grid, which has been initialized
      // in particular g.dim = problem spatial dimn
      // copy only the spatial axes, leaving the rest to
      // the sampler constructors
      init_default_grid(&g);
      int dim = (((w->getStateArray())[0])->model).g.dim;      
      int check = 0;
      for (int i=0;i<RARR_MAX_NDIM;i++) {
	int id = (((w->getStateArray())[0])->model).g.axes[i].id; 
	if ( (id> -1) && (id < dim) ) {
	  copy_axis(&(g.axes[id]),&((((w->getStateArray())[0])->model).g.axes[i]));
	  check++;
	}
      }
      if (check != dim) {
	RVLException e;
	e<<"Error: IWaveSim constructor\n";
	e<<"  grid in model indices funged up\n";
	throw e;
      }
      // trash all axes above spatial dimn
      for (int i=dim;i<RARR_MAX_NDIM;i++) {
	g.axes[i].n=0;
	g.axes[i].d=0.0;
	g.axes[i].o=0.0;
	g.axes[i].id=-1;
      }
      g.dim =  dim;
      g.gdim=g.dim;
      // cerr<<"---------------\n";
      // cerr<<"IWaveSim: sim grid g after initial construction from model.g:\n";
      //      fprint_grid(stderr,g);
      // cerr<<"---------------\n";

      // axis order is rigidly:
      // z, x, y (for 3d) or z, x (for 2d);
      // axis index g.dim is time
      // axis indices above g.dim are extended
      // note that this ordering has nothing to do with actual 
      // storage order in rsf data structure - correspondence is
      // by axis id
      
      // cerr<<"step 2a: build checkpoint structure if required\n";
      // sanity
      snaps = valparse<int>(pars,"nsnaps",0);
      if (!fwd && ((snaps<=0) && (order>0))) {
	RVLException e;
	e<<"Error: IWaveSim constructor\n";
	e<<"  adjoint mode (fwd=false):"<<"\n";
	e<<"  must provide positive number of \n";
	e<<"  storage units for checkpoints (argument \"nsnaps\")\n";
	e<<"  however snaps = "<<snaps<<"\n";
	throw e;
      }

      // work out number of dynamic arrays
      ndyn = 0;
      for (int i=0;i<RDOM_MAX_NARR;i++) {
	// cerr<<"i="<<i<<" ndyn="<<ndyn<<"\n";
	if (fd_isarr(i,(w->getStateArray())[0]->model,ic) && fd_isdyn(i,ic)) ndyn++; 
      }
      if (ndyn==0) {
	RVLException e;
	e<<"Error: IWaveSim constructor\n";
	e<<"  count of dynamic fields = 0\n";
	throw e;
      }

      // cerr<<"step 2b: in adjoint case, allocate checkpoint buffers\n";
      // IWAVEs 0,...,order-1 contain the reference data for the 
      // adjoint computation, so need order * snaps * ndyn RARRAYs.
      // NOTE ADDED 22.03.16: needed only if adjoint is computed for nonlinear
      // params, i.e. order>0
      if (!fwd && order>0) {
	typedef RARR* RP;
	typedef RARR** RPP;
	narr = snaps*pow2(order-1)*ndyn;
	RARR * cpstmp = new RARR[narr];
	RARR ** cpstmp2 = new RP[snaps*pow2(order-1)];
	for (size_t i=0;i<snaps*pow2(order-1);i++) cpstmp2[i]=&(cpstmp[i*ndyn]);
	cps = new RPP[snaps];
	for (int i=0;i<snaps;i++) cps[i]=&(cpstmp2[i*pow2(order-1)]);
	int l=0;
	for (int k=0;k<(w->getRefStateArray())[0]->model.ld_a.narr;k++) {
	  if (fd_isarr(k,(w->getRefStateArray())[0]->model,ic) && fd_isdyn(k,ic)) {
	    // pull out gs, ge for kth rarr in w->getRefStateArray[0] 
	    IPNT gs;
	    IPNT ge;
	    ra_a_gse(&((w->getRefStateArray())[0]->model.ld_a._s[k]),gs,ge);
	    //	    int ndim = (w->getRefStateArray())[0]->model.ld_a._s[k].ndim;
	    for (int i=0;i<snaps;i++) {
	      for (size_t j=0;j<pow2(order-1);j++) {
		// ra_create cps[i][j][l] -- check l
		if (l<0 || l>ndyn-1) {
		  RVLException e;
		  e<<"Error: IWaveSim: construct checkpoint array\n";
		  e<<"  bad news - index into rarray "<<l<<" out of range [0,"<<ndyn<<"\n";
		  throw e;
		}
		if (int err=ra_create(&(cps[i][j][l]),gs,ge)) {
		  RVLException e;
		  e<<"Error: IWaveSim from ra_create - checkpoint array\n";
		  e<<"  snap # = "<<i<<"\n";
		  e<<"  branch # = "<<j<<"\n";
		  e<<"  dyn arr # = "<<l<<"\n";
		  e<<"  rarr # = "<<k<<"\n";
		  e<<"  gs: ";
		  for (int i=0; i<RARR_MAX_NDIM; i++) e<<gs[i]<<" ";
		  e<<"\n";
	       	  e<<"  ge: ";
		  for (int i=0; i<RARR_MAX_NDIM; i++) e<<ge[i]<<" ";
		  e<<"\n";
		  e<<"  err = "<<err<<"\n";
		  throw e;
		}
	      }
	    }
	    l++;
	  }
	}
      }

      // cerr<<"step 3: construct list of samplers, axes\n";
      // note sampler exposed to first IWAVE - grid, time step info
      s.clear();
      for (size_t i=0; i<t.size(); i++) {
	IWaveSampler * tmp = NULL;
#ifdef IWAVE_VERBOSE
	fprintf(stream,"construct sampler %d on keyword %s\n",i,(t[i]->keyword).c_str());
#endif
	tmp = new IWaveSampler(w->getStateArray()[0], t[i]->keyword, pars, stream);
	// mod of 07.12.13: if current job does not provide source/sink for 
	// this i/o task, then sampler returns no axes. in that case, set
	// this pointer in s to NULL
	if (tmp->getNumAxes() == 0) {
#ifdef IWAVE_VERBOSE
	  cstr<<"sampler "<<i<<" on keyword "<<t[i]->keyword<<" is null\n";
	  cstr.flush();
#endif
	  delete tmp;
	  tmp = NULL;
	}
	else {
#ifdef IWAVE_VERBOSE
	  cstr<<"build sampler "<<i<<" on keyword "<<t[i]->keyword<<"\n";
	  cstr.flush();
#endif
	  /////// loop through sampler axes, add each to grid, then on
	  /////// success add sampler address to sim list
	  //	  cerr<<"sampler "<<i<<" on keyword "<<t[i]->keyword<<" is live\n";
	  // note these are internal simulation axes, not archival, 
	  // which will require some re-arranging in traceio
	  // grid_union(grid *, std::vector<axis *>)
	  //	  cerr<<"in sampler loop: sampler index = "<<i<<" keyword = "<<t[i]->keyword<<"\n";
	  //	  for (int k=0;k<tmp->getNumAxes(); k++) 
	  //	    fprint_axis(stderr,tmp->getAxis(k));
	  //	  fprintf(stderr,"  add to grid:\n");
	  //	  fprint_grid(stderr,g);
	  
	  // IMPORTANT MOD 11.12.14: since physical domain is spec'd by FIELD ID 0,
	  // only add spatial axes for FIELD ID 0. These axes were already fixed in 
	  // step 1, by initial construction of grid. So skip all subsequent axes
	  // with id < dim
	  
	  for (int j=0; j<tmp->getNumAxes();j++) {
	    // IMPORTANT CHANGE 07.12.13: the time axis is special, because
	    // the internal step is already fixed by the IWAVE constructor.
	    // therefore, change this axis to have the simulation time step
	    // before registering it via grid_union.
	    // Change 08.12.13: copy axis to workspace - only const access
	    // now allowed
	    axis * a = new axis;
	    init_default_axis(a);
	    copy_axis(a, &(tmp->getAxis(j)));
	    float tntmp = 0.0f;
	    if (a->id == g.dim) {     // signals time axis!!
	      if (ProtectedDivision<float>((a->n * a->d),
					   (w->getStateArray()[0]->model.tsind.dt),
					   tntmp)) {
		RVLException e;
		e<<"Error: IWaveSim constructor\n";
		e<<"  zerodivide by dt stored in IWaveTree model struct\n";
		e<<"  dt = "<<w->getStateArray()[0]->model.tsind.dt<<"\n";
		throw e;
	      }
	      // adjust number, spacing
	      a->n = iwave_max(1, (int)(tntmp+0.1f));
	      a->d = w->getStateArray()[0]->model.tsind.dt;
	      // adjust origin so that in all cases 0 is a (virtual) grid point
	      // this assures that all time axes are compatible. Adjust o and n
	      // so that new grid interval contains old.
	      float tmpo = a->o;
	      a->o =
		//		((int)(((a->o)/(a->d))-0.1))
		((int)((a->o)/(a->d)))
		*(a->d);
	      if (tmpo<a->o) {
		a->o -= a->d;
		a->n ++;
	      }
	      // notice that now simulation time axis has been DETERMINED, so 
	      // must be used in trace sampler!!!!
	    }
	    
	    if ((a->id > g.dim-1) && (!grid_union(&g,a))) {
	      //	      cerr<<"Error: IWaveSim constructor from grid_union\n";
	      RVLException e;
	      e<<"Error: IWaveSim constructor from grid_union\n";
	      fprintf(stream,"Error: IWaveSim constructor from grid_union\n");
	      fprintf(stream,"  failed to add these axes:\n");
	      for (int k=0;k<tmp->getNumAxes(); k++) 
		fprint_axis(stream,tmp->getAxis(k));
	      fprintf(stream,"  to grid:\n");
	      fprint_grid(stream,g);
	      throw e;
	    }
	    // having added to grid, or not, trash it
	    delete a;
	    
	  }
	  //      	  fprintf(stderr,"  result grid:\n");
	  //       	  fprint_grid(stderr,g);
	}
#ifdef IWAVE_VERBOSE
	cstr<<"finished sampler "<<i<<" on keyword "<<t[i]->keyword<<"\n";
	cstr.flush();
#endif	  
	s.push_back(tmp);
      }
#ifdef IWAVE_VERBOSE
      cstr<<"finished sampler construction loop\n";
      cstr<<"simulation grid constructed:\n";
      fprint_grid(stream,g);      
      cstr.flush();
#endif	  	
      // 08.01.14: panelindex construction makes this unnecessary
      //      fprint_grid(stderr,g);
      // step 4: modify grid if necessary
      // no-op if no mpi
      //#ifdef IWAVE_USE_MPI
      //      mpi_update_grid(&g);
      //#endif
      // cerr<<"exit iwavesim constructor\n";

#ifdef IWAVE_SIM_TIME
      cstr<<"IWaveSim constructor time = "<<((double)(clock()-t_start))/CLOCKS_PER_SEC
	  <<"\n";
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from IWaveSim constructor\n";
      throw e;
    }
  }
  
  IWaveSim::~IWaveSim() {
    //    cerr<<"destructor\n";
    for (size_t i=0; i<t.size(); i++) {
      //cerr<<"destroy sampler "<<i<<"\n";
      if (s.at(i)) {
	//	cerr<<"destroy sampler "<<i<<"\n";
	delete s.at(i);
      }
    }
    //    cerr<<"destroy state\n";
    if (w) delete w;
    for (size_t i=0; i<t.size(); i++) {
      //      cerr<<"destroy iotask "<<i<<"\n";
      if (t.at(i)) delete t.at(i); 
    }
    if (cps) {
      //      cerr<<"destroy checkpoints\n";
      for (int i=0;i<snaps;i++) {
	for (size_t j=0;j<pow2(order-1);j++) {
	  for (int k=0;k<ndyn;k++) {
	    ra_destroy(&(cps[i][j][k]));
	  }
	}
      }
      delete [] &(cps[0][0][0]);
      delete [] &(cps[0][0]);
      delete [] &(cps[0]);
    }
  }

  void IWaveSim::run() {
    try {
#ifdef IWAVE_VERBOSE
      cstr<<"IwaveSim::run\n";
      cstr.flush();
#endif
#ifdef IWAVE_SIM_TIME
      clock_t t_start = clock();
#endif
      // initialize step
      IPNT step;
      IPNT start;
      IPNT stop;
      get_gs(start, g);
      get_ge(stop, g);
      // initially, assign step = start
      IASN(step,start);
      // next, work out first and last panels for extended
      // axes, this process
      int first;
      int last;
      int panelnum=1;
      int panelindex=0;
      for (int i=g.dim+1; i<g.gdim; i++) {
	panelnum *= stop[i]-start[i]+1;
      }
      calc_group(&first, &last, panelnum, stream);
#ifdef IWAVE_VERBOSE
      cstr<<"panelnum="<<panelnum<<" first="<<first<<" last="<<last<<"\n";
      cstr.flush();      
#endif
      if ((first > last) && (fwd==0)) {
	RVLException e;
	e<<"ERROR: IWaveSim::run, adjoint simulation\n";
	e<<"  number of MPI process groups = "<<retrieveNumGroups()<<"\n";
	e<<"  number of data gathers = "<<panelnum<<"\n";
	e<<"  because adjoint simulation requires reduction,\n";
	e<<"  all process groups must actively participate,\n";
	e<<"  so number of data gathers cannot be less than\n";
	e<<"  number of process groups\n";
	e<<"  Try again with fewer MPI processes\n";
	throw e;
      }
      
      // increment step until at first panel
      while (panelindex < first) {
	// update step array
	for (int i=g.dim+1;i<g.gdim;i++) {
	  if (step[i]<stop[i]) step[i]++;
	  else {
	    if (i<g.gdim-1) step[i]=start[i];
	  }
	}
	panelindex++;
      }

      // pull out fdpars for use in time step - same in every
      // step, and for every RDOM, so do it once here and get 
      // from root RDOM
      void * fdm = NULL;
      if (!dryrun) {
	fdm = (w->getStateArray()[0]->model.specs);
      }

      // time axis - index=dim: last data point if fwd = false (reverse)
      
      // step 5: initial sample - includes any post-contruction
      // initialization
      
      // step 6: sim loop - includes time, record,...
      // time step function knows what to do based 
      // merely on size. use next_step to increment
      // extended axes - time step loop is explicit

      // here more = more records
      for (panelindex=first; panelindex<=last; panelindex++) {
#ifdef IWAVE_VERBOSE
	cstr<<"panelindex="<<panelindex<<" first="<<first<<" last="<<last<<"\n";
	cstr.flush();
#endif
#ifdef IWAVE_SIM_TIME
	if (panelindex==first) {
	  cstr<<"IWaveSim::run: panelindex="<<panelindex<<"time="
	      <<((double)(clock()-t_start))/CLOCKS_PER_SEC<<"\n";
	}
	else {
	  cstr<<"IWaveSim::run: panelindex="<<panelindex<<"\n";
	}
	cstr.flush();
	t_start=clock();
#endif
	/*
	  for (int i=0;i<w->getStateArray().size();i++) {
	  if (int err=rd_a_zero(&((w->getStateArray()[i]->model).ld_a))) {
	  RVLException e;
	  e<<"Error: IWaveTree main constructor, IWAVE["<<i<<"]\n";
	  e<<"  returning error "<<err<<" from rd_a_zero\n";
	  throw e;
	  }
	  }
	*/
	// compute panelindex = step through external extended indices
	if (dryrun) {
	  if (panelindex==first) {
	    drystr<<"\n*** IWaveSim: model="<<ic.iwave_model<<" fwd="<<fwd<<" deriv="<<order<<"\n";
	  }
	  drystr<<"    panel="<<panelindex<<" in range=["<<first<<", "<<last<<"], rank="<<retrieveGlobalRank()<<"\n";
	}
	else {
	  if (printact > -1) {
	    if (panelindex==first && retrieveGlobalRank()==0) {
	      //	      fprintf(stream,"\n*** IWaveSim: model=%s  fwd=%d deriv=%d\n",ic.iwave_model.c_str(), fwd, order);
	      fprintf(stderr,"\n*** IWaveSim: model=%s  fwd=%d deriv=%d\n",ic.iwave_model.c_str(), fwd, order);
	    }
	    //	    fprintf(stream,"    panel=%d in range=[%d, %d], rank=%d\n", panelindex,first,last,retrieveGlobalRank());
	    fprintf(stderr,"    panel=%d in range=[%d, %d], rank=%d\n", panelindex,first,last,retrieveGlobalRank());
	  }

#ifdef IWAVE_VERBOSE
	  cstr<<"simulation grid:\n";
	  for (int i=0;i<g.gdim;i++) {
	    if (i<=g.dim) step[i]=start[i];
	    cstr<<"axis "<<i<<" n="<<g.axes[i].n<<" o="<<g.axes[i].o<<" d="<<g.axes[i].d<<" id="<<g.axes[i].id<<" gs="<<start[i]<<" ge="<<stop[i]<<" step="<<step[i]<<"\n";
	  }
	  cstr.flush();
#endif
	}
	if (dryrun) {
	  drystr<<"\nIWaveSim::run - initialize dynamic arrays"<<"\n";
	  drystr<<"simulation grid = \n";
	  for (int i=0;i<g.gdim;i++) {
	    if (i<=g.dim) step[i]=start[i];
	    drystr<<"axis "<<i<<" n="<<g.axes[i].n<<" o="<<g.axes[i].o<<" d="<<g.axes[i].d<<" id="<<g.axes[i].id<<" gs="<<start[i]<<" ge="<<stop[i]<<" step="<<step[i]<<"\n";
	  }
	  drystr<<"\n\n";
	}

	/******* FORWARD LOOP *********/

#ifdef IWAVE_SIM_TIME
	cstr<<"IWaveSim::run - setup time = "
	    <<((double)(clock()-t_start))/CLOCKS_PER_SEC<<"\n";
	t_start=clock();
#endif
	if (fwd) {

	  

	  // cerr<<"\nIWaveSim::run - FORWARD WET RUN\n\n";
	  if (dryrun) {
	    drystr<<"\nIWaveSim::run - FORWARD DRY RUN\n\n";
	  }
#ifdef IWAVE_VERBOSE
	  cstr<<"IWaveSim::run - dynamic init\n";
	  cstr.flush();
#endif
	  for (size_t i=0;i<w->getStateArray().size();i++) {
	    iwave_dynamic_init(w->getStateArray()[i],start[g.dim],ic);
	  }

#ifdef IWAVE_SIM_TIME
	  clock_t t_step = 0.0;
	  clock_t t_samp = 0.0;
#endif
	  //	  for (int it=start[g.dim]; it<stop[g.dim]; it++) {
	  for (int it=0; it<stop[g.dim]-start[g.dim]; it++) {

	    float dt = g.axes[g.dim].d;
	    float ot = g.axes[g.dim].o;
	    //	    float tmax = dt*(g.axes[g.dim].n);

#ifdef IWAVE_VERBOSE
	    cstr<<"IWaveSim::run - it="<<it<<" t="<<ot+it*dt<<"\n";
	    cstr.flush();
#endif
	    if (printact > 0)
	      fprintf(stream,"    it=%d t=%12.4e\n",it,ot+it*dt);
	  
	    if (dryrun) drystr<<"\n";
	    step[g.dim]=it+start[g.dim];

#ifdef IWAVE_SIM_TIME
	    clock_t dt_samp=clock();
#endif
	    for (size_t i=0; i<t.size(); i++) {
	      if (s[i]) {
#ifdef IWAVE_VERBOSE
		cstr<<"IWaveSim::run - rk="<<retrieveGlobalRank()<<" -> sampler["<<i<<"]"<<"g.dim="<<g.dim<<" step="<<step[g.dim]<<"\n";
		cstr.flush();		
#endif
		s[i]->sample(g,step,ic.get_iwave_fields()[t[i]->rarrindex].gtype,fwd,
			     t[i]->input,w->getStateArray()[t[i]->iwaveindex],
			     t[i]->rarrindex,t[i]->iwaveindex,stream,
			     dryrun,drystr);
#ifdef IWAVE_VERBOSE
		cstr<<"IWaveSim::run - rk="<<retrieveGlobalRank()<<" <- sampler["<<i<<"]\n";
		cstr.flush();
#endif
	      }

	      // refine loop def - use info from linear (source) inputs only
	      /*
		if (!dryrun) {
		if ((it==start[g.dim]) &&
		(ic.get_property_iokeys(t[i]->keyword,"active") == ACTIVE_LINEAR)) {
		ic.get_loopdef()(s[i]->get_gmin(),s[i]->get_gmax(),tmax,t[i]->input,stream,fdm);
		}
		}
	      */
	    }
#ifdef IWAVE_SIM_TIME
	    t_samp += clock()-dt_samp;
#endif
	    
	    // first time through, call check 
#ifdef IWAVE_VERBOSE
	    cstr<<"IWaveSim::run - call check\n";
	    cstr.flush();	    
#endif
            if (!dryrun) {
	      //	      if (it==start[g.dim]) {
	      if (it==0) {
		ic.get_check()(w->getRDOMArray()[0],fdm,stream);
	      }
            }
	    if (dryrun) {
	      drystr<<"\nIWaveSim::run fwd step "<<it<<" -> "<<it+1<<"\n";
	    }
	    else { 
#ifdef IWAVE_VERBOSE
	      cstr<<"IWaveSim::run - rk="<<retrieveGlobalRank()<<": timestep\n";
	      cstr.flush();
#endif
	      for (int iv=0;iv<fd_numsubsteps(ic);iv++) {
		//		cerr<<"timestep\n";
#ifdef IWAVE_USE_MPI
		//		MPI_Barrier(retrieveComm());
#endif
		if (printact>5) {
		  fprintf(stream,"iwsim::run timestep it=%d iv=%d\n",
			  it,iv);
		  fflush(stream);
		}
#ifdef IWAVE_SIM_TIME
		clock_t dt_step = clock();
#endif
		
		ic.get_timestep()(w->getRDOMArray(),fwd,iv,fdm);
#ifdef IWAVE_SIM_TIME
		t_step += clock()-dt_step;
#endif
		if (printact>5) {
		  fprintf(stream,"iwsim::run return from timestep\n");
		  fflush(stream);
		}
		// in fwd loop, synch ALL dynamic arrays
#ifdef IWAVE_USE_MPI
		//		MPI_Barrier(retrieveComm());
#endif
		for (size_t k=0; k<w->getStateArray().size(); k++) {
		  if (w->getStateArray()[k]->printact > 5) {
		    fprintf(stream,"iwsim::run: before synch k=%lu it=%d iv=%d\n"
			    ,k,it,iv);
		    fflush(stream);
		  }
#ifdef IWAVE_VERBOSE
		  cstr<<"IWaveSim::run -synch substep="<<iv<<" iwdx="<<k<<"\n";
		  cstr.flush();
#endif
		  synch(w->getStateArray()[k],fwd,it,iv,ic,stream);
		  if (printact>5) {
		    fprintf(stream,"iwsim::run: after synch k=%lu it=%d iv=%d\n",
			    k,it,iv);
		    fflush(stream);
		  }
		}
	      }
	    }
	  }
	  
	  if (dryrun) drystr<<"\n";
	  step[g.dim]=stop[g.dim];
	  for (size_t i=0; i<t.size(); i++) {
	    if (s[i]) {
#ifdef IWAVE_VERBOSE
	      cstr<<"IWaveSim::run - rk="<<retrieveGlobalRank()<<" sampler["<<i<<"]\n";
	      cstr.flush();
#endif
	      s[i]->sample(g,step,ic.get_iwave_fields()[t[i]->rarrindex].gtype,fwd,t[i]->input,
			   w->getStateArray()[t[i]->iwaveindex],
			   t[i]->rarrindex,t[i]->iwaveindex,stream,
			   dryrun,drystr);
	    }
	  }
#ifdef IWAVE_SIM_TIME
	  cstr<<"IWaveSim::run: loop time = "
	      <<((double)(clock()-t_start))/CLOCKS_PER_SEC<<"\n";
	  cstr<<"               samp time = "
	      <<((double)(t_samp))/CLOCKS_PER_SEC<<"\n";
	  cstr<<"               step time = "
	      <<((double)(t_step))/CLOCKS_PER_SEC<<"\n";
	  t_start=clock();
#endif
	}

	/******* ADJOINT LOOP, ORDER=0 - SOURCE CASE *********/

	else if ((fwd==0) && (order==0)) {

	  // cerr<<"\nIWaveSim::run - BACKWARD WET RUN\n\n";
	  if (dryrun) {
	    drystr<<"\nIWaveSim::run - BACKWARD DRY RUN\n\n";
	  }
#ifdef IWAVE_VERBOSE
	  cstr<<"IWaveSim::run - dynamic init\n";
	  cstr.flush();
#endif
	  for (size_t i=0;i<w->getStateArray().size();i++) {
	    iwave_dynamic_init(w->getStateArray()[i],start[g.dim],ic);
	  }

	  for (int it=stop[g.dim]; it>=start[g.dim]; it--) {

	    float dt = g.axes[g.dim].d;
	    float ot = g.axes[g.dim].o;

#ifdef IWAVE_VERBOSE
	    cstr<<"IWaveSim::run - dt="<<dt<<" ot="<<ot<<"\n";
	    cstr<<"IWaveSim::run - it="<<it<<" t="<<ot+it*dt<<"\n";
	    cstr.flush();
#endif
	    if (printact > 0)
	      fprintf(stream,"    it=%d t=%12.4e\n",it,ot+it*dt);
	  
	    if (dryrun) drystr<<"\n";
	    step[g.dim]=it;
  
	    for (size_t i=0; i<t.size(); i++) {
	      if (s[i]) {
#ifdef IWAVE_VERBOSE
		cstr<<"IWaveSim::run - rk="<<retrieveGlobalRank()<<" sampler["<<i<<"]\n";
		cstr.flush();
#endif
		s[i]->sample(g,step,ic.get_iwave_fields()[t[i]->rarrindex].gtype,fwd,
			     t[i]->input,w->getStateArray()[t[i]->iwaveindex],
			     t[i]->rarrindex,t[i]->iwaveindex,stream,
			     dryrun,drystr);
	      }
	    }
	    
	    // first time through, call check
#ifdef IWAVE_VERBOSE
	    cstr<<"IWaveSim::run - call check\n";
	    cstr.flush();
#endif
            if (!dryrun) {
	      if (it==start[g.dim]) ic.get_check()(w->getRDOMArray()[0],fdm,stream);
            }
	    if (dryrun) {
	      drystr<<"\nIWaveSim::run fwd step "<<it<<" -> "<<it+1<<"\n";
	    }
	    else { 
#ifdef IWAVE_VERBOSE
	      cstr<<"IWaveSim::run - rk="<<retrieveGlobalRank()<<": timestep\n";
	      cstr.flush();
#endif
	      //// NOTE - BACKWARDS!!!
	      //	      for (int iv=fd_numsubsteps(ic)-1;iv>-1;iv--) {
	      for (int iv=0; iv<fd_numsubsteps(ic);iv++) {
		//		cerr<<"timestep\n";
		ic.get_timestep()(w->getRDOMArray(),fwd,iv,fdm);
		// in fwd loop, synch ALL dynamic arrays
		
		for (size_t k=0; k<w->getStateArray().size(); k++) {
		  //		  if (w->getStateArray()[k]->printact > 5) 
		  //		    fprintf(stream,"\n*** SYNCH: iwdx = %d\n\n",k);
#ifdef IWAVE_VERBOSE
		  cstr<<"IWaveSim::run -synch substep="<<iv<<" iwdx="<<k<<"\n";
		  cstr.flush();
#endif
		  synch(w->getStateArray()[k],fwd,it,iv,ic,stream);
		  //		  cerr<<"  k="<<k<<"\n";
		}
	      }
	    }
	  }
	  
	  if (dryrun) drystr<<"\n";
	  step[g.dim]=stop[g.dim];
	  for (size_t i=0; i<t.size(); i++) {
	    if (s[i]) {
#ifdef IWAVE_VERBOSE
	      cstr<<"IWaveSim::run - rk="<<retrieveGlobalRank()<<" sampler["<<i<<"]\n";
	      cstr.flush();
#endif
	      s[i]->sample(g,step,ic.get_iwave_fields()[t[i]->rarrindex].gtype,fwd,
			   t[i]->input,w->getStateArray()[t[i]->iwaveindex],
			   t[i]->rarrindex,t[i]->iwaveindex,stream,
			   dryrun,drystr);
	    }
	  }
	}

	/******** ONE STEP CASE OF ADJOINT *******/

	else if (stop[g.dim]-start[g.dim]==1) {
	  
	  for (size_t i=0;i<w->getStateArray().size();i++) {
	    iwave_dynamic_init(w->getStateArray()[i],start[g.dim],ic);
	  }

	  int it = start[g.dim]; // ref state index
	  int at = stop[g.dim];  // pert state index
	  
	  if (dryrun) {
	    drystr<<"\nIWaveSim::run - ADJOINT DRY RUN\n";
	    drystr<<"  SINGLE STEP CASE\n\n";
	  }	 
#ifdef IWAVE_VERBOSE
	  cstr<<"\nIWaveSim::run - ADJOINT DRY RUN\n";
	  cstr<<"  SINGLE STEP CASE\n\n";
	  cstr.flush();
#endif
	  // load reference data
	  step[g.dim]=it;
	  bool reffwd = true;
	  for (size_t i=0; i<t.size(); i++) {
	    // sample data for reference
	    if (t[i]->iwaveindex < (int)pow2(order-1) && s[i]) {
	      s[i]->sample(g,step,ic.get_iwave_fields()[t[i]->rarrindex].gtype,reffwd,
			   t[i]->input,w->getStateArray()[t[i]->iwaveindex],
			   t[i]->rarrindex,t[i]->iwaveindex,stream,
			   dryrun,drystr);
	    }
	  }
	  // check
          if (!dryrun) {
	    ic.get_check()(w->getRDOMArray()[0],fdm,stream);
          }

	  // load adjoint data 
	  step[g.dim]=at;
	  for (size_t i=0; i<t.size(); i++) {
	    // pert sample
	    if (t[i]->iwaveindex >= (int)pow2(order-1) && s[i]) { 
	      s[i]->sample(g,step,ic.get_iwave_fields()[t[i]->rarrindex].gtype,fwd,
			   t[i]->input,w->getStateArray()[t[i]->iwaveindex],
			   t[i]->rarrindex,t[i]->iwaveindex,stream,
			   dryrun,drystr);
	    }
	  }
	  if (!dryrun) {
	    // backwards step - only need to synch top-order pert arrays
	    for (int iv=fd_numsubsteps(ic)-1; iv>-1; iv--) {
	      ic.get_timestep()(w->getRDOMArray(),fwd,iv,fdm);
	      for (size_t k=pow2(order-1); k<pow2(order); k++) {
		//		if (w->getStateArray()[k]->printact > 5) 
		//		  fprintf(stream,"\n*** SYNCH: iwdx = %d\n\n",k);
		synch(w->getStateArray()[k],fwd,at,iv,ic,stream);
	      }
	    }	
	  }
	  else {
	    drystr<<"IWaveSim::run adj - first adjoint step: at="<<at;
	  }
	  at--;
	  if (dryrun) {
	    drystr<<"->"<<at<<"; ref step = "<<it<<"\n";
	  }
	  step[g.dim]=at;

	  for (size_t i=0; i<t.size(); i++) {
	    // pert sample
	    if (t[i]->iwaveindex >= (int)pow2(order-1) && s[i]) { 
	      s[i]->sample(g,step,ic.get_iwave_fields()[t[i]->rarrindex].gtype,fwd,
			   t[i]->input,w->getStateArray()[t[i]->iwaveindex],
			   t[i]->rarrindex,t[i]->iwaveindex,stream,
			   dryrun,drystr);
	    }
	  }

	}

	/****** ADJOINT STATE LOOP *******/
	else {
	  // flag for step, sample on fwd steps
	  bool reffwd = true;

	  // checkpointing algorithm - should be abstracted!!!!
	  
	  // create checkpoint schedule
	  ostringstream ostr;
	  Revolve r(stop[g.dim]-start[g.dim],snaps,ostr);

	  for (int i=0;i<snaps;i++) {
	    for (int j=0;j<(int)pow2(order-1);j++) {
	      for (int l=0;l<ndyn;l++) {
		if (int err=ra_a_zero(&(cps[i][j][l]))) {
		  RVLException e;
		  e<<"Error: IWaveSim from ra_create - zero rarray\n"; 
		  e<<"  snapindex="<<i<<" iwaveindex="<<j<<" rarrindex="<<l<<"\n";
		  e<<"  err = "<<err<<"\n";
		}
	      }
	    }
	  }

	  std::vector<int> cplist(snaps);
	  int it = start[g.dim]; // ref state index
	  int at = stop[g.dim];  // pert state index

#ifdef IWAVE_VERBOSE
	  cstr<<"\nIWaveSim::run adj - load initial data time step "<<it<<"\n";
	  cstr.flush();
#endif
	  step[g.dim]=it;

	  for (size_t i=0;i<w->getStateArray().size();i++) {
	    iwave_dynamic_init(w->getStateArray()[i],start[g.dim],ic);
	  }
	  
	  for (size_t i=0; i<t.size(); i++) {
	    // sample data for reference
	    if (t[i]->iwaveindex < (int)pow2(order-1) && s[i]) {
	      s[i]->sample(g,step,ic.get_iwave_fields()[t[i]->rarrindex].gtype,reffwd,
			   t[i]->input,w->getStateArray()[t[i]->iwaveindex],
			   t[i]->rarrindex,t[i]->iwaveindex,stream,
			   dryrun,drystr);
	    }
	  }

	  // check if it=start
          if (!dryrun) {
	    if (it==start[g.dim]) ic.get_check()(w->getRDOMArray()[0],fdm,stream);
          }

	  ACTION::action whatodo;

	  if (dryrun) {
	    drystr<<"\nIWaveSim::run - ADJOINT DRY RUN\n";
	    drystr<<"  FORWARD TIME LOOP\n\n";
	  }
	  do {
#ifdef IWAVE_VERBOSE
	    cstr<<"top of adjoint loop, it="<<it<<" at="<<at<<"\n";
	    cstr.flush();
#endif
	    whatodo = r.revolve();
#ifdef IWAVE_VERBOSE
	    cstr<<"whatodo = "<<whatodo<<"\n";
	    cstr.flush();
#endif
	    if (whatodo == ACTION::takeshot) {
	      int cp = r.getcheck();
#ifdef IWAVE_VERBOSE
	      cstr<<"\nIWaveSim::run adj - store step "<<it<<" in checkpoint "<<cp<<"\n";
	      cstr.flush();
#endif
	      for (int j=0;j<(int)pow2(order-1);j++) {
		int l = 0;
		for (int k=0;k<RDOM_MAX_NARR;k++) {
		  if (fd_isarr(k,w->getStateArray()[0]->model,ic) && fd_isdyn(k,ic)) {		  
		    if (ra_a_copy(&(cps[cp][j][l]),&(((w->getRefRDOMArray())[j])->_s[k]))) {
		      RVLException e;
		      e<<"Error: IWaveSim::run\n";
		      e<<"attempt to store checkpoint "<<cp
		       <<" encountered error from ra_a_copy\n";
		      e<<"IWAVE index = "<<j<<" checkpoint RARR index = "
		       <<l<<" RDOM RARR index = "<<k<<"\n";
		      throw e;
		    }
#ifdef IWAVE_VERBOSE
		    float tmp=0.0;
		    ra_a_inner(&(((w->getRefRDOMArray())[j])->_s[k]),
			       &(((w->getRefRDOMArray())[j])->_s[k]),
			       &tmp);
		    cstr<<"  level="<<j<<" cpindx="<<l<<" ridx="<<k<<" normsq="<<tmp<<"\n";
		    cstr.flush();
#endif
		    l++;
		  }
		}
	      }
	      cplist.at(cp)=it;
	      if (dryrun) {
		drystr<<"\nIWaveSim::run adj - stored step "<<it<<" in checkpoint "<<cp<<"\n";
	      }
	    }
	    if (whatodo == ACTION::advance) { 
	      for (int j = r.getoldcapo(); j < r.getcapo(); j++) {
		if (dryrun) drystr<<"\n";
		if (!dryrun) {
		  for (int iv=0;iv<fd_numsubsteps(ic);iv++) {
		    ic.get_timestep()(w->getRefRDOMArray(),reffwd,iv,fdm);
		    for (size_t k=0; k<pow2(order-1); k++) {
		      //		      if (w->getStateArray()[k]->printact > 5) 
		      //			fprintf(stream,"\n*** SYNCH: iwdx = %d\n\n",k);
		      synch(w->getRefStateArray()[k],reffwd,it,iv,ic,stream);
		    }
		  }
#ifdef IWAVE_VERBOSE	
		  cstr<<"IWaveSim::run adj - fwd step "<<it;
		  cstr.flush();
#endif
		}
		else {
		  drystr<<"IWaveSim::run adj - fwd step "<<it;
		}
		it++;
#ifdef IWAVE_VERBOSE
		cstr<<"->"<<it<<"\n";
		cstr.flush();
#endif
		if (dryrun) {
		  drystr<<"->"<<it<<"\n";
		}
		step[g.dim]=it;
		for (size_t i=0; i<t.size(); i++) {
		  // sample data for reference
		  if (t[i]->iwaveindex < (int)pow2(order-1) && s[i]) {
		    s[i]->sample(g,step,ic.get_iwave_fields()[t[i]->rarrindex].gtype,reffwd,
				 t[i]->input,w->getStateArray()[t[i]->iwaveindex],
				 t[i]->rarrindex,t[i]->iwaveindex,stream,
				 dryrun,drystr);
		  }
		}
	      }
	    }
	    if (whatodo == ACTION::firsturn) {
#ifdef IWAVE_VERBOSE
	      cstr<<"\nIWaveSim::run - ADJOINT WET RUN\n";
	      cstr<<"  REVERSE TIME LOOP\n\n";
	      cstr.flush();
#endif
	      if (dryrun) {
		drystr<<"\nIWaveSim::run - ADJOINT DRY RUN\n";
		drystr<<"  REVERSE TIME LOOP\n\n";
	      }	      
	      step[g.dim]=at;
	      for (size_t i=0; i<t.size(); i++) {
		// need to clean output arrays at outset of reverse time loop
		if (t[i]->iwaveindex >= (int)pow2(order-1)) {
		  // iwave_dynamic_init(w->getStateArray()[t[i]->iwaveindex],at,ic);
		  //		  if (!(t[i]->input) )
		  //		    ra_a_zero(&(w->getStateArray()[t[i]->iwaveindex]->model.ld_a._s[t[i]->rarrindex]));
		}
		// pert sample 
		if (t[i]->iwaveindex >= (int)pow2(order-1) && s[i]) { 
		  s[i]->sample(g,step,ic.get_iwave_fields()[t[i]->rarrindex].gtype,fwd,
			       t[i]->input,w->getStateArray()[t[i]->iwaveindex],
			       t[i]->rarrindex,t[i]->iwaveindex,stream,
			       dryrun,drystr);
		}
	      }
	      if (!dryrun) {
		for (int iv=fd_numsubsteps(ic)-1; iv>-1; iv--) {
		  ic.get_timestep()(w->getRDOMArray(),fwd,iv,fdm);
		  for (size_t k=pow2(order-1); k<pow2(order); k++) {
		    //		    if (w->getStateArray()[k]->printact > 5) 
		    //		      fprintf(stream,"\n*** SYNCH: iwdx = %d\n\n",k);
		    synch(w->getStateArray()[k],fwd,at,iv,ic,stream);
		  }
		}	
#ifdef IWAVE_VERBOSE
		cstr<<"IWaveSim::run adj - first adjoint step: at="<<at;
		cstr.flush();
#endif
	      }
	      else {
		drystr<<"IWaveSim::run adj - first adjoint step: at="<<at;
	      }
	      at--;
#ifdef IWAVE_VERBOSE
	      cstr<<"->"<<at<<"; ref step = "<<it<<"\n";
	      cstr.flush();
#endif
	      if (dryrun) {
		drystr<<"->"<<at<<"; ref step = "<<it<<"\n";
	      }
	    }
	    if (whatodo == ACTION::youturn) {
#ifdef IWAVE_VERBOSE
	      cstr<<"\n";
	      cstr.flush();
#endif
	      if (dryrun) drystr<<"\n";
	      step[g.dim]=at;
	      for (size_t i=0; i<t.size(); i++) {
		// pert sample
		if (t[i]->iwaveindex >= (int)pow2(order-1) && s[i]) {
#ifdef IWAVE_VERBOSE
		  cstr<<"backwards step "<<at<<": sample index=" << i <<"\n";
		  cstr<<" order="<< order <<" input="<< t[i]->input
		      <<" iwdx="<< t[i]->iwaveindex 
		      <<" ridx="<< t[i]->rarrindex <<"\n";
		  if (!(t[i]->input)) {
		    float tmp;
		    RARR const * arr = &(((w->getStateArray()[t[i]->iwaveindex])->model).ld_a._s[t[i]->rarrindex]);
		    ra_a_inner(arr,arr,&tmp);
		    cstr<<" l2 norm-squared of output = "<<tmp<<"\n";
		  }
		  cstr.flush();
#endif

		  s[i]->sample(g,step,ic.get_iwave_fields()[t[i]->rarrindex].gtype,fwd,
			       t[i]->input,w->getStateArray()[t[i]->iwaveindex],
			       t[i]->rarrindex,t[i]->iwaveindex,stream,
			       dryrun,drystr);
		}
	      }
	      if (!dryrun) {
		for (int iv=fd_numsubsteps(ic)-1; iv>-1; iv--) {
		  ic.get_timestep()(w->getRDOMArray(),fwd,iv,fdm);
		  for (size_t k=pow2(order-1); k<pow2(order); k++) {
		    //		    if (w->getStateArray()[k]->printact > 5) 
		    //		      fprintf(stream,"\n*** SYNCH: iwdx = %d\n\n",k);
		    synch(w->getStateArray()[k],fwd,at,iv,ic,stream);
		  }
		}
#ifdef IWAVE_VERBOSE
		cstr<<"IWaveSim::run adj step "<<at;
		cstr.flush();
#endif
	      }
	      else {
		drystr<<"IWaveSim::run adj step "<<at;
	      }
	      at--;
#ifdef IWAVE_VERBOSE
	      cstr<<"->"<<at<<"; ref step ="<<it<<"\n";
	      cstr.flush();
#endif
	      if (dryrun) {
		drystr<<"->"<<at<<"; ref step ="<<it<<"\n";
	      }
	    }
	    if (whatodo == ACTION::restore) {
	      int cp = r.getcheck();
	      it = cplist.at(cp);

#ifdef IWAVE_VERBOSE
	      cstr<<"\nIWaveSim::run adj - restore step "<<it<<" from checkpoint "<<cp<<"\n";
	      cstr.flush();
#endif	      
	      for (int j=0;j<(int)pow2(order-1);j++) {
		int l = 0;
		for (int k=0;k<RDOM_MAX_NARR;k++) {
		  if (fd_isarr(k,w->getStateArray()[0]->model,ic) && fd_isdyn(k,ic)) {		  
		    if (ra_a_copy(&(((w->getRefRDOMArray())[j])->_s[k]),&(cps[cp][j][l]))) {
		      RVLException e;
		      e<<"Error: IWaveSim::run\n";
		      e<<"attempt to restore checkpoint "<<cp
		       <<" encountered error from ra_a_copy\n";
		      e<<"IWAVE index = "<<j<<" checkpoint RARR index = "
		       <<l<<" RDOM RARR index = "<<k<<"\n";
		      throw e;
		    }
#ifdef IWAVE_VERBOSE
		    float tmp=0.0;
		    ra_a_inner(&(((w->getRefRDOMArray())[j])->_s[k]),
			       &(((w->getRefRDOMArray())[j])->_s[k]),
			       &tmp);
		    cstr<<"  level="<<j<<" cpindx="<<l<<" ridx="<<k<<" normsq="<<tmp<<"\n";
		    cstr.flush();
#endif
		    l++;
		  }
		}
	      }

	      if (dryrun) {
		drystr<<"\nIWaveSim::run adj - restored step "<<it<<" from checkpoint "<<cp<<"\n";
	      }
	    }
	    if (whatodo == ACTION::error) {
	      RVLException e;
	      e<<"Error: IWaveSim::run, adjoint mode\n";
	      e<<"  irregular termination of revolve\n";
	      throw e;
	    }
	  }
	  while ((whatodo != ACTION::terminate) && (whatodo != ACTION::error));

#ifdef IWAVE_VERBOSE
	  cstr<<"IWaveSim: exited time loop\n";
	  cstr<<"final sample at="<<at<<"\n";
	  cstr.flush();
#endif
	  if (dryrun) drystr<<"\n";
	  step[g.dim]=at;
	  for (size_t i=0; i<t.size(); i++) {
	    // pert sample
	    if (t[i]->iwaveindex >= (int)pow2(order-1) && s[i]) {
#ifdef IWAVE_VERBOSE
	      cstr<<"begin sample index "<<i<<" on keyword "<<t[i]->keyword<<"\n";
	      cstr.flush();
#endif
	      s[i]->sample(g,step,ic.get_iwave_fields()[t[i]->rarrindex].gtype,fwd,
			   t[i]->input,w->getStateArray()[t[i]->iwaveindex],
			   t[i]->rarrindex,t[i]->iwaveindex,stream,
			   dryrun,drystr);
#ifdef IWAVE_VERBOSE
	      cstr<<"end sample index "<<i<<" on keyword "<<t[i]->keyword<<"\n";
	      cstr.flush();
#endif	      
	    }
	  }
#ifdef IWAVE_VERBOSE
	  cstr<<"revolve report\n";
	  cstr.flush();
#endif
	  // for dry run, dump revolve output to
	  // dry run file. Else to sim output.
	  if (dryrun) drystr<<"\nRevolve report:\n"<<ostr.str()<<"\n";
	  else {
	    if (retrieveGlobalRank()==0)
	      (fprintf(stream,"\nRevolve report:\n%s\n",ostr.str().c_str()));
	  
	  }
	}
	if (g.gdim > g.dim+2) {
	  RVLException e;
	  e<<"Error: IWaveSim::run\n";
	  e<<"  current data model is 1D array of panels\n";
	  e<<"  so simulation grid should have exactly two\n";
	  e<<"  axes beyond the spatial: time and panel\n";
	  e<<"  however g.dim="<<g.dim<<" and g.gdim="<<g.gdim<<"\n";
	  throw e;
	}
	step[g.dim+1]++;
	if ((step[g.dim+1] >= g.axes[g.dim+1].n) &&
	    (panelindex != panelnum-1)) {
	  //	  if (!(next_step(g,step)) && (panelindex != panelnum-1)) {
	  RVLException e;
	  e<<"Error: IWaveSim::run\n";
	  e<<"  step array at last element but \n";
	  e<<"  panelindex = "<<panelindex<<" not =\n";
	  e<<"  panelnum   = "<<panelnum<<"\n";
	  throw e;
	}
	if (dryrun) {
	  for (int i=g.dim+1;i<g.gdim;i++) {
	    drystr<<"NEXT: step["<<i<<"]="<<step[i]<<"\n";
	  }
	}
#ifdef IWAVE_VERBOSE
	cstr<<"IWaveSim: after step update, end of panelindex loop: panelindex="<<panelindex<<"\n";
	for (int idim=0;idim<g.gdim;idim++) {
	  cstr<<"  dim="<<idim<<" step="<<step[idim]<<"\n";
	}
#endif
      }
#ifdef IWAVE_VERBOSE
      cstr<<"IWaveSim::run - exit\n";
      cstr.flush();
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from IWaveSim::run\n";
      throw e;
    }
  }

  ostream & IWaveSim::write(ostream & str) const {
    str<<"IWaveSim ";
    if (fwd) str<<"forward simulator\n";
    else str<<"adjoint simulator\n";
    str<<"  derivative order = "<<order<<"\n";
    str<<"  number of dynamic arrays in each RDOM = "<<ndyn<<"\n";
    int nio = 0;
    for (size_t i=0;i<t.size();i++) {
      if (s[i]) nio++;
    }
    str<<"  number of i/o tasks = "<<nio<<"\n";
    if (!fwd) {
      str<<"  number of checkpoint states = "<<snaps<<"\n";
      str<<"  number of checkpoint RARRs  = "<<narr<<"\n";
    }
    return str;
  } 
  
  void IWaveSim::printgrid(FILE * fp) const {
    fprint_grid(fp, g);
  }

  void IWaveApply(int argc, char ** argv) {
    try {
      // set up environment
      PARARRAY * par = NULL;
      FILE * stream = NULL;
#ifdef IWAVE_APPLY_VERBOSE
      cerr<<"IWaveApply -> IWaveEnvironment\n";
#endif
      IWaveEnvironment(argc, argv, 0, &par, &stream);

#ifdef IWAVE_USE_MPI
      if (retrieveGroupID() == MPI_UNDEFINED) {
#ifdef IWAVE_APPLY_VERBOSE
	fprintf(stderr,"NOTE: idle rank=%d, finalize MPI, cleanup, exit\n",retrieveGlobalRank());
#endif	
	fprintf(stream,"NOTE: finalize MPI, cleanup, exit\n");
	fflush(stream);
      }
      else {
#endif
#ifdef IWAVE_APPLY_VERBOSE
	fprintf(stderr,"NOTE: working rank=%d proceed with sim\n",retrieveGlobalRank());
#endif	
	int dryrun = 0;
	parse(*par,"dryrun",dryrun);
	ofstream * drystr = NULL;
	if (dryrun) {
	  stringstream dryfile;
	  dryfile <<"dryrun";
	  dryfile <<retrieveGroupID();
	  drystr = new ofstream(dryfile.str().c_str());
	}

	// read basic job params from parfile
	// these are mandatory, more or less
	int deriv=0;
	if (!parse(*par, "deriv", deriv)) {
	  RVLException e;
	  e<<"Error: IWaveApply: failed to read parameter deriv\n";
	  e<<"  (int) = order of derivative to compute\n";
	  e<<"  values: 0 = forward map;\n";
	  e<<"          1 = first derivative or adjoint first derivative\n";
	  e<<"          2 = second derivative or adjoint second derivative\n";
	  e<<"          ...\n";
	  e<<"  check parameter file\n";
	  throw e;
	}
#ifdef IWAVE_APPLY_VERBOSE
	cerr<<"IWaveApply: deriv = "<<deriv<<"\n";
#endif

	bool fwd = true;
	int adj = 0;
	if (!parse(*par, "adjoint", adj) && deriv>0) {
	  RVLException e;
	  e<<"Error: IWaveApply: failed to read parameter adj\n";
	  e<<"  required for application of derivative of any positive order\n";
	  e<<"  (int) = [0 = apply derivative, 1 = apply adjoint derivative]\n";
	  e<<"  check parameter file\n";
	  throw e;
	}
	if (adj == 0) fwd = true;
	else fwd = false;
#ifdef IWAVE_APPLY_VERBOSE
	cerr<<"IWaveApply: fwd = "<<fwd<<"\n";
#endif
	int printact = valparse<int>(*par,"printact",0);
#ifdef IWAVE_APPLY_VERBOSE
	cerr<<"IWaveApply: printact = "<<printact<<"\n";
#endif
	IWaveInfo ic;
#ifdef IWAVE_APPLY_VERBOSE      
	cerr<<"IWaveApply -> IWaveSim constructor\n";
	cerr<<"IWaveApply: fwd="<<fwd<<" deriv="<<deriv<<"\n";
#endif
	if (dryrun) {
	  IWaveSim sim(deriv,fwd,*par,stream,ic,printact,dryrun,*drystr);
#ifdef IWAVE_APPLY_VERBOSE
	  cerr<<"IWaveApply -> IWaveSim::run\n";
#endif
	  sim.run();
	}
	else {
	  IWaveSim sim(deriv,fwd,*par,stream,ic,printact,dryrun,cerr);
#ifdef IWAVE_APPLY_VERBOSE
	  cerr<<"IWaveApply -> IWaveSim::run\n";
#endif
	  sim.run();
	}
#ifdef IWAVE_APPLY_VERBOSE
	cerr<<"IWaveApply -> exit\n";
#endif
	if (drystr) {
	  drystr->close();
	  delete drystr;
	}
#ifdef IWAVE_USE_MPI
	/* end nontriv comm branch */
      }
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e<<"\ncalled from IWaveApply\n";
      throw e;
    }
  }

} // end namespace TSOpt

