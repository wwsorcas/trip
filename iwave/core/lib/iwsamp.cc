#include "iwsamp.hh"

//#define IWAVE_VERBOSE

namespace TSOpt {

  IWaveSampler::IWaveSampler(IWAVE * state, string key, PARARRAY & pars, FILE * stream)
    : axes(0), prev_panelindex(-1), tg(NULL), dump_term(0), increment(true) {

    // initialize receiver range indices to global grid
    get_gs(gmin,(state->model).g);
    get_ge(gmax,(state->model).g);
    
    samplekey=key;
    pname="";
    if (parse(pars,key,pname)) {
      string gname = "";

      // extract proto name - use preferentially to actual file name,
      // particularly to handle temp filenames, which do not have suffixes.
      // note that only geometrical info will be extracted here, which can
      // as well be found from the proto files
      const char * rname = NULL;

      // because current implementation of out-of-core classes does all 
      // file handling on rk 0, must extract prototype file name there.
      // this block of code should be abtracted eventually to hide 
      // implementation details against future changes eg. to distrtibuted
      // or network database. As is, will work properly so long as prototype
      // file is used in iwave_fopen on rk=0.

#ifdef IWAVE_USE_MPI
      int cprotolen=0;
      if (retrieveGlobalRank()==0) {
#endif
	if ((rname = iwave_getproto(pname.c_str()))) gname = rname;
	else gname = pname;
#ifdef IWAVE_USE_MPI
	// C string length 
	cprotolen = gname.size()+1;
      }

      MPI_Bcast(&cprotolen,1,MPI_INT,0,retrieveGlobalComm());
      char * cproto = new char[cprotolen];
      if (retrieveGlobalRank()==0) strcpy(cproto,gname.c_str());
      MPI_Bcast(cproto,cprotolen,MPI_CHAR,0,retrieveGlobalComm());
      gname=cproto;
      delete [] cproto;
#endif

      // check for file type
      /*
      size_t pos = gname.find_last_of(".");
      if (pos==std::string::npos || pos >= gname.size()-1) {
	RVLException e;
	e<<"Error: IWaveSampler constructor - filename "<<gname<<" has no suffix\n";
	e<<"gname = "<<gname<<"\n";
	e<<"pname = "<<pname<<"\n";
	e<<"rname = "<<rname<<"\n";
	throw e;
      }
      size_t net = gname.size()-pos-1;
      suffix = gname.substr(pos+1,net);
      */
      suffix = findsuf(gname);
      
#ifdef IWAVE_VERBOSE      
      fprintf(stream,"IWaveSampler: gname=%s pname=%s rank=%d\n",gname.c_str(),pname.c_str(),retrieveGlobalRank());
      fflush(stream);
#endif

      if (suffix=="rsf") {
	grid g;
	if (read_grid(&g,gname.c_str(),stream)) {
	  RVLException e;
	  e<<"Error: IWaveSampler constructor - rsf case, failed to read grid\n";
	  e<<"  key = "<<samplekey<<" proto filename = "<<gname<<"\n";
	  throw e;
	}
	// copy in rigid order - 0 - dim-1 are spatial axes in order (z, x, y); 
	// dim is time, hence skipped; all others are extended - identify by axis id
	// also check that the index mapping is bijective
	for (int i=0; i<g.gdim; i++) {
	  for (int j=0;j<i;j++) {
	    if (g.axes[i].id == g.axes[j].id) {
	      RVLException e;
	      e<<"Error: IWaveSampler constructor\n";
	      e<<"  repetition of axis ids\n";
	      throw e;
	    }
	  }
	}
	for (int k=0; k<g.gdim; k++) {
	  // don't include internal extended axes! they do not
	  // contribute to simulation grid - they are ** internal **
	  if (g.axes[k].id < EXTINT) {
	    //	    cerr<<"  copy axis "<<i<<"\n";
	    axis * a = new axis;
	    copy_axis(a,&(g.axes[k]));
	    axes.push_back(a);
	  }
	  else {
	    //	    cerr<<"  do not copy axis "<<i<<"\n";
	  }
	}

	// determine whether the grid has spatial axes - these would have id's 
	// between 0 and spatial dim
	has_Spatial_Axes = true;
	for (int j=0; j<(state->model).g.dim; j++) 
	  has_Spatial_Axes = has_Spatial_Axes && has_Axis(j);

	//      cerr<<"rsf case: axis vector = \n";
	//      for (int i=0;i<axes.size();i++) {
	//	fprint_axis(stderr,*(axes[i]));
	//      }
	//      cerr<<"exit IWaveSampler constructor, rsf case\n";
	tracestart=0;
	tracestop=0;
      }
      else if (suffix=="su") {

	// by default, not used to sample spatial fields - 
	// only trace time series
	has_Spatial_Axes = false;

#ifdef IWAVE_VERBOSE
	fprintf(stream,"IWaveSamp constructor -> construct_tracegeom, fin=%s\n",pname.c_str());
	fflush(stream);
#endif
	
	/* construct trace geometry object */
	tg = new tracegeom;
	int err=construct_tracegeom(tg,
				    pname.c_str(),
				    (state->model).tsind.dt,
				    SRC_TOL,
				    stream);
  
	if (err) {
	  fprintf(stream,"ERROR: IWaveSampler\n");
	  fprintf(stream,"return from construct_tracegeom with err=%d\n",err);
	  RVLException e;
	  e<<"ERROR: IWaveSampler\n";
	  e<<"  return from construct_tracegeom with err="<<err<<"\n";
	  throw e;
	}
#ifdef IWAVE_VERBOSE
	fprintf(stream,"IWaveSamp constructor <- construct_tracegeom, fin=%s\n",pname.c_str());
	fflush(stream);
#endif
	dump_term = valparse<int>(pars,"dump_term",0);
	taperwidth= valparse<int>(pars,"taperwidth",0);
	timewidth = valparse<int>(pars,"timewidth",0);
        muteslope = valparse<float>(pars,"muteslope",0.0f);
        mutezotime= valparse<float>(pars,"mutezotime",0.0f);
        mutewidth = valparse<float>(pars,"mutewidth",0.0f);
	/* set up axes */
	axis * at = new axis;
	at->n = tg->nt;
	at->o = tg->t0;
	at->d = (state->model).tsind.dt;
	at->id = (state->model).g.dim;
	axes.push_back(at);
	/*
	  if (retrieveGlobalRank()==1) {
	  cerr<<"IWaveSampler - su case, time axis: n="<<at->n<<" d="<<at->d<<" o="<<at->o<<" id="<<at->id<<endl;
	  cerr<<"IWaveSampler - su case, model grid =\n";
	  fprint_grid(stderr,(state->model).g);
	  cerr<<"IWaveSampler - su case, new axis =\n";
	  fprint_axis(stderr,*at);
	  cerr<<"IWaveSampler - su case, recorded axis "<<axes.size()-1<<"=\n";
	  fprint_axis(stderr,*(axes[axes.size()-1]));
	  cerr<<"IWaveSampler - su case, nrec ="<<tg->nrec<<"\n";
	  }
	*/

	tracestart=(int)(floor(((at->o)/(at->d))+0.1));
	tracestop=tracestart + at->n - 1;

	// create record axis
	if (tg->nrec > 1) {
	  axis * ae = new axis;
	  ae->n = tg->nrec;
	  ae->o = 0.0;
	  ae->d = 1.0;
	  ae->id = (state->model).g.dim + 1;
	  axes.push_back(ae);
	}
	
	/* read sampling order */
	sampord=0;
	parse(pars,"sampord",sampord);

	// 2016.10.09 WWS: optionally unset increment, triggering
	// overwrite mode
	std::string ikey = "increment_"+key;
	increment = valparse<bool>(pars,ikey,true);

#ifdef IWAVE_VERBOSE
	fprintf(stream,"IWaveSamp constructor exit\n");
	fflush(stream);
#endif
      }
      else {
	RVLException e;
	e<<"Error: IWaveSampler constructor - suffix "<<suffix<<" represents unknown data structure\n";
	e<<"  rk="<<retrieveGlobalRank()<<"\n";
	iwave_fprintall(stderr);
	throw e;
      }
    }

  }

  IWaveSampler::~IWaveSampler() {
    if (tg) {
      destroy_tracegeom(tg);
      delete tg;
    }
    for (int i=0;i<(int)axes.size(); i++) {
      delete axes[i];
    }
  }

  axis const & IWaveSampler::getAxis(int i) const { 
    try {
      return *(axes.at(i)); 
    }
    catch (out_of_range) {
      RVLException e;
      e<<"Error: IWaveSampler::getAxis\n";
      e<<"  axis index "<<i<<" out of range = [0, "<<axes.size()<<"\n";
      throw e;
    }
  }

  bool IWaveSampler::has_Axis(int id) const {
    bool has = false;
    for (int i=0; i<(int)axes.size(); i++) 
      has = has || (axes[i]->id == id);
    return has;
  }

  axis const * IWaveSampler::getTimeAxis(int dim) const {
    for (int i=0; i<(int)axes.size();i++) {
      if (axes[i]->id == dim) return axes[i];
    }
    return NULL;
  }

  ireal IWaveSampler::getCellVol() const {
    ireal cv = REAL_ONE;
    for (int i=0; i<(int)axes.size();i++) 
      cv *= axes[i]->d;
    return cv;
  }

  ireal IWaveSampler::getRecipCellVol() const {
    ireal cv = REAL_ONE;
    ireal rcv = REAL_ONE;
    for (int i=0; i<(int)axes.size();i++) 
      cv *= axes[i]->d;
    if (ProtectedDivision<ireal>(REAL_ONE,cv,rcv)) {
      RVLException e;
      e<<"Error: IWaveSampler::getRecipCellVol\n";
      e<<"  samplekey = "<<samplekey<<"\n";
      e<<"  zerodivide by cell volume - check data geometry\n";
      throw e;
    } 
    return rcv;
  }
    
  void IWaveSampler::sample(grid g, IPNT step, IPNT gtype, bool fwd, bool input,
			    IWAVE * state, int ridx, int iwdx,
			    FILE * stream,
			    bool dryrun, ostream & drystr) {
    try {

      if (this->axes.size() > 0) {

	bool load = false;       // load data from disk
	bool save = false;       // save data to disk
	bool init = false;       // initialize internal buffers and transfer info, if needed
	bool sample = false;     // sample to internal or external buffer

	// compute istart
	IPNT istart;
	IPNT istop;

	get_gs(istart,g);
	get_ge(istop,g);

#ifdef IWAVE_VERBOSE
	fprintf(stream,"IWaveSampler::sample\n");
	fprintf(stream,"  sim grid: dim=%d gdim=%d\n",g.dim,g.gdim);
	for (int idim=0; idim<g.gdim; idim++)
	  fprintf(stream,"  dim=%d istart=%d istop=%d istep=%d\n",idim,istart[idim],istop[idim],step[idim]);
#endif

	// detect time axis
	axis const * a = this->getTimeAxis(g.dim);

	// compute panelindex = step through external extended indices
	// load - if at BEGIN of all axes that you DON'T HAVE, and at begin of time axis
	// save - if at END of all axes that you DON'T HAVE, and at end of time axis
	// for determining "begin" and "end" of time axis, sense is reversed when
	// fwd=false

	int panelnum=1;
	int panelindex=0;

	// nontrivial extended axes
	if (g.gdim > g.dim+1) {
	  // multipanel movies not permitted
	  /*
	  if (has_Spatial_Axes && a) {
	    RVLException e;
	    e<<"Error: IWaveSampler::sample, samplekey="<<samplekey<<"\n";
	    e<<"  time axis i/o of full spatial grid\n";
	    e<<"  not allowed in conjunction with multiple\n";
	    e<<"  simulations (no multishot movies)\n";
	    throw e;
	  }
	  */
	  for (int i=g.dim+1; i<g.gdim; i++) {
	    panelindex += (step[i]-istart[i])*panelnum;
	    panelnum *= istop[i]-istart[i]+1;
	    //	    cerr<<"axis "<<i<<" has="<<has_Axis(g.axes[i].id) << " input="<<input<<"\n";
	  }
	}
	else {
	  panelnum=1;
	  panelindex=0;
	}

	// compute group panelindex limits
	int first;
	int last;
	calc_group(&first, &last, panelnum,stream);
	if ((panelindex < first) ||
	    (panelindex > last)) {
	  RVLException e;
	  e<<"Error: IWaveSampler::sample\n";
	  e<<"  panelindex "<<panelindex<<" out of range ["<<first<<", "<<last<<"]\n";
	  throw e;
	}
#ifdef IWAVE_VERBOSE
	fprintf(stream,"IWaveSampler::sample: samplekey=%s\n",samplekey.c_str());
	fprintf(stream,"  first=%d last=%d\n",first,last);
	fprintf(stream,"  panelindex=%d\n",panelindex);
	fprintf(stream,"  fwd=%d\n",fwd);
	fprintf(stream,"  input=%d\n",input);
	fflush(stream);
#endif


	// now compute sampler-specific panel index

	int this_panelnum=1;
	int this_panelindex=prev_panelindex;	

	// case 1: extended axes, excluding time axis
	if (dryrun)
	  drystr<<"\nSampler decision loop: samplekey="<<samplekey<<"\n";
	std::vector<int> myaxes(0);
	for (int i=g.dim+1; i<g.gdim; i++) {
	  if (has_Axis(g.axes[i].id)) {
	    myaxes.push_back(i);
	  }
	}
	if (myaxes.size()>0) this_panelindex=0;
	for (size_t i=0; i< myaxes.size(); i++) {
	  this_panelindex += (step[myaxes[i]]-istart[myaxes[i]])*this_panelnum;
	  this_panelnum *= istop[myaxes[i]]-istart[myaxes[i]]+1;
	}
	if (dryrun) 
	  drystr<<myaxes.size()<<" extended axes, sampling enabled\n";
	// no extended axes
	if (myaxes.size()==0) {
	  if ((input && (panelindex==first)) ||
	      (!input && (panelindex==last))) {
	    this_panelindex++;
	    this_panelindex=iwave_min(this_panelindex, this_panelnum-1);
	    if (dryrun) 
	      drystr<<"no extended axes - sampling enabled, update this_panelindex"<<endl;
	  }
	  else {
	    this_panelindex=prev_panelindex;
	    if (dryrun)
	      drystr<<"no extended axes - sampling disabled, do not update this_panelindex"<<endl;
	  }	  
	}

	if (a && has_Spatial_Axes) {
	  // now account for time panels, if time axis is part of this sampler's
	  // grid and it samples a full spatial field - then sampling and load/save
	  // are assumed to be the same
	  this_panelnum = a->n;
	  // compute timeslices for sampling
	  float t = g.axes[g.dim].d*(step[g.dim]-istart[g.dim]) + g.axes[g.dim].o;
	  if (dryrun) drystr<<"  sampler on rarr "<<ridx<<" iwave "<<iwdx<<" time = "<<t<<"\n";
	  // warn if sampling will be forced to duplicate
	  if ((istop[g.dim]-istart[g.dim]+1 < (int)(a->n)) && (step[g.dim]==istart[g.dim])) {
	    cerr<<"NOTE: IWaveSampler::sample\n";
	    cerr<<"  number of time steps = "<<istop[g.dim]-istart[g.dim]+1;
	    cerr<<" less than requested number of time samples = "<<a->n<<"\n";
	    cerr<<"  some time samples will be duplicated\n";
	  }
	  // time before begin of sample range
	  if (t < a->o - 0.01*a->d || t > a->o + (a->n-1+0.01)*(a->d)) this_panelindex = prev_panelindex;
	  else {
	    // compute relative time on axis with positive margin of one-half internal time step
	    float curr_time_rel   = ((t - a->o + 0.5*g.axes[g.dim].d)/(a->d));
	    if (dryrun) drystr<<"curr_time_rel="<<(int)(floor(curr_time_rel))<<"\n";
	    // get panel index
	    this_panelindex = iwave_max(0,iwave_min((int)(floor(curr_time_rel)),this_panelnum-1));
	  }
	}	    

	// time condition for load/save
	// note that input flags input for forward mode - has
	// opposite meaning for reverse mode
	bool buffered = (a && !has_Spatial_Axes);
	bool atbeg = (fwd && (step[g.dim]==istart[g.dim])) ||
	  (!fwd && (step[g.dim]==istop[g.dim]));
	bool atend = (!fwd && (step[g.dim]==istart[g.dim])) ||
	  (fwd && (step[g.dim]==istop[g.dim]));

	if (dryrun && (this_panelindex==prev_panelindex)) {
	  drystr<<"this_panelindex unchanged - no sampling:\n";
	  drystr<<"this_panelindex = "<<this_panelindex<<" this_panelnum="<<this_panelnum<<endl; 
	}
	if (this_panelindex != prev_panelindex) {

	  if (dryrun) {
	    drystr<<"this_panelindex has changed - sampling is possible\n";
	    drystr<<"this_panelnum="<<this_panelnum<<" this_panelindex="<<this_panelindex<<" prev_panelindex="<<prev_panelindex<<"\n";
	    drystr<<"input="<<input<<" buffered="<<buffered<<" atbeg="<<atbeg<<" atend="<<atend<<"\n";
	  }

	  if (buffered) {
	    // for buffered i/o, load only at start, save only at 
	    // end of time axis
	    if (atbeg && input) {
	      load=true;
	      prev_panelindex = this_panelindex;
	    }	   
	    if (atend && !input) {
	      save=true;
	      prev_panelindex = this_panelindex;
	    }
	  }
	  // for unbuffered i/o, load/save whenever panelindex changes, if
	  // you have a time axis, otherwise only at ends of time axis
	  // also zero target array at ends of time axis, not otherwise
	  else if (a) {
	    load = input;
	    save = !input;
	    prev_panelindex = this_panelindex;	    
	  }
	  else {
	    load = input && atbeg;
	    save = !input && atend;
	    if (load || save)  
	      prev_panelindex = this_panelindex;
	    if (load) {
#ifdef IWAVE_VERBOSE
	      fprintf(stream,"zero buffer\n");
#endif
	      if (int err = ra_a_zero(&((state->model).ld_a._s[ridx]))) {
		RVLException e;
		e<<"Error: IWaveSampler::sample from ra_a_zero\n";
		e<<"err="<<err<<"\n";
		throw e;
	      }
	    }
	  }
	}

	// for buffered i/o, sample at every time step, 
	// initialize buffer at begin of time axis
	if (buffered) {
	  sample=true;
	  init = atbeg;
	}
	// otherwise, sampling is synonymous with load/save
	else {
	  sample = (input && load) || (!input && save);
	}

#ifdef IWAVE_VERBOSE
	fprintf(stream,"  panelnum=%d\n",panelnum);
	fprintf(stream,"  panelindex=%d\n",panelindex);
	fprintf(stream,"  this_panelnum=%d\n",this_panelnum);
	fprintf(stream,"  this_panelindex=%d\n",this_panelindex);
	fprintf(stream,"  buffered=%d\n",buffered);
	fprintf(stream,"  atbeg=%d\n",atbeg);
	fprintf(stream,"  atend=%d\n",atend);
	fprintf(stream,"  load=%d\n",load);
	fprintf(stream,"  save=%d\n",save);
	fprintf(stream,"  sample=%d\n",sample);
	fprintf(stream,"  tracestart=%d\n",tracestart);
	fprintf(stream,"  tracestop=%d\n",tracestop);
	fprintf(stream,"  time step=%d\n",step[g.dim]);
	fflush(stream);
#endif
	
	// internal grid info
	int dim;
	IPNT n0;
	IPNT gs0;
	IPNT ge0;
	IPNT n;
	IPNT gs;
	IPNT ge;
	RDOM * dom = &((state->model).ld_a) ;
	if (rd_ndim(dom,ridx,&dim) ||
	    rd_size(dom,ridx,n0) ||
	    rd_gse(dom,ridx,gs0,ge0)) {
	  RVLException e;
	  e<<"Error: IWaveSample::sample\n";
	  e<<"  failure on size query, allocated array\n";
	  e<<"  iwdx="<<iwdx<<" ridx="<<ridx<<"\n";
	  throw e;
	}
	if (rd_size(&((state->model).ld_c),ridx,n) ||
	    rd_gse(&((state->model).ld_c),ridx,gs,ge)) {
	  RVLException e;
	  e<<"Error: IWaveSample::sample\n";
	  e<<"  failure on size query, computational array\n";
	  e<<"  iwdx="<<iwdx<<" ridx="<<ridx<<"\n";
	  throw e;
	}

	// data buffer for current RARR
	float * data = dom->_s[ridx]._s0;
	  
	// scale for adjoint linearized i/o
	//   - input - by cell volume 
	//   - output - by reciprocal cell volume
	// detect linearized i/o by iwdx - if > 0 this i/o op
	// reads or writes a perturbation
	ireal scale = REAL_ONE;
	if (!fwd && input && (iwdx>0)) scale = this->getCellVol();
	if (!fwd && !input && (iwdx>0)) scale = this->getRecipCellVol();
	
	if (sample) {

	  if (dryrun) {
	    
	    drystr<<"\nIWaveSample: samplekey="<<samplekey<<"\n";
	    drystr<<"  panelnum        = "<<panelnum<<"\n";
	    drystr<<"  panelindex      = "<<panelindex<<"\n";
	    drystr<<"  this_panelnum   = "<<this_panelnum<<"\n";
	    drystr<<"  this_panelindex = "<<this_panelindex<<"\n";
	    drystr<<"  input           = "<<input<<"\n";
	    drystr<<"  fwd             = "<<fwd<<"\n";
	    drystr<<"  load            = "<<load<<"\n";
	    drystr<<"  save            = "<<save<<"\n";
	    drystr<<"  init            = "<<init<<"\n";
	    if (fwd) {
	      drystr<<"  fwd time step    = "<<step[g.dim]<<"\n";
	    }
	    else {
	      drystr<<"  bwd time step    = "<<step[g.dim]<<"\n";
	    }
	    if (step[g.dim]==istart[g.dim]) drystr<<"  at start of time loop\n";
	    if (step[g.dim]==istop[g.dim])  drystr<<"  at end of time loop\n";
	    if (load) drystr<<"  read data\n";
	    if (save) drystr<<"  write data\n";
	    for (int i=g.dim+1;i<g.gdim;i++) {
	      if (this->has_Axis(g.axes[i].id)) {
		drystr<<"  extd axis    = "<<i<<" id = "<<g.axes[i].id<<" step = "<<step[i]<<"\n";
	      }
	    }
	    
	    if (this_panelnum > 1) {
	      drystr<<"  this_panelnum     = "<<this_panelnum<<"\n";
	      drystr<<"  this_panelindex   = "<<this_panelindex<<"\n";
	    }
	    drystr<<"  iwave index  = "<<iwdx<<"\n";
	    drystr<<"  rarray index = "<<ridx<<"\n";
	    if (a) drystr <<"  has time axis\n";
	    if (has_Spatial_Axes) drystr<<"  samples spatial fields\n";
	  }
	  
	  else {
	    if (suffix=="rsf") {

	      if (load) {

#ifdef IWAVE_VERBOSE
		fprintf(stream,"IWaveSampler::sample: rsf load on samplekey=%s iwdx=%d ridx=%d\n",this->samplekey.c_str(),iwdx,ridx);
		fprintf(stream,"  rsfread pname=%s gs0[0]=%d gs0[1]=%d gs0[2]=%d\n",pname.c_str(),gs0[0],gs0[1],gs0[2]);
		fprintf(stream,"  rsfread n0[0]=%d n0[1]=%d n0[2]=%d\n",n0[0],n0[1],n0[2]);
		fprintf(stream,"  rsfread scale=%e panelindex=%d\n",scale,this_panelindex);
		fflush(stream);
#endif

		// 03.02.14: read into buffer; axpy via loop - hoist all scaling
		// out of rsfread, which is used in overwrite mode so no need to
		// independently initialize buf
		//
		int err=rsfread(data,gs0,n0,pname.c_str(),scale,stream,this_panelindex);
		if (err) {
		  RVLException e;
		  e<<"Error: IWaveSample::sample, rsf case\n";
		  e<<"  failed to read panel "<<panelindex<<" from "<<pname<<"\n";
		  throw e;
		}	      
		// 03.02.14: temporarily disable scaling for unit conversion - scale
		// only for adjoint input, update array by saxpy
		
	      }  
	      if (save) {
#ifdef IWAVE_VERBOSE
		fprintf(stream,"rsf save on samplekey=%s panelnum=%d iwdx=%d ridx=%d\n",
			this->samplekey.c_str(),this_panelnum,iwdx,ridx);
		fflush(stream);
#endif
		int err=0;

		// stack case 
		if (this_panelnum == 1) {

#ifdef IWAVE_USE_MPI
		  // if no time axis, then local stack
		  MPI_Barrier(MPI_COMM_WORLD);
		  size_t ntot=1;
		  for (int i=0;i<dim;i++) ntot *= n0[i];
		  float * buf = new float[ntot];
		  // reduce along remote comm axis - same domain in every group
		  int err=MPI_Reduce(data,buf,ntot,IWAVE_MPI_REAL,MPI_SUM,0,retrieveRemComm());
		  if (err) {
		    RVLException e;
		    e<<"Error: IWaveSampler::sample from MPI_Reduce\n";
		    e<<"  iwdx = "<<iwdx<<" ridx = "<<ridx<<"\n";
		    throw e;
		  }
		  // on group 0, copy back to alloca domain
		  if (retrieveGroupID()==0) {
		    memcpy(data,buf,ntot*sizeof(ireal));
		  }
		  delete [] buf;
		  if (retrieveGroupID()==0) {
#endif
		  err=rsfwrite(data,gs0,n0,pname.c_str(),scale,stream,0);
		  if (err) {
		    RVLException e;
		    e<<"Error: IWaveSample::sample, rsf single panel / time sample case\n";
		    e<<"  failed to write panel global index = "<<panelindex<<"\n";
		    e<<"  local index = "<<this_panelindex<<" from "<<pname<<"\n";
		    e<<"  return value from rsfwrite = "<<err<<"\n";
		    throw e;
		  }
#ifdef IWAVE_USE_MPI
		  }
#endif
		}
		else {
		  err=rsfwrite(data,gs0,n0,pname.c_str(),scale,stream,this_panelindex);
		  if (err) {
		    RVLException e;
		    e<<"Error: IWaveSample::sample, rsf multiple panel case\n";
		    e<<"  failed to write panel global index = "<<panelindex<<"\n";
		    e<<"  local index = "<<this_panelindex<<" from "<<pname<<"\n";
		    e<<"  return value from rsfwrite = "<<err<<"\n";
		    throw e;
		  }
		  // if more than one panel index for this sampler, then
		  // zero out accumulation in connected rarray, in preparation
		  // for computing the next panel
		  // probably suffices for most instances
		  // BUT don't do this if sampling on the time axis!!!

		  if (!a) {
#ifdef IWAVE_VERBOSE
		    fprintf(stream,"IWaveSample::sample multipanel case\n");
		    fprintf(stream,"  zero buffer for %s\n",
			    this->samplekey.c_str());
#endif
		    ra_a_zero(&(dom->_s[ridx]));
		  }
		}
	      }
	    }

	    else if (suffix=="su") {
	      
	      /* extract grid params - space and time */
	      IPNT nl;                       /* axis lengths, local grid */
	      RPNT ol;                       /* axis origins, local grid */
	      RPNT dl;                       /* axis steps, local grid */
	      IPNT ng;                       /* axis lengths, global grid */
	      RPNT og;                       /* axis origins, global grid */
	      RPNT dg;                       /* axis steps, global grid */
	      IPNT axord;                    /* axis order array */
	      
	      get_n(nl,(state->model).gl);
	      get_o(ol,(state->model).gl);
	      get_d(dl,(state->model).gl);
	      get_n(ng,(state->model).g);
	      get_o(og,(state->model).g);
	      get_d(dg,(state->model).g);
	      get_ord(axord,(state->model).g);

	      // shift axis origin if necessary for staggering - this info
	      // should already be part of grid, which should be part of
	      // rarray.
	      // NOTE: ASSUMED THAT PRIMARY GRID DOES NOT CONTAIN BDRY, so
	      // first staggered grid point is to left of first primary grid
	      // point in simulation grid
	      for (int i=0;i<(state->model).g.dim;i++) {
		ol[i] += 0.5*gtype[i]*dl[i];
		og[i] += 0.5*gtype[i]*dg[i];
	      }
	      // int flag for initializing trace buffer
	      // note that "load" is a subcase of "init" for buffered i/o
	      int initbuf = 0;
	      /*
	      if (!fwd && load) initbuf = -1; // -1 -> adjoint cubic interpolation
	      if ( fwd && load) initbuf =  1; //  1 -> cubic interpolation
	      */
	      // WWS 2016.10.05: always adjoint-interpolate on load
	      if (load && increment) initbuf=-1;
	      if (load && !increment) initbuf=1;
#ifdef IWAVE_VERBOSE
	      if (load) fprintf(stream,"IWaveSampler::sample: su load on samplekey=%s\n",samplekey.c_str());
	      if (save) fprintf(stream,"IWaveSampler::sample: su save on samplekey=%s\n",samplekey.c_str());
#endif

	      if (init) {
		//		cerr<<"samplekey="<<samplekey<<" rk="<<retrieveGlobalRank()<<" init_tracegeom\n";
		/* initialize tracegeom - note data is UNSCALED */
		int err=init_tracegeom(tg,
				       this_panelindex,
				       og,nl,dl,ol,axord,
				       sampord,
				       (state->model).g.dim,
				       initbuf,
				       stream);
		if (err) {
		  RVLException e;
		  e<<"Error: IWaveSampler::sample from init_tracegeom\n";
		  e<<"  err = "<<err<<"\n";
		  throw e;
		}

		// fish out min and max s/r points
		indexrange(tg,gmin,gmax);
		
		if (dump_term) fprint_tracegeom(tg,stream);

	      }

#ifdef IWAVE_VERBOSE
	      fprintf(stream,"IWaveSampler::sample: step=%d tracestart=%d tracestop=%d\n",step[g.dim],tracestart,tracestop);
#endif

	      // sample if in axis range
	      if ((step[g.dim] >= tracestart) &&
		  (step[g.dim] <= tracestop)) {

		// read/write trace switch		
		int traceinput = 0;
		if (input) {
		  traceinput = 1;
		  if (!increment) {
		    traceinput=2;
		  }
		  // overwrite option 2016.10.09 WWS
		  // this assumes that for reference state
		  // that (i) all trace input is RHS (questionable),
		  // and (ii) all trace input should be treated as
		  // defining a numerical delta (safe)
		  if (iwdx==0) {
		    scale=1.0f;
		    for (int kk=0;kk<(state->model).g.dim;kk++) 
		      scale/=(state->model).g.axes[kk].d;
		    if (increment) scale*= (state->model).tsind.rhs;
		  }
		    //		    cerr<<"samplekey="<<samplekey<<" traceinput="<<traceinput<<" rhs="<<(state->model).tsind.rhs<<" iwdx="<<iwdx<<" scale="<<scale<<endl;
		}
		/*
		if (traceinput && (samplekey.compare("data")==0)) {
	            tapermutetraces(tg,step[g.dim] - tracestart, dl[1], taperwidth, timewidth, muteslope, mutezotime, mutewidth);
		}
		*/
		
		sampletraces(tg,
			     sampord,
			     traceinput,
			     step[g.dim] - tracestart,
			     n0, gs0,
			     n, gs,
			     data,
			     scale);

		if (!traceinput && (samplekey.compare("data")==0)) {
	            tapermutetraces(tg,step[g.dim] - tracestart, dl[1], taperwidth, timewidth, muteslope, mutezotime, mutewidth);
		}

	      }

	      if (save) {
#ifdef IWAVE_VERBOSE
		fprintf(stream,"IWaveSampler::sample: samplekey=%s -> writetraces\n",samplekey.c_str());
#endif
		int err=writetraces(tg,dg,og,stream);
		if (err) {
		  RVLException e;
		  e<<"Error: IWaveSampler::sample from writetraces\n";
		  e<<"  err = "<<err<<"\n";
		}
	      }
	    }
	  }
	}

      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from IWaveSampler::sample\n";
      throw e;
    }
  }
    
}

