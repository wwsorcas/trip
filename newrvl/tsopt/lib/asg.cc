#include "asg.hh"
#include "asgfcns.hh"

namespace RVL {

  std::map<std::string,int> make_asg_indices(int dim) {
    std::map<std::string,int> ind;
    ind["bulk"]=0;
    ind["buoy"]=1;
    ind["p0"]=0;
    ind["p1"]=1;
    ind["p2"]=2;
    ind["v0"]=dim;
    ind["v1"]=dim+1;
    ind["v2"]=dim+2;
    return ind;
  }

  std::vector<Grid *> make_asg_ctrllist(Grid const & phys,
					int maxoff) {
    try {
      std::vector<Grid *> glist;
      glist.push_back(new Grid(phys));
      glist.push_back(new Grid(phys));
      return glist;
    }
    catch (RVLException & e) {
      e<<"Error: asg_ctrllist\n";
      throw e;
    }
  }
      
  std::vector<Grid *> make_asg_statelist(Grid const & phys,
					 int maxoff) {
    try {
      std::vector<Grid *> glist;

      // p grids - set up for pml
      for (int i=0; i< phys.getDimension(); i++) {
	Grid * g = new Grid(phys.getDimension());
	for (int j=0; j< phys.getDimension(); j++) {
	  // cerr<<"p: i="<<i<<" j="<<j<<endl;
	  if (i==j) {
	    int n=phys.getAxis(j).getLength()+2*(maxoff-1);
	    float d=phys.getAxis(j).getStep();
	    float o=phys.getAxis(j).getOrigin()-(maxoff-1)*(phys.getAxis(j)).getStep();
	    int id=phys.getAxis(j).getID();
	    // cerr<<"make\n";

	    Axis a(n,d,o,id);
	    //	    a.write(cerr);
	    a.setSubAxisStartIndex(1);
	    a.setSubAxisEndIndex(phys.getAxis(j).getAxisEndIndex()-1);
	    g->addAxis(a);
	  }
	  else {
	    Axis a(phys.getAxis(j).getLength()-2,
		   phys.getAxis(j).getStep(),
		   phys.getAxis(j).getOrigin()+phys.getAxis(j).getStep(),
		   phys.getAxis(j).getID());
	    a.setSubAxisStartIndex(1);
	    a.setSubAxisEndIndex(phys.getAxis(j).getAxisEndIndex()-1);
	    g->addAxis(a);
	  }
	}
	glist.push_back(g);
      }
      // v grids are not
      for (int i=0; i< phys.getDimension(); i++) {
	Grid * g = new Grid(phys.getDimension());
	for (int j=0; j< phys.getDimension(); j++) {
	  // cerr<<"v: i="<<i<<" j="<<j<<endl;
	  if (i==j) {
	    Axis a(phys.getAxis(j).getLength()+2*(maxoff-1)-1,
		   phys.getAxis(j).getStep(),
		   phys.getAxis(j).getOrigin()-(maxoff-0.5)*phys.getAxis(j).getStep(),
		   phys.getAxis(j).getID());
	    a.setSubAxisStartIndex(0);
	    a.setSubAxisEndIndex(phys.getAxis(j).getAxisEndIndex()-1);
	    g->addAxis(a);
	  }
	  else {
	    Axis a(phys.getAxis(j).getLength()-2,
		   phys.getAxis(j).getStep(),
		   phys.getAxis(j).getOrigin()+phys.getAxis(j).getStep(),
		   phys.getAxis(j).getID());
	    a.setSubAxisStartIndex(1);
	    a.setSubAxisEndIndex(phys.getAxis(j).getAxisEndIndex()-1);
	    g->addAxis(a);
	  }
	}
	glist.push_back(g);
      }
      return glist;
    }
    catch (RVLException & e) {
      e<<"Error: asg_statelist\n";
      throw e;
    }
  }

  void asg_pmlaxis(int n0, int nl, int nr,
		   float amp, float dt, int gtype,
		   float ** ep, float ** epp) {

    float t_one = ScalarFieldTraits<float>::One();
    
    if (*ep || *epp) return;
    nl=iwave_max(0,nl);
    nr=iwave_max(0,nr);
    
    //    cerr<<"pmlaxis: gtype="<<gtype<<" n0="<<n0<<" nr="<<nr<<" nl="<<nl<<"\n";
    
    *ep  = (float*)usermalloc_(sizeof(float)*n0);
    *epp = (float*)usermalloc_(sizeof(float)*n0);
    for (int i=0;i<nl;i++) {
      float p = (float(nl-i)-0.5*gtype)/float(nl);
      p = amp*fabs(p*p*p);
      (*ep)[i] = t_one - 0.5*dt*p;
      (*epp)[i] = t_one/(t_one + 0.5*dt*p);
    }
    for (int i=nl; i<n0-nr;i++) {
      (*ep)[i] = t_one;
      (*epp)[i]= t_one;
    }
    for (int i=n0-nr; i<n0;i++) {
      float p =  (float(i-(n0-nr))+1.0-0.5*gtype)/float(nr);
      p = amp*fabs(p*p*p);
      (*ep)[i] = t_one - 0.5*dt*p;
      (*epp)[i] = t_one/(t_one + 0.5*dt*p);
      //        cerr<<"i="<<i<<" p="<<p<<" ep="<<(*ep)[i]<<" epp="<<(*epp)[i]<<endl;
    }
  }

  void ASGapplyFO::operator()(TSDC & y) const {
    try {
      TSDC & cs = dynamic_cast<TSDC &>(y);
      if (cs.getSize() != 2) {
	RVLException e;
	e<<"Error: ASG_applyFO\n";
	e<<"  input not size 2 TSDC\n";
	using RVL::Writeable;
	y.RVL::Writeable::write(e);
	throw e;
      }
      TSDC & ctrl = dynamic_cast<TSDC &>(cs[0]);
      TSDC & state = dynamic_cast<TSDC &>(cs[1]);
      GridDC & bulkdc = dynamic_cast<GridDC &>(ctrl[aux.gdom.get_ind().at("bulk")]);
      GridDC & buoydc = dynamic_cast<GridDC &>(ctrl[aux.gdom.get_ind().at("buoy")]);	
      GridDC & p0dc = dynamic_cast<GridDC &>(state[aux.gdom.get_ind().at("p0")]);
      GridDC & p1dc = dynamic_cast<GridDC &>(state[aux.gdom.get_ind().at("p1")]);      
      GridDC & v0dc = dynamic_cast<GridDC &>(state[aux.gdom.get_ind().at("v0")]);
      GridDC & v1dc = dynamic_cast<GridDC &>(state[aux.gdom.get_ind().at("v1")]);
      
      if (integerstep) {
	for (int i1=(aux.gdom.get_gsc()[aux.gdom.get_ind().at("p0")])[1]+aux.nls[1];
	     i1<=aux.gdom.get_gec()[aux.gdom.get_ind().at("p0")][1]-aux.nrs[1];i1++) {
	  for (int i0=(aux.gdom.get_gsc()[aux.gdom.get_ind().at("p0")])[0]+aux.nls[0];
	       i0<=(aux.gdom.get_gec()[aux.gdom.get_ind().at("p0")])[0]-aux.nrs[0];i0++) {
	    p1dc.getData2DGlobal()[i1][i0] = p0dc.getData2DGlobal()[i1][i0]; 
	  }
	}
	asg_pstep2d(bulkdc.getData2DGlobal(),
		    p0dc.getData2DGlobal(),
		    p1dc.getData2DGlobal(),
		    v0dc.getData2DGlobal(),
		    v1dc.getData2DGlobal(),
		    aux.ep, aux.epp,
		    aux.vdiv,
		    aux.gdom.get_gsc()[aux.gdom.get_ind().at("p0")],
		    aux.gdom.get_gec()[aux.gdom.get_ind().at("p0")],
		    aux.lbc,aux.rbc,
		    aux.maxoff,aux.coeffs);
      }
      else {
	asg_vstep2d(buoydc.getData2DGlobal(),
		    p0dc.getData2DGlobal(),
		    p1dc.getData2DGlobal(),
		    v0dc.getData2DGlobal(),
		    v1dc.getData2DGlobal(),
		    aux.ev, aux.evp,
		    aux.pgrad,
		    aux.gdom.get_gsc()[aux.gdom.get_ind().at("v0")],
		    aux.gdom.get_gec()[aux.gdom.get_ind().at("v0")],
		    aux.gdom.get_gsc()[aux.gdom.get_ind().at("v1")],
		    aux.gdom.get_gec()[aux.gdom.get_ind().at("v1")],
		    aux.lbc,aux.rbc,
		    aux.maxoff,aux.coeffs);
	
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: ASGapplyFO::operator()\n";
      e<<"  either argument not TSDC or components not GridDCs\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from ASGapplyFO::operator()\n";
      throw e;
    }
  }

  void ASGapplyAdjFO::operator()(TSDC & y) const {}

  void ASGapplyTangentFO::operator()(TSDC & y) const {}

  void ASGapplyAdjTangentFO::operator()(TSDC & y) const {}

  void ASGaux::initialize(GridDomain const & gdom,
			  int maxoff, float dt,
			  std::vector<int> nls,
			  std::vector<int> nrs,
			  float amp) {
    try {

      // required for kernel calls - fake call, must be replaced for DD
      for (int i=0; i< gdom.getDimension(); i++) {
	lbc[i]=1; rbc[i]=1;
      }
      
      // fast index workspace - div v
      vdiv_alloc
	= (float *)usermalloc_(((gdom.get_gec()[gdom.get_ind()["p0"]])[0]
				-(gdom.get_gsc()[gdom.get_ind()["p0"]])[0]+1)
			       *sizeof(float));
      vdiv = &(vdiv_alloc[-(gdom.get_gsc()[gdom.get_ind()["p0"]])[0]]);

      // fast index workspace - grad p
      pgrad_alloc = (float **)usermalloc_(gdom.getDimension()*sizeof(float*));
      pgrad = (float **)usermalloc_(gdom.getDimension()*sizeof(float*));
      for (int i=0; i< gdom.getDimension(); i++) {
	std::stringstream vstr;
	vstr<<"v"<<i;
	pgrad_alloc[i] = (float *)usermalloc_(((gdom.get_gec()[gdom.get_ind()[vstr.str()]])[0]
					      -(gdom.get_gsc()[gdom.get_ind()[vstr.str()]])[0]
					      +1)*sizeof(float));
	pgrad[i] = &(pgrad_alloc[i][-(gdom.get_gsc()[gdom.get_ind()[vstr.str()]])[0]]);
      }

      ep_alloc = (float **)usermalloc_(gdom.getDimension()*sizeof(float*));
      ep = (float **)usermalloc_(gdom.getDimension()*sizeof(float*));
      epp_alloc = (float **)usermalloc_(gdom.getDimension()*sizeof(float*));
      epp = (float **)usermalloc_(gdom.getDimension()*sizeof(float*));
      ev_alloc = (float **)usermalloc_(gdom.getDimension()*sizeof(float*));
      ev = (float **)usermalloc_(gdom.getDimension()*sizeof(float*));
      evp_alloc = (float **)usermalloc_(gdom.getDimension()*sizeof(float*));
      evp = (float **)usermalloc_(gdom.getDimension()*sizeof(float*));
	
      for (int idim=0; idim<gdom.getDimension(); idim++) {
	std::stringstream ips;
	ips<<"p"<<idim;
	ep_alloc[idim]=NULL;
	epp_alloc[idim]=NULL;
	asg_pmlaxis((gdom.get_gea()[gdom.get_ind()[ips.str()]])[idim]
		    -(gdom.get_gsa()[gdom.get_ind()[ips.str()]])[idim]
		    +1,
		    nls[idim], nrs[idim],
		    amp, dt, 0,
		    &((ep_alloc)[idim]), &((epp_alloc)[idim]));
	ep[idim]=ep_alloc[idim] - (gdom.get_gsa()[gdom.get_ind()[ips.str()]])[idim];
	epp[idim]=epp_alloc[idim] - (gdom.get_gsa()[gdom.get_ind()[ips.str()]])[idim];

	std::stringstream ivs;
	ivs<<"v"<<idim;
	ev_alloc[idim]=NULL;
	evp_alloc[idim]=NULL;
	asg_pmlaxis((gdom.get_gea()[gdom.get_ind()[ivs.str()]])[idim]
		    -(gdom.get_gsa()[gdom.get_ind()[ivs.str()]])[idim]
		    +1,
		    nls[idim], nrs[idim],
		    amp, dt, 1,
		    &((ev_alloc)[idim]), &((evp_alloc)[idim]));
	ev[idim]=ev_alloc[idim] - (gdom.get_gsa()[gdom.get_ind()[ivs.str()]])[idim];
	evp[idim]=evp_alloc[idim] - (gdom.get_gsa()[gdom.get_ind()[ivs.str()]])[idim];
      }
	
      // fd coeffs
      /* initialize scaled Courant arrays, pml layers */
      coeffs = (float **)usermalloc_(gdom.getDimension()*sizeof(float*));
      for (int idim = 0;idim < gdom.getDimension();idim ++) {
	if (gdom.getSteps()[idim] <= 0.0) {
	  RVLException e;
	  e<<"ERROR: asg_modelinit\n";
	  e<<"  bad input: wrong grid space step, dim="<<idim<<", step="<<gdom.getSteps()[idim]<<"\n";
	  throw e;
	}
	float lam = dt / gdom.getSteps()[idim];
	coeffs[idim] = sgcoeffs(maxoff);
	for (int i=0;i<maxoff;i++) coeffs[idim][i]*=lam;
      }
    }
            
    catch (RVLException & e) {
      e<<"\ncalled from ASGaux::initialize\n";
      throw e;
    }
  }

  ASGaux::ASGaux(GridDomain const & _gdom,
		 int _maxoff, float _dt,
		 std::vector<int> _nls,
		 std::vector<int> _nrs,		    
		 float _amp) 
    : gdom(_gdom), maxoff(_maxoff), dt(_dt),
      nls(_nls), nrs(_nrs),
      amp(_amp)
  { 
    try {
      initialize(gdom,maxoff,dt,nls,nrs,amp);
    }
    catch (RVLException & e) {
      e<<"\ncalled from ASGaux constructor\n";
      throw e;
    }
  }

  ASGaux::ASGaux(ASGaux const & aux) 
    : gdom(aux.gdom), maxoff(aux.maxoff), dt(aux.dt),
      nls(aux.nls), nrs(aux.nrs), amp(aux.amp) {
    try {
      initialize(gdom,maxoff,dt,nls,nrs,amp);
    }
    catch (RVLException & e) {
      e<<"\ncalled from ASGaux copy constructor\n";
      throw e;
    }
  }

  ASGaux::~ASGaux() {
    for (int i=0; i< gdom.getDimension(); i++) {
      userfree_(pgrad_alloc[i]);
      userfree_(ep_alloc[i]);
      userfree_(epp_alloc[i]);
      userfree_(ev_alloc[i]);
      userfree_(evp_alloc[i]);
      userfree_(coeffs[i]);
    }
    userfree_(vdiv_alloc);
    userfree_(pgrad_alloc);
    userfree_(ep);
    userfree_(ep_alloc);
    userfree_(epp);
    userfree_(epp_alloc);
    userfree_(ev);
    userfree_(ev_alloc);
    userfree_(evp);
    userfree_(evp_alloc);
    userfree_(coeffs);
  }

  ostream & ASGaux::write(ostream & str) const {
    str<<"ASG auxiliary data struct \n";
    return str;
  }

  /*
    ASGStateSpace::ASGStateSpace(Grid const & phys, int maxoff)
    : StdProductSpace<float>(2*phys.getDimension()) {
    try {
    vector<Grid *> glist = asg_gridlist(phys,maxoff);
    for (size_t i=0; i<glist.size(); i++) {
    GridSpace sp(*(glist[i]));
    this->set(sp,i);
    }
    }
    catch (RVLException & e) {
    e<<"\ncalled from ASGStateSpace constructor\n";
    throw e;
    }
    }

    // x=[x0,x1]; x0 = control; x1 = state; y = state; usually x1=y.
    void ASG_applyFO::operator()(DataContainer & y,
    DataContainer const & x0,
    DataContainer const & x1) const {
    try {
    // first check to see if this has been called with same input
    // and output if not, set all reference, output dc ptrs to NULL -
    // they will be read from input in that case. Otherwise, skip that step.
    if (refptr != &x0) {
    refdc[0]=NULL; refdc[1]=NULL;
    }
    if (outptr != &y) {
    for (int i=0; i< 2*dim; i++) outdc[i]=NULL;
    }
    // extract control DC
    ProductDataContainer<float> const & px0 =
    dynamic_cast<ProductDataContainer<float> const &>(x0);
    for (int i=0; i< 2; i++) {
    if (!(refdc[i]))
    refdc[i] = dynamic_cast<GridDC const *>(&(px0[i]));
    if (!(refdc[i])) {
    RVLException e;
    e<<"Error: ASG_applyFO::operator()\n";
    e<<"  failed to extract component "<<i<<" as GridDC from input x0\n";
    throw e;
    }
    }

    refdc[0]->getMetaData().getSubAxisLimits(
    if (dim==1) {
    bulk1D = refdc[0]->getData1DGlobal();
    buoy1D = refdc[1]->getData1DGlobal();
    }	
    if (dim==2) {
    bulk2D = refdc[0]->getData2DGlobal();
    buoy2D = refdc[1]->getData2DGlobal();
    }
    else if (dim==3) {
    bulk3D = refdc[0]->getData3DGlobal();
    buoy3D = refdc[1]->getData3DGlobal();
    }
      
    // extract state DC - input and output may be aliased or not
    ProductDataContainer<float> & py =
    dynamic_cast<ProductDataContainer<float> &>(y);
    for (int i=0; i<2*dim; i++) {
    if (!outdc[i])
    outdc[i] = dynamic_cast<GridDC const *>(py[i]);
    if (!outdc[i]) {
    RVLException e;
    e<<"Error: ASG_applyFO::operator()\n";
    e<<"  failed to extract component "<<i<<" as GridDC from input y\n";
    throw e;
    }
    }

      
    if (dim==1) {
    p01D = outdc[0]->getData();
    v01D = outdc[1]->getData();
	
    gsc_p[0] = outdc[0]->getAxes(0).getSubAxisStart();
    gec_p[0] = outdc[0]->getAxes(0).getSubAxisEnd();
    gsc_v[0][0] = outdc[1]->getAxes(0).getSubAxisStart();
    gec_v[0][0] = outdc[1]->getAxes(0).getSubAxisEnd();
    }
    else if (dim==2) {
    p02D = outdc[0]->getData2D();
    p12D = outdc[1]->getData2D();
    v02D = outdc[2]->getData2D();
    v12D = outdc[3]->getData2D();
    IPNT gsc_p, gec_p;
    outdc[0]->getMetadata().getSubAxisLimits(
    gsc_p[0] = outdc[0]->getAxes(0).getSubAxisStart();
    gsc_p[1] = outdc[0]->getAxes(1).getSubAxisStart();
    gec_p[0] = outdc[0]->getAxes(0).getSubAxisEnd();
    gec_p[1] = outdc[0]->getAxes(1).getSubAxisEnd();
    gsc_v[0][0] = outdc[2]->getAxes(0).getSubAxisStart();
    gsc_v[0][1] = outdc[2]->getAxes(1).getSubAxisStart();
    gec_v[0][0] = outdc[2]->getAxes(0).getSubAxisEnd();
    gec_v[0][1] = outdc[2]->getAxes(1).getSubAxisEnd();
    gsc_v[1][0] = outdc[3]->getAxes(0).getSubAxisStart();
    gsc_v[1][1] = outdc[3]->getAxes(1).getSubAxisStart();
    gec_v[1][0] = outdc[3]->getAxes(0).getSubAxisEnd();
    gec_v[1][1] = outdc[3]->getAxes(1).getSubAxisEnd();	  
    }
    else if (dim==3) {
    p03D = outdc[0]->getData3D();
    p13D = outdc[1]->getData3D();
    p23D = outdc[2]->getData3D();
    v03D = outdc[3]->getData3D();
    v13D = outdc[4]->getData3D();
    v23D = outdc[5]->getData3D();
    gsc_p[0] = outdc[0]->getAxes(0).getSubAxisStart();
    gsc_p[1] = outdc[0]->getAxes(1).getSubAxisStart();
    gsc_p[2] = outdc[0]->getAxes(2).getSubAxisStart();
    gec_p[0] = outdc[0]->getAxes(0).getSubAxisEnd();
    gec_p[1] = outdc[0]->getAxes(1).getSubAxisEnd();
    gec_p[2] = outdc[0]->getAxes(2).getSubAxisEnd();
    gsc_v[0][0] = outdc[3]->getAxes(0).getSubAxisStart();
    gsc_v[0][1] = outdc[3]->getAxes(1).getSubAxisStart();
    gsc_v[0][2] = outdc[3]->getAxes(2).getSubAxisStart();
    gec_v[0][0] = outdc[3]->getAxes(0).getSubAxisEnd();
    gec_v[0][1] = outdc[3]->getAxes(1).getSubAxisEnd();
    gec_v[0][2] = outdc[3]->getAxes(2).getSubAxisEnd();
    gsc_v[1][0] = outdc[4]->getAxes(0).getSubAxisStart();
    gsc_v[1][1] = outdc[4]->getAxes(1).getSubAxisStart();
    gsc_v[1][2] = outdc[4]->getAxes(2).getSubAxisStart();
    gec_v[1][0] = outdc[4]->getAxes(0).getSubAxisEnd();
    gec_v[1][1] = outdc[4]->getAxes(1).getSubAxisEnd();	  
    gec_v[1][2] = outdc[4]->getAxes(2).getSubAxisEnd();
    gsc_v[2][0] = outdc[5]->getAxes(0).getSubAxisStart();
    gsc_v[2][1] = outdc[5]->getAxes(1).getSubAxisStart();
    gsc_v[2][2] = outdc[5]->getAxes(2).getSubAxisStart();
    gec_v[2][0] = outdc[5]->getAxes(0).getSubAxisEnd();
    gec_v[2][1] = outdc[5]->getAxes(1).getSubAxisEnd();	  
    gec_v[2][2] = outdc[5]->getAxes(2).getSubAxisEnd();	  	  
    }
    else {
    RVLException e;
    e<<"Error: ASG_applyFO::operator()\n";
    e<<"  dim="<<dim<<" outside of available range {1,2,3}\n";
    throw e;
    } 
    }
    if (&x1 != &y) {
    RVLCopy<float> cp;
    ProductDataContainer<float> const & px1 =
    dynamic_cast<ProductDataContainer<float> const &>(x1);
    for (int i=0; i<2*dim; i++) {
    indc[i] = dynamic_cast<GridDC const *>(px1[i]);
    if (!indc[i]) {
    RVLException e;
    e<<"Error: ASG_applyFO::operator()\n";
    e<<"  failed to extract component "<<i<<" as GridDC from input x1\n";
    throw e;
    }
    }
    // the non-aliased variant is provided only to allow testing via
    // standard LinearOp functions, so can afford some inefficient data
    // motion
    cp(*(outdc[i]),*(indc[i]));
    }
    if (ndim == 2) {
    register float ** restrict bulk2 = xdc[0].getData2DGlobal();
    register ireal ** restrict buoy2 = xdc[1].getData2DGlobal();
    register ireal ** restrict p02 = ydc[0].getData2DGlobal();
    register ireal ** restrict p12 = ydc[1].getData2DGlobal();
    register ireal ** restrict v02 = ydc[2].getData2DGlobal();
    register ireal ** restrict v12 = ydc[3].getData2DGlobal();

    IPNT gsa_p; IPNT gsc_p; IPNT gea_p; IPNT gec_p;
    ydc[0].getAxisIndices(gsa_p,gsc_p,gea_p,gec_p);
    ydc[2].getAxisIndices(gsa_v0,gsc_v0,gea_v0,gec_v0);
    ydc[3].getAxisIndices(gsa_v1,gsc_v1,gea_v1,gec_v1);
      
    if (iv == 0) {
    ireal * sdiv_alloc = (ireal *)usermalloc_((gec_p[0]-gsc_p[0]+1)*sizeof(ireal));
    ireal * sdiv = &(sdiv_alloc[-gsc_p[0]]);
    for (int i1=gsc_p[1]+nls[1]; i1<=gec_p[1]-nrs[1];i1++) {
    for (int i0=gsc_p[0]+nls[0]; i0<=gec_p[0]-nrs[0];i0++) {
    p12[i1][i0]=p02[i1][i0];
    }
    }
    asg_pstep2d(bulk2,
    p02,p12,
    v02,v12,
    ep, epp,
    sdiv,
    gsc_p,gec_p,
    lbc,rbc,
    k,coeffs);
    userfree_(sdiv_alloc);
    }
    if (iv == 1) {
    #ifdef VERBOSE
    fprintf(asgpars->stream,"asg_timestep iv=1 branch\n");
    fflush(asgpars->stream);
    #endif
    ireal * gradp[RARR_MAX_NDIM];
    ireal * gradp_alloc[RARR_MAX_NDIM];
    gradp_alloc[0] = (ireal *)usermalloc_((gec_v[0][0]-gsc_v[0][0]+1)*sizeof(ireal));
    gradp_alloc[1] = (ireal *)usermalloc_((gec_v[1][0]-gsc_v[1][0]+1)*sizeof(ireal));
    gradp[0]=&(gradp_alloc[0][-gsc_v[0][0]]);
    gradp[1]=&(gradp_alloc[1][-gsc_v[1][0]]);
    #ifdef VERBOSE
    fprintf(asgpars->stream,"asg_timestep -> asv_vstep2d\n");
    fflush(asgpars->stream);
    #endif
    asg_vstep2d(buoy2,
    p02,p12,
    v02,v12,
    ev, evp,
    gradp,
    &(gsc_v[0][0]),&(gec_v[0][0]),
    &(gsc_v[1][0]),&(gec_v[1][0]),
    asgpars->lbc,asgpars->rbc,
    asgpars->k,asgpars->coeffs,
    asgpars->stream);
    #ifdef VERBOSE
    fprintf(asgpars->stream,"asg_timestep <- asv_vstep2d\n");
    fflush(asgpars->stream);
    #endif
    userfree_(gradp_alloc[0]);
    userfree_(gradp_alloc[1]);
    }
    }
    if (ndim == 3) {
    register ireal *** restrict bulk3 = (dom[0]->_s)[D_BULK]._s3;
    register ireal *** restrict buoy3 = (dom[0]->_s)[D_BUOY]._s3;
    register ireal *** restrict p03 = (dom[0]->_s)[i_p[0]]._s3;
    register ireal *** restrict p13 = (dom[0]->_s)[i_p[1]]._s3;
    register ireal *** restrict p23 = (dom[0]->_s)[i_p[2]]._s3;
    register ireal *** restrict v03 = (dom[0]->_s)[i_v[0]]._s3;
    register ireal *** restrict v13 = (dom[0]->_s)[i_v[1]]._s3;
    register ireal *** restrict v23 = (dom[0]->_s)[i_v[2]]._s3;

                
    if (iv == 0) {
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
    }
    if (iv == 1) {
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
    }
    }
    }
  */   
}
