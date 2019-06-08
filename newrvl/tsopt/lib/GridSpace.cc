#include "GridSpace.hh"

namespace RVL {

  Axis::Axis(size_t _n,
	     float _d,
	     float _o,
	     int _id,
	     std::string _dunit)
    : n(_n), d(_d), o(_o), id(_id), dunit(_dunit) {
    if (d<TOLFAC*std::numeric_limits<float>::epsilon()) {
      RVLException e;
      e<<"Error: Axis constructor\n";
      e<<"  step d="<<d<<" too small\n";
      throw e;
    }
    if (o/d < 0) 
      gs0 = (o/d)-TOLFAC*std::numeric_limits<float>::epsilon();
    else
      gs0 = (o/d)+TOLFAC*std::numeric_limits<float>::epsilon();
    ge0 = gs0+n-1;
    gs=gs0;
    ge=ge0;
    if (abs(o-gs0*d) < TOLFAC*std::numeric_limits<float>::epsilon())
      atype=0;
    else if (abs(o-(0.5+gs0)*d) < TOLFAC*std::numeric_limits<float>::epsilon())
      atype=1;
    else {
      atype=-1; // signifies neither primal nor dual
    }
  }
  Axis::Axis(Axis const & a)
    : n(a.n), d(a.d), o(a.o), id(a.id), dunit(a.dunit), atype(a.atype),
      gs0(a.gs0), ge0(a.ge0), gs(a.gs), ge(a.ge) {}
    
  bool Axis::operator==(Axis const & a) const {
    float tol=std::numeric_limits<float>::epsilon();
    tol *=TOLFAC;
    bool err=true;
    err = err && (n == a.n);
    err = err && (fabs(d-a.d) < tol*min(d,a.d));
    err = err && (fabs(o-a.o) < tol*min(d,a.d));
    err = err && (id == a.id);
    err = err && (dunit == a.dunit); // this should be replaced with an intelligent
    // equivalent-unit test
    return err;
  }

  void Axis::operator=(Axis const & a) {
    n=a.n; d=a.d; o=a.o; id=a.id; atype=a.atype;
    gs0=a.gs0; ge0=a.ge0; gs=a.gs; ge=a.ge;
  }

  void Axis::setSubAxisStartIndex(int _gs) const { 
    if (_gs<gs0) {
      RVLException e;
      e<<"Error: Axis::setSubAxisStartIndex\n";
      e<<"  gs must be >= gs0;\n";
      e<<"  gs="<<_gs<<" gs0="<<gs0<<"\n";
      throw e;
    }
    gs=_gs;
  }
  void Axis::setSubAxisEndIndex(int _ge) const { 
    if (_ge>ge0) {
      RVLException e;
      e<<"Error: Axis::setSubAxisEndIndex\n";
      e<<"  ge must be <= ge0;\n";
      e<<"  ge="<<_ge<<" ge0="<<ge0<<"\n";
      throw e;
    }
    ge=_ge;
  }

  // used in implementation of getMetaSize helper function for ContentPackage class
  size_t Axis::getAxisSize() const {
    size_t axsize =
      sizeof(size_t) +
      sizeof(float) +
      sizeof(float) +
      6*sizeof(int) +
      dunit.size()+1;
    return axsize;
  }

  ostream & Axis::write(ostream & str) const {
    str<<"n="<<n<<" d="<<d<<" o="<<o<<" id="<<id<<" unit="<<dunit<<" atype="<<atype
       <<" gs0="<<gs0<<" ge0="<<ge0<<" gs="<<gs<<" ge="<<ge<<"\n";
    return str;
  }

  /** returns true if 
      1. a.d = b.d
      2. a.o - b.o is divisible by a.d
      These conditions imply that a and b are sub-axes of a global axis
  */
  bool areCompatible(Axis const & a, Axis const & b) {
    if (abs(a.getStep()-b.getStep()) > TOLFAC*std::numeric_limits<float>::epsilon() * a.getStep())
      return false;
    if ((abs(a.getStep()*((int)((b.getOrigin()-a.getOrigin()+0.1*a.getStep())/a.getStep())) - (b.getOrigin()-a.getOrigin()))) > TOL*std::numeric_limits<float>::epsilon() * a.getStep())
      return false;
    return true;
  }
  Grid::Grid(size_t _gdim): gdim(_gdim), init(false) {}
  Grid::Grid(Grid const & g)
    : gdim(g.gdim), init(g.init) {
    for (size_t i=0; i<g.axes.size();i++) {
      Axis a(g.axes[i]);
      axes.push_back(a);
    }
  }
  void Grid::addAxis(Axis const & a) {
    if (axes.size() < gdim) axes.push_back(a);
    else {
      RVLException e;
      e<<"Error: Grid constructor\n";
      e<<"  gdim = "<<gdim<<" axes already assigned\n";
      e<<"  cannat assign more\n";
      throw e;
    }
    //      cerr<<"Grid::addAxis - axes.size = "<<axes.size()<<"\n";
    if (axes.size()==gdim) init=true;
  }
  Axis const & Grid::getAxis(size_t i) const {
    try {
      return axes.at(i);
    }
    catch (out_of_range) {
      RVLException e;
      e<<"Error: Grid::getAxis\n";
      e<<"  index out of range of initialized axes = [0, "<<axes.size()<<"]\n";
      throw e;
    }
  }
  float Grid::getCellVolume() const {
    if (init) {
      float vol=1.0f;
      for (int i=0; i<gdim; i++) vol *= getAxis(i).getStep();
      return vol;
    }
    else {
      RVLException e;
      e<<"Error: Grid::getCellVolume\n";
      e<<"  grid not initialized\n";
      throw e;
    }
  }
  void Grid::getAxisLimits(IPNT gs0, IPNT ge0) const {
    if (init) {
      for (int i=0;i<gdim;i++) {
	gs0[i]=getAxis(i).getAxisStartIndex();
	ge0[i]=getAxis(i).getAxisEndIndex();
      }
    }
    else {
      RVLException e;
      e<<"Error: Grid::getAxisLimits\n";
      e<<"  grid not initialized\n";
      throw e;
    }
  }
  void Grid::getSubAxisLimits(IPNT gs, IPNT ge) const {
    if (init) {
      for (int i=0;i<gdim;i++) {
	gs[i]=getAxis(i).getSubAxisStartIndex();
	ge[i]=getAxis(i).getSubAxisEndIndex();
      }
    }
    else {
      RVLException e;
      e<<"Error: Grid::getSubAxisLimits\n";
      e<<"  grid not initialized\n";
      throw e;
    }
  }
  void Grid::setSubAxisLimits(IPNT gs, IPNT ge) const {
    try {
      if (init) {
	for (int i=0;i<gdim;i++) {
	  getAxis(i).setSubAxisStartIndex(gs[i]);
	  getAxis(i).setSubAxisEndIndex(ge[i]);
	}
      }
      else {
	RVLException e;
	e<<"Error: Grid::getSubAxisLimits\n";
	e<<"  grid not initialized\n";
	throw e;
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from Grid::setSubAxisLimits\n";
      throw e;
    }
  }
    
  bool Grid::operator==(Grid const & g) const {
    if (g.gdim != gdim) return false;
    if (g.axes.size() != axes.size()) return false;
    bool ret=true;
    for (size_t i=0; i<gdim; i++) ret = ret && (g.axes[i] == axes[i]);
    return ret;
  }

  ostream & Grid::write(ostream & str) const {
    str<<"Grid, dimension="<<gdim<<", "<<axes.size()<<" defined\n";
    for (size_t i=0; i<axes.size(); i++) {
      str<<"Axis "<<i<<": ";
      axes[i].write(str);
    }
  }   
 

  /** returns true if all axes are compatible, pairwise */
  bool areCompatible(Grid const & g, Grid const & h) {
    if (g.getDimension() != h.getDimension()) return false;
    bool res=true;
    for (int idim=0; idim<g.getDimension(); idim++)
      res = res && areCompatible(g.getAxis(idim), h.getAxis(idim));
    return res;
  }

  template<> ostream & writeMeta<Grid>(Grid const & g, ostream & e) {
    g.write(e);
  }
  
  template<> size_t getDataSize<Grid>(Grid const & g) {
    if (g.isInit()) {
      size_t nnn = 1;
      for (size_t i=0; i<g.getDimension(); i++) 
	nnn *= g.getAxis(i).getLength();
      return nnn;
    }
    else {
      RVLException e;
      e<<"Error: Grid::getDataSize()\n";
      e<<"  object not initialized\n";
      throw e;
    }
  }

  template<> size_t getMetaSize<Grid>(Grid const & g) {
    if (g.isInit()) {
      size_t nnn = 0;
      for (int i=0; i<g.getDimension(); i++)
	nnn += g.getAxis(i).getAxisSize();
      return nnn;
    }
    else {
      RVLException e;
      e<<"Error: Grid::getMetaSize()\n";
      e<<"  object not initialized\n";
      throw e;
    }
  }

  template<> char * serialize<Axis>(Axis const & a, size_t & len) {
    len = a.getAxisSize();
    size_t loc=0;
    size_t n = a.getLength();
    float d = a.getStep();
    float o = a.getOrigin();
    int id = a.getID();
    std::string dunit = a.getUnit();
    char * cbuf = new char[len];
    memcpy((&(cbuf[loc])),(char *)(&len),sizeof(size_t)); loc += sizeof(size_t);
    memcpy((&(cbuf[loc])),(char *)(&n),sizeof(size_t)); loc += sizeof(size_t);
    memcpy((&(cbuf[loc])),(char *)(&d),sizeof(float)); loc += sizeof(float);
    memcpy((&(cbuf[loc])),(char *)(&o),sizeof(float)); loc += sizeof(float);
    memcpy((&(cbuf[loc])),(char *)(&id),sizeof(int)); loc += sizeof(int);
    memcpy((&(cbuf[loc])),(char *)(dunit.c_str()),dunit.size()+1);
    return cbuf;
  }

  template<> Axis * deserialize<Axis>(char * cbuf, size_t len) {
    if (len < sizeof(size_t)) {
      RVLException e;
      e<<"Error: deserialize<Axis>\n";
      e<<"  2nd param must hold at least one size_t\n";
      e<<"  you passed "<<len<<"\n";
      throw e;
    }
    size_t loc=0;
    size_t alen;
    memcpy((char *)(&alen),(&(cbuf[loc])),sizeof(size_t)); loc += sizeof(size_t);
    if (len < alen) {
      RVLException e;
      e<<"Error: deserialize<Axis>\n";
      e<<"  2nd param must = axis size = "<<alen<<"\n";
      e<<"  you passed "<<len<<"\n";
      throw e;
    }
    size_t n;
    float d;
    float o;
    int id;
    memcpy((char *)(&n),(&(cbuf[loc])),sizeof(size_t)); loc += sizeof(size_t);
    memcpy((char *)(&d),(&(cbuf[loc])),sizeof(float)); loc += sizeof(float);
    memcpy((char *)(&o),(&(cbuf[loc])),sizeof(float)); loc += sizeof(float);
    memcpy((char *)(&id),(&(cbuf[loc])),sizeof(int)); loc += sizeof(int);
    char * tmp = new char[len-loc];
    memcpy(tmp,(&(cbuf[loc])),(len-loc)*sizeof(char));
    std::string dunit=tmp;
    Axis * a = new Axis(n,d,o,id,dunit);
    return a;
  }

  template<> char * serialize<Grid>(Grid const & g, size_t & len) {
    len = getMetaSize<Grid>(g);
    char * c = new char[len];
    size_t loc=0;
    size_t tmp;
    for (size_t i=0; i< g.getDimension(); i++) {
      memcpy(&(c[loc]),serialize<Axis>(g.getAxis(i),tmp),g.getAxis(i).getAxisSize());
      if (tmp != g.getAxis(i).getAxisSize()) {
	RVLException e;
	e<<"Error: serialize<Grid>\n";
	e<<"  wrong number of bytes read on axis "<<i<<"\n";
	throw e;
      }
      loc += g.getAxis(i).getAxisSize();
    }
  }
  
  template<> Grid * deserialize<Grid>(char * cbuf, size_t len) {
    // first check that the length is divisible by axis length, and compute
    // number of axes

  }

  std::shared_ptr<Grid> make_padded_grid(Grid const & phys,
					 std::vector<int> nlsloc,
					 std::vector<int> nrsloc) {
    try {
      auto g = std::make_shared<Grid>(phys.getDimension());
      for (int j=0; j< phys.getDimension(); j++) {
	int n=phys.getAxis(j).getLength() + nlsloc[j] + nrsloc[j];
	float d=phys.getAxis(j).getStep();
	float o=phys.getAxis(j).getOrigin()-nlsloc[j]*(phys.getAxis(j)).getStep();
	int id=phys.getAxis(j).getID();
	Axis a(n,d,o,id);
	g->addAxis(a);
      }
      return g;
    }
    catch (RVLException & e) {
      e<<"\ncalled from make_padded_grid\n";
      throw e;
    }
  }

  GridDC::GridDC(Grid const & g)
    : ContentPackage<float,Grid>() {
    this->initialize(g);
    data1ds=NULL;
    data1dsc=NULL;
    data2d=NULL;
    data2ds=NULL;
    data2dsc=NULL;
    data3d=NULL;
    data3ds=NULL;
    data3dsc=NULL;
    data4d=NULL;
    data5d=NULL;
  }
  GridDC::GridDC(GridDC const & d)
    : ContentPackage<float,Grid>(d) {
    data1ds=NULL;
    data1dsc=NULL;
    data2d=NULL;
    data2ds=NULL;
    data2dsc=NULL;
    data3d=NULL;
    data3ds=NULL;
    data3dsc=NULL;
    data4d=NULL;
    data5d=NULL;
  }
  GridDC::~GridDC() {
    IPNT gs0;
    IPNT ge0;
    Grid const & g = this->getMetadata();
    g.getAxisLimits(gs0,ge0);
    if (data4d) userfree_(data4d);
    if (data5d) userfree_(data5d);
    if (data3d) {
      //cerr<<"3d\n";
      userfree_(data3d); }
    if (data3ds) {
      //cerr<<"3ds\n";
      userfree_(data3ds+gs0[2]);}
    if (data3dsc) {
      //cerr<<"3dsc\n";
      userfree_(data3dsc+gs0[2]);}
    if (data2d) {
      //cerr <<"2d\n";
      userfree_(data2d);}
    if (data2ds) {
      //cerr<<"2ds\n";
      userfree_(data2ds+gs0[1]); }
    if (data2dsc) {
      //cerr<<"2dsc\n";
      userfree_(data2dsc+gs0[1]); }
  }
  
  float * GridDC::getData1D() { return this->getData(); }
  float const * GridDC::getData1D() const { return this->getData1D(); }
  
  float * GridDC::getData1DGlobal() {
    if (!data1ds) {
      Grid const & g = this->getMetadata();
      if (g.getDimension() < 1) {
	RVLException e;
	e<<"Error: GridDC::getData2DGlobal\n";
	e<<"  dimension of grid = "<<g.getDimension()<<" must be at least 1\n";
	throw e;
      }
      IPNT gs0;
      IPNT ge0;
      g.getAxisLimits(gs0,ge0);
      data1ds=this->getData() - gs0[0];
    }
    return data1ds;
  }
  //  float const * GridDC::getData1DGlobal() const { return this->getData1DGlobalPrep(); }
    float const * GridDC::getData1DGlobal() const {
      if (!data1dsc) {
	Grid const & g = this->getMetadata();
	if (g.getDimension() < 1) {
	  RVLException e;
	  e<<"Error: GridDC::getData2DGlobal\n";
	  e<<"  dimension of grid = "<<g.getDimension()<<" must be at least 1\n";
	  throw e;
	}
	IPNT gs0;
	IPNT ge0;
	g.getAxisLimits(gs0,ge0);
	data1dsc=this->getData() - gs0[0];
      }
      return data1dsc;
    }
  
  float ** GridDC::getData2D() {
    if (!data2d) {
      Grid const & g = this->getMetadata();
      if (g.getDimension() < 2) {
	RVLException e;
	e<<"Error: GridDC::getData2D\n";
	e<<"  dimension of grid = "<<g.getDimension()<<" must be at least 2\n";
	throw e;
      }
      size_t n = getDataSize<Grid>(g);
      data2d=dimView(this->getData(),n,g.getAxis(0).getLength());
    }
    return data2d;
  }
  float const ** GridDC::getData2D() const { return this->getData2D();}
  
  float ** GridDC::getData2DGlobal() {
    if (!data2ds) {
      Grid const & g = dynamic_cast<Grid const &>(this->getMetadata());
      if (g.getDimension() < 2) {
	RVLException e;
	e<<"Error: GridDC::getData2DGlobal\n";
	e<<"  dimension of grid = "<<g.getDimension()<<" must be at least 2\n";
	throw e;
      }
      size_t n = getDataSize<Grid>(g);
      data2ds=dimView(this->getData1DGlobal(),n,g.getAxis(0).getLength());
      IPNT gs0;
      IPNT ge0;
      g.getAxisLimits(gs0,ge0);
      data2ds -= gs0[1];
    }
    return data2ds;
  }
  
  float const ** GridDC::getData2DGlobal() const {
    if (!data2dsc) {
      Grid const & g = dynamic_cast<Grid const &>(this->getMetadata());
      if (g.getDimension() < 2) {
	RVLException e;
	e<<"Error: GridDC::getData2DGlobal\n";
	e<<"  dimension of grid = "<<g.getDimension()<<" must be at least 2\n";
	throw e;
      }
      size_t n = getDataSize<Grid>(g);
      data2dsc=dimView(this->getData1DGlobal(),n,g.getAxis(0).getLength());
      IPNT gs0;
      IPNT ge0;
      g.getAxisLimits(gs0,ge0);
      data2dsc -= gs0[1];
    }
    return data2dsc;
  }
  
  float *** GridDC::getData3D() {
    if (!data3d) {
      Grid const & g = this->getMetadata();
      if (g.getDimension() < 3) {
	RVLException e;
	e<<"Error: GridDC::getData3D\n";
	e<<"  dimension of grid = "<<g.getDimension()<<" must be at least 3\n";
	throw e;
      }
      size_t n = getDataSize<Grid>(g)/
	g.getAxis(0).getLength();
      data3d=dimView(this->getData2D(),n,g.getAxis(1).getLength());
    }
    return data3d;
  }
  float const *** GridDC::getData3D() const { return this->getData3D(); }

  float *** GridDC::getData3DGlobal() {
    if (!data3ds) {
      Grid const & g = this->getMetadata();
      if (g.getDimension() < 3) {
	RVLException e;
	e<<"Error: GridDC::getData3DGlobal\n";
	e<<"  dimension of grid = "<<g.getDimension()<<" must be at least 3\n";
	throw e;
      }
      size_t n = getDataSize<Grid>(g)/
	g.getAxis(0).getLength();
      data3ds=dimView(this->getData2DGlobal(),n,g.getAxis(1).getLength());
      IPNT gs0;
      IPNT ge0;
      g.getAxisLimits(gs0,ge0);
      data3ds -= gs0[2];
    }
    return data3ds;
  }
  float const *** GridDC::getData3DGlobal() const {
    if (!data3dsc) {
      Grid const & g = this->getMetadata();
      if (g.getDimension() < 3) {
	RVLException e;
	e<<"Error: GridDC::getData3DGlobal\n";
	e<<"  dimension of grid = "<<g.getDimension()<<" must be at least 3\n";
	throw e;
      }
      size_t n = getDataSize<Grid>(g)/
	g.getAxis(0).getLength();
      data3dsc=dimView(this->getData2DGlobal(),n,g.getAxis(1).getLength());
      IPNT gs0;
      IPNT ge0;
      g.getAxisLimits(gs0,ge0);
      data3dsc -= gs0[2];
    }
    return data3dsc;
  }

  
  float **** GridDC::getData4D() {
    if (!data4d) {
      Grid const & g = this->getMetadata();
      if (g.getDimension() < 4) {
	RVLException e;
	e<<"Error: GridDC::getData4D\n";
	e<<"  dimension of grid = "<<g.getDimension()<<" must be at least 4\n";
	throw e;
      }
      size_t n = getDataSize<Grid>(g)/
	(g.getAxis(1).getLength()*
	 g.getAxis(0).getLength());
      data4d=dimView(this->getData3D(),n,g.getAxis(2).getLength());
    }
    return data4d;
  }
  float const **** GridDC::getData4D() const { return this->getData4D(); }
  
  float ***** GridDC::getData5D() {
    if (!data5d) {
      Grid const & g = this->getMetadata();
      if (g.getDimension() < 5) {
	RVLException e;
	e<<"Error: GridDC::getData5D\n";
	e<<"  dimension of grid = "<<g.getDimension()<<" must be at least 5\n";
	throw e;
      }
      size_t n = getDataSize<Grid>(g)/
	(g.getAxis(1).getLength()*
	 g.getAxis(0).getLength()*
	 g.getAxis(2).getLength());
      data5d=dimView(this->getData4D(),n,g.getAxis(3).getLength());
    }
    return data5d;
  }
  float const ***** GridDC::getData5D() const { return this->getData5D(); }
  
  GridDomain::GridDomain(std::vector<Grid *> ctrlgrids,
			 std::vector<Grid *> stategrids,
			 std::map<std::string,int> _ind)
    : TSSpace<float>(2),
    gsa(stategrids.size()),
    gea(stategrids.size()),
    gsc(stategrids.size()),
    gec(stategrids.size()),
    ind(_ind),
    gdim((*(ctrlgrids[0])).getDimension()),
    dxs((*(ctrlgrids[0])).getDimension()) {
      
    TSSpace<float> csp(ctrlgrids.size());
    for (size_t i=0; i<ctrlgrids.size(); i++) {
      GridSpace gsp(*(ctrlgrids[i]));
      csp.set(gsp,i);
    }
    this->set(csp,0);
    TSSpace<float> ssp(stategrids.size());
    for (size_t i=0; i<stategrids.size(); i++) {
      GridSpace gsp(*(stategrids[i]));
      ssp.set(gsp,i);
      gsp.getGrid().getAxisLimits(gsa[i],gea[i]);
      gsp.getGrid().getSubAxisLimits(gsc[i],gec[i]);
    }
    this->set(ssp,1);
      
    for (int i=0; i< gdim; i++) {
      dxs[i] = (*(ctrlgrids[0])).getAxis(i).getStep();
    }
  }
    
  GridDomain::GridDomain(GridDomain const & gd)
    : TSSpace<float>(gd),
    gsa(gd.gsc.size()),
    gea(gd.gec.size()),
    gsc(gd.gsc.size()),
    gec(gd.gec.size()),
    ind(gd.ind),
    gdim(gd.gdim),
    dxs(gd.gdim) {
    try {
      for (int i=0; i<gd.gsc.size(); i++) {
	IASN(gsc[i],gd.gsc[i]);
	IASN(gec[i],gd.gec[i]);
	IASN(gsa[i],gd.gsa[i]);
	IASN(gea[i],gd.gea[i]);
      }
      TSSpace<float> const & ctrlsp = dynamic_cast<TSSpace<float> const &>(gd[0]);
      GridSpace const & primsp = dynamic_cast<GridSpace const &>(ctrlsp[0]);
	
      for (int i=0; i< gdim; i++) {
	dxs[i] = primsp.getGrid().getAxis(i).getStep();
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridDomain copy constructor\n";
      e<<"  structural error - comp 0 not TSSpace or its comp 0  not GridSpace\n";
      throw e;
    }
  }

  ostream & GridDomain::write(ostream & str) const {
    str<<"\n++++++++++++ begin GridDomain +++++++++++++++\n";
    str<<"GridDomain specialization of TSSpace\n";
    str<<"+++++++++++++++++++++++++++++++++++++++++++++\n";      
    str<<"  all final leaves are GridSpaces\n";
    str<<"  control subspace:\n";
    (*this)[0].write(str);
    str<<"+++++++++++++++++++++++++++++++++++++++++++++\n";
    str<<"  state subspace:\n";
    (*this)[1].write(str);
    str<<"+++++++++++++ end GridDomain ++++++++++++++++\n\n";
    return str;
  }

  std::shared_ptr<Grid> GridfromFile(std::string fname) {
    try {
      PARARRAY * par = ps_new();
      if (ps_createfile(par, fname.c_str())) {
	RVLException e;
	e<<"Error: GridLoad - failed to create param table from file "<<fname<<"\n";
	throw e;
      }
      // figger out dimn - first nx that returns 0
      std::vector<int> nvec;
      int gdim=0;
      std::stringstream nstr;
      nstr<<"n"<<gdim+1;
      int ntmp;	float otmp; float dtmp; int idtmp; std::string unittmp;

      while (((ntmp=valparse<int>(*par,nstr.str(),0)) > 0) && 
	     (gdim < GRID_MAXDIM)) {
	gdim++;
	nstr.str(""); nstr<<"n"<<gdim+1;
      }

      auto g = std::make_shared<Grid>(gdim);
      
      std::stringstream ostr;
      std::stringstream dstr;
      std::stringstream idstr;
      std::stringstream unitstr;
      for (int idim=0; idim<gdim; idim++) {
	nstr.str(""); nstr<<"n"<<idim+1;
	ostr.str(""); ostr<<"o"<<idim+1;
	dstr.str(""); dstr<<"d"<<idim+1;
	idstr.str(""); idstr<<"id"<<idim+1;
	unitstr.str(""); unitstr<<"unit"<<idim+1;
	ntmp=valparse<int>(*par,nstr.str(),0);
	otmp=valparse<float>(*par,ostr.str(),0.0f);
	dtmp=valparse<float>(*par,dstr.str(),1.0f);
	idtmp=valparse<int>(*par,idstr.str(),idim);
	unittmp=valparse<std::string>(*par,unitstr.str(),"");
	Axis a(ntmp,dtmp,otmp,idtmp,unittmp);
	g->addAxis(a);
      }

      return g;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridfromFile\n";
      throw e;
    }
  }
  
  void GridDCIOFO::operator()(LocalDataContainer<float> & x) {
    try {
      GridDC & gx = dynamic_cast<GridDC &>(x);
      if (gx.getGrid() != *g) {
	RVLException e;
	e<<"Error: GridDCIOFO::operator()\n";
	e<<"  target, file grids not same\n";
	e<<"  grid from file: \n";
	g->Writeable::write(e);
	e<<"  target LDC grid:\n";
	gx.getGrid().Writeable::write(e);
	throw e;
      }
      // loademup
      PARARRAY * par = ps_new();
      if (ps_createfile(par, fname.c_str())) { 
	RVLException e;
	e<<"Error: GridDCIOFO::operator()\n";
	e<<"   failed to create param table from file "<<fname<<"\n";
	throw e;
      }
      
      FILE * fp = NULL;
      if (!(fp=fopen((valparse<std::string>(*par,"in")).c_str(),"r+"))) {
	RVLException e;
	e<<"Error: GridDCIOFO::operator()\n";
	e<<"   failed to open data file "<<valparse<std::string>(*par,"in")<<"\n";
	e<<"   read from header file "<<fname<<"\n";
	throw e;
      }		     

      size_t nread = getDataSize<Grid>(*g);
      bool done = true;
      if (load) done = (nread == fread(gx.getData(),sizeof(float),nread,fp));
      else done = (nread == fwrite(gx.getData(),sizeof(float),nread,fp));
      if (!load) fflush(fp);
      fclose(fp);
      if (!done) {
	RVLException e;
	e<<"Error: GridDCIOFO::operator(), load="<<load<<"\n";
	if (load) 
	  e<<"  failed to read "<<nread<<" floats from file "<<fname<<"\n";
	else
	  e<<"  failed to write "<<nread<<" floats to file "<<fname<<"\n";	  
	throw e;
      }
      
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridDCIOFO load = "<<load<<"\n";
      e<<"  DC not a GridDC cannot read to / write from "<<fname<<"\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridDCIOFO::operator()\n";
      throw e;
    }
  }
}
