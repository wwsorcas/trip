#ifndef __TSOPT_GRIDPP_OPS__
#define __TSOPT_GRIDPP_OPS__

#include "rn.hh"
#include "op.hh"
#include "productspace.hh"
#include "mpiserialfo.hh"
#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#else
#include "gridpp.hh"
#endif

using RVL::ScalarFieldTraits;
using RVL::SpaceTest;
using RVL::Operator;
using RVL::LinearOp;
using RVL::Space;
using RVL::ProductSpace;
using RVL::Vector;
using RVL::Components;
using RVL::ProtectedDivision;
using RVL::RnArray;
using RVL::RVLScale;

/** piecewise linear taper function. no sanity check - calling unit
    responsible for assuring that w > 0 */
template<typename Scalar>
inline Scalar window(Scalar a, Scalar b, Scalar w, Scalar x) {
  return min(min(ScalarFieldTraits<Scalar>::One(),max(ScalarFieldTraits<Scalar>::Zero(),(x-a)/w)),min(ScalarFieldTraits<Scalar>::One(),max(ScalarFieldTraits<Scalar>::Zero(),(b-x)/w)));
}

namespace TSOpt {

#ifdef IWAVE_USE_MPI
  typedef MPIGridSpace myGridSpace;
#else
  typedef GridSpace myGridSpace;
#endif


  using RVL::BinaryLocalFunctionObject;
  using RVL::RVLException;
  using RVL::ContentPackage;
  using RVL::LocalDataContainer;

    
  /** function object for masking a 1D grid function, in the form of
      a CP with Grid metadata. Scalar arrays ffset between input and
      output grids computed by calling unit
  */
    
  class GridMaskFO: public BinaryLocalFunctionObject<float> {
        
  private:
        
    IPNT siw;     // mask start points in gridpoints
    IPNT eiw;     // mask end points in gridpoints
    RPNT width; 
    bool bias;
    GridMaskFO();
        
  public:
        
    GridMaskFO(IPNT const & _siw, IPNT const & _eiw, RPNT const & _width, bool _bias=false)
      : bias(_bias){
      IASN(siw,_siw);
      IASN(eiw,_eiw);
      RASN(width,_width);
    }
        
    GridMaskFO(GridMaskFO const & f)
    : bias(f.bias){
      IASN(siw,f.siw);
      IASN(eiw,f.eiw);
      RASN(width,f.width);
    }
        
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> &,
		    LocalDataContainer<float> const &);
        
    string getName() const { string tmp = "GridMaskFO"; return tmp; }
        
  };
    
  /** Mask operator for grid objects. Apply method outputs
      masked version of background Vector data member: thus
     
      \f$ y = x outside of mask, or 0 inside mask\f$
     
      Derivative and adjoint derivative are implement standard linear injection and
      projection operators.
  */
  class GridMaskOp: public Operator<float> {
        
  private:
        
    Space<float> const & dom;
    Vector<float> const & bg;
    IPNT siw;
    IPNT eiw;
    RPNT width;
    GridMaskOp();
        
  protected:
        
    void apply(Vector<float> const &,
	       Vector<float> &) const;
    void applyDeriv(Vector<float> const &,
		    Vector<float> const &,
		    Vector<float> &) const;
    void applyAdjDeriv(Vector<float> const &,
		       Vector<float> const &,
		       Vector<float> &) const;
    void applyDeriv2(const Vector<float> &,
		     const Vector<float> &,
		     const Vector<float> &,
		     Vector<float> & dy) const { dy.zero(); }
    void applyAdjDeriv2(const Vector<float> &,
			const Vector<float> &,
			const Vector<float> &,
			Vector<float> & dx1) const { dx1.zero(); }
        
    Operator<float> * clone() const { return new GridMaskOp(*this); }
        
  public:
        
    /** main constructor: takes Grid defining increment window, and width
	parameter defining zone of smooth taper - same on all sides.
    */
    GridMaskOp(Space<float> const & _dom,
	       Vector<float> const & _bg,
	       RPNT const & sw = RPNT_0, RPNT const & ew = RPNT_0, 
               RPNT const & width = RPNT_1);
        
    /** Copy constructor - memberwise */
    GridMaskOp(GridMaskOp const & op)
    : dom(op.dom), bg(op.bg) {
      IASN(siw,op.siw);
      IASN(eiw,op.eiw);
      RASN(width,op.width);
    }
        
    ~GridMaskOp() {}
        
    Space<float> const & getDomain() const { return dom; }
    Space<float> const & getRange() const { return bg.getSpace(); }
        
    ostream & write(ostream & str) const;
  };
    
  /** function object for tapering a 1D grid function, in the form of
      a CP with Grid metadata. Scalar arrays ffset between input and
      output grids computed by calling unit
  */

  class GridWindowFO: public BinaryLocalFunctionObject<float> {

  private:

    IPNT iw;     // window width in gridpoints
    bool bias;   // add windowed values if set, else overwrite
    
    GridWindowFO();
    
  public:
    
    GridWindowFO(IPNT const & _iw, bool _bias=false) 
      : bias(_bias) { 
      IASN(iw,_iw);
    }

    GridWindowFO(GridWindowFO const & f) 
    : bias(f.bias) {
      IASN(iw,f.iw);
    }

    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> & x,
		    LocalDataContainer<float> const & y);
    
    string getName() const { string tmp = "GridWindoFO"; return tmp; }
    
  };

  /** Affine window operator for grid objects. Apply method outputs 
      windowed increment of background Vector data member: thus 

      \f$ y = x_{bg} + \phi x\f$

      Derivative and adjoint derivative are independent of
      \f$x_{bg}\f$ and implement standard linear injection and
      projection operators.
  */
  class GridWindowOp: public Operator<float> {

  private:

    std::shared_ptr<Space<float> > dom;
    std::shared_ptr<Vector<float> > bg;
    IPNT iw;
    void initialize(RPNT const w);
    GridWindowOp();

  protected:

    void apply(Vector<float> const &,
	       Vector<float> &) const;
    void applyDeriv(Vector<float> const &,
		    Vector<float> const &,
		    Vector<float> &) const;
    void applyAdjDeriv(Vector<float> const &,
		       Vector<float> const &,
		       Vector<float> &) const;
    void applyDeriv2(const Vector<float> &,
		     const Vector<float> &,
		     const Vector<float> &,
		     Vector<float> & dy) const { dy.zero(); }
    void applyAdjDeriv2(const Vector<float> &,
			const Vector<float> &,
			const Vector<float> &,
			Vector<float> & dx1) const { dx1.zero(); }

    Operator<float> * clone() const { return new GridWindowOp(*this); }
    
  public:
    
    /** main constructor -- old form, deprecated, forces Vector copy
    */
    GridWindowOp(Space<float> const & _dom,
		 Vector<float> const & _bg,
      		 RPNT const sw = RPNT_0);

    /** new style main constructor: 
	- space as const reference, lightweight object to be cloned
	- vector as shared_ptr, heavyweight object to be referenced
	- parameters defining taper on the various axes
    */
    GridWindowOp(Space<float> const & _dom,
		 shared_ptr<RVL::Vector<float> > const _bg,
      		 RPNT const sw = RPNT_0);

    /** Copy constructor - memberwise */
    GridWindowOp(GridWindowOp const & op) 
      : dom(op.dom), bg(op.bg) {
      IASN(iw,op.iw);
    }

    ~GridWindowOp() {}
    
    Space<float> const & getDomain() const { return *dom; }
    Space<float> const & getRange() const { return bg->getSpace(); }

    ostream & write(ostream & str) const;
  };

  class GridFwdDerivFO: public BinaryLocalFunctionObject<float> {

  private:

    int dir;
    float fac; 

    GridFwdDerivFO();

  public:
    GridFwdDerivFO(int _dir, float _fac)
      : dir(_dir), fac(_fac) {}
    GridFwdDerivFO(GridFwdDerivFO const & a)
    : dir(a.dir), fac(a.fac) {}
    
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> & x,
		    LocalDataContainer<float> const & y);
    
    string getName() const { string tmp = "GridFwdDerivFO"; return tmp; }
  };

  class GridAdjDerivFO: public BinaryLocalFunctionObject<float> {

  private:

    int dir;
    float fac; 

    GridAdjDerivFO();

  public:
    GridAdjDerivFO(int _dir, float _fac)
      : dir(_dir), fac(_fac) {}
    GridAdjDerivFO(GridAdjDerivFO const & a)
    : dir(a.dir), fac(a.fac) {}
    
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> & x,
		    LocalDataContainer<float> const & y);
    
    string getName() const { string tmp = "GridAdjDerivFO"; return tmp; }
  };

  class GridDerivOp: public LinearOp<float> {
    
  private:
    
    int dir;
    std::vector<float> fac;
    Space<float> const & dom;

    GridDerivOp();

  protected:

    void apply(Vector<float> const & x,
	       Vector<float> & y) const;
    void applyAdj(Vector<float> const & x,
		  Vector<float> & y) const;

    Operator<float> * clone() const { return new GridDerivOp(*this); }
    
  public:
    
    /** main constructor: takes Grid defining increment window, and width
	parameter defining zone of smooth taper - same on all sides. 
    */
    GridDerivOp(Space<float> const & _dom, int _dir, float scale = REAL_ONE);

    /** Copy constructor - memberwise */
    GridDerivOp(GridDerivOp const & op);

    ~GridDerivOp() {}
    
    Space<float> const & getDomain() const { return dom; }
    Space<float> const & getRange() const { return dom; }

    ostream & write(ostream & str) const;
  };	

  class GridFwdExtendFO: public BinaryLocalFunctionObject<float> {

  private:

    int n_ext;
    bool ext;
    float fac; 

    GridFwdExtendFO();

  public:
    GridFwdExtendFO(int _n_ext, bool _ext, float _fac)
      : n_ext(_n_ext), ext(_ext), fac(_fac) {}
    GridFwdExtendFO(GridFwdExtendFO const & a)
    : n_ext(a.n_ext), ext(a.ext), fac(a.fac) {}
    
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> & x,
		    LocalDataContainer<float> const & y);
    
    string getName() const { string tmp = "GridFwdExtendFO"; return tmp; }
  };

  class GridAdjExtendFO: public BinaryLocalFunctionObject<float> {

  private:

    int n_ext;
    bool ext;
    float fac; 

    GridAdjExtendFO();

  public:

    GridAdjExtendFO(int _n_ext, bool _ext, float _fac)
      : n_ext(_n_ext), ext(_ext), fac(_fac) {}
    GridAdjExtendFO(GridAdjExtendFO const & a)
    : n_ext(a.n_ext), ext(a.ext), fac(a.fac) {}
    
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> & x,
		    LocalDataContainer<float> const & y);
    
    string getName() const { string tmp = "GridAdjExtendFO"; return tmp; }
  };

  class GridExtendOp: public LinearOp<float> {
    
  private:
    
    Space<float> const & dom;
    Space<float> const & rng;
    std::vector<int> n_ext;
    std::vector<bool> ext;
    std::vector<float> fac;

    GridExtendOp();

  protected:

    void apply(Vector<float> const & x,
	       Vector<float> & y) const;
    void applyAdj(Vector<float> const & x,
		  Vector<float> & y) const;

    Operator<float> * clone() const { return new GridExtendOp(*this); }
    
  public:
    
    /** main constructor: takes Grid defining increment window, and width
	parameter defining zone of smooth taper - same on all sides. 
    */
    GridExtendOp(Space<float> const & _dom, Space<float> const & _rng);

    /** Copy constructor - memberwise */
    GridExtendOp(GridExtendOp const & op);

    ~GridExtendOp() {}
    
    Space<float> const & getDomain() const { return dom; }
    Space<float> const & getRange() const { return rng; }

    ostream & write(ostream & str) const;
  };

  /** target = source in overlap. if axis between 0 and dim-1 specified, 
      result extended by constant in axis direction. */
  class GridOverSprayFO: public BinaryLocalFunctionObject<float> {

  private:

    int axis;

    GridOverSprayFO();

  public:
    GridOverSprayFO(int _axis)
      : axis(_axis) {}
    GridOverSprayFO(GridOverSprayFO const & a)
    : axis(a.axis) {}
    
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> & x,
		    LocalDataContainer<float> const & y);
    
    string getName() const { string tmp = "GridOverSprayFO"; return tmp; }
  };

  /** target scaled with factor and power of source, pointwise. Source may be 
      defined on intersection with lower-dimensinoal hyperplane - 
      implicitly extended by constant in orthogonal directions. Version
      of dot-star
  */
      
  class GridScaleFO: public BinaryLocalFunctionObject<float> {

  private:

    float fac;
    float pwr;

    GridScaleFO();
    
  public:

    GridScaleFO(float _fac, float _pwr)
      : fac(_fac), pwr(_pwr) {}

    GridScaleFO(GridScaleFO const & g)
      : fac(g.fac), pwr(g.pwr) {}
    
    void operator()(LocalDataContainer<float> & x,
		    LocalDataContainer<float> const & y);
    std::string getName() const {
      std::string tmp = "GridScaleFO\n"; return tmp; }
    
  };

  class GridFwdZDerivFO: public BinaryLocalFunctionObject<float> {

  private:

    float fac;

    GridFwdZDerivFO();
    
  public:

    GridFwdZDerivFO(float _fac)
      : fac(_fac) {}

    GridFwdZDerivFO(GridFwdZDerivFO const & g)
      : fac(g.fac) {}
    
    void operator()(LocalDataContainer<float> & x,
		    LocalDataContainer<float> const & y);

    std::string getName() const {
      std::string tmp = "GridFwdZDerivFO\n"; return tmp; }
    
  };

  class GridZDerivOp: public LinearOp<float> {
    
  private:
    
    Space<float> const & dom;
    float fac;

    GridZDerivOp();

  protected:

    void apply(Vector<float> const & x,
	       Vector<float> & y) const;
    void applyAdj(Vector<float> const & x,
		  Vector<float> & y) const;

    Operator<float> * clone() const { return new GridZDerivOp(*this); }
    
  public:
    
    /** main constructor: takes Grid defining increment window, and width
	parameter defining zone of smooth taper - same on all sides. 
    */
    GridZDerivOp(Space<float> const & _dom);

    /** Copy constructor - memberwise */
    GridZDerivOp(GridZDerivOp const & op);

    ~GridZDerivOp() {}
    
    Space<float> const & getDomain() const { return dom; }
    Space<float> const & getRange() const { return dom; }

    ostream & write(ostream & str) const;
  };	

  class GridZTaperFO: public RVL::UnaryLocalFunctionObject<float> {
  private:
    int itop;
    int ibot;
    
    GridZTaperFO();

  public:
    GridZTaperFO(int _itop, int _ibot): itop(_itop), ibot(_ibot) {}
    GridZTaperFO(GridZTaperFO const & fo): itop(fo.itop), ibot(fo.ibot) {}
    void operator()(LocalDataContainer<float> & x);

    std::string getName() const {
      std::string tmp = "GridZTaperFO\n"; return tmp; }
    
  };

  class GridZTaperOp: public LinearOp<float> {
    
  private:
    
   Space<float> const & dom;
   int itop;
   int ibot;
   
   GridZTaperOp();
   
 protected:
   
   void apply(Vector<float> const & x,
	      Vector<float> & y) const;
   void applyAdj(Vector<float> const & x,
		 Vector<float> & y) const;
   
   Operator<float> * clone() const { return new GridZTaperOp(*this); }
   
  public:
    
    GridZTaperOp(Space<float> const & _dom, float top, float bot);

    /** Copy constructor - memberwise */
    GridZTaperOp(GridZTaperOp const & op);

    ~GridZTaperOp() {}
    
    Space<float> const & getDomain() const { return dom; }
    Space<float> const & getRange() const { return dom; }

    ostream & write(ostream & str) const;
  };	

  class GridGardnerFO: public BinaryLocalFunctionObject<float> {

  public:

    void operator()(LocalDataContainer<float> & x,
		    LocalDataContainer<float> const & y);

    std::string getName() const {
      std::string tmp = "GridGardnerFO\n"; return tmp; }
    
  };

  class GridCoordScaleFO: public RVL::UnaryLocalFunctionObject<float> {
  private:
    std::vector<float> a;
    float b;
    RPNT o;
    RPNT d;
    
    GridCoordScaleFO();

  public:
    GridCoordScaleFO(std::vector<float> _a,
		     float _b,
		     RPNT const _o,
		     RPNT const _d)
      : a(_a), b(_b) {
      RASN(o,_o);
      RASN(d,_d);
    }
    GridCoordScaleFO(GridCoordScaleFO const & fo)
      : a(fo.a), b(fo.b) {
      RASN(o,fo.o);
      RASN(d,fo.d);
    }
    
    void operator()(LocalDataContainer<float> & x);

    std::string getName() const {
      std::string tmp = "GridCoordScaleFO\n"; return tmp; }
    
  };

  class GridCoordScaleOp: public LinearOp<float> {
    
  private:
    
    Space<float> const & dom;
    std::vector<float> a;
    float b;
    RPNT o;
    RPNT d;
    
    GridCoordScaleOp();
   
 protected:
   
   void apply(Vector<float> const & x,
	      Vector<float> & y) const;
   void applyAdj(Vector<float> const & x,
		 Vector<float> & y) const;
   
   Operator<float> * clone() const { return new GridCoordScaleOp(*this); }
   
  public:
    
    GridCoordScaleOp(Space<float> const & _dom,
		     std::vector<float> _a, float _b);

    /** Copy constructor - memberwise */
    GridCoordScaleOp(GridCoordScaleOp const & op)
      : dom(op.dom), a(op.a), b(op.b) {
      RASN(o,op.o);
      RASN(d,op.d);
    }
    
    ~GridCoordScaleOp() {}
    
    Space<float> const & getDomain() const { return dom; }
    Space<float> const & getRange() const { return dom; }

    ostream & write(ostream & str) const;
  };	

  class GridRadialScaleFO: public RVL::UnaryLocalFunctionObject<float> {
  private:
    float alpha;
    RPNT o;
    RPNT d;
    
    GridRadialScaleFO();

  public:
    GridRadialScaleFO(float _alpha,
		      RPNT const _o,
		      RPNT const _d)
      : alpha(_alpha) {
      RASN(o,_o);
      RASN(d,_d);
    }
    GridRadialScaleFO(GridRadialScaleFO const & fo)
      : alpha(fo.alpha) {
      RASN(o,fo.o);
      RASN(d,fo.d);
    }
    
    void operator()(LocalDataContainer<float> & x);

    std::string getName() const {
      std::string tmp = "GridRadialScaleFO\n"; return tmp; }
    
  };

  class GridRadialScaleOp: public LinearOp<float> {
    
  private:
    
    Space<float> const & dom;
    float alpha;
    RPNT o;
    RPNT d;
    
    GridRadialScaleOp();
   
 protected:
   
   void apply(Vector<float> const & x,
	      Vector<float> & y) const;
   void applyAdj(Vector<float> const & x,
		 Vector<float> & y) const;
   
   Operator<float> * clone() const { return new GridRadialScaleOp(*this); }
   
  public:
    
    GridRadialScaleOp(Space<float> const & _dom,
		      float _alpha);

    /** Copy constructor - memberwise */
    GridRadialScaleOp(GridRadialScaleOp const & op)
      : dom(op.dom), alpha(op.alpha) {
      RASN(o,op.o);
      RASN(d,op.d);
    }
    
    ~GridRadialScaleOp() {}
    
    Space<float> const & getDomain() const { return dom; }
    Space<float> const & getRange() const { return dom; }

    ostream & write(ostream & str) const;
  };



  /** version of box smoothing in 1, 2, and 3 dimensions. constructor
      arguments are dimension and vector of nonnegative smoothing radii. 
  */
  class GridMAFO: public RVL::BinaryLocalFunctionObject<float> {
  private:

    int rep;
    IPNT rad;
    
    GridMAFO();

  public:

    GridMAFO(int _rep, IPNT _rad): rep(_rep) {
      IASN(rad,_rad);
    }
    
    GridMAFO(GridMAFO const & fo) {
      rep=fo.rep;
      IASN(rad,fo.rad);
    }
    
    void operator()(LocalDataContainer<float> & x,
		    LocalDataContainer<float> const & y);

    std::string getName() const {
      std::string tmp = "GridMAFO\n"; return tmp; }

    std::ostream & write(std::ostream & str) const {
      str<<"GridMAFO: function object implementing moving average on grid\n";
      str<<"  Syntax and effect identical to sf_boxsmooth\n";
      str<<"  applicable to functions on grids of dimension <= 3\n";
      for (int i=0; i<3; i++) {
	stringstream ss;
	ss<<"radius on axis "<<i;
	str<<ss.str()<<" = "<<rad[i]<<"\n";
      }
      str<<"repititions of moving average = "<<rep<<"\n";
      return str;
    }
    
  };

  // should be sufficiently general to back up most
  // preconditioner constructions
  class GridOpData {
  private:
    GridOpData();
  public:
    PARARRAY * pars;
    RVL::Space<float> const & dom;
    
    GridOpData(PARARRAY * _pars,
	       RVL::Space<float> const & _dom)
      : pars(_pars), dom(_dom) {}
  
    GridOpData(GridOpData const & d)
      : pars(d.pars), dom(d.dom) {}

    ~GridOpData() {}
  };
  
  class GridMAOpBuilder: public RVL::Builder<RVL::LinearOp<float>, GridOpData> {

  private:

    mutable std::shared_ptr<GridMAFO> f;
    
  public:

    std::shared_ptr<RVL::LinearOp<float> >
    build(shared_ptr<GridOpData> d) const;
  
    std::ostream & write(std::ostream & str) const;
  };


  
  /* lenwork must be > 6*n1*n2+3*max(n2,2*n1)+21 */

  /*
    
  class HelmFO: public BinaryLocalFunctionObject<float> {
    
  private:
    float scale1, scale2;
    float power, datum;
    int DirichletSides;
    IPNT n_arr;
    RPNT d_arr;
    HelmFO();
        
  public:
    HelmFO(IPNT const & _narr,
	   RPNT const & _darr,
	   float _scale1=1.0f,
	   float _scale2=1.0f,
	   float _power=0.0f,
	   float _datum=0.0f,
	   int _DirichletSides=0)
      : scale1(_scale1),scale2(_scale2),power(_power), datum(_datum), DirichletSides(_DirichletSides){
      IASN(n_arr,_narr);
      RASN(d_arr,_darr);
    }
        
    HelmFO(HelmFO const & f)
    : scale1(f.scale1), scale2(f.scale2), power(f.power), datum(f.datum), DirichletSides(f.DirichletSides){
      IASN(n_arr,f.n_arr);
      RASN(d_arr,f.d_arr);
    }
        
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> & x,
		    LocalDataContainer<float> const & y);
        
    string getName() const { string tmp = "HelmFO"; return tmp; }
        
  };
  
  class GridHelmOp: public LinearOp<float> {
  private:
        
    Space<float> const & dom;
    RPNT weights;
    float power, datum;
    int DirichletSides;
        
    // default construction disabled
    GridHelmOp();
        
  protected:
        
    void apply(const Vector<float> & x,
	       Vector<float> & y) const;
        
    void applyAdj(const Vector<float> & x,
		  Vector<float> & y) const;
        
  public:
        
    GridHelmOp(GridHelmOp const & A)
    : dom(A.dom),
      power(A.power), datum(A.datum), DirichletSides(A.DirichletSides) {
      RASN(weights,A.weights);
    }
        
    GridHelmOp(Space<float> const & _dom,
	       RPNT _weights,
	       float _power=0.0f,
	       float _datum=0.0f,
	       int _DirichletSides=0):
      dom(_dom),
      power(_power), datum(_datum), DirichletSides(_DirichletSides){
      try{
	RASN(weights,_weights);
      }
      catch (RVLException & e) {
	e<<"\ncalled from GridHelmOp constructor\n";
	throw e;
      }
    }
        
    ~GridHelmOp() {}
        
    // this class is considered terminal, with no overrides foreseen,
    // so clone method is not virtual
    LinearOp<float> * clone() const { return new GridHelmOp(*this); }
        
    // access to domain, range
    const Space<float> & getDomain() const { return dom; }
    const Space<float> & getRange() const { return dom; }
        
    ostream & write(ostream & str) const {
      str<<"GridHelmOp\n";
      return str;
    }
        
  };
  */
}
#endif
