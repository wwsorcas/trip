#ifndef __TSOPT__SEGY__OPS__
#define __TSOPT__SEGY__OPS__

#include "taperpars.hh"
#include "local.hh"
#include "linop.hh"
#include "mpiserialfo.hh"
#include "functions.hh"

#ifdef IWAVE_USE_MPI
#include "mpisegypp.hh"
#else
#include "segypp.hh"
#endif

#ifdef IWAVE_USE_MPI
using TSOpt::MPISEGYSpace;
typedef TSOpt::MPISEGYSpace gsp;
#else
using TSOpt::SEGYSpace;
typedef TSOpt::SEGYSpace gsp;
#endif

#include "fftops.hh"

namespace TSOpt {
    
  using RVL::BinaryLocalFunctionObject;
  using RVL::TernaryLocalFunctionObject;
  using RVL::LocalDataContainer;


  void conv(int shift, int nout, int nin, int nker,
	    float * restrict out,
	    const float * restrict in,
	    const float * restrict ker,
	    float scal=1.0);

  void corr(int shift, int nout, int nin, int nker,
	    float * restrict out,
	    const float * restrict in,
	    const float * restrict ker,
	    float scal=1.0);



  // performs either convolution (adj == 0) or correlation (adj != 0)
  // possible future optimization: assuming that all calls will involve
  // same tracec geometries, move initializations of header-dep quantities
  // into setadj, make them mutable member data
  class SEGYConvolve: public TernaryLocalFunctionObject<float> {
  private:
    mutable int adj;
  public:
    SEGYConvolve() { setadj(0); }
    void setadj(int newadj) { adj=newadj; }
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> &,
		    LocalDataContainer<float> const &,
    		    LocalDataContainer<float> const &);
    string getName() const { string tmp = "SEGYConvolve"; return tmp; }
  };

  /*
  class SEGYConvolveOp: public RVL::LinearOp<float> {
  private:
    RVL::Space<float> const & dom;
    RVL::Space<float> const & rng;
    TSOpt::SEGYSpace const & ker;
    RVL::Vector<float> const & k;
    
  protected:
    
    void apply(RVL::Vector<float> const & x,
	       RVL::Vector<float> & y) const {
      try {
	SEGYConvolve f();
	RVL::MPISerialFunctionObject<float> mpif(f);	
	y.eval(mpif,x,k);
      }
      catch (RVL::RVLException & e) {
	e<<"\ncalled from SEGYConvolveOp::apply\n";
	throw e;
      }
    }
    
    void applyAdj(RVL::Vector<float> const & x,
		  RVL::Vector<float> & y) const {
      try {
	SEGYConvolve f();
	f.setadj(1);
	RVL::MPISerialFunctionObject<float> mpif(f);		
	y.eval(mpif,x,k);
      }
      catch (RVL::RVLException & e) {
	e<<"\ncalled from SEGYConvolveOp::applyAdj\n";
	throw e;
      }
    }

    RVL::Operator<float> * clone() const { return new SEGYConvolveOp(*this); }
    
  public:
    
    SEGYConvolveOp(RVL::Space<float> const & _dom,
	           RVL::Space<float> const & _rng,
		   std::string kname)
      : dom(_dom), rng(_rng), ker(kname,"none"), k(ker) {
      try {
	segy tr;
	ker.initialize(tr);
	FILE * fp = NULL;
	if (!(fp=iwave_const_fopen(kname.c_str(),"r+",NULL,stderr))) {
	  RVL::RVLException e;
	  e<<"Error: SEGYTraceConvolveOp constructor\n";
	  e<<"  failed to open file = "<<kname<<"\n";
	  throw e;
	}
	if (!fgettr(fp,&(ker.getMetadata()))) {
	  RVL::RVLException e;
	  e<<"Error: SEGYTraceConvolveOp constructor\n";
	  e<<"  failed to read trace from file = "<<kname<<"\n";
	  throw e;
	}
	iwave_fclose(fp);
	
	gsp const & gdom = dynamic_cast<gsp const &>(dom);
	gsp const & grng = dynamic_cast<gsp const &>(rng);
	
	Value val;
	float dt;
	string dtstr="dt";
	
	gethdval(&(ker.getMetadata()),  (char*)(dtstr.c_str()), &val);
	dt=vtof(hdtype((char*)(dtstr.c_str())), val);
	dt*=0.001;
	
	if((fabs(dt-gdom.getDt()) > 0.0001*dt) ||
	   (fabs(dt-grng.getDt()) > 0.0001*dt)) {
	  RVLException e;
	  e<<"Error: SEGYTraceConvolveOp constructor\n";
	  e<<"  input, output, and kernel do not share same time step\n";
	  throw e;
	}
      }
      catch (bad_cast) {
	RVLException e;
	e<<"Error SEGYTraceConvolveOp constructor\n";
	e<<"  either domain or range is not SEGY space\n";
	throw e;
      }

    }
    
  
    SEGYTraceConvolveOp(SEGYTraceConvolveOp const & op)
      : dom(op.dom), rng(op.rng), ker(op.ker) {}
    
    ~SEGYTraceConvolveOp() {}
    
    RVL::Space<float> const & getDomain() const { return dom; }
    RVL::Space<float> const & getRange() const { return rng; }

    ostream & write(ostream & str) const {
      str<<"SEGYTraceConvolveOp\n";
      return str;
    }
  };	
  
  */
  
  class SEGYTraceConvolve: public BinaryLocalFunctionObject<float> {
  private:
    mutable int adj;
    segytrace const & ker;
    SEGYConvolve f;
    SEGYTraceConvolve();
  public:
    SEGYTraceConvolve(segytrace const & _ker)
      : ker(_ker) { f.setadj(0); }
    SEGYTraceConvolve(SEGYTraceConvolve const & g)
      : ker(g.ker), f(g.f) {}
    void setadj(int newadj) { f.setadj(newadj); }
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> & x,
		    LocalDataContainer<float> const & y) {
      f(x,y,ker);
    }
    string getName() const { string tmp = "SEGYTraceConvolve"; return tmp; }
  };

  class SEGYTraceConvolveOp: public RVL::LinearOp<float> {
  private:
    RVL::Space<float> const & dom;
    RVL::Space<float> const & rng;
    segytrace ker;
  protected:
    void apply(RVL::Vector<float> const & x,
	       RVL::Vector<float> & y) const {
      try {
	SEGYTraceConvolve f(ker);
	RVL::MPISerialFunctionObject<float> mpif(f);	
	y.eval(mpif,x);
      }
      catch (RVL::RVLException & e) {
	e<<"\ncalled from SEGYTraceConvolveOp::apply\n";
	throw e;
      }
    }
    void applyAdj(RVL::Vector<float> const & x,
		  RVL::Vector<float> & y) const {
      try {
	SEGYTraceConvolve f(ker);
	f.setadj(1);
	RVL::MPISerialFunctionObject<float> mpif(f);		
	y.eval(mpif,x);
      }
      catch (RVL::RVLException & e) {
	e<<"\ncalled from SEGYTraceConvolveOp::applyAdj\n";
	throw e;
      }
    }

    RVL::Operator<float> * clone() const { return new SEGYTraceConvolveOp(*this); }
    
  public:
    
    SEGYTraceConvolveOp(RVL::Space<float> const & _dom,
			RVL::Space<float> const & _rng,
			std::string kname)
      : dom(_dom), rng(_rng) {
      try {
	segy tr;
	ker.initialize(tr);
	FILE * fp = NULL;
	if (!(fp=iwave_const_fopen(kname.c_str(),"r+",NULL,stderr))) {
	  RVL::RVLException e;
	  e<<"Error: SEGYTraceConvolveOp constructor\n";
	  e<<"  failed to open file = "<<kname<<"\n";
	  throw e;
	}
	if (!fgettr(fp,&(ker.getMetadata()))) {
	  RVL::RVLException e;
	  e<<"Error: SEGYTraceConvolveOp constructor\n";
	  e<<"  failed to read trace from file = "<<kname<<"\n";
	  throw e;
	}
	iwave_fclose(fp);
	
	gsp const & gdom = dynamic_cast<gsp const &>(dom);
	gsp const & grng = dynamic_cast<gsp const &>(rng);
	
	Value val;
	float dt;
	string dtstr="dt";
	
	gethdval(&(ker.getMetadata()),  (char*)(dtstr.c_str()), &val);
	dt=vtof(hdtype((char*)(dtstr.c_str())), val);
	dt*=0.001;
	
	if((fabs(dt-gdom.getDt()) > 0.0001*dt) ||
	   (fabs(dt-grng.getDt()) > 0.0001*dt)) {
	  RVLException e;
	  e<<"Error: SEGYTraceConvolveOp constructor\n";
	  e<<"  input, output, and kernel do not share same time step\n";
	  throw e;
	}
      }
      catch (bad_cast) {
	RVLException e;
	e<<"Error SEGYTraceConvolveOp constructor\n";
	e<<"  either domain or range is not SEGY space\n";
	throw e;
      }

    }
    
    SEGYTraceConvolveOp(SEGYTraceConvolveOp const & op)
      : dom(op.dom), rng(op.rng), ker(op.ker) {}
    
    ~SEGYTraceConvolveOp() {}
    
    RVL::Space<float> const & getDomain() const { return dom; }
    RVL::Space<float> const & getRange() const { return rng; }

    ostream & write(ostream & str) const {
      str<<"SEGYTraceConvolveOp\n";
      return str;
    }
  };	

  class SEGYTaper: public BinaryLocalFunctionObject<float> {
  private:
    TaperPars tap;
  public:
    SEGYTaper(std::string td): tap(td) {
      //      tap.write(cerr);
    }
    SEGYTaper(SEGYTaper const & x): tap(x.tap) {}
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> &,
    		    LocalDataContainer<float> const &);      
    string getName() const { string tmp = "SEGYTaper"; return tmp; }
  };
  
  class SEGYLinMute: public BinaryLocalFunctionObject<float> {
        
  private:
    float s;   // type 0: slope of mute (dt/dx)
    // type 1: minimum value of gx where taper starts (amp=0.)
        
    float tm;  // type 0: mute onset at zero offset, time AFTER first sample
    // type 1: maximum value of gx where taper starts (amp=0.)
        
    float w;   // width of mute zone
        
    int mute_type; // muting type: 0, conic-like (point source); 1, rectangle (plane-wave src)

    int mode; // 0, zero for t<mute; 1, zero for t>mute
        
    float gxbeg;
    float gxend;
  public:
        
    SEGYLinMute(float _s=0.0f, float _tm=0.0f, float _w=0.0f, int _mode=0, int _type = 0)
      : s(_s),tm(_tm),w(_w),mode(_mode),mute_type(_type) {}
        
    SEGYLinMute(SEGYLinMute const & m)
      : s(m.s),tm(m.tm),w(m.w),mode(m.mode),mute_type(m.mute_type),gxbeg(m.gxbeg),gxend(m.gxend) {}
    void set(float _s, float _tm, float _w, int _type = 0){
      mute_type = _type;
      if (!mute_type) {
	s = _s;
	tm = _tm;
	w = _w;
      }
      else {
	s = _s; //gxbeg
	tm = _tm; //gxend
	w = _w;
      }
    }
        
    ~SEGYLinMute() {}
        
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> &,
		    LocalDataContainer<float> const &);
        
    string getName() const { string tmp = "SEGYLinMute"; return tmp; }
        
  };
    
  
  class SEGYFwdDiffIP: public RVL::UnaryLocalFunctionObject<float> {
  private:
        
    int nint;  // number of trace integrations
        
  public:
        
    SEGYFwdDiffIP(int _nint): nint(_nint) {}
        
    SEGYFwdDiffIP(SEGYFwdDiffIP const & m): nint(m.nint) {}
        
    ~SEGYFwdDiffIP() {}
        
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> &);
        
    string getName() const { string tmp = "SEGYFwdDiffIP"; return tmp; }
        
  };

  class SEGYFwdIntIP: public RVL::UnaryLocalFunctionObject<float> {
        
  private:
        
    int nint;  // number of trace integrations
        
  public:
        
    SEGYFwdIntIP(int _nint): nint(_nint) {}
        
    SEGYFwdIntIP(SEGYFwdIntIP const & m): nint(m.nint) {}
        
    ~SEGYFwdIntIP() {}
        
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> &);
        
    string getName() const { string tmp = "SEGYFwdIntIP"; return tmp; }
        
  };
    
  class SEGYAdjIntIP: public RVL::UnaryLocalFunctionObject<float> {
        
  private:
        
    int nint;  // number of trace integrations
        
  public:
        
    SEGYAdjIntIP(int _nint): nint(_nint) {}
        
    SEGYAdjIntIP(SEGYAdjIntIP const & m): nint(m.nint) {}
        
    ~SEGYAdjIntIP() {}
        
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> &);
        
    string getName() const { string tmp = "SEGYAdjIntIP"; return tmp; }
        
  };

  class SEGYFwdInt: public BinaryLocalFunctionObject<float> {
        
  private:
        
    int nint;  // number of trace integrations
        
  public:
        
    SEGYFwdInt(int _nint): nint(_nint) {}
        
    SEGYFwdInt(SEGYFwdInt const & m): nint(m.nint) {}
        
    ~SEGYFwdInt() {}
        
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> & out,
		    LocalDataContainer<float> const & in) {
      try {
	SEGYFwdIntIP iop(nint);
	RVL::RVLCopy<float> cp;
	cp(out,in);
	iop(out);
      }
      catch (RVL::RVLException e) {
	e<<"\ncalled from SEGYAdjInt::operator()\n";
	throw e;
      }
    }      
        
    string getName() const { string tmp = "SEGYFwdInt"; return tmp; }
        
  };
    
  class SEGYAdjInt: public BinaryLocalFunctionObject<float> {
        
  private:
        
    int nint;  // number of trace integrations
        
  public:
        
    SEGYAdjInt(int _nint): nint(_nint) {}
        
    SEGYAdjInt(SEGYAdjInt const & m): nint(m.nint) {}
        
    ~SEGYAdjInt() {}
        
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> & out,
		    LocalDataContainer<float> const & in) {
      try {
	SEGYAdjIntIP iop(nint);
	RVL::RVLCopy<float> cp;
	cp(out,in);
	iop(out);
      }
      catch (RVL::RVLException e) {
	e<<"\ncalled from SEGYAdjInt::operator()\n";
	throw e;
      }
    }
        
    string getName() const { string tmp = "SEGYAdjInt"; return tmp; }
        
  };

  class SEGYTaperMute: public BinaryLocalFunctionObject<float> {
        
  private:
    float s;   // type 0: slope of mute (dt/dx)
    // type 1: minimum value of gx where taper starts (amp=0.)
        
    float tm;  // type 0: mute onset at zero offset, time AFTER first sample
    // type 1: maximum value of gx where taper starts (amp=0.)
        
    float w;   // width of mute zone
        
    int mute_type;     // muting type: 0, conic-like (point source); 1, rectangle (plane-wave src)
        
    float taper_min;   // minimum value of keyword where taper starts
    float taper_max;   // maximum value of keyword where taper ends
    float width;       // width of taper zone
    float tw;          // taper of end time width, unit(ms)
    int taper_type;    // taper type: 0, geophone position; 1, offset

    // taper parameters for source location
    float sx_min;      // minimum value of source location
    float sx_max;      // maximum value of source location
    float sx_width;    // taper width

  public:
        
    SEGYTaperMute(float _s=0.0f, float _tm=0.0f, float _w=0.0f, int _type = 0, 
		  float _taper_min=0.0f, float _taper_max=0.0f, float _width=0.0f, 
		  int _tapertype=0, float _tw=0.0f,
		  float _sxmin=0.0f, float _sxmax=numeric_limits<float>::max(),
		  float _sxw=0.0f)
      : s(_s),tm(_tm),w(_w),mute_type(_type), 
	taper_min(_taper_min), taper_max(_taper_max), 
	width(_width), tw(_tw), taper_type(_tapertype),
	sx_min(_sxmin), sx_max(_sxmax), sx_width(_sxw) {}
        
    SEGYTaperMute(SEGYTaperMute const & m)
      : s(m.s),tm(m.tm),w(m.w),mute_type(m.mute_type),
	taper_min(m.taper_min),taper_max(m.taper_max),
	width(m.width), tw(m.tw), taper_type(m.taper_type),
	sx_min(m.sx_min), sx_max(m.sx_max), sx_width(m.sx_width) {}
 
    void set(float _s, float _tm, float _w, int _type = 0, 
	     float _taper_min =0.0f , float _taper_max = numeric_limits<float>::max(), 
	     float _width=0.0f, int _tapertype=0, float _tw=0.0f,
	     float _sxmin=0.0f, float _sxmax=numeric_limits<float>::max(),
	     float _sxw=0.0f)
    {
      mute_type = _type;
      s  = _s;
      tm = _tm;
      w  = _w;
            
      taper_min = _taper_min;
      taper_max = _taper_max;
      width     = _width;
      taper_type = _tapertype;
            
      tw = _tw;
            
      sx_min = _sxmin;
      sx_max = _sxmax;
      sx_width = _sxw; 
    }
        
    ~SEGYTaperMute() {}
        
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> &,
		    LocalDataContainer<float> const &);
        
    string getName() const { string tmp = "SEGYTaperMute"; return tmp; }
        
  };

  class readTraceSampleFO: public RVL::UnaryLocalFunctionObjectConstEval<float> {
  private:
    LocalDataContainer<float> & v;
    mutable float dtr;
    mutable int n;
    mutable int beg;
    size_t idx;
    bool init;
    readTraceSampleFO();
    readTraceSampleFO(readTraceSampleFO const & f);
  public:
    readTraceSampleFO(LocalDataContainer<float> & _v)
      : v(_v), dtr(0.0), n(0), beg(0), idx(0), init(false) {}
    void operator()(LocalDataContainer<float> const & x);
    std::string getName() const { std::string tmp = "readTraceSampleFO\n"; return tmp; }

  };

  class writeTraceSampleFO: public RVL::UnaryLocalFunctionObject<float> {
  private:
    LocalDataContainer<float> & v;
    mutable float dtr;
    mutable int n;
    mutable int beg;
    size_t idx;
    bool init;
    writeTraceSampleFO();
    writeTraceSampleFO(writeTraceSampleFO const & f);
  public:
    writeTraceSampleFO(LocalDataContainer<float> & _v)
      : v(_v), dtr(0.0), n(0), beg(0), idx(0), init(false) {}
    void operator()(LocalDataContainer<float> & x);
    std::string getName() const { std::string tmp = "readTraceSampleFO\n"; return tmp; }
  };

  /** scale by uniform factor and optionally a header value. If header
      key is default ("") then only uniform factor is applied. Checks
      that key is legal, but not that it is defined - if you see a
      zero section come out, likely the key was not defined in your
      trace headers as gethdval will then return a zero value. */
  class TraceScaleFO: public RVL::BinaryLocalFunctionObject<float> {
  private:
    float fac;
    std::string key;
    float (*f)(float);
    
  public:
    TraceScaleFO(float _fac, std::string _key="", float (*_f)(float) = NULL)
      : fac(_fac), key(_key), f(_f) {
      try {
	if ((key.size() > 0) &&
	    (getindex(key.c_str()) < 0)) {
	  RVL::RVLException e;
	  e<<"Error: TraceScaleFO constructor\n";
	  e<<"  input header string "<<key<<" is not valid SEGY key\n";
	  throw e;
	}
      }
      catch (RVL::RVLException & e) {
	e<<"\ncalled from TraceScaleFO constructor\n";
	throw e;
      }
    }
    TraceScaleFO(TraceScaleFO const & fo)
      : fac(fo.fac), key(fo.key), f(fo.f) {}
    using RVL::BinaryLocalEvaluation<float>::operator();
    void operator()(RVL::LocalDataContainer<float> & xout,
		    RVL::LocalDataContainer<float> const & xin);
    std::string getName() const {
      std::string tmp = "TraceScaleFO"; return tmp;
    }
  };

  /** assumption: traces in each gather increasing gy, increasing gx
      over gathers: increasing sy, increasing sx */
  class SRCellVolFO
    : public RVL::UnaryLocalFunctionObjectScalarRedn<float,float> {
  private:
    mutable bool init;
    mutable float dgx; mutable float dgy;
    mutable float dsx; mutable float dsy;
    mutable float gxlast; mutable float gylast;
    mutable float sxlast; mutable float sylast;
  public:
    SRCellVolFO(float vol=1.0f) 
      : RVL::UnaryLocalFunctionObjectScalarRedn<float,float>(vol),
      init(false), dgx(1.0f), dgy(1.0f), dsx(1.0f), dsy(1.0f),
      gxlast(0.0f), gylast(0.0f), sxlast(0.0f), sylast(0.0f) {}
    SRCellVolFO(const SRCellVolFO & m)
      : RVL::UnaryLocalFunctionObjectScalarRedn<float,float>(m),
      init(m.init), dgx(m.dgx), dgy(m.dgy), dsx(m.dsx), dsy(m.dsy),
      gxlast(m.gxlast), gylast(m.gylast), sxlast(m.sxlast), sylast(m.sylast) {}
    ~SRCellVolFO() {}
    
    // override base class
    void setValue() { RVL::ScalarRedn<float>::setValue(1.0f); }
    using RVL::UnaryLocalConstEval<float>::operator();
    void operator() (RVL::LocalDataContainer<float> const & x);
    std::string getName() const  { std::string tmp="SRCellVolFO"; return tmp; } 
  };
  
  /** 
      FO to shift data from one set of traces to another also scaling by 
      1/(4*z_s*z_r). Fwd assumes z_s, z_r to be read from input trace,
      Adj from output trace. Other trace should have zero src/rec depths -
      checked. The private data factor should be initialized to cell volume 
      for regularly spaced traces.
  */

  class OpZFwdFO: public RVL::BinaryLocalFunctionObject<float> {

  private:

    float fac;
    
  public:
        
    using RVL::LocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> & x,
		    LocalDataContainer<float> const & y);
    std::string getName() const { std::string tmp = "OpZFwdFO"; return tmp; }
  };

  class OpZAdjFO: public RVL::BinaryLocalFunctionObject<float> {
    
  public:
        
    using RVL::LocalEvaluation<float>::operator();    
    void operator()(LocalDataContainer<float> & x,
		    LocalDataContainer<float> const & y);
    std::string getName() const { std::string tmp = "OpZAdjFO"; return tmp; }
  };

  class SEGYTFTScaleFO: public RVL::UnaryLocalFunctionObject<float> {
  private:
    float power;
    int ab;
    float locut;
    float lopas;
    float hipas;
    float hicut;
    
  public:
    SEGYTFTScaleFO(float _power = 0.0f,
		   int _ab = 0,
		   float _locut = 0.0f,
		   float _lopas = 0.0f,
		   float _hipas = std::numeric_limits<float>::max(),
		   float _hicut = std::numeric_limits<float>::max())
      : power(_power), ab(_ab),
	locut(_locut), lopas(_lopas),
	hipas(_hipas), hicut(_hicut) {
      //      cerr<<"SEGYTFTScaleFO:\n";
      //      cerr<<"  power="<<power<<"\n";
      //      cerr<<"  ab="<<ab<<"\n";
      //      cerr<<"  locut="<<locut<<"\n";
      //      cerr<<"  lop="<<lopas<<"\n";
      //      cerr<<"  power="<<power<<"\n";
    }
    SEGYTFTScaleFO(SEGYTFTScaleFO & fo)
      : power(fo.power), ab(fo.ab),
	locut(fo.locut), lopas(fo.lopas),
	hipas(fo.hipas), hicut(fo.hicut) {}
    using RVL::LocalEvaluation<float>::operator();    
    void operator()(LocalDataContainer<float> & x);
    std::string getName() const { std::string tmp = "SEGYTFTScaleFO"; return tmp; }
  };

  class SSEScaleFO: public RVL::BinaryLocalFunctionObject<float> {
  private:
    float shift;
    float alpha;
    float p;
    
  public:
    SSEScaleFO(float _shift, float _alpha, float _p)
      : shift(_shift), alpha(_alpha), p(_p)  {}
    SSEScaleFO(SSEScaleFO const & fo)
      : shift(fo.shift), alpha(fo.alpha), p(fo.p) {}
    using RVL::BinaryLocalEvaluation<float>::operator();
    void operator()(RVL::LocalDataContainer<float> & xout,
		    RVL::LocalDataContainer<float> const & xin);
    std::string getName() const {
      std::string tmp = "SSEScaleFO"; return tmp;
    }
  };
}

#endif
