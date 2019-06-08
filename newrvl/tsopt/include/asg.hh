#include "GridSpace.hh"

namespace RVL {

  float * sgcoeffs(int k);
  std::vector<Grid *> make_asg_ctrllist(Grid const & phys, int maxoff);
  std::vector<Grid *> make_asg_statelist(Grid const & phys, int maxoff);
  std::map<std::string,int> make_asg_indices(int dim);
  void asg_pmlaxis(int n0, int nl, int nr,
		   float amp, float dt, int gtype,
		   float ** ep, float ** epp);
  float * sgcoeffs(int maxoff);
  
  class ASGaux: public Writeable {
  private:
    void initialize(GridDomain const & gdom,
		    int maxoff, float dt,
		    std::vector<int> nls,
		    std::vector<int> nrs,		    
		    float amp);
    ASGaux();
  public:
    GridDomain const & gdom;
    float dt; // required for copy construct
    float amp; // required for copy construct
    std::vector<int> nls;
    std::vector<int> nrs;
    IPNT lbc; // required as int * in call to kernel fcns
    IPNT rbc;
    float * vdiv_alloc;
    float * vdiv;
    float ** pgrad_alloc;
    float ** pgrad;
    float ** ep;
    float ** epp;
    float ** ev;
    float ** evp;
    float ** ep_alloc;
    float ** epp_alloc;
    float ** ev_alloc;
    float ** evp_alloc;
    float ** coeffs;
    int maxoff;

    ASGaux(GridDomain const & gdom,
	   int maxoff, float dt,
	   std::vector<int> nls,
	   std::vector<int> nrs,
	   float amp);
    ASGaux(ASGaux const & aux);
    ~ASGaux();

    ostream & write(ostream & str) const;
    
  };
    
  class ASGapplyFO: public TSFO {
  private:
    ASGaux const & aux;
    bool const & integerstep;
  public:
    ASGapplyFO(ASGaux const & _aux,
	       bool const & _integerstep)
      : aux(_aux),
	integerstep(_integerstep) {}
    ASGapplyFO(ASGapplyFO const & f)
      : aux(f.aux), integerstep(f.integerstep) {}
    ~ASGapplyFO() {}

    void operator()(TSDC & y) const;

    std::string getName() const { std::string tmp = "ASGapplyFO"; return tmp;}
  };
    
  class ASGapplyAdjFO: public TSFO {
  private:
    ASGaux const & aux;
    bool const & integerstep;
  public:
    ASGapplyAdjFO(ASGaux const & _aux,
	       bool const & _integerstep)
      : aux(_aux),
	integerstep(_integerstep) {}
    ASGapplyAdjFO(ASGapplyAdjFO const & f)
      : aux(f.aux), integerstep(f.integerstep) {}
    ~ASGapplyAdjFO() {}

    void operator()(TSDC & y) const;

    std::string getName() const { std::string tmp = "ASGapplyAdjFO"; return tmp;}
  };
  
  class ASGapplyTangentFO: public TSFO {
  private:
    ASGaux const & aux;
    bool const & integerstep;
  public:
    ASGapplyTangentFO(ASGaux const & _aux,
	       bool const & _integerstep)
      : aux(_aux),
	integerstep(_integerstep) {}
    ASGapplyTangentFO(ASGapplyTangentFO const & f)
      : aux(f.aux), integerstep(f.integerstep) {}
    ~ASGapplyTangentFO() {}

    void operator()(TSDC & y) const;

    std::string getName() const { std::string tmp = "ASGapplyTangentFO"; return tmp;}
  };

  class ASGapplyAdjTangentFO: public TSFO{
  private:
    ASGaux const & aux;
    bool const & integerstep;
  public:
    ASGapplyAdjTangentFO(ASGaux const & _aux,
	       bool const & _integerstep)
      : aux(_aux),
	integerstep(_integerstep) {}
    ASGapplyAdjTangentFO(ASGapplyAdjTangentFO const & f)
      : aux(f.aux), integerstep(f.integerstep) {}
    ~ASGapplyAdjTangentFO() {}

    void operator()(TSDC & y) const;

    std::string getName() const { std::string tmp = "ASGapplyAdjTangentFO"; return tmp;}
  };
  
  class ASGStep: public TSStep<float> {

  private:

    std::shared_ptr<Grid> gpad; // made in constructor
    std::vector<Grid *> ctrllist;
    std::vector<Grid *> statelist;
    std::map<std::string,int> indexmap;
    GridDomain gdom; // made in constructor
    float dt; // copied in constructor
    float tmin;
    float amp; // copied in constructor
    std::vector<int> nls; // copied in constructor
    std::vector<int> nrs; // copied in constructor

    ASGaux aux;
    
    mutable bool integerstep;
    
    mutable ASGapplyFO app;
    mutable ASGapplyAdjFO appAdj;
    mutable ASGapplyTangentFO appTan;
    mutable ASGapplyAdjTangentFO appAdjTan;

    TSFO & applyFO() const { return app; }
    TSFO & applyAdjFO() const { return appAdj; }
    TSFO & applyTangentFO() const { return appTan; }
    TSFO & applyAdjTangentFO() const { return appAdjTan; }

  protected:

    void stepTimeFwd() const {
      this->setTime(this->getTime()+0.5*dt);
      integerstep = !integerstep;
    } 
    void stepTimeBwd() const {
      this->setTime(this->getTime()-0.5*dt);
      integerstep = !integerstep;
    } 
    
    Operator<float> * clone() const { return new ASGStep(*this); }
    
  public:

    ASGStep(ASGStep const & stp)
      : gpad(stp.gpad),     
	gdom(stp.gdom),
	tmin(stp.tmin), dt(stp.dt),
	nls(stp.nls), nrs(stp.nrs),
	aux(stp.aux),
	app(aux,integerstep),
	appAdj(aux,integerstep),
	appTan(aux,integerstep),
	appAdjTan(aux,integerstep) {

      try {
	this->setTime(stp.tmin);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ASGStep constructor\n";
	throw e;
      }
    }

    ASGStep(Grid const & phys, int maxoff,
	    float _tmin, float _dt,
	    std::vector<int> _nls,
	    std::vector<int> _nrs,		    
	    float amp)
      : gpad(make_padded_grid(phys,_nls,_nrs)),
	gdom(make_asg_ctrllist(*gpad,maxoff),
	     make_asg_statelist(*gpad,maxoff),
	     make_asg_indices(gpad->getDimension())),
	tmin(_tmin), dt(_dt),
	nls(_nls), nrs(_nrs),
	aux(gdom,maxoff,dt,_nls,_nrs,amp),
	integerstep(true),
	app(aux,integerstep),
	appAdj(aux,integerstep),
	appTan(aux,integerstep),
	appAdjTan(aux,integerstep) {

      try {
	this->setTime(tmin);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ASGStep constructor\n";
	throw e;
      }
    }

    float getTimeStep() const {return dt; }
    Space<float> const & getDomain() const { return gdom; }
    Space<float> const & getRange() const { return gdom; }

    ostream & write(ostream & str) const {
      str<<"ASGStep\n";
      str<<"  tmin="<<tmin<<" dt="<<dt<<"\n";
      str<<"  padding:\n";
      for (int i=0;i<gpad->getDimension(); i++)
	str<<"  axis="<<i<<" left="<<nls[i]<<" right="<<nrs[i]<<"\n";
      str<<"  padded physical grid:\n";
      gpad->write(str);
      str<<"  GridDomain:\n";
      gdom.write(str);
      return str;
    }
  };

}
