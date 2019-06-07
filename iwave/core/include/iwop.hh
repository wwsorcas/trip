#ifndef __IWAVE_OP
#define __IWAVE_OP

#define DEFAULT_SNAPS 10

//#include "alg.hh"
#include "op.hh"
#include "productspace.hh"
#include "blockop.hh"
#include "ocdc.hh"
#include "gridpp.hh"
#include "segypp.hh"
#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#include "mpisegypp.hh"
#endif
#include "gridops.hh"
#include "logistic.hh"
#include "iwsim.hh"

namespace TSOpt {

  //using namespace RVLAlg;
  using RVLAlg::ListAlg;
  using RVL::DataContainer;
  using RVL::ProductDataContainer;
  using RVL::StdProductDataContainer;
  using RVL::Space;
  using RVL::SpaceDCF;
  using RVL::ProductSpace;
  using RVL::StdProductSpace;
  using RVL::ConstContainer;
  using RVL::STRING_PAIR;
  using RVL::Vector;
  using RVL::Components;
  using RVL::FunctionObject;
  using RVL::Operator;
  using RVL::LinOpValOp;
  using RVL::Writeable;
  using RVL::AssignParams;
  
  class IWaveSpace: public ProductSpace<float> {

  private:

    /** vector of const pointers to Space */
    std::vector< Space<float> const * > _s;
    std::vector< std::string > _keys;

  protected:

    Space<float> * clone() const {
      return new IWaveSpace(*this);
    }
    
  public:

    IWaveSpace(PARARRAY const & par, 
	       IWaveInfo const & ic,
	       bool input,
	       int lflag = 1,
	       ostream & outfile = cerr);

    IWaveSpace(IWaveSpace const & sp);
    
    ~IWaveSpace();

    /** implements virtual DataContainer constructor via
	StdProductDataContainer class. */
    DataContainer * buildDataContainer() const;

    size_t getSize() const;
    Space<float> const & operator[](size_t i) const;
    std::vector<std::string> getKeys() const;
  };

  // convenience filename transfer - weak sanity check, presume 
  // that membership is already established
  // change 25.03.16: hoisted out of IWaveOp class,
  // membership test added
  void param_set(RVL::Vector<float> const & x, 
		 PARARRAY & pars, 
		 IWaveSpace const & sp,
		 std::string const & suf,
		 FILE * stream=stderr);

  class IWaveOp: public Operator<float>  {
      
  private:

    IWaveInfo ic;
    IWaveSpace dom;
    IWaveSpace rng;
    mutable FILE * stream;              /* output stream            */
    PARARRAY * pars;            /* parameter array ref copy */

    // verbosity control
    int dump_steps;
    int dump_pars;
    int dump_term;

    // other verbosity control handled within iwave code
    
    // dry run option
    bool dryrun;
    ostream & drystr;
    
    // verbose announcements
    ostream & announce;

    IWaveOp();
      
  protected:
      
    void apply(const Vector<float> & x, 
	       Vector<float> & y) const;

    void applyDeriv(const Vector<float> & x, 
		    const Vector<float> & dx,
		    Vector<float> & dy) const;
      
    void applyAdjDeriv(const Vector<float> & x, 
		       const Vector<float> & dy,
		       Vector<float> & dx) const;
      
    void applyDeriv2(const Vector<float> & x, 
		     const Vector<float> & dx0,
		     const Vector<float> & dx1,
		     Vector<float> & dy) const;

    void applyAdjDeriv2(const Vector<float> & x, 
			const Vector<float> & dx0,
			const Vector<float> & dy,
			Vector<float> & dx1) const;

    Operator<float> * clone() const;
    
  public:
    
    IWaveOp(PARARRAY _pars, FILE * _stream,
	    bool _dryrun=false,
	    ostream & _drystr=cerr,
	    ostream & _announce=cerr);
    
    IWaveOp(IWaveOp const & x);

    ~IWaveOp();

    const IWaveSpace & getIWaveDomain() const { return dom; } 
    const IWaveSpace & getIWaveRange() const { return rng; } 
    const Space<float> & getDomain() const { return getIWaveDomain(); }
    const Space<float> & getRange() const { return getIWaveRange(); }

    // added 23.06.10 to facilitate using source as variable
    // without admitting that it's part of domain
    PARARRAY & getPar();
    PARARRAY const & getPar() const;

    std::string getModel() const { return ic.get_iwave_model(); }

    ostream & write(ostream & str) const;
  };

  class IWaveLOVOp: public LinOpValOp<float>  {
      
  private:

    IWaveInfo ic;
    StdProductSpace<float> dom;
    IWaveSpace rng;
    mutable FILE * stream;      /* output stream            */
    PARARRAY * pars;            /* parameter array ref copy */

    // verbosity control
    int dump_steps;
    int dump_pars;
    int dump_term;

    // other verbosity control handled within iwave code
    
    // dry run option
    bool dryrun;
    ostream & drystr;
    
    // verbose announcements
    ostream & announce;

    IWaveLOVOp();
      
  protected:
      
    virtual void apply0(const Vector<float> & x0,
			const Vector<float> & x1, 
			Vector<float> & y) const;
    
    void applyAdj0(const Vector<float> & x0,
		   const Vector<float> & y, 
		   Vector<float> & x1) const;
    
    void applyPartialDeriv0(const Vector<float> & x0,
			    const Vector<float> & x1,
			    const Vector<float> & dx0,
			    Vector<float> & dy) const;
    
    void applyAdjPartialDeriv0(const Vector<float> & x0,
			       const Vector<float> & x1,
			       const Vector<float> & dy,
			       Vector<float> & dx0) const;
    
    void applyPartialDeriv20(const Vector<float> & x0,
			     const Vector<float> & x1,
			     const Vector<float> & dx00,
			     const Vector<float> & dx01,
			     Vector<float> & dy) const;
    
    void applyAdjPartialDeriv20(const Vector<float> & x0,
				const Vector<float> & x1,
				const Vector<float> & dx00,
				const Vector<float> & dy,
				Vector<float> & dx01) const;

    RVL::OperatorProductDomain<float> * clonePD() const {
      return new IWaveLOVOp(*this);
    }
    
  public:
    
    IWaveLOVOp(PARARRAY _pars, FILE * _stream,
	       bool _dryrun=false,
	       ostream & _drystr=cerr,
	       ostream & _announce=cerr);
    
    IWaveLOVOp(IWaveLOVOp const & x);

    ~IWaveLOVOp();

    // convenience functions
    const IWaveSpace & getNonLinDomain() const {
      try {
	IWaveSpace const & iwspace =
	  dynamic_cast<IWaveSpace const &>(this->getProductDomain()[0]);
	return iwspace;
      }
      catch (RVLException & e) {
	e<<"\ncalled from IWaveLOVOp::getNonLinDomain\n";
	throw e;
      }
      catch (bad_cast) {
	RVLException e;
	e<<"ERROR: IWaveLOVOp::getNonLinDomain\n";
	e<<"  who knows how - first component of domain not iwspace\n";
	throw e;
      }
    }
    const IWaveSpace & getLinDomain() const {
      try {
	IWaveSpace const & iwspace =
	  dynamic_cast<IWaveSpace const &>(this->getProductDomain()[1]);
	return iwspace;
      }
      catch (RVLException & e) {
	e<<"\ncalled from IWaveLOVOp::getLinDomain\n";
	throw e;
      }
      catch (bad_cast) {
	RVLException e;
	e<<"ERROR: IWaveLOVOp::getLinDomain\n";
	e<<"  who knows how - first component of domain not iwspace\n";
	throw e;
      }
    }
    const IWaveSpace & getIWaveRange() const { return rng; }

    // mandatory
    const ProductSpace<float> & getProductDomain() const { return dom; }
    const Space<float> & getRange() const { return rng; }

    // added 23.06.10 to facilitate using source as variable
    // without admitting that it's part of domain
    PARARRAY & getPar();
    PARARRAY const & getPar() const;

    std::string getModel() const { return ic.get_iwave_model(); }

    ostream & write(ostream & str) const;
  };

}

#endif
