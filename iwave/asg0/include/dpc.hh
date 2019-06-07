#ifndef __ASG_DPC
#define __ASG_DPC

#include "iwop.hh"

#ifdef IWAVE_USE_MPI
#include "mpisegypp.hh"
#else
#include "segypp.hh"
#endif

namespace TSOpt {

#ifdef IWAVE_USE_MPI
  typedef MPISEGYSpace mySEGYSpace;
#else
  typedef SEGYSpace mySEGYSpace;
#endif

  class Parametrix: public RVL::LinOpValOp<float> {
  
  private:

    IWaveLOVOp iwop0;
    PARARRAY * modpars;
    IWaveLOVOp iwop;
    RVL::Vector<float> m;
    RVL::RestrictOp<float> rop;
    RVL::TangentMap<float> top;
    
    Parametrix();

  protected:

    void apply0(const RVL::Vector<float> & x0,
		const RVL::Vector<float> & x1, 
		RVL::Vector<float> & y) const;
    
    void applyAdj0(const RVL::Vector<float> & x0,
		   const RVL::Vector<float> & y, 
		   RVL::Vector<float> & x1) const;
    
    void applyPartialDeriv0(const RVL::Vector<float> & x0,
			    const RVL::Vector<float> & x1,
			    const RVL::Vector<float> & dx0,
			    RVL::Vector<float> & dy) const;
    
    void applyAdjPartialDeriv0(const RVL::Vector<float> & x0,
			       const RVL::Vector<float> & x1,
			       const RVL::Vector<float> & dy,
			       RVL::Vector<float> & dx0) const;
    
    void applyPartialDeriv20(const RVL::Vector<float> & x0,
			     const RVL::Vector<float> & x1,
			     const RVL::Vector<float> & dx00,
			     const RVL::Vector<float> & dx01,
			     RVL::Vector<float> & dy) const {
      RVL::RVLException e;
      e<<"Error: Parametrix::applyPartialDeriv20 not implemented\n";
      throw e;
    }
    
    void applyAdjPartialDeriv20(const RVL::Vector<float> & x0,
				const RVL::Vector<float> & x1,
				const RVL::Vector<float> & dx00,
				const RVL::Vector<float> & dy,
				RVL::Vector<float> & dx01) const {
      RVL::RVLException e;
      e<<"Error: Parametrix::applyAdjPartialDeriv20 not implemented\n";
      throw e;
    }

    RVL::OperatorProductDomain<float> * clonePD() const {
      return new Parametrix(*this);
    }
    
  public:

    Parametrix(PARARRAY pars, FILE * stream);
    Parametrix(const Parametrix &);
    ~Parametrix() { ps_delete(&modpars); }

       // convenience functions
    const IWaveSpace & getNonLinDomain() const {
      try {
	// nonlin domain of iwop = (bulkmod, buoyancy)
	return iwop.getNonLinDomain();
      }
      catch (RVLException & e) {
	e<<"\ncalled from Parametrix::getLinDomain\n";
	throw e;
      }
    }
    const IWaveSpace & getLinDomain() const {
      try {
	return iwop.getNonLinDomain();
      }
      catch (RVLException & e) {
	e<<"\ncalled from Parametrix::getLinDomain\n";
	throw e;
      }
    }
    const IWaveSpace & getIWaveRange() const {
      try {
	return iwop0.getIWaveRange();
      }
      catch (RVLException & e) {
	e<<"\ncalled from Parametrix::getLinDomain\n";
	throw e;
      }
    }
    
    // note - domain is nonlinear subspace of *** modified *** op (domain padded,
    // range to be abs surface data
    const RVL::ProductSpace<float> & getProductDomain() const {
      return iwop.getNonLinDomain(); }
    // whereas range is free surface data from unmodified op
    const RVL::Space<float> & getRange() const {
      return iwop0.getRange();
    }

    PARARRAY & getPar() { return iwop.getPar(); }
    PARARRAY const & getPar() const { return iwop.getPar(); }

    ostream & write(ostream & str) const;
  };

}

#endif
