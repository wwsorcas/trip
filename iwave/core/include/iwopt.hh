#ifndef __IWAVE_OPT
#define __IWAVE_OPT

#include "functional.hh"
#include "ls.hh"
#include "cgnealg.hh"
#include "TRGNAlg.hh"
#include "LBFGSBT.hh"
#include "iwop.hh"

namespace TSOpt {

#ifdef IWAVE_USE_MPI
  typedef MPIGridSpace myGridSpace;
#else
  typedef GridSpace myGridSpace;
#endif
  
  void IWaveLoad(PARARRAY,
		 Vector<float> &,
		 std::vector<std::string>);

  /**
     parameters for iwavefwi:

     all parameters for underlying modeling operator

     for each component to be updated: name should take form 
     [name from keys]_est, with [name from keys] being the keyword
     associated with that component

     grid spec'd in file [name from keys]_est is update grid, to be
     injected into grid spec'd by [name_from_keys] = archival inverted
     data (initial on call) - so former must be subgrid of latter

     upper and lower bounds are either
     [name_from_keys]_ub, [name_from_keys]_lb (files for envelope, 
     must have same grid as _est file, or
     [name_from_keys]_max, [name_from_keys]_min (floats)

     tapers on grid windows (same for all): taper1, taper2, taper3

   */
  
  class IWaveFWIOp: public Operator<float> {

  private:

    IWaveOp iwop;
    std::shared_ptr<Vector<float> > m0;
    std::shared_ptr<Vector<float> > m;
    std::vector<std::string> est;
    std::vector<std::string> estkeys;
    std::vector<int> estidx;
    std::shared_ptr<Space<float> > dom;
    std::shared_ptr<Space<float> > med;
    std::shared_ptr<Vector<float> > ub;
    std::shared_ptr<Vector<float> > lb;
    std::shared_ptr<Vector<float> > lm0;
    std::shared_ptr<Operator<float> > dwop;
    std::shared_ptr<Operator<float> > dlop;
    std::shared_ptr<Operator<float> > injop;
    std::shared_ptr<Operator<float> > transop;
    std::shared_ptr<Operator<float> > modelop;
    std::shared_ptr<Operator<float> > op;

    IWaveFWIOp();

  protected:

    void apply(Vector<float> const & x,
	       Vector<float> & y) const;
    
    void applyDeriv(Vector<float> const & x,
		    Vector<float> const & dx,
		    Vector<float> & dy) const;
    
    void applyAdjDeriv(Vector<float> const & x,
		       Vector<float> const & dy,
		       Vector<float> & dx) const;
    
    void applyDeriv2(const Vector<float> & x,
		     const Vector<float> & dx1,
		     const Vector<float> & dx2,
		     Vector<float> & dy) const;
    
    void applyAdjDeriv2(const Vector<float> & x,
			const Vector<float> & dy,
			const Vector<float> & dx2,
			Vector<float> & dx1) const;

    Operator<float> * clone() const { return new IWaveFWIOp(*this); }

  public:

    IWaveFWIOp(PARARRAY par, FILE * stream);

    IWaveFWIOp(IWaveFWIOp const & x);

    ~IWaveFWIOp() {}

    Space<float> const & getDomain() const { return op->getDomain(); }
    Space<float> const & getRange() const { return op->getRange(); }

    Operator<float> const & getTransformOp() const { return *transop; }
    Operator<float> const & getModelingOp() const { return *modelop; }

    std::vector<std::string> getDomainKeys() const { return estkeys; }
    std::vector<std::string> getModelKeys() const {
      IWaveSpace const & iwdom =
	dynamic_cast<IWaveSpace const &>(iwop.getDomain());
      return iwdom.getKeys();
    }
    std::vector<std::string> getRangeKeys() const {
      IWaveSpace const & iwrng =
	dynamic_cast<IWaveSpace const &>(iwop.getRange());
      return iwrng.getKeys();
    }
    
    ostream & write(ostream & str) const;

  };
  
  void IWaveOpt(PARARRAY pars,
		Operator<float> const & op,
		std::shared_ptr<RVL::LinearOp<float> const> prec,		
		Vector<float> const & d,
		Vector<float> & m);

  void IWaveResMod(PARARRAY pars,
		   Operator<float> const & op,
		   Vector<float> const & d,
		   Vector<float> const & m);

  int IWaveOpApply(int, char **,
		   void (*invfcn)(PARARRAY,
				  RVL::Builder<RVL::LinearOp<float>, GridOpData> &,
				  FILE *),
		   RVL::Builder<RVL::LinearOp<float>, GridOpData> &);

  void StraightOLS(PARARRAY, RVL::Builder<RVL::LinearOp<float>, GridOpData> &, FILE*);
  
}

#endif
