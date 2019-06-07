#ifndef __DTON_IWLOVOP__
#define __DTON_IWLOVOP__

#include "iwop.hh"
#include "segyops.hh"
#include "gridops.hh"

#ifdef IWAVE_USE_MPI
typedef TSOpt::MPISEGYSpace mySEGYSpace;
typedef TSOpt::MPIGridSpace myGridSpace;
#else
typedef TSOpt::SEGYSpace mySEGYSpace;
typedef TSOpt::GridSpace myGridSpace;  
#endif

namespace TSOpt {

  /**
     Class designed for linear acoustics to update bulk modulus (that is, velocity) by surface source extension principle. Active variables: bulkmod (nonlinear input), pressure source (linear input), pressure data (linear output). Uses internal IWaveLOVOp that adds velocity data (linear output) to implement D-to-N map preconditioning. 
   */
  
  class DToNIWaveLOVOp: public RVL::LinOpValOp<float> {
    
  private:

    IWaveInfo ic;
    RVL::StdProductSpace<float> dom;
    std::shared_ptr<TSOpt::IWaveLOVOp> iwop;
    int appinv;
    
  protected:

    void apply0(const Vector<float> & x0,
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
    
    RVL::OperatorProductDomain<float> * clonePD() const {
      return new DToNIWaveLOVOp(*this);
    }
    
  public:
    
    DToNIWaveLOVOp(PARARRAY pars,
		   FILE * stream);

    DToNIWaveLOVOp(DToNIWaveLOVOp const & op);


    ProductSpace<float> const & getProductDomain() const {
      return dom;
    }
    
    Space<float> const & getRange() const {
      return iwop->getIWaveRange()[0];
    }

    std::shared_ptr<TSOpt::IWaveLOVOp> const getIWOP() const {
      return iwop;
    }

    ostream & write(ostream & str) const;

  };
}

#endif
