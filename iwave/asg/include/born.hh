#ifndef __BORN_IWLOVOP__
#define __BORN_IWLOVOP__

//#define __INTEGRAL__
#undef __INTEGRAL__

#include "iwop.hh"
#include "segyops.hh"
#include "gridops.hh"
#include "gridfftops.hh"
#include "appinv_utils.hh"

#ifdef IWAVE_USE_MPI
typedef TSOpt::MPISEGYSpace mySEGYSpace;
typedef TSOpt::MPIGridSpace myGridSpace;
#else
typedef TSOpt::SEGYSpace mySEGYSpace;
typedef TSOpt::GridSpace myGridSpace;  
#endif


namespace TSOpt {
  
  class BornIWaveLOVOp: public RVL::LinOpValOp<float> {
    
  private:

    IWaveInfo ic;
    std::vector<std::string> bkeys;
    std::shared_ptr<IWaveLOVOp> iwopfree;
    std::shared_ptr<IWaveLOVOp> iwopabsb;
    RVL::StdProductSpace<float> freedom;
    RVL::StdProductSpace<float> absbdom;
    std::vector<size_t> idx;
    int freechoice;
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
      return new BornIWaveLOVOp(*this);
    }
    
  public:
    
    BornIWaveLOVOp(std::vector<std::string> _bkeys,
		   PARARRAY pars,
		   FILE * stream);

    BornIWaveLOVOp(BornIWaveLOVOp const & bwop);

    std::vector<std::string> getRangeKeys() const {
      // keys are some in both cases
      return iwopfree->getIWaveRange().getKeys();
    }

    ProductSpace<float> const & getProductDomain() const {
      if (freechoice) return freedom;
      else return absbdom;
    }
    Space<float> const & getRange() const {
      if (freechoice) return iwopfree->getRange();
      else return iwopabsb->getRange();
    }

    IWaveLOVOp const & getIWOP() const {
      if (freechoice) return *iwopfree;
      else return *iwopabsb;
    }

    std::vector<size_t> getIDX() const { return idx; }
    
    ostream & write(ostream & str) const;

  };
}

#endif
