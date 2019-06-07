#ifndef __ASG_APPINV_UTILS__
#define __ASG_APPINV_UTILS__

#include "linop.hh"
#include "mpiserialfo.hh"
#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#include "mpisegypp.hh"
typedef TSOpt::MPIGridSpace myGridSpace;
typedef TSOpt::MPISEGYSpace mySEGYSpace;
#else
#include "gridpp.hh"
#include "segypp.hh"
typedef TSOpt::GridSpace myGridSpace;
typedef TSOpt::SEGYSpace mySEGYSpace;
#endif
#include "gridops.hh"
#include "gridfftops.hh"

namespace TSOpt {

  class ASGModelWeightOp: public RVL::LinearOpWithInverse<float> {

  private:
    RVL::Vector<float> const & bulk;
    RVL::Vector<float> const & buoy;
    RVL::Vector<float> const & extbulk;
    float * band;
    float dz;
    int symm;

  protected:
    void apply(RVL::Vector<float> const & x,
	       RVL::Vector<float> & y) const;

    void applyAdj(RVL::Vector<float> const & x,
		  RVL::Vector<float> & y) const;

    void applyInv(RVL::Vector<float> const & x,
		  RVL::Vector<float> & y) const;

    void applyInvAdj(RVL::Vector<float> const & x,
		     RVL::Vector<float> & y) const;

    RVL::Operator<float> * clone() const {
      return new ASGModelWeightOp(*this);
    }
      
  public:
    ASGModelWeightOp(RVL::Vector<float> const & bulk,
		     RVL::Vector<float> const & buoy,
		     RVL::Vector<float> const & extbulk,
		     float locut, float lopas, float hipas, float hicut,
		     int symm);

    ASGModelWeightOp(ASGModelWeightOp const & op);

    ~ASGModelWeightOp() { delete [] band; }

    RVL::Space<float> const & getDomain() const { return extbulk.getSpace(); }
    RVL::Space<float> const & getRange() const { return extbulk.getSpace(); }

    ostream & write(ostream & str) const;

  };

}

#endif
