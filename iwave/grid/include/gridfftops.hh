#ifndef __TSOPT_GRIDFFT_OPS__
#define __TSOPT_GRIDFFT_OPS__

#include "op.hh"
#include "productspace.hh"
#include "mpiserialfo.hh"
#include "fftops.hh"

#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
typedef TSOpt::MPIGridSpace myGridSpace;
#else
#include "gridpp.hh"
typedef TSOpt::GridSpace myGridSpace;
#endif

namespace TSOpt {

  class GridZFTScaleFO: public RVL::UnaryLocalFunctionObject<float> {
  private:
    float dz;
    float pow;
    int ab;
    float * band;
    int adj;
  public:
    GridZFTScaleFO(float _dz,
		   float _pow,
		   int _ab,
		   float * _band,
		   int _adj) 
      : dz(_dz), pow(_pow), ab(_ab),
	band(new float[4]), adj(_adj) {
      band[0]=_band[0];
      band[1]=_band[1];
      band[2]=_band[2];
      band[3]=_band[3];
    }
    GridZFTScaleFO(GridZFTScaleFO const & fo)
      : dz(fo.dz), pow(fo.pow), ab(fo.ab),
	band(new float[4]), adj(fo.adj) {
      band[0]=fo.band[0];
      band[1]=fo.band[1];
      band[2]=fo.band[2];
      band[3]=fo.band[3];
    }      
    ~GridZFTScaleFO() { delete [] band; }
    using RVL::LocalEvaluation<float>::operator();
    void operator()(RVL::LocalDataContainer<float> &);
    string getName() const { string tmp = "GridZFTScaleFO"; return tmp; }
  };
    
  class GridZFTScaleOp: public RVL::LinearOp<float> {
  
  private:
    
    RVL::Space<float> const & dom;
    float pow;
    int ab;
    float * band;

    float dz;
    
    GridZFTScaleOp();

  protected:

    void apply(RVL::Vector<float> const & x,
	       RVL::Vector<float> & y) const;
    void applyAdj(RVL::Vector<float> const & x,
		  RVL::Vector<float> & y) const;

    RVL::Operator<float> * clone() const { return new GridZFTScaleOp(*this); }
    
  public:
    
    /** main constructor: takes Grid defining increment window, and width
	parameter defining zone of smooth taper - same on all sides. 
    */
    GridZFTScaleOp(RVL::Space<float> const & _dom,
		 float _pow,
		 int ab,
		 float _locut,
		 float _lopas,
		 float _hipas,
		 float _hicut);

    /** Copy constructor - memberwise */
    GridZFTScaleOp(GridZFTScaleOp const & op);

    ~GridZFTScaleOp() {}
    
    RVL::Space<float> const & getDomain() const { return dom; }
    RVL::Space<float> const & getRange() const { return dom; }

    ostream & write(ostream & str) const;
  };	
  
}
#endif
