#ifndef __IWAVE_ANG_FFTOPS__
#define __IWAVE_ANG_FFTOPS__

#include <cmath>
#include <complex>
#include <memory.h>
#include "local.hh"
#include "utils.h"

using namespace std;

typedef complex<float> float_complex;

namespace TSOpt {

  void FTScale(float *Trace, int Nt, float dt,
	       int ab, int power, float *band);

  /** 1D scale by frequency or absolute frequency power, followed
      by bandpass filter  - in-place version */
  class FTScaleFO: public RVL::UnaryLocalFunctionObject<float> {
  private:
    float dt;
    float power;
    int ab;
    float * band;
  public:
    FTScaleFO(float _dt = 1.0f,
	      float _power = 0.0f,
	      int _ab = 0,
	      float locut = 0.0f,
	      float lopas = 0.0f,
	      float hipas = std::numeric_limits<float>::max(),
	      float hicut = std::numeric_limits<float>::max())
      : dt(_dt), power(_power), ab(_ab),
	band(new float[4]) {
      band[0]=locut;
      band[1]=lopas;
      band[2]=hipas;
      band[3]=hicut;
    }
    FTScaleFO(FTScaleFO const & fo)
      : dt(fo.dt), power(fo.power), ab(fo.ab),
	band(new float[4]) {
      band[0]=fo.band[0];
      band[1]=fo.band[1];
      band[2]=fo.band[2];
      band[3]=fo.band[3];
    }
    ~FTScaleFO() { delete [] band; }
    using RVL::LocalEvaluation<float>::operator();
    void operator()(RVL::LocalDataContainer<float> &);
    string getName() const { string tmp = "FTScaleFO"; return tmp; }
  };    
    
}

#endif
