#ifndef __RVL_UMIN_ALPHAUPDATE__
#define __RVL_UMIN_ALPHAUPDATE__

#include "space.hh"

namespace RVLUmin {

  float alphaupdate(RVL::Vector<float> & d,
		    RVL::Vector<float> & e,
		    RVL::Vector<float> & p,
		    float alphacurr,
		    float ubnd,
		    float lbnd,
		    float decr);

}

#endif
