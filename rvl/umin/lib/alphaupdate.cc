#include "alphaupdate.hh"

namespace RVLUmin {

  float alphaupdate(RVL::Vector<float> & d,
		    RVL::Vector<float> & e,
		    RVL::Vector<float> & p,
		    float alphacurr,
		    float ubnd,
		    float lbnd,
		    float decr) {
    try {
      if ((ubnd > 1.0) || (lbnd < 0.0) ||
	  (ubnd <= lbnd)) {
	RVL::RVLException e;
	e<<"Error: RVLUmin::alphaupdate\n";
	e<<"  lbnd, ubnd=["<<lbnd<<", "<<ubnd<<"] do no obey necessary \n";
	e<<"  relation 0 < lbnd < ubnd < 1\n";
	throw e;
      }	
      if ((decr >= 1.0f) ||
	  (decr <= 0.0f)) {
	RVL::RVLException e;
	e<<"Error: RVLUmin::alphaupdate\n";
	e<<"  decrease factor "<<decr<<" outside of range (0,1)\n";
	throw e;
      }
      if (alphacurr < 0.0f) {
	RVL::RVLException e;
	e<<"Error: RVLUmin::alphaupdate\n";
	e<<"  current penalty weight "<<alphacurr<<" negative\n";
	throw e;
      }
      //      if (alphacurr < numeric_limits<float>::epsilon()) {
      float dn = d.normsq();
      float en = e.normsq();
      float pn = p.normsq();
      float rn=0.0f;
      float sn=0.0f;
      if (RVL::ProtectedDivision<float>(en,dn,rn) ||
	  RVL::ProtectedDivision<float>(pn,dn,sn)) {
	RVL::RVLException e;
	e<<"Error: RVLUmin::alphaupdate\n";
	e<<"  input reference vector effectively zero! cannot proceed\n";
	throw e;
      }
      if (ubnd < rn) {
	if (alphacurr < numeric_limits<float>::epsilon()) {
	  /*
	  RVL::RVLException e;
	  e<<"Error: RVLUmin::alphaupdate\n";
	  e<<"  current penalty weight "<<alphacurr<<" too small to decrease\n";
	  e<<"  either residual (1st arg) not reduced sufficiently or\n";
	  e<<"  data cannot be fit arbitrarily accurately\n";
	  throw e;
	  */
	  cerr<<"Error: RVLUmin::alphaupdate\n";
	  cerr<<"  current penalty weight "<<alphacurr<<" too small to decrease further\n";
	  cerr<<"  either residual (1st arg) not reduced sufficiently or\n";
	  cerr<<"  data cannot be fit arbitrarily accurately\n";
	  
	  return alphacurr;
	}
	return decr*alphacurr;
      }
      // on first step, with alpha=0, always update
      else if ((lbnd > rn) ||
	       (alphacurr*sn < numeric_limits<float>::epsilon())) {
	return (alphacurr + 0.5*(ubnd*dn - en)/pn);
      }
      else {
	return alphacurr;
      }
    }
    catch (RVL::RVLException & e) {
      e<<"\ncalled from RVLUmin::alphaupdate\n";
      throw e;
    }
  }
}

