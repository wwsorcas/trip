#ifndef __IW_SAMP__
#define __IW_SAMP__

#include "except.hh"
#include "parser.h"
#include "parserdecl.hh"
#include "iwave.h"
#include "grid.h"
#include "traceio.h"
#include "findsuffix.hh"

namespace TSOpt {
  using RVL::parse;
  using RVL::valparse;
  using RVL::RVLException;
  using RVL::ProtectedDivision;

  class IWaveSampler {
  private:
    string samplekey;
    string suffix;
    string pname;
    // FILE * stream; 
    std::vector<axis *> axes; 
    axis const * getTimeAxis(int dim) const;
    bool has_Axis(int i) const;
    bool has_Spatial_Axes;
    mutable int prev_panelindex;
    bool increment;

    // raw info for refining loop limits
    IPNT gmin;
    IPNT gmax;
    
    // this object really should be private data
    // of a trace sampler subclass or independent object
    tracegeom * tg;
    int sampord;
    int tracestart;
    int tracestop;
    int taperwidth;   // width of taper for each side for one shot in (# of traces)
    int timewidth;    // width of time to taper at the end of simulation in (ms)

    float muteslope;  // slope of mute (dt/dx) (ms/m)
    float mutezotime; // mute onset at zero offset, time AFTER first sample (ms)
    float mutewidth;  // width of mute zone (ms)

    int dump_term;

  public:
    /** main constructor: associates sampler with grid structure 
	of file = value for key in pars. Does not define sampling - 
	that would require IOTASK defn, happens on call to sample
    */
    IWaveSampler(IWAVE * state, string key, PARARRAY & pars, FILE * stream);
    virtual ~IWaveSampler();
    int getNumAxes() const { return axes.size(); }
    axis const & getAxis(int i) const;
    ireal getCellVol() const;
    ireal getRecipCellVol() const;
    int const * get_gmin() const { return gmin; }
    int const * get_gmax() const { return gmax; }
    /** main method - does far too much work in this version 
	@param g = simulation grid
	@param step = simulation step vector step[dim] is time step index
	@param gtype = stagger index array - should be data of rarray
	@param fwd = fwd/adj flag
	@param input = set if input to simulation, else output
	@param state = collection of rdoms defining simulation
	@param ridx = rarray index within rdom
	@param iwdx = rdom index within state
	@param stream = verbose output
     */
    void sample(grid g, IPNT step, IPNT gtype, bool fwd, bool input, 
		IWAVE * state, int ridx, int iwdx, FILE * stream,
		bool dryrun=false, ostream & drystr=cerr);
  };
	  
}

#endif
