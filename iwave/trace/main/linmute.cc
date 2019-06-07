#include "parser.h"
#include "segypp.hh"
#include "segyops.hh"
#include "op.hh"

using RVL::valparse;
using RVL::RVLException;
using RVL::Vector;
using RVL::Components;
using RVL::LinearOpFO;
using RVL::AssignFilename;

using TSOpt::SEGYLinMute;
using TSOpt::SEGYSpace;
typedef TSOpt::SEGYSpace tsp;

const char * sdoc[] = {
  " IWAVE linear mute command",
  " arguments:",
  "   input      = filename of input SEGY data",
  "  output      = filename of output SEGY data",
  " mute_slope   = dtime/d|offset| in ms/m - typical value = 0.67 corresponds to H2O Vp",
  " mute_zotime  = mute onset time at zero offset",
  " mute_width   = width of linear taper zone",
  NULL};

int xargc;
char **xargv;

int main(int argc, char ** argv) {

  try {

    if (retrieveGlobalRank()==0 && argc<2) {
      pagedoc();
      exit(0);
    }

    PARARRAY * pars = ps_new();
    
    if (ps_createargs(pars,argc-1,&(argv[1]))) {
      RVLException e;
      e<<"ERROR: mute from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++) 
	e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }
    // since the product of grid spaces is not really an 
    // out-of-core structure, this driver operates on single
    // grid spaces
    string inp = valparse<string>(*pars,"input");
    string outp = valparse<string>(*pars,"output");

    tsp dom(inp,"data");
        
    Vector<ireal> ddin(dom);
    Vector<ireal> ddout(dom);
        
    AssignFilename ddinfn(inp);
    Components<ireal> cddin(ddin);
    cddin[0].eval(ddinfn);
        
    AssignFilename ddoutfn(outp);
    Components<ireal> cddout(ddout);
    cddout[0].eval(ddoutfn);

    SEGYLinMute mute(valparse<float>(*pars,"mute_slope",0.0f),
		     valparse<float>(*pars,"mute_zotime",0.0f),
		     valparse<float>(*pars,"mute_width",0.0f));

    LinearOpFO<float> muteop(dom,dom,mute,mute);
    muteop.applyOp(ddin,ddout);
  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
