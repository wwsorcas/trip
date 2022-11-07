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
using TSOpt::SSEScaleFO;
using TSOpt::SEGYSpace;
typedef TSOpt::SEGYSpace tsp;

const char * sdoc[] = {
  " Trace scaling by (shift+alpha^2*offset^2)^p",
  " arguments:",
  "   input      = filename of input SEGY data (string)",
  "  output      = filename of output SEGY data (string)",
  "   shift      = shift (float)",
  "   alpha      = weight (float)",
  "       p      = exponent (float)",
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

    string inp = valparse<string>(*pars,"input");
    string outp = valparse<string>(*pars,"output");
    float shift = valparse<float>(*pars,"shift");
    float alpha = valparse<float>(*pars,"alpha");
    float p = valparse<float>(*pars,"p");

    tsp dom(inp,"data");
        
    Vector<ireal> ddin(dom);
    Vector<ireal> ddout(dom);
        
    AssignFilename ddinfn(inp);
    Components<ireal> cddin(ddin);
    cddin[0].eval(ddinfn);
        
    AssignFilename ddoutfn(outp);
    Components<ireal> cddout(ddout);
    cddout[0].eval(ddoutfn);

    SSEScaleFO fo(shift,alpha,p);
    LinearOpFO<float> tsop(dom,dom,fo,fo);
    tsop.applyOp(ddin,ddout);
  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
