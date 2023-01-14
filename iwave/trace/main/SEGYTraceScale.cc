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

using TSOpt::TraceScaleFO;
using TSOpt::SEGYSpace;
typedef TSOpt::SEGYSpace tsp;

const char * sdoc[] = {
  " Trace scaling by uniform factor and optional header word",
  " arguments:",
  "   input      = filename of input SEGY data (string)",
  "  output      = filename of output SEGY data (string)",
  "     fac      = uniform factor (float)",
  "     key      = header key (string)",
  " NOTE: applies scaling ONLY to first factor of product data set",
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

    tsp dom(inp,"data");
        
    Vector<ireal> ddin(dom);
    Vector<ireal> ddout(dom);
        
    AssignFilename ddinfn(inp);
    Components<ireal> cddin(ddin);
    cddin[0].eval(ddinfn);
        
    AssignFilename ddoutfn(outp);
    Components<ireal> cddout(ddout);
    cddout[0].eval(ddoutfn);

    float fac = valparse<float>(*pars,"fac");
    std::string key = valparse<std::string>(*pars,"key","");
    TraceScaleFO fo(fac,key);
    LinearOpFO<float> tsop(dom,dom,fo,fo);
    tsop.applyOp(ddin,ddout);

    ps_delete(&pars);
    iwave_fdestroy();
    exit(0);
  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
