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

using TSOpt::SEGYTaper;
using TSOpt::SEGYSpace;
typedef TSOpt::SEGYSpace tsp;

int xargc;
char **xargv;

int main(int argc, char ** argv) {

  try {

    PARARRAY * pars = ps_new();
    
    if (ps_createargs(pars,argc-1,&(argv[1]))) {
      RVLException e;
      e<<"ERROR: taper from ps_creatargs \n";
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
    ddin.eval(ddinfn);
        
    AssignFilename ddoutfn(outp);
    ddout.eval(ddoutfn);

    SEGYTaper fo(valparse<std::string>(*pars,"taperpars"));

    LinearOpFO<float> op(dom,dom,fo,fo);
    op.applyOp(ddin,ddout);
  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
