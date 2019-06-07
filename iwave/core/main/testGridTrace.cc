#include "parser.h"
#include "segypp.hh"
#include "gridpp.hh"
#include "gridtrace.hh"

int xargc;
char ** xargv;

int main(int argc, char ** argv) {

  try {

    PARARRAY * pars = ps_new();
    ps_createargs(pars,argc-1,argv+1);
    
    FILE * stream = fopen("cout0.txt","w");
    
    // grid = input
    std::string inp = RVL::valparse<string>(*pars,"input");
    std::string outp = RVL::valparse<string>(*pars,"output");
    std::string datap = RVL::valparse<string>(*pars,"data");
    std::string cwproot = getenv("CWPROOT");
    createExplReflTraces(inp,outp,datap,cwproot,stream);

    ps_delete(&pars);
    
    fflush(stream);
    fclose(stream);
  }
  catch (RVL::RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
