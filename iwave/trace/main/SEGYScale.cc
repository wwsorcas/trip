#include "parser.h"
#ifdef IWAVE_USE_MPI
#include "mpisegypp.hh"
#else
#include "segypp.hh"
#endif
#include "par.h"

using RVL::valparse;
using RVL::RVLException;
using RVL::Vector;
using RVL::AssignFilename;

#ifdef IWAVE_USE_MPI
using TSOpt::MPISEGYSpace;
typedef TSOpt::MPISEGYSpace gsp;
#else
using TSOpt::SEGYSpace;
typedef TSOpt::SEGYSpace gsp;
#endif

int xargc;
char **xargv;

const char * sdoc[] = {
  "usage: SEGYScale.x in=<string> scale=<float>",
  "string is name of SU file. scale is factor applied to every trace sample",
  NULL};

int main(int argc, char ** argv) {

  try {
#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);   
    storeGlobalComm(MPI_COMM_WORLD);
#endif

    if (retrieveGlobalRank()==0 && argc<2) {
      pagedoc();
      exit(0);
    }    
    
    PARARRAY * pars = ps_new();
    
    if (ps_createargs(pars,argc-1,&(argv[1]))) {
      RVLException e;
      e<<"ERROR: SEGYNorm from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++) 
	e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }

    string in = valparse<string>(*pars,"in");
    float a = valparse<float>(*pars,"scale");

    gsp sp(in,"notype"
#ifdef IWAVE_USE_MPI 
	   , retrieveGlobalComm()
#endif
	   );
    Vector<float> vec(sp);
    AssignFilename af(in);
    vec.eval(af);
    vec.scale(a);
    ps_delete(&pars);
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
  
}


