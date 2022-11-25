#include "parser.h"
#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#else
#include "gridpp.hh"
#endif
#include "gridops.hh"
#include "par.h"

using RVL::valparse;
using RVL::RVLException;
using RVL::Vector;
using RVL::LinearOp;
using RVL::AssignFilename;

using TSOpt::GridDerivOp;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
typedef TSOpt::MPIGridSpace gsp;
#else
using TSOpt::GridSpace;
typedef TSOpt::GridSpace gsp;
#endif

int xargc;
char **xargv;

const char * sdoc[] = {
  "usage: GridRand.x in=<string> out=<string>",
  "in is name of RSF file pair defining GridSpace. out is name of  ",
  "RSF file pair storing random vector in this space.",
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
      e<<"ERROR: GridRand from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++) 
	e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }
    // since the product of grid spaces is not really an 
    // out-of-core structure, this driver operates on single
    // grid spaces
    string in = valparse<string>(*pars,"in");
    string out = valparse<string>(*pars,"out");

    gsp sp(in,"notype",true
#ifdef IWAVE_USE_MPI 
	   , retrieveGlobalComm()
#endif
	   );
    Vector<float> vec(sp);
    AssignFilename af(out);
    vec.eval(af);
    RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
    vec.eval(rnd);
    
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


