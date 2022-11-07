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
  "usage: GridLinComb.x a=<float> in1=<string> b=<float> in2=<string>",
  "strings are names of RSF file pairs. The grid geometries defined",
  "defined by the header files must be the same. The command modifies",
  "the second vector argument: in2 = a*in1 + b*in2.",
  NULL};

int main(int argc, char ** argv) {

  try {
#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);   
    storeGlobalComm(MPI_COMM_WORLD);
#endif

    PARARRAY * pars = ps_new();

    if (retrieveGlobalRank()==0 && argc<2) {
      pagedoc();
      exit(0);
    }
    
    if (ps_createargs(pars,argc-1,&(argv[1]))) {
      RVLException e;
      e<<"ERROR: GridLinComb from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++) 
	e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }
    // since the product of grid spaces is not really an 
    // out-of-core structure, this driver operates on single
    // grid spaces
    float a = valparse<float>(*pars,"a");
    string in1 = valparse<string>(*pars,"in1");
    float b = valparse<float>(*pars,"b",1.0f);
    string in2 = valparse<string>(*pars,"in2");

    gsp sp(in1,"notype",true
#ifdef IWAVE_USE_MPI 
	   , retrieveGlobalComm()
#endif
	   );
    Vector<float> vec1(sp);
    Vector<float> vec2(sp);
    AssignFilename af1(in1);
    AssignFilename af2(in2);
    vec1.eval(af1);
    vec2.eval(af2);
    vec2.linComb(a, vec1, b);
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


