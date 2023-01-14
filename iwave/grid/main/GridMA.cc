#include "parser.h"
#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#else
#include "gridpp.hh"
#endif
#include "gridops.hh"

using RVL::valparse;
using RVL::RVLException;
using RVL::Vector;
using RVL::LinearOp;
using RVL::AssignFilename;

using TSOpt::GridMAFO;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
typedef TSOpt::MPIGridSpace gsp;
#else
using TSOpt::GridSpace;
typedef TSOpt::GridSpace gsp;
#endif

const char * sdoc[] = {
  "Moving average in dimension 1, 2, or 3. Input averaging radius",
  "for each axis - default is 0 (i.e. no averaging). Similar in function ",
  "to sfsmooth.",
  NULL};
  
int xargc;
char **xargv;

int main(int argc, char ** argv) {

  try {
#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);   
    storeGlobalComm(MPI_COMM_WORLD);
#endif

    PARARRAY * pars = ps_new();

    if (ps_createargs(pars,argc-1,&(argv[1]))) {
      RVLException e;
      e<<"ERROR: GridDot from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++) 
	e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }
    // since the product of grid spaces is not really an 
    // out-of-core structure, this driver operates on single
    // grid spaces
    string inp = valparse<string>(*pars,"inp");
    string outp = valparse<string>(*pars,"outp");
    
    gsp sp(inp,"notype",true
#ifdef IWAVE_USE_MPI 
	   , retrieveGlobalComm()
#endif
	   );
    Vector<float> vecin(sp);
    Vector<float> vecout(sp);
    AssignFilename afin(inp);
    AssignFilename afout(outp);
    vecin.eval(afin);
    vecout.eval(afout);
    
    IPNT rad;
    for (int i=1; i < RARR_MAX_NDIM+1; i++) {
      std::stringstream foo;
      foo<<i;
      std::string barf = "rad";
      rad[i-1] = valparse<int>(*pars,barf+foo.str(),0);
    }

    int rep = valparse<int>(*pars,"repeat",1);
    
    GridMAFO f(rep,rad);
    vecout.eval(f,vecin);
    
    ps_delete(&pars);
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
    exit(0);
  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
  
}


