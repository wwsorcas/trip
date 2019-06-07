#include "parser.h"
#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#else
#include "gridpp.hh"
#endif
#include "gridops.hh"
#include "gridfftops.hh"
#include "adjtest.hh"
#include "functions.hh"
#include <par.h>

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
  "GridZFTScale: scale columns (depth traces) of multi-d grid by power of  ",
  "frequency, followed by cosine shoulder bandpass filter.",
  "inputs:",
  "in       = filename of input rsf file",
  "out      = filename of output rsf file",
  "power    = power of frequency or absolute frequency",
  "absfreq  = 0 for power of frequency, 1 for power of absolute frequency",
  "locut    = low cut frequency for bandpass",
  "lopas    = low pass frequency for bandpass",
  "hipas    = high pass frequency for bandpass",
  "hicut    = hich cut frequency for bandpass",
  " ",
  "Note: frequency should be given in appropriate units. For example, if",
  "depth is in meters, frequency is in cycles/meter, generally in the range",
  "0.001-0.1",
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
      RVL::RVLException e;
      e<<"ERROR: GridDerivOp from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++) 
	e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }
    // since the product of grid spaces is not really an 
    // out-of-core structure, this driver operates on single
    // grid spaces
    string inp    = RVL::valparse<string>(*pars,"in");
    string outp   = RVL::valparse<string>(*pars,"out");
    int ab        = RVL::valparse<int>(*pars,"absfreq");
    float pow     = RVL::valparse<float>(*pars,"power");
    float locut   = RVL::valparse<float>(*pars,"locut");
    float lopas   = RVL::valparse<float>(*pars,"lopas");
    float hipas   = RVL::valparse<float>(*pars,"hipas");
    float hicut   = RVL::valparse<float>(*pars,"hicut");

    gsp sp(inp,"notype",true
#ifdef IWAVE_USE_MPI 
	   , retrieveGlobalComm()
#endif
	   );
    TSOpt::GridZFTScaleOp op(sp,pow,ab,locut,lopas,hipas,hicut);
    RVL::Vector<float> invec(sp);
    RVL::Vector<float> outvec(sp);
    RVL::AssignFilename afin(inp);
    RVL::AssignFilename afout(outp);
    invec.eval(afin);
    outvec.eval(afout);
    op.applyOp(invec,outvec);
    ps_delete(&pars);
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
  }
  catch (RVL::RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
  
}


