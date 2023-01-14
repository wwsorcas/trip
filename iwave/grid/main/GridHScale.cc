#include "parser.h"
#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#else
#include "gridpp.hh"
#endif
#include "gridops.hh"
#include "adjtest.hh"
#include "functions.hh"
#include "par.h"

using RVL::valparse;
using RVL::RVLException;
using RVL::Vector;
using RVL::LinearOp;
using RVL::AssignFilename;

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
  "scales input by distance along internal extended axes (standard",
  "DSO annihilator).",
  " ",
  "written so that axes dim,...,gdim-1 are assumed to be internal ",
  "extended, and outputs out[i], i=0,...gdim-dim-1 need to be specified",
  "in case gdim=dim+1, out is the only output key",
  " ",
  "Thus works for 2D or 3D equally.",
  " ",
  "Arguments:",
  "  in = rsf filename for input image volume",
  "  out or out0, out1, ... = rsf filename for output volume (2D)",
  "  or volumes (3 or larger D)",
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
      e<<"ERROR: GridHScale from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++) 
	e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }
    string inp = valparse<string>(*pars,"in");
    gsp sp(inp,"notype",true
#ifdef IWAVE_USE_MPI 
	   , retrieveGlobalComm()
#endif
	   );

    int dim = sp.getGrid().dim;
    if (sp.getGrid().dim >= sp.getGrid().gdim) {
      RVLException e;
      e<<"ERROR: GridHScale\n";
      e<<"  input grid from file "<<inp<<" has gdim <= dim\n";
      throw e;
    }
    std::vector<std::string> outp;
    std::vector<float> a(sp.getGrid().gdim, 0.0f);
    if (sp.getGrid().gdim == dim + 1) {
      outp.push_back(valparse<string>(*pars,"out"));
    }
    else {
      for (int i=0; i< sp.getGrid().gdim - dim; i++) {
	std::stringstream nm;
	nm<<"out"<<i;
	outp.push_back(valparse<string>(*pars,nm.str()));
      }
    }

    Vector<float> invec(sp);
    AssignFilename afin(inp);
    invec.eval(afin);
    
    for (int i=0; i< sp.getGrid().gdim - dim; i++) {
      a[dim+i]=1.0f;      
      TSOpt::GridCoordScaleOp op(sp,a,0.0f);
      Vector<float> outvec(sp);
      AssignFilename afout(outp[i]);
      outvec.eval(afout);
      op.applyOp(invec,outvec);
      a[dim+i]=1.0f;
    }
    ps_delete(&pars);
    iwave_fdestroy();
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
    exit(0);
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


