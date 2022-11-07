#include "parser.h"
#ifdef IWAVE_USE_MPI
#include "mpisegypp.hh"
#else
#include "segypp.hh"
#endif

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
  "lines 8-10 of PCG alg",
  "on call:",
  "  rc = current data residual",
  "  rp = next data residual",  
  "  gc = current gradient",
  "  gp = next gradient",
  "  pc = current search dir",
  "  pp = next search dir",
  "on return: pp updated",
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
      e<<"ERROR: SEGYCGL8to10 from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++) 
	e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }

    std::string rc = valparse<std::string>(*pars,"rc");
    std::string rp = valparse<std::string>(*pars,"rp");
    std::string gc = valparse<std::string>(*pars,"gc");
    std::string gp = valparse<std::string>(*pars,"gp");
    std::string pc = valparse<std::string>(*pars,"pc");
    std::string pp = valparse<std::string>(*pars,"pp");

    gsp sp(rc,"notype"
#ifdef IWAVE_USE_MPI 
	   , retrieveGlobalComm()
#endif
	   );

    Vector<float> rcvec(sp);
    Vector<float> rpvec(sp);
    Vector<float> gcvec(sp);
    Vector<float> gpvec(sp);
    Vector<float> pcvec(sp);
    Vector<float> ppvec(sp);
    AssignFilename afrc(rc);
    AssignFilename afrp(rp);
    AssignFilename afgc(gc);
    AssignFilename afgp(gp);
    AssignFilename afpc(pc);
    AssignFilename afpp(pp);
    rcvec.eval(afrc);
    rpvec.eval(afrp);
    gcvec.eval(afgc);
    gpvec.eval(afgp);
    pcvec.eval(afpc);
    ppvec.eval(afpp);

    float beta = gpvec.inner(rpvec)/gcvec.inner(rcvec);
    ppvec.copy(gpvec);
    ppvec.linComb(beta,pcvec);
    
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


