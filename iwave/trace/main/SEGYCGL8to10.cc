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
  "  hc  = current iterate",
  "  hp = next iterate",
  "  rc  = current normal residual",
  "  rp = next normal residual",
  "  dc = current data residual",
  "  dp = next data residual",
  "  g = current gradient",
  "  p = current search dir",
  "  q = current p image",
  "on return: hp, rp, dp updated",
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
    // since the product of grid spaces is not really an 
    // out-of-core structure, this driver operates on single
    // grid spaces
    std::string hc = valparse<std::string>(*pars,"hc");
    std::string rc = valparse<std::string>(*pars,"rc");
    std::string dc = valparse<std::string>(*pars,"dc");
    std::string hp = valparse<std::string>(*pars,"hp");
    std::string rp = valparse<std::string>(*pars,"rp");
    std::string dp = valparse<std::string>(*pars,"dp");
    std::string g = valparse<std::string>(*pars,"g");
    std::string p = valparse<std::string>(*pars,"p");
    std::string s = valparse<std::string>(*pars,"s");
    std::string q = valparse<std::string>(*pars,"q");

    gsp dom(hc,"notype"
#ifdef IWAVE_USE_MPI 
	   , retrieveGlobalComm()
#endif
	   );
    gsp rng(dc,"notype"
#ifdef IWAVE_USE_MPI 
	   , retrieveGlobalComm()
#endif
	   );
    rng.write(cerr);

    Vector<float> hcvec(dom);
    Vector<float> rcvec(dom);
    Vector<float> dcvec(rng);    
    Vector<float> hpvec(dom);
    Vector<float> rpvec(dom);
    Vector<float> dpvec(rng);
    Vector<float> gvec(dom);
    Vector<float> pvec(dom);
    Vector<float> svec(rng);
    Vector<float> qvec(dom);
    AssignFilename afhc(hc);
    AssignFilename afrc(rc);
    AssignFilename afdc(dc);
    AssignFilename afhp(hp);
    AssignFilename afrp(rp);
    AssignFilename afdp(dp);
    AssignFilename afg(g);
    AssignFilename afp(p);
    AssignFilename afs(s);
    AssignFilename afq(q);    
    hcvec.eval(afhc);
    hpvec.eval(afhp);
    rcvec.eval(afrc);
    rpvec.eval(afrp);
    dcvec.eval(afdc);
    dpvec.eval(afdp);
    gvec.eval(afg);    
    pvec.eval(afp);
    svec.eval(afs);
    qvec.eval(afq);    

    float alpha = gvec.inner(rcvec)/pvec.inner(qvec);
    hpvec.copy(hcvec);
    hpvec.linComb(alpha,pvec);
    rpvec.copy(rcvec);
    rpvec.linComb(-alpha,qvec);
    dpvec.copy(dcvec);
    dpvec.linComb(-alpha,svec);
    
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


