// acd_scal_Jinv_towed.cc
// Author: Mario J. Bencomo
// last modified: 04/02/18

#include "acd_defn.hh"
#include "grid.h"
#include "iwop.hh"
#include "iwopt.hh"
#include "gridpp.hh"
#include "functions.hh"
#include "cgnealg.hh"
#include "RedOp.hh"
//#include "RedOp_mod.hh"
#include "linop_apps.hh"
#include "adjtest.hh"
#include "segyops.hh"
#include "MPS_iwop.hh"
#include "MPS_Space_Examples.hh"
#include "MPS_frac_cal.hh"
#include "MPS_spread.hh"


IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csq",    0, 1, 1 },
  {"data",   1, 0, 2 },
  {"source", 1, 1, 2 },
  {"",       0, 0, 0 }
};

using RVL::parse;
using RVL::RVLException;
using RVL::Vector;
using RVL::Components;
using RVL::Operator;
using RVL::OperatorEvaluation;
using RVL::LinearOp;
using RVL::LinearOpFO;
using RVL::OpComp;
using RVL::SymmetricBilinearOp;
using RVL::AssignFilename;
using RVL::AssignParams;
using RVL::RVLRandomize;
using RVL::AdjointTest;
using RVL::LinCompLOVOp;
using RVL::CompLinearOp;

using TSOpt::IWaveEnvironment;
using TSOpt::IWaveTree;
using TSOpt::IWaveSampler;
using TSOpt::IWaveSim;
using TSOpt::TASK_RELN;
using TSOpt::IOTask;
using TSOpt::IWaveOp;
using TSOpt::IWaveLOVOp;
using TSOpt::SEGYLinMute;
using TSOpt::MPS_IWaveLOVOp;
using TSOpt::Scal_MPS_Space;
using TSOpt::MPS_frac_cal;
using TSOpt::MPS_spread;
using TSOpt::MPS_KEYS;
using TSOpt::MPS_to_RHS;

#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
using TSOpt::MPISEGYSpace;
#else
using TSOpt::GridSpace;
using TSOpt::SEGYSpace;
#endif

using RVLUmin::CGNEAlg;
using RVLUmin::CGNEPolicy;
using RVLUmin::CGNEPolicyData;
//using RVLUmin::ReduceOp_mod;
using RVLUmin::ReduceOp;

const char *sdoc[] = {
  " ============================================================================",
  " acd_trgn.x ",
  " ============================================================================",
  " Authors: Mario J. Bencomo",
  " ",
  " Medium inversion driver with acoustic constant density modeling and single ",
  " scalar multipole sources and pressure data output. This version allows for ",
  " multiple shots from a towed multipole source, i.e., many sources having the",
  " same source waveforms (see MPS_spread operator).",
  " ",
  " Trust region Gauss-Newton is used for the underlying nonlinear optimization",
  " problem, with conjugate gradient for the undelying linear solves.",
  " ",
  " Forward map F, F(m) w = d",
  "       m = medium parameters (bulk modulus)",
  "       w = MPS coefficients (pressure source)",
  "       d = data (sampled pressure field)",
  " FWI formulation for medium parameter estimation:",
  "       min_{m} 0.5*|F(m)w - d|^2",
  " ",
  " Typical parameter list. May be copied, edited, and used for input: either",
  " include parameters on command line (for example in Flow), or place",
  " in file <foo> and include \"par=<foo>\" on command line. Any parameter",
  " included explicitly in command line overrides parameter with same key",
  " in par file.",
  " ",
  "  Invoke single threaded execution by ",
  " \"acd_trgn.x [parameters] [standalone install]\"",
  " ",
  " or multi-threaded execution using interactive or batch MPI (for which",
  " code must be compiled with MPI enabled).",
  " ",
  " Given values are defaults; non-optional values indicated by corner brackets.",
  " ",
  " --------------------------- begin parameters ---------------------------",
  " FD info:",
  " ",
  "          order = 2           spatial half-order",
  "            cfl = 0.75        proportion of max dt/dx",
  "           cmin = 1.0         min permitted velocity (m/ms) - sanity check",
  "           cmax = 4.5         max permitted velocity (m/ms) - used in dt comp",
  " ",
  " ------------------------------------------------------------------------",
  " Source info:",
  " ",
  "         source=empty         key required related to RHS sources for acd driver,",
  "                              set to empty to generate filename from code",
  " ",
  " ------------------------------------------------------------------------",
  " Trace info:",
  " ",
  "          data  = <path>      data file, SU format - must exist ",
  "                              headers establish acquisition geometry",
  "                              output for forward modeling (including all",
  "                              orders of derivative), input for adjoint",
  " ",
  " ------------------------------------------------------------------------",
  " Model info:",
  " ",
  "        csq = <path>          input velocitu squared model, ",
  "                              determines simulation spatial grid,",
  " ",
  " ------------------------------------------------------------------------",
  " MPI info:",
  " ",
  "        mpi_np1 = 1           number of subdomains along axis 1",
  "        mpi_np2 = 1           number of subdomains along axis 2",
  "        mpi_np3 = 1           number of subdomains along axis 3",
  "        partask = 1           number of shots to execute in parallel",
  " ",
  " ------------------------------------------------------------------------",
  " Output info:",
  " ",
  " FD ouput - written to coutxxx.txt on rank xxx",
  "       printact = 0           monitor output:",
  "                            < 0 - none",
  "                              0 - announce each shot simulation",
  "                              1 - print time step index",
  "                              2 - diagnostic messages from main routines",
  "                              > 2 - much more, useful only for serious ",
  "                                 debugging",
  "        dump_pi = 0           dump parallel/dom. decomp info",
  "       dump_lda = 0           dump grid data for allocated arrays",
  "       dump_ldc = 0           dump grid data for computational arrays",
  "       dump_ldr = 0           dump grid data for receive buffers",
  "       dump_lds = 0           dump grid data for send buffers",
  "      dump_term = 0           dump trace header data",
  " ",
  " ------------------------------------------------------------------------",
  " MPS info:",
  " ",
  "      appx_ord = 4          singular source approximation order ",
  "     MPS_ord_0 = 0          MPS order in z-direction ",
  "     MPS_ord_1 = 0          MPS order in x-direction ",
  "     MPS_ord_2 = 0          MPS order in y-direction ",
  "     MPS_towed = <path>     SU file containing MPS coefficients of a towed ",
  "                            source, SU format",
  "                            Must exist and contain source locations.",
  "     MPS_delta = <path>     MPS delta file, SU format - must exist",
  "      data_p_g = <path>     Green's function file or 'Green data', SU format - must exist",
  " ------------------------------------------------------------------------",
  " TRGN parameters for medium inversion:",
  " ",
  "         GN_MaxIter = 5             maximum number of GN iterations",
  "     GN_ResidualTol = 1e3*eps_mach  residual tolerance",
  "      GN_AbsGradTol = 1e3*eps_mach  absolute gradient residual tolerance",
  "      GN_RelGradTol = 0.01          relative gradient residual tolerance",
  "     GN_MinDecrease = 0.1           lower G-A parameter",
  "    GN_GoodDecrease = 0.9           upper G-A parameter",
  "  GN_StepDecrFactor = 0.5           trust region reduction factor",
  "  GN_StepIncrFactor = 1.8           trust region expansion factor",
  " GN_InitTrustRadius = 1e-6          initial trust region radius",
  " ",
  " CG parameters for medium inversion:",
  " ",
  "            CG_RTol = 1e3*eps_mach residual tolerance",
  "            CG_GTol = 1e3*eps_mach  gradient (normal residual) tolerance",
  "        CG_TrustRad = 1.0           trust region radius",
  "         CG_MaxIter = 10            maximum number of CG iterations",
  "         CG_Verbose = 0             verbose output flag",
  " ",
  " ---------------------------end parameters ------------------------------",
  NULL };


int xargc;
char **xargv;

int main(int argc, char ** argv) {
  stringstream sstream;

  try {

#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
#endif

    PARARRAY * pars = NULL;
    FILE * stream = NULL;
    IWaveEnvironment(argc, argv, 0, &pars, &stream);

    cerr << "\n////////////////////////////////////////////////////////////\n"
	 << "Running acd_trgn.x\n\n";

    // Setting MPS keys
    MPS_KEYS mks;
    mks.appx_ord  = "appx_ord";
    mks.grid_file = "csq";
    mks.RHS_files.push_back("source");
    mks.hdr_files.push_back("data");
    mks.MPS_file = "MPS_towed";

    // Initializing MPS_to_RHS op
    Scal_MPS_Space mps_sp(mks,*pars);
    MPS_to_RHS m2r(mks,*pars,mps_sp);
    
    // Initializing the forward map F
    IWaveLOVOp lovop(m2r.get_pars(),stream);
    LinCompLOVOp<float> F(m2r,lovop);    
      
    // Pairing vectors with file names
    Vector<float> m0(F.getProductDomain()[0]);
    Vector<float>  m(F.getProductDomain()[0]);
    Vector<float>  w(F.getProductDomain()[1]);
    Vector<float>  d(F.getRange());
    AssignFilename af_csq(valparse<string>(*pars,"csq"));
    AssignFilename af_csq_inv(valparse<string>(*pars,"csq_inv"));
    AssignFilename af_mps(valparse<string>(*pars,"MPS_towed"));
    AssignFilename af_data(valparse<string>(*pars,"data"));
    m0.eval(af_csq);
    m.eval(af_csq_inv);
    w.eval(af_mps);
    d.eval(af_data);

    Vector<float> x(F.getDomain());
    Components<float> cx(x);
    cx[0].copy(m);
    cx[1].copy(w);

    RestrictOp<float> F0(F,x,0); //op with fixed correct source


    // Initializing TRGN algorithm for medium inversion

    int   GN_MaxI  = valparse<int>(*pars,"GN_MaxIter",5); 
    float GN_RTol  = valparse<float>(*pars,"GN_ResidualTol",1000.0*numeric_limits<float>::epsilon());
    float GN_AGTol = valparse<float>(*pars,"GN_AbsGradTol",1000.0*numeric_limits<float>::epsilon());
    float GN_RGTol = valparse<float>(*pars,"GN_RelGradTol",1.0e-2);
    float GN_Eta1  = valparse<float>(*pars,"GN_MinDecrease",0.1f);   
    float GN_Eta2  = valparse<float>(*pars,"GN_GoodDecrease",0.9f);
    float GN_Gamma1= valparse<float>(*pars,"GN_StepDecrFactor",0.5f);
    float GN_Gamma2= valparse<float>(*pars,"GN_StepIncrFactor",1.8f);
    float GN_TRad0 = valparse<float>(*pars,"GN_InitTrustRadius",1.e-06);

    float CG_RTol = valparse<float>(*pars,"CG_RTol",1000.0*numeric_limits<float>::epsilon());
    float CG_GTol = valparse<float>(*pars,"CG_GTol",0.001);
    float CG_TRad = valparse<float>(*pars,"CG_TrustRad",1.f);
    int   CG_MaxI = valparse<int>(*pars,"CG_MaxIter",10);
    int   CG_Verb = valparse<int>(*pars,"CG_Verbose",0);

    sstream<<scientific;
    sstream<<"\n*******************************************************\n"
	   <<"* Acoustic constant density (acd) medium inversion via   \n"
	   <<"* trust region Gauss-Newton with conjugate gradient for  \n"
	   <<"* normal eqns algorithm for the underlying linear solves.\n"
	   <<"*                                                        \n"
	   <<"* TRGN parameters for medium inversion:                  \n"
	   <<"*   max itertations         = "<< GN_MaxI   <<endl
	   <<"*   residual tolerance      = "<< GN_RTol   <<endl
	   <<"*   absolute normal res tol = "<< GN_AGTol  <<endl
	   <<"*   relative normal res tol = "<< GN_RGTol  <<endl
	   <<"*   minimum decrease        = "<< GN_Eta1   <<endl
	   <<"*   good decrease           = "<< GN_Eta2   <<endl
	   <<"*   trust region reduction  = "<< GN_Gamma1 <<endl
	   <<"*   trust region expansion  = "<< GN_Gamma2 <<endl
	   <<"*   initial trust radius    = "<< GN_TRad0  <<endl
	   <<"*                                                        \n"
	   <<"* CGNE parameters for medium inversion:                  \n"
	   <<"*   max iterations       = "<< CG_MaxI  <<endl
	   <<"*   residual tolerance   = "<< CG_RTol  <<endl
	   <<"*   normal res tolerance = "<< CG_GTol  <<endl
	   <<"*   trust radius         = "<< CG_TRad  <<endl
	   <<"*******************************************************\n";

    RVLAlg::Algorithm * alg = NULL;
    RVLUmin::TRGNAlg<float, RVLUmin::CGNEPolicy<float> > * tralg = 
      new RVLUmin::TRGNAlg<float, RVLUmin::CGNEPolicy<float> >
      (F0, m,
       GN_MaxI,
       GN_RTol,
       GN_AGTol,
       GN_RGTol,
       GN_Eta1,
       GN_Eta2,
       GN_Gamma1,
       GN_Gamma2,
       GN_TRad0,
       sstream);

    tralg->assign(CG_RTol,
		  CG_GTol,
		  CG_TRad,
		  CG_MaxI,
		  CG_Verb);

    alg=tralg;
    alg->run();

    // recomputing data residual
    Vector<float> r(F0.getRange());
    AssignFilename af_res(valparse<string>(*pars,"datares"));
    r.eval(af_res);

    OperatorEvaluation<float> opeval(F0,m);
    r.copy(opeval.getValue());


#ifdef IWAVE_USE_MPI
    if (retrieveGlobalRank()==0) {
#endif
      string outfile = valparse<string>(*pars,"outfile","");
      if (outfile.size()>0) {
        ofstream outf(outfile.c_str());
        outf<<sstream.str();
        outf.close();
      }
      else {
        cout<<sstream.str();
      }
#ifdef IWAVE_USE_MPI
    }
#endif

    //clean up
    ps_delete(&pars);
    iwave_fdestroy();
    
    cerr << "\nFinishing acd_trgn.x\n"
	 << "////////////////////////////////////////////////////////////\n";

#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
  }
  catch (bad_cast) {
    RVLException e;
    e<<"Error: acd_trgn.x\n";
    e<<"  op returned from opeval is not ReduceOp - WTF?\n";
    throw e;
  }
  catch (RVLException & e) {
    e << "Exiting with error!\n"
      << "Dumping output stream:\n";
#ifdef IWAVE_USE_MPI
    if( retrieveGlobalRank()==0 ){
#endif
      e << sstream.str();
      e.write(cerr);
#ifdef IWAVE_USE_MPI
    }
    MPI_Barrier(retrieveGlobalComm());
    MPI_Abort(retrieveGlobalComm(),0);
#endif
    exit(1);
  }
}
