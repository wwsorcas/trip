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
#include "adjtest.hh"
#include "segyops.hh"
#include "MPS_iwop.hh"
#include "MPS_Space_Examples.hh"
#include "MPS_frac_cal.hh"
#include "MPS_spread.hh"


IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csq_inv", 0, 1, 1 },
  {"data",    1, 0, 2 },
  {"source",  1, 1, 2 },
  {"",        0, 0, 0 }
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

using TSOpt::IWaveEnvironment;
using TSOpt::IWaveTree;
using TSOpt::IWaveSampler;
using TSOpt::IWaveSim;
using TSOpt::TASK_RELN;
using TSOpt::IOTask;
using TSOpt::IWaveOp;
using TSOpt::SEGYLinMute;
using TSOpt::MPS_IWaveLOVOp;
using TSOpt::Scal_MPS_Space;
using TSOpt::MPS_frac_cal;
using TSOpt::MPS_spread;
using TSOpt::MPS_KEYS;

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
using RVLUmin::ReduceOp;

const char *sdoc[] = {
  " ============================================================================",
  " acd_scal_Jinv_towed.x ",
  " ============================================================================",
  " Authors: Mario J. Bencomo",
  " ",
  " Joint source-medium inversion driver with acoustic constant density modeling",
  " single scalar multipole sources and pressure data output. Forward modeling  ",
  " is done via convolution with precomputed Green's functions (see MPS_IWaveLOVOp",
  " operator). This version allows for multiple shots from a towed multipole source",
  " ,i.e., many sources having the same source waveforms (see MPS_spread operator).",
  " ",
  " Forward map F(m), ",
  "       F(m)w=d",
  "       m = medium parameters (bulk modulus)",
  "       w = MPS coefficients (pressure source)",
  "       d = data (sampled pressure field)",
  " FWI formulation of source estimation:",
  "       min_{m,w} 0.5*|F(m)w - d|^2",
  " ",
  " The variable projection method is applied thus reducing the forward map F to a",
  " function of medium parameters only. In particular, the reduced forward map F_red",
  " is defined as"
  "       F_red(m) = F(m)w(m) ",
  " where w(m) = argmin_{w} 0.5*|F(m)w - d|^2",
  " The medium parameters are then determined by solving the LS problem with the",
  " reduced forward map:",
  "       min_{m} F_red(m)",
  " The variable projection method effectively eliminates the (linear) source",
  " parameters in the inversion. Equivalently, we implement variable projection",
  " as a nested optimization scheme where the source and medium parameters",
  " are updated in the inner and outer optimization loops respectively. More",
  " precisely, the inner optimization loop is ran for every evaluation of the",
  " reduced forward map if the medium parameter has been modified."
  " Trust region Gauss-Newton is used for the (nonlinear) medium parameter",
  " estimation subproblem, with conjugate gradient for the underlying linear solves.",
  " Conjugate gradient is also used to solve normal equations related to source estimation.",
  " Fractional calculus operators are available for preconditioning.",
  " ",
  " ",
  " Typical parameter list. May be copied, edited, and used for input: either",
  " include parameters on command line (for example in Flow), or place",
  " in file <foo> and include \"par=<foo>\" on command line. Any parameter",
  " included explicitly in command line overrides parameter with same key",
  " in par file.",
  " ",
  "  Invoke single threaded execution by ",
  " \"acd_scal_Jinv_towed.x [parameters] [standalone install]\"",
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
  "        csq_inv = <path>      output estimated velocity squared model, ",
  "                              also used as initial estimate",
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
  "  MPS_file_inv = <path>     SU file containing MPS coefficients and source locations",
  "                            (if not at towed source).",
  "     MPS_towed = <path>     SU file containing MPS coefficients of a towed ",
  "                            source specified by MPS_file, SU format",
  "                            Must exist and contain source locations, but need",
  "                            not contain actual data. MPS_file and MPS_towed ",
  "                            must be compatible (spanned by same MPS basis with",
  "                            time axes.",
  "     MPS_delta = <path>     MPS delta file, SU format - must exist",
  "      data_p_g = <path>     Green's function file or 'Green data', SU format - must exist",
  " ------------------------------------------------------------------------",
  " CG parameters for MPS inversion:",
  " ",
  "           MPS_RTol = 1e3*eps_mach  residual tolerance",
  "           MPS_GTol = 1e3*eps_mach  gradient (normal residual) tolerance",
  "       MPS_TrustRad = eps_max       trust region radius",
  "        MPS_MaxIter = 10            maximum number of CG iterations",
  "        MPS_Verbose = 0             verbose output flag",
  " ",
  " Preconditioner info:",
  " ",
  "            order_0 = 0  preconditioner base order",
  "            order_d = 0  preconditioner incerement per MPS order",
  " ",
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
	 << "Running acd_scal_Jinv_towed.x\n\n";

    // Initializing MPS Spaces and MPS spread op
    MPS_KEYS mk_sp;
    mk_sp.MPS_file  = "MPS_file_inv";
    mk_sp.grid_file = "csq_inv";
    Scal_MPS_Space MPS_SP(mk_sp,*pars);    

    MPS_KEYS mk_tw;
    mk_tw.MPS_file  = "MPS_towed";
    mk_tw.grid_file = "csq_inv"; 
    Scal_MPS_Space MPS_TW(mk_tw,*pars);

    MPS_spread MPS_SPR(MPS_SP,MPS_TW);

    // Initializing forward map F
    MPS_KEYS mk;
    mk.MPS_file   = "MPS_towed";
    mk.grid_file  = "csq_inv";     
    mk.appx_ord   = "appx_ord";
    mk.delta_file = "MPS_delta";
    mk.RHS_files.push_back("source");
    mk.hdr_files.push_back("data");
    mk.G_files.push_back("data_g");    
    MPS_IWaveLOVOp<Scal_MPS_Space> MPS_LOVO(mk,*pars,stream);
    LinCompLOVOp<float> F(MPS_SPR,MPS_LOVO);


    // Pairing input/output vectors with file names
    Vector<float>  m(F.getProductDomain()[0]);
    Vector<float>  w(F.getProductDomain()[1]);
    Vector<float>  d(F.getRange());
    Vector<float>  r(F.getRange());

    AssignFilename af_csq_inv(valparse<string>(*pars,"csq_inv"));
    AssignFilename af_mps(valparse<string>(*pars,"MPS_file_inv"));
    AssignFilename af_data(valparse<string>(*pars,"data"));
    AssignFilename af_res(valparse<string>(*pars,"datares"));

    m.eval(af_csq_inv);
    w.eval(af_mps);
    d.eval(af_data);
    r.eval(af_res);

    /*
    // Initializing preconditioner op M
    float c       = valparse<float>(*pars,"cmax",1.0f);
    float order_0 = valparse<float>(*pars,"order_0",0);
    float order_d = valparse<float>(*pars,"order_d",0);
    bool precond  = (order_0!=0) || (order_d!=0);
   
    MPS_frac_cal Q_inv(*pars,
                       mps_space,
		       c,
		       order_0,
		       order_d,
		       true );
    
    AdjLinearOp<float> Q_inv_adj(Q_inv);
    NormalLinearOp<float> M(Q_inv_adj);
    */


    // Initializing CGNE policy for source inversion 
    // and TRGN algorithm for medium inversion

    float MPS_RTol = valparse<float>(*pars,"MPS_RTol",1000.0*numeric_limits<float>::epsilon());
    float MPS_GTol = valparse<float>(*pars,"MPS_GTol",1000.0*numeric_limits<float>::epsilon());
    float MPS_TRad = valparse<float>(*pars,"MPS_TrustRad",numeric_limits<float>::max());
    int   MPS_MaxI = valparse<int>(*pars,"MPS_MaxIter",10);
    int   MPS_Verb = valparse<int>(*pars,"MPS_Verbose",0);

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
	   <<"* Acoustic constant density (acd) joint source-medium    \n"
	   <<"* inversion via trust region Gauss-Newton for medium     \n"
	   <<"* parameters with conjugate gradient for normal eqns     \n"
	   <<"* algorithm for the underlying linear solves. Conjugate  \n"
	   <<"* gradient for normal eqns algorithm is also used for    \n"
	   <<"* inner inversion of source parameters.                  \n"
	   <<"*                                                        \n"
	   <<"* CGNE parameters for source inversion:                  \n"
	   <<"*   max iterations       = "<< MPS_MaxI  <<endl
	   <<"*   residual tolerance   = "<< MPS_RTol  <<endl
	   <<"*   normal res tolerance = "<< MPS_GTol  <<endl
	   <<"*   trust radius         = "<< MPS_TRad  <<endl
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

    CGNEPolicyData<float> cgne(MPS_RTol,
                               MPS_GTol,
                               MPS_TRad,
                               MPS_MaxI,
                               MPS_Verb);
    
    //ReduceOp_mod<float,CGNEPolicy<float>,CGNEPolicyData<float> > rop(F,M,d,cgne,sstream);
    ReduceOp<float,CGNEPolicy<float>,CGNEPolicyData<float> > R_red(F,d,cgne,sstream);
    
    RVLAlg::Algorithm * alg = NULL;
    RVLUmin::TRGNAlg<float, RVLUmin::CGNEPolicy<float> > * tralg = 
      new RVLUmin::TRGNAlg<float, RVLUmin::CGNEPolicy<float> >
      (R_red, m,
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

    // recomputing data residual and source
    OperatorEvaluation<float> opeval(R_red,m);
    r.copy(opeval.getValue());
    ReduceOp<float,CGNEPolicy<float>,CGNEPolicyData<float> > const & rop
      = dynamic_cast< ReduceOp<float,CGNEPolicy<float>,CGNEPolicyData<float> > const &>(opeval.getOp());
    w.copy(rop.getLSSoln());


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
    
    cerr << "\nFinishing acd_scal_Jinv_towed.x\n"
	 << "////////////////////////////////////////////////////////////\n";

#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
  }
  catch (bad_cast) {
    RVLException e;
    e<<"Error: acd_scal_Jinv_towed.x\n";
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
