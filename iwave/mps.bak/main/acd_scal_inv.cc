// acd_scal_inv.cc
// Author: Mario J. Bencomo
// last modified: 04/10/18

#include "acd_defn.hh"
#include "acd.hh"
#include "functions.hh"
#include "cgnealg.hh"
#include "MPS_iwop.hh"
#include "MPS_Space_Examples.hh"
#include "MPS_frac_cal.hh"
#include "MPS_spread.hh"

//#define VERBOSE_MJB

using RVL::valparse;
using RVL::RVLException;
using RVL::Vector;
using RVL::Components;
using RVL::OperatorEvaluation;
using RVL::FunctionalEvaluation;
using RVL::AssignFilename;
using RVL::NormalLinearOp;
using RVL::LinCompLOVOp;
using RVL::AdjLinearOp;
using RVL::InvLinearOp;

using RVLUmin::CGNEPolicy;
using RVLUmin::CGNEPolicyData;

using TSOpt::IWaveEnvironment;
using TSOpt::MPS_IWaveLOVOp;
using TSOpt::CanScal_MPS_Space;
using TSOpt::CanVec_MPS_Space;
using TSOpt::ExVec_MPS_Space;
using TSOpt::Scal_MPS_Space;
using TSOpt::MPS_frac_cal;
using TSOpt::MPS_spread;
using TSOpt::MPS_KEYS;

const char *sdoc[] = {
  " ============================================================================",
  " acd_scal_inv.x ",
  " ============================================================================",
  " Authors: Mario J. Bencomo",
  " ",
  " Source inversion driver with acoustic constant density centered difference",
  " modelling for single scalar multipole sources given pressure data output.",
  " Forward modeling is done via convolution with precomputed Green's functions,",
  " See MPS_IWaveLOVOp operator.",
  "",
  " Conjugate gradient is used to solve normal equations related to source estimation. ",
  " Fractional calculus operators are available for preconditioning.",
  " ",
  " Forward map F(m), ",
  "       F(m)w = d",
  "       m = medium parameters (bulk modulus)",
  "       w = MPS coefficients (pressure source)",
  "       d = data (sampled pressure field)",
  " FWI formulation of source estimation:",
  "       min_{w} 0.5*|F(m)w - d|^2",
  " Normal equations:",
  "       F(m)^T ( F(m)w - d ) = 0",
  " Preconditioned normal equations:",
  "     M F(m)^T ( F(m)w - d ) = 0",
  "     M = (Q^T Q)^{-1} with frac calculus op Q",
  " ",
  " Typical parameter list. May be copied, edited, and used for input: either",
  " include parameters on command line (for example in Flow), or place",
  " in file <foo> and include \"par=<foo>\" on command line. Any parameter",
  " included explicitly in command line overrides parameter with same key",
  " in par file.",
  " ",
  "  Invoke single threaded execution by ",
  " \"acd_scal_inv.x [parameters] [standalone install]\"",
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
  "  MPS_file_inv = <path>     SU file containing MPS coefficients and source locations",
  "                            for estimates source.",
  "     MPS_delta = <path>     SU file containing filtered delta waveform to compute Green's functions.",
  "      data_p_g = <path>     SU file for formating Green's function output, or 'Green data'.",
  " ------------------------------------------------------------------------",
  " CG info:",
  " ",
  "   CG_MaxIter = 10            Max number of CG iterations",
  "      CG_RTol = 1e3*eps_mach  Residual tolerance",
  "      CG_GTol = 1e3*eps_mach  Gradient (normal residual) tolerance",
  "  CG_TrustRad = eps_max       Trust region radius",
  "      outfile = cerr          Text file to dump CG iterate info", 
  " ",
  " Preconditioner info:",
  "          order_0 = (optional)  preconditioner base order",
  "          order_d = (optional)  preconditioner incerement per MPS order",
  "",
  " ---------------------------end parameters ------------------------------",
  NULL };

IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csq",    0, true,  1 },
  {"data",   1, false, 2 },
  {"source", 1, true,  2 },
  {"",       0, false, 0 }
};


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

#ifdef IWAVE_USE_MPI
    if(retrieveGlobalRank()==0){
#endif

      if(argc<2){
	pagedoc();
	exit(0);
      }

      cerr << "\n////////////////////////////////////////////////////////////\n"
	   << "Running acd_scal_inv.x\n\n";

#ifdef IWAVE_USE_MPI
    }
#endif
    
    // Initializing forward map F
    MPS_KEYS mk;
    mk.MPS_file   = "MPS_file_inv";
    mk.grid_file  = "csq";     
    mk.appx_ord   = "appx_ord";
    mk.delta_file = "MPS_delta";
    mk.RHS_files.push_back("source");
    mk.hdr_files.push_back("data");
    mk.G_files.push_back("data_g");
    MPS_IWaveLOVOp<Scal_MPS_Space> F(mk,*pars,stream);
    
    //Initializing input/output vectors
    Vector<float> m(F.getProductDomain()[0]); //medium parameter vector
    Vector<float> w(F.getProductDomain()[1]); //estimated source parameter vector
    Vector<float> d(F.getRange());            //data vector
    Vector<float> r(F.getRange());            //residual vector, r=F(m)w-d

    AssignFilename af_csq(valparse<string>(*pars,"csq"));
    AssignFilename af_mps(valparse<string>(*pars,"MPS_file_inv"));
    AssignFilename af_data(valparse<string>(*pars,"data"));    
    AssignFilename af_datares(valparse<string>(*pars,"datares"));    

    m.eval(af_csq);
    w.eval(af_mps);
    d.eval(af_data);
    r.eval(af_datares);

    //Extracting linear restriction of F
    LinearRestrictOp<float> F_lin(F,m);
    
    //Initializing preconditioner op M
    float c       = valparse<float>(*pars,"cmax",1.0f);
    float order_0 = valparse<float>(*pars,"order_0",0);
    float order_d = valparse<float>(*pars,"order_d",0);
    bool precond  = (order_0!=0) || (order_d!=0);
    
    MPS_frac_cal Q(*pars,
		    F.getLinDomain(),
		    c,
		    order_0,
		    order_d );
    InvLinearOp<ireal> Q_inv(Q);
    AdjLinearOp<ireal> Q_inv_adj(Q_inv);
    NormalLinearOp<ireal> M(Q_inv_adj);

    // Initializing CG algorithm
    float rtol     = valparse<float>(*pars,"CG_RTol",1000.0*numeric_limits<float>::epsilon());
    float gtol     = valparse<float>(*pars,"CG_GTol",1000.0*numeric_limits<float>::epsilon());
    float trustrad = valparse<float>(*pars,"CG_TrustRad",numeric_limits<float>::max());
    int   maxiter  = valparse<int>(*pars,"CG_MaxIter",10);

    sstream<<scientific;
    sstream<<"\n*******************************************************\n"
	   <<"* Acoustic (acd) MPS Inversion via                       \n"
	   <<"* Conjugate Gradient Algorithm for Normal Eqns           \n"
	   <<"* max iterations       = "<< maxiter   <<"\n"
	   <<"* residual tolerance   = "<< rtol      <<"\n"
	   <<"* normal res tolerance = "<< gtol      <<"\n"
	   <<"* trust radius         = "<< trustrad <<"\n"
	   <<"*******************************************************\n";
      
    float rnorm, nrnorm;
    RVLAlg::Algorithm * alg = NULL;

    if( precond ){
      RVLUmin::CGNEAlg<float> *cgalg = 
	new RVLUmin::CGNEAlg<float> 
	(w,
	 F_lin,
	 M,
	 d,
	 rnorm,
	 nrnorm,
	 rtol,
	 gtol,
	 maxiter,
	 trustrad,
	 sstream);
      
      alg=cgalg;
    }
    else{
      RVLUmin::CGNEAlg<float> *cgalg = 
	new RVLUmin::CGNEAlg<float> 
	(w,
	 F_lin,
	 d,
	 rnorm,
	 nrnorm,
	 rtol,
	 gtol,
	 maxiter,
	 trustrad,
	 sstream);
      
      alg=cgalg;
    }

    //Running CG
    alg->run();
    sstream << "Objective functional: |F[m]w-d|^2 = "<< (rnorm*rnorm) << "\n";
    
    //Recomputing residual
    F_lin.applyOp(w,r);

    
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

#ifdef IWAVE_USE_MPI
    if(retrieveGlobalRank()==0){
#endif

      cerr << "\nFinishing acd_scal_inv.x\n"
	   << "////////////////////////////////////////////////////////////\n";

#ifdef IWAVE_USE_MPI
    }
    MPI_Finalize();
#endif

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
