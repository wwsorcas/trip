// asg_canscal_p_inv.cc
// Author: Mario J. Bencomo
// last modified: 04/07/18

#include "asg_defn.hh"
#include "asg.hh"
#include "cgnealg_mod.hh"
#include "MPS_iwop.hh"
#include "MPS_Space_Examples.hh"
#include "MPS_frac_cal.hh"

using RVL::valparse;
using RVL::RVLException;
using RVL::Vector;
using RVL::Components;
using RVL::AssignFilename;
using RVL::NormalLinearOp;
using RVL::AdjLinearOp;
using RVL::InvLinearOp;
using RVL::LinearRestrictOp;

using TSOpt::MPS_IWaveLOVOp;
using TSOpt::CanScal_MPS_Space;
using TSOpt::MPS_frac_cal;
using TSOpt::MPS_KEYS;

const char *sdoc[] = {
  " ============================================================================",
  " asg_canscal_p_inv.x ",
  " ============================================================================",
  " Authors: Mario J. Bencomo",
  " ",
  " Source inversion driver with acoustic variable density centered difference ",
  " modelling for canonical scalar multipole sources and pressure data output,",
  " see asg_p_canscal.cc for more info on the forward map.",
  "", 
  " Conjugate gradient is used to solve normal equations related to source estimation. ",
  " Fractional calculus operators are available for preconditioning.",
  " ",
  " Forward map F(m), ",
  "       F(m)w = d",
  "       m = model parameters",
  "       w = MPS coefficients",
  "       d = data (pressure)",
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
  " \"asg_p_canscal_inv.x [parameters] [standalone install]\"",
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
  "           dmin = 1.0         min permitted density (g/cm^3) - sanity check",
  "           dmax = 1.0         max permitted density (g/cm^3) - sanity check",
  " ",
  " ------------------------------------------------------------------------",
  " Source info:",
  " ",
  "         source_p=empty       dummy key required for asg driver.",
  " ",
  " ------------------------------------------------------------------------",
  " Trace info:",
  " ",
  "          data_p  = <path>    data file, SU format - must exist ",
  "                              headers establish acquisition geometry",
  "                              output for forward modeling (including all",
  "                              orders of derivative), input for adjoint",
  " ",
  " ------------------------------------------------------------------------",
  " Model info:",
  " ",
  "        bulkmod = <path>      input bulk modulus, ",
  "                              determines simulation spatial grid,",
  "        buoyancy= <path>      input buoyancy",
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
  "      appx_ord = 2        singular source approximation order ",
  "                          should match spatial order of FD method, i.e., 2*order",
  "       MPS_ord = 0        MPS order ",
  "      MPS_file = <path>   MPS file containing MPS coefficients, SU format - must exist",
  "                          Must contain source locations.",
  "     MPS_delta = <path>   MPS file containing delta traces, SU format - must exist",
  "      data_p_g = <path>   file for storing Green's functions, SU format - must exist",
  " ------------------------------------------------------------------------",
  " CG info:",
  " ",
  "          Maxiter = 10             Max number of CG iterations",
  "      ResidualTol = 1000*eps_mach  Residual tolerance",
  "      GradientTol = 1000*eps_mach  Gradient tolerance",
  "         TrustRad = infty          Trust region radius",
  "          outfile = cerr           Text file to dump CG iterate info",
  " ",
  " Preconditioner info:",
  "          order_0 = 0  preconditioner base order",
  "          order_d = 0  preconditioner incerement per MPS order",
  "",
  " Stacking over sources:",
  "          stack = 0    flag for stacking inverted sources, 1=true, 0=false",
  " ---------------------------end parameters ------------------------------",
  NULL };


IOKEY IWaveInfo::iwave_iokeys[]
= {
  { "bulkmod",   0, true,  1 },
  { "buoyancy",  1, true,  1 },
  { "source_p",  2, true,  2 },
  { "data_p",    2, false, 2 },
  {"",           0, false, 0 }
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
    TSOpt::IWaveEnvironment(argc, argv, 0, &pars, &stream);

#ifdef IWAVE_USE_MPI
    if(retrieveGlobalRank()==0){
#endif

      if(argc<2){
	pagedoc();
	exit(0);
      }

      cerr << "\n////////////////////////////////////////////////////////////\n"
	   << "Running asg_canscal_p_inv.x\n";

#ifdef IWAVE_USE_MPI
    }
#endif


    //Initializing forward map F
    TSOpt::MPS_KEYS mks;
    mks.appx_ord   = "appx_ord";
    mks.MPS_file   = "MPS_file";
    mks.delta_file = "MPS_delta";
    mks.grid_file  = "bulkmod";
    mks.RHS_files.push_back("source_p");
    mks.hdr_files.push_back("data_p");
    mks.G_files.push_back("data_p_g");
    MPS_IWaveLOVOp<CanScal_MPS_Space> F(mks,*pars,stream);

    //Setting up input/output vectors
    Vector<ireal> m(F.getProductDomain()[0]); //medium parameters
    Components<ireal> cm(m);
    Vector<ireal> w(F.getProductDomain()[1]); //estimated source parameters
    Vector<ireal> d(F.getRange()); //data vector
    
    AssignFilename af_bmd(valparse<string>(*pars,"bulkmod"));
    AssignFilename af_buo(valparse<string>(*pars,"buoyancy"));
    AssignFilename af_mps(valparse<string>(*pars,"MPS_file"));
    AssignFilename af_data(valparse<string>(*pars,"data_p"));

    cm[0].eval(af_bmd);
    cm[1].eval(af_buo);
    w.eval(af_mps);
    d.eval(af_data);
    
    Vector<ireal> w0(F.getProductDomain()[1]);
    AssignFilename af_mps_true(valparse<string>(*pars,"MPS_file_true",
						valparse<string>(*pars,"MPS_file")));
    w0.eval(af_mps_true);
    
    //Extracting linear restriction of F
    LinearRestrictOp<float> F_lin(F,m);

    //Initializing preconditioner op M
    float c       = valparse<float>(*pars,"PC_c",1.0f);
    float order_0 = valparse<float>(*pars,"PC_ord0",0);
    float order_d = valparse<float>(*pars,"PC_ordd",0);
    bool  precond = (order_0!=0) || (order_d!=0);

    MPS_frac_cal Q(*pars,
		   F.getLinDomain(),
		   c,
		   order_0,
		   order_d );
    
    //Initializing CG algorithm
    float rtol     = valparse<float>(*pars,"ResidualTol",1000.0*numeric_limits<float>::epsilon());
    float gtol     = valparse<float>(*pars,"GradientTol",1000.0*numeric_limits<float>::epsilon());
    float xtol     = valparse<float>(*pars,"xTol",-numeric_limits<float>::max());    
    float Lxtol    = valparse<float>(*pars,"LxTol",-numeric_limits<float>::max());    
    float trustrad = valparse<float>(*pars,"TrustRad",numeric_limits<float>::max());
    int   maxiter  = valparse<int>  (*pars,"MaxIter",10);
    float mu       = valparse<float>(*pars,"Tikhonov_reg",0.0);

    if( mu<0 ){
      RVLException e;
      e << "Error, Tikhonov regularization parameter must be >=0!\n"
	<< "  mu = "<<mu<<"\n";
      throw e;
    }

    sstream <<scientific;
    sstream <<"\n*******************************************************\n"
	    <<"* Acoustic (asg) MPS Inversion via \n"
	    <<"* Conjugate Gradient Algorithm for Normal Eqns \n"
	    <<"* max iterations        = "<< maxiter  <<"\n"
	    <<"* residual tolerance    = "<< rtol     <<"\n"
	    <<"* normal res tolerance  = "<< gtol    <<"\n"
	    <<"* solution max size     = "<< xtol    <<"\n"
	    <<"* weighted sol max size = "<< Lxtol    <<"\n"
	    <<"* trust radius          = "<< trustrad <<"\n"
	    <<"*******************************************************\n";

    float rnorm, nrnorm, xnorm, Lxnorm, xerrs, Lxerrs;
    RVLAlg::Algorithm * alg = NULL;


    //with preconditioning
    if( precond ){
      //with regularization
      if( mu!=0.0 ){
	RVLUmin::CGNEAlg_mod<float> *cgalg = 
	  new RVLUmin::CGNEAlg_mod<float> 
	  (w,
	   F_lin,
	   Q,
	   d,
	   w0,
	   mu,
	   rnorm,
	   nrnorm,
	   xnorm,
	   Lxnorm,
	   xerrs,
	   Lxerrs,
	   rtol,
	   gtol,
	   xtol,
	   Lxtol,
	   maxiter,
	   trustrad,
	   sstream);
      
	alg=cgalg;
      }
      //no regularization
      else{
	RVLUmin::CGNEAlg_mod<float> *cgalg = 
	  new RVLUmin::CGNEAlg_mod<float> 
	  (w,
	   F_lin,
	   Q,
	   d,
	   w0,
	   rnorm,
	   nrnorm,
	   xnorm,
	   Lxnorm,
	   xerrs,
	   Lxerrs,
	   rtol,
	   gtol,
	   xtol,
	   Lxtol,
	   maxiter,
	   trustrad,
	   sstream);
      
	alg=cgalg;
      }
    }
    //no preconditioning
    else{
      //with regularization
      if( mu!=0 ){
	RVLUmin::CGNEAlg_mod<float> *cgalg = 
	  new RVLUmin::CGNEAlg_mod<float> 
	  (w,
	   F_lin,
	   d,
	   w0,
	   mu,
	   rnorm,
	   nrnorm,
	   xnorm,
	   xerrs,
	   rtol,
	   gtol,
	   xtol,
	   maxiter,
	   trustrad,
	   sstream);
      
	alg=cgalg;
      }
      //no regularization
      else{
	RVLUmin::CGNEAlg_mod<float> *cgalg = 
	  new RVLUmin::CGNEAlg_mod<float> 
	  (w,
	   F_lin,
	   d,
	   w0,
	   rnorm,
	   nrnorm,
	   xnorm,
	   xerrs,
	   rtol,
	   gtol,
	   xtol,
	   maxiter,
	   trustrad,
	   sstream);
      
	alg=cgalg;
      }
    }

    //Running CG
    alg->run();
    sstream << "Final value of objective functional: |F(m)w-d|^2 = "<<(rnorm*rnorm)<<"\n";

    //Recomputing residual
    //F_lin.applyOp(w,r);


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

      cerr << "\nFinishing asg_canscal_p_inv.x\n"
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
