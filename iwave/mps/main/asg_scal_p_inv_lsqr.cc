// asg_scal_p_inv_lsqr.cc
// Author: Mario J. Bencomo
// last modified: 04/01/18

#include "asg_defn.hh"
#include "asg.hh"
#include "functions.hh"
#include "lsqr.hh"
#include "VPM.hh"
#include "VPM_mod.hh"
#include "MPS_iwop.hh"
#include "MPS_Space_Examples.hh"
#include "MPS_frac_cal.hh"

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

using RVLUmin::LSQRPolicy;
using RVLUmin::LSQRPolicyData;
using RVLUmin::VPM;
using RVLUmin::VPM_mod;

using TSOpt::IWaveEnvironment;
using TSOpt::MPS_IWaveLOVOp;
using TSOpt::CanScal_MPS_Space;
using TSOpt::CanVec_MPS_Space;
using TSOpt::ExVec_MPS_Space;
using TSOpt::Scal_MPS_Space;
using TSOpt::MPS_frac_cal;

const char *sdoc[] = {
  " ============================================================================",
  " asg_scal_p_inv_lsqr.x ",
  " ============================================================================",
  " Authors: Mario J. Bencomo",
  " ",
  " Source inversion driver with acoustic variable density centered difference ",
  " modelling for single scalar multipole sources and pressure data output,", 

  " LSQR is used to solve normal equations related to source estimation. ",
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
  " \"asg_p_scal_inv.x [parameters] [standalone install]\"",
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
  "         source_p=empty       key required related to RHS sources for asg driver,",
  "                              set to empty to generate filename from code",
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
  "      appx_ord = 4        singular source approximation order ",
  "     MPS_ord_0 = 0        MPS order in z-direction ",
  "     MPS_ord_1 = 0        MPS order in x-direction ",
  "     MPS_ord_2 = 0        MPS order in y-direction ",
  "      MPS_file = <path>   MPS file containing MPS coefficients, SU format - must exist",
  "                          Must contain source locations.",
  "     MPS_delta = <path>   MPS delta file, SU format - must exist",
  "      data_p_g = <path>   Green's function file or 'Green data', SU format - must exist",
  " ------------------------------------------------------------------------",
  " LSQR info:",
  " ",
  "          Maxiter = (optional)  Max number of CG iterations",
  "      ResidualTol = (optional)  Residual tolerance",
  "      GradientTol = (optional)  Gradient tolerance",
  "          MaxStep = (optional)  Trust region radius",
  "          outfile = (optional)  Text file to dump CG iterate info",
  " ",
  " Preconditioner info:",
  "          order_0 = (optional)  preconditioner base order",
  "          order_d = (optional)  preconditioner incerement per MPS order",
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

  std::stringstream str;

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
	   << "Running asg_scal_p_inv_lsqr.x\n";

#ifdef IWAVE_USE_MPI
    }
#endif

    TSOpt::MPS_KEYS mks;
    mks.appx_ord   = "appx_ord";
    mks.MPS_file   = "MPS_file";
    mks.delta_file = "MPS_delta";
    mks.grid_file  = "bulkmod";
    mks.RHS_files.push_back("source_p");
    mks.hdr_files.push_back("data_p");
    mks.G_files.push_back("data_p_g");
    

    // Initializing forward map F
    MPS_IWaveLOVOp<Scal_MPS_Space> F(mks,*pars,stream);

    Vector<ireal> m(F.getProductDomain()[0]);
    Vector<ireal> w(F.getProductDomain()[1]);
    Components<ireal> cm(m);
    
    AssignFilename af_bmd(valparse<string>(*pars,"bulkmod"));
    AssignFilename af_buo(valparse<string>(*pars,"buoyancy"));
    AssignFilename af_mps(valparse<string>(*pars,"MPS_file"));

    cm[0].eval(af_bmd);
    cm[1].eval(af_buo);
    w.eval(af_mps);


    // Initializing preconditioner op M
    float c = valparse<float>(*pars,"cmax",1.0f);
    float order_0 = valparse<float>(*pars,"order_0",0);
    float order_d = valparse<float>(*pars,"order_d",0);
    bool precond = (order_0!=0) || (order_d!=0);
    
    /*
    MPS_frac_cal Q_inv(*pars,
		       F.getLinDomain(),
		       c,
		       order_0,
		       order_d );
    
    AdjLinearOp<ireal> Q_inv_adj(Q_inv);
    NormalLinearOp<ireal> M(Q_inv_adj);
    */
    //cerr << "---> Preconditioner:\n ";
    //M.write(cerr);

    // Initializing data y
    Vector<ireal> y(F.getRange());
    AssignFilename af_data(valparse<string>(*pars,"data_p"));


    // Initializing CG policy
    float rtol=valparse<float>(*pars,"ResidualTol",100.0*numeric_limits<float>::epsilon());
    float nrtol=valparse<float>(*pars,"GradientTol",100.0*numeric_limits<float>::epsilon());
    int maxcount=valparse<int>(*pars,"MaxIter",10);
    float maxstep=valparse<float>(*pars,"MaxStep",numeric_limits<float>::max());
    std::stringstream res;
    res<<scientific;
    
    res<<endl<<"*******************************************************"<<endl;
    res<<"* Acoustic (asg) MPS Inversion via";
    res<<"* LSQR Algorithm for Normal Eqns"<<endl;
    res<<"* max iterations       = "<<maxcount<<endl;
    res<<"* residual tolerance   = "<<rtol<<endl;
    res<<"* normal res tolerance = "<<nrtol<<endl;
    res<<"* trust radius         = "<<maxstep<<endl;
    res<<"*******************************************************"<<endl;

    LSQRPolicyData<ireal> lsqr(rtol,
			       nrtol,
			       maxstep,
			       maxcount,
			       true);


    VPM< ireal,LSQRPolicy<ireal>,LSQRPolicyData<ireal> > F_red(F,y,lsqr,res);
    FunctionalEvaluation<ireal> feval(F_red,m);
    float val = feval.getValue();
    res << "VPM value = 0.5*|F[m]w-d|^2 = "<< val << "\n";

    VPM< ireal,LSQRPolicy<ireal>,LSQRPolicyData<ireal> > const &fcp =
      dynamic_cast< VPM< ireal,LSQRPolicy<ireal>,LSQRPolicyData<ireal> > const &>(feval.getFunctional());
    w.copy(fcp.getLSSoln());     

    /*
    // Initializing and running VPM
    if( precond ){

      VPM_mod< ireal,CGNEPolicy<ireal>,CGNEPolicyData<ireal> > vpmfun_pc(lovop,M,y,cgnepd,res);
      FunctionalEvaluation<ireal> feval_pc(vpmfun_pc,cx[0]);
      float val_pc = feval_pc.getValue();
      res << "VPM with PC value = 0.5*||F[x]-d||^2 = "<< val_pc << "\n";

      VPM_mod< ireal,CGNEPolicy<ireal>,CGNEPolicyData<ireal> > const &fcp_pc =
	dynamic_cast< VPM_mod< ireal,CGNEPolicy<ireal>,CGNEPolicyData<ireal> > const &>(feval_pc.getFunctional());
      cx[1].copy(fcp_pc.getLSSoln());     

    }
    else{

      VPM< ireal,CGNEPolicy<ireal>,CGNEPolicyData<ireal> > vpmfun(lovop,y,cgnepd,res);
      FunctionalEvaluation<ireal> feval(vpmfun,cx[0]);
      float val    = feval.getValue();
      res << "VPM value = 0.5*||F[x]-d||^2 = "<< val << "\n";
    

      VPM< ireal,CGNEPolicy<ireal>,CGNEPolicyData<ireal> > const &fcp =
	dynamic_cast< VPM< ireal,CGNEPolicy<ireal>,CGNEPolicyData<ireal> > const &>(feval.getFunctional());
      cx[1].copy(fcp.getLSSoln());
    }
    */

    //writing out CG info
    lsqr.write(res);
    str << res.str();
    
#ifdef IWAVE_USE_MPI
    if (retrieveGlobalRank()==0) {
#endif
      string outfile = valparse<string>(*pars,"outfile","");
      if (outfile.size()>0) {
	ofstream outf(outfile.c_str());
	outf<<str.str();
	outf.close();
      }
      else {
	cout<<str.str();
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

      cerr << "\nFinishing asg_scal_p_inv_lsqr.x\n"
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
      e << str.str();
      e.write(cerr); 
#ifdef IWAVE_USE_MPI
    }
    MPI_Barrier(retrieveGlobalComm());
    MPI_Abort(retrieveGlobalComm(),0);
#endif 
    exit(1);
  }
}
