// acd_scal.cc
// Author: Mario J. Bencomo
// last modified: 03/31/18

#include "acd_defn.hh"
#include "acd.hh"
#include "MPS_includes.hh"
#include "MPS_Space_Examples.hh"
#include "MPS_to_RHS.hh"

//#define VERBOSE_MJB


using RVL::valparse;
using RVL::AssignFilename;
using RVL::OperatorEvaluation;
using RVL::LinCompLOVOp;

using TSOpt::Scal_MPS_Space;
using TSOpt::MPS_to_RHS;
using TSOpt::add_to_pars;
using TSOpt::MPS_KEYS;
using TSOpt::IWaveLOVOp;


const char *sdoc[] = {
  " ============================================================================",
  " acd_scal.x ",
  " ============================================================================",
  " Authors: Mario J. Bencomo",
  " ",
  " Acoustic constant density centered difference modeling for single scalar",
  " multipole sources.",
  " Forward map from source-medium parameters to data is implemented by operator",
  " composition of MPS_to_RHS operator with IWave LOVOp."
  " ",
  " PDE:",
  "   (d/dt)^2 p - c^2 nabla^2 p = f(x,t)",
  " Source term:",
  "   2-D, f(x,t) = w(t) (d/dx_0)^{s_0} (d/dx_1)^{s_1} delta(x-x^*)",
  "   3-D, f(x,t) = w(t) (d/dx_0)^{s_0} (d/dx_1)^{s_1} (d/dx_2)^{s_2} delta(x-x^*)",
  " ",
  " Typical parameter list. May be copied, edited, and used for input: either",
  " include parameters on command line (for example in Flow), or place",
  " in file <foo> and include \"par=<foo>\" on command line. Any parameter",
  " included explicitly in command line overrides parameter with same key",
  " in par file.",
  " ",
  "  Invoke single threaded execution by ",
  " \"acd_scal.x [parameters] [standalone install]\"",
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
  "         source=empty         dummy key required for acd driver.",
  " ",
  " ------------------------------------------------------------------------",
  " Trace info:",
  " ",
  "          data  = <path>    data file, SU format - must exist ",
  "                            headers establish acquisition geometry",
  "                            output for forward modeling (including all",
  "                            orders of derivative), input for adjoint",
  " ",
  " ------------------------------------------------------------------------",
  " Model info:",
  " ",
  "        csq = <path>      input velocity squared model, ",
  "                          determines simulation spatial grid,",
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
  "        appx_ord = 4       singular source approximation order",
  "                           (should match spatial order of FD method, i.e., 2*order)",
  "       MPS_ord_0 = 0       MPS order in z-direction, s_0",
  "       MPS_ord_1 = 0       MPS order in x-direction, s_1",
  "       MPS_ord_2 = 0       MPS order in y-direction, s_2",
  "        MPS_file = <path>  SU file containing MPS coefficients and source location.",
  " ---------------------------end parameters ------------------------------",
  NULL };


IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csq",      0, true,  1 },
  {"data",     1, false, 2 },
  {"source",   1, true,  2 },
  {"",         0, false, 0 }
};


int xargc;
char **xargv;

int main(int argc, char ** argv) {

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
	   << "Running acd_scal.x\n";

#ifdef VERBOSE_MJB
      ps_printall(*pars,stderr);
#endif

#ifdef IWAVE_USE_MPI
    }
#endif

    //initializing MPS space
    MPS_KEYS mk_sp;
    mk_sp.MPS_file  = "MPS_file";
    mk_sp.grid_file = "csq";
    Scal_MPS_Space MPS_SP(mk_sp,*pars);
    
    //initializing MPS_to_RHS operator
    MPS_KEYS mk_m2r;
    mk_m2r.MPS_file  = "MPS_file";
    mk_m2r.grid_file = "csq";
    mk_m2r.appx_ord  = "appx_ord";
    mk_m2r.RHS_files.push_back("source");
    MPS_to_RHS M2R(mk_m2r,*pars,MPS_SP);
    
    //constructing forward map
    IWaveLOVOp IWOP(M2R.get_pars(),stream);
    LinCompLOVOp<float> F(M2R,IWOP);

    //setting up domain vector
    Vector<float> x(F.getDomain());
    Components<float> cx(x);
    AssignFilename af_csq(valparse<string>(*pars,"csq"));
    cx[0].eval(af_csq);
    AssignFilename af_mps(valparse<string>(*pars,"MPS_file"));
    cx[1].eval(af_mps);

    //setting up range vector
    Vector<float> d(F.getRange());
    AssignFilename af_data(valparse<string>(*pars,"data"));
    d.eval(af_data);
    
    //applying op
    OperatorEvaluation<float> eval(F,x);
    d.copy(eval.getValue());

    //clean up
    ps_delete(&pars);
    iwave_fdestroy();


#ifdef IWAVE_USE_MPI
    if(retrieveGlobalRank()==0){
#endif

      cerr << "\nFinishing acd_scal.x\n"
	   << "////////////////////////////////////////////////////////////\n";

#ifdef IWAVE_USE_MPI
    }
    MPI_Finalize();
#endif

  }
  catch(RVLException &e){
    e << "Exiting with error!\n";
#ifdef IWAVE_USE_MPI
    if( retrieveGlobalRank()==0 ){
#endif
      e.write(cerr);
#ifdef IWAVE_USE_MPI
    }
    MPI_Barrier(retrieveGlobalComm());
    MPI_Abort(retrieveGlobalComm(),0);
#endif 
    exit(1);
  }
}
