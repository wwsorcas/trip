// acd_scal_conv.cc
// Author: Mario J. Bencomo
// last modified: 03/31/18

#include "acd_defn.hh"
#include "acd.hh"
#include "MPS_includes.hh"
#include "MPS_Space_Examples.hh"
#include "MPS_to_RHS.hh"
#include "MPS_iwop.hh"

using RVL::valparse;
using RVL::AssignFilename;
using RVL::OperatorEvaluation;

using TSOpt::Scal_MPS_Space;
using TSOpt::MPS_to_RHS;
using TSOpt::add_to_pars;
using TSOpt::MPS_IWaveLOVOp;
using TSOpt::MPS_KEYS;

const char *sdoc[] = {
  " ============================================================================",
  " acd_scal_conv.x ",
  " ============================================================================",
  " Authors: Mario J. Bencomo",
  " ",
  " Acoustic constant density centered difference modeling for single scalar",
  " multipole sources.",
  " Forward map from source-medium parameters to data is implemented via the",
  " MPS_IWaveLOVOp operator, which computes outputs by convolution of input",
  " MPS waveforms with corresponding precomputed Green's functions (computed",
  " for a given medium parameter).",
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
  " \"acd_scal_conv.x [parameters] [standalone install]\"",
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
  "         source=empty         dummy key required for asg driver.",
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
  "        appx_ord = 4       singular source approximation order,",
  "                           should match spatial order of FD method, i.e., 2*order",
  "       MPS_ord_0 = 0       MPS order in z-direction, s_0",
  "       MPS_ord_1 = 0       MPS order in x-direction, s_1",
  "       MPS_ord_2 = 0       MPS order in y-direction, s_2",
  "        MPS_file = <path>  SU file containing MPS coefficients and source location.",
  "       MPS_delta = <path>  SU file containing filtered delta wavelet for computing Green's functions.",
  "          data_g = <path>  SU file for formating computed Green's functions.",
  "",
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
    ps_printall(*pars,stderr);

#ifdef IWAVE_USE_MPI
    if(retrieveGlobalRank()==0){
#endif

      if(argc<2){
	pagedoc();
	exit(0);
      }

      cerr << "\n////////////////////////////////////////////////////////////\n"
	   << "Running acd_scal_conv.x\n";

#ifdef IWAVE_USE_MPI
    }
#endif

    MPS_KEYS mks;
    mks.appx_ord = "appx_ord";
    mks.MPS_file = "MPS_file";
    mks.delta_file = "MPS_delta";
    mks.grid_file = "csq";
    mks.RHS_files.push_back("source");
    mks.hdr_files.push_back("data");
    mks.G_files.push_back("data_g");

    //initializing forward map F
    MPS_IWaveLOVOp<Scal_MPS_Space> F(mks,*pars,stream);

    //setting up domain vector
    Vector<ireal> x(F.getDomain());
    Components<ireal> cx(x);

    AssignFilename af_csq(valparse<string>(*pars,"csq"));
    cx[0].eval(af_csq);
    AssignFilename af_mps(valparse<string>(*pars,"MPS_file"));
    cx[1].eval(af_mps);

    //setting up range vector
    Vector<ireal> y(F.getRange());
    AssignFilename af_data(valparse<string>(*pars,"data"));
    y.eval(af_data);
    
    //applying op
    OperatorEvaluation<ireal> eval(F,x);
    y.copy(eval.getValue());

    //clean up
    ps_delete(&pars);
    iwave_fdestroy();


#ifdef IWAVE_USE_MPI
    if(retrieveGlobalRank()==0){
#endif

      cerr << "\nFinishing acd_scal_conv.x\n"
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
