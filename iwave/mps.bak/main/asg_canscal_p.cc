// asg_canscal_p.cc
// Author: Mario J. Bencomo
// last modified: 04/07/18

#include "asg_defn.hh"
#include "asg.hh"
#include "MPS_includes.hh"
#include "MPS_Space_Examples.hh"
#include "MPS_to_RHS.hh"

using RVL::valparse;
using RVL::AssignFilename;
using RVL::OperatorEvaluation;
using TSOpt::CanScal_MPS_Space;
using TSOpt::MPS_to_RHS;
using TSOpt::add_to_pars;

const char *sdoc[] = {
  " ============================================================================",
  " asg_canscal_p.x ",
  " ============================================================================",
  " Authors: Mario J. Bencomo",
  " ",
  " Acoustic variable density centered difference modeling for canonical scalar",
  " multipole sources and pressure data output.", 
  "",
  " PDE:",
  "   dp/dt + bulk{ div v } = f(t,x)",
  "   dv/dt + buoy grad p = 0",
  "",
  " Source term:",
  "   MPS_order=0,",
  "      f(x,t) = w_0(t) delta(x-x^*)",
  "   MPS_order=1, 2-D,",
  "      f(x,t) = w_0(t) delta(x-x^*)",
  "             + w_1(t) d/dx_0 delta(x-x^*)",
  "             + w_2(t) d/dx_1 delta(x-x^*)",
  "   MPS_order=1, 3-D,",
  "      f(x,t) = w_0(t) delta(x-x^*)",
  "             + w_1(t) d/dx_0 delta(x-x^*)",
  "             + w_2(t) d/dx_1 delta(x-x^*)",
  "             + w_3(t) d/dx_2 delta(x-x^*)",
  " ",
  " Typical parameter list. May be copied, edited, and used for input: either",
  " include parameters on command line (for example in Flow), or place",
  " in file <foo> and include \"par=<foo>\" on command line. Any parameter",
  " included explicitly in command line overrides parameter with same key",
  " in par file.",
  " ",
  "  Invoke single threaded execution by ",
  " \"asg_p_canscal.x [parameters] [standalone install]\"",
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
  "      appx_ord = 4        singular source approximation order,",
  "                          should match spatial order of FD method, i.e., 2*order",
  "       MPS_ord = 0        MPS order ",
  "      MPS_file = <path>   MPS file containing MPS coefficients, SU format - must exist",
  "                          Must contain source locations.",
  " ---------------------------end parameters ------------------------------",
  NULL };

IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"bulkmod",  0, true,  1 },
  {"buoyancy", 1, true,  1 },
  {"source_p", 2, true,  2 },
  {"data_p",   2, false, 2 },
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
	   << "Running asg_canscal_p.x\n";

#ifdef IWAVE_USE_MPI
    }
#endif

    TSOpt::MPS_KEYS mks;
    mks.appx_ord = "appx_ord";
    mks.MPS_file = "MPS_file";
    mks.grid_file = "bulkmod";
    mks.RHS_files.push_back("source_p");

    //initializing MPS space
    CanScal_MPS_Space mps_sp(mks,*pars);
    
    //initializing MPS_to_RHS
    MPS_to_RHS m2r(mks,*pars,mps_sp);
    
    //constructing forward map
    TSOpt::IWaveLOVOp lovop(m2r.get_pars(),stream);
    RVL::LinCompLOVOp<float> F(m2r,lovop);

    //setting up domain vector
    Vector<ireal> x(F.getDomain());
    Components<ireal> cx(x);

    Components<ireal> cx_m(cx[0]); //model params
    RVL::AssignFilename af_bmd(valparse<string>(*pars,"bulkmod"));
    cx_m[0].eval(af_bmd);
    RVL::AssignFilename af_buo(valparse<string>(*pars,"buoyancy"));
    cx_m[1].eval(af_buo);
    Components<ireal> cx_w(cx[1]); //source params
    RVL::AssignFilename af_mps(valparse<string>(*pars,"MPS_file"));
    cx_w[0].eval(af_mps);

    //setting up range vector
    Vector<ireal> y(F.getRange());
    RVL::AssignFilename af_data(RVL::valparse<string>(*pars,"data_p"));
    y.eval(af_data);
    
    //applying op
    RVL::OperatorEvaluation<ireal> eval(F,x);
    y.copy(eval.getValue());

    //clean up
    ps_delete(&pars);
    iwave_fdestroy();


#ifdef IWAVE_USE_MPI
    if(retrieveGlobalRank()==0){
#endif

      cerr << "\nFinishing acd_canscal_p.x\n"
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
