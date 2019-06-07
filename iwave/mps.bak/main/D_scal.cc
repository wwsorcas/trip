// D_scal.cc
// Author: Mario J. Bencomo
// last modified: 11/10/16

#include "acd_defn.hh"
#include "acd.hh"
#include "MPS_includes.hh"
#include "MPS_Space_Examples.hh"
#include "MPS_frac_cal.hh"

using RVL::valparse;
using RVL::AssignFilename;
using RVL::OperatorEvaluation;

using TSOpt::Scal_MPS_Space;
using TSOpt::add_to_pars;
using TSOpt::MPS_frac_cal;

const char *sdoc[] = {
  " ============================================================================",
  " D_scal.x ",
  " ============================================================================",
  " Authors: Mario J. Bencomo",
  " ",
  " Applying fractional calculus operator to MPS coefficients of a scalar MPS,",
  " of the form:",
  "   2-D, f(x,t) = f_0(t) (d/dx_0)^{s_0} (d/dx_1)^{s_1} delta(x-x^*)",
  "   3-D, f(x,t) = f_0(t) (d/dx_0)^{s_0} (d/dx_1)^{s_1} (d/dx_2)^{s_2} delta(x-x^*)",
  " where {s_0,s_1,s_2} are real valued.",
  " Note, negative derivative powers s_i refer to fractional integral operators.",
  "",
  " The fractional calculus operator applied to an MPS with multi-index s={s_0,s_1,s_2}",
  " is given as follows:",
  "   c^{-q}(d/dt)^{q_0 + q} f_0(t) ",
  " where q_0 is the starting order, q_0=order_0, and q is the incremental order, ",
  " q=order_d*max{s_0,s_1,s_2}.",
  " ",
  " Typical parameter list. May be copied, edited, and used for input: either",
  " include parameters on command line (for example in Flow), or place",
  " in file <foo> and include \"par=<foo>\" on command line. Any parameter",
  " included explicitly in command line overrides parameter with same key",
  " in par file.",
  " ",
  "  Invoke single threaded execution by ",
  " \"D_scal.x [parameters] [standalone install]\"",
  " ",
  " --------------------------- begin parameters ---------------------------",
  "      order_0 = 0.0    Starting order of frac calculus op.",
  "      order_d = 1.0    Incremental of order of frac calculus with MPS order.",
  "    MPS_ord_0 = 0.0    MPS order in z-direction, s_0",
  "    MPS_ord_1 = 0.0    MPS order in x-direction, s_1",
  "    MPS_ord_2 = 0.0    MPS order in y-direction, s_2",
  "        c_max = 3.0    Maximum velocity in km/s, used as scaling factor in op",
  "       D_file = <path> Output MPS coefficients SU file - must exist",
  "     MPS_file = <path> Input MPS coefficients SU file - must exist",
  "         grid = <path> RSF file determining spatial grid info - must exist",
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
	   << "Running D_acd_scal.x\n";
      //ps_printall(*pars,stderr);
#ifdef IWAVE_USE_MPI
    }
#endif


    TSOpt::MPS_KEYS mks;
    mks.appx_ord = "appx_ord";
    mks.MPS_file = "MPS_file";
    mks.grid_file = "grid";    

    float c = valparse<float>(*pars,"c_max");

    Scal_MPS_Space mps_sp(mks,*pars);

    MPS_frac_cal D( *pars, 
		    mps_sp,
		    c,
		    valparse<float>(*pars,"order_0"),
		    valparse<float>(*pars,"order_d") );

    Vector<ireal> x(mps_sp);
    AssignFilename af_x(valparse<string>(*pars,"MPS_file"));
    x.eval(af_x);

    Vector<ireal> y(mps_sp);
    AssignFilename af_y(valparse<string>(*pars,"D_file"));
    y.eval(af_y);
    
    OperatorEvaluation<ireal> eval_D(D,x);
    y.copy(eval_D.getValue());
    

    //clean up
    ps_delete(&pars);
    iwave_fdestroy();


#ifdef IWAVE_USE_MPI
    if(retrieveGlobalRank()==0){
#endif

      cerr << "\nFinishing D_scal.x\n"
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
