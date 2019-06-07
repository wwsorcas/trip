// MPS_howto.cc
// Author: Mario J. Bencomo
// last modified: 03/31/18

#define MJB_VERBOSE

#include "asg_defn.hh"
#include "asg.hh"
#include "MPS_howto.hh"
#include "MPS_to_RHS.hh"

using RVL::AssignFilename;
using RVL::valparse;
using RVL::LinCompLOVOp;
using RVL::OperatorEvaluation;

using TSOpt::add_to_pars;
using TSOpt::ExScal_MPS_Space;
using TSOpt::MPS_to_RHS;
using TSOpt::MPS_KEYS;
using TSOpt::IWaveLOVOp;


const char *sdoc[] = {
  " ============================================================================",
  " MPS_howto.x ",
  " ============================================================================",
  " Authors: Mario J. Bencomo",
  " ",
  " 2D/3D acoustic variable density centered difference modeling (asg) with an example",
  " scalar (pressure) multipole source."
  "",
  " PDE:",
  "   dp/dt + beta div(v) = w_0(t) delta(x-x^*) + w_1(t) {d/dx_0 + d/dx_1} delta(x-x^*)",
  "   dv/dt + kappa grad(p) = 0",
  "",
  " Input:",
  "   kappa = bulk modulus (units GPa)",
  "   beta = buoyancy (units cm^3/g)",
  "   {w_0,w_1} = multipole coeffiecients ",
  "               (units GPa*m^2/s for w_0, units GPa*m^3/s for w_1)",   
  "",
  " Output:",
  "   p = pressure data (units GPa)",
  " ",
  " Typical parameter list. May be copied, edited, and used for input: either",
  " include parameters on command line (for example in Flow), or place",
  " in file <foo> and include \"par=<foo>\" on command line. Any parameter",
  " included explicitly in command line overrides parameter with same key",
  " in par file.",
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
  "      source_p = <path>    RHS source file, SU format - need filename.",
  "                           Currently driver assumes file does not exist,"
  "                           thus will generate such file as part of the"
  "                           initialization of MPS_to_RHS operators.",
  "      appx_ord = <2*order> Singular source approximation order.",
  "                           Should ideally be equal to the FD spatial order,",
  "                           i.e., appx_ord=2*order",
  "      MPS_file = <path>    MPS file containing MPS coefficients, SU format - must exist.",
  "                           For this particular MPS space the MPS file must contain"
  "                           two traces, and source locations.",
  " ",
  " ------------------------------------------------------------------------",
  " Data info:",
  " ",
  "          data_p = <path>    data file, SU format - must exist ",
  "                             headers establish acquisition geometry",
  "                             output for forward modeling",
  " ",
  " ------------------------------------------------------------------------",
  " Model info:",
  " ",
  "        bulkmod = <path>      input bulk modulus, ",
  "                              determines simulation spatial grid,",
  "       buoyancy = <path>      input buoyancy",
  " ",
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
    if( retrieveGlobalRank()==0 ){
#endif

      if(argc<2){
	pagedoc();
	exit(0);
      }

      cerr << "\n============================================================\n" 
	   << "Running MPS_howto.x\n";

      //ps_printall(*pars,stderr);
#ifdef IWAVE_USE_MPI
    }
#endif

    //Initializing MPS space
    MPS_KEYS mk_sp;
    mk_sp.MPS_file  = "MPS_file"; 
    mk_sp.grid_file = "bulkmod";  
    ExScal_MPS_Space MPS_SP(mk_sp,*pars);

#ifdef MJB_VERBOSE
    str << "------------------------------\n"
	<< "Printing out info of MPS space\n"
	<< "------------------------------\n";
    MPS_SP.write(str);
#endif

    //Initializing MPS_to_RHS op
    MPS_KEYS mk_m2r;
    mk_m2r.MPS_file  = "MPS_file";
    mk_m2r.grid_file = "bulkmod";
    mk_m2r.appx_ord  = "appx_ord";
    mk_m2r.RHS_files.push_back("source_p"); 
    MPS_to_RHS M2R(mk_m2r,*pars,MPS_SP);

#ifdef MJB_VERBOSE
    str << "----------------------------------\n"
	<< "Printing out info of MPS_to_RHS op\n"
	<< "----------------------------------\n";
    M2R.write(str);
#endif

    //setting up MPS and RHS vectors
    Vector<float> mps(M2R.getDomain());
    AssignFilename af_mps(valparse<string>(*pars,"MPS_file"));
    mps.eval(af_mps);
    Vector<float> rhs(M2R.getRange());
    AssignFilename af_rhs(valparse<string>(*pars,"source_p"));
    rhs.eval(af_rhs);
    
    //applying MPS-to-RHS op
    M2R.applyOp(mps,rhs);


    //operator composition to build forward map F
    IWaveLOVOp IWOP(M2R.get_pars(),stream);
    LinCompLOVOp<float> F(M2R,IWOP);

#ifdef MJB_VERBOSE
    str << "-------------------------\n"
	<< "Printing out forward map \n"
	<< "-------------------------\n";
    F.write(str);
#endif

    //setting up domain vector
    Vector<float> x(F.getDomain());
    Components<float> cx(x);
    Components<float> m(cx[0]); //medium params
    AssignFilename af_bmod(valparse<string>(*pars,"bulkmod"));
    AssignFilename af_buoy(valparse<string>(*pars,"buoyancy"));
    m[0].eval(af_bmod);
    m[1].eval(af_buoy);
    cx[1].eval(af_mps);

    //setting up range vector
    Vector<float> y(F.getRange());
    AssignFilename af_data(valparse<string>(*pars,"data_p"));
    y.eval(af_data);
    
    //applying op
    OperatorEvaluation<float> opeval(F,x);
    y.copy(opeval.getValue());
    
    //writing out
    string outfile = valparse<string>(*pars,"outfile","");

#ifdef IWAVE_USE_MPI
    for( int r=0; r<retrieveGlobalSize(); r++ ){
      if( retrieveGlobalRank()==r ){
	cerr << "rank="<<r<<"\n";
#endif
	if(outfile.size()>0){
	  ofstream outf;
	  outf.open(outfile.c_str(),ios::app);
	  outf << str.str();
	  outf.close();
	}
	else{
	  cerr << str.str();
	}
#ifdef IWAVE_USE_MPI
      }
      MPI_Barrier(retrieveGlobalComm());
    }
#endif

    //clean up
    ps_delete(&pars);
    iwave_fdestroy();

#ifdef IWAVE_USE_MPI    
    if( retrieveGlobalRank()==0 ){
#endif            
      cerr << "\nExiting MSP_howto.x\n"
	   << "============================================================\n\n";
#ifdef IWAVE_USE_MPI    
    }
    MPI_Finalize();
#endif

  }
  catch(RVLException &e){
    e << "Exiting with error!\n"
      << "Dumping output stream:\n"
      << str.str();

#ifdef IWAVE_USE_MPI
    if( retrieveGlobalRank()==0 ){
#endif
      e.write(cerr);
#ifdef IWAVE_USE_MPI
    }   
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
}
