// W_res.cc
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

using TSOpt::CanScal_MPS_Space;
using TSOpt::CanVec_MPS_Space;
using TSOpt::Scal_MPS_Space;
using TSOpt::add_to_pars;
using TSOpt::MPS_frac_cal;

const char *sdoc[] = {
  " ============================================================================",
  " W_res.x ",
  " ============================================================================",
  " Authors: Mario J. Bencomo",
  " ",
  " Applying fractional calculus operator to compute normalized weighted residuals ",
  " of two multipole sources, MPS_true and MPS_est:",
  "      |W*(MPS_true-MPS_est)|/|W*MPS_true| ",
  " ",
  " Typical parameter list. May be copied, edited, and used for input: either",
  " include parameters on command line (for example in Flow), or place",
  " in file <foo> and include \"par=<foo>\" on command line. Any parameter",
  " included explicitly in command line overrides parameter with same key",
  " in par file.",
  " ",
  "  Invoke single threaded execution by ",
  " \"W_res.x [parameters] [standalone install]\"",
  " ",
  " --------------------------- begin parameters ---------------------------",
  "     appx_ord = 2      Singular source approx order",
  "      order_0 = 0.0    Starting order of frac calculus op.",
  "      order_d = 1.0    Incremental of order of frac calculus with MPS order.",
  "    MPS_ord_0 = 0.0    MPS order in z-direction",
  "    MPS_ord_1 = 0.0    MPS order in x-direction",
  "    MPS_ord_2 = 0.0    MPS order in y-direction",
  "        c_max = 3.0    Maximum velocity in km/s, used as scaling factor in op",
  "     MPS_file = <path> Input MPS SU file for building MPS space and ops.",
  "     MPS_true = <path> Input/true MPS coefficients SU file"
  "      MPS_est = <path> Input/estimate MPS coefficients SU file"
  "         grid = <path> RSF file determining spatial grid info",
  "     MPS_type = scal   Type of MPS space, options: CanScal, CanVec, Scal",
  "      outfile = <path> Output text file to dump results",
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
	   << "Running W_res.x\n";
      //ps_printall(*pars,stderr);
#ifdef IWAVE_USE_MPI
    }
#endif

    std::stringstream str;

    TSOpt::MPS_KEYS mks;
    mks.appx_ord = "appx_ord";
    mks.MPS_file = "MPS_file";
    mks.grid_file = "grid";    

    string MPS_type = valparse<string>(*pars,"MPS_type");
    string MPS_true = valparse<string>(*pars,"MPS_true");
    string MPS_est  = valparse<string>(*pars,"MPS_est");
    float c = valparse<float>(*pars,"c_max");

    if( MPS_type.compare("CanScal")==0 ){

      CanScal_MPS_Space mps_sp(mks,*pars);
      MPS_frac_cal D( *pars, 
		      mps_sp,
		      c,
		      valparse<float>(*pars,"order_0"),
		      valparse<float>(*pars,"order_d") );
      
      Vector<ireal> x_true(mps_sp);
      AssignFilename af_true(valparse<string>(*pars,"MPS_true"));
      x_true.eval(af_true);

      Vector<ireal> x_est(mps_sp);
      AssignFilename af_est(valparse<string>(*pars,"MPS_est"));
      x_est.eval(af_est);

      Vector<ireal> x_res(x_est);
      x_res.linComb(-1.0,x_true);

      float nres_norm = x_res.norm()/x_true.norm();
      
      Vector<ireal> Dx_res(mps_sp);
      Vector<ireal> Dx_true(mps_sp);

      OperatorEvaluation<ireal> eval_Dx_res(D,x_res);
      Dx_res.copy(eval_Dx_res.getValue());      
      OperatorEvaluation<ireal> eval_Dx_true(D,x_true);
      Dx_true.copy(eval_Dx_true.getValue());

      float Wnres_norm = Dx_res.norm()/Dx_true.norm();
      
      str << "\n" 
	  << "|x_est - x_true|/|x_true| = "<< nres_norm <<"\n"
	  << "|x_true| = "<< x_true.norm() <<"\n"
	  << "|x_est - x_true| = "<< x_res.norm() <<"\n";

      str << "\n" 
	  << "|W(x_est - x_true)|/|W(x_true)| = "<< Wnres_norm <<"\n"
	  << "|W(x_true)| = "<< Dx_true.norm() <<"\n"
	  << "|W(x_est - x_true)| = "<< Dx_res.norm() <<"\n";

    }
    else if( MPS_type.compare("CanVec")==0 ){

      CanVec_MPS_Space mps_sp(mks,*pars);
      MPS_frac_cal D( *pars, 
		      mps_sp,
		      c,
		      valparse<float>(*pars,"order_0"),
		      valparse<float>(*pars,"order_d") );
      
      Vector<ireal> x_true(mps_sp);
      AssignFilename af_true(valparse<string>(*pars,"MPS_true"));
      x_true.eval(af_true);

      Vector<ireal> x_est(mps_sp);
      AssignFilename af_est(valparse<string>(*pars,"MPS_est"));
      x_est.eval(af_est);

      Vector<ireal> x_res(x_est);
      x_res.linComb(-1.0,x_true);

      float nres_norm = x_res.norm()/x_true.norm();
      
      Vector<ireal> Dx_res(mps_sp);
      Vector<ireal> Dx_true(mps_sp);

      OperatorEvaluation<ireal> eval_Dx_res(D,x_res);
      Dx_res.copy(eval_Dx_res.getValue());      
      OperatorEvaluation<ireal> eval_Dx_true(D,x_true);
      Dx_true.copy(eval_Dx_true.getValue());

      float Wnres_norm = Dx_res.norm()/Dx_true.norm();
      
      str << "\n" 
	  << "|x_est - x_true|/|x_true| = "<< nres_norm <<"\n"
	  << "|x_true| = "<< x_true.norm() <<"\n"
	  << "|x_est - x_true| = "<< x_res.norm() <<"\n";

      str << "\n" 
	  << "|W(x_est - x_true)|/|W(x_true)| = "<< Wnres_norm <<"\n"
	  << "|W(x_true)| = "<< Dx_true.norm() <<"\n"
	  << "|W(x_est - x_true)| = "<< Dx_res.norm() <<"\n";

      
    }
    else if( MPS_type.compare("Scal")==0 ){
      
      Scal_MPS_Space mps_sp(mks,*pars);
      MPS_frac_cal D( *pars, 
		      mps_sp,
		      c,
		      valparse<float>(*pars,"order_0"),
		      valparse<float>(*pars,"order_d") );

      Vector<ireal> x_true(mps_sp);
      AssignFilename af_true(valparse<string>(*pars,"MPS_true"));
      x_true.eval(af_true);

      Vector<ireal> x_est(mps_sp);
      AssignFilename af_est(valparse<string>(*pars,"MPS_est"));
      x_est.eval(af_est);

      Vector<ireal> x_res(x_est);
      x_res.linComb(-1.0,x_true);

      float nres_norm = x_res.norm()/x_true.norm();
      
      Vector<ireal> Dx_res(mps_sp);
      Vector<ireal> Dx_true(mps_sp);

      OperatorEvaluation<ireal> eval_Dx_res(D,x_res);
      Dx_res.copy(eval_Dx_res.getValue());      
      OperatorEvaluation<ireal> eval_Dx_true(D,x_true);
      Dx_true.copy(eval_Dx_true.getValue());

      float Wnres_norm = Dx_res.norm()/Dx_true.norm();
      
      str << "\n" 
	  << "|x_est - x_true|/|x_true| = "<< nres_norm <<"\n"
	  << "|x_true| = "<< x_true.norm() <<"\n"
	  << "|x_est - x_true| = "<< x_res.norm() <<"\n";

      str << "\n" 
	  << "|W(x_est - x_true)|/|W(x_true)| = "<< Wnres_norm <<"\n"
	  << "|W(x_true)| = "<< Dx_true.norm() <<"\n"
	  << "|W(x_est - x_true)| = "<< Dx_res.norm() <<"\n";
    }
    


#ifdef IWAVE_USE_MPI
    if(retrieveGlobalRank()==0){
#endif

      string outfile = valparse<string>(*pars,"outfile","");
      if(outfile.size()>0){
	ofstream outf;
	outf.open(outfile, std::ofstream::out | std::ofstream::app);	
	outf<<str.str();
	outf.close();
      }
      else{
	cout<<str.str();
      }

      cerr << "\nFinishing W_res.x\n"
	   << "////////////////////////////////////////////////////////////\n";

#ifdef IWAVE_USE_MPI
    }
#endif 
    
    //clean up
    ps_delete(&pars);
    iwave_fdestroy();
    
#ifdef IWAVE_USE_MPI
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
