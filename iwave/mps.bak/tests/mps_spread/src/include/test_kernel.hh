#include "adjtest.hh"

#include "MPS_Space.hh"
#include "MPS_Space_Examples.hh"
#include "MPS_spread.hh"
#include "MPS_iwop.hh"

#define CSQ "csq"
#define BULKMOD "bulkmod"
#define BUOYANCY "buoyancy"
#define SOURCE_P "source_p"
#define SOURCE_V0 "source_v0"
#define SOURCE_V1 "source_v1"
#define DATA_P "data_p"
#define DATA_V0 "data_v0"
#define DATA_V1 "data_v1"

#ifndef MY_TOL
#define MY_TOL 1e-3
#endif

#define VERBOSE_MJB

using RVL::valparse;
using RVL::AssignFilename;
using RVL::OperatorEvaluation;

using TSOpt::CanScal_MPS_Space;
using TSOpt::CanVec_MPS_Space;
using TSOpt::ExVec_MPS_Space;
using TSOpt::Scal_MPS_Space;
using TSOpt::MPS_spread;


template< class T_MPS_Space >
void test_kernel(PARARRAY &pars, FILE &stream, string test_info){

  std::stringstream str;
  string message;

  string jobname = valparse<string>(pars,"jobname");

  try {
    str << "\n============================================================\n"
	<< "Test name: "<< jobname <<"\n"
	<< "Runnning test for MPS_spread.\n"
	<< test_info <<"\n";
#ifdef IWAVE_USE_MPI
    str << "rank = "<< retrieveGlobalRank() <<"\n";
#endif
    str << "============================================================\n";  
    
    //ps_printall(pars,stderr);

    string driver_type = valparse<string>(pars,"driver");
    string src_type = valparse<string>(pars,"src_type");
    string data_type = valparse<string>(pars,"data_type");

    if( driver_type.compare("acd")!=0 &&
        driver_type.compare("asg")!=0 ){
      RVLException e;
      e << "Driver ("<< driver_type <<") is of unexpected type!\n";
      throw e;
    }
    if( src_type.compare("pressure")!=0 &&
	src_type.compare("velocity")!=0 ){
      RVLException e;
      e << "Source ("<< src_type <<") is of unexpected type!\n";
      throw e;
    }
    if( data_type.compare("pressure")!=0 &&
	data_type.compare("velocity")!=0 ){
      RVLException e;
      e << "Data ("<< data_type <<") is of unexpected type!\n";
      throw e;
    }


    TSOpt::MPS_KEYS mks;
    mks.appx_ord   = "appx_ord";
    
    //setting mks.grid_file
    if( driver_type.compare("acd")==0 ){
      mks.grid_file = CSQ;
    }
    else{
      mks.grid_file = BULKMOD;
    }
    //setting mks.RHS_files
    if( src_type.compare("pressure")==0 ){
      mks.RHS_files.push_back(SOURCE_P);
    }
    else{
      mks.RHS_files.push_back(SOURCE_V0);
      mks.RHS_files.push_back(SOURCE_V1);
    }
    //setting mks.hdr_files and mks.G_files
    if( data_type.compare("pressure")==0 ){
      mks.hdr_files.push_back(DATA_P);
      mks.G_files.push_back("data_p_g");
    }
    else{
      mks.hdr_files.push_back(DATA_V0);
      mks.hdr_files.push_back(DATA_V1);
      mks.G_files.push_back("data_v0_g");
      mks.G_files.push_back("data_v1_g");
    }


    //////////////////////
    // Running operator //
    //////////////////////
#ifdef VERBOSE_MJB
    str << "\n[] Making data:\n";
#endif

    //initializing
    mks.MPS_file = "MPS_file";
    T_MPS_Space mps_sp(mks,pars);
    
    mks.MPS_file = "MPS_towed";
    T_MPS_Space towed_mps_sp(mks,pars);

    MPS_spread mps_spread(mps_sp,towed_mps_sp);

    //setting input and output vectors
    Vector<ireal> x(mps_spread.getMPSDomain());
    AssignFilename af_x(valparse<string>(pars,"MPS_file"));
    x.eval(af_x);

    Vector<ireal> y(mps_spread.getMPSRange());
    AssignFilename af_y(valparse<string>(pars,"MPS_towed"));
    y.eval(af_y);

    //evaluating operator
    OperatorEvaluation<ireal> opeval(mps_spread,x);
    y.copy(opeval.getValue());


    //////////////////
    // Adjoint test //
    //////////////////
#ifdef VERBOSE_MJB
    str << "\n[] Running adjoint test on MPS_spread op:\n\n";
#endif

    std::stringstream res;
    RVL::RVLRandomize<float> rnd(getpid(),-1.0,1.0);
    bool adj_test=RVL::AdjointTest<float>(mps_spread,rnd,res);
    
    if(adj_test) message="\n\033[1;32m[PASSED]\033[0m";
    else message="\n\033[1;31m[FAILED]\033[0m";
    
    message += " Adjoint test for MPS_spread op\n";
    str  << message;
    if(!adj_test) str << res.str();

#ifdef VERBOSE_MJB
    str << res.str();
#endif

#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif



    /*
    // Initializing forwad map F
    TSOpt::MPS_IWaveLOVOp<Scal_MPS_Space> lovop(mks,pars,stream);
    RVL::LinCompLOVOp<ireal> towed_lovop(mps_spread,lovop);
#ifdef VERBOSE_MJB
    cerr << "After LOVOP and MPS_spread op comp,";
    towed_lovop.write(cerr);
#endif

    Vector<ireal> X(towed_lovop.getDomain());
    Components<ireal> cX(X);
    Components<ireal> cX0(cX[0]);
    
    AssignFilename af_csq(valparse<string>(pars,"csq"));
    AssignFilename af_mps(valparse<string>(pars,"MPS_file"));
    cX0[0].eval(af_csq);
    cX[1].eval(af_mps);    

    Vector<ireal> Y(towed_lovop.getRange());
    
    AssignFilename af_data(valparse<string>(pars,"data_p"));    
    Y.eval(af_data);

    //evaluating operator
    OperatorEvaluation<ireal> opeval_lovop(towed_lovop,X);
    Y.copy(opeval_lovop.getValue());

    */


    // writing out
    string outfile = valparse<string>(pars,"outfile","");
      
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
	  cout << str.str();
	}
#ifdef IWAVE_USE_MPI
      }
      MPI_Barrier(retrieveGlobalComm());
    }
#endif
  }
  catch(RVLException &e){
    e << "ERROR from test_kernel()!\n"
      << "Dumping output stream:\n"
      << str.str();
    throw e;
  }
}
