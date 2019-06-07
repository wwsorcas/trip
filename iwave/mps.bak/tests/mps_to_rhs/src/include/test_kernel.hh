#include "adjtest.hh"

#include "MPS_Space_Examples.hh"
#include "MPS_to_RHS.hh"

#define CSQ "csq"
#define BULKMOD "bulkmod"
#define BUOYANCY "buoyancy"
#define SOURCE_P "source_p"
#define SOURCE_V0 "source_v0"
#define SOURCE_V1 "source_v1"
#define DATA_P "data_p"
#define DATA_V0 "data_v0"
#define DATA_V1 "data_v1"

#define VERBOSE_MJB

using RVL::valparse;
using TSOpt::IWaveSpace;
using TSOpt::add_to_pars;
using TSOpt::MPS_to_RHS;
using TSOpt::CanScal_MPS_Space;
using TSOpt::CanVec_MPS_Space;
using TSOpt::ExVec_MPS_Space;
using TSOpt::Scal_MPS_Space;

template< class T_MPS_Space >
void test_kernel(PARARRAY &pars, FILE &stream, string test_info){

  std::stringstream str;
  string message;

  try {

    string jobname = valparse<string>(pars,"jobname");
    str << "\n============================================================\n"
	<< "Test name: "<< jobname <<"\n"
	<< "Runnning test for MPS_to_RHS.\n"
	<< test_info <<"\n";
#ifdef IWAVE_USE_MPI
    str << "rank = "<< retrieveGlobalRank() <<"\n";
#endif
    str << "============================================================\n";

    vector<string> RHS_files(0);
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
    mks.appx_ord  = "appx_ord";
    mks.MPS_file  = "MPS_file";

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


    //Initializing MPS space
    T_MPS_Space mps_sp(mks,pars);    


    //Constructor with make_RHS=true case and given RHS filename
#ifdef VERBOSE_MJB
    str << "\n[] Testing constructor with make_RHS=true w/ given RHS filenames\n\n";
#endif
    
    if( src_type.compare("pressure")==0 ){
      RHS_files.push_back(valparse<string>(pars,"RHS_p"));
      add_to_pars(pars,mks.RHS_files,RHS_files);
      RHS_files.clear();
    }
    else{
      RHS_files.push_back(valparse<string>(pars,"RHS_v0"));
      RHS_files.push_back(valparse<string>(pars,"RHS_v1"));
      add_to_pars(pars,mks.RHS_files,RHS_files);
      RHS_files.clear();      
    }

    MPS_to_RHS m2r_1(mks,pars,mps_sp);

#ifdef VERBOSE_MJB
    m2r_1.write(str);
    m2r_1.print_info(str);
#endif

    //Constructor with make_RHS=true case and w/o given RHS filename
#ifdef VERBOSE_MJB
    str << "\n[] Testing constructor with make_RHS=true w/o given RHS filenames\n\n";
#endif

    if( src_type.compare("pressure")==0 ){
      RHS_files.push_back("empty");
      add_to_pars(pars,mks.RHS_files,RHS_files);
      RHS_files.clear();
    }
    else{
      RHS_files.push_back("empty");
      RHS_files.push_back("empty");
      add_to_pars(pars,mks.RHS_files,RHS_files);
      RHS_files.clear();      
    }

    MPS_to_RHS m2r_2(mks,pars,mps_sp);

#ifdef VERBOSE_MJB
    m2r_2.write(str);
    m2r_2.print_info(str);
#endif


    //Adjoint test
#ifdef VERBOSE_MJB
    str << "\n[] Running adjoint test:\n";    
#endif
    std::stringstream res;
    RVL::RVLRandomize<float> rnd(getpid(),-1.0,1.0);
    bool adj_test;
    adj_test = RVL::AdjointTest<float>(m2r_1,rnd,res);
    
    if(adj_test) message="\n\033[1;32m[PASSED]\033[0m";
    else message="\n\033[1;31m[FAILED]\033[0m";
      
    message += " Adjoint test for M_to_R";
    str  << message << "\n";
    if(!adj_test) str << res.str();
    //cerr << message << "; rank="<<retrieveGlobalRank()<<"\n";
#ifdef VERBOSE_MJB
    str << res.str();
#endif

    // writing out
    string outfile = valparse<string>(pars,"outfile","");
#ifdef IWAVE_USE_MPI
    for( int r=0; r<retrieveGlobalSize(); r++ ){
      if( retrieveGlobalRank()==r ){
	cerr << message << "; "
	     << "rank="<<r<<"\n";
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
      << "Dumping output stream: \n"
      << str.str();
    throw e;
  }
}
