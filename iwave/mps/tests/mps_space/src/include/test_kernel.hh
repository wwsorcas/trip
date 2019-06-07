#include "MPS_Space_Examples.hh"

#define CSQ "csq"
#define BULKMOD "bulkmod"
#define BUOYANCY "buoyancy"
#define SOURCE_P "source_p"
#define SOURCE_V0 "source_v0"
#define SOURCE_V1 "source_v1"
#define DATA_P "data_p"
#define DATA_V0 "data_v0"
#define DATA_V1 "data_v1"

using RVL::valparse;
using TSOpt::IWaveSpace;
using TSOpt::add_to_pars;
using TSOpt::CanScal_MPS_Space;
using TSOpt::CanVec_MPS_Space;
using TSOpt::ExVec_MPS_Space;
using TSOpt::Scal_MPS_Space;

template< class T_MPS_Space >
void test_kernel(PARARRAY &pars, FILE &stream){

  std::stringstream str;

  try {

    str << "\n////////////////////////////////////////////////////////////\n"
	<< "////////////////////////////////////////////////////////////\n"
	<< "Testing MPS_Space\n";
#ifdef IWAVE_USE_MPI
    str << "rank = "<< retrieveGlobalRank() <<"\n";
#endif
    str << "////////////////////////////////////////////////////////////\n"
	<< "////////////////////////////////////////////////////////////\n";

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
    mks.MPS_ref   = "MPS_ref";

    //setting mks.grid_file
    if( driver_type.compare("acd")==0 ){
      mks.grid_file = CSQ;
    }
    else{
      mks.grid_file = BULKMOD;
    }
    //setting mks.hdr_files and mks.G_files
    if( data_type.compare("pressure")==0 ){
      mks.hdr_files.push_back(DATA_P);
    }
    else{
      mks.hdr_files.push_back(DATA_V0);
      mks.hdr_files.push_back(DATA_V1);
    }

    //extracting RHS_files
    vector<string> RHS_keys;
    if( src_type.compare("pressure")==0 ){
      RHS_keys.push_back(SOURCE_P);
    }
    else{
      RHS_keys.push_back(SOURCE_V0);
      RHS_keys.push_back(SOURCE_V1);
    }

    //Empty constructor case
    str << "\n[] Testing empty constructor.\n\n";
    T_MPS_Space sp0;
    sp0.write(str);


    //Constructor with make=false case
    str << "\n[] Testing constructor with make=false\n\n";
    add_to_pars(pars,mks.MPS_file,valparse<string>(pars,"MPS_file_1"));
    T_MPS_Space sp1(mks,pars,false);
    sp1.write(str);
    sp1.print_info(str);
   

    //Constructor with make=true case
    str << "\n[] Testing constructor with make=true\n\n";
    add_to_pars(pars,mks.MPS_file,valparse<string>(pars,"MPS_file_2"));
    T_MPS_Space sp2(mks,pars,true);
    sp2.write(str);
    sp2.print_info(str);

    /*    
    //Testing gen_RHS_Space
    str << "\n[] Testing gen_RHS_Space method\n\n";
    IWaveInfo ic;
    int a_ord = valparse<int>(pars,mks.appx_ord);
    vector<string> RHS_files(0);
    for( int i=0; i<RHS_keys.size(); i++){
      RHS_files.push_back( valparse<string>(pars,RHS_keys[i]) );
    }
    std::shared_ptr<IWaveSpace> rhs_ptr= sp1.gen_RHS_Space(a_ord,
							   RHS_keys,
							   RHS_files,
							   ic);    
    rhs_ptr->write(str);
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
