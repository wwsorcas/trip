#include "iwop.hh"
#include "blockop.hh"
#include "adjtest.hh"

#include "MPS_to_RHS.hh"
#include "MPS_Space_Examples.hh"
#include "MPS_conv.hh"

#define CSQ "csq"
#define BULKMOD "bulkmod"
#define BUOYANCY "buoyancy"
#define SOURCE_P "source_p"
#define SOURCE_V0 "source_v0"
#define SOURCE_V1 "source_v1"
#define DATA_P "data_p"
#define DATA_V0 "data_v0"
#define DATA_V1 "data_v1"
#define DATA_P_G "data_p_g"
#define DATA_V0_G "data_v0_g"
#define DATA_V1_G "data_v1_g"
#define DATA_P_RE "data_p_re"
#define DATA_V0_RE "data_v0_re"
#define DATA_V1_RE "data_v1_re"
#define MPS_ADJ "MPS_adj"
#define MPS_ADJ_RE "MPS_adj_re"
#define SOURCE_P_G "source_p_g"
#define SOURCE_V0_G "source_v0_g"
#define SOURCE_V1_G "source_v1_g"

#ifndef MY_TOL
#define MY_TOL 1e-3
#endif

//#define VERBOSE_TEST_KERNEL_HH

using RVL::valparse;
using RVL::AssignFilename;
using RVL::OperatorEvaluation;
using RVL::LinCompLOVOp;
using RVL::LinearRestrictOp;
using RVL::RVLRandomize;
using RVL::AdjLinearOp;
using RVL::LinearOpAdjEvaluation;

using TSOpt::IWaveLOVOp;
using TSOpt::IWaveSpace;
using TSOpt::CanScal_MPS_Space;
using TSOpt::CanVec_MPS_Space;
using TSOpt::ExVec_MPS_Space;
using TSOpt::Scal_MPS_Space;
using TSOpt::MPS_to_RHS;
using TSOpt::MPS_conv;
using TSOpt::MPS_KEYS;

template< class T_MPS_Space >
void test_kernel(PARARRAY &pars, FILE &stream, string test_info){

  stringstream str;
  string message;
  
  string jobname = valparse<string>(pars,"jobname");

  try {
    str << "\n============================================================\n"
	<< "Test name: "<< jobname <<"\n"
	<< "Runnning test for MPS_conv.\n"
	<< test_info <<"\n";
#ifdef IWAVE_USE_MPI
    str << "rank = "<< retrieveGlobalRank() <<"\n";
#endif
    str << "============================================================\n";

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
  
    //setting mps keys
    MPS_KEYS mks;
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
    //setting mks.hdr_files and mks.G_files
    if( data_type.compare("pressure")==0 ){
      mks.hdr_files.push_back(DATA_P);
      mks.G_files.push_back(DATA_P_G);
    }
    else{
      mks.hdr_files.push_back(DATA_V0);
      mks.hdr_files.push_back(DATA_V1);
      mks.G_files.push_back(DATA_V0_G);
      mks.G_files.push_back(DATA_V1_G);
    }
  

    ////////////////////////////////
    // Making data via IWaveLOVOp //
    ////////////////////////////////
#ifdef VERBOSE_TEST_KERNEL_HH
    str << "\n[] Making data via IWaveLOVOp:\n";
#endif

    //Initializing
    T_MPS_Space mps_sp(mks,pars);
    MPS_to_RHS M2R(mks,pars,mps_sp);
    IWaveLOVOp IWOP(M2R.get_pars(),&stream);
    LinCompLOVOp<float> LOVOP(M2R,IWOP);

    //setting inputs
    Vector<ireal> x(LOVOP.getDomain());
    Components<ireal> cx(x);
    Components<ireal> cx0(cx[0]);
    Components<ireal> cx1(cx[1]);

    if( driver_type.compare("acd")==0 ){
      AssignFilename af_csq(valparse<string>(pars,CSQ));
      cx0[0].eval(af_csq);
    }
    else{
      AssignFilename af_bmod(valparse<string>(pars,BULKMOD));
      AssignFilename af_buoy(valparse<string>(pars,BUOYANCY));
      cx0[0].eval(af_bmod);
      cx0[1].eval(af_buoy);
    }
    AssignFilename af_src(valparse<string>(pars,mks.MPS_file));
    cx1[0].eval(af_src);

    //setting outputs
    Vector<ireal> y(LOVOP.getRange());
    Components<ireal> cy(y);

    if(data_type.compare("pressure")==0){
      AssignFilename af_data_p(valparse<string>(pars,DATA_P));
      cy[0].eval(af_data_p);
    }
    else{
      AssignFilename af_data_v0(valparse<string>(pars,DATA_V0));
      AssignFilename af_data_v1(valparse<string>(pars,DATA_V1));
      cy[0].eval(af_data_v0);
      cy[1].eval(af_data_v1);
    }

    //computing data via IWaveLOVOp
    OperatorEvaluation<ireal> opeval(LOVOP,x);
    y.copy(opeval.getValue());

#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif
 
  
    ///////////////////////
    // making green data //
    ///////////////////////
#ifdef VERBOSE_TEST_KERNEL_HH
    str << "\n[] Making Green data:\n";
#endif

    PARARRAY *pars_2 = ps_new();
    ps_copy(&pars_2,pars);
    
    //adding key-val pairs for Green data simulation
    vector<string> keys(0);
    vector<string> vals(0);

    if(data_type.compare("pressure")==0){
      keys.push_back(DATA_P);
      vals.push_back(valparse<string>(pars,DATA_P_G));
    }
    else{
      keys.push_back(DATA_V0);
      vals.push_back(valparse<string>(pars,DATA_V0_G));
      keys.push_back(DATA_V1);
      vals.push_back(valparse<string>(pars,DATA_V1_G));
    }
    if(src_type.compare("pressure")==0){
      keys.push_back(SOURCE_P);
      vals.push_back(valparse<string>(pars,SOURCE_P_G));
    }
    else{
      keys.push_back(SOURCE_V0);
      vals.push_back(valparse<string>(pars,SOURCE_V0_G));
      keys.push_back(SOURCE_V1);
      vals.push_back(valparse<string>(pars,SOURCE_V1_G));
    }
    keys.push_back("MPS_file");
    vals.push_back(valparse<string>(pars,"MPS_file_g"));

    TSOpt::add_to_pars(*pars_2,keys,vals);
  
    //Initializing 
    T_MPS_Space mps_sp_g(mks,*pars_2);
    MPS_to_RHS M2R_g(mks,*pars_2,mps_sp_g);
    IWaveLOVOp IWOP_g(M2R_g.get_pars(),&stream);
    LinCompLOVOp<float> LOVOP_g(M2R_g,IWOP_g);
    
    //setting inputs
    Vector<ireal> x_g(LOVOP_g.getDomain());
    Components<ireal> cx_g(x_g);
    Components<ireal> cx0_g(cx_g[0]);
    Components<ireal> cx1_g(cx_g[1]);

    if( driver_type.compare("acd")==0 ){
      AssignFilename af_csq(valparse<string>(pars,CSQ));
      cx0_g[0].eval(af_csq);      
    }
    else{
      AssignFilename af_bmod(valparse<string>(pars,BULKMOD));
      AssignFilename af_buoy(valparse<string>(pars,BUOYANCY));
      cx0_g[0].eval(af_bmod);
      cx0_g[1].eval(af_buoy);
    }
    AssignFilename af_src_g(valparse<string>(*pars_2,"MPS_file_g"));
    cx1_g[0].eval(af_src_g);

    //setting outputs
    Vector<ireal> y_g(LOVOP_g.getRange());
    Components<ireal> cy_g(y_g);

    if( data_type.compare("pressure")==0 ){
      AssignFilename af_data_p_g(valparse<string>(*pars_2,DATA_P_G));
      cy_g[0].eval(af_data_p_g);
    }
    else{
      AssignFilename af_data_v0_g(valparse<string>(*pars_2,DATA_V0_G));
      AssignFilename af_data_v1_g(valparse<string>(*pars_2,DATA_V1_G));
      cy_g[0].eval(af_data_v0_g);
      cy_g[1].eval(af_data_v1_g);
    }

    //computing green data
    OperatorEvaluation<ireal> opeval_g(LOVOP_g,x_g);
    y_g.copy(opeval_g.getValue());

#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif
 

    /////////////////////////////////
    // making data via convolution //
    /////////////////////////////////
#ifdef VERBOSE_TEST_KERNEL_HH
    str << "\n[] Making data via convolution:\n\n";
#endif
    //initialize
    MPS_conv CONV( pars,
		   y_g,
		   M2R.get_MPS_Domain(),
		   IWOP.getIWaveRange() );
    
    //setting outputs
    Vector<ireal> y_re(CONV.getRange());
    Components<ireal> cy_re(y_re);
   
    if( data_type.compare("pressure")==0 ){
      AssignFilename af_data_p_re(valparse<string>(pars,DATA_P_RE));
      cy_re[0].eval(af_data_p_re);
    }
    else{
      AssignFilename af_data_v0_re(valparse<string>(pars,DATA_V0_RE));
      AssignFilename af_data_v1_re(valparse<string>(pars,DATA_V1_RE));
      cy_re[0].eval(af_data_v0_re);
      cy_re[1].eval(af_data_v1_re);
    }

    //computing data via convolution 
    int MPS_idx = valparse<int>(pars,"MPS_idx");
    CONV.set_MPS_idx(MPS_idx);
    OperatorEvaluation<ireal> opeval_conv(CONV,cx[1]);
    y_re.copy(opeval_conv.getValue());

#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
    for( int r=0; r<retrieveGlobalSize(); r++ ){
      if( retrieveGlobalRank()==r ){
	cerr <<"rank="<<r<<", after conv eval\n";
      }
      MPI_Barrier(retrieveGlobalComm());
    } 
#endif
  

    //////////////////
    // FORWARD TEST //
    //////////////////
#ifdef VERBOSE_TEST_KERNEL_HH
    str << "\n[] Running forward residual test:\n\n";
#endif

    Vector<ireal> fwd_res(y_re);
    fwd_res.copy(y_re);
    fwd_res.linComb(-1,y);
    float norm_fwd_res = fwd_res.norm() / y.norm();

    bool fwd_test = (fabs(norm_fwd_res)<MY_TOL);
    if( fwd_test ) message="\n\033[1;32m[PASSED]\033[0m";
    else message="\n\033[1;31m[FAILED]\033[0m";
    
    message += " Forward residual test\n";
    str  << message;

#ifdef VERBOSE_TEST_KERNEL_HH
    str << "|data - G*f|/|data| = " << norm_fwd_res << "\n"
	<< "|data| = " << y.norm() << "\n"
	<< "|G*f| = " << y_re.norm() << "\n";
    
    if(data_type.compare("pressure")==0){
      str << "See IWaveLOVOp data file = "
	  << valparse<string>(pars,DATA_P) << "\n";
      str << "See MPS_IWaveLOVOp data file = "
	  << valparse<string>(pars,DATA_P_RE) << "\n";
    }
    else{
      str << "See IWaveLOVOp data files = "
	  << valparse<string>(pars,DATA_V0) << ", "
	  << valparse<string>(pars,DATA_V1) << "\n";
      str << "See MPS_IWaveLOVOp data files = "
	  << valparse<string>(pars,DATA_V0_RE) << ", "
	  << valparse<string>(pars,DATA_V1_RE) << "\n";
    }
#else
    if(!fwd_test){
      str << "|data - G*f|/|data| = " << norm_fwd_res << "\n"
	  << "|data| = " << y.norm() << "\n"
	  << "|G*f| = " << y_re.norm() << "\n";
      
      if(data_type.compare("pressure")==0){
	str << "See IWaveLOVOp data file = "
	    << valparse<string>(pars,DATA_P) << "\n";
	str << "See MPS_IWaveLOVOp data file = "
	    << valparse<string>(pars,DATA_P_RE) << "\n";
      }
      else{
	str << "See IWaveLOVOp data files = "
	    << valparse<string>(pars,DATA_V0) << ", "
	    << valparse<string>(pars,DATA_V1) << "\n";
	str << "See MPS_IWaveLOVOp data files = "
	    << valparse<string>(pars,DATA_V0_RE) << ", "
	    << valparse<string>(pars,DATA_V1_RE) << "\n";
      }
    }
#endif

  
#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif


    /////////////////////////////////////
    // Adjoint test for convolution op //
    /////////////////////////////////////
#ifdef VERBOSE_TEST_KERNEL_HH
    str << "\n[] Running adjoint test on convolution op:\n\n";
#endif

    stringstream res;
    RVLRandomize<float> rnd(getpid(),-1.0,1.0);
    bool adj_test = RVL::AdjointTest<float>(CONV,rnd,res);

    if(adj_test) message="\n\033[1;32m[PASSED]\033[0m";
    else message="\n\033[1;31m[FAILED]\033[0m";
    
    message += " Adjoint test for MPS_conv op\n";
    str  << message;
  
#ifdef VERBOSE_TEST_KERNEL_HH
    str << res.str();
#else
    if(!adj_test) str << res.str();
#endif


#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif

    /* Doesn't make sense to run adjoint residual tests
       at this point. This should be done for MPS_IWaveLOVOp.

    // Adjoint operations only available for acd currently
    if( driver_type.compare("acd")==0 ){

      ////////////////////////////
      // Adjoint test for LOVOp //
      ////////////////////////////
#ifdef VERBOSE_TEST_KERNEL_HH
      str << "\n[] Running adjoint test for IWaveLOVOp:\n";
#endif
      LinearRestrictOp<float> LOP(LOVOP,cx[0]);
      stringstream res_iw;
      RVLRandomize<float> rnd_iw(getpid(),-1.0,1.0);
      bool adj_test_iw=RVL::AdjointTest<float>(LOP,rnd_iw,res_iw);

      if(adj_test_iw) message="\n\033[1;32m[PASSED]\033[0m";
      else message="\n\033[1;31m[FAILED]\033[0m";

      message += " Adjoint test for IWaveLOVOp\n";
      str  << message;
        
#ifdef VERBOSE_TEST_KERNEL_HH
      str << res_iw.str();
#else
      if(!adj_test_iw) str << res_iw.str();
#endif

#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif
        
    
      ///////////////////////////
      // Adjoint residual test //
      ///////////////////////////
#ifdef VERBOSE_TEST_KERNEL_HH
      str << "\n[] Making adjoint data via IWaveLOVOp:\n";
#endif
      //setting output
      Vector<float> x1_adj(LOP.getDomain());
      AssignFilename af_src_adj(valparse<string>(pars,MPS_ADJ));
      x1_adj.eval(af_src_adj);
      
      //computing adjoint data via IWaveLOVOp
      LinearOpAdjEvaluation<ireal> LOP_adj_eval(LOP,y);
      x1_adj.copy(LOP_adj_eval.getValue());

#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif
      

#ifdef VERBOSE_TEST_KERNEL_HH
      str << "\n[] Making adjoint data via convolution:\n\n";
#endif
      //setting output
      Vector<float> x1_adj_re(LOP.getDomain());
      AssignFilename af_src_adj_re(valparse<string>(pars,MPS_ADJ_RE));    
      x1_adj_re.eval(af_src_adj_re);
      
      //computing adjoint data via convolution
      LinearOpAdjEvaluation<ireal> conv_adj_eval(CONV,y);
      x1_adj_re.copy(conv_adj_eval.getValue());
      
      //residual
      Vector<ireal> adj_res(x1_adj_re.getSpace());
      adj_res.copy(x1_adj_re);
      adj_res.linComb(-1,x1_adj);
      float norm_adj_res = adj_res.norm() / x1_adj.norm();

      //adjoint residual test
      bool aadj_test = (fabs(norm_adj_res)<MY_TOL);
      if( aadj_test ) message="\n\033[1;32m[PASSED]\033[0m";
      else message="\n\033[1;31m[FAILED]\033[0m";

      message += " Adjoint residual test\n";
      str  << message;

#ifdef VERBOSE_TEST_KERNEL_HH
      str << "|src_adj - G^T*data|/|src_adj| = " << norm_adj_res << "\n"
	  << "|src_adj| = " << x1_adj.norm()    << "\n"
	  << "|G^T*data| = "  << x1_adj_re.norm() << "\n";
      
      str << "See IWaveLOVOp adjoint data file = "
	  << valparse<string>(pars,"MPS_adj") << "\n";
      str << "See MPS_IWaveLOVOp adjoint data file = "
	  << valparse<string>(pars,"MPS_adj_re") << "\n";
#else
      if(!aadj_test){
	str << "|src_adj - G^T*data|/|src_adj| = " << norm_adj_res << "\n"
	    << "|src_adj| = " << x1_adj.norm()    << "\n"
	    << "|G^T*data| = "  << x1_adj_re.norm() << "\n";
	
	str << "See IWaveLOVOp adjoint data file = "
	    << valparse<string>(pars,"MPS_adj") << "\n";
	str << "See MPS_IWaveLOVOp adjoint data file = "
	    << valparse<string>(pars,"MPS_adj_re") << "\n";
      }
#endif
    }
    */

#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif

    ps_delete(&pars_2);

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
	  cerr << str.str();
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
