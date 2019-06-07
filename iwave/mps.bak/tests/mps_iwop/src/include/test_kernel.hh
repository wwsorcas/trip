#include "iwop.hh"
#include "adjtest.hh"

#include "MPS_to_RHS.hh"
#include "MPS_iwop.hh"
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
#define DATA_P_G "data_p_g"
#define DATA_V0_G "data_v0_g"
#define DATA_V1_G "data_v1_g"
#define DATA_P_RE "data_p_re"
#define DATA_V0_RE "data_v0_re"
#define DATA_V1_RE "data_v1_re"
#define MPS_ADJ "MPS_adj"
#define MPS_ADJ_RE "MPS_adj_re"

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
using RVL::LinearOpAdjEvaluation;

using TSOpt::IWaveLOVOp;
using TSOpt::IWaveSpace;
using TSOpt::CanScal_MPS_Space;
using TSOpt::CanVec_MPS_Space;
using TSOpt::ExVec_MPS_Space;
using TSOpt::Scal_MPS_Space;
using TSOpt::MPS_to_RHS;
using TSOpt::MPS_IWaveLOVOp;
using TSOpt::MPS_KEYS;

template< class T_MPS_Space >
void test_kernel(PARARRAY &pars, FILE &stream, string test_info){

  stringstream str;
  string message;

  string jobname = valparse<string>(pars,"jobname");

  try {
    str << "\n============================================================\n"
	<< "Test name: "<< jobname <<"\n"
	<< "Runnning test for MPS_iwop.\n"
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

    //setting mps keys
    MPS_KEYS mks;
    mks.appx_ord   = "appx_ord";
    mks.MPS_file   = "MPS_file";
    mks.delta_file = "MPS_delta";

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
      mks.G_files.push_back(DATA_V0_G);
    }


    //////////////////////
    // Forward map test //
    //////////////////////

    //initializing ops and spaces
    T_MPS_Space mps_sp(mks,pars);
    MPS_to_RHS M2R(mks,pars,mps_sp);
    IWaveLOVOp IWOP(M2R.get_pars(),&stream);
    LinCompLOVOp<float> LOVOP(M2R,IWOP);
    MPS_IWaveLOVOp<T_MPS_Space> MPS_LOVOP(mks,pars,&stream);   

#ifdef VERBOSE_TEST_KERNEL_HH
    str << "After initializing ops\n";
    LOVOP.write(str);
    MPS_LOVOP.write(str);
#endif

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
    Vector<ireal> y_re(MPS_LOVOP.getRange());
    Components<ireal> cy_re(y_re);

    if(data_type.compare("pressure")==0){
      AssignFilename af_data_p(valparse<string>(pars,DATA_P));
      cy[0].eval(af_data_p);
      AssignFilename af_data_p_re(valparse<string>(pars,DATA_P_RE));
      cy_re[0].eval(af_data_p_re);
    }
    else{
      AssignFilename af_data_v0(valparse<string>(pars,DATA_V0));
      AssignFilename af_data_v1(valparse<string>(pars,DATA_V1));
      cy[0].eval(af_data_v0);
      cy[1].eval(af_data_v1);
      AssignFilename af_data_v0_re(valparse<string>(pars,DATA_V0_RE));
      AssignFilename af_data_v1_re(valparse<string>(pars,DATA_V1_RE));
      cy_re[0].eval(af_data_v0_re);
      cy_re[1].eval(af_data_v1_re);
    }

#ifdef VERBOSE_TEST_KERNEL_HH
    str << "\n[] Making data via IWaveLOVOp:\n";
#endif 
    //computing data via IWaveLOVOp
    OperatorEvaluation<ireal> opeval(LOVOP,x);
    y.copy(opeval.getValue());

#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif


#ifdef VERBOSE_TEST_KERNEL_HH
    str << "\n[] Making data via MPS_IWaveLOVOp (convolution):\n\n";
#endif
    //computing data via MPS_IWaveLOVOp (convolution)
    OperatorEvaluation<ireal> opeval_re(MPS_LOVOP,x);
    y_re.copy(opeval_re.getValue());

#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif
    

#ifdef VERBOSE_TEST_KERNEL_HH
    str << "\n[] Forward residual test:\n\n";
#endif
    Vector<ireal> fwd_res(y_re);
    fwd_res.copy(y_re);
    fwd_res.linComb(-1,y);
    float norm_fwd_res = fwd_res.norm() / y.norm();

    //running fwd residual test
    bool fwd_test = (fabs(norm_fwd_res)<MY_TOL);
    if( fwd_test ) message="\n\033[1;32m[PASSED]\033[0m";
    else message="\n\033[1;31m[FAILED]\033[0m";

    message += " Forward residual test\n";
    str  << message;

#ifdef VERBOSE_TEST_KERNEL_HH
    str << "|data - G*f|/|data| = " << norm_fwd_res << "\n"
	<< "|data| = " << y.norm()    << "\n"
	<< "|G*f| = "  << y_re.norm() << "\n";
    
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
	  << "|data| = " << y.norm()    << "\n"
	  << "|G*f| = "  << y_re.norm() << "\n";
      
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
    // Adjoint test for MPS_IWaveLOVOp //
    /////////////////////////////////////
#ifdef VERBOSE_TEST_KERNEL_HH
    str << "\n[] Running adjoint test for MPS_IWaveLOVOp:\n\n";
#endif
    LinearRestrictOp<float> MPS_LOP(MPS_LOVOP,cx[0]);
    stringstream res;
    RVLRandomize<float> rnd(getpid(),-1.0,1.0);
    bool adj_test=RVL::AdjointTest<float>(MPS_LOP,rnd,res);

    if(adj_test) message="\n\033[1;32m[PASSED]\033[0m";
    else message="\n\033[1;31m[FAILED]\033[0m";    
    message += " Adjoint test for MPS_IWaveLOVOp\n";
    str  << message;

#ifdef VERBOSE_TEST_KERNEL_HH
    str << res.str();
#else
    if(!adj_test) str << res.str();
#endif


#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif


    // Adjoint of IWaveLOVOp only available for acd currently
    if( driver_type.compare("acd")==0 ){

      /////////////////////////////////
      // Adjoint test for IWaveLOVOp //
      /////////////////////////////////
#ifdef VERBOSE_TEST_KERNEL_HH
      str << "\n[] Running adjoint test for IWaveLOVOp:\n\n";
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
      Vector<ireal> x1_adj(LOP.getDomain());
      AssignFilename af_src_adj(valparse<string>(pars,MPS_ADJ));
      x1_adj.eval(af_src_adj);
      
      //computing data via IWaveLOVOp
      LinearOpAdjEvaluation<ireal> LOP_adj_eval(LOP,y);
      x1_adj.copy(LOP_adj_eval.getValue());
      
#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif


#ifdef VERBOSE_TEST_KERNEL_HH
      str << "\n[] Making adjoint data via MPS_IWaveLOVOp (convolution):\n\n";
#endif      
      //setting output
      Vector<ireal> x1_adj_re(MPS_LOP.getDomain());
      AssignFilename af_src_adj_re(valparse<string>(pars,MPS_ADJ_RE));
      x1_adj_re.eval(af_src_adj_re);
      
      //computing data via MPS_IWaveLOVOp (convolution)
      LinearOpAdjEvaluation<ireal> MPS_LOP_adj_eval(MPS_LOP,y);
      x1_adj_re.copy(MPS_LOP_adj_eval.getValue());


#ifdef VERBOSE_TEST_KERNEL_HH
      str << "\n[] Adjoint residual test\n\n";
#endif
      Vector<ireal> adj_res(x1_adj_re);
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

      
#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
      
    }      

    //writing out
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
