#include "adjtest.hh"
#include "op.hh"

#include "MPS_iwop.hh"
#include "MPS_Space_Examples.hh"
#include "MPS_frac_cal.hh"

#define CSQ "csq"
#define BULKMOD "bulkmod"
#define BUOYANCY "buoyancy"
#define SOURCE_P "source_p"
#define SOURCE_V0 "source_v0"
#define SOURCE_V1 "source_v1"
#define DATA_P "data_p"
#define DATA_V0 "data_v0"
#define DATA_V1 "data_v1"

#ifndef TOL
#define TOL 1e-5
#endif

//#define VERBOSE_TEST_KERNEL_HH

using RVL::valparse;
using RVL::AssignFilename;
using RVL::OperatorEvaluation;
using RVL::NormalLinearOp;
using RVL::AdjLinearOp;

using TSOpt::CanScal_MPS_Space;
using TSOpt::CanVec_MPS_Space;
using TSOpt::ExVec_MPS_Space;
using TSOpt::Scal_MPS_Space;
using TSOpt::MPS_IWaveLOVOp;
using TSOpt::MPS_frac_cal;


template< class T_MPS_Space >
void test_kernel(PARARRAY &pars, FILE &stream, string test_info) {

  std::stringstream str;
  string message;

  string jobname = valparse<string>(pars,"jobname");

  try {

    str << "\n============================================================\n"
	<< "Test name: "<< jobname <<"\n"
	<< "Runnning test for MPS_frac_cal.\n"
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

    TSOpt::MPS_KEYS mks;
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
      mks.G_files.push_back("data_p_g");
    }
    else{
      mks.hdr_files.push_back(DATA_V0);
      mks.hdr_files.push_back(DATA_V1);
      mks.G_files.push_back("data_v0_g");
      mks.G_files.push_back("data_v1_g");
    }


    //Initializing
#ifdef VERBOSE_TEST_KERNEL_HH
    str << "\n[] Initializing:\n";
#endif

    MPS_IWaveLOVOp<T_MPS_Space> mps_lovop(mks,pars,&stream);
    T_MPS_Space mps_sp(mps_lovop.getLinDomain());

    float order_0 = valparse<float>(pars,"order_0");
    float order_d = valparse<float>(pars,"order_d");
    float c = valparse<float>(pars,"c_max",1.0f);

    MPS_frac_cal D( pars,
		    mps_sp,
		    c,
		    order_0,
		    order_d );

    RVL::InvLinearOp<ireal> Dinv(D);

    /*
    NormalLinearOp<ireal> Minv(D);
    AdjLinearOp<ireal> Dinv_adj(Dinv);
    NormalLinearOp<ireal> M(Dinv_adj);
    */

#ifdef VERBOSE_TEST_KERNEL_HH
    str << "\n- writing D op: ";
    D.write(str);
    str << "\n- writing Dinv op: ";
    Dinv.write(str);
    str << "\n";
#endif

#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif

    // Testing diffint ops
#ifdef VERBOSE_TEST_KERNEL_HH
    str << "\n[] Testing MPS_frac_cal ops:\n";
#endif

    Vector<ireal> x(mps_sp);
    AssignFilename af_x(valparse<string>(pars,mks.MPS_file));
    x.eval(af_x);   

    Vector<ireal> y_D(mps_sp);
    AssignFilename af_y_D(valparse<string>(pars,"y_D"));
    y_D.eval(af_y_D);    
    Vector<ireal> y_D_Dinv(mps_sp);
    AssignFilename af_y_D_Dinv(valparse<string>(pars,"y_D_Dinv"));
    y_D_Dinv.eval(af_y_D_Dinv);    
    Vector<ireal> y_Dinv(mps_sp);
    AssignFilename af_y_Dinv(valparse<string>(pars,"y_Dinv"));
    y_Dinv.eval(af_y_Dinv);    
    Vector<ireal> y_Dinv_D(mps_sp);
    AssignFilename af_y_Dinv_D(valparse<string>(pars,"y_Dinv_D"));
    y_Dinv_D.eval(af_y_Dinv_D);    

    OperatorEvaluation<ireal> eval_D(D,x);
    y_D.copy(eval_D.getValue());
    OperatorEvaluation<ireal> eval_D_Dinv(Dinv,y_D);
    y_D_Dinv.copy(eval_D_Dinv.getValue());
    
    Vector<ireal> residual(x);
    residual.copy(x);
    residual.linComb(-1,y_D_Dinv);
    float norm_res = residual.norm() / x.norm();
    
    bool inv_test = (fabs(norm_res)<TOL);
    if( inv_test ) message="\n\033[1;32m[PASSED]\033[0m";
    else message="\n\033[1;31m[FAILED]\033[0m";

    message += " Dinv*D test\n";    
    str  << message;
    if(!inv_test){
      str << "\n"
	  << "|x - Dinv*D*x|/|x| = " << norm_res << "\n"
	  << "|x| = " << x.norm() << "\n"
	  << "|Dinv*D*x| = " << y_D_Dinv.norm() << "\n";
    }


#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif

    Vector<ireal> residual_inv(x);
    residual_inv.copy(x);
    residual_inv.linComb(-1,y_Dinv_D);
    float norm_res_inv = residual_inv.norm() / x.norm();

    inv_test = (fabs(norm_res_inv)<TOL);
    if( inv_test ) message="\n\033[1;32m[PASSED]\033[0m";
    else message="\n\033[1;31m[FAILED]\033[0m";

    message += " D*Dinv test\n";
    str  << message;
    if(!inv_test){
      str  << "\n"
	   << "|x - D*Dinv*x|/|x| = " << norm_res_inv << "\n"
	   << "|x| = " << x.norm() << "\n"
	   << "|D*Dinv*x| = " << y_Dinv_D.norm() << "\n";
    }


#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif

    // Adjoint Test
#ifdef VERBOSE_TEST_KERNEL_HH
    str << "\n[] Running adjoint test:\n";
#endif

    bool adj_test;

    std::stringstream res;
    RVL::RVLRandomize<float> rnd(getpid(),-1.0,1.0);
    adj_test = RVL::AdjointTest<float>(D,rnd,res);
    
    if(adj_test) message="\n\033[1;32m[PASSED]\033[0m";
    else message="\n\033[1;31m[FAILED]\033[0m";

    message += " Adjoint test for D\n";
    str << message;
    if(!adj_test) str << res.str();

#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif

    std::stringstream res_inv;
    RVL::RVLRandomize<float> rnd_inv(getpid(),-1.0,1.0);
    adj_test = RVL::AdjointTest<float>(Dinv,rnd_inv,res_inv);

    if(adj_test) message="\n\033[1;32m[PASSED]\033[0m";
    else message="\n\033[1;31m[FAILED]\033[0m";

    message += " Adjoint test for Dinv\n";
    str  << message;
    if(!adj_test) str << res.str();
    
#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif

    /*
    // Testing M and M_inv
#ifdef VERBOSE_TEST_KERNEL_HH
    str << "\n[] Testing M and Minv\n";
#endif

    Vector<ireal> y_Minv(mps_sp);
    AssignFilename af_y_Minv(valparse<string>(pars,"y_Minv"));
    y_Minv.eval(af_y_Minv);    
    Vector<ireal> y_Minv_M(mps_sp);
    AssignFilename af_y_Minv_M(valparse<string>(pars,"y_Minv_M"));
    y_Minv_M.eval(af_y_Minv_M);

    OperatorEvaluation<ireal> eval_Minv(Minv,x);
    y_Minv.copy(eval_Minv.getValue());
    OperatorEvaluation<ireal> eval_Minv_M(M,y_Minv);
    y_Minv_M.copy(eval_Minv_M.getValue());

    Vector<ireal> residual_M(x);
    residual_M.copy(x);
    residual_M.linComb(-1,y_Minv_M);
    float norm_res_M = residual_M.norm() / x.norm();
    
    bool M_test = ( fabs(norm_res_M)<TOL*1e+3 );
    if( M_test ) message="\n\033[1;32m[PASSED]\033[0m";
    else message="\n\033[1;31m[FAILED]\033[0m";
    
    message += " M*Minv test\n";
    str  << message;
    if(!M_test){
      str << "\n"
	  << "|x - M*Minv*x|/|x| = " << norm_res_M << "\n"
	  << "|x| = " << x.norm() << "\n"
	  << "|M*Minv*x| = " << y_Minv_M.norm() << "\n";
    }
    */

#ifdef IWAVE_USE_MPI
    MPI_Barrier(retrieveGlobalComm());
#endif


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
