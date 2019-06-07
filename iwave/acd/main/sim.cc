#include "acd_defn.hh"
#include "grid.h"
#include "gridpp.hh"
#include "gridops.hh"
#include "iwop.hh"
#include "functions.hh"
#include "op.hh"
#include "linop_apps.hh"
#include "ls.hh"
#include "blockop.hh"
#include "cgnealg.hh"
#include "LBFGSBT.hh"
#include "acd_selfdoc.h"
#include <par.h>
#include "segyops.hh"
#include "adjtest.hh"
#include <omp.h>

//#define GTEST_VERBOSE
IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csq",    0, 1, 1 },
  {"data",   1, 0, 2 },
  {"source", 1, 1, 2 },
  {"movie",  1, 0, 2 },
  {"",       0, 0, 0}
};

using RVL::parse;
using RVL::RVLException;
using RVL::Vector;
using RVL::Components;
using RVL::Operator;
using RVL::OperatorEvaluation;
using RVL::LinearOp;
using RVL::LinearOpFO;
using RVL::RestrictOp;
using RVL::TangentMap;
using RVL::LinearRestrictOp;
using RVL::SymmetricBilinearOp;
using RVL::LinearBilinearOp;
using RVL::AssignFilename;
using RVL::AssignParams;
using RVL::RVLRandomize;
using RVL::ScaleOpFwd;
using RVL::TensorOp;
using RVL::CompLOVOp;
using RVL::FunctionObject;
using RVL::AdjointTest;
using RVLUmin::CGNEPolicy;
using RVLUmin::CGNEPolicyData;
using RVLUmin::LBFGSBT;
using RVLUmin::CGNEAlg;
using TSOpt::IWaveEnvironment;
using TSOpt::IWaveTree;
using TSOpt::IWaveSampler;
using TSOpt::IWaveSim;
using TSOpt::TASK_RELN;
using TSOpt::IOTask;
using TSOpt::IWaveOp;
using TSOpt::IWaveLOVOp;
using TSOpt::SEGYTaper;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
using TSOpt::MPISEGYSpace;
using RVL::MPISerialFunctionObject;
#else
using TSOpt::GridSpace;
using TSOpt::SEGYSpace;
#endif
using TSOpt::GridExtendOp;

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
    // internal tests to ensure valid pointers returned.
    IWaveEnvironment(argc, argv, 0, &pars, &stream);

    if (retrieveGlobalRank()==0 && argc<2) {
      pagedoc();
      exit(0);
    }

#ifdef _OPENMP
    omp_set_num_threads(valparse<int>(*pars,"num_threads",1));
    std::cout << "Number of OMP threads in use: "
	      << omp_get_max_threads() << std::endl;
#endif
    
#ifdef IWAVE_USE_MPI
    if (retrieveGroupID() == MPI_UNDEFINED) {
      fprintf(stream,"NOTE: finalize MPI, cleanup, exit\n");
      fflush(stream);
    }
    else {
#endif
      // cerr<<"001\n";
      
      /* This driver offers physical (non-extended) simulation,
         extended simulation, and the mixed case of linearized
         extended simulation about a physical background. To detect
         this last situation, must unpack the reference and
         perturbation
      */

      int deriv = valparse<int>(*pars,"deriv",0);
      int adjoint = valparse<int>(*pars,"adjoint",0);

      // key for source
      std::string src = "";
      
      // the following keys refer to square-velocity and its
      // perturbations
      std::string refr = "";
      std::string extdrefr = "";
      std::string pert = "";
      std::string pert2 = "";
      std::string extdpert2 = "";
      std::string refrkey = "csq";
      std::string srckey = "source";
      std::string pertkey = "csq_d1";
      std::string pert2key = "csq_d2";
      
      if (deriv==1 && adjoint) pertkey="csq_b1";
      if (deriv==2 && adjoint) pert2key="csq_b2";

      // name of background csq file
      refr = valparse<std::string>(*pars,refrkey);
      // name of first pert csq file
      if (deriv > 0) pert  = valparse<std::string>(*pars,pertkey);
      // name of second pert csq file
      if (deriv > 1) pert2 = valparse<std::string>(*pars,pert2key);
      // name of source file
      src = valparse<std::string>(*pars,srckey);

      std::stringstream err;

      // the mixed case is significant only for deriv > 0
      // mixed = flag set, else either physical or completely extended
      int physbg = 0; 
      if (deriv && (retrieveGlobalRank()==0)) {
	// construct reference, perturbational grids
	grid grefr; init_default_grid(&grefr);
	read_grid(&grefr,refr.c_str(),stream);
	grid gpert; init_default_grid(&gpert);
	read_grid(&gpert,pert.c_str(),stream);
	// sanity check
	if (grefr.dim != gpert.dim) {
	  RVLException e;
	  e<<"ERROR: acd, deriv>0\n";
	  e<<"  ref, pert grids have different physical dimensions\n";
	  throw e;
	}
	// set physical bg flag if necessary
	if (grefr.gdim < gpert.gdim) physbg=1;
      }
#ifdef IWAVE_USE_MPI
      // broadcast flag
      MPI_Bcast(&physbg,1,MPI_INT,0,retrieveGlobalComm());
#endif
      // record physical background flag in param 
      // table for transfer to iwave
      if (ps_slint(*pars,"physbg",physbg)) {
	RVLException e;
	e<<"ERROR: acd, deriv>0\n"; 
	e<<"  failed to add physbg=1 to par table\n";
	throw e;
      }

#ifdef IWAVE_USE_MPI
      // wait for everyone to catch up to rk 0
      MPI_Barrier(retrieveGlobalComm());
#endif

      // allocate an OpComp as handle
      OpComp<float> op;

      // if mixed case, reset refr fields to pert fields so
      // IWaveLOVOp treats them all as extended - but need to
      // do this only for deriv>0 case
      if (physbg && (deriv>0)) {
	if (ps_slcstring(*pars,"csq",pert.c_str())) {
	  RVLException e;
	  e<<"ERROR: acd, deriv>0\n"; 
	  e<<"  failed to add physbg=1 to par table\n";
	  throw e;
	}
      }

      IWaveLOVOp iwop(*pars,stream);
      op.setNext(iwop);

      // possible taper operator
      std::string taperpars = valparse<std::string>(*pars,"taperpars","");
      std::shared_ptr<RVL::FunctionObject> foptr;
      std::shared_ptr<RVL::Operator<float> > tapop;
#ifdef IWAVE_USE_MPI
      std::shared_ptr<RVL::MPISerialFunctionObject<float> > mpifoptr;
#endif
      if (taperpars.size()>0) {
	std::shared_ptr<RVL::FunctionObject> aptr(new SEGYTaper(taperpars));
	foptr = aptr; 
#ifdef IWAVE_USE_MPI
	std::shared_ptr<RVL::MPISerialFunctionObject<float> > 
	  bptr(new  MPISerialFunctionObject<float>(*foptr));
	mpifoptr = bptr;
	LinearOpFO<float> taplfo(iwop.getRange(),iwop.getRange(),*mpifoptr,*mpifoptr);
#else
	LinearOpFO<float> taplfo(iwop.getRange(),iwop.getRange(),*foptr,*foptr);
#endif
	tapop = RVL::Operator<float>::clonePtr(taplfo);
      }
      /*      
#ifdef IWAVE_USE_MPI
      int rk = retrieveGlobalRank();
      if (rk==0) {
	// cerr<<"rk="<<rk<<endl;
#endif
	op.write(cerr);
#ifdef IWAVE_USE_MPI
      }
#endif
      */

      // link to data objects common to all simulations

      // model vector - note that if model has been extended, then
      // the file accessed here is the result of the extension, and
      // the domain of the op is the extended model space
      AssignFilename m0fn(refr);
      AssignFilename m1fn(src);
      
      // data vector
      Vector<float> d(op.getRange());
      AssignFilename dfn(valparse<std::string>(*pars,"data"));
      Components<float> cd(d);
      cd[0].eval(dfn);

      // forward map
      if (deriv==0) {
	Vector<float> m(iwop.getDomain());
	Components<float> cm(m);
	cm[0].eval(m0fn);
	cm[1].eval(m1fn);
	RVL::CompLinearOp<float> clop;
	LinearRestrictOp<float> rlop(iwop,cm[0]);
	clop.setNext(rlop);
	if (taperpars.size()>0) {
	  LinearOp<float> const & tapl = dynamic_cast<LinearOp<float> const &>(*tapop);
	  clop.setNext(tapl);
	}
	if (valparse<int>(*pars,"adjtest",0)) {
	  RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
	  RVL::AdjointTest(clop,rnd,cerr);
	}
	else {
	  if (adjoint==0) {
	    clop.applyOp(cm[1],d);
	  }
	  else {
	    clop.applyAdjOp(d,cm[1]);
	  }
	}
      }
      
      // derivatives
      else {

	// cerr<<"01\n";
	// entire construction from here on based on restriction to
	// 0th factor
	Vector<float> m(iwop.getDomain());
	Components<float> cm(m);
	//	cm[0].eval(m0fn); not initialized until after
	// branch on extension
	cm[1].eval(m1fn); // initializes source
	RestrictOp<float> rop(iwop,m,0);
	// here is a different LOVOp
	TangentMap<float> tm(rop);
	std::shared_ptr<RVL::Operator<float> > gextop;
	std::shared_ptr<RVL::Operator<float> > top;
	std::shared_ptr<RVL::LinOpValOp<float> > toplovo;
	std::shared_ptr<RVL::Operator<float> > tmp = RVL::Operator<float>::clonePtr(tm);
	std::shared_ptr<RVL::LinOpValOp<float> > tmplovo
	  = dynamic_pointer_cast<RVL::LinOpValOp<float> >(tmp);
	// mixed case: compose on first factor with tangent map
	if (physbg) {
	  	// cerr<<"02\n";
	  // if necessary create grid extension operator
#ifdef IWAVE_USE_MPI
	  MPIGridSpace refsp(refr,"notype",true);
#else
	  GridSpace refsp(refr,"notype",true);
#endif    
	  GridExtendOp extop(refsp,iwop.getProductDomain()[0]);
	  gextop = RVL::Operator<float>::clonePtr(extop);
	  // a lot of this would be easier with a construction like
	  // top = make_shared<CompLOVOp<float> >(*gextop,*tm);
	  CompLOVOp<float> clo(*gextop,*tmplovo);
	  top=RVL::Operator<float>::clonePtr(clo);
	  toplovo= dynamic_pointer_cast<RVL::LinOpValOp<float> >(top);
	}
	else {
	  // cerr<<"03\n";
	  toplovo=tmplovo;
	}
	// pert is csq_d1 if adjoint=0, csq_b1 if adjoint=1
	// cerr<<"04\n";
	AssignFilename m01fn(pert);
	// cerr<<"05\n";
	Vector<float> m01(toplovo->getDomain());
	// cerr<<"06\n";
	Components<float> cm01(m01);
	cm01[0].eval(m0fn);  // background csq - non-extended if phsybg
	cm01[1].eval(m01fn); // pert csq
	//	cerr<<"07\n";
	
	//	m01.write(cerr);
	// first derivatives
	if (deriv==1) {
	  RVL::CompLinearOp<float> clop;
	  LinearRestrictOp<float> lop(*toplovo,cm01[0]);
	  clop.setNext(lop);
	  if (taperpars.size()>0) {
	    LinearOp<float> const & tapl = dynamic_cast<LinearOp<float> const &>(*tapop);
	    clop.setNext(tapl);
	  }
	  int adjtest = valparse<int>(*pars,"adjtest",0);
	  //	  cerr<<"adjtest="<<adjtest<<endl;
	  if (adjtest) {
	    RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
	    RVL::AdjointTest(clop,rnd,cerr);
	  }
	  else {
	    if (adjoint==0) clop.applyOp(cm01[1],d);
	    else clop.applyAdjOp(d,cm01[1]);
	  }
	}
	
	// second deriv 
	else if (deriv==2) {
	  // cerr<<"11\n";
	  // another perturbation vect
	  Vector<float> m02(toplovo->getProductDomain()[0]);
	  // cerr<<"12\n";
	  // this is a PHYSICAL perturbation in either physical or mixed cases
	  pert2 = valparse<std::string>(*pars,pert2key);
	  // cerr<<"pert2="<<pert2<<"\n";
	  AssignFilename m02fn(pert2);
	  m02.eval(m02fn);
	  // cerr<<"2\n";
	  // fix 1st pert - so deriv is 2nd deriv of original map
	  RestrictOp<float> rop2(*toplovo,m01,0);
	  	  // cerr<<"3\n";
	  OperatorEvaluation<float> rop2eval(rop2,cm01[0]);
	  	  // cerr<<"4\n";
	  LinearOp<float> const & lop = rop2eval.getDeriv();
	  	  // cerr<<"6\n";
	  RVL::CompLinearOp<float> clop;
	  clop.setNext(lop);
	  	  // cerr<<"7\n";
	  if (taperpars.size()>0) {
	    LinearOp<float> const & tapl = dynamic_cast<LinearOp<float> const &>(*tapop);
	    clop.setNext(tapl);
	    	  // cerr<<"8\n";
	  }
	  if (valparse<int>(*pars,"adjtest",0)) {
	    RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
	    	  // cerr<<"9\n";
	    RVL::AdjointTest(clop,rnd,cerr);
	    	  // cerr<<"10\n";
	  }
	  else {
	    	  // cerr<<"11\n";
	    if (adjoint==0) clop.applyOp(m02,d);
	    else clop.applyAdjOp(d,m02);
	    	  // cerr<<"12\n";
	  }
	}
	
	else {
	  RVLException e;
	  e<<"ERROR: acd\n";
	  e<<"  currently only 0th, 1st, and 2nd derivatives implemented\n";
	  throw e;
	}
      }
#ifdef IWAVE_USE_MPI
    }
#endif

    ps_delete(&pars);
    iwave_fdestroy();
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
  
}
