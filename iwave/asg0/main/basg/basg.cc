#include "iwop.hh"
#include "functions.hh"
#include "linop_apps.hh"
#include <par.h>
#include "segyops.hh"
#include "adjtest.hh"
#include "basg_defn.hh"
#include "sim_selfdoc.h"
#include "parser.h"
#include <omp.h>

// field list
enum{
  D_BULK=0,
  D_BUOY=1,
  D_P0=2,
  D_P1=3,
  D_P2=4,
  D_V0=5,
  D_V1=6,
  D_V2=7,
  D_DBULK=8,
  D_DBUOY=9,
  D_DP0=10,
  D_DP1=11,
  D_DP2=12,
  D_DV0=13,
  D_DV1=14,
  D_DV2=15,
};

IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"bulkmod",    D_BULK,  1, 1},
  {"dbulkmod",   D_DBULK, 1, 2},  
  {"buoyancy",   D_BUOY,  1, 1},
  {"dbuoyancy",  D_DBUOY, 1, 2},  
  {"source_p",   D_P0,    1, 0},
  {"data_p",     D_DP0,   0, 2},
  {"data_v0",    D_DV0,   0, 2},
  {"data_v1",    D_DV1,   0, 2},
  {"data_v2",    D_DV2,   0, 2},
  {"",           0,       0, 0}
};

int main(int argc, char ** argv) {

  try {
        
#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);
#endif
        
    PARARRAY * pars = NULL;
    FILE * stream = NULL;
    TSOpt::IWaveEnvironment(argc, argv, 0, &pars, &stream);

    if (retrieveGlobalRank()==0 && argc<2) {
      pagedoc();
      exit(0);
    }

    int num_threads = RVL::valparse<int>(*pars,"num_threads",1);
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#else
    num_threads=1;
#endif
    
#ifdef _OPENMP
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

      // basic op
      TSOpt::IWaveLOVOp iwop(*pars,stream);

      // link to data objects common to all simulations

      // model vector -
      // bulk, buoy, dbuoy are physical
      // dbulk is (possibly) extended
      RVL::Vector<float> m(iwop.getProductDomain());
      RVL::Vector<float> d(iwop.getIWaveRange());
      RVL::Components<float> cm(m);
      // components of nonlinear params (bulk, buoy)
      RVL::Components<float> cm0(cm[0]);
      // components of linear params (dbulk, dbuoy)
      RVL::Components<float> cm1(cm[1]);
      // components of output
      RVL::Components<float> dm(d);
      for (int j=0; j< iwop.getNonLinDomain().getSize(); j++) {
#ifdef IWAVE_USE_MPI
	if (retrieveRank()==0) 
#endif
	  cerr<<"load nonlin input param "<<j<<" key = "<< iwop.getNonLinDomain().getKeys()[j]<<"\n";
	std::string fn =
	  RVL::valparse<std::string>(*pars,
				iwop.getNonLinDomain().getKeys()[j]);
	RVL::AssignFilename af(fn);
	cm0[j].eval(af);
      }
      for (int j=0; j<iwop.getLinDomain().getSize(); j++) {
#ifdef IWAVE_USE_MPI
	if (retrieveRank()==0) 
#endif
	  cerr<<"load lin input param "<<j<<" key = "<< iwop.getLinDomain().getKeys()[j]<<"\n";	
	std::string fn =
	  RVL::valparse<std::string>(*pars,
				iwop.getLinDomain().getKeys()[j]);
	RVL::AssignFilename af(fn);
	cm1[j].eval(af);
      }
      for (int j=0; j< iwop.getIWaveRange().getSize(); j++) {
#ifdef IWAVE_USE_MPI
	if (retrieveRank()==0) 
#endif
	  cerr<<"load output param "<<j<<" key = "<< iwop.getIWaveRange().getKeys()[j]<<"\n";	
	std::string fn =
	  RVL::valparse<std::string>(*pars,
				iwop.getIWaveRange().getKeys()[j]);
	RVL::AssignFilename af(fn);
	dm[j].eval(af);
      }

      int deriv = RVL::valparse<int>(*pars,"deriv",0);
      int adjoint = RVL::valparse<int>(*pars,"adjoint",0);
      
      // forward map
      if (deriv==0) {
	RVL::LinearRestrictOp<float> rlop(iwop,cm[0]);
	if (RVL::valparse<int>(*pars,"adjtest",0)) {
	  RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
	  RVL::AdjointTest(rlop,rnd,cerr);
	}
	else {
	  if (adjoint==0) {
	    rlop.applyOp(cm[1],d);
	  }
	  else {
	    rlop.applyAdjOp(d,cm[1]);
	    cerr<<"output adjoint norm = "<<cm[1].norm()<<endl;
	  }
	}
      }
      else if (deriv==1) {
	// extract perturbational input/output
	std::string psymb="";
	if (adjoint==0) psymb = "_d1";
	else psymb = "_b1";
	// pert vector in dom[0]
	RVL::Vector<float> pm(iwop.getNonLinDomain());
	RVL::Components<float> cpm(pm);
	for (int j=0; j< iwop.getNonLinDomain().getSize(); j++) {
	  std::string fn =
	    RVL::valparse<std::string>(*pars,
				       iwop.getNonLinDomain().getKeys()[j]
				       + psymb);
	  RVL::AssignFilename af(fn);
	  cpm[j].eval(af);
	}
	// restrict to first component, evaluate
	RVL::RestrictOp<float> rop(iwop,m,0);
	RVL::OperatorEvaluation<float> ropeval(rop,cm[0]);
	// optional adjoint test
	if (RVL::valparse<int>(*pars,"adjtest",0)) {
	  RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
	  RVL::AdjointTest(ropeval.getDeriv(),rnd,cerr);
	}
	// application of derivative
	else {
	  if (adjoint==0) {
	    ropeval.getDeriv().applyOp(pm,d);
	  }
	  else {
	    ropeval.getDeriv().applyAdjOp(d,pm);
	  }
	}	
      }
      else {
	RVL::RVLException e;
	e<<"Error: asg\n";
	e<<"  only deriv=0 and deriv=1 implemented\n";
	throw e;
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
  catch (RVL::RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
}
