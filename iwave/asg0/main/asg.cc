#include "iwop.hh"
#include "functions.hh"
#include "linop_apps.hh"
#include <par.h>
#include "segyops.hh"
#include "adjtest.hh"
#include "asg_defn.hh"
#include "sim_selfdoc.h"
#include "parser.h"
#include <omp.h>

enum{
  D_BULK=0,
  D_BUOY=1,
  D_P0=2,
  D_P1=3,
  D_P2=4,
  D_V0=5,
  D_V1=6,
  D_V2=7
};

IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"bulkmod",    D_BULK, 1, 1},
  {"buoyancy",   D_BUOY, 1, 1},
  {"source_p",   D_P0, 1, 2},
  {"data_p",     D_P0, 0, 2},
  {"data_v0",    D_V0, 0, 2},
  {"data_v1",    D_V1, 0, 2},
  {"data_v2",    D_V2, 0, 2},
  {"",           0, 0, 0}
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

#ifdef _OPENMP
    omp_set_num_threads(RVL::valparse<int>(*pars,"num_threads",1));
    //    std::cout << "Number of OMP threads in use: "
    //	      << omp_get_max_threads() << std::endl;
    if (retrieveGlobalRank()==0) {
      fprintf(stream,"Number of OMP threads in use: %d\n",omp_get_max_threads());
      fflush(stream);
    }
#endif    
    
#ifdef IWAVE_USE_MPI
    fprintf(stream,"rk=%d groupID=%d\n",retrieveGlobalRank(),retrieveGroupID());
    if (retrieveGroupID() == MPI_UNDEFINED) {
      fprintf(stream,"NOTE: idle process - finalize MPI, cleanup, exit\n");
      fflush(stream);
    }
    else {
#endif

      // basic op
      TSOpt::IWaveLOVOp iwop(*pars,stream);

      // link to data objects common to all simulations

      // model vector - note that if model has been extended, then
      // the file accessed here is the result of the extension, and
      // the domain of the op is the extended model space
      RVL::Vector<float> m(iwop.getProductDomain());
      RVL::Vector<float> d(iwop.getIWaveRange());
      RVL::Components<float> cm(m);
      RVL::Components<float> cm0(cm[0]);
      RVL::Components<float> cm1(cm[1]);
      RVL::Components<float> dm(d);
      for (int j=0; j< iwop.getNonLinDomain().getSize(); j++) {
	std::string fn =
	  RVL::valparse<std::string>(*pars,
				iwop.getNonLinDomain().getKeys()[j]);
	RVL::AssignFilename af(fn);
	cm0[j].eval(af);
      }
      for (int j=0; j<iwop.getLinDomain().getSize(); j++) {
	std::string fn =
	  RVL::valparse<std::string>(*pars,
				iwop.getLinDomain().getKeys()[j]);
	RVL::AssignFilename af(fn);
	cm1[j].eval(af);
      }
      for (int j=0; j< iwop.getIWaveRange().getSize(); j++) {
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
	  std::stringstream stuff;
	  RVL::AdjointTest(rlop,rnd,stuff);
	  fprintf(stream,stuff.str().c_str());
	}
	else {
	  if (adjoint==0) {
	    rlop.applyOp(cm[1],d);
	  }
	  else {
	    rlop.applyAdjOp(d,cm[1]);
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
	  std::stringstream stuff;
	  RVL::AdjointTest(ropeval.getDeriv(),rnd,stuff);
	  fprintf(stream,stuff.str().c_str());	  
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
      MPI_Barrier(MPI_COMM_WORLD);
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
