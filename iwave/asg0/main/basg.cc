#include "iwop.hh"
#include "functions.hh"
#include "linop_apps.hh"
#include <par.h>
#include "segyops.hh"
#include "gridops.hh"
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

#ifdef IWAVE_USE_MPI
typedef TSOpt::MPISEGYSpace mySEGYSpace;
typedef TSOpt::MPIGridSpace myGridSpace;
#else
typedef TSOpt::SEGYSpace mySEGYSpace;
typedef TSOpt::GridSpace myGridSpace;  
#endif

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

      // must modify parameters to use extended bulk modulus as
      // reference in construction of IWaveOp

      int deriv = RVL::valparse<int>(*pars,"deriv",0);
      if (deriv==0) {
	RVL::RVLException e;
	e<<"Error: basg driver\n";
	e<<"  only first and second derivatives of underlying nonlinear\n";
	e<<"  modeling operator are permitted\n";
	throw e;
      }

      int adjoint = RVL::valparse<int>(*pars,"adjoint",0);
      std::string psymb="";
      if (adjoint==0) psymb = "_d1";
      else psymb = "_b1";
      std::string dbulkmodname = "bulkmod" + psymb;
      std::string dbulkmod = RVL::valparse<std::string>(*pars,dbulkmodname);
      // create modifed par
      PARARRAY * modpars = ps_new();
      ps_copy(&modpars,*pars);
      // add bulkmod = same file as pert
      ps_slcstring(*modpars,"bulkmod",dbulkmod.c_str());
      // basic op which now has extended bulkmod grid (assuming that perturbation had it)
      TSOpt::IWaveLOVOp iwop(*modpars,stream);

      // create input bulkmod space, load input bulkmod
      myGridSpace bulksp(RVL::valparse<string>(*pars,"bulkmod"),"notype",true);      
      RVL::Vector<float> bulk(bulksp);
      RVL::AssignFilename afbulk(RVL::valparse<string>(*pars,"bulkmod"));
      bulk.eval(afbulk);

      // 0.0 0 component of cm[0] is bulk
      // physical
      RVL::Vector<float> m(iwop.getProductDomain());
      RVL::Components<float> cm(m);
      RVL::Components<float> cm0(cm[0]);
      TSOpt::GridExtendOp gop(bulksp,cm0[0].getSpace());
      gop.applyOp(bulk,cm0[0]);

      // 0.1 1 component of cm[0] is buoy assumed physical
      RVL::AssignFilename afbuoy(RVL::valparse<std::string>(*pars,"buoyancy"));
      cm0[1].eval(afbuoy);

      // perturbation vectors
      RVL::Vector<float> pm(iwop.getNonLinDomain());
      RVL::Components<float> cpm(pm);
      
      // 1.0 0 component of cpm[0] is bulk (possibly extended)
      RVL::AssignFilename afdbulk(RVL::valparse<std::string>(*pars,"bulkmod"+psymb));
      cpm[0].eval(afdbulk);

      // 1.1 1 component of cpm[0] is buoy
      RVL::AssignFilename afdbuoy(RVL::valparse<std::string>(*pars,"buoyancy"+psymb));
      cpm[1].eval(afdbuoy);

      // load the rest
      RVL::Components<float> cm1(cm[1]);
      RVL::Vector<float> d(iwop.getIWaveRange());      
      RVL::Components<float> dm(d);
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

      if (deriv==1) {

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
