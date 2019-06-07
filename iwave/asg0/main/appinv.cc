#include "asg_defn.hh"
#include "born.hh"
#include <par.h>
#include "appinv_selfdoc.h"
#include "adjtest.hh"

//#define IWAVE_VERBOSE
#undef IWAVE_VERBOSE

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
    if (retrieveGlobalRank()==0) {
      std::cout << "Number of OMP threads in use: "
		<< omp_get_max_threads() << std::endl;
    }
#endif    
    
#ifdef IWAVE_USE_MPI
    if (retrieveGroupID() == MPI_UNDEFINED) {
      fprintf(stream,"NOTE: finalize MPI, cleanup, exit\n");
      fflush(stream);
    }
    else {
#endif

      // for appinv, must alter par array to point to abs surf data at z=0
      // before creating iwop, since change of geometry is involved
      int deriv   = RVL::valparse<int>(*pars,"deriv",0);
      int adjoint = RVL::valparse<int>(*pars,"adjoint",0);
      int freechoice = RVL::valparse<int>(*pars,"freesurface");
      string prefix;
      if (freechoice) prefix="free";
      else prefix="absb";

      std::vector<std::string> bkeys;
      bkeys.push_back("bulkmod");
      TSOpt::BornIWaveLOVOp biwop(bkeys,*pars,stream);

#ifdef IWAVE_VERBOSE
      if (retrieveGlobalRank()==0) {
	//	CStream cstr(stream);
	biwop.write(cerr);
	//	cstr.flush();
      }
#endif
      
      RVL::Vector<float> m(biwop.getProductDomain());
      RVL::Vector<float> d(biwop.getRange());
      RVL::Components<float> cm(m);
      RVL::Components<float> cm0(cm[0]);
      RVL::Components<float> cm1(cm[1]);
      RVL::Components<float> cd(d);
      for (int j=0; j<bkeys.size(); j++) {
	std::string key = prefix + bkeys[j];
	std::string fn = RVL::valparse<std::string>(*pars, key);
	RVL::AssignFilename af(fn);
	cm0[j].eval(af);
#ifdef IWAVE_VERBOSE
	cerr<<"loaded "<<fn<<" on key "<<key<<" into BornIWOP ref component "<<j<<"\n";
#endif
      }

      for (int j=0; j<bkeys.size(); j++) {
	std::string key = prefix + bkeys[j]+"pert";
	std::string fn = RVL::valparse<std::string>(*pars, key);
	RVL::AssignFilename af(fn);
	cm1[j].eval(af);
#ifdef IWAVE_VERBOSE
	cerr<<"loaded "<<fn<<" on key "<<key<<" into BornIWOP pert component "<<j<<"\n";
#endif
      }
      
      for (int j=0; j< biwop.getRangeKeys().size(); j++) {
	std::string key = prefix + biwop.getRangeKeys()[j];
	std::string fn = RVL::valparse<std::string>(*pars, key);
	RVL::AssignFilename af(fn);
	cd[j].eval(af);
#ifdef IWAVE_VERBOSE	
	cerr<<"loaded "<<fn<<" on key "<<key<<" into BornIWOP range component "<<j<<"\n";
#endif
      }

#ifdef IWAVE_VERBOSE
      cerr<<"deriv = "<<deriv<<" adjoint = "<<adjoint<<"\n";
#endif
      if (deriv==0) {
	RVL::LinearRestrictOp<float> lrop(biwop,cm[0]);
	if (RVL::valparse<int>(*pars,"adjtest",0)) {
	  RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
	  std::stringstream stuff;
	  RVL::AdjointTest(lrop,rnd,stuff);
	  if (retrieveGlobalRank()==0) {
	    fprintf(stream,stuff.str().c_str());
	  }
	}	
	if (adjoint==0) {
#ifdef IWAVE_VERBOSE
	  cerr<<"forward Born map\n";
#endif
	  lrop.applyOp(cm[1],d);
	}
	else {
#ifdef IWAVE_VERBOSE	  
	  cerr<<"adjoint Born map\n";
#endif
	  lrop.applyAdjOp(d,cm[1]);
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
  catch (RVL::RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
}
