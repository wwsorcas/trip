#include "asg_defn.hh"
//#include "params.h"
#include "grid.h"
#include "iwop.hh"
#include "gridpp.hh"
#include "functions.hh"
#include "cgnealg.hh"
#include "segyops.hh"
#include "born.hh"

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

using RVL::parse;
using RVL::RVLException;
using RVL::Vector;
using RVL::Components;
using RVL::Operator;
using RVL::OperatorEvaluation;
using RVL::LinearOp;
using RVL::LinearOpFO;
using RVL::OpComp;
using RVL::SymmetricBilinearOp;
using RVL::AssignFilename;
using RVL::AssignParams;
using RVL::RVLRandomize;
using TSOpt::IWaveEnvironment;
using TSOpt::IWaveTree;
using TSOpt::IWaveSampler;
using TSOpt::IWaveSim;
using TSOpt::TASK_RELN;
using TSOpt::IOTask;
using TSOpt::IWaveOp;
using TSOpt::SEGYLinMute;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
using TSOpt::MPISEGYSpace;
typedef TSOpt::MPIGridSpace myGridSpace;
typedef TSOpt::MPISEGYSpace mySEGYSpace;
#else
using TSOpt::GridSpace;
using TSOpt::SEGYSpace;
typedef TSOpt::GridSpace myGridSpace;
typedef TSOpt::SEGYSpace mySEGYSpace;
#endif

using RVLUmin::CGNEAlg;

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
    IWaveEnvironment(argc, argv, 0, &pars, &stream);

#ifdef IWAVE_USE_MPI
    if (retrieveGroupID() == MPI_UNDEFINED) {
      fprintf(stream,"NOTE: finalize MPI, cleanup, exit\n");
      fflush(stream);
    }
    else {
#endif    
    
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
	std::string key = prefix + bkeys[j]+"inv";
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
    
      float rtol=valparse<float>(*pars,"ResidualTol",100.0*numeric_limits<float>::epsilon());
      float nrtol=valparse<float>(*pars,"GradientTol",100.0*numeric_limits<float>::epsilon());
      int maxcount=valparse<int>(*pars,"MaxIter",10);
      float maxstep=valparse<float>(*pars,"MaxStep",numeric_limits<float>::max());

      /* output stream */
      std::stringstream outstr;
    
      std::ostream * optr = NULL;
      std::string outfile = RVL::valparse<std::string>(*pars,"outfile","");
      if (retrieveRank()==0) {
	if (outfile.size()==0) optr = &cerr;
	else {
	  optr = new std::ofstream(outfile.c_str());
	  (*optr)<<scientific;	  
	}
      }
      else {
	optr = new std::stringstream;
      }
      ostream & res = *optr;
    
      res<<endl<<"*******************************************************"<<endl;
      res<<"* Acoustic Bulk Modulus Only Linearized Inversion via";
      res<<"* Conjugate Gradient Algorithm for Normal Eqns"<<endl;
      res<<"* max iterations       = "<<maxcount<<endl;
      res<<"* residual tolerance   = "<<rtol<<endl;
      res<<"* normal res tolerance = "<<nrtol<<endl;
      res<<"* trust radius         = "<<maxstep<<endl;
      res<<"*******************************************************"<<endl;
    
      /* create CGNE object */
      float rnorm;
      float nrnorm;
      RVL::LinearRestrictOp<float> lrop(biwop,cm[0]);
      cm[1].zero();
      CGNEAlg<float> alg(cm[1],lrop,d,
			 rnorm, nrnorm, rtol, nrtol, maxcount, maxstep, res);
    
      alg.run();
    
      std::string dataest = valparse<std::string>(*pars,"dataest","");
      std::string datares = valparse<std::string>(*pars,"datares","");
      if (dataest.size()>0) {
	Vector<float> est(lrop.getRange());
	AssignFilename estfn(dataest);
	est.eval(estfn);
	lrop.applyOp(cm[1],est);
	if (datares.size()>0) {
	  Vector<float> dres(lrop.getRange());
	  AssignFilename resfn(datares);
	  dres.eval(resfn);
	  dres.copy(d);
	  dres.linComb(-1.0f,est);
	} 
      }

      if ((outfile.size()>0) && (retrieveRank()==0) && (optr)) {
	std::ofstream * ofptr = NULL;
	ofptr = dynamic_cast<std::ofstream *>(optr);
	if (ofptr) {
	  ofptr->flush();
	  ofptr->close();
	}
	delete optr;
      }
    
#ifdef IWAVE_USE_MPI
      MPI_Finalize();
    }
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
