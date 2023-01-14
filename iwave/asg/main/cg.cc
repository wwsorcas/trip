#include "asg_defn.hh"
#include "lin_selfdoc.h"
#include "grid.h"
#include "iwop.hh"
#include "gridpp.hh"
#include "functions.hh"
#include "cgnealg.hh"
#include "segyops.hh"

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
#else
using TSOpt::GridSpace;
using TSOpt::SEGYSpace;
#endif

using RVLUmin::CGNEAlg;

int xargc;
char **xargv;

int main(int argc, char ** argv) {

  try {

    if (retrieveGlobalRank()==0 && argc<2) {
      pagedoc();
      exit(0);
    }
    
#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
#endif

    PARARRAY * pars = NULL;
    FILE * stream = NULL;
    IWaveEnvironment(argc, argv, 0, &pars, &stream);
    
    // the Op
    TSOpt::IWaveLOVOp iwop(*pars,stream);
    /*      
    SEGYLinMute mute(valparse<float>(*pars,"mute_slope",0.0f),
                     valparse<float>(*pars,"mute_zotime",0.0f),
                     valparse<float>(*pars,"mute_width",0.0f));
      
    LinearOpFO<float> muteop(iwop.getRange(),iwop.getRange(),mute,mute);
    */
    Vector<ireal> m(iwop.getProductDomain());
    Components<ireal> cm(m);
    // 0 component is bulkmod buoyancy
    Components<ireal> cm0(cm[0]);
    // 1 component is source_p
    Components<ireal> cm1(cm[1]);
    
    Vector<ireal> dm(iwop.getProductDomain()[0]);
    Components<ireal> cdm(dm);

    Vector<ireal> dd(iwop.getRange());
    Vector<ireal> mdd(iwop.getRange());
    
    AssignFilename mfn0(valparse<std::string>(*pars,"bulkmod"));
    cm0[0].eval(mfn0);
    AssignFilename mfn1(valparse<std::string>(*pars,"buoyancy"));
    cm0[1].eval(mfn1);

    AssignFilename dmfn0(valparse<std::string>(*pars,"bulkmod_est"));
    cdm[0].eval(dmfn0);
    AssignFilename dmfn1(valparse<std::string>(*pars,"buoyancy_est"));
    cdm[1].eval(dmfn1);
    dm.zero();

    AssignFilename ffn(valparse<std::string>(*pars,"source_p"));
    cm1[0].eval(ffn);
    
    AssignFilename ddfn(valparse<std::string>(*pars,"data_p"));
    dd.eval(ddfn);
    /*
    std::string mddnm = valparse<std::string>(*pars,"datamut","");
    if (mddnm.size()>0) {
      AssignFilename mddfn(mddnm);
      mdd.eval(mddfn);
    }
    muteop.applyOp(dd,mdd);
    */
    
    RVL::RestrictOp<ireal> op(iwop,m,0);
    //    OpComp<float> op(rwop,muteop);
    
    float rtol=valparse<float>(*pars,"ResidualTol",100.0*numeric_limits<float>::epsilon());
    float nrtol=valparse<float>(*pars,"GradientTol",100.0*numeric_limits<float>::epsilon());
    int maxcount=valparse<int>(*pars,"MaxIter",10);
    float maxstep=valparse<float>(*pars,"MaxStep",numeric_limits<float>::max());

    RVLRandomize<float> rnd(getpid(),-1.0,1.0);

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
    res<<"* Acoustic Variable Density Linearized Inversion via";
    res<<"* Conjugate Gradient Algorithm for Normal Eqns"<<endl;
    res<<"* max iterations       = "<<maxcount<<endl;
    res<<"* residual tolerance   = "<<rtol<<endl;
    res<<"* normal res tolerance = "<<nrtol<<endl;
    res<<"* trust radius         = "<<maxstep<<endl;
    res<<"*******************************************************"<<endl;
    
    /* create CGNE object */
    float rnorm;
    float nrnorm;
    OperatorEvaluation<ireal> opeval(op,cm[0]);
    CGNEAlg<float> alg(dm,opeval.getDeriv(),dd,
		       rnorm, nrnorm, rtol, nrtol, maxcount, maxstep, res);
    float nrnorm0=nrnorm;
    float rnorm0=rnorm;
    
    alg.run();
    
    // display results
    res<<"\n ******* summary ********  "<<endl;
    res<<"initial residual norm      = "<<rnorm0<<endl;
    res<<"residual norm              = "<<rnorm<<endl;
    res<<"residual redn              = "<<rnorm/rnorm0<<endl;
    res<<"initial gradient norm      = "<<nrnorm0<<endl;
    res<<"gradient norm              = "<<nrnorm<<endl;
    res<<"gradient redn              = "<<nrnorm/nrnorm0<<endl;
    
    std::string dataest = valparse<std::string>(*pars,"dataest","");
    std::string datares = valparse<std::string>(*pars,"datares","");
    if (dataest.size()>0) {
      Vector<float> est(op.getRange());
      AssignFilename estfn(dataest);
      est.eval(estfn);
      opeval.getDeriv().applyOp(dm,est);
      if (datares.size()>0) {
	Vector<float> dres(op.getRange());
	AssignFilename resfn(datares);
	dres.eval(resfn);
	dres.copy(mdd);
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
#endif
    exit(0);
  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
  
}
