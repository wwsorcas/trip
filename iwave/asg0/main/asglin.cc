#include "asg_defn.hh"
#include "lin_selfdoc.h"
#include "grid.h"
#include "iwop.hh"
#include "gridpp.hh"
#include "functions.hh"
#include "cgnealg.hh"
#include "adjtest.hh"
#include "segyops.hh"
#include "linop_apps.hh"

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
  {"source_v0",  D_V0, 1, 2},
  {"source_v1",  D_V1, 1, 2},
  {"source_v2",  D_V2, 1, 2},
  {"data_p",     D_P0, 0, 2},
  {"data_v0",    D_V0, 0, 2},
  {"data_v1",    D_V1, 0, 2},
  {"data_v2",    D_V2, 0, 2},
  {"",           0, 0, 0}
};

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
    TSOpt::IWaveEnvironment(argc, argv, 0, &pars, &stream);

    if (retrieveGlobalRank()==0 && argc<2) {
      pagedoc();
      exit(0);
    }

#ifdef IWAVE_USE_MPI
    if (retrieveGroupID() == MPI_UNDEFINED) {
      fprintf(stream,"NOTE: finalize MPI, cleanup, exit\n");
      fflush(stream);
    }
    else {
#endif

      // basic op
      TSOpt::IWaveLOVOp iwop(*pars,stream);

      // model vector - note that if model has been extended, then
      // the file accessed here is the result of the extension, and
      // the domain of the op is the extended model space
      RVL::Vector<float> m(iwop.getProductDomain());
      RVL::Components<float> cm(m);
      
      RVL::Vector<float> d(iwop.getIWaveRange());

      // load background model components
      TSOpt::IWaveSpace const & iwdom0
	= dynamic_cast<TSOpt::IWaveSpace const &>(iwop.getProductDomain()[0]);
      TSOpt::IWaveSpace const & iwdom1
	= dynamic_cast<TSOpt::IWaveSpace const &>(iwop.getProductDomain()[1]); 
      auto m0 = std::make_shared<RVL::Vector<float> >(iwdom0);
      auto m1 = std::make_shared<RVL::Vector<float> >(iwdom1);
      cerr<<"IWaveLoad on factor 0:\n";
      cerr<<"  keys:\n";
      for (int i=0; i< iwdom0.getKeys().size(); i++) 
	cerr<<"keys["<<i<<"] = "<<iwdom0.getKeys()[i]<<"\n";
      TSOpt::IWaveLoad(*pars, *m0, iwdom0.getKeys());
      cerr<<"IWaveLoad on factor 1:\n";
      cerr<<"  keys:\n";
      for (int i=0; i< iwdom1.getKeys().size(); i++) 
	cerr<<"keys["<<i<<"] = "<<iwdom1.getKeys()[i]<<"\n";
      TSOpt::IWaveLoad(*pars, *m1, iwdom1.getKeys());

      // select parameters to invert
      std::vector<std::string> est0;
      std::vector<std::string> est1;
      std::vector<std::string> estkeys0;
      std::vector<std::string> estkeys1;
      std::vector<int> estidx0;
      std::vector<int> estidx1;
      
      for (int i=0; i<iwdom0.getKeys().size(); i++) {
	std::string tmp = valparse<std::string>(*pars,
						iwdom0.getKeys()[i] + "_est",
						"");
	if (tmp.size()>0) {
	  est0.push_back(tmp);
	  estidx0.push_back(i);
	  estkeys0.push_back(iwdom0.getKeys()[i]+"_est");
	}
      }
      
      for (int i=0; i<iwdom1.getKeys().size(); i++) {
	std::string tmp = valparse<std::string>(*pars,
						iwdom1.getKeys()[i] + "_est",
						"");
	if (tmp.size()>0) {
	  est1.push_back(tmp);
	  estidx1.push_back(i);
	  estkeys1.push_back(iwdom1.getKeys()[i]+"_est");
	}
      }
      
      // not enough...
      if (est0.size() + est1.size() ==0) {
	RVL::RVLException e;
	e<<"ERROR: IWaveCGNE constructor - no ";
	e<<"  inversion targets specified (_est suffix)\n";
	throw e;
      }

      // too much...
      if ((est0.size() > 0) && (est1.size() > 0)) {
	RVL::RVLException e;
	e<<"ERROR: IWaveCGNE constructor";
	e<<"  inversion targets specified (_est suffix) for both\n";
	e<<"  linearized nonlinear and linear domain components\n";
	e<<"  Linear solver can't do both at once - that's a nonlinear problem\n";
	throw e;
      }	

      // parameter for CG
      float rtol=valparse<float>(*pars,"ResidualTol",100.0*numeric_limits<float>::epsilon());
      float nrtol=valparse<float>(*pars,"GradientTol",100.0*numeric_limits<float>::epsilon());
      int maxcount=valparse<int>(*pars,"MaxIter",10);
      float maxstep=valparse<float>(*pars,"MaxStep",numeric_limits<float>::max());
      
      // output stream 
      std::stringstream outstr;
      std::ostream * optr = NULL;
      std::string outfile = RVL::valparse<std::string>(*pars,"outfile","");
      if (retrieveRank()==0) {
	if (outfile.size()==0) optr = &cout;
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
      res<<"* Linear Inversion via CG Algorithm for Normal Eqns"<<endl;
      if (est0.size()>0) 
	res<<"* applied to Born approximation for model = "
	   <<iwop.getModel()<<endl;
      else
	res<<"* applied to estimate source params of model = "
	   <<iwop.getModel()<<endl;
      res<<"* max iterations       = "<<maxcount<<endl;
      res<<"* residual tolerance   = "<<rtol<<endl;
      res<<"* normal res tolerance = "<<nrtol<<endl;
      res<<"* trust radius         = "<<maxstep<<endl;
      res<<"*******************************************************"<<endl;

      // declare op for use in CG
      RVL::CompLinearOp<float> op;
      float rnorm;
      float nrnorm;
      float nrnorm0;
      float rnorm0;

      // source inversion
      if (est1.size() > 0) {
	RVL::InjectOp<float> injop(m1,estidx1);
	// would not be necessary if injop exposed its linear member
	RVL::Vector<float> fake(injop.getDomain());
	fake.zero();	
	RVL::OperatorEvaluation<float> injopeval(injop,fake);
	RVL::LinearRestrictOp<float> rlop(iwop,*m0);
	op.setNext(injopeval.getDeriv());
	op.setNext(rlop);
	RVL::Vector<float> dm(injop.getDomain());
	TSOpt::IWaveLoad(*pars, dm, estkeys1);
	dm.zero();
	RVLUmin::CGNEAlg<float> alg(dm,op,d,
				    rnorm, nrnorm,
				    rtol, nrtol,
				    maxcount, maxstep,
				    res);
	nrnorm0=nrnorm;
	rnorm0=rnorm;
	alg.run();
      }    
      // born inversion
      else {
	RVL::InjectOp<float> injop(m0,estidx0);
	RVL::Vector<float> fake(injop.getDomain());
	fake.zero();	
	RVL::OperatorEvaluation<float> injopeval(injop,fake);
	cm[0].copy(*m0);
	cm[1].copy(*m1);
	RVL::RestrictOp<float> rlop(iwop,m,1);
	RVL::OperatorEvaluation<float> rlopeval(rlop,*m0);
	op.setNext(injopeval.getDeriv());
	op.setNext(rlopeval.getDeriv());
	RVL::Vector<float> dm(injop.getDomain());
	TSOpt::IWaveLoad(*pars, dm, estkeys0);
	dm.zero();
	RVLUmin::CGNEAlg<float> alg(dm,op,d,
				    rnorm, nrnorm,
				    rtol, nrtol,
				    maxcount, maxstep,
				    res);
	nrnorm0=nrnorm;
	rnorm0=rnorm;
	alg.run();
      }
      
      // display results
      res<<"\n ******* summary ********  "<<endl;
      res<<"initial residual norm      = "<<rnorm0<<endl;
      res<<"residual norm              = "<<rnorm<<endl;
      res<<"residual redn              = "<<rnorm/rnorm0<<endl;
      res<<"initial gradient norm      = "<<nrnorm0<<endl;
      res<<"gradient norm              = "<<nrnorm<<endl;
      res<<"gradient redn              = "<<nrnorm/nrnorm0<<endl;
      
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
    }
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
