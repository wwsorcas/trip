#include "parser.h"
#ifdef IWAVE_USE_MPI
#include "mpisegypp.hh"
#else
#include "segypp.hh"
#endif
#include "segyops.hh"
#include "cgnealg.hh"
#include "adjtest.hh"

#ifdef IWAVE_USE_MPI
using TSOpt::MPISEGYSpace;
typedef TSOpt::MPISEGYSpace gsp;
#else
using TSOpt::SEGYSpace;
typedef TSOpt::SEGYSpace gsp;
#endif

int xargc;
char **xargv;

const char * sdoc[] = {
  "SEGYCGDecon: convolve, cross-correlate, or deconvolve a single trace signature ",
  "with/from a collection of traces. Specify su files containing source and target",
  "traces and signature trace. Files must exist on call, as they specify the domain",
  "and range of the convolution operator. Relation is",
  " ",
  "    target = source [convolve] signature (mode=fwd or mode=solve)",
  "    source = target [adjoint-convolve] signature (mode=adj)",
  " ",
  "required inputs:",
  "data     = name of su file containing target traces: input for mode=adj or",
  "           mode=solve, output for mode=fwd",
  "soln     = name of su file containing source traces: input for mode=fwd,",
  "           output for mode=adj or mode=solve",
  "kernel   = name of file containing signature trace", 
  "mode     = fwd -   convolve signature with soln, output data",
  "           adj -   adjoint-convolve (cross-correlate) signature with data, output soln",
  "           test -  apply adjoint test to convolution operator, with domain and range "
  "                   defined by soln and data files",
  "           solve - perform deconvolution via CG iteration",
  "additional inputs needed for mode=solve [default value]:",
  "ResidualTol      = stopping threshhold for residual norm [100*macheps]",
  "GradientTol      = stopping threshhold for gradient (normal residual) ",
  "                   norm [100*macheps]",
  "MaxIter          = max number of CG iterations permitted [10]",
  NULL};

int main(int argc, char ** argv) {

  try {

    if (retrieveGlobalRank()==0 && argc<2) {
      pagedoc();
      exit(0);
    }
    
#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);   
    storeGlobalComm(MPI_COMM_WORLD);
#endif

    PARARRAY * pars = ps_new();
    
    if (ps_createargs(pars,argc-1,&(argv[1]))) {
      RVL::RVLException e;
      e<<"ERROR: SEGYCGDecon from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++) 
	e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }
    // solves x*k=d for x given k and d
    std::string dname = RVL::valparse<string>(*pars,"data");
    std::string xname = RVL::valparse<string>(*pars,"soln");
    std::string kname = RVL::valparse<string>(*pars,"kernel");

    cerr<<"construct domain on "<<xname<<endl;
    gsp dom(xname,"notype"
#ifdef IWAVE_USE_MPI 
	   , retrieveGlobalComm()
#endif
	   );
    dom.write(cerr);
    cerr<<"construct range on "<<dname<<endl;
    gsp rng(dname,"notype"
#ifdef IWAVE_USE_MPI 
	   , retrieveGlobalComm()
#endif
	   );
    rng.write(cerr);
    TSOpt::SEGYTraceConvolveOp op(dom,rng,kname);
    cerr<<"dom ="<<xname<<endl;
    op.getDomain().write(cerr);
    cerr<<"rng ="<<dname<<endl;
    op.getRange().write(cerr);
    
    
    RVL::Vector<float> d(rng);
    RVL::AssignFilename daf(dname);
    d.eval(daf);

    RVL::Vector<float> x(dom);
    RVL::AssignFilename xaf(xname);
    x.eval(xaf);

    std::string mode = RVL::valparse<string>(*pars,"mode");

    if (mode == "fwd") op.applyOp(x,d);
    if (mode == "adj") op.applyAdjOp(d,x);
    if (mode == "test") {
      RVL::RVLRandomize<float> rnd(getpid(),-1.0,1.0);
      RVL::AdjointTest(op,rnd,cerr);
    }
    if (mode == "solve") {
      float rtol=RVL::valparse<float>(*pars,"ResidualTol",100.0*numeric_limits<float>::epsilon());
      float nrtol=RVL::valparse<float>(*pars,"GradientTol",100.0*numeric_limits<float>::epsilon());
      int maxcount=RVL::valparse<int>(*pars,"MaxIter",10);
      float maxstep=RVL::valparse<float>(*pars,"MaxStep",numeric_limits<float>::max());

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
      res<<"* Signature Deconvolution via Least Squares";
      res<<"* Conjugate Gradient Algorithm for Normal Eqns"<<endl;
      res<<"* max iterations       = "<<maxcount<<endl;
      res<<"* residual tolerance   = "<<rtol<<endl;
      res<<"* normal res tolerance = "<<nrtol<<endl;
      res<<"* trust radius         = "<<maxstep<<endl;
      res<<"*******************************************************"<<endl;
    
      /* create CGNE object */
      float rnorm;
      float nrnorm;
      x.zero();
      RVLUmin::CGNEAlg<float> alg(x,op,d,
				  rnorm, nrnorm, rtol, nrtol, maxcount, maxstep, res);
      //      float nrnorm0=nrnorm;
      //float rnorm0=rnorm;
      
      alg.run();
      
      // display results
      /*
      res<<"\n ******* summary ********  "<<endl;
      res<<"initial residual norm      = "<<rnorm0<<endl;
      res<<"residual norm              = "<<rnorm<<endl;
      res<<"residual redn              = "<<rnorm/rnorm0<<endl;
      res<<"initial gradient norm      = "<<nrnorm0<<endl;
      res<<"gradient norm              = "<<nrnorm<<endl;
      res<<"gradient redn              = "<<nrnorm/nrnorm0<<endl;
      */
      std::string dataest = valparse<std::string>(*pars,"dataest","");
      std::string datares = valparse<std::string>(*pars,"datares","");
      if (dataest.size()>0) {
	RVL::Vector<float> est(op.getRange());
	RVL::AssignFilename estfn(dataest);
	est.eval(estfn);
	op.applyOp(x,est);
	if (datares.size()>0) {
	  RVL::Vector<float> dres(op.getRange());
	  RVL::AssignFilename resfn(datares);
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
    }
    
    ps_delete(&pars);
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


