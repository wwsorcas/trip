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
#include "gridtrace.hh"
#ifdef IWAVE_USE_MPI
#include "mpiserialfo.hh"
#endif

enum{
    
  D_BULK=0,
  D_BUOY=1,
  D_P0=2,
  D_P1=3,
  D_P2=4,
  D_V0=5,
  D_V1=6,
  D_V2=7,
  D_PS=8    
};

IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"bulkmod",    D_BULK, 1, 1},
  {"buoyancy",   D_BUOY, 1, 1},
  {"source_p",   D_PS, 1, 2},
  {"data_p",     D_P0, 0, 2},
  {"data_v0",    D_V0, 0, 2},
  {"data_v1",    D_V1, 0, 2},
  {"data_v2",    D_V2, 0, 2},
  {"movie_p",    D_P0, 0, 0},
  {"",           0, 0, 0}
};

int xargc;
char **xargv;

namespace TSOpt {

  class readGridFO: public RVL::UnaryLocalFunctionObjectConstEval<float> {

  private:
    RVL::LocalDataContainer<float> & v;

    readGridFO();
    readGridFO(readGridFO const &);
    
  public:

    readGridFO(RVL::LocalDataContainer<float> & _v)
      : v(_v) {}
    
    void operator()(RVL::LocalDataContainer<float> const & x) {

      if (x.getSize() != v.getSize()) {
	RVL::RVLException e;
	e<<"Error: readGridFO::operator()\n";
	e<<"  input, stored vector sizes differ\n";
	throw e;
      }

      for (size_t i = 0; i < v.getSize(); i++ ) 
	v.getData()[i]=x.getData()[i];

    }

    std::string getName() const { std::string tmp = "readGridFO"; return tmp; }
  };

  class writeGridFO: public RVL::UnaryLocalFunctionObject<float> {

  private:
    RVL::LocalDataContainer<float> & v;
    float fac;

    writeGridFO();
    writeGridFO(readGridFO const &);
    
  public:

    writeGridFO(RVL::LocalDataContainer<float> & _v,
	       float _fac)
      : v(_v), fac(_fac) {}
    
    void operator()(RVL::LocalDataContainer<float> & x) {

      if (x.getSize() != v.getSize()) {
	RVL::RVLException e;
	e<<"Error: writeGridFO::operator()\n";
	e<<"  input, stored vector sizes differ\n";
	throw e;
      }

      for (size_t i = 0; i < v.getSize(); i++ ) 
	x.getData()[i]=fac*v.getData()[i];

    }

    std::string getName() const { std::string tmp = "writeGridFO"; return tmp; }
  };

  class GridTraceOp: public RVL::LinearOp<float> {

  private:

#ifdef IWAVE_USE_MPI
    MPIGridSpace const & dom;
#else
    GridSpace const & dom;
#endif
    IWaveSpace const & rng;
    float fac;
    
  protected:

    void apply(RVL::Vector<float> const & x,
	       RVL::Vector<float> & y) const {
      try {
	// sanity - no need, handled by construction and Op methods
	// GridSpace const & gsp = dynamic_cast<GridSpace const &>(x.getSpace());
	// SEGYSpace const & tsp = dynamic_cast<GridSpace const &>(y.getSpace());

	RVL::RnArray<float> buf(get_datasize_grid(dom.getGrid()));
	y.zero();
	/*
#ifdef IWAVE_USE_MPI
	readGridFO gf(buf);
	MPISerialFunctionObject<float> mpigf(gf);
	x.eval(mpigf);
	writeTraceSampleFO zs(buf);
	MPISerialFunctionObject<float> mpizs(zs);
	Components<float> cy(y);
	cy[0].eval(mpizs);
#else
	*/
	readGridFO gf(buf);
	x.eval(gf);
	writeTraceSampleFO zs(buf);
	Components<float> cy(y);
	cy[0].eval(zs);
	/*
#endif
	*/
      }
      catch (RVL::RVLException & e) {
	e<<"\ncalled from GridTraceOp::apply\n";
	throw e;
      }
    }

    void applyAdj(RVL::Vector<float> const & x,
		  RVL::Vector<float> & y) const {
      try {

	RVL::RnArray<float> buf(get_datasize_grid(dom.getGrid()));
	y.zero();
	readTraceSampleFO zs(buf);
	Components<float> cx(x);
	cx[0].eval(zs);
	writeGridFO gf(buf,fac);
	y.eval(gf);
      }
      catch (RVL::RVLException & e) {
	e<<"\ncalled from GridTraceOp::applyAdj\n";
	throw e;
      }
    }

    RVL::LinearOp<float> * clone() const { return new GridTraceOp(*this); }

  public:

#ifdef IWAVE_USE_MPI
    GridTraceOp(MPIGridSpace const & _dom,
#else
    GridTraceOp(GridSpace const & _dom,
#endif
		IWaveSpace const & _rng)
      : dom(_dom), rng(_rng) {
      // sanity checks
      if (RVL::ProtectedDivision<float>(1.0f,get_cellvol_grid(dom.getGrid()),fac)) {
	RVL::RVLException e;
	e<<"Error: GridTraceOp constructor\n";
	e<<"  zerodivide by cell vol\n";
	throw e;
      }
      if (rng.getSize() > 1) {
	RVL::RVLException e;
	e<<"Error: GridTraceOp constructor\n";
	e<<"  for exploding reflector modeling, range should be single\n";
	e<<"  component, sampling p\n";
	throw e;
      }
#ifdef IWAVE_USE_MPI
      MPISEGYSpace const * tmp = NULL;
      tmp = dynamic_cast<MPISEGYSpace const *>(&(rng[0]));
#else
      SEGYSpace const * tmp = NULL;
      tmp = dynamic_cast<SEGYSpace const *>(&(rng[0]));
#endif
      if (!tmp) {
	RVL::RVLException e;
	e<<"Error: GridTraceOp constructor\n";
	e<<"  for exploding reflector modeling, single component of range\n";
	e<<"  IWaveSpace should be SEGYSpace\n";
	throw e;
      }
    }
    
    GridTraceOp(GridTraceOp const & gt)
      : dom(gt.dom), rng(gt.rng), fac(gt.fac) {}

    RVL::Space<float> const & getDomain() const { return dom; }
    RVL::Space<float> const & getRange() const { return rng; }

    ostream & write(ostream & str) const {
      str<<"GridTraceOp: map between gridded data and traces\n";
      str<<"  suitable for exploding reflector sources\n";
      str<<"  domain:\n";
      dom.write(str);
      str<<"  range:\n";
      rng.write(str);
      return str;
    }
  };
}

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

      // exploding reflector model as grid function 
      // workspace for source traces
      std::string tracebuf = "tracebuf.su";
      // create source trace space for use in defining IWaveLOVOp. Add
      // to parameter list under source_p key - becomes lin domain for
      // iwop, by construction consistent with explref space
      createExplReflTraces(RVL::valparse<std::string>(*pars,"explref"),
			   tracebuf,
			   RVL::valparse<std::string>(*pars,"data_p"),
			   RVL::valparse<std::string>(*pars,"CWPROOT"),
			   stream);
      
      ps_slcstring(*pars,"source_p",tracebuf.c_str());
		   
      // basic op
      TSOpt::IWaveLOVOp iwop(*pars,stream);

      // sanity check
      if ((iwop.getLinDomain().getSize() != 1) &&
	  (iwop.getLinDomain().getKeys()[0] != "source_p")) {
	RVLException e;
	e<<"Error: explref\n";
	e<<"  lin domain of operator should consists of only one\n";
	e<<"  component with key source_p\n";
	throw e;
      }

      // model vector - note that if model has been extended, then
      // the file accessed here is the result of the extension, and
      // the domain of the op is the extended model space
      RVL::Vector<float> m(iwop.getNonLinDomain());
      RVL::Components<float> cm(m);
      {
	RVL::AssignFilename af(RVL::valparse<std::string>(*pars,"bulkmod"));
	cm[0].eval(af);
      }
      {
	RVL::AssignFilename af(RVL::valparse<std::string>(*pars,"buoyancy"));
	cm[1].eval(af);
      }
      
      RVL::Vector<float> d(iwop.getIWaveRange());
      /*
      // data - first check on key
      if ((iwop.getIWaveRange().getSize() == 1) &&
	  (iwop.getIWaveRange().getKeys()[0] == "data_p")) {
	{
	  RVL::AssignFilename af(RVL::valparse<std::string>(*pars,"data_p"));
	  RVL::Components<float> cd(d);
	  cd[0].eval(af);
	}
      }
      else {
	RVLException e;
	e<<"Error: explref\n";
	e<<"  Range space not singleton consisting of SEGYSpace with key = data_p\n";
	throw e;
      }
      */
      RVL::Components<float> cd(d);
      for (int i=0; i<iwop.getIWaveRange().getSize(); i++) {
	{
	  cerr<<"data key "<<i<<" = "<<iwop.getIWaveRange().getKeys()[i]<<"\n";
	  std::string dname =
	    RVL::valparse<std::string>(*pars,
				       iwop.getIWaveRange().getKeys()[i]);
	  cerr<<"file = "<<dname<<"\n";
	  RVL::AssignFilename af(dname);
	  cd[i].eval(af);
	}
      }
	
#ifdef IWAVE_USE_MPI
      TSOpt::MPIGridSpace gsp(RVL::valparse<std::string>(*pars,"explref"),"data");
#else
      TSOpt::GridSpace gsp(RVL::valparse<std::string>(*pars,"explref"),"data");
#endif

      // grid to traces - dom, range consistent by construction
      TSOpt::GridTraceOp gtr(gsp,iwop.getLinDomain());

      // linear op from iwop
      RVL::LinearRestrictOp<float> rlop(iwop,m);      

      // compose
      RVL::CompLinearOp<float> op;
      op.setNext(gtr);
      op.setNext(rlop);

      // branch

      // ask for inversion output as usual, by explref_inv key
      std::string explref_inv = RVL::valparse<std::string>(*pars,"explref_inv","");

      // set up workspace for exploding reflector model, connect either
      // to inversion output file or model file
      RVL::Vector<float> r(gsp);
      if (explref_inv.size() > 0) {
	RVL::AssignFilename af(explref_inv);
	r.eval(af);
      }
      else {
	RVL::AssignFilename af(RVL::valparse<std::string>(*pars,"explref"));
	r.eval(af);
      }
      
      // ask for adjoint tests by flag
      // for GridTraceOp
      int adjtest0 = RVL::valparse<int>(*pars,"adjtest0",0);
      // for fwd map
      int adjtest1 = RVL::valparse<int>(*pars,"adjtest1",0);
      // for composite
      int adjtest = RVL::valparse<int>(*pars,"adjtest",0);
      
      // ask for forward or adjoint map by flag
      int adjoint = RVL::valparse<int>(*pars,"adjoint",0);
      
      // adjoint tests
      RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
      if (adjtest0) RVL::AdjointTest(gtr,rnd,cerr);
      else if (adjtest1) RVL::AdjointTest(rlop,rnd,cerr);
      else if (adjtest) {
     	RVL::AdjointTest(op,rnd,cerr);
	/*
	RVL::Vector<float> x1(gtr.getDomain());
	RVL::Vector<float> x2(gtr.getRange());
	RVL::Vector<float> x3(rlop.getRange());
	RVL::Vector<float> y3(rlop.getRange());
	RVL::Vector<float> y2(gtr.getRange());
	RVL::Vector<float> y1(gtr.getDomain());
	x1.eval(rnd);
	y3.eval(rnd);
	gtr.applyOp(x1,x2);
	rlop.applyOp(x2,x3);
	rlop.applyAdjOp(y3,y2);
	gtr.applyAdjOp(y2,y1);

	cerr<<"x1.y1="<<x1.inner(y1)<<" x2.y2="<<x2.inner(y2)<<" x3.y3="<<x3.inner(y3)<<endl;
	y2.eval(rnd);
	gtr.applyAdjOp(y2,y1);
	cerr<<"x1.y1="<<x1.inner(y1)<<" x2.y2="<<x2.inner(y2)<<endl;
	*/
      }
      // exploding reflector inversion
      else if (explref_inv.size() > 0) {
	// parameter for CG
	float rtol=RVL::valparse<float>(*pars,"ResidualTol",100.0*numeric_limits<float>::epsilon());
	float nrtol=RVL::valparse<float>(*pars,"GradientTol",100.0*numeric_limits<float>::epsilon());
	int maxcount=RVL::valparse<int>(*pars,"MaxIter",10);
	float maxstep=RVL::valparse<float>(*pars,"MaxStep",numeric_limits<float>::max());
	
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
	
	res<<endl
	   <<"*********************************************************"<<endl;
	res<<"* Exploding Reflector Linear Inversion via CGNE Algorithm "<<endl;
	res<<"* max iterations       = "<<maxcount<<endl;
	res<<"* residual tolerance   = "<<rtol<<endl;
	res<<"* normal res tolerance = "<<nrtol<<endl;
	res<<"* trust radius         = "<<maxstep<<endl;
	res<<"*********************************************************"<<endl;

	r.zero();
	
	float rnorm;
	float nrnorm;
	float nrnorm0;
	float rnorm0;

	RVLUmin::CGNEAlg<float> alg(r,op,d,
				    rnorm, nrnorm,
				    rtol, nrtol,
				    maxcount, maxstep,
				    res);
	nrnorm0=nrnorm;
	rnorm0=rnorm;
	alg.run();

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
      }

      // exploding reflector modeling
      else if (adjoint==0) {
        op.applyOp(r,d);
      }
      // post-stack RTM
      else {
	cerr<<"data norm = "<<d.norm()<<endl;
	op.applyAdjOp(d,r);
	cerr<<"image norm = "<<r.norm()<<endl;
      }

      // clean up
      //      unlink(tracebuf.c_str());
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
