#include "iwopt.hh"

//#define VERBOSE

namespace TSOpt {

  /** builder class for preconditioners, functions as container
      to pass build instructions through IWaveOpApply interface.
      implemented base class with virtual build function, returning
      a null shared_ptr.
  */
  
  
  void IWaveLoad(PARARRAY pars,
		 RVL::Vector<float> & v,
		 std::vector<std::string> keys) {
    try {
      RVL::Components<ireal> cv(v);
      if (cv.getSize() != keys.size()) {
	RVL::RVLException e;
	e<<"ERROR: IWaveLoad\n";
	e<<"  number of vector components differs from number of keys\n";
	e<<"  aaaaaiiiiiieeeeeeeeeeeeeee!!!!!!!!!!!!!\n";
	e<<"  keys:\n";
	for (int i=0; i<keys.size(); i++) e<<"keys["<<i<<"] = "<<keys[i]<<"\n";
	e<<"  RVL::Vector:\n";
	v.RVL::Writeable::write(e);
	throw e;
      }	 
      for (int i=0; i<keys.size(); i++) {
	std::string fs = valparse<std::string>(pars,keys[i]);
#ifdef VERBOSE
	cerr<<"IWAVE load: "<<endl;
	cerr<<"  assigning file = "<<fs<<endl;
	cerr<<"  to component index = "<<i<<" key = "<<keys[i]<<endl;
#endif
	RVL::AssignFilename fn(fs);
	cv[i].eval(fn);
      }
    }
    catch (RVL::RVLException & e) {
      e<<"\ncalled from IWaveLoad\n";
      throw e;
    }
  }

  IWaveFWIOp::IWaveFWIOp(PARARRAY pars, FILE * stream)
    : iwop(pars, stream) {
          
    try {

#ifdef VERBOSE
      cerr<<"begin IWaveFWIOp\n";
      cerr<<"load bg model\n";
#endif
      
      IWaveSpace const & iwdom = dynamic_cast<IWaveSpace const &>(iwop.getDomain());
      m0 = RVL::Vector<float>::newPtr(iwop.getDomain());
      IWaveLoad(pars, *m0,iwdom.getKeys());
      RVL::Components<float> cm0(*m0);


      /* select parameters to invert 
       */
      for (int i=0; i<iwdom.getKeys().size(); i++) {
	std::string tmp = valparse<std::string>(pars,
						iwdom.getKeys()[i] + "_est",
						"");
	if (tmp.size()>0) {
	  est.push_back(tmp);
	  estidx.push_back(i);
	  estkeys.push_back(iwdom.getKeys()[i]+"_est");
	}
      }
      
      /* detect null inversion */
      if (est.size()==0) {
	RVL::RVLException e;
	e<<"ERROR: FWI constructor - no ";
	e<<"  inversion targets specified (_est suffix)\n";
	throw e;
      }

      /* create domain, range spaces for diagonal window op */
      // might as well use clone since at this point nothing much
      // is allocated
#ifdef VERBOSE
      cerr<<"create domain, intermed spaces\n";
#endif
      {
#ifdef VERBOSE
	cerr<<"StdProductSpaces of size "<<est.size()<<endl;
#endif
	RVL::StdProductSpace<float> _dom(est.size());
	RVL::StdProductSpace<float> _med(est.size());
	dom = RVL::Space<float>::clonePtr(_dom);
	med = RVL::Space<float>::clonePtr(_med);
      }
      RVL::StdProductSpace<float> & domref = dynamic_cast<RVL::StdProductSpace<float> &>(*dom);
      RVL::StdProductSpace<float> & medref = dynamic_cast<RVL::StdProductSpace<float> &>(*med);
      for (int i=0; i<est.size(); i++) {
#ifdef IWAVE_USE_MPI
	TSOpt::MPIGridSpace dsp(est[i],"notype",true);
#else
	TSOpt::GridSpace dsp(est[i],"notype",true);
#endif
#ifdef VERBOSE
	cerr<<"set space factor "<<i<<" constructed on "<<est[i]<<endl;
	cerr<<"space index "<<estidx[i]<<endl;
#endif
	domref.set(dsp,i);
	medref.set(iwdom[estidx[i]],i);
      }

#ifdef VERBOSE
      cerr<<"upper, lower bounds\n";
#endif
      
      /* upper and lower bounds - note that these are required only
	 for parameters to be inverted */
      ub = RVL::Vector<float>::newPtr(*med);
      lb = RVL::Vector<float>::newPtr(*med);
      RVL::Components<ireal> cub(*ub);
      RVL::Components<ireal> clb(*lb);
      for (int i=0; i<est.size(); i++) {
#ifdef VERBOSE
	cerr<<"upper bound component "<<i<<"\n";
#endif
	std::string ubn = valparse<std::string>(pars,
						iwdom.getKeys()[estidx[i]] + "_ub",
						"");
        if (ubn.size() > 0) {
	  RVL::AssignFilename fn(ubn);
	  cub[i].eval(fn);
	}
	else {
	  try {
	    float val = valparse<float>(pars,
					iwdom.getKeys()[estidx[i]] + "_max");
	    RVL::RVLAssignConst<float> ac(val);
	    cub[i].eval(ac);
#ifdef VERBOSE
	    cerr<<"assigned upper bound for param "<<i<<" = "<<val<<"\n";
#endif
	  }
	  catch (RVL::RVLException & e) {
	    e<<"ERROR: FWI\n";
	    e<<"  neither upper envelope vector nor scalar max given\n";
	    e<<"  for parameter "<<iwdom.getKeys()[estidx[i]]<<"\n";
	    throw e;
	  }
	}
	  
#ifdef VERBOSE
	cerr<<"lower bound component "<<i<<"\n";
#endif
	std::string lbn = valparse<std::string>(pars,
						iwdom.getKeys()[estidx[i]] + "_lb",
						"");
        if (lbn.size() > 0) {
	  RVL::AssignFilename fn(lbn);
	  clb[i].eval(fn);
	}
	else {
	  try {
	    float val = valparse<float>(pars,
					iwdom.getKeys()[estidx[i]] + "_min");
	    RVL::RVLAssignConst<float> ac(val);
	    clb[i].eval(ac);
#ifdef VERBOSE
	    cerr<<"assigned lower bound for param "<<i<<" = "<<val<<"\n";
#endif
	  }
	  catch (RVL::RVLException & e) {
	    e<<"ERROR: FWI\n";
	    e<<"  neither lower envelope vector nor scalar min given\n";
	    e<<"  for parameter "<<iwdom.getKeys()[estidx[i]]<<"\n";
	    throw e;
	  }
	}
      }
#ifdef VERBOSE
      cerr<<"logistic repn\n";
#endif

      /* logistic representation of background */
      lm0 = RVL::Vector<float>::newPtr(*med);
      RVL::Components<float> clm0(*lm0);
      RVL::RVLVectorLogisticInverse<float> vli;
      for (int i=0; i<est.size(); i++) {
	std::vector<RVL::Vector<float> const *> vv(3);
	vv[0]=&(cm0[estidx[i]]);
	vv[1]=&(clb[i]);
	vv[2]=&(cub[i]);
	clm0[i].eval(vli,vv);
#ifdef VERBOSE
#ifdef IWAVE_USE_MPI
	RVL::RVLMin<float> svmin;
	RVL::RVLMax<float> svmax;
	RVL::MPISerialFunctionObjectRedn<float,float> vmin(svmin);
	RVL::MPISerialFunctionObjectRedn<float,float> vmax(svmax);
	clm0[i].eval(vmin);
	clm0[i].eval(vmax);
	cerr<<"*** background model comp "<<i<<" min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
#else
	RVL::RVLMin<float> vmin;
	RVL::RVLMax<float> vmax;
	clm0[i].eval(vmin);
	clm0[i].eval(vmax);
	cerr<<"*** background model comp "<<i<<" min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
#endif
#endif
      }

#ifdef VERBOSE
      cerr<<"diag window op\n";
#endif
      // create diag window op with domain and range
      // tapers same for all components
      RPNT sw;
      sw[0]=valparse<float>(pars,"taper1",0.0f);
      sw[1]=valparse<float>(pars,"taper2",0.0f);
      sw[2]=valparse<float>(pars,"taper3",0.0f);
      {
	RVL::DiagOp<float> _dwop(est.size());
	dwop = RVL::Operator<float>::clonePtr(_dwop);
      }
      RVL::DiagOp<float> & dwopref = dynamic_cast<RVL::DiagOp<float> &>(*dwop);
      for (int i=0; i<est.size(); i++) {
#ifdef IWAVE_USE_MPI
	TSOpt::MPIGridSpace dsp(est[i],"notype",true);
#else
	TSOpt::GridSpace dsp(est[i],"notype",true);
#endif
	TSOpt::GridWindowOp wop(domref[i],clm0[i],sw);
	dwopref.set(wop,i);
      }

#ifdef VERBOSE
      cerr<<"logistic transform\n";
#endif
      
      /* Logistic transform - acts on range of window op */
      RVL::VectorLogisticOp<float> _dlop(*lb,*ub);
      dlop = RVL::Operator<float>::clonePtr(_dlop);
      
#ifdef VERBOSE
      cerr<<"injection op\n";
#endif
      /* inject range of window op into domain of IWaveOp -
	 alter components in range, leave other at bg values */
      {
	RVL::InjectOp<float> _injop(m0,estidx);
	injop = RVL::Operator<float>::clonePtr(_injop);
#ifdef VERBOSE
	injop->write(cerr);
#endif
      }

#ifdef VERBOSE
      cerr<<"opcomps\n";
#endif
      
      /* build the ops */
      {
	// subgrid logistic params to grid logistic params
	// to grid model params, est model comps only
#ifdef VERBOSE
	cerr<<"transop = diagonal window op then diagonal logistic op\n";
#endif
	RVL::OpComp<float> _transop(*dwop,*dlop);
	// est model comps to all model comps to data
#ifdef VERBOSE
	cerr<<"modelop = inject op then iwop\n";
#endif
	RVL::OpComp<float> _modelop(*injop,iwop);
	// subgrid est logistic params to data
	RVL::OpComp<float> _op(_transop,_modelop);
	transop = RVL::Operator<float>::clonePtr(_transop);
	modelop = RVL::Operator<float>::clonePtr(_modelop);
	op = RVL::Operator<float>::clonePtr(_op);
      }

#ifdef VERBOSE
      cerr<<"end constructor\n";
#endif
      
    }
    catch (RVL::RVLException & e) {
      e<<"\ncalled from IWaveFWIOp constructor\n";
      throw e;
    }
  }

  IWaveFWIOp::IWaveFWIOp(IWaveFWIOp const & x) 
    : iwop(x.iwop), m0(x.m0),
      est(x.est), estidx(x.estidx),
      dom(x.dom), med(x.med),
      ub(x.ub), lb(x.lb), lm0(x.lm0),
      dwop(x.dwop), dlop(x.dlop),
      injop(x.injop), transop(x.transop),
      modelop(x.modelop), op(x.op) {}

  // transop = dwop, dlop
  // modelop = injop, iwop
  // op = transop, modelop
  
  void IWaveFWIOp::apply(RVL::Vector<float> const & x,
			 RVL::Vector<float> & y) const {
    try {
      y.zero();

#ifndef VERBOSE
      
      RVL::Operator<float>::export_apply(*op,x,y);
      
#else 
      cerr<<"IWaveFWIOp::apply\n";

      RVL::Vector<float> z1(dwop->getRange());
      RVL::Vector<float> z2(dlop->getRange());
      RVL::Vector<float> z3(injop->getRange());

      RVL::Operator<float>::export_apply(*dwop,x,z1);
      RVL::Operator<float>::export_apply(*dlop,z1,z2);
      RVL::Operator<float>::export_apply(*injop,z2,z3);
      RVL::Operator<float>::export_apply(iwop,z3,y);
      {
#ifdef IWAVE_USE_MPI
	{
	  RVL::RVLMin<float> svmin;
	  RVL::RVLMax<float> svmax;
	  RVL::MPISerialFunctionObjectRedn<float,float> vmin(svmin);
	  RVL::MPISerialFunctionObjectRedn<float,float> vmax(svmax);
	  x.eval(vmin);
	  x.eval(vmax);
	  cerr<<"*** input min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
	}
	{
	  RVL::RVLMin<float> svmin;
	  RVL::RVLMax<float> svmax;
	  RVL::MPISerialFunctionObjectRedn<float,float> vmin(svmin);
	  RVL::MPISerialFunctionObjectRedn<float,float> vmax(svmax);
	  z1.eval(vmin);
	  z1.eval(vmax);
	  cerr<<"*** z1 min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
	}
	{
	  RVL::RVLMin<float> svmin;
	  RVL::RVLMax<float> svmax;
	  RVL::MPISerialFunctionObjectRedn<float,float> vmin(svmin);
	  RVL::MPISerialFunctionObjectRedn<float,float> vmax(svmax);
	  z2.eval(vmin);
	  z2.eval(vmax);
	  cerr<<"*** z2 min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
	}
	{
	  RVL::RVLMin<float> svmin;
	  RVL::RVLMax<float> svmax;
	  RVL::MPISerialFunctionObjectRedn<float,float> vmin(svmin);
	  RVL::MPISerialFunctionObjectRedn<float,float> vmax(svmax);	  
	  z3.eval(vmin);
	  z3.eval(vmax);
	  cerr<<"*** z3 min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";		}
	{
	  RVL::RVLMin<float> svmin;
	  RVL::RVLMax<float> svmax;
	  RVL::MPISerialFunctionObjectRedn<float,float> vmin(svmin);
	  RVL::MPISerialFunctionObjectRedn<float,float> vmax(svmax);
	  y.eval(vmin);
	  y.eval(vmax);
	  cerr<<"*** output min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
	}
#else
	{
	  RVL::RVLMin<float> vmin;
	  RVL::RVLMax<float> vmax;
	  x.eval(vmin);
	  x.eval(vmax);
	  cerr<<"*** input min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
	}
	{
	  RVL::RVLMin<float> vmin;
	  RVL::RVLMax<float> vmax;
	  z1.eval(vmin);
	  z1.eval(vmax);
	  cerr<<"*** z1 min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
	}
	{
	  RVL::RVLMin<float> vmin;
	  RVL::RVLMax<float> vmax;
	  z2.eval(vmin);
	  z2.eval(vmax);
	  cerr<<"*** z2 min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
	}
	RVL::Components<float> cz3(z3);
	for (int i=0;i<cz3.getSize(); i++) {
	  RVL::RVLMin<float> vmin;
	  RVL::RVLMax<float> vmax;
	  cz3[i].eval(vmin);
	  cz3[i].eval(vmax);
	  cerr<<"*** z3 comp "<<i<<" min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
	}	
	{
	  RVL::RVLMin<float> vmin;
	  RVL::RVLMax<float> vmax;
	  y.eval(vmin);
	  y.eval(vmax);
	  cerr<<"*** output min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";	}
#endif
      }
#endif
    }
    catch (RVL::RVLException e) {
      e<<"\ncalled from IWaveFWIOp::apply\n";
      throw e;
    }
  }
    
  void IWaveFWIOp::applyDeriv(RVL::Vector<float> const & x,
			      RVL::Vector<float> const & dx,
			      RVL::Vector<float> & dy) const {
    try {
      RVL::Operator<float>::export_applyDeriv(*op,x,dx,dy);
    }
    catch (RVL::RVLException e) {
      e<<"\ncalled from IWaveFWIOp::applyDeriv\n";
      throw e;
    }
  }
    
  void IWaveFWIOp::applyAdjDeriv(RVL::Vector<float> const & x,
				 RVL::Vector<float> const & dy,
				 RVL::Vector<float> & dx) const {
    try {
      RVL::Operator<float>::export_applyAdjDeriv(*op,x,dy,dx);
    }
    catch (RVL::RVLException e) {
      e<<"\ncalled from IWaveFWIOp::applyAdjDeriv\n";
      throw e;
    }
  }
    
  void IWaveFWIOp::applyDeriv2(const RVL::Vector<float> & x,
			       const RVL::Vector<float> & dx1,
			       const RVL::Vector<float> & dx2,
			       RVL::Vector<float> & dy) const {
    try {
      RVL::Operator<float>::export_applyDeriv2(*op,x,dx1,dx2,dy);
    }
    catch (RVL::RVLException e) {
      e<<"\ncalled from IWaveFWIOp::applyDeriv2\n";
      throw e;
    }
  }
    
  void IWaveFWIOp::applyAdjDeriv2(const RVL::Vector<float> & x,
				  const RVL::Vector<float> & dy,
				  const RVL::Vector<float> & dx2,
				  RVL::Vector<float> & dx1) const {
    try {
      RVL::Operator<float>::export_applyAdjDeriv2(*op,x,dy,dx2,dx1);
    }
    catch (RVL::RVLException e) {
      e<<"\ncalled from IWaveFWIOp::applyAdjDeriv2\n";
      throw e;
    }
  }

  ostream & IWaveFWIOp::write(ostream & str) const {
    str<<"IWaveFWIOp::write\n";
    return str;
  }

  void IWaveOpt(PARARRAY pars,
		RVL::Operator<float> const & op,
		std::shared_ptr<RVL::LinearOp<float> const> prec,
		RVL::Vector<float> const & d,
		RVL::Vector<float> & m) {
    try {
      /*
      Components<float> cm(m);

      for (int i=0;i<cm.getSize();i++) {
#ifdef IWAVE_USE_MPI
	RVL::RVLMin<float> svmin;
	RVL::RVLMax<float> svmax;
	RVL::MPISerialFunctionObjectRedn<float,float> vmin(svmin);
	RVL::MPISerialFunctionObjectRedn<float,float> vmax(svmax);
	cm[i].eval(vmin);
	cm[i].eval(vmax);
	cerr<<"StdOLS bfore \n";
	cerr<<"  background model comp "<<i<<" min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
	int nnn;
	cin>>nnn;
#else
	RVL::RVLMin<float> vmin;
	RVL::RVLMax<float> vmax;
	cm[i].eval(vmin);
	cm[i].eval(vmax);
	cerr<<"  background model comp "<<i<<" min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
#endif
      }
      */
      /* output stream */
      std::stringstream outstr;

      std::ostream * optr = NULL;
      std::string outfile = RVL::valparse<std::string>(pars,"outfile","");
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
      if (!optr) {
	RVL::RVLException e;
	e<<"Error: IWaveOpt\n";
	e<<"  failed to initialize output pointer\n";
	throw e;
      }

      /* switch on method */
      RVLAlg::Algorithm * alg = NULL;
      std::string optmethod = RVL::valparse<std::string>(pars,"OptMethod","lbfgs");
#ifdef VERBOSE
      cerr<<"IWaveOpt: optmethod="<<optmethod<<"\n";
#endif
      if (optmethod == "lbfgs") {
	RVL::StdLeastSquaresFcnlGN<float> f(op,d);
	alg = new RVLUmin::LBFGSBT<float>
	  (f, m,
	   RVL::valparse<float>(pars,"InvHessianScale",1.0f),
	   RVL::valparse<int>(pars,"MaxInvHessianUpdates",5),
	   RVL::valparse<int>(pars,"MaxSubSteps",10),
	   RVL::valparse<bool>(pars,"VerboseDisplay",true), 
	   RVL::valparse<float>(pars,"InitStepBound",1.0f), 
	   RVL::valparse<float>(pars,"MinDecrease",0.1f),
	   RVL::valparse<float>(pars,"GoodDecrease",0.9f), 
	   RVL::valparse<float>(pars,"StepDecrFactor",0.5f), 
	   RVL::valparse<float>(pars,"StepIncrFactor",1.8f),
	   RVL::valparse<float>(pars,"MaxFracDistToBdry",1.0), 
	   RVL::valparse<float>(pars,"MinStepTol",1.e-06),     
	   RVL::valparse<int>(pars,"MaxSteps",10), 
	   RVL::valparse<float>(pars,"AbsGradThresh",0.0), 
	   RVL::valparse<float>(pars,"RelGradThresh",1.e-2), 
	   *optr);
      }
      else if (optmethod == "trcg") {
	RVL::ResidualOperator<float> rop(op,d);      
	RVLUmin::TRGNAlg<float, RVLUmin::CGNEPolicy<float> > * tralg = new RVLUmin::TRGNAlg<float, RVLUmin::CGNEPolicy<float> >
	  (rop, prec, m,
	   RVL::valparse<int>(pars,"MaxSteps",10),             // _maxcount,
	   RVL::valparse<float>(pars,"ResidualTol",0.0f),       // _jtol,
	   RVL::valparse<float>(pars,"AbsGradThresh",0.0f),    // _agtol,
	   RVL::valparse<float>(pars,"RelGradThresh",1.0e-2),  // _rgtol,
	   RVL::valparse<float>(pars,"MinDecrease",0.1f),      // _eta1
	   RVL::valparse<float>(pars,"GoodDecrease",0.9f),     // _eta2
	   RVL::valparse<float>(pars,"StepDecrFactor",0.5f),   // _gamma1
	   RVL::valparse<float>(pars,"StepIncrFactor",1.8f),   // _gamma2
	   RVL::valparse<float>(pars,"MinStepTol",1.e-06),     // min step as frac
	   *optr);
	
	// assign CG params
	tralg->assign
	  (RVL::valparse<float>(pars,"CGNE_ResTol",0.0f),      // rtol
	   RVL::valparse<float>(pars,"CGNE_GradTol",0.001f),   // nrtol,
	   RVL::valparse<float>(pars,"InitStepBound",1.0f),    // Delta
	   RVL::valparse<int>(pars,"MaxSubSteps",10),          // maxcount
	   RVL::valparse<bool>(pars,"VerboseDisplay",true));   // verbose
	alg=tralg;
      }
      else {
	RVL::RVLException e;
	e<<"ERROR: IWaveOpt\n";
	e<<"  no optimization algorithm selected, or algoritm specified is not\n";
	e<<"  supported. Currently supported algorithms:\n";
	e<<"  OptMethod = lbfgs (Limited Memory Broyden-Fletcher-Goldfarb-Shanno)\n";
	e<<"  OptMethod = trcg (Trust Region Gauss-Newton with CG linear solver\n)";
	throw e;
      }

#ifdef VERBOSE
      cerr<<"IWaveOpt -> run\n";
#endif
      alg->run();
#ifdef VERBOSE
      cerr<<"IWaveOpt -> report\n";
#endif
      delete alg;
      /*
      for (int i=0; i<cm.getSize();i++) {
#ifdef IWAVE_USE_MPI
	RVL::RVLMin<float> svmin;
	RVL::RVLMax<float> svmax;
	RVL::MPISerialFunctionObjectRedn<float,float> vmin(svmin);
	RVL::MPISerialFunctionObjectRedn<float,float> vmax(svmax);
	cm[i].eval(vmin);
	cm[i].eval(vmax);
	cerr<<"StdOLS after\n";
	cerr<<"  background model comp "<<i<<" min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
	int nnn;
	cin>>nnn;
#else
	RVL::RVLMin<float> vmin;
	RVL::RVLMax<float> vmax;
	cm[i].eval(vmin);
	cm[i].eval(vmax);
	cerr<<"  background model comp "<<i<<" min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
#endif	
      }
      */
    
      if ((outfile.size()>0) && (retrieveRank()==0) && (optr)) {
	std::ofstream * ofptr = NULL;
	ofptr = dynamic_cast<std::ofstream *>(optr);
	if (ofptr) {
	  ofptr->flush();
	  ofptr->close();
	}
	delete optr;
      }

      return;

    }
    catch (RVL::RVLException & e) {
      e<<"\ncalled from IWaveOpt\n";
      throw e;

    }
  }

  void IWaveResMod(PARARRAY pars,
		   IWaveFWIOp const & op,
		   RVL::Vector<float> const d,
		   RVL::Vector<float> const r) {

    try {

      // create model vector, optionally assign components to
      // archival files
      RVL::Vector<float> m(op.getTransformOp().getRange());
      RVL::Components<float> cm(m);
      std::vector<std::string> mkeys = op.getDomainKeys();

      size_t flag = 0;
      for (int i=0; i<mkeys.size(); i++) {
	// chop off "_est"
	std::string mdl = mkeys[i].substr(0,mkeys[i].size()-4);
	// tack on "_inv", check to see if file is specified
	std::string ename = RVL::valparse<std::string>(pars,mdl+"_est","");
	std::string fname = RVL::valparse<std::string>(pars,mdl+"_inv","");
	//	cerr<<"IWaveResMod:\n";
	//	cerr<<"  model index="<<i<<" key="<<mkeys[i]<<" root="<<mdl<<" inv="<<fname<<"\n";
	// if so, connect to vector
	if ((ename.size() > 0) && (fname.size() > 0)) {
	  RVL::AssignFilename fn(fname);
	  cm[i].eval(fn);
	  flag += fname.size();
	  /*
#ifdef IWAVE_USE_MPI
	  RVL::RVLMin<float> svmin;
	  RVL::RVLMax<float> svmax;
	  RVL::MPISerialFunctionObjectRedn<float,float> vmin(svmin);
	  RVL::MPISerialFunctionObjectRedn<float,float> vmax(svmax);
	  cm[i].eval(vmin);
	  cm[i].eval(vmax);
	  cerr<<"  background model comp "<<i<<" min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
	  int nnn;
	  cin>>nnn;
#else
	  RVL::RVLMin<float> vmin;
	  RVL::RVLMax<float> vmax;
	  cm[i].eval(vmin);
	  cm[i].eval(vmax);
	  cerr<<"  background model comp "<<i<<" min="<<vmin.getValue()<<" max="<<vmax.getValue()<<"\n";
#endif	
	  */
	}
      }
      // don't bother if the user hasn't asked for any model files
      if (flag > 0) {
	RVL::OperatorEvaluation<float> transeval(op.getTransformOp(),r);
	m.copy(transeval.getValue());
      }

      // create residual vector, optionally assign components to
      // archival files
      RVL::Vector<float> res(op.getModelingOp().getRange());
      RVL::Components<float> cres(res);
      std::vector<std::string> ekeys = op.getRangeKeys();

      size_t elag = 0;
      for (int i=0; i<ekeys.size(); i++) {
	std::string fname = RVL::valparse<std::string>(pars,ekeys[i]+"_res","");
	if (fname.size() > 0) {
	  RVL::AssignFilename fn(fname);
	  cres[i].eval(fn);
	  elag += fname.size();
	}
      }
      // don't bother if the user hasn't asked for any residual files
      if (elag > 0) {
	if (flag > 0) {
	  RVL::OperatorEvaluation<float> modeleval(op.getModelingOp(),m);
	  res.copy(modeleval.getValue());
	}
	else {
	  RVL::OperatorEvaluation<float> modeleval(op,r);
	  res.copy(modeleval.getValue());
	}
      }
      res.linComb(-1.0f,d);
    }
    catch (RVL::RVLException & e) {
      e<<"\ncalled from IWaveResMod\n";
      throw e;
    }
  }

  int IWaveOpApply(int argc, char ** argv,
		   void (*invfcn)(PARARRAY,
				  RVL::Builder<RVL::LinearOp<float>, GridOpData> &,
				  FILE *),
		   RVL::Builder<RVL::LinearOp<float>, GridOpData> & b) {

    try {

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

	invfcn(*pars,b,stream);
	
#ifdef IWAVE_USE_MPI
	/* end nontriv comm branch */
      }
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      ps_delete(&pars);
      fclose(stream);
      
      return 0;
    }
    catch (RVL::RVLException & e) {
      e<<"\ncalled from IWaveOpApply\n";
      throw e;
    }
  }

  void StraightOLS(PARARRAY pars,
		   RVL::Builder<RVL::LinearOp<float>, GridOpData> & b,
		   FILE * stream) {
    try {
#ifdef VERBOSE
      cerr<<"StraightOLS -> IWaveFWIOp constructor\n";
#endif
      TSOpt::IWaveFWIOp op(pars, stream);
#ifdef VERBOSE
      cerr<<"StraightOLS -> preconditioner factory\n";
#endif
      std::shared_ptr<GridOpData> h;
      h = make_shared<GridOpData>(&pars,op.getDomain());
      shared_ptr<LinearOp<float> const> prec = b.build(h);
#ifdef VERBOSE
      cerr<<"StraightOLS -> domain vector\n";
#endif
      RVL::Vector<float> r(op.getDomain());
#ifdef VERBOSE
      cerr<<"StraightOLS -> range vector\n";
#endif
      RVL::Vector<ireal> d(op.getRange());
#ifdef VERBOSE
      for (int i=0;i<op.getDomainKeys().size(); i++)
	cerr<<"StraightOLS: domain key "<<i<<" = "<<op.getDomainKeys()[i]<<"\n";
      cerr<<"StraightOLS->load domain vector\n";
#endif
      TSOpt::IWaveLoad(pars,r,op.getDomainKeys());
#ifdef VERBOSE
      for (int i=0;i<op.getRangeKeys().size(); i++)
	cerr<<"StraightOLS: range key "<<i<<" = "<<op.getRangeKeys()[i]<<"\n";
      cerr<<"StraightOLS->load range vector\n";
#endif      
      TSOpt::IWaveLoad(pars,d,op.getRangeKeys());
#ifdef VERBOSE
      cerr<<"StraightOLS: optimize\n";
#endif
      TSOpt::IWaveOpt(pars, op, prec, d, r);
#ifdef VERBOSE
      cerr<<"StraightOLS: cleanup\n";
#endif      
      TSOpt::IWaveResMod(pars,op,d,r);
#ifdef VERBOSE
      cerr<<"StraightOLS exit\n";
#endif
    }    
    catch (RVL::RVLException & e) {
      e<<"\ncalled from TSOpt::StraightOLS\n";
      throw e;
    }
  }
}
