#include "iwop.hh"
#include "functions.hh"
#include "linop_apps.hh"
#include <par.h>
#include "segyops.hh"
#include "gridops.hh"
#include "adjtest.hh"
#include "asg_defn.hh"
#include "sse_selfdoc.h"
#include "parser.h"
#include "blockop.hh"
#include "cgnealg.hh"
#include "alphaupdate.hh"
//#include <omp.h>

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
      /*
	for (int j=0; j< iwop.getIWaveRange().getSize(); j++) {
	std::string fn =
	RVL::valparse<std::string>(*pars,
	iwop.getIWaveRange().getKeys()[j]);
	RVL::AssignFilename af(fn);
	dm[j].eval(af);
	}
      */
      // mute
      TSOpt::SEGYLinMute mute(RVL::valparse<float>(*pars,"mute_slope",0.0f),
			      RVL::valparse<float>(*pars,"mute_zotime",0.0f),
			      RVL::valparse<float>(*pars,"mute_width",0.0f),
			      RVL::valparse<int>(*pars,"mute_mode",0));
      
      RVL::LinearOpFO<float> muteop(iwop.getIWaveRange(),iwop.getIWaveRange(),mute,mute);

      // taper
      TSOpt::SEGYTaper fo(RVL::valparse<std::string>(*pars,"taperpars",""));
      RVL::LinearOpFO<float> taperop(iwop.getRange(),
				     iwop.getRange(),
				     fo,fo);
    
      //      RVL::LinearRestrictOp<float> rlop(iwop,cm[0]);
      //      RVL::CompLinearOp<float> mrlop(rlop,muteop);
      //      RVL::CompLinearOp<float> tmrlop(mrlop,taperop);

      // alternative construction: first compose with succesion of
      // linear ops on output side, then restrict to LOP

      RVL::LOVCompLinOp<float> mrop(iwop,muteop);
      RVL::LOVCompLinOp<float> tmrop(mrop,taperop);
      RVL::LinearRestrictOp<float> tmrlop(tmrop,cm[0]);      

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
      std::ostream & res = *optr;
    
      // CG parameters
      float rtol=RVL::valparse<float>(*pars,"ResidualTol",
				      100.0*numeric_limits<float>::epsilon());
      float nrtol=RVL::valparse<float>(*pars,"GradientTol",
				       100.0*numeric_limits<float>::epsilon());
      int maxcount=RVL::valparse<int>(*pars,"MaxIter",10);
      float maxstep=RVL::valparse<float>(*pars,"MaxStep",
					 numeric_limits<float>::max());



      // discrepancy principle params
      float alphatol=RVL::valparse<float>(*pars,"AlphaUpdateTol",0.01);
      int maxalpha=RVL::valparse<int>(*pars,"MaxAlphaUpdate",1);
      float ubnd = RVL::valparse<float>(*pars,"RelResUpperBnd",0.12);
      float lbnd = RVL::valparse<float>(*pars,"RelResLowerBnd",0.08);
      float decr = RVL::valparse<float>(*pars,"AlphaDecrFac",0.5);

      // optional outputs
      std::string dataest = valparse<std::string>(*pars,"est_data_p","");
      std::string datares = valparse<std::string>(*pars,"res_data_p","");
      std::string dataann = valparse<std::string>(*pars,"res_anni_p","");

      // load data
      RVL::StdProductSpace<float> tsp(tmrlop.getRange(),tmrlop.getDomain());
      RVL::Vector<float> dd(tsp);
      RVL::Components<float> cdd(dd);
      RVL::AssignFilename ddfn(RVL::valparse<std::string>(*pars,"data_p"));
      cdd[0].eval(ddfn);
      cdd[1].zero();

      // remove part of data predicted from initial input model
      RVL::Vector<float> ldd(tsp);
      ldd.copy(dd);
      RVL::Vector<float> pd(tmrlop.getRange());
      tmrlop.applyOp(cm[1],pd);
      Components<float> cldd(ldd);
      cldd[0].linComb(-1.0,pd);
      
      // vector workspace, may be archived
      RVL::StdProductSpace<float> rsp(tmrlop.getRange(),tmrlop.getDomain());
      RVL::Vector<float> est(rsp);
      RVL::Components<float> cest(est);
      RVL::Vector<float> dres(tmrlop.getRange());

      // archive if filenames provided
      if (dataest.size()>0) {
	AssignFilename estfn(dataest);
	cest[0].eval(estfn);
      }
      if (dataann.size()>0) {
	AssignFilename annfn(dataann);
	cest[1].eval(annfn);
      }
      if (datares.size()>0) {
	RVL::AssignFilename resfn(datares);
	dres.eval(resfn);
      }
            
      // initialize alpha
      float alpha =0.0f;
      std::string iwnm = RVL::valparse<std::string>(*pars, "weightin","");
      if (iwnm.size() > 0) {
	ifstream iwstr(iwnm);
	iwstr >> alpha;
	iwstr.close();
      }
      else {
	alpha = RVL::valparse<float>(*pars, "weight", 0.0);
      }
      // if set, precondition by (1 + alpha^2 offset^2)^{-1}
      int precond=RVL::valparse<int>(*pars,"precond",0);

      // return of resid, normal resid - must persist
      float rnorm;
      float nrnorm;

      // workspace for CG - initially 0 in each call, on
      // success used to update cm[1]
      RVL::Vector<float> x(tmrlop.getDomain());

      for (int ialpha=0; ialpha<maxalpha; ialpha++) {

	// trace scaling by weighted offset
	float alphasq=alpha*alpha;
	TSOpt::TraceScaleFO ts(alpha,"offset");
	RVL::LinearOpFO<float> tsop(tmrlop.getDomain(),tmrlop.getDomain(),ts,ts);
	// block operator (rlop, tsop)^T
	RVL::TensorLinearOp<float> top(tmrlop,tsop);

	res<<endl<<"*******************************************************"<<endl;
	res<<"* Acoustic Variable Density Source Inversion via";
	res<<"* Conjugate Gradient Algorithm for Normal Eqns"<<endl;
	res<<"* max iterations       = "<<maxcount<<endl;
	res<<"* residual tolerance   = "<<rtol<<endl;
	res<<"* normal res tolerance = "<<nrtol<<endl;
	res<<"* trust radius         = "<<maxstep<<endl;
	res<<"*******************************************************"<<endl;

	if (precond) {
	  float shift=1.0f;
	  float pwr=-1.0f;
	  TSOpt::SSEScaleFO pfo(shift,alpha,pwr);
	  RVL::LinearOpFO<float> prec(top.getDomain(),top.getDomain(),pfo,pfo);
	  RVLUmin::CGNEAlg<float> alg(x,top,prec,ldd,
	  			      rnorm, nrnorm, rtol, nrtol,
				      maxcount, maxstep, res);
	  alg.run();
	}
	else {
	  RVLUmin::CGNEAlg<float> alg(x,top,ldd,
	  			      rnorm, nrnorm, rtol, nrtol,
				      maxcount, maxstep, res);
	  alg.run();
	}


	// post-optimization computations
	// compute estimated data
	top.applyOp(x,est);
	// compute residual - S[c]f-d rather than other way round
	dres.copy(cest[0]);
	dres.linComb(-1.0f,cldd[0]);
	// compute penalty output (unweighted)
	TSOpt::TraceScaleFO ts1(1.0f,"offset");
	RVL::LinearOpFO<float> tsop1(tmrlop.getDomain(),
				     tmrlop.getDomain(),ts1,ts1);
	tsop1.applyOp(x,cest[1]);
    
	if (maxalpha > 1) {
	  float alphasqnext =
	    RVLUmin::alphaupdate(cldd[0],dres,cest[1],
				 alphasq,ubnd*ubnd,lbnd*lbnd,decr*decr);
	  float alphanext = sqrt(alphasqnext);
	  res<<"current alpha              = "<<alpha<<"\n";
	  float dn = cldd[0].norm();
	  float en = dres.norm();
	  res<<"data norm                  = "<<dn<<"\n";
	  res<<"res norm                   = "<<en<<"\n";
	  res<<"upper bnd res              = "<<dn*ubnd<<"\n";
	  res<<"lower bnd res              = "<<dn*lbnd<<"\n";
	  res<<"pen norm                   = "<<cest[1].norm()<<"\n";
	  res<<"updated alpha              = "<<alphanext<<"\n";

	  cerr<<"fabs(alpha-alphanext)="<<fabs(alpha-alphanext)<<" alphatol*alphanext="<<alphatol*alphanext<<endl;
	  // check whether penalty is updated
	  if ((alphanext < std::numeric_limits<float>::epsilon()*alpha) ||
	      (fabs(alpha-alphanext) < alphatol*alphanext)) {
	    break;
	  }
	  else {
	    alpha=alphanext;
	  }
      
	}

      }

      // update cm[1] with result of cg loop
      cm[1].linComb(1.0f,x);
      
      // record weight in file to transfer to Python
      std::string ownm = RVL::valparse<std::string>(*pars,"weightout","");
      if (ownm.size() > 0) {
	ofstream wstr(ownm);
	wstr<<alpha<<"\n";
	wstr.flush();
	wstr.close();
      }
      
      // apply Born adjoint to data residual to compute VPM gradient
      // switch on if any gradient filenames are provided - name of
      // nonlin param + "_grad"
      std::vector<std::string> gradlist;
      size_t tot=0;
      for (int j=0; j< iwop.getNonLinDomain().getSize(); j++) {
	std::string fn =
	  RVL::valparse<std::string>(*pars,
				     iwop.getNonLinDomain().getKeys()[j]
				     + "_grad","");
	cerr<<iwop.getNonLinDomain().getKeys()[j]<<"_grad = " <<fn<<endl;
	gradlist.push_back(fn);
	tot += fn.size();
      }
      
      // switch on gradient computation if any gradient files indicated
      RVL::Vector<float> pm(iwop.getNonLinDomain());
      RVL::Components<float> cpm(pm);
      RVL::Vector<float> gpm(iwop.getNonLinDomain());
      RVL::Components<float> cgpm(gpm);
      
      if (tot) {
	for (int j=0; j< iwop.getNonLinDomain().getSize(); j++) {	  
	  if (gradlist[j].size()>0) {
	    RVL::AssignFilename af(gradlist[j]);
	    cpm[j].eval(af);
	    res<<" *** computing L2 gradient component "<<gradlist[j]<<"\n";
	  }
	}	  
	// restrict to first component, evaluate
	// note that dres is full residual, even though
	// it was computed with reference data removed
	pm.zero();
	RVL::RestrictOp<float> rtmrop(tmrop,m,0);
	RVL::OperatorEvaluation<float> rtmropeval(rtmrop,cm[0]);
	rtmropeval.getDeriv().applyAdjOp(dres,pm);
	  
	// pm is now L2 gradient. if precond is provided, apply -
	// use in FO form to avoid having to create full op interface
	IPNT rad;
	for (int i=1; i < RARR_MAX_NDIM+1; i++) {
	  std::stringstream foo;
	  foo<<i;
	  std::string barf = "rect";
	  rad[i-1] = RVL::valparse<int>(*pars,barf+foo.str(),0);
	  if (rad[i-1]<0) {
	    RVL::RVLException e;
	    e<<"Error: asgiva application of preconditioner\n";
	    e<<"  one of the axis box lengths rect["<<i<<"] < 0\n";
	    throw e;
	  }
	}
	  
	int rep = RVL::valparse<int>(*pars,"repeat",1);
	if (rep<1) {
	  RVL::RVLException e;
	  e<<"Error: asgiva application of preconditioner\n";
	  e<<"  repetitions of moving average given as < 1\n";
	  throw e;
	}
	  
	// look for filenames for precond grad components, count
	std::vector<std::string> pgradlist;
	size_t ptot=0;
	for (int j=0; j< iwop.getNonLinDomain().getSize(); j++) {
	  std::string fn =
	    RVL::valparse<std::string>(*pars,
				       iwop.getNonLinDomain().getKeys()[j]
				       + "_pgrad","");
	  pgradlist.push_back(fn);
	  ptot += fn.size();
	}	  
	if (ptot) {
	  TSOpt::GridMAFO f(rep,rad);
	  for (int j=0; j< iwop.getNonLinDomain().getSize(); j++) {	  
	    if (pgradlist[j].size()>0) {
	      RVL::AssignFilename af(pgradlist[j]);
	      cgpm[j].eval(af);
	      res<<" *** computing L2 gradient component "<<pgradlist[j]<<"\n";		
	    }
	  }
	    
	  gpm.eval(f,pm);
	}	  	    
      }

      // update bulkmod - input (a) first step as fraction of distance
      // to boundary of feasible set, (b) max permitted fraction of
      // distance to boundary, (c) Goldstein-Armijo params, (d)
      // backtrack fraction, (e) max bactracks/ internal doublings
      // these names are from standard Umin interface
      /*
  InitStepBound = 1.0f         max length of first step
    MinDecrease = 0.1f         minimum fractional decrease (1st G-A param)
   GoodDecrease = 0.9f         decrease threshhold for step increase
                                   (2nd G-A param)
 StepDecrFactor = 0.5f         decrease step or trust radius by this factor
                                   if fractional decrease less than MinDecrease
 StepIncrFactor = 1.8f         increase step or trust radius by this factor
                                   if fractional decrease more than GoodDecrease
 MaxFracDistToBdry 
                = 1.0          do not allow update to exceed this fractional
                                   step to boundary of feasible set (limits step
                                   or trust radius)
      */

      std::string bulkupd =
	RVL::valparse<std::string>(*pars,"bulkmod_upd","");
      if (bulkupd.size()>0) {
	
	int maxback = RVL::valparse<int>(*pars,"MaxSubSteps",10);
	float firststep = RVL::valparse<float>(*pars,"InitStepBound",0.1);
	//x	float maxfrac = RVL::valparse<float>(*pars,"MaxFracDistToBdry",0.5);
	float mindecr = RVL::valparse<float>(*pars,"MinDecrease",0.1);
	float maxdecr = RVL::valparse<float>(*pars,"GoodDecrease",0.7);
	float decrfac = RVL::valparse<float>(*pars,"StepDecrFactor",0.5);
	float incrfac = RVL::valparse<float>(*pars,"StepIncrFactor",1.8);

	// compute distance to boundary
	// simple since bounds are const and density not changed so
	// bulkmod bound is equivalent to velocity bound
	float cmin = RVL::valparse<float>(*pars,"cmin");
	float cmax = RVL::valparse<float>(*pars,"cmax");      
	float dmin = RVL::valparse<float>(*pars,"dmin");
	float dmax = RVL::valparse<float>(*pars,"dmax");

	float kmin = dmin*cmin*cmin;
	float kmax = dmax*cmax*cmax;      

	// distance of current iterate to boundary
#ifdef IWAVE_USE_MPI
	RVL::RVLMin<float> svmin;
	RVL::RVLMax<float> svmax;
	RVL::MPISerialFunctionObjectRedn<float,float> vmin(svmin);
	RVL::MPISerialFunctionObjectRedn<float,float> vmax(svmax);
#else
	RVL::RVLMin<float> vmin;
	RVL::RVLMax<float> vmax;
#endif

#ifdef IWAVE_USE_MPI
	RVL::RVLMin<float> dsvmin;
	RVL::RVLMax<float> dsvmax;
	RVL::MPISerialFunctionObjectRedn<float,float> dvmin(dsvmin);
	RVL::MPISerialFunctionObjectRedn<float,float> dvmax(dsvmax);
#else
	RVL::RVLMin<float> dvmin;
	RVL::RVLMax<float> dvmax;
#endif     
	// this is specific - 0th comp is bulk mod
	cm0[0].eval(vmin);
	cm0[0].eval(vmax);
	cgpm[0].eval(dvmin);
	cgpm[0].eval(dvmax);

	// scale gradient to have unit L^infty length, to avoid other issues.
	RVL::Vector<float> sgrad(iwop.getNonLinDomain()[0]);
	float grecip;
	if (RVL::ProtectedDivision<float>(1.0f,
					  iwave_max(fabs(dvmin.getValue()),
						    fabs(dvmax.getValue())),
					  grecip)) {
	  RVL::RVLException e;
	  e<<"Error: asgiva\n";
	  e<<"  L^infty norm of gradient too small\n";
	  throw e;
	}
	sgrad.scale(grecip,cgpm[0]);
      
	float dist = iwave_min(vmin.getValue()-kmin, kmax-vmax.getValue());
	cerr<<"distance to bdry = "<<dist<<endl;

	// test vector
	RVL::Vector<float> x(iwop.getDomain());
	RVL::Components<float> cx(x);
	RVL::Components<float> cx0(cx[0]);
	if (bulkupd.size()>0) {
	  RVL::AssignFilename af(bulkupd);
	  cx0[0].eval(af);
	}

	int upd=1;
	int stepback=0;
	float step = -dist*firststep;

	// compute rate of decrease
	float rate = -cpm[0].inner(sgrad);
      
	while (upd && (stepback < maxback)) { 
      
	  // copy base for this step
	  x.copy(m);
	  // update bulk mod component
	  cx0[0].linComb(step,sgrad);
	
	  // create, run CG
	  // first reserve last values of residual, normal residual
	  float rnorm0=rnorm;
	
	  // restrict op, block op
	  RVL::LinearRestrictOp<float> xmrlop(tmrop,cx[0]);
	  TSOpt::TraceScaleFO ts(alpha,"offset");
	  RVL::LinearOpFO<float> tsop(xmrlop.getDomain(),xmrlop.getDomain(),ts,ts);
	  RVL::TensorLinearOp<float> xop(xmrlop,tsop);      
	  RVLUmin::CGNEAlg<float> alg(cx[1],xop,dd,
				      rnorm, nrnorm, rtol, nrtol,
				      maxcount, maxstep, res);
	  alg.run();

	  res<<"prior objective value      = "<<rnorm0<<"\n";
	  res<<"new objective value        = "<<rnorm<<"\n";
	  res<<"step                       = "<<step<<"\n";
	  res<<"rate                       = "<<rate<<"\n";
	  // check GA conditions
	  float predred = -step*rate;
	  float actred = 0.5*(rnorm-rnorm0);
	  res<<"predicted decrease         = "<<predred<<"\n";
	  res<<"actual decrease            = "<<actred<<"\n";
	  if (actred>mindecr*predred) {
	    step *= decrfac; upd=1; stepback++;
	    res<<"insufficient decrease - backtrack\n";
	    res<<"step                       = "<<step<<"\n";
	  }
	  else if ((actred<maxdecr*predred) &&
		   (stepback==0)) {
	    step *= incrfac; upd=1; stepback++;
	    res<<"large decrease - internal doubling\n";
	    res<<"step                       = "<<step<<"\n";	  
	  }
	  else {
	    res<<"accepted step\n";
	    upd=0;
	  }
	}
      }
      
      // finish report
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
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
    ps_delete(&pars);
    iwave_fdestroy();
#ifdef IWAVE_USE_MPI    
    MPI_Finalize();
#endif
    exit(0);
  }
  catch (RVL::RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
}
