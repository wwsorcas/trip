#include "asg_defn.hh"
#include "born.hh"
#include <par.h>
#include "appinv_selfdoc.h"
#include "appinv_utils.hh"
#include "adjtest.hh"
#include "cgnealg.hh"
#include "linop_apps.hh"
#include "blockop.hh"

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
      int freechoice = RVL::valparse<int>(*pars,"freesurface");
      string prefix;
      if (freechoice) prefix="free";
      else prefix="absb";

// precond flag: if set, implies use of appinv
      int precond = RVL::valparse<int>(*pars,"precond",0);
      if (precond) {
	ps_slint(*pars,"appinv",1);
	ps_slint(*pars,"ft",1);
	ps_slint(*pars,"ab",1);
      }
      else {
	ps_slint(*pars,"appinv",0);
	ps_slint(*pars,"ft",0);
	ps_slint(*pars,"ab",0);
      }
      
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
	if (retrieveGlobalRank()==0) cerr<<"loaded "<<fn<<" on key "<<key<<" into BornIWOP ref component "<<j<<"\n";
#endif
      }

      for (int j=0; j<bkeys.size(); j++) {
	std::string key = prefix + bkeys[j]+"pert";
	std::string fn = RVL::valparse<std::string>(*pars, key);
	RVL::AssignFilename af(fn);
	cm1[j].eval(af);
#ifdef IWAVE_VERBOSE
	if (retrieveGlobalRank()==0) cerr<<"loaded "<<fn<<" on key "<<key<<" into BornIWOP pert component "<<j<<"\n";
#endif
      }
      
      for (int j=0; j< biwop.getRangeKeys().size(); j++) {
	std::string key = prefix + biwop.getRangeKeys()[j];
	std::string fn = RVL::valparse<std::string>(*pars, key);
	RVL::AssignFilename af(fn);
	cd[j].eval(af);
#ifdef IWAVE_VERBOSE	
	if (retrieveGlobalRank()==0) cerr<<"loaded "<<fn<<" on key "<<key<<" into BornIWOP range component "<<j<<"\n";
#endif
      }

#ifdef IWAVE_VERBOSE
      if (retrieveGlobalRank()==0) cerr<<"deriv = "<<deriv<<" adjoint = "<<adjoint<<"\n";
#endif

      RVL::LinearRestrictOp<float> lrop(biwop,cm[0]);

      if (RVL::valparse<int>(*pars,"adjtest",0)) {
	RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
	std::stringstream stuff;
	RVL::AdjointTest(lrop,rnd,stuff);
	if (retrieveGlobalRank()==0) {
	  fprintf(stream,stuff.str().c_str());
	}
      }

      // mute
      TSOpt::SEGYLinMute mute(RVL::valparse<float>(*pars,"mute_slope",0.0f),
			      RVL::valparse<float>(*pars,"mute_zotime",0.0f),
			      RVL::valparse<float>(*pars,"mute_width",0.0f));
      
      RVL::LinearOpFO<float> muteop(lrop.getRange(),lrop.getRange(),mute,mute);

      TSOpt::GridZTaperOp top(lrop.getDomain(),
			      RVL::valparse<float>(*pars,"top",0.0f),
			      RVL::valparse<float>(*pars,"bot",0.0f));

      RVL::CompLinearOp<float> cop;
      cop.setNext(top);
      cop.setNext(lrop);
      cop.setNext(muteop);

      TSOpt::GridRadialScaleOp hop(lrop.getDomain(),
				   RVL::valparse<float>(*pars,"alpha",0.0f));
      
      /* output stream */
      std::stringstream outstr;
      std::shared_ptr<std::ostream> optr;
      std::string outfile = RVL::valparse<std::string>(*pars,"outfile","");
      if ((retrieveGlobalRank()==0) && (outfile.size() > 0)) {
	optr = make_shared<std::ofstream>(outfile.c_str());
	(*optr)<<scientific;
	fprintf(stream,"output on file %s\n",outfile.c_str());
      }
      else {
	optr = make_shared<std::stringstream>();
      }

      // detect number of iterations desired
      int maxcount=RVL::valparse<int>(*pars,"MaxIter",0);

      // deriv, adjoint
      int deriv = RVL::valparse<int>(*pars,"deriv",0);
      if (deriv > 0) {
	RVL::RVLException e;
	e<<"Error: borninv\n";
	e<<"  deriv > 0 not implemented yet\n";
	throw e;
      }
      int adjoint = RVL::valparse<int>(*pars,"adjoint",1);

      if (precond) {
	if (deriv > 0) {
	  RVL::RVLException e;
	  e<<"Error: borninv, preconditioned mode (precond=1)\n";
	  e<<"  derivative or adjoint derivative (tomo op) incompatible with precond\n";
	  throw e;
	}
	if (adjoint == 0) {
	  RVL::RVLException e;
	  e<<"Error: borninv, preconditioned mode (precond=1)\n";
	  e<<"  only compatible with adjoint mode (adjoint=1)\n";
	  throw e;
	}
	// set up model weight op 
	float maxmaxspf = RVL::valparse<float>(*pars,"maxfreq")/
	  RVL::valparse<float>(*pars,"cmin");
	float minmaxspf = RVL::valparse<float>(*pars,"maxfreq")/
	  RVL::valparse<float>(*pars,"cmax");
#ifdef IWAVE_VERBOSE
	if (retrieveGlobalRank()==0) 
	  cerr<<"model space weight op constructor\n";
#endif
	RVL::Vector<float> & bulkvec = cm0[0];
	//	    RVL::Vector<float> buoyvec(biwop.getIWOP().getNonLinDomain()[1]);
	// a crude hack, ignores unit conflict
	RVL::Vector<float> buoyvec(bulkvec.getSpace());	    
	std::string fn = RVL::valparse<std::string>(*pars,prefix + "buoyancy");
	RVL::AssignFilename af(fn);
	buoyvec.eval(af);
	
	RVL::Vector<float> & extdvec = cm[1];
	
	TSOpt::ASGModelWeightOp wm(bulkvec,buoyvec,extdvec,
				   RVL::valparse<float>(*pars,"Zlocut", 0.0625*minmaxspf),
				   RVL::valparse<float>(*pars,"Zlopas", 0.125*minmaxspf),
				   RVL::valparse<float>(*pars,"Zhipas", 0.9375*maxmaxspf),
				   RVL::valparse<float>(*pars,"Zhicut", maxmaxspf),
				   RVL::valparse<int>(*pars,"symm",1));
		
	if (maxcount==0) {
	  float dnorm=d.norm();
	  if (retrieveGlobalRank()==0) cerr<<" appinv1: norm of input data = "<<dnorm<<endl;
	  RVL::Vector<float> tmp(cm[1].getSpace());
	  //	  cop.applyAdjOp(d,tmp);
	  RVL::Vector<float> dd(lrop.getRange());
	  muteop.applyAdjOp(d,dd);
	  float ddnorm=dd.norm();
	  if (retrieveGlobalRank()==0) cerr<<" appinv1: norm of muted data = "<<ddnorm<<endl;
	  lrop.applyAdjOp(dd,tmp);
	  float tmpnorm=tmp.norm();
	  if (retrieveGlobalRank()==0) cerr<<" appinv1: norm of adjoint output = "<<tmpnorm<<endl;
	  wm.applyInvOp(tmp,cm[1]);
	  float cmnorm=cm[1].norm();
	  if (retrieveGlobalRank()==0) cerr<<" appinv1: norm of wm output = "<<cmnorm<<endl;
	}
	else {
	  (*optr) << "start pcg\n";
	  RVL::InvLinearOp<float> iwm(wm);
	  TSOpt::GridRadialScaleOp hmult(lrop.getDomain(),valparse<float>(*pars,"alpha",0.0f));
	  RVL::RWAdjLinearOp<float> rwhmult(hmult,wm);
	  (*optr) << "build tensop\n";
	  //	  RVL::TensorLinearOp<float> tensop(cop,rwhmult);
	  RVL::TensorLinearOp<float> tensop(cop,hmult);
	  (*optr) << "build rhs\n";
	  RVL::Vector<float> td(tensop.getRange());
	  Components<float> ctd(td);
	  ctd[0].copy(d);
	  ctd[1].zero();
	  (*optr) <<"rhs norm = "<<td.norm()<<"\n";
	  //	  optr->flush();
	  float rtol=RVL::valparse<float>(*pars,"ResidualTol",100.0*numeric_limits<float>::epsilon());
	  float nrtol=RVL::valparse<float>(*pars,"GradientTol",100.0*numeric_limits<float>::epsilon());
	  float maxstep=RVL::valparse<float>(*pars,"MaxStep",numeric_limits<float>::max());
	  float rnorm;
	  float nrnorm;

	  (*optr) << "set up pcg\n";
	  RVLUmin::CGNEAlg<float> alg(cm[1],tensop,iwm,td,
				      rnorm, nrnorm, rtol, nrtol, maxcount, maxstep, *optr);
	  //	  optr->flush();
	  /*
	  (*optr) << "create local output\n";
	  RVL::Vector<float> r(tensop.getDomain());
	  RVL::AssignFilename afr("graw.rsf");
	  r.eval(afr);
	  RVL::Vector<float> g(tensop.getDomain());
	  RVL::AssignFilename afg("gpre.rsf");
	  g.eval(afg);
	  (*optr) <<"apply adjoint modeling op\n";
	  tensop.applyAdjOp(td,r);
	  (*optr) <<"apply precond op\n";
	  iwm.applyOp(r,g);
	  
	  (*optr)<<"res norm      = "<<td.norm()<<"\n";
	  (*optr)<<"raw grad norm = "<<r.norm()<<"\n";
	  (*optr)<<"pre grad norm = "<<g.inner(r)<<"\n";
	  */
	  (*optr)<<"run pcg\n";
	  alg.run();
	  // create, record residual components
	  RVL::Vector<float> atd(tensop.getRange());
	  Components<float> catd(atd);
	  std::string resnm = RVL::valparse<std::string>(*pars,"datares","");
	  if (resnm.size()>0) {
	    RVL::AssignFilename afres(resnm);
	    catd[0].eval(afres);
	  }
	  tensop.applyOp(cm[1],atd);
	  catd[0].linComb(-1.0f,d);
	  float e = catd[0].normsq();
	  RVL::Vector<float> schr(catd[1].getSpace());
	  RVL::Vector<float> hr(catd[1].getSpace());
	  //	  hmult.applyOp(cm[1],hr);
	  //	  wm.applyOp(hr,schr);
	  //	  float p = hr.inner(schr);
	  float p = catd[1].normsq();
	  // use xplus = (0.125*|d|)^2
	  float alpha = valparse<float>(*pars,"alpha",0.0f);
	  float alphaplussq = alpha*alpha + ((0.125*(d.norm()))*(0.125*(d.norm())) - e)/(2.0*p);
	  (*optr)<<"data residual = "<<catd[0].norm()<<"\n";
	  (*optr)<<"norm of h*r   = "<<catd[1].norm()<<"\n";
	  (*optr)<<"alpha         = "<<RVL::valparse<float>(*pars,"alpha",0.0f)<<"\n";
	  (*optr)<<"alpha_+^2     = "<<alphaplussq<<"\n";
	  if ((outfile.size()==0) && (retrieveGlobalRank()==0)) {
	    std::shared_ptr<std::stringstream> fish = dynamic_pointer_cast<std::stringstream>(optr);
	    fprintf(stream,fish->str().c_str());
	  }
	}
      }
      else{
	if (adjoint==0) {
	  cop.applyOp(cm[1],d);
	}
	else {
	  if (maxcount==0) {
	    cop.applyAdjOp(d,cm[1]);
	    cerr<<"migration output norm = "<<cm[1].norm()<<endl;
	  }
	  else {
	    TSOpt::GridRadialScaleOp hmult(lrop.getDomain(),valparse<float>(*pars,"alpha",0.0f));
	    RVL::TensorLinearOp<float> tensop(cop,hmult);	    
	    RVL::Vector<float> td(tensop.getRange());
	    Components<float> ctd(td);
	    ctd[0].copy(d);
	    ctd[1].zero();
	    float rtol=valparse<float>(*pars,"ResidualTol",100.0*numeric_limits<float>::epsilon());
	    float nrtol=valparse<float>(*pars,"GradientTol",100.0*numeric_limits<float>::epsilon());
	    float maxstep=valparse<float>(*pars,"MaxStep",numeric_limits<float>::max());
	    float rnorm;
	    float nrnorm;
	    RVLUmin::CGNEAlg<float> alg(cm[1],tensop,td,
					rnorm, nrnorm, rtol, nrtol, maxcount, maxstep, *optr);
	    alg.run();
	    if ((outfile.size()==0) && (retrieveGlobalRank()==0)) {
	      std::shared_ptr<std::stringstream> fish = dynamic_pointer_cast<std::stringstream>(optr);
	      fprintf(stream,fish->str().c_str());
	    }	  
	  }
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
