#include "iwop.hh"
#include "functions.hh"
#include "linop_apps.hh"
#include <par.h>
#include "segyops.hh"
#include "adjtest.hh"
#include "asg_defn.hh"
#include "sim_selfdoc.h"
#include "parser.h"
#include "blockop.hh"
#include "cgnealg.hh"
#include "alphaupdate.hh"
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
			      RVL::valparse<float>(*pars,"mute_width",0.0f));
      
      RVL::LinearOpFO<float> muteop(iwop.getIWaveRange(),iwop.getIWaveRange(),mute,mute);

      // taper
      TSOpt::SEGYTaper fo(RVL::valparse<std::string>(*pars,"taperpars",""));
      RVL::LinearOpFO<float> taperop(iwop.getRange(),
				     iwop.getRange(),
				     fo,fo);
    
      RVL::LinearRestrictOp<float> rlop(iwop,cm[0]);
      RVL::CompLinearOp<float> mrlop(rlop,muteop);
      RVL::CompLinearOp<float> tmrlop(mrlop,taperop);

      // trace scaling by weighted offset
      float alpha=RVL::valparse<float>(*pars,"weight",0.0f);
      float alphasq=alpha*alpha;
      TSOpt::TraceScaleFO ts(alpha,"offset");
      RVL::LinearOpFO<float> tsop(tmrlop.getDomain(),tmrlop.getDomain(),ts,ts);

      // block operator (rlop, tsop)^T
      RVL::TensorLinearOp<float> top(tmrlop,tsop);
      
      RVL::Vector<float> dd(top.getRange());
      RVL::Components<float> cdd(dd);

      // load data
      RVL::AssignFilename ddfn(RVL::valparse<std::string>(*pars,"data_p"));
      cdd[0].eval(ddfn);
      cdd[1].zero();

      //      top.applyAdjOp(dd,cm[1]);

      //      cerr<<"data norm="<<cdd[0].norm()<<endl;
      //      cerr<<"adj output norm="<<cm[1].norm()<<endl;

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

    res<<endl<<"*******************************************************"<<endl;
    res<<"* Acoustic Variable Density Source Inversion via";
    res<<"* Conjugate Gradient Algorithm for Normal Eqns"<<endl;
    res<<"* max iterations       = "<<maxcount<<endl;
    res<<"* residual tolerance   = "<<rtol<<endl;
    res<<"* normal res tolerance = "<<nrtol<<endl;
    res<<"* trust radius         = "<<maxstep<<endl;
    res<<"*******************************************************"<<endl;
    
    /* create CGNE object */
    float rnorm;
    float nrnorm;
    RVLUmin::CGNEAlg<float> alg(cm[1],top,dd,
				rnorm, nrnorm, rtol, nrtol, maxcount, maxstep, res);
    //    float nrnorm0=nrnorm;
    //    float rnorm0=rnorm;

    // run CG
    alg.run();
      
    std::string dataest = valparse<std::string>(*pars,"est_data_p","");
    std::string datares = valparse<std::string>(*pars,"res_data_p","");
    std::string dataann = valparse<std::string>(*pars,"res_anni_p","");
    if (dataest.size()>0) {
      Vector<float> est(top.getRange());
      Components<float> cest(est);
      AssignFilename estfn(dataest);
      cest[0].eval(estfn);
      if (dataann.size()>0) {
	AssignFilename annfn(dataann);
	cest[1].eval(annfn);
      }
      top.applyOp(cm[1],est);
      RVL::Vector<float> dres(rlop.getRange());
      if (datares.size()>0) {
	RVL::AssignFilename resfn(datares);
	dres.eval(resfn);
      }
      // S[c]f-d rather than other way round
      dres.copy(cest[0]);
      dres.linComb(-1.0f,cdd[0]);
      if (dataann.size()>0) {
	TSOpt::TraceScaleFO ts1(1.0f,"offset");
	RVL::LinearOpFO<float> tsop1(rlop.getDomain(),rlop.getDomain(),ts1,ts1);
	tsop1.applyOp(cm[1],cest[1]);
      }
      if (RVL::valparse<int>(*pars,"alphaupdate",0)) {
	float ubnd = RVL::valparse<float>(*pars,"upperbnd");
	float lbnd = RVL::valparse<float>(*pars,"lowerbnd");
	float decr = RVL::valparse<float>(*pars,"decrfac");
	float alphasqnext =
	  RVLUmin::alphaupdate(cdd[0],dres,cest[1],
			       alphasq,ubnd*ubnd,lbnd*lbnd,decr*decr);
	res<<"current alpha = "<<alpha<<"\n";
	float dn = cdd[0].norm();
	float en = dres.norm();
	res<<"data norm     = "<<dn<<"\n";
	res<<"res norm      = "<<en<<"\n";
	res<<"upper bnd res = "<<dn*ubnd<<"\n";
	res<<"lower bnd res = "<<dn*lbnd<<"\n";
	res<<"pen norm      = "<<cest[1].norm()<<"\n";
	res<<"updated alpha = "<<sqrt(alphasqnext)<<"\n";
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
      MPI_Barrier(MPI_COMM_WORLD);
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
