#include "born.hh"
//#define IWAVE_VERBOSE
//#define RECORD
namespace TSOpt {

  void BornIWaveLOVOp::apply0(const Vector<float> & x0,
			      const Vector<float> & x1, 
			      Vector<float> & y) const {

    try {

      // choice - free or absorbing surface
      std::shared_ptr<TSOpt::IWaveLOVOp> iwop;
      std::string prefix;
      if (freechoice) {
	iwop = iwopfree;
	prefix = "free";
      }
      else {
	iwop = iwopabsb;
	prefix = "absb";
      }
	
      // Full ref vector - workspace for pert
      RVL::Vector<float> m(iwop->getProductDomain()); // domain vector - all components
      RVL::Vector<float> dm(iwop->getNonLinDomain()); // pert vector - all components
      dm.zero();
      Components<float> cdm(dm);  // pert vector - all components
      Components<float> cx0(x0);  // Born ref vector - active components
      Components<float> cx1(x1);  // Born pert vector
      Components<float> cm(m);   // components of stored Full ref
      Components<float> cm0(cm[0]); // nonlinear components of stored Full ref
      Components<float> cm1(cm[1]); // nonlinear components of stored Full ref
      int ib = 0; // index for active components
      for (int j=0; j< iwop->getNonLinDomain().getSize(); j++) {
	std::string key = iwop->getNonLinDomain().getKeys()[j];
	if (j == idx[ib]) {
	  myGridSpace const & gsp = dynamic_cast<myGridSpace const &>(cx0[ib].getSpace());
	  myGridSpace const & egsp = dynamic_cast<myGridSpace const &>(cm0[idx[ib]].getSpace());
	  GridExtendOp gop(gsp,egsp);
	  gop.applyOp(cx0[ib],cm0[idx[ib]]);
#ifdef IWAVE_VERBOSE
	  if (retrieveGlobalRank()==0) {
	    cerr<<"extended nonlin input comp "<<ib<<" = iwop domain component "<<idx[ib]<<"\n";
	    cerr<<"   with key "<<key<<"\n";
	  }
#endif
	  ib++;
	}
	else {
	  std::string fn = RVL::valparse<std::string>(iwop->getPar(),key);
	  RVL::AssignFilename af(fn);
	  cm0[j].eval(af);
#ifdef IWAVE_VERBOSE
	  if (retrieveGlobalRank()==0) 
	    cerr<<"loaded "<<fn<<" on key "<<key<<" into iwop vector, nonlin component "<<j<<"\n";
#endif
	}
      }      

      // load source parameters - must use prefix!!!
      for (int j=0; j<iwop->getLinDomain().getSize(); j++) {
	std::string key = iwop->getLinDomain().getKeys()[j];
	std::string fn = RVL::valparse<std::string>(iwop->getPar(),prefix+key);
	RVL::AssignFilename af(fn);
	cm1[j].eval(af);
#ifdef IWAVE_VERBOSE
	if (retrieveGlobalRank()==0) 
	  cerr<<"loaded "<<fn<<" on key "<<prefix+key<<" into iwop vector, linear component "<<j<<"\n";
#endif
      }

#ifdef IWAVE_VERBOSE
      if (retrieveGlobalRank()==0) {
	cerr<<"\n ACTUAL PAR ARRAY FOR FORWARD SIMULATION\n";
	ps_printall(iwop->getPar(),cerr);
	cerr<<"\n NONLINEAR INPUTS:\n";
	cm[0].write(cerr);
	cerr<<"\n LINEAR INPUTS:\n";
	cm[1].write(cerr);
	cerr<<"\n PERT NONLIN INPUTS:\n";
	dm.write(cerr);
	cerr<<"\n OUTPUTS:\n";
	y.write(cerr);
      }
#endif
	  
      for (size_t i=0; i<bkeys.size(); i++) {	
	cdm[idx[i]].copy(cx1[i]);
      }
      RVL::LinOpValOp<float>::export_applyPartialDeriv0(*iwop,cm[0],cm[1],dm,y);

      // WWS 2017.12.15 shift integrations to adjoint
      // two trace integrations in appinv + free surface case
      //      if (appinv && freechoice) {
      //	SEGYFwdIntIP intop(2);
      //	y.eval(intop);
      //      }

    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: BornIWaveLOVOp::apply0\n";
      e<<"  failure to cast one of the domain spaces as a grid space\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from BornIWaveLOVOp::apply0\n";
      throw e;
    }

  }

  void BornIWaveLOVOp::applyAdj0(const Vector<float> & x0,
				 const Vector<float> & y, 
				 Vector<float> & x1) const {

    try {
      // choice - free or absorbing surface
#ifdef IWAVE_VERBOSE
      if (retrieveGlobalRank()==0) {
	cerr<<"BornIWaveLOVOp::applyAdj0\n";
	cerr<<"freechoice = "<<freechoice<<" appinv = "<<appinv<<"\n";
      }
#endif
      /*** free surface ***/
      if (freechoice) {

	/*** appinv ***/
	if (appinv) {
	  // 1. map data to zero-depth data - OpZFwd (includes scaling by s/r cellvol)
	  // 2. inegrate in time - WWS 2017.12.15 three times!
	  // 3. adjoint for absb
	  // 4. z-deriv
	  // 5. scale, map to free-surface model
#ifdef IWAVE_VERBOSE
	  float ynorm = y.norm();
	  if (retrieveGlobalRank()==0) {
	    cerr<<"biwop::applyAdj0 norm of input free surface data = "<<ynorm<<endl;
	    cerr<<"create zero-depth traces\n";
	  }
#endif
	  RVL::Vector<float> zdtraces(iwopabsb->getRange());

#ifdef RECORD
	  RVL::AssignFilename afzd("zdtraces.su");
	  zdtraces.eval(afzd);
#endif
	  
	  // 1.1 cell volume
	  SRCellVolFO sr;
	  RVL::MPISerialFunctionObjectRedn<float,float> mpisr(sr);
	  y.eval(mpisr);
#ifdef IWAVE_VERBOSE	  
	  if (retrieveGlobalRank()==0) 	  
	    cerr<<"derived s/r cell volume = "<<sr.getValue()<<"\n";
#endif
	  // 1.2 map data to zero-depth data
	  OpZFwdFO ozfwd;
	  RVL::MPISerialFunctionObject<float> mpiozfwd(ozfwd);	  
	  zdtraces.eval(mpiozfwd,y);
#ifdef IWAVE_VERBOSE
	  float zdnorm1=zdtraces.norm();
	  if (retrieveGlobalRank()==0)
	    cerr<<"biwop::applyAdj0 norm of unscaled zd data = "<<zdnorm1<<endl;
#endif
	  // 1.3 scale by MINUS (20170206) cell volume
	  TraceScaleFO trsc(-sr.getValue());
	  RVL::MPISerialFunctionObject<float> mpitrsc(trsc);	  	  
	  zdtraces.eval(mpitrsc);
#ifdef IWAVE_VERBOSE
	  float zdnorm2=zdtraces.norm();
	  if (retrieveGlobalRank()==0)
	    cerr<<"biwop::applyAdj0 norm of scaled zd data = "<<zdnorm2<<endl;
#endif
#ifdef IWAVE_VERBOSE
	  if (retrieveGlobalRank()==0) 	  
	    cerr<<"time integrations\n";
#endif

	  //WWS 20170206: hard-wire FT computations
	  //	  if (RVL::valparse<int>(iwopabsb->getPar(),"ft")) {
	  float maxfreq = RVL::valparse<float>(iwopabsb->getPar(),"maxfreq");
	  float power = -3.0f;
	  int ab = 1;
	  SEGYTFTScaleFO ftscale(power, ab,
				 RVL::valparse<float>(iwopabsb->getPar(),"Tlocut",
						      0.0625*maxfreq),
				 RVL::valparse<float>(iwopabsb->getPar(),"Tlopas",
						      0.125*maxfreq),
				 RVL::valparse<float>(iwopabsb->getPar(),"Thipas",
						      0.9375*maxfreq),
				 RVL::valparse<float>(iwopabsb->getPar(),"Thicut",
						      maxfreq));
	  RVL::MPISerialFunctionObject<float> mpiftscale(ftscale);	  	  
	  zdtraces.eval(mpiftscale);
	  
#ifdef IWAVE_VERBOSE
	  float zdnorm3=zdtraces.norm();
	  if (retrieveGlobalRank()==0)
	    cerr<<"biwop::applyAdj0 norm of integrated zd data = "<<zdnorm3<<endl;
#endif
	  /*	
		else {
		#ifdef IWAVE_VERBOSE	    
		if (retrieveGlobalRank()==0)
		cerr<<"ft = 0 ab = 0\n";
		#endif
		SEGYFwdIntIP intfwd(3);
		RVL::MPISerialFunctionObject<float> mpiintfwd(intfwd);	  	  	    
		zdtraces.eval(mpiintfwd);
		}
	  */
#ifdef IWAVE_VERBOSE
	  if (retrieveGlobalRank()==0) 	  
	    cerr<<"IWLOVOP call\n";
#endif
	  // 3.1 set up mapping by opposite iwop
	  // Full ref vector - workspace for pert
	  RVL::Vector<float> m(iwopabsb->getProductDomain()); // domain vector - all components
	  RVL::Vector<float> dm(iwopabsb->getNonLinDomain()); // pert vector - all components
	  dm.zero();
	  Components<float> cdm(dm);  // pert vector - all components
	  Components<float> cx0(x0);  // Born ref vector - active components
	  Components<float> cx1(x1);  // Born pert vector
	  Components<float> cm(m);    // components of stored Full ref
	  Components<float> cm0(cm[0]); // nonlinear components of stored Full ref
	  Components<float> cm1(cm[1]); // linear components of stored Full ref (sources)

	  // save vector for data weighted output
	  RVL::Vector<float> cx1save(cx1[0].getSpace());
	  RVL::AssignFilename afs("cx1save.rsf");
	  cx1save.eval(afs);
	  
	  int ib = 0; // index for active components
	  for (int j=0; j< iwopabsb->getNonLinDomain().getSize(); j++) {
	    std::string key = iwopabsb->getNonLinDomain().getKeys()[j];
	    if (j == idx[ib]) {
	      // How to get extension to work: first create intermediate extension using
	      // pert component as model - already extended, should be same extension
	      // the spray extension along 0 axis to obtain component of m[0] - assumption
	      // that difference is only spray along 0 axis
	      myGridSpace const & gsp = dynamic_cast<myGridSpace const &>(cx0[ib].getSpace());
	      myGridSpace const & egsp = dynamic_cast<myGridSpace const &>(cx1[ib].getSpace());
	      RVL::Vector<float> tmp(egsp);
	      GridExtendOp gop(gsp,egsp);
	      gop.applyOp(cx0[ib],tmp);
	      GridOverSprayFO glp(0);
	      RVL::MPISerialFunctionObject<float> mpiglp(glp);	  	  	    
	      cm0[idx[ib]].eval(mpiglp,tmp);
#ifdef IWAVE_VERBOSE
	      if (retrieveGlobalRank()==0) {
		cerr<<"extended nonlin input comp "<<ib<<" = iwop domain component "<<idx[ib]<<"\n";
		cerr<<"   with key "<<key<<"\n";
	      }
#endif
	      ib++;
	    }
	    else {
	      std::string fn = RVL::valparse<std::string>(iwopabsb->getPar(),key);
	      RVL::AssignFilename af(fn);
	      cm0[j].eval(af);
#ifdef IWAVE_VERBOSE
	      if (retrieveGlobalRank()==0) 
		cerr<<"loaded "<<fn<<" on key "<<key<<" into iwop vector, nonlin component "<<j<<"\n";
#endif
	    }
	  }
	  for (int j=0; j<iwopabsb->getLinDomain().getSize(); j++) {
	    std::string key = iwopabsb->getLinDomain().getKeys()[j];
	    std::string fn = RVL::valparse<std::string>(iwopabsb->getPar(),key);
	    RVL::AssignFilename af(fn);
	    cm1[j].eval(af);
#ifdef IWAVE_VERBOSE
	    if (retrieveGlobalRank()==0) 
	      cerr<<"loaded "<<fn<<" on key "<<key<<" into iwop vector, linear component "<<j<<"\n";
#endif
	  }
#ifdef IWAVE_VERBOSE
	  if (retrieveGlobalRank()==0) 
	    cerr<<"-->export_apply\n";
#endif
	  RVL::LinOpValOp<float>::export_applyAdjPartialDeriv0(*iwopabsb,cm[0],cm[1],zdtraces,dm);
#ifdef IWAVE_VERBOSE
	  if (retrieveGlobalRank()==0) 
	    cerr<<"<--export_apply\n";
#endif
	  
	  // 4.1 z deriv after map absb -> free model perts

	  GridOverSprayFO glp(0); // makes target equal to source on overlap, zero else
	  RVL::MPISerialFunctionObject<float> mpiglp(glp);	  	  	    	      
	  ib = 0; // index for active components
	  for (int j=0; j< iwopabsb->getNonLinDomain().getSize(); j++) {
	    if ((ib<idx.size()) && (j == idx[ib])) {
	      // also sprays, axis=arg
#ifdef IWAVE_VERBOSE
	      if (retrieveGlobalRank()==0) {	      
		cerr<<"ib="<<ib<<" j="<<j<<" rank="<<retrieveGlobalRank()<<"\n";	    	      
		cerr<<"-->overspray component "<<ib<<" j="<<j<<" rank="<<retrieveGlobalRank()<<"\n";
	      }
#endif
	      
	      //#ifndef __INTEGRAL__		  
	      cx1[ib].eval(mpiglp,cdm[j]);
#ifdef IWAVE_VERBOSE
	      if (retrieveGlobalRank()==0) 	      
		cerr<<"<--overspray component "<<ib<<" j="<<j<<" rank="<<retrieveGlobalRank()<<"\n";
#endif
	      if (ib==0) {
		cx1save.copy(cx1[ib]);
	      }
	      /*	      
#else	  
	      RVL::Vector<float> tmp(cx1[ib].getSpace());
	      tmp.zero();
	      tmp.eval(mpiglp,cdm[j]);
	      if (ib==0) {
		cx1save.copy(tmp);
	      }
#ifdef IWAVE_VERBOSE
	      if (retrieveGlobalRank()==0) 		
		cerr<<"differentiate sprayed adj output, store in x1 comp\n";
#endif
	      //	      if (RVL::valparse<int>(iwopabsb->getPar(),"ft")) {
	      float power = 1.0;
	      float maxmaxspf = RVL::valparse<float>(iwopabsb->getPar(),"maxfreq")/
		RVL::valparse<float>(iwopabsb->getPar(),"cmin");
	      float minmaxspf = RVL::valparse<float>(iwopabsb->getPar(),"maxfreq")/
		RVL::valparse<float>(iwopabsb->getPar(),"cmax");
	      TSOpt::GridZFTScaleOp gop(cx1[ib].getSpace(),
					power,
					RVL::valparse<int>(iwopabsb->getPar(),"ab"),
					RVL::valparse<float>(iwopabsb->getPar(),"Zlocut",
							     0.0625*minmaxspf),
					RVL::valparse<float>(iwopabsb->getPar(),"Zlopas",
							     0.125*minmaxspf),
					RVL::valparse<float>(iwopabsb->getPar(),"Zhipas",
							     0.9375*maxmaxspf),
					RVL::valparse<float>(iwopabsb->getPar(),"Zhicut",
							     maxmaxspf));
#ifdef IWAVE_VERBOSE
	      if (retrieveGlobalRank()==0) 				  
		cerr  <<" locut="<< RVL::valparse<float>(iwopabsb->getPar(),"Zlocut",
							 0.0625*minmaxspf)
		      <<" lopas="<< RVL::valparse<float>(iwopabsb->getPar(),"Zlopas",
							 0.125*minmaxspf)
		      <<" hipas="<< RVL::valparse<float>(iwopabsb->getPar(),"Zhipas",
							 0.9375*maxmaxspf)
		      <<" hicut="<< RVL::valparse<float>(iwopabsb->getPar(),"Zhicut",
							 maxmaxspf) <<"\n";
#endif
	      gop.applyOp(tmp,cx1[ib]);

	      else {
		TSOpt::GridZDerivOp gop(cx1[ib].getSpace());
		gop.applyOp(tmp,cx1[ib]);
	      }

	      // 5.1 weight by powers of coeffs - active from x0, passive from m
	      // purely acoustic choices - maybe should be made elsewhere
	      std::vector<float> pwr(iwopabsb->getNonLinDomain().getSize());
	      std::vector<float> fac(iwopabsb->getNonLinDomain().getSize());
	      // this is the explicitly acoustic choice
	      // bulkmod 2.5 pwr - gives rho^2.5 v^5
	      // then 1.5 pwr of beta gives rho v^5 = kappa v^3
	      // factor of 8 included in second scale FO
	      pwr[0]=3; pwr[1]=1;
	      fac[0]=1.0; fac[1]=32.0;
	      // the rest of this code is general, assuming that the correct
	      // scaling is a product of powers of coeffs
	      int ib1=0;
	      for (int k=0; k< iwopabsb->getNonLinDomain().getSize(); k++) {
		GridScaleFO gsc(fac[k],pwr[k]);
		RVL::MPISerialFunctionObject<float> mpigsc(gsc);
#ifdef IWAVE_VERBOSE
		if (retrieveGlobalRank()==0)            
		  cerr<<"scale output "<<ib<<" by coefficient field"<<k<<"\n";
#endif
		if ((ib1<idx.size()) && (k == idx[ib1])) {
		  cx1[ib].eval(mpigsc,cx0[ib1]);
		  ib1++;
		}
		else {
		  // different grid, have to spray first - re-use tmp
		  // grid is same as cx0[ib] grid - since cx1[ib] grid is
		  // extension - this avoids loading the free surface grids
		  // explicitly
		  RVL::Vector<float> tmp(cx0[ib].getSpace());
		  tmp.eval(mpiglp,cm0[k]);
		  cx1[ib].eval(mpigsc,tmp);
		}
	      }
#endif
	      */
	      ib++;
	    }
	  }
#ifdef IWAVE_VERBOSE
	  if (retrieveGlobalRank()==0) 	  
	    cerr<<"end free surface appinv\n";
#endif
	}

	/*** end appinv ***/

	/*** migration ***/
	else {
	  
	  // Full ref vector - workspace for pert
	  RVL::Vector<float> m(iwopfree->getProductDomain()); // domain vector - all components
	  RVL::Vector<float> dm(iwopfree->getNonLinDomain()); // pert vector - all components
	  dm.zero();
	  Components<float> cdm(dm);  // pert vector - all components
	  Components<float> cx0(x0);  // Born ref vector - active components
	  Components<float> cx1(x1);  // Born pert vector
	  Components<float> cm(m);    // components of stored Full ref
	  Components<float> cm0(cm[0]); // nonlinear components of stored Full ref
	  Components<float> cm1(cm[1]); // linear components of stored Full ref (sources)
	  int ib = 0; // index for active components
	  for (int j=0; j< iwopfree->getNonLinDomain().getSize(); j++) {
	    std::string key = iwopfree->getNonLinDomain().getKeys()[j];
	    if (j == idx[ib]) {
	      myGridSpace const & gsp = dynamic_cast<myGridSpace const &>(cx0[ib].getSpace());
	      myGridSpace const & egsp = dynamic_cast<myGridSpace const &>(cm0[idx[ib]].getSpace());
	      GridExtendOp gop(gsp,egsp);
	      gop.applyOp(cx0[ib],cm0[idx[ib]]);
#ifdef IWAVE_VERBOSE
	      cerr<<"extended nonlin input comp "<<ib<<" = iwop domain component "<<idx[ib]<<"\n";
	      cerr<<"   with key "<<key<<"\n";
#endif
	      ib++;
	    }
	    else {
	      std::string fn = RVL::valparse<std::string>(iwopfree->getPar(),key);
	      RVL::AssignFilename af(fn);
	      cm0[j].eval(af);
#ifdef IWAVE_VERBOSE
	      cerr<<"loaded "<<fn<<" on key "<<key<<" into iwop vector, nonlin component "<<j<<"\n";
#endif
	    }
	  }
	  for (int j=0; j<iwopfree->getLinDomain().getSize(); j++) {
	    std::string key = iwopfree->getLinDomain().getKeys()[j];
	    std::string fn = RVL::valparse<std::string>(iwopfree->getPar(),key);
	    RVL::AssignFilename af(fn);
	    cm1[j].eval(af);
#ifdef IWAVE_VERBOSE	  
	    cerr<<"loaded "<<fn<<" on key "<<key<<" into iwop vector, linear component "<<j<<"\n";
#endif
	  }
	    
	  RVL::LinOpValOp<float>::export_applyAdjPartialDeriv0(*iwopfree,cm[0],cm[1],y,dm);
	    
	  // copy works here - same grid
	    
	  for (size_t i=0; i<bkeys.size(); i++) {	
	    cx1[i].copy(cdm[idx[i]]);
	  }
#ifdef IWAVE_VERBOSE
	  cerr<<"end free surface migration\n";
#endif
	}
	/*** end migration ***/
#ifdef IWAVE_VERBOSE
	cerr<<"end free surface case\n";
#endif
      }
      /*** end free surface ***/
    
      /*** absorbing boundary ***/ 
      else {

	/*** appinv ***/
	if (appinv) {
	  // 1. integrate in time
	  // 2. map data to tow-depth data - OpZAdj 
	  // 3. adjoint for free - includes two time integrals!!
	  // 4. z-deriv
	  // 5. scale, map to abs-surface model
	  // for convenience communte time integration and map to tow depth
#ifdef IWAVE_VERBOSE
	  if (retrieveGlobalRank()==0) 
	    cerr<<"create tow-depth traces\n";
#endif
	  RVL::Vector<float> tdtraces(iwopfree->getRange());
#ifdef RECORD	  
	  RVL::Vector<float> stage1(iwopfree->getRange());
	  RVL::AssignFilename fn1("stage1.su");
	  stage1.eval(fn1);
	  stage1.copy(tdtraces);
#endif	  
	  // 1.2 map data to tow-depth data
	  OpZAdjFO ozadj;
	  RVL::MPISerialFunctionObject<float> mpiozadj(ozadj);	  	  
	  tdtraces.eval(mpiozadj,y);
#ifdef RECORD	  
	  stage1.copy(tdtraces);
#endif	  
	  // 1.1 cell volume
	  SRCellVolFO sr;
	  RVL::MPISerialFunctionObjectRedn<float,float> mpisr(sr);
	  y.eval(mpisr);
#ifdef IWAVE_VERBOSE
	  if (retrieveGlobalRank()==0) 	  
	    cerr<<"derived s/r cell volume = "<<sr.getValue()<<"\n";
#endif
	  // 1.3 scale by MINUS cell volume (2017.02.06)
	  TraceScaleFO trsc(-sr.getValue());
	  RVL::MPISerialFunctionObject<float> mpitrsc(trsc);	  	  
	  tdtraces.eval(mpitrsc);
#ifdef RECORD
	  cerr<<"scaled by cell volume\n";
	  RVL::Vector<float> stage2(iwopfree->getRange());
	  RVL::AssignFilename fn2("stage2.su");
	  stage2.eval(fn2);
	  stage2.copy(tdtraces);
#endif
#ifdef IWAVE_VERBOSE
	  if (retrieveGlobalRank()==0) 	  
	    cerr<<"time integration - includes two for F_1\n";
#endif
	  //	  if (RVL::valparse<int>(iwopfree->getPar(),"ft")) {
	  float maxfreq = RVL::valparse<float>(iwopfree->getPar(),"maxfreq");
	  float power = -3.0f;
	  int ab=1;
	  SEGYTFTScaleFO ftscale(power, ab,
				 RVL::valparse<float>(iwopfree->getPar(),"Tlocut",
						      0.0625*maxfreq),
				 RVL::valparse<float>(iwopfree->getPar(),"Tlopas",
						      0.125*maxfreq),
				 RVL::valparse<float>(iwopfree->getPar(),"Thipas",
						      0.9375*maxfreq),
				 RVL::valparse<float>(iwopfree->getPar(),"Thicut",
						      maxfreq));
	  RVL::MPISerialFunctionObject<float> mpiftscale(ftscale);	  	    	    
	  tdtraces.eval(mpiftscale);
	  /*
	  }
	  else {	  
	    SEGYFwdIntIP intfwd(3);
	    RVL::MPISerialFunctionObject<float> mpiintfwd(intfwd);	  	    	    
	    tdtraces.eval(mpiintfwd);
	  }
	  */
#ifdef RECORD
	  RVL::Vector<float> stage3(iwopfree->getRange());
	  RVL::AssignFilename fn3("stage3.su");
	  stage3.eval(fn3);
	  stage3.copy(tdtraces);
#endif
#ifdef IWAVE_VERBOSE
	  if (retrieveGlobalRank()==0) 	  
	    cerr<<"IWLOVOP call\n";
#endif
	  // 3.1 set up mapping by opposite iwop
	  // Full ref vector - workspace for pert
	  RVL::Vector<float> m(iwopfree->getProductDomain()); // domain vector - all components
	  RVL::Vector<float> dm(iwopfree->getNonLinDomain()); // pert vector - all components
	  dm.zero();
	  Components<float> cdm(dm);  // pert vector - all components
	  Components<float> cx0(x0);  // Born ref vector - active components
	  Components<float> cx1(x1);  // Born pert vector
	  Components<float> cm(m);    // components of stored Full ref
	  Components<float> cm0(cm[0]); // nonlinear components of stored Full ref
	  Components<float> cm1(cm[1]); // linear components of stored Full ref (sources)
	  int ib = 0; // index for active components
	  GridOverSprayFO glp(0);
	  RVL::MPISerialFunctionObject<float> mpiglp(glp);	  	    
	  for (int j=0; j< iwopfree->getNonLinDomain().getSize(); j++) {
	    std::string key = iwopfree->getNonLinDomain().getKeys()[j];
	    if (j == idx[ib]) {
	      // How to get extension to work: first create intermediate extension using
	      // pert component as model - already extended, should be same extension
	      // the spray extension along 0 axis to obtain component of m[0] - assumption
	      // that difference is only spray along 0 axis
	      myGridSpace const & gsp = dynamic_cast<myGridSpace const &>(cx0[ib].getSpace());
	      myGridSpace const & egsp = dynamic_cast<myGridSpace const &>(cx1[ib].getSpace());
	      RVL::Vector<float> tmp(egsp);
	      GridExtendOp gop(gsp,egsp);
	      gop.applyOp(cx0[ib],tmp);
	      cm0[idx[ib]].eval(mpiglp,tmp);
#ifdef IWAVE_VERBOSE
	      cerr<<"extended nonlin input comp "<<ib<<" = iwop domain component "<<idx[ib]<<"\n";
	      cerr<<"   with key "<<key<<"\n";
#endif
	      ib++;
	    }
	    else {
	      std::string fn = RVL::valparse<std::string>(iwopfree->getPar(),key);
	      RVL::AssignFilename af(fn);
	      cm0[j].eval(af);
#ifdef IWAVE_VERBOSE
	      cerr<<"loaded "<<fn<<" on key "<<key<<" into iwop vector, nonlin component "<<j<<"\n";
#endif
	    }
	  }
	  for (int j=0; j<iwopfree->getLinDomain().getSize(); j++) {
	    std::string key = iwopfree->getLinDomain().getKeys()[j];
	    std::string fn = RVL::valparse<std::string>(iwopfree->getPar(),key);
	    RVL::AssignFilename af(fn);
	    cm1[j].eval(af);
#ifdef IWAVE_VERBOSE	  
	    cerr<<"loaded "<<fn<<" on key "<<key<<" into iwop vector, linear component "<<j<<"\n";
#endif
	  }
	    
	  RVL::LinOpValOp<float>::export_applyAdjPartialDeriv0(*iwopfree,cm[0],cm[1],tdtraces,dm);

	  // 4.1 z deriv after map free -> absb model perts

#ifdef RECORD
	  {
	    RVL::Vector<float> stage4(cdm[0].getSpace());
	    RVL::AssignFilename fn4("stage4.rsf");
	    stage4.eval(fn4);
	    stage4.copy(cdm[0]);
	  }
#endif
	  
	  ib = 0; // index for active components
	  for (int j=0; j< iwopfree->getNonLinDomain().getSize(); j++) {
	    
	    //#ifndef __INTEGRAL__		  
	    if ((ib<idx.size()) && (j == idx[ib])) {
	      GridOverSprayFO glp(0); // makes target equal to source on overlap, zero else
	      // also sprays, axis=arg
#ifdef IWAVE_VERBOSE
	      if (retrieveGlobalRank()==0) 		
		cerr<<"overspray component "<<ib<<"\n";
#endif
	      cx1[ib].eval(mpiglp,cdm[j]);
	      ib++;
	    }

	    /*	    
#else	  		
	    if ((ib<idx.size()) && (j == idx[ib])) {
	      GridOverSprayFO glp(0); // makes target equal to source on overlap, zero else
	      // also sprays, axis=arg
#ifdef IWAVE_VERBOSE
	      if (retrieveGlobalRank()==0) 		
		cerr<<"overspray component "<<ib<<"\n";	      
#endif
	      RVL::Vector<float> tmp(cx1[ib].getSpace());	      
	      tmp.zero();
	      tmp.eval(mpiglp,cdm[j]);
#ifdef IWAVE_VERBOSE
	      if (retrieveGlobalRank()==0) 		
		cerr<<"differentiate sprayed adj output, store in x1 comp\n";
#endif
	      if (RVL::valparse<int>(iwopfree->getPar(),"ft")) {
		float power = 1.0;
		float maxspf = RVL::valparse<float>(iwopfree->getPar(),"maxfreq")/
		  RVL::valparse<float>(iwopfree->getPar(),"cmax");
		TSOpt::GridZFTScaleOp gop(cx1[ib].getSpace(),
					  power,
					  RVL::valparse<int>(iwopfree->getPar(),"ab"),
					  RVL::valparse<float>(iwopfree->getPar(),"Zlocut",
							       0.0625*maxspf),
					  RVL::valparse<float>(iwopfree->getPar(),"Zlopas",
							       0.125*maxspf),
					  RVL::valparse<float>(iwopfree->getPar(),"Zhipas",
							       0.9375*maxspf),
					  RVL::valparse<float>(iwopfree->getPar(),"Zhicut",
							       maxspf));
		gop.applyOp(tmp,cx1[ib]);
	      }
	      else {		
		TSOpt::GridZDerivOp gop(cx1[ib].getSpace());
		gop.applyOp(tmp,cx1[ib]);
	      }
#ifdef RECORD
	      {
		RVL::Vector<float> stage5(cx1[ib].getSpace());		
		RVL::AssignFilename fn5("stage5.rsf");
		stage5.eval(fn5);
		stage5.copy(cx1[ib]);
	      }
#endif
	      // 5.1 weight by powers of coeffs - active from x0, passive from m
	      // purely acoustic choices - maybe should be made elsewhere
	      std::vector<float> pwr(iwopfree->getNonLinDomain().getSize());
	      std::vector<float> fac(iwopfree->getNonLinDomain().getSize());
	      // this is the explicitly acoustic choice
	      // bulkmod 2.5 pwr - gives rho^2.5 v^5
	      // then 1.5 pwr of beta gives rho v^5 = kappa v^3
	      // factor of 8 included in second scale FO
	      pwr[0]=3.0; pwr[1]=1.0;
	      fac[0]=1.0; fac[1]=32.0;
	      // the rest of this code is general, assuming that the correct
	      // scaling is a product of powers of coeffs
	      int ib1=0;
	      for (int k=0; k< iwopfree->getNonLinDomain().getSize(); k++) {
		GridScaleFO gsc(fac[k],pwr[k]);
		GridOverSprayFO glp(0);
		RVL::MPISerialFunctionObject<float> mpigsc(gsc);
#ifdef IWAVE_VERBOSE
		if (retrieveGlobalRank()==0) 		
		  cerr<<"scale output "<<ib<<" by coefficient field"<<k<<"\n";
#endif
		if ((ib1<idx.size()) && (k == idx[ib1])) {
		  cx1[ib].eval(mpigsc,cx0[ib1]);
		  ib1++;
		}
		else {
		  // different grid, have to spray first - re-use tmp
		  // grid is same as cx0[ib] grid - since cx1[ib] grid is
		  // extension - this avoids loading the free surface grids
		  // explicitly
		  RVL::Vector<float> tmp(cx0[ib].getSpace());
		  tmp.eval(mpiglp,cm0[k]);
		  cx1[ib].eval(mpigsc,tmp);
		}
	      }
	      ib++;	      
	    }
#endif
	    */
	  }
	}
	/*** end appinv ***/

	/*** migration ***/
	else {
	  
	  // Full ref vector - workspace for pert
	  RVL::Vector<float> m(iwopabsb->getProductDomain()); // domain vector - all components
	  RVL::Vector<float> dm(iwopabsb->getNonLinDomain()); // pert vector - all components
	  dm.zero();
	  Components<float> cdm(dm);  // pert vector - all components
	  Components<float> cx0(x0);  // Born ref vector - active components
	  Components<float> cx1(x1);  // Born pert vector
	  Components<float> cm(m);    // components of stored Full ref
	  Components<float> cm0(cm[0]); // nonlinear components of stored Full ref
	  Components<float> cm1(cm[1]); // linear components of stored Full ref (sources)
	  int ib = 0; // index for active components
	  for (int j=0; j< iwopabsb->getNonLinDomain().getSize(); j++) {
	    std::string key = iwopabsb->getNonLinDomain().getKeys()[j];
	    if (j == idx[ib]) {
	      myGridSpace const & gsp = dynamic_cast<myGridSpace const &>(cx0[ib].getSpace());
	      myGridSpace const & egsp = dynamic_cast<myGridSpace const &>(cm0[idx[ib]].getSpace());
	      GridExtendOp gop(gsp,egsp);
	      gop.applyOp(cx0[ib],cm0[idx[ib]]);
#ifdef IWAVE_VERBOSE
	      cerr<<"extended nonlin input comp "<<ib<<" = iwop domain component "<<idx[ib]<<"\n";
	      cerr<<"   with key "<<key<<"\n";
#endif
	      ib++;
	    }
	    else {
	      std::string fn = RVL::valparse<std::string>(iwopabsb->getPar(),key);
	      RVL::AssignFilename af(fn);
	      cm0[j].eval(af);
#ifdef IWAVE_VERBOSE
	      cerr<<"loaded "<<fn<<" on key "<<key<<" into iwop vector, nonlin component "<<j<<"\n";
#endif
	    }
	  }
	  for (int j=0; j<iwopabsb->getLinDomain().getSize(); j++) {
	    std::string key = iwopabsb->getLinDomain().getKeys()[j];
	    std::string fn = RVL::valparse<std::string>(iwopabsb->getPar(),key);
	    RVL::AssignFilename af(fn);
	    cm1[j].eval(af);
#ifdef IWAVE_VERBOSE	  
	    cerr<<"loaded "<<fn<<" on key "<<key<<" into iwop vector, linear component "<<j<<"\n";
#endif
	  }
	  
	  RVL::LinOpValOp<float>::export_applyAdjPartialDeriv0(*iwopabsb,cm[0],cm[1],y,dm);
	  
	  // copy works here - same grid
	  
	  for (size_t i=0; i<bkeys.size(); i++) {	
	    cx1[i].copy(cdm[idx[i]]);
	  }
	}
	/*** end migration ***/
	
      }
      /*** end absorbing surface ***/
#ifdef IWAVE_VERBOSE
      cerr<<"end appinv code\n";
#endif
    }
  
    catch (bad_cast) {
      RVLException e;
      e<<"Error: BornIWaveLOVOp::applyAdj0\n";
      e<<"  failure to cast one of the domain spaces as a grid space\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from BornIWaveLOVOp::applyAdj0\n";
      throw e;
    }
  }
      
  void BornIWaveLOVOp::applyPartialDeriv0(const Vector<float> & x0,
					  const Vector<float> & x1,
					  const Vector<float> & dx0,
					  Vector<float> & dy) const {} 
    
  void BornIWaveLOVOp::applyAdjPartialDeriv0(const Vector<float> & x0,
					     const Vector<float> & x1,
					     const Vector<float> & dy,
					     Vector<float> & dx0) const {}
    
  BornIWaveLOVOp::BornIWaveLOVOp(std::vector<std::string> _bkeys, 
				 PARARRAY pars,
				 FILE * stream)
    : bkeys(_bkeys), freedom(2), absbdom(2), idx(0) {

    // free vs. absb surface for primary problem
    freechoice = valparse<int>(pars,"freesurface");
    std::string prefix;
    if (freechoice) prefix="free";
    else prefix="absb";

    // appinv vs. Euclidean IP in defns of operators
    appinv =  valparse<int>(pars,"appinv");

    // modify par array to create extended domain on chosen keys
    PARARRAY * freepars = ps_new();
    PARARRAY * absbpars = ps_new();      
    ps_copy(&freepars, pars);
    ps_copy(&absbpars, pars);

    // pml parameters
    float pmlwidth = RVL::valparse<float>(pars,"pmlwidth");
    ps_slfloat(*freepars,"nl1",0.0f);      
    ps_slfloat(*freepars,"nr1",pmlwidth);            
    ps_slfloat(*freepars,"nl2",pmlwidth);      
    ps_slfloat(*freepars,"nr2",pmlwidth);            
    ps_slfloat(*freepars,"nl3",pmlwidth);      
    ps_slfloat(*freepars,"nr3",pmlwidth);            
    ps_slfloat(*absbpars,"nl1",pmlwidth);      
    ps_slfloat(*absbpars,"nr1",pmlwidth);            
    ps_slfloat(*absbpars,"nl2",pmlwidth);      
    ps_slfloat(*absbpars,"nr2",pmlwidth);            
    ps_slfloat(*absbpars,"nl3",pmlwidth);      
    ps_slfloat(*absbpars,"nr3",pmlwidth);

    // key lists
    //    IWaveInfo ic;
    // find nonlin keywords - all essential
    std::vector<std::string> pkeys;
    for (int i=0; i< ic.get_num_iokeys(); i++) {
      if ((ic.iwave_iokeys[i].input > 0) &&
	  (ic.iwave_iokeys[i].active == 1) ) {
	pkeys.push_back(ic.iwave_iokeys[i].keyword);
      }
    }

    // find data (output active linear) keywords present in par file
    std::vector<std::string> dkeys;
    for (int i=0; i< ic.get_num_iokeys(); i++) {
      if ((ic.iwave_iokeys[i].input == 0) &&
	  (ic.iwave_iokeys[i].active == 2) ) {
	std::string tmp = RVL::valparse<std::string>(pars,prefix+ic.iwave_iokeys[i].keyword,"");
	if (tmp.size()>0) {
	  dkeys.push_back(ic.iwave_iokeys[i].keyword);
	}
      }
    }
    if (dkeys.size()<1) {
      RVLException e;
      e<<"Error: dpc\n";
      e<<"  must have some data keys\n";
      throw e;
    }
      
    // find source (input active linear) keywords present in par file
    std::vector<std::string> skeys;
    for (int i=0; i< ic.get_num_iokeys(); i++) {
      if ((ic.iwave_iokeys[i].input == 1) &&
	  (ic.iwave_iokeys[i].active == 2) ) {
	std::string tmp = RVL::valparse<std::string>(pars,prefix+ic.iwave_iokeys[i].keyword,"");
	if (tmp.size()>0) {
	  skeys.push_back(ic.iwave_iokeys[i].keyword);
	}
      }
    }
    if (skeys.size()<1) {
      RVLException e;
      e<<"Error: dpc\n";
      e<<"  must have some source keys\n";
      throw e;
    }
            
    // background parameters - replace active params by corresponding pert param key = bc+key+pert
    // as proxy for extended model, replace inactive params by corresponding params in each
    // bc-dependent pararray
    std::string key;
    std::string val;
    int ib=0;
    for (size_t i=0; i<pkeys.size(); i++) {
      if (pkeys[i]==bkeys[ib]) {
	key = "free"+bkeys[ib]+"pert";
	val = RVL::valparse<std::string>(pars,key);
	ps_slcstring(*freepars,bkeys[ib].c_str(),val.c_str());
#ifdef IWAVE_VERBOSE
	if (retrieveGlobalRank()==0) 	
	  cerr<<"assigned "<<val<<" to key "<<bkeys[ib]<<" in free surface IWOP par array\n"<<endl;
#endif
	key = "absb"+bkeys[ib]+"pert";
	val = RVL::valparse<std::string>(pars,key);
	ps_slcstring(*absbpars,bkeys[ib].c_str(),val.c_str());
#ifdef IWAVE_VERBOSE
	if (retrieveGlobalRank()==0) 
	  cerr<<"assigned "<<val<<" to key "<<bkeys[ib]<<" in absb surface IWOP par array\n"<<endl;
#endif
	ib++;
      }
      else {
	key = "free"+pkeys[i];
	val = RVL::valparse<std::string>(pars,key);
	ps_slcstring(*freepars,pkeys[i].c_str(),val.c_str());
#ifdef IWAVE_VERBOSE
	if (retrieveGlobalRank()==0) 
	  cerr<<"assigned "<<val<<" to key "<<pkeys[i]<<" in free surface IWOP par array\n"<<endl;
#endif
	key = "absb"+pkeys[i];
	val = RVL::valparse<std::string>(pars,key);
	ps_slcstring(*absbpars,pkeys[i].c_str(),val.c_str());
#ifdef IWAVE_VERBOSE
	if (retrieveGlobalRank()==0) 
	  cerr<<"assigned "<<val<<" to key "<<pkeys[i]<<" in absb surface IWOP par array\n"<<endl;
#endif
      }
    }
      
    // data parameters
    for (size_t i=0; i<dkeys.size(); i++) {
      key = "free"+dkeys[i];
      val = RVL::valparse<std::string>(pars,key);
      ps_slcstring(*freepars,dkeys[i].c_str(),val.c_str());
#ifdef IWAVE_VERBOSE
      if (retrieveGlobalRank()==0) 
	cerr<<"assigned "<<val<<" to key "<<dkeys[i]<<" in free surface IWOP par array\n"<<endl;
#endif
      key = "absb"+dkeys[i];
      val = RVL::valparse<std::string>(pars,key);
      ps_slcstring(*absbpars,dkeys[i].c_str(),val.c_str());
#ifdef IWAVE_VERBOSE
      if (retrieveGlobalRank()==0) 
	cerr<<"assigned "<<val<<" to key "<<dkeys[i]<<" in absb surface IWOP par array\n"<<endl;
#endif	
    }

    // source parameters
    for (size_t i=0; i<skeys.size(); i++) {
      key = "free"+skeys[i];
      val = RVL::valparse<std::string>(pars,key);
      ps_slcstring(*freepars,skeys[i].c_str(),val.c_str());
#ifdef IWAVE_VERBOSE
      if (retrieveGlobalRank()==0) 
	cerr<<"assigned "<<val<<" to key "<<skeys[i]<<" in free surface IWOP par array\n"<<endl;
#endif
      key = "absb"+skeys[i];
      val = RVL::valparse<std::string>(pars,key);
      ps_slcstring(*absbpars,skeys[i].c_str(),val.c_str());
#ifdef IWAVE_VERBOSE
      if (retrieveGlobalRank()==0) 
	cerr<<"assigned "<<val<<" to key "<<skeys[i]<<" in absb surface IWOP par array\n"<<endl;
#endif	
    }      
      
    // basic op which now has extended grids (assuming that perturbation had it)
    iwopfree = make_shared<TSOpt::IWaveLOVOp>(*freepars,stream);
    iwopabsb = make_shared<TSOpt::IWaveLOVOp>(*absbpars,stream);

    // create domain vector, initialize from modpars - non-extended params treated
    // as const
    // presumption is that original pars has filenames for non-extended params in
    // the extended list, defining non-extended geometry
    StdProductSpace<float> freenldom(bkeys.size());
    StdProductSpace<float> freelndom(bkeys.size());
    StdProductSpace<float> absbnldom(bkeys.size());
    StdProductSpace<float> absblndom(bkeys.size());
    for (int j=0; j< iwopfree->getNonLinDomain().getSize(); j++) {
      key = iwopfree->getNonLinDomain().getKeys()[j];
      // IWOP nonlin keys now assigned to extended pert fields, use to define
      // linear factor of domain
      std::string fn = RVL::valparse<std::string>(*freepars,key);
      for (size_t i=0; i<bkeys.size(); i++) {
	if (key == bkeys[i]) {
	  myGridSpace egsp(fn,"notype",true);
	  freelndom.set(egsp,i);
#ifdef IWAVE_VERBOSE
	  if (retrieveGlobalRank()==0) 	  
	    cerr<<"  constructed factor "<<i<<" in free surface linear dom space from file "<<fn<<"\n";
#endif
	  // nonlin keys with prefix point to correct prototypes for
	  // nonlinear factor of domain
	  fn=RVL::valparse<std::string>(pars,"free"+key);
	  myGridSpace gsp(fn,"notype",true);
	  freenldom.set(gsp,i);
#ifdef IWAVE_VERBOSE
	  if (retrieveGlobalRank()==0) 	  
	    cerr<<"  constructed factor "<<i<<" in free surface nonlin dom space from file "<<fn<<"\n";
#endif
	  idx.push_back(i);
#ifdef IWAVE_VERBOSE
	  if (retrieveGlobalRank()==0) 	  
	    cerr<<"  domain factor "<<i<<" -> iwop domain factor "<<j<<"\n";
#endif
	}
      }
    }
    for (int j=0; j< iwopabsb->getNonLinDomain().getSize(); j++) {
      key = iwopabsb->getNonLinDomain().getKeys()[j];
      std::string fn = RVL::valparse<std::string>(*absbpars,key);
      for (size_t i=0; i<bkeys.size(); i++) {
	if (key == bkeys[i]) {
	  myGridSpace egsp(fn,"notype",true);
	  absblndom.set(egsp,i);
#ifdef IWAVE_VERBOSE
	  cerr<<"  constructed factor "<<i<<" in absb surface linear dom space from file "<<fn<<"\n";
#endif
	  fn=RVL::valparse<std::string>(pars,"absb"+key);
	  myGridSpace gsp(fn,"notype",true);
	  absbnldom.set(gsp,i);
#ifdef IWAVE_VERBOSE
	  cerr<<"  constructed factor "<<i<<" in absb surface nonlin dom space from file "<<fn<<"\n";
#endif
	  //	    idx.push_back(i);
#ifdef IWAVE_VERBOSE
	  cerr<<"  domain factor "<<i<<" -> iwop domain factor "<<j<<"\n";
#endif
	}
      }
    }

    if (idx.size() != bkeys.size()) {
      RVLException e;
      e<<"Error: BornIWaveLOVOp constructor\n";
      e<<"  at least one extended parameter key not found in IWOP extended \n";
      e<<"  param list.\n";
      throw e;
    }

    // dom is product of non-extended, extended spaces
    freedom.set(freenldom,0);
    freedom.set(freelndom,1);
    absbdom.set(absbnldom,0);      
    absbdom.set(absblndom,1);

    ps_delete(&freepars);
    ps_delete(&absbpars);
  }

  BornIWaveLOVOp::BornIWaveLOVOp(BornIWaveLOVOp const & vop)
    : bkeys(vop.bkeys),
      iwopfree(make_shared<TSOpt::IWaveLOVOp>(*(vop.iwopfree))),
      iwopabsb(make_shared<TSOpt::IWaveLOVOp>(*(vop.iwopabsb))),
      freedom(vop.freedom), absbdom(vop.absbdom), idx(vop.idx) {}

  ostream & BornIWaveLOVOp::write(ostream & str) const {
    str<<"IWave-based Born LinOpValOp\n";
    str<<"  \n******* domain:\n\n";
    this->getDomain().write(str);
    str<<"  \n******* range:\n\n";
    this->getRange().write(str);
    str<<"  \n******* free param array\n\n";
    ps_printall(iwopfree->getPar(),str);
    str<<"  \n******* absb param array\n\n";
    ps_printall(iwopabsb->getPar(),str);
    return str;
  }

}
