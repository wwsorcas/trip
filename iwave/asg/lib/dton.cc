#include "dton.hh"
//#define IWAVE_VERBOSE
//#define RECORD
namespace TSOpt {

  void DToNIWaveLOVOp::apply0(const Vector<float> & x0,
			      const Vector<float> & x1, 
			      Vector<float> & y) const {

    try {
      RVL::Vector<float> inp(iwop->getLinDomain());
      Components<float> cinp(inp);
      cinp[0].copy(x1);
      cinp[1].zero();
      RVL::Vector<float> outp(iwop->getRange());
      Components<float> coutp(outp);
      export_apply0(*iwop,x0,inp,outp);
      if (this->appinv) {
	y.copy(coutp[1]);
      }
      else {
	y.copy(coutp[0]);
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: DToNIWaveLOVOp::apply0\n";
      e<<"  failure to cast one of the domain spaces as a grid space\n";
      throw e;
    }    
    catch (RVLException & e) {
      e<<"\ncalled from DToNIWaveLOVOp::apply0\n";
      throw e;
    }

  }

  void DToNIWaveLOVOp::applyAdj0(const Vector<float> & x0,
				 const Vector<float> & y, 
				 Vector<float> & x1) const {

    try {
    }
  
    catch (bad_cast) {
      RVLException e;
      e<<"Error: DToNIWaveLOVOp::applyAdj0\n";
      e<<"  failure to cast one of the domain spaces as a grid space\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from DToNIWaveLOVOp::applyAdj0\n";
      throw e;
    }
  }
      
  void DToNIWaveLOVOp::applyPartialDeriv0(const Vector<float> & x0,
					  const Vector<float> & x1,
					  const Vector<float> & dx0,
					  Vector<float> & dy) const {} 
    
  void DToNIWaveLOVOp::applyAdjPartialDeriv0(const Vector<float> & x0,
					     const Vector<float> & x1,
					     const Vector<float> & dy,
					     Vector<float> & dx0) const {}
    
  DToNIWaveLOVOp::DToNIWaveLOVOp(PARARRAY pars,
				 FILE * stream)
    : dom(2) {

    // appinv vs. Euclidean IP in defns of operators
    appinv =  valparse<int>(pars,"appinv");

    std::string tmp;

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
    if ((pkeys.size()<2) ||
	(pkeys[0] != "bulkmod") ||
	(pkeys[1] != "buoyancy")) {
      RVLException e;
      e<<"Error: dton\n";
      e<<"  IWaveInfo should define 2 nonlinear input parameters\n";
      e<<"  key[0]=bulkmod key[1]=buoyancy\n";
      throw e;
    }

    // find data (output active linear) keywords present in par file
    std::vector<std::string> dkeys;
    for (int i=0; i< ic.get_num_iokeys(); i++) {
      if ((ic.iwave_iokeys[i].input == 0) &&
	  (ic.iwave_iokeys[i].active == 2) ) {
	tmp = RVL::valparse<std::string>(pars,ic.iwave_iokeys[i].keyword,"");
	if (tmp.size()>0) {
	  dkeys.push_back(ic.iwave_iokeys[i].keyword);
	}
      }
    }
    if ((dkeys.size() < 2) ||
	(dkeys[0] != "data_p") ||
	(dkeys[1] != "data_v0")) {
      RVLException e;
      e<<"Error: dton\n";
      e<<"  IWaveInfo should define >= 2 linear output parameters\n";
      e<<"  key[0]=data_p key[1]=data_v0 \n";
      throw e;
    }
      
    // find source (input active linear) keywords present in par file
    std::vector<std::string> skeys;
    for (int i=0; i< ic.get_num_iokeys(); i++) {
      if ((ic.iwave_iokeys[i].input == 1) &&
	  (ic.iwave_iokeys[i].active == 2) ) {
	tmp = RVL::valparse<std::string>(pars,ic.iwave_iokeys[i].keyword,"");
	if (tmp.size()>0) {
	  skeys.push_back(ic.iwave_iokeys[i].keyword);
	}
      }
    }
    if ((skeys.size()<2) ||
	(skeys[0] != "source_p") ||
	(skeys[1] != "source_v0")) {
      RVLException e;
      e<<"Error: dton\n";
      e<<"  IWaveInfo should define >= 2 linear output parameters\n";
      e<<"  key[0]=source_p key[1]=source_v0 \n";
      throw e;
    }

    // modify par array to create domain on chosen keys
    PARARRAY * iwpars = ps_new();
    ps_copy(&iwpars, pars);

    // pml parameters
    float pmlwidth = RVL::valparse<float>(pars,"pmlwidth");
    ps_slfloat(*iwpars,"nl1",pmlwidth);      
    ps_slfloat(*iwpars,"nr1",pmlwidth);            
    ps_slfloat(*iwpars,"nl2",pmlwidth);      
    ps_slfloat(*iwpars,"nr2",pmlwidth);            
    ps_slfloat(*iwpars,"nl3",pmlwidth);      
    ps_slfloat(*iwpars,"nr3",pmlwidth);

    // check that neither source_v0 nor data_v0 are present in pars
    if ((RVL::parse(pars,"source_v0",tmp)) ||
	(RVL::parse(pars,"data_v0",tmp))) {
      RVLException e;
      e<<"Error: dton\n";
      e<<"  input par array should define neither source_v0 nor data_v0\n";
      throw e;
    }

    // add copies of source_p and data_p to iwpars, to initialize
    // source_vz and data_vz spaces
    if (!(RVL::parse(pars,"source_p",tmp))) {
            RVLException e;
      e<<"Error: dton\n";
      e<<"  input par array must define source_p\n";
      throw e;
    }
    ps_slcstring(*iwpars,"source_v0",tmp.c_str());

    if (!(RVL::parse(pars,"data_p",tmp))) {
            RVLException e;
      e<<"Error: dton\n";
      e<<"  input par array must define data_p\n";
      throw e;
    }
    ps_slcstring(*iwpars,"data_v0",tmp.c_str());
    
    // basic op which includes velocity source and data
    iwop = make_shared<TSOpt::IWaveLOVOp>(*iwpars,stream);

    dom.set(iwop->getProductDomain()[0],0);

    RVL::ProductSpace<float> const & iwldom =
      dynamic_cast<RVL::ProductSpace<float> const &>(iwop->getProductDomain()[1]);

    RVL::StdProductSpace<float> ldom(iwldom[0]);

    dom.set(ldom,1);

    ps_delete(&iwpars);
  }

  DToNIWaveLOVOp::DToNIWaveLOVOp(DToNIWaveLOVOp const & op)
    : dom(op.dom), iwop(op.iwop), appinv(op.appinv) {}


  ostream & DToNIWaveLOVOp::write(ostream & str) const {
    str<<"IWave-based DToN LinOpValOp\n";
    str<<"  \n******* domain:\n\n";
    this->getDomain().write(str);
    str<<"  \n******* range:\n\n";
    this->getRange().write(str);
    str<<"  \n******* iwop param array\n\n";
    ps_printall(iwop->getPar(),str);
    return str;
  }

}
