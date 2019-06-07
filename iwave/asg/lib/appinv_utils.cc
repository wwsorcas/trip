#include "appinv_utils.hh"

namespace TSOpt {

  void ASGModelWeightOp::apply(RVL::Vector<float> const & x,
			       RVL::Vector<float> & y) const {
    try {

      Components<float> cx(x);
      Components<float> cy(y);
      
      if (symm) {
	
	float power   = -1.0f;
	float bulkpwr = -1.5f;
	float buoypwr = -0.5f;
	float fac     = 1.0f/32.0f;
	
	GridScaleFO bulksc(1.0f,bulkpwr);
	GridScaleFO buoysc(1.0f,buoypwr);
	GridZFTScaleFO ftscale(dz,power,1,band,0);
	RVL::MPISerialFunctionObject<float> mpibulksc(bulksc);
	RVL::MPISerialFunctionObject<float> mpibuoysc(buoysc);
	RVL::MPISerialFunctionObject<float> mpiftscale(ftscale);
	
	cy[0].copy(cx[0]);
	cy[0].eval(mpibuoysc,buoy);
	cy[0].eval(mpibulksc,bulk);
	cy[0].eval(mpiftscale);
	cy[0].eval(mpibulksc,bulk);	
	cy[0].eval(mpibuoysc,buoy);
	cy[0].scale(fac);
      }
      else {
	RVL::RVLException e;
	e<<"Error: ASGModelWeightOp::apply\n";
	e<<"  non available in non-symm case\n";
      }      
    }
    catch (RVL::RVLException & e) {
      e<<"\ncalled from ASGModelWeightOp::apply\n";
      throw e;
    }
  }

  void ASGModelWeightOp::applyAdj(RVL::Vector<float> const & x,
				  RVL::Vector<float> & y) const {
    try {
      if (symm) this->apply(x,y);
      else {
	RVL::RVLException e;
	e<<"Error: ASGModelWeightOp::applyAdj\n";
	e<<"  not available in non-symm case\n";
      }	
    }
    catch (RVL::RVLException & e) {
      e<<"\ncalled from ASGModelWeightOp::applyAdj\n";
      throw e;
    }
  }      

  void ASGModelWeightOp::applyInv(RVL::Vector<float> const & x,
				  RVL::Vector<float> & y) const {
    try {
      float power   = 1.0f;
      float bulkpwr = 1.5f;
      float buoypwr = 0.5f;
      float fac     = 32.0f;

      Components<float> cx(x);
      Components<float> cy(y);
      
      if (symm) {
	
	GridScaleFO bulksc(1.0f,bulkpwr);
	GridScaleFO buoysc(1.0f,buoypwr);
	GridZFTScaleFO ftscale(dz,power,1,band,0);
	RVL::MPISerialFunctionObject<float> mpibulksc(bulksc);
	RVL::MPISerialFunctionObject<float> mpibuoysc(buoysc);
	RVL::MPISerialFunctionObject<float> mpiftscale(ftscale);
	
	cy[0].copy(cx[0]);
	cy[0].eval(mpibuoysc,buoy);
	cy[0].eval(mpibulksc,bulk);
	cy[0].eval(mpiftscale);
	cy[0].eval(mpibulksc,bulk);	
	cy[0].eval(mpibuoysc,buoy);
	cy[0].scale(fac);
      }
      else {
	GridScaleFO bulksc(1.0f,2.0f*bulkpwr);
	GridScaleFO buoysc(1.0f,2.0f*buoypwr);
	GridFwdDerivFO fwdderiv(0,1.0f/dz);
	RVL::MPISerialFunctionObject<float> mpibulksc(bulksc);
	RVL::MPISerialFunctionObject<float> mpibuoysc(buoysc);
	RVL::MPISerialFunctionObject<float> mpifwdderiv(fwdderiv);
	
	cy[0].eval(mpifwdderiv,cx[0]);
	cy[0].eval(mpibuoysc,buoy);
	cy[0].eval(mpibulksc,bulk);
	cy[0].scale(fac);
      }
      
    }
    catch (RVL::RVLException & e) {
      e<<"\ncalled from ASGModelWeightOp::apply\n";
      throw e;
    }
  }
  
  void ASGModelWeightOp::applyInvAdj(RVL::Vector<float> const & x,
				     RVL::Vector<float> & y) const {
    try {
      if (symm) this->applyInv(x,y);
      else {
	RVL::RVLException e;
	e<<"Error: ASGModelWeightOp::applyInvAdj\n";
	e<<"  non available in non-symm case\n";
      }
    }
    catch (RVL::RVLException & e) {
      e<<"\ncalled from ASGModelWeightOp::applyInvAdj\n";
      throw e;
    }
  }      
  
  ASGModelWeightOp::ASGModelWeightOp(RVL::Vector<float> const & _bulk,
				     RVL::Vector<float> const & _buoy,
				     RVL::Vector<float> const & _extbulk,
				     float locut, float lopas, float hipas, float hicut,
				     int _symm)
    : bulk(_bulk), buoy(_buoy), extbulk(_extbulk), band(new float[4]), dz(1.0f), symm(_symm) {
    try {
      if (bulk.getSpace() != buoy.getSpace()) {
	RVL::RVLException e;
	e<<"Error: ASGModelWeightOp constructor\n";
	e<<"  bulk, buoy grids not same\n";
	e<<"  bulk space:\n";
	bulk.getSpace().write(e);
	e<<"  buoy space:\n";
	buoy.getSpace().write(e);
	throw e;
      }
      band[0]=locut;
      band[1]=lopas;
      band[2]=hipas;
      band[3]=hicut;
      
      // find dz, allowing for nontrivial axis assignment
      // note that dz stays = 1.0 on all ranks other than 0
      myGridSpace const & gsp = dynamic_cast<myGridSpace const &>(bulk.getSpace());
      if (retrieveGlobalRank()==0) {
	RPNT d;
	IPNT id;
	get_d(d,gsp.getGrid());
	get_id(id,gsp.getGrid());
	int idx=-1;
	for (int i=0;i<gsp.getGrid().dim;i++) 
	  if (id[i]==0) idx=i;
	// default choice is axis not assigned - assume idx=0
	if (idx<0) dz=d[0];
	else dz=d[idx];
      }
    }
    catch (bad_cast) {
      RVL::RVLException e;
      e<<"Error: ASGModelWeightOp constructor\n";
      e<<"  bulk mod space not grid space\n";
      throw e;
    }
    catch (RVL::RVLException & e) {
      e<<"\ncalled from ASGModelWeightOp constructor\n";
      throw e;
    }
  }
  ASGModelWeightOp::ASGModelWeightOp(ASGModelWeightOp const & op)
    : bulk(op.bulk), buoy(op.buoy), extbulk(op.extbulk), symm(op.symm), band(new float[4]) {
    band[0]=op.band[0];
    band[1]=op.band[1];
    band[2]=op.band[2];
    band[3]=op.band[3];
    dz=op.dz;
  }

  ostream & ASGModelWeightOp::write(ostream & str) const {
    str<<"ASGModelWeightOp\n";
    str<<"  bulk modulus:\n";
    bulk.write(str);
    str<<"  buoyancy:\n";
    buoy.write(str);
    str<<"  extd bulk modulus:\n";
    extbulk.write(str);
    str<<"locut="<<band[0]<<"\n";
    str<<"lopas="<<band[1]<<"\n";
    str<<"hipas="<<band[2]<<"\n";
    str<<"hicut="<<band[3]<<"\n";
    str<<"symmetry flag = "<<symm<<"\n";
    return str;
  }
}


