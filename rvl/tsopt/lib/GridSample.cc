#include "GridSample.hh"

namespace RVL {
     
  void RVL::GridCopyOverlapFO::operator()(LocalDataContainer<float> & inside,
					  LocalDataContainer<float> const & outside) {
    try {
      
      GridDC & indc = dynamic_cast<GridDC &>(inside);
      GridDC const & outdc = dynamic_cast<GridDC const &>(outside);
      if (!(areCompatible(indc.getGrid(),outdc.getGrid()))) {
	RVLException e;
	e<<"GridCopyOverlapFO::operator()\n";
	e<<"  args not compatible\n";
	e<<"  src:\n";
	outdc.getGrid().RVL::Writeable::write(e);
	e<<"  tgt\n";
	indc.getGrid().RVL::Writeable::write(e);
	throw e;
      }
      //      cerr<<"1.1\n";
      int dim = indc.getGrid().getDimension();
      IPNT gsin; IPNT gein; IPNT gsout; IPNT geout; IPNT gs; IPNT ge;
      indc.getGrid().getAxisLimits(gsin,gein);
      outdc.getGrid().getAxisLimits(gsout,geout);
      for (int i=0;i<dim;i++) {
	gs[i]=max(gsin[i],gsout[i]); //cerr<<"gs["<<i<<"]="<<gs[i]<<"\n";
	ge[i]=min(gein[i],geout[i]); //cerr<<"ge["<<i<<"]="<<ge[i]<<"\n";
      }
      if (dim==1) {
	float * in = indc.getData1DGlobal();
	float const * out = outdc.getData1DGlobal();
	if (plus) for (int i0=gs[0]; i0<=ge[0]; i0++) in[i0]+=out[i0];
        else for (int i0=gs[0]; i0<=ge[0]; i0++) in[i0]=out[i0];
      }
      else if (dim==2) {
	//	cerr<<"dim=2\n";
	float ** in = indc.getData2DGlobal();
	//	cerr<<"in ptr\n";
	float const ** out = outdc.getData2DGlobal();
	//	cerr<<"out ptr\n";
	for (int i1=gs[1]; i1<=ge[1]; i1++) {
	  //	  cerr<<"i1="<<i1<<endl;
	  if (plus) for (int i0=gs[0]; i0<=ge[0]; i0++) in[i1][i0]+=out[i1][i0];
	  else for (int i0=gs[0]; i0<=ge[0]; i0++) in[i1][i0]=out[i1][i0];
	}
	//        cerr<<"end loop"<<endl;
      }
      else if (dim==3) {
	float *** in = indc.getData3DGlobal();
	float const *** out = outdc.getData3DGlobal();
	for (int i2=gs[2]; i2<=ge[2]; i2++) {
	  for (int i1=gs[1]; i1<=ge[1]; i1++) {
	    if (plus) 
	      for (int i0=gs[0]; i0<=ge[0]; i0++) {
		in[i2][i1][i0]+=out[i2][i1][i0];
	      }
	    else 
	      for (int i0=gs[0]; i0<=ge[0]; i0++) {
		in[i2][i1][i0]=out[i2][i1][i0];
	      }
	  }
	}
      }
      else {
	RVLException e;
	e<<"GridCopyOverlapFO::operator()\n";
	e<<"  dim="<<dim<<" > 3 not supported\n";
	e<<"  if this does not work for you - did we just hear you volunteer?\n";
	throw e;
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridCopyOverlapFO::operator()\n";
      e<<"  either first argument not TSDC or\n";
      e<<"  second is not GridDC:\n";
      e<<"  first arg:\n";
      inside.write(e);
      e<<"  second arg:\n";
      outside.write(e);
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridCopyOverlapFO::operator()\n";
      throw e;
    }
  }
  
  void GridExtendFO::operator()(LocalDataContainer<float> & x) {
    try {
      
      GridDC & gx = dynamic_cast<GridDC &>(x);
      if (!(areCompatible(gx.getGrid(),g))) {
	RVLException e;
	e<<"GridExtendFO::operator()\n";
	e<<"  args not compatible\n";
	e<<"  fixed grid:\n";
	g.RVL::Writeable::write(e);
	e<<"  extend grid\n";
	gx.getGrid().RVL::Writeable::write(e);
	throw e;
      }
	
      int dim = gx.getGrid().getDimension();
      IPNT gsin; IPNT gein; IPNT gsout; IPNT geout; IPNT gs; IPNT ge; IPNT idx;
      gx.getGrid().getAxisLimits(gsin,gein);
      g.getAxisLimits(gsout,geout);
      for (int i=0;i<dim;i++) {
	gs[i]=max(gsin[i],gsout[i]);
	ge[i]=min(gein[i],geout[i]);
      }
      
      if (dim == 3) {
	for (idx[2]=gsin[2]; idx[2]<gs[2]; idx[2]++) {
	  for (idx[1]=gs[1]; idx[1]<=ge[1]; idx[1]++) {
	    for (idx[0]=gs[0]; idx[0]<=ge[0]; idx[0]++) {
	      gx.getData3DGlobal()[idx[2]][idx[1]][idx[0]]
		= gx.getData3DGlobal()[gs[2]][idx[1]][idx[0]];
	    }
	  }
	}
	for (idx[2]=ge[2]+1; idx[2]<=gein[2]; idx[2]++) {
	  for (idx[1]=gs[1]; idx[1]<=ge[1]; idx[1]++) {
	    for (idx[0]=gs[0]; idx[0]<=ge[0]; idx[0]++) {
	      gx.getData3DGlobal()[idx[2]][idx[1]][idx[0]]
		= gx.getData3DGlobal()[ge[2]][idx[1]][idx[0]];
	    }
	  }
	}
	for (idx[2]=gsin[2]; idx[2]<=gein[2]; idx[2]++) {
	  for (idx[1]=gsin[1]; idx[1]<gs[1]; idx[1]++) {
	    for (idx[0]=gs[0]; idx[0]<=ge[0]; idx[0]++) {
	      gx.getData3DGlobal()[idx[2]][idx[1]][idx[0]]
		= gx.getData3DGlobal()[idx[2]][gs[1]][idx[0]];
	    }
	  }
	}
	for (idx[2]=gsin[2]; idx[2]<=gein[2]; idx[2]++) {
	  for (idx[1]=ge[1]+1; idx[1]<=gein[1]; idx[1]++) {
	    for (idx[0]=gs[0]; idx[0]<=ge[0]; idx[0]++) {
	      gx.getData3DGlobal()[idx[2]][idx[1]][idx[0]]
		= gx.getData3DGlobal()[idx[2]][ge[1]][idx[0]];
	    }
	  }
	}
	for (idx[2]=gsin[2]; idx[2]<=gein[2]; idx[2]++) {
	  for (idx[1]=gsin[1]; idx[1]<=gein[1]; idx[1]++) {
	    for (idx[0]=gsin[0]; idx[0]<gs[0]; idx[0]++) {
	      gx.getData3DGlobal()[idx[2]][idx[1]][idx[0]]
		= gx.getData3DGlobal()[idx[2]][idx[1]][gs[0]];
	    }
	  }
	}
	for (idx[2]=gsin[2]; idx[2]<=gein[2]; idx[2]++) {
	  for (idx[1]=gsin[1]; idx[1]<=gein[1]; idx[1]++) {
	    for (idx[0]=ge[0]+1; idx[0]<=gein[0]; idx[0]++) {
	      gx.getData3DGlobal()[idx[2]][idx[1]][idx[0]]
		= gx.getData3DGlobal()[idx[2]][idx[1]][ge[0]];
	    }
	  }
	}		
      }
      else if (dim == 2) {
	for (idx[1]=gsin[1]; idx[1]<gs[1]; idx[1]++) {
	  for (idx[0]=gs[0]; idx[0]<=ge[0]; idx[0]++) {
	    gx.getData2DGlobal()[idx[1]][idx[0]]
	      = gx.getData2DGlobal()[gs[1]][idx[0]];
	  }
	}
	for (idx[1]=ge[1]+1; idx[1]<=gein[1]; idx[1]++) {
	  for (idx[0]=gs[0]; idx[0]<=ge[0]; idx[0]++) {
	    gx.getData2DGlobal()[idx[1]][idx[0]]
	      = gx.getData2DGlobal()[ge[1]][idx[0]];
	  }
	}
	for (idx[1]=gsin[1]; idx[1]<=gein[1]; idx[1]++) {
	  for (idx[0]=gsin[0]; idx[0]<gs[0]; idx[0]++) {
	    gx.getData2DGlobal()[idx[1]][idx[0]]
	      = gx.getData2DGlobal()[idx[1]][gs[0]];
	  }
	}
	for (idx[1]=gsin[1]; idx[1]<=gein[1]; idx[1]++) {
	  for (idx[0]=ge[0]+1; idx[0]<=gein[0]; idx[0]++) {
	    gx.getData2DGlobal()[idx[1]][idx[0]]
	      = gx.getData2DGlobal()[idx[1]][ge[0]];
	  }
	}
      }
      else if (dim == 1) {
	for (idx[0]=gsin[0]; idx[0]<gs[0]; idx[0]++) {
	  gx.getData1DGlobal()[idx[0]]
	    = gx.getData1DGlobal()[gs[0]];
	}
	for (idx[0]=ge[0]+1; idx[0]<=gein[0]; idx[0]++) {
	  gx.getData1DGlobal()[idx[0]]
	    = gx.getData1DGlobal()[ge[0]];
	}
      }
    }	   
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridExtendFO::operator()\n";
      e<<"  argument not GridDC:\n";
      x.write(e);
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridExtendFO::operator()\n";
      throw e;
    }
  }

  // for testing - otherwise inefficient
  void GridtoTSOp::apply(Vector<float> const & x,
			 Vector<float> & y) const {
    try {
      if (abs(this->getTime() - tsamp) < tol) {
	y.zero();
	GridCopyOverlapFO f;
	// control or state
	Components<float> cy(y);
	// which component
	Components<float> ccy(cy[ind[0]]);
	ccy[ind[1]].eval(f,x);
	if (extend) {
	  GridExtendFO g(dom.getGrid());
	  ccy[ind[1]].eval(g);
	}
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridtoTSOp::apply\n";
      throw e;
    }
  }

  void GridtoTSOp::applyPlus(Vector<float> const & x,
			     Vector<float> & y) const {
    try {
      if (abs(this->getTime() - tsamp) < tol) {
	GridCopyOverlapFO f(true);
	// control or state
	Components<float> cy(y);
	// which component
	Components<float> ccy(cy[ind[0]]);
	ccy[ind[1]].eval(f,x);
	if (extend) {
	  GridExtendFO g(dom.getGrid());
	  ccy[ind[1]].eval(g);
	}
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridtoTSOp::applyPlus\n";
      throw e;
    }
  }

  void GridtoTSOp::applyAdj(Vector<float> const & x,
			    Vector<float> & y) const {
    try {
      if (abs(this->getTime() - tsamp) < tol) {
	GridCopyOverlapFO f;
	// control or state
	Components<float> cx(x);
	// which component
	Components<float> ccx(cx[ind[0]]);
	y.eval(f,ccx[ind[1]]);
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridtoTSOp::applyAdj\n";
      throw e;
    }
  }

  void GridtoTSOp::applyAdjPlus(Vector<float> const & x,
			    Vector<float> & y) const {
    try {
      if (abs(this->getTime() - tsamp) < tol) {
	GridCopyOverlapFO f(true);
	// control or state
	Components<float> cx(x);
	// which component
	Components<float> ccx(cx[ind[0]]);
	y.eval(f,ccx[ind[1]]);
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridtoTSOp::apply\n";
      throw e;
    }
  }

  GridtoTSOp::GridtoTSOp(float _tsamp,
			 std::vector<size_t> _ind,
			 GridSpace const & _gsp,
			 Space<float> const & _tsp,
			 bool _extend,
			 float _tol)
    : tsamp(_tsamp), ind(_ind), dom(_gsp), rng(_tsp), extend(_extend), tol(_tol) {

    try {
      this->setMinTime(tsamp);
      this->setMaxTime(tsamp);
      // rangs should be GridDomain
      GridDomain const & grng = dynamic_cast<GridDomain const &>(rng);
      if (ind.size()<2) {
	RVLException e;
	e<<"Error: GridtoTSOpt constructor\n";
	e<<"  index array must have size at least 2 to indicate indices\n";
	e<<"  in control or state spaces\n";
	throw e;
      }
      if ((ind[0]!=0) && (ind[0]!=1)) {
	RVLException e;
	e<<"Error: GridtoTSOpt constructor\n";
	e<<"  ind[0] must = 0 (control) or = 1 (state)\n";
	throw e;
      }
      if (ind[0]==0) {
	TSSpace<float> const & ctrl = dynamic_cast<TSSpace<float> const & >(grng[0]);
	if (ind[1]>ctrl.getSize()) {
	  RVLException e;
	  e<<"Error: GridtoTSOpt constructor\n";
	  e<<"  ind[1]="<<ind[1]<<" but must be less than number of \n";
	  e<<"  control components = "<<ctrl.getSize()<<"\n";
	  throw e;
	}
      }
      else if (ind[0]==1) {
	TSSpace<float> const & state = dynamic_cast<TSSpace<float> const & >(grng[1]);	
	if (ind[1]>state.getSize()) {
	  RVLException e;
	  e<<"Error: GridtoTSOpt constructor\n";
	  e<<"  ind[1]="<<ind[1]<<" but must be less than number of \n";
	  e<<"  state components = "<<state.getSize()<<"\n";
	  throw e;
	}
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridtoTSOp constructor\n";
      e<<"  range space not GridDomain\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridtoTSOp constructor\n";
      throw e;
    }
  }

  void GridtoTSOp::applyPlusOp(Vector<float> const & x,
			       Vector<float> & y) const {
    try {
      sanity(x,y);
      this->applyPlus(x,y);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridtoTSOp::applyPlusOp\n";
      throw e;
    }
  }

  void GridtoTSOp::applyPlusAdjOp(Vector<float> const & x,
				  Vector<float> & y) const {
    try {
      sanity(y,x);
      this->applyAdjPlus(x,y);
    }
    catch (RVLException & e) {
      e<<"\ncalled from LinearOp::applyOp (with lincomb) \n";
      throw e;
    }
  }
    
}
