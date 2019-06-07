#include "gridops.hh"

namespace TSOpt {

  using RVL::ScalarFieldTraits;
  using RVL::SpaceTest;
  using RVL::Operator;
  using RVL::LinearOp;
  using RVL::Space;
  using RVL::ProductSpace;
  using RVL::Vector;
  using RVL::Components;
  using RVL::ProtectedDivision;
  using RVL::RnArray;
  using RVL::RVLScale;
  using RVL::BinaryLocalFunctionObject;
  using RVL::RVLException;
  using RVL::ContentPackage;
  using RVL::LocalDataContainer;
  using RVL::MPISerialFunctionObject;
  using RVL::MPISerialFunctionObjectRedn;

  void GridMaskFO::operator()(LocalDataContainer<float> & x,
			      LocalDataContainer<float> const & y) {
    try {
      // cerr<<"GridWindowFO::operator() begin\n";
      ContentPackage< float, RARR > const & gy =
	dynamic_cast<ContentPackage< float, RARR > const &>(y);
            
      ContentPackage< float, RARR > & gx =
	dynamic_cast<ContentPackage< float, RARR > &>(x);
            
      // precondition - metadata are same dimn
      RARR & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();
      int dimx; int dimy;
      ra_ndim(&rax,&dimx);
      ra_ndim(&ray,&dimy);
      if (dimx != dimy) {
	RVLException e;
	e<<"Error: GridMaskFO::operator()\n";
	e<<"arguments have different dims:\n";
	e<<"dimx="<<dimx<<" dimy="<<dimy<<"\n";
	throw e;
      }
      // compute grid params
      IPNT gsx; IPNT gex;
      IPNT gsy; IPNT gey;
      IPNT s; IPNT e;
      ra_a_gse(&rax,gsx,gex);
      ra_a_gse(&ray,gsy,gey);
      // calculate grid overlap
      for (int ii=0;ii<dimx;ii++)  {
	s[ii]=max(gsy[ii],gsx[ii]);
	e[ii]=min(gey[ii],gex[ii]);
      }
            
      IPNT i;
      RPNT fac;
      RASN(fac,RPNT_1);
#if RARR_MAX_NDIM > 0
      if (dimx==1) {
#pragma ivdep
	for (i[0]=s[0]+siw[0];i[0]<=e[0]-eiw[0];i[0]++) {
          if (i[0] < width[0]+s[0]+siw[0]) 
            fac[0] = cosfun((width[0]+s[0]+siw[0]-i[0])/width[0]);
          else if (i[0] > e[0]-eiw[0]-width[0]) 
	    fac[0] = cosfun((i[0]-e[0]+eiw[0]+width[0])/width[0]);
	  else fac[0] = 1.0f; //iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(i[0]-s[0]-siw[0]+1))/float(width[0]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(e[0]-eiw[0]+1-i[0]))/float(width[0]))));
	  if (bias) {
	    rax._s1[i[0]]+=ray._s1[i[0]]*fac[0];
	  }
	  else{
	    rax._s1[i[0]]=ray._s1[i[0]]*fac[0];
	  }
	}
      }
#endif
#if RARR_MAX_NDIM > 1
      if (dimx==2) {
	for (i[1]=s[1]+siw[1];i[1]<=e[1]-eiw[1];i[1]++) {
          if (i[1] < width[1]+s[1]+siw[1]) 
            fac[1] = cosfun((width[1]+s[1]+siw[1]-i[1])/width[1]);
          else if (i[1] > e[1]-eiw[1]-width[1]) 
	    fac[1] = cosfun((i[1]-e[1]+eiw[1]+width[1])/width[1]);
	  else fac[1] = 1.0f; // iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(i[1]-s[1]-siw[1]+1))/float(width[1]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(e[1]-eiw[1]+1-i[1]))/float(width[1]))));
#pragma ivdep
	  for (i[0]=s[0]+siw[0];i[0]<=e[0]-eiw[0];i[0]++) {
	    if (i[0] < width[0]+s[0]+siw[0]) 
	      fac[0] = cosfun((width[0]+s[0]+siw[0]-i[0])/width[0]);
	    else if (i[0] > e[0]-eiw[0]-width[0]) 
	      fac[0] = cosfun((i[0]-e[0]+eiw[0]+width[0])/width[0]);
	    else fac[0] = 1.0f; //iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(i[0]-s[0]-siw[0]+1))/float(width[0]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(e[0]-eiw[0]+1-i[0]))/float(width[0]))));
	    if (bias) {
	      rax._s2[i[1]][i[0]]+=fac[0]*fac[1]*ray._s2[i[1]][i[0]];
	    }
	    else{
	      rax._s2[i[1]][i[0]]=fac[0]*fac[1]*ray._s2[i[1]][i[0]];
	    }
	  }
	}
      }
#endif
#if RARR_MAX_NDIM > 2
      if (dimx==3) {
	for (i[2]=s[2]+siw[2];i[2]<=e[2]-eiw[2];i[2]++) {
          if (i[2] < width[2]+s[2]+siw[2]) 
            fac[2] = cosfun((width[2]+s[2]+siw[2]-i[2])/width[2]);
          else if (i[2] > e[2]-eiw[2]-width[2]) 
	    fac[2] = cosfun((i[2]-e[2]+eiw[2]+width[2])/width[2]);
	  else fac[2] = 1.0f; // iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(i[2]-s[2]-siw[2]+1))/float(width[2]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(e[2]-eiw[2]+1-i[2]))/float(width[2]))));
	  for (i[1]=s[1]+siw[1];i[1]<=e[1]-eiw[1];i[1]++) {
	    if (i[1] < width[1]+s[1]+siw[1]) 
	      fac[1] = cosfun((width[1]+s[1]+siw[1]-i[1])/width[1]);
	    else if (i[1] > e[1]-eiw[1]-width[1]) 
	      fac[1] = cosfun((i[1]-e[1]+eiw[1]+width[1])/width[1]);
	    else fac[1] = 1.0f; // iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(i[1]-s[1]-siw[1]+1))/float(width[1]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(e[1]-eiw[1]+1-i[1]))/float(width[1]))));
#pragma ivdep
	    for (i[0]=s[0]+siw[0];i[0]<=e[0]-eiw[0];i[0]++) {
	      if (i[0] < width[0]+s[0]+siw[0]) 
		fac[0] = cosfun((width[0]+s[0]+siw[0]-i[0])/width[0]);
	      else if (i[0] > e[0]-eiw[0]-width[0]) 
		fac[0] = cosfun((i[0]-e[0]+eiw[0]+width[0])/width[0]);
	      else fac[0] = 1.0f; //iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(i[0]-s[0]-siw[0]+1))/float(width[0]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(e[0]-eiw[0]+1-i[0]))/float(width[0]))));
	      if (bias) {
		rax._s3[i[2]][i[1]][i[0]]+=fac[0]*fac[1]*ray._s3[i[2]][i[1]][i[0]];
	      }
	      else{
		rax._s3[i[2]][i[1]][i[0]]=fac[0]*fac[1]*ray._s3[i[2]][i[1]][i[0]];
	      }
	    }
	  }
	}
      }
#endif
      if (dimx<1 || dimx>3) {
	RVLException e;
	e<<"Error: GridMaskFO::operator()\n";
	e<<"dim = "<<dimx<<" outside of admissible set {1, 2, 3}\n";
	throw e;
      }
      // cerr<<"GridWindowFO::operator() end\n";
    }
    catch (bad_cast) {
      RVLException e;
      e<<"\nError: GridMaskFO::operator()\n";
      e<<"at least one arg is not ContentPackage<float,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridMaskFO::operator()\n";
      throw e;
    }
  }
    
   
#ifdef IWAVE_USE_MPI
  typedef MPIGridSpace myGridSpace;
#else
  typedef GridSpace myGridSpace;
#endif
    
  GridMaskOp::GridMaskOp(Space<float> const & _dom,
			 Vector<float> const & _bg,
			 RPNT const & sw, RPNT const & ew, RPNT const & _width)
    : dom(_dom), bg(_bg) {
    try {
      // generic initialization of iw
      IASN(siw,IPNT_0);
      IASN(eiw,IPNT_0);
      RASN(width,_width);
  
      // branch on product structure - unfortunately required
      ProductSpace<float> const * pdom = dynamic_cast<ProductSpace<float> const *>(&dom);
      ProductSpace<float> const * prng = dynamic_cast<ProductSpace<float> const *>(&(bg.getSpace()));
      if (pdom && prng) {
	if (pdom->getSize() != prng->getSize()) {
	  RVLException e;
	  e<<"Error GridMaskOp constructor\n";
	  e<<"  domain and range are product spaces with different numbers of components\n";
	  throw e;
	}
	// check that component gridspaces are pairwise compatible and compatible
	// with 0th domain component - then they are all compatible
                
	myGridSpace const & gref = dynamic_cast<myGridSpace const &>((*pdom)[0]);
	for (int j=0; j<(int)pdom->getSize(); j++) {
	  myGridSpace const & gdom = dynamic_cast<myGridSpace const &>((*pdom)[j]);
	  myGridSpace const & grng = dynamic_cast<myGridSpace const &>((*prng)[j]);
                    
	  if (retrieveGlobalRank()==0) {
	    if (compatible_grid(gdom.getGrid(),grng.getGrid()) ||
		compatible_grid(gref.getGrid(),grng.getGrid())) {
	      RVLException e;
	      e<<"Error: GridMaskOp constructor\n";
	      e<<"  domain, range defined on incompatible grids\n";
	      e<<"  product case, component = "<<j<<"\n";
	      e<<"  domain:\n";
	      for (int i=0;i<gdom.getGrid().gdim;i++)
		e<<"    axis "<<i<<" d="<<gdom.getGrid().axes[i].d<<" o="<<gdom.getGrid().axes[i].o<<"\n";
	      e<<"  range:\n";
	      for (int i=0;i<grng.getGrid().gdim;i++)
		e<<"    axis "<<i<<" d="<<grng.getGrid().axes[i].d<<" o="<<grng.getGrid().axes[i].o<<"\n";
	      throw e;
	    }
	  }
	}
	if (retrieveGlobalRank()==0) {
	  grid const & g = gref.getGrid();
	  for (int i=0; i< g.dim; i++) {
	    siw[i]=(int) (sw[i]/(g.axes[i].d) + 0.1);
	    eiw[i]=(int) (ew[i]/(g.axes[i].d) + 0.1);
            
	    siw[i]=iwave_max(siw[i],1);
	    eiw[i]=iwave_max(eiw[i],1);
            // cerr << "g.axes[" << i << "].d=" << g.axes[i].d << endl;
	  }
	}
      }
      else {
	myGridSpace const & gdom = dynamic_cast<myGridSpace const &> (dom);
	myGridSpace const & grng = dynamic_cast<myGridSpace const &>(bg.getSpace());
	if (retrieveGlobalRank()==0) {
	  if (compatible_grid(gdom.getGrid(),grng.getGrid())) {
	    RVLException e;
	    e<<"Error: GridMaskOp constructor\n";
	    e<<"  domain, range defined on incompatible grids\n";
	    e<<"  domain:\n";
	    for (int i=0;i<gdom.getGrid().gdim;i++) 
	      e<<"    axis "<<i<<" d="<<gdom.getGrid().axes[i].d<<" o="<<gdom.getGrid().axes[i].o<<"\n";
	    e<<"  range:\n";
	    for (int i=0;i<grng.getGrid().gdim;i++) 
	      e<<"    axis "<<i<<" d="<<grng.getGrid().axes[i].d<<" o="<<grng.getGrid().axes[i].o<<"\n";
	    throw e;
	  }
	  grid const & g = gdom.getGrid();
	  for (int i=0; i< g.dim; i++) {
	    siw[i]=(int) (sw[i]/(g.axes[i].d) + 0.1);
	    eiw[i]=(int) (ew[i]/(g.axes[i].d) + 0.1);

	    siw[i]=iwave_max(siw[i],0);
	    eiw[i]=iwave_max(eiw[i],0);
	  }
	}
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridMaskOp constructor\n";
      e<<"  either domain or range is neither product nor a GridSpace,\n";
      e<<"  or some component is not a GridSpace\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridMaskOp constructor\n";
      throw e;
    }
  }
    
  void GridMaskOp::apply(Vector<float> const & x,
			 Vector<float> & y) const {
    try {
      GridMaskFO op(siw,eiw,width,true);
      MPISerialFunctionObject<float> mpiop(op);
      y.copy(bg);
      y.eval(mpiop,x);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridMaskOp::apply\n";
      throw e;
    }
  }
    
  void GridMaskOp::applyDeriv(Vector<float> const & x,
			      Vector<float> const & dx,
			      Vector<float> & dy) const {
    try {
      GridMaskFO op(siw,eiw,width,false);
      MPISerialFunctionObject<float> mpiop(op);
      dy.zero();
      dy.eval(mpiop,dx);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridMaskOp::applyDeriv\n";
      throw e;
    }
  }
    
  void GridMaskOp::applyAdjDeriv(Vector<float> const & x,
				 Vector<float> const & dy,
				 Vector<float> & dx) const {
    try {
      GridMaskFO op(siw,eiw,width,false);
      MPISerialFunctionObject<float> mpiop(op);
      dx.zero();
      dx.eval(mpiop,dy);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridMaskOp::applyAdjDeriv\n";
      throw e;
    }
        
  }
  ostream & GridMaskOp::write(ostream & str) const {
    if (!retrieveGlobalRank()) {
      str<<"GridMaskOp\n";
    }
    return str;
  }
    
  void GridWindowFO::operator()(LocalDataContainer<float> & x,
				LocalDataContainer<float> const & y) {
    try {
      // cerr<<"GridWindowFO::operator() begin\n";
      ContentPackage< float, RARR > const & gy = 
	dynamic_cast<ContentPackage< float, RARR > const &>(y);

      ContentPackage< float, RARR > & gx = 
	dynamic_cast<ContentPackage< float, RARR > &>(x);
      
      // precondition - metadata are same dimn
      RARR & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();
      int dimx; int dimy;
      ra_ndim(&rax,&dimx);
      ra_ndim(&ray,&dimy);
      if (dimx != dimy) {
	RVLException e;
	e<<"Error: GridWindow::operator()\n";
	e<<"arguments have different dims:\n";
	e<<"dimx="<<dimx<<" dimy="<<dimy<<"\n";
	throw e;
      }

      // compute grid params
      IPNT gsx; IPNT gex;
      IPNT gsy; IPNT gey;
      IPNT s; IPNT e; 
      ra_a_gse(&rax,gsx,gex);
      ra_a_gse(&ray,gsy,gey);
      // calculate grid overlap
      for (int ii=0;ii<dimx;ii++)  {
	s[ii]=max(gsy[ii],gsx[ii]);
	e[ii]=min(gey[ii],gex[ii]);
      }	

      IPNT i;
      RPNT fac;
      RASN(fac,RPNT_1);
#if RARR_MAX_NDIM > 0
      if (dimx==1) {
#pragma ivdep
	for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	  fac[0] = iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(i[0]-s[0]+1))/float(iw[0]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(e[0]+1-i[0]))/float(iw[0]))));
	  if (bias) {
	    rax._s1[i[0]]+=fac[0]*ray._s1[i[0]];
	  }
	  else {
	    rax._s1[i[0]] =fac[0]*ray._s1[i[0]];
	  }
	}
      }
#endif
#if RARR_MAX_NDIM > 1 
      if (dimx==2) {
	for (i[1]=s[1];i[1]<=e[1];i[1]++) {
	  fac[1] = iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(i[1]-s[1]+1))/float(iw[1]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(e[1]+1-i[1]))/float(iw[1]))));
#pragma ivdep
	  for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	    fac[0] = fac[1]*iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(i[0]-s[0]+1))/float(iw[0]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(e[0]+1-i[0]))/float(iw[0]))));
	    if (bias) {
	      rax._s2[i[1]][i[0]]+=fac[0]*ray._s2[i[1]][i[0]];	  
	    }
	    else {
	      rax._s2[i[1]][i[0]] =fac[0]*ray._s2[i[1]][i[0]];	  
	    }		
	  }
	}
      }
#endif
#if RARR_MAX_NDIM > 2 
      if (dimx==3) {
	for (i[2]=s[2];i[2]<=e[2];i[2]++) {
	  fac[2] = iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(i[2]-s[2]+1))/float(iw[2]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(e[2]+1-i[2]))/float(iw[2]))));
	  for (i[1]=s[1];i[1]<=e[1];i[1]++) {
	    fac[1] = fac[2]*iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(i[1]-s[1]+1))/float(iw[1]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(e[1]+1-i[1]))/float(iw[1]))));
#pragma ivdep
	    for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	      fac[0] = fac[1]*iwave_min(iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(i[0]-s[0]+1))/float(iw[0]))),iwave_min(REAL_ONE,iwave_max(REAL_ZERO,(float(e[0]+1-i[0]))/float(iw[0]))));
	      if (bias) {
		rax._s3[i[2]][i[1]][i[0]]+=fac[0]*ray._s3[i[2]][i[1]][i[0]];	  
	      }
	      else {
		rax._s3[i[2]][i[1]][i[0]] =fac[0]*ray._s3[i[2]][i[1]][i[0]];	  
	      }
	    }
	  }
	}
      }
#endif
      if (dimx<1 || dimx>3) {
	RVLException e;
	e<<"Error: GridWindowFO::operator()\n";
	e<<"dim = "<<dimx<<" outside of admissible set {1, 2, 3}\n";
	throw e;
      }
      // cerr<<"GridWindowFO::operator() end\n";
    }
    catch (bad_cast) {
      RVLException e;
      e<<"\nError: GridWindowFO::operator()\n";
      e<<"at least one arg is not ContentPackage<float,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridWindowFO::operator()\n";
      throw e;
    }
  }



  void GridWindowOp::initialize(RPNT const w) {

    try {

      // generic initialization of iw
      IASN(iw,IPNT_0);

      // branch on product structure - unfortunately required
      ProductSpace<float> const * pdom = dynamic_cast<ProductSpace<float> const *>(dom.get());
      ProductSpace<float> const * prng = dynamic_cast<ProductSpace<float> const *>(&(bg->getSpace()));
      if (pdom && prng) {
	if (pdom->getSize() != prng->getSize()) {
	  RVLException e;
	  e<<"Error GridWindowOp constructor\n";
	  e<<"  domain and range are product spaces with different numbers of components\n";
	  throw e;
	}
	// check that component gridspaces are pairwise compatible and compatible
	// with 0th domain component - then they are all compatible

	myGridSpace const & gref = dynamic_cast<myGridSpace const &>((*pdom)[0]);
	for (int j=0; j<(int)pdom->getSize(); j++) {
	  myGridSpace const & gdom = dynamic_cast<myGridSpace const &>((*pdom)[j]);
	  myGridSpace const & grng = dynamic_cast<myGridSpace const &>((*prng)[j]);
	  
	  if (retrieveGlobalRank()==0) {
	    if (compatible_grid(gdom.getGrid(),grng.getGrid()) ||
		compatible_grid(gref.getGrid(),grng.getGrid())) {
	      RVLException e;
	      e<<"Error: GridWindowOp constructor\n";
	      e<<"  domain, range defined on incompatible grids\n";
	      e<<"  product case, component = "<<j<<"\n";
	      e<<"  domain:\n";
	      for (int i=0;i<gdom.getGrid().gdim;i++) 
		e<<"    axis "<<i<<" d="<<gdom.getGrid().axes[i].d<<" o="<<gdom.getGrid().axes[i].o<<"\n";
	      e<<"  range:\n";
	      for (int i=0;i<grng.getGrid().gdim;i++) 
		e<<"    axis "<<i<<" d="<<grng.getGrid().axes[i].d<<" o="<<grng.getGrid().axes[i].o<<"\n";
	      throw e;
	    }
	  }
	}
	if (retrieveGlobalRank()==0) {
	  grid const & g = gref.getGrid();
	  for (int i=0; i< g.dim; i++) {
	    iw[i]=(int) (w[i]/(g.axes[i].d) + 0.1);
	    iw[i]=iwave_max(iw[i],1);
	  }
	}
      }
      else {
	myGridSpace const * gdom = dynamic_cast<myGridSpace const *> (dom.get());
	myGridSpace const * grng = dynamic_cast<myGridSpace const *>(&(bg->getSpace()));
	if (!gdom) {
	  RVLException e;
	  e<<"Error: GridWindowOp constructor\n";
	  e<<"  domain is neither product nor a GridSpace,\n";
	  e<<"  or some component is not a GridSpace\n";
	  e<<"  DOMAIN:\n";
	  dom->write(e);
	  throw e;
	}
	if (!grng) {
	  RVLException e;
	  e<<"Error: GridWindowOp constructor\n";
	  e<<"  range is neither product nor a GridSpace,\n";
	  e<<"  or some component is not a GridSpace\n";
	  e<<"  RANGE:\n";
	  bg->getSpace().write(e);
	  throw e;
	}
	if (retrieveGlobalRank()==0) {
	  if (compatible_grid(gdom->getGrid(),grng->getGrid())) {
	    RVLException e;
	    e<<"Error: GridWindowOp constructor\n";
	    e<<"  domain, range defined on incompatible grids\n";
	    e<<"  domain:\n";
	    for (int i=0;i<gdom->getGrid().gdim;i++) 
	      e<<"    axis "<<i<<" d="<<gdom->getGrid().axes[i].d<<" o="<<gdom->getGrid().axes[i].o<<"\n";
	    e<<"  range:\n";
	    for (int i=0;i<grng->getGrid().gdim;i++) 
	      e<<"    axis "<<i<<" d="<<grng->getGrid().axes[i].d<<" o="<<grng->getGrid().axes[i].o<<"\n";
	    throw e;
	  }
	  grid const & g = gdom->getGrid();
	  for (int i=0; i< g.dim; i++) {
	    iw[i]=(int) (w[i]/(g.axes[i].d) + 0.1);
	    iw[i]=iwave_max(iw[i],0);
	  }
	}
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridWindowOp::initialize()\n";
      throw e;
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridWindowOp constructor\n";
      e<<"  either domain or range is neither product nor a GridSpace,\n";
      e<<"  or some component is not a GridSpace\n";
      e<<"  DOMAIN:\n";
      dom->write(e);
      e<<"  RANGE:\n";
      bg->getSpace().write(e);
      throw e;
    }
  }  

  GridWindowOp::GridWindowOp(Space<float> const & _dom,
			     Vector<float> const & _bg,
			     RPNT const w) 
    : dom(RVL::Space<float>::clonePtr(_dom)),
      bg(RVL::Vector<float>::newPtr(_bg.getSpace())) {
    try {

      bg->copy(_bg);
      
      this->initialize(w);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridWindowOp constructor\n";
      throw e;
    }
  }

  GridWindowOp::GridWindowOp(Space<float> const & _dom,
			     std::shared_ptr< Vector<float> > const _bg,
			     RPNT const w) 
    : dom(RVL::Space<float>::clonePtr(_dom)), bg(_bg) {
    try {

      this->initialize(w);
      
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridWindowOp constructor\n";
      throw e;
    }
  }

  void GridWindowOp::apply(Vector<float> const & x,
			   Vector<float> & y) const {
    try {
      GridWindowFO op(iw,true);
      MPISerialFunctionObject<float> mpiop(op);
      y.copy(*bg);
      y.eval(mpiop,x);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridWindowOp::apply\n";
      throw e;
    }
  }

  void GridWindowOp::applyDeriv(Vector<float> const & x,
				Vector<float> const & dx,
				Vector<float> & dy) const {
    try {
      GridWindowFO op(iw,false);
      MPISerialFunctionObject<float> mpiop(op);
      dy.zero();
      dy.eval(mpiop,dx);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridWindowOp::applyDeriv\n";
      throw e;
    }
  }

  void GridWindowOp::applyAdjDeriv(Vector<float> const & x,
				   Vector<float> const & dy,
				   Vector<float> & dx) const {
    try {
      GridWindowFO op(iw,false);
      MPISerialFunctionObject<float> mpiop(op);    
      dx.zero();
      dx.eval(mpiop,dy);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridWindowOp::applyAdjDeriv\n";
      throw e;
    }

  }

  ostream & GridWindowOp::write(ostream & str) const {
    if (!retrieveGlobalRank()) {
      str<<"GridWindowOp\n";
    }
    return str;
  }

  void GridFwdDerivFO::operator()(LocalDataContainer<float> & x,
				  LocalDataContainer<float> const & y) {
    try {
      // cerr<<"GridFwdDerivFO::operator() begin\n";
      ContentPackage< float, RARR > const & gy = 
	dynamic_cast<ContentPackage< float, RARR > const &>(y);

      ContentPackage< float, RARR > & gx = 
	dynamic_cast<ContentPackage< float, RARR > &>(x);

      RARR & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();	  

      // sanity test 
      if (ra_compare_meta(&rax,&ray)) {
	RVLException e;
	e<<"Error: GridFwdDerivFO::operator()\n";
	e<<"  incompatible input, output LDC\n";
	e<<"  input RARR:\n";
	throw e;
      }

      IPNT s0, e0, i;
      ra_a_gse(&rax,s0,e0);

      //      for (int ii=0; ii<RARR_MAX_NDIM; ii++) {
      //	cerr<<"GridDeriv: i="<<ii<<" s="<<s0[ii]<<" e="<<e0[ii]<<endl;
      //      }
#if RARR_MAX_NDIM > 0

      if (dir==0 && rax.ndim==1) {
	rax._s1[e0[0]] = REAL_ZERO;
#pragma ivdep
	for (i[0]=s0[0];i[0]<e0[0];i[0]++) 
	  rax._s1[i[0]] = (ray._s1[i[0]+1]-ray._s1[i[0]])*fac;
      }

#if RARR_MAX_NDIM > 1

      else if (dir==0 && rax.ndim==2) {
	for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	  rax._s2[i[1]][e0[0]] = REAL_ZERO;
#pragma ivdep
	  for (i[0]=s0[0];i[0]<e0[0];i[0]++) 
	    rax._s2[i[1]][i[0]] = (ray._s2[i[1]][i[0]+1]-ray._s2[i[1]][i[0]])*fac;
	}
      }

      else if (dir==1 && rax.ndim==2) {
	for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	  rax._s2[e0[1]][i[0]] = REAL_ZERO;
	  for (i[1]=s0[1];i[1]<e0[1];i[1]++) {
	    rax._s2[i[1]][i[0]] = (ray._s2[i[1]+1][i[0]]-ray._s2[i[1]][i[0]])*fac;
	  }
	}
      }

#if RARR_MAX_NDIM > 2

      else if (dir==0 && rax.ndim==3) {
	for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	    rax._s3[i[2]][i[1]][e0[0]] = REAL_ZERO;
#pragma ivdep
	    for (i[0]=s0[0];i[0]<e0[0];i[0]++) 
	      rax._s3[i[2]][i[1]][i[0]] = (ray._s3[i[2]][i[1]][i[0]+1]-ray._s3[i[2]][i[1]][i[0]])*fac;
	  }
	}
      }
      else if (dir==1 && rax.ndim==3) {
	for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
	    rax._s3[i[2]][e0[1]][i[0]] = REAL_ZERO;
#pragma ivdep
	    for (i[1]=s0[1];i[1]<e0[1];i[1]++) {
	      rax._s3[i[2]][i[1]][i[0]] = (ray._s3[i[2]][i[1]+1][i[0]]-ray._s3[i[2]][i[1]][i[0]])*fac;
	    }	
	  }
	}
      }

      else if (dir==2 && rax.ndim==3) {
	for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	    rax._s3[e0[2]][i[1]][i[0]] = REAL_ZERO;
#pragma ivdep
	    for (i[2]=s0[2];i[2]<e0[2];i[2]++) {
	      rax._s3[i[2]][i[1]][i[0]] = (ray._s3[i[2]+1][i[1]][i[0]]-ray._s3[i[2]][i[1]][i[0]])*fac;
	    }
	  }
	}
      }

#if RARR_MAX_NDIM > 3
	
      else if (dir==0 && rax.ndim==4) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      rax._s4[i[3]][i[2]][i[1]][e0[0]] = REAL_ZERO;
#pragma ivdep
	      for (i[0]=s0[0];i[0]<e0[0];i[0]++) 
		rax._s4[i[3]][i[2]][i[1]][i[0]] = (ray._s4[i[3]][i[2]][i[1]][i[0]+1]-
						   ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;
	    }
	  }
	}
      }
      else if (dir==1 && rax.ndim==4) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
	      rax._s4[i[3]][i[2]][e0[1]][i[0]] = REAL_ZERO;
#pragma ivdep
	      for (i[1]=s0[1];i[1]<e0[1];i[1]++) {
		rax._s4[i[3]][i[2]][i[1]][i[0]] = (ray._s4[i[3]][i[2]][i[1]+1][i[0]]-
						   ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;
	      }	
	    }
	  }
	}
      }
      else if (dir==2 && rax.ndim==4) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
	      rax._s4[i[3]][e0[2]][i[1]][i[0]] = REAL_ZERO;
#pragma ivdep
	      for (i[2]=s0[2];i[2]<e0[2];i[2]++) {
		rax._s4[i[3]][i[2]][i[1]][i[0]] = (ray._s4[i[3]][i[2]+1][i[1]][i[0]]-
						   ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;
	      }	
	    }
	  }
	}
      }

      else if (dir==3 && rax.ndim==4) {
	for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
	      rax._s4[e0[3]][i[2]][i[1]][i[0]] = REAL_ZERO;
#pragma ivdep
	      for (i[3]=s0[3];i[3]<e0[3];i[3]++) {
		rax._s4[i[3]][i[2]][i[1]][i[0]] = (ray._s4[i[3]+1][i[2]][i[1]][i[0]]-
						   ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;
	      }	
	    }
	  }
	}
      }

#if RARR_MAX_NDIM > 4

      else if (dir==0 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
		rax._s5[i[4]][i[3]][i[2]][i[1]][e0[0]] = REAL_ZERO;
#pragma ivdep
		for (i[0]=s0[0];i[0]<e0[0];i[0]++) 
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = (ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]+1]-
							   ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;
	      }
	    }
	  }
	}
      }
      else if (dir==1 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
		rax._s5[i[4]][i[3]][i[2]][e0[1]][i[0]] = REAL_ZERO;
#pragma ivdep
		for (i[1]=s0[1];i[1]<e0[1];i[1]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = (ray._s5[i[4]][i[3]][i[2]][i[1]+1][i[0]]-
							   ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;
		}	
	      }
	    }
	  }
	}
      }
      else if (dir==2 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
		rax._s5[i[4]][i[3]][e0[2]][i[1]][i[0]] = REAL_ZERO;
#pragma ivdep
		for (i[2]=s0[2];i[2]<e0[2];i[2]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = (ray._s5[i[4]][i[3]][i[2]+1][i[1]][i[0]]-
							   ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;
		}	
	      }
	    }
	  }
	}
      }
      else if (dir==3 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
		rax._s5[i[4]][e0[3]][i[2]][i[1]][i[0]] = REAL_ZERO;
#pragma ivdep
		for (i[3]=s0[3];i[3]<e0[3];i[3]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = (ray._s5[i[4]][i[3]+1][i[2]][i[1]][i[0]]-
							   ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;
		}	
	      }
	    }
	  }
	}
      }

      else if (dir==4 && rax.ndim==5) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) { 
		rax._s5[e0[4]][i[3]][i[2]][i[1]][i[0]] = REAL_ZERO;
#pragma ivdep
		for (i[4]=s0[4];i[4]<e0[4];i[4]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = (ray._s5[i[4]+1][i[3]][i[2]][i[1]][i[0]]-
							   ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;
		}	
	      }
	    }
	  }
	}
      }

#endif // RARR_MAX_NDIM > 4
#endif // > 3
#endif // > 2
#endif // > 1

      else {
	RVLException e;
	e<<"Error: GridFwdDerivFO::operator()\n";
	e<<"  attempt to apply divided diff in direction "<<dir<<" on array of dim "<<rax.ndim<<"\n";
	throw e;
      }
#endif // > 0
      // cerr<<"GridFwdDerivFO::operator() end\n";
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridFwdDerivFO::operator()\n";
      e<<"input type error - not CP<float,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridFwdDerivFO::operator()\n";
      throw e;
    }
  }

  void GridAdjDerivFO::operator()(LocalDataContainer<float> & x,
				  LocalDataContainer<float> const & y) {
    try {
      // cerr<<"GridAdjDerivFO::operator() begin\n";
      ContentPackage< float, RARR > const & gy = 
	dynamic_cast<ContentPackage< float, RARR > const &>(y);

      ContentPackage< float, RARR > & gx = 
	dynamic_cast<ContentPackage< float, RARR > &>(x);

      RARR & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();	  

      // sanity test 
      if (ra_compare_meta(&rax,&ray)) {
	RVLException e;
	e<<"Error: GridAdjDerivFO::operator()\n";
	e<<"  incompatible input, output LDC\n";
	throw e;
      }

      IPNT s0, e0, i;
      ra_a_gse(&rax,s0,e0);
      
#if RARR_MAX_NDIM > 0

      if (dir==0 && rax.ndim==1) {
	rax._s1[s0[0]] = -ray._s1[s0[0]]*fac;
	rax._s1[e0[0]] = ray._s1[e0[0]-1]*fac;
#pragma ivdep
	for (i[0]=s0[0]+1;i[0]<e0[0];i[0]++) 
	  rax._s1[i[0]]=(ray._s1[i[0]-1]-ray._s1[i[0]])*fac;      }

#if RARR_MAX_NDIM > 1

      else if (dir==0 && rax.ndim==2) {
	for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	  rax._s2[i[1]][s0[0]] = -ray._s2[i[1]][s0[0]]*fac;
	  rax._s2[i[1]][e0[0]] = ray._s2[i[1]][e0[0]-1]*fac;
#pragma ivdep
	  for (i[0]=s0[0]+1;i[0]<e0[0];i[0]++) {
	    rax._s2[i[1]][i[0]]=(ray._s2[i[1]][i[0]-1]-ray._s2[i[1]][i[0]])*fac;      
	  }
	}
      }
      
      else if (dir==1 && rax.ndim==2) {
	for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	  rax._s2[s0[1]][i[0]] = -ray._s2[s0[1]][i[0]]*fac;
	  rax._s2[e0[1]][i[0]] = ray._s2[e0[1]-1][i[0]]*fac;
	  for (i[1]=s0[1]+1;i[1]<e0[1];i[1]++) {
	    rax._s2[i[1]][i[0]]=(ray._s2[i[1]-1][i[0]]-ray._s2[i[1]][i[0]])*fac;      
	  }
	}
      }

#if RARR_MAX_NDIM > 2

      else if (dir==0 && rax.ndim==3) {
	for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	    rax._s3[i[2]][i[1]][s0[0]] = -ray._s3[i[2]][i[1]][s0[0]]*fac;
	    rax._s3[i[2]][i[1]][e0[0]] = ray._s3[i[2]][i[1]][e0[0]-1]*fac;
#pragma ivdep
	    for (i[0]=s0[0]+1;i[0]<e0[0];i[0]++) {
	      rax._s3[i[2]][i[1]][i[0]]=(ray._s3[i[2]][i[1]][i[0]-1]-ray._s3[i[2]][i[1]][i[0]])*fac;      
	    }
	  }
	}
      }

      else if (dir==1 && rax.ndim==3) {
	for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	    rax._s3[i[2]][s0[1]][i[0]] = -ray._s3[i[2]][s0[1]][i[0]]*fac;
	    rax._s3[i[2]][e0[1]][i[0]] = ray._s3[i[2]][e0[1]-1][i[0]]*fac;
#pragma ivdep
	    for (i[1]=s0[1]+1;i[1]<e0[1];i[1]++) {
	      rax._s3[i[2]][i[1]][i[0]]=(ray._s3[i[2]][i[1]-1][i[0]]-ray._s3[i[2]][i[1]][i[0]])*fac;      
	    }
	  }
	}
      }

      else if (dir==2 && rax.ndim == 3) {
	for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	    rax._s3[s0[2]][i[1]][i[0]] = -ray._s3[s0[2]][i[1]][i[0]]*fac;
	    rax._s3[e0[2]][i[1]][i[0]] = ray._s3[e0[2]-1][i[1]][i[0]]*fac;
#pragma ivdep
	    for (i[2]=s0[2]+1;i[2]<e0[2];i[2]++) {
	      rax._s3[i[2]][i[1]][i[0]]=(ray._s3[i[2]-1][i[1]][i[0]]-ray._s3[i[2]][i[1]][i[0]])*fac;      
	    }
	  }
	}
      }

#if RARR_MAX_NDIM > 3
	
      else if (dir==0 && rax.ndim==4) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      rax._s4[i[3]][i[2]][i[1]][s0[0]] = -ray._s4[i[3]][i[2]][i[1]][s0[0]]*fac;
	      rax._s4[i[3]][i[2]][i[1]][e0[0]] = ray._s4[i[3]][i[2]][i[1]][e0[0]-1]*fac;
#pragma ivdep
	      for (i[0]=s0[0]+1;i[0]<e0[0];i[0]++) {
		rax._s4[i[3]][i[2]][i[1]][i[0]]=(ray._s4[i[3]][i[2]][i[1]][i[0]-1]-
						 ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;      
	      }
	    }
	  }
	}
      }
      else if (dir==1 && rax.ndim==4) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      rax._s4[i[3]][i[2]][s0[1]][i[0]] = -ray._s4[i[3]][i[2]][s0[1]][i[0]]*fac;
	      rax._s4[i[3]][i[2]][e0[1]][i[0]] = ray._s4[i[3]][i[2]][e0[1]-1][i[0]]*fac;
#pragma ivdep
	      for (i[1]=s0[1]+1;i[1]<e0[1];i[1]++) {
		rax._s4[i[3]][i[2]][i[1]][i[0]]=(ray._s4[i[3]][i[2]][i[1]-1][i[0]]-
						 ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;      
	      }
	    }
	  }
	}
      }
      else if (dir==2 && rax.ndim==4) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      rax._s4[i[3]][s0[2]][i[1]][i[0]] = -ray._s4[i[3]][s0[2]][i[1]][i[0]]*fac;
	      rax._s4[i[3]][e0[2]][i[1]][i[0]] = ray._s4[i[3]][e0[2]-1][i[1]][i[0]]*fac;
#pragma ivdep
	      for (i[2]=s0[2]+1;i[2]<e0[2];i[2]++) {
		rax._s4[i[3]][i[2]][i[1]][i[0]]=(ray._s4[i[3]][i[2]-1][i[1]][i[0]]-
						 ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;      
	      }
	    }
	  }
	}
      }

      else if (dir==3 && rax.ndim==4) {
	for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      rax._s4[s0[3]][i[2]][i[1]][i[0]] = -ray._s4[s0[3]][i[2]][i[1]][i[0]]*fac;
	      rax._s4[e0[3]][i[2]][i[1]][i[0]] = ray._s4[e0[3]-1][i[2]][i[1]][i[0]]*fac;
#pragma ivdep
	      for (i[3]=s0[3]+1;i[3]<e0[3];i[3]++) {
		rax._s4[i[3]][i[2]][i[1]][i[0]]=(ray._s4[i[3]-1][i[2]][i[1]][i[0]]-
						 ray._s4[i[3]][i[2]][i[1]][i[0]])*fac;      
	      }
	    }
	  }
	}
      }

#if RARR_MAX_NDIM > 4

      else if (dir==0 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
		rax._s5[i[4]][i[3]][i[2]][i[1]][s0[0]] = -ray._s5[i[4]][i[3]][i[2]][i[1]][s0[0]]*fac;
		rax._s5[i[4]][i[3]][i[2]][i[1]][e0[0]] = ray._s5[i[4]][i[3]][i[2]][i[1]][e0[0]-1]*fac;
#pragma ivdep
		for (i[0]=s0[0]+1;i[0]<e0[0];i[0]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]]=(ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]-1]-
							 ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;      
		}
	      }
	    }
	  }
	}
      }
      else if (dir==1 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s5[i[4]][i[3]][i[2]][s0[1]][i[0]] = -ray._s5[i[4]][i[3]][i[2]][s0[1]][i[0]]*fac;
		rax._s5[i[4]][i[3]][i[2]][e0[1]][i[0]] = ray._s5[i[4]][i[3]][i[2]][e0[1]-1][i[0]]*fac;
#pragma ivdep
		for (i[1]=s0[1]+1;i[1]<e0[1];i[1]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]]=(ray._s5[i[4]][i[3]][i[2]][i[1]-1][i[0]]-
							 ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;      
		}
	      }
	    }
	  }
	}
      }
      else if (dir==2 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s5[i[4]][i[3]][s0[2]][i[1]][i[0]] = -ray._s5[i[4]][i[3]][s0[2]][i[1]][i[0]]*fac;
		rax._s5[i[4]][i[3]][e0[2]][i[1]][i[0]] = ray._s5[i[4]][i[3]][e0[2]-1][i[1]][i[0]]*fac;
#pragma ivdep
		for (i[2]=s0[2]+1;i[2]<e0[2];i[2]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]]=(ray._s5[i[4]][i[3]][i[2]-1][i[1]][i[0]]-
							 ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;      
		}
	      }
	    }
	  }
	}
      }
      else if (dir==3 && rax.ndim==5) {
	for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s5[i[4]][s0[3]][i[2]][i[1]][i[0]] = -ray._s5[i[4]][s0[3]][i[2]][i[1]][i[0]]*fac;
		rax._s5[i[4]][e0[3]][i[2]][i[1]][i[0]] = ray._s5[i[4]][e0[3]-1][i[2]][i[1]][i[0]]*fac;
#pragma ivdep
		for (i[3]=s0[3]+1;i[3]<e0[3];i[3]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]]=(ray._s5[i[4]][i[3]-1][i[2]][i[1]][i[0]]-
							 ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;      
		}
	      }
	    }
	  }
	}
      }

      else if (dir==4 && rax.ndim==5) {
	for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	  for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s5[s0[4]][i[3]][i[2]][i[1]][i[0]] = -ray._s5[s0[4]][i[3]][i[2]][i[1]][i[0]]*fac;
		rax._s5[e0[4]][i[3]][i[2]][i[1]][i[0]] = ray._s5[e0[4]-1][i[3]][i[2]][i[1]][i[0]]*fac;
#pragma ivdep
		for (i[4]=s0[4]+1;i[4]<e0[4];i[4]++) {
		  rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]]=(ray._s5[i[4]-1][i[3]][i[2]][i[1]][i[0]]-
							 ray._s5[i[4]][i[3]][i[2]][i[1]][i[0]])*fac;      
		}
	      }
	    }
	  }
	}	
      }

#endif // RARR_MAX_NDIM > 4
#endif // > 3
#endif // > 2
#endif // > 1

      else {
	RVLException e;
	e<<"Error: GridAdjDerivFO::operator()\n";
	e<<"  attempt to apply divided diff in direction "<<dir<<" on array of dim "<<rax.ndim<<"\n";
	throw e;
      }
#endif // > 0
      // cerr<<"GridAdjDerivFO::operator() end\n";
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridAdjDerivFO::operator()\n";
      e<<"input type error - not CP<float,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridAdjDerivFO::operator()\n";
      throw e;
    }
  }

  GridDerivOp::GridDerivOp(Space<float> const & _dom,
			   int _dir, float scale)
    : dir(_dir), fac(0), dom(_dom) {
    try {
      ProductSpace<float> const * pdom = NULL;
      pdom = dynamic_cast<ProductSpace<float> const *>(&dom);
      int n_fac=1;
      if (pdom) n_fac=pdom->getSize();
      Space<float> const * sp = NULL;
      for (int j=0; j<n_fac; j++) {
	if (pdom) sp = &((*pdom)[j]);
	else sp = &dom;
	myGridSpace const * gdom = dynamic_cast<myGridSpace const *>(sp);
	if (!gdom) {
	  RVLException e;
	  e<<"Error: GridDerivOp constructor\n";
	  e<<"  factor "<<j<<" of input space is not a GridSpace\n";
	  e<<"  description:\n";
	  sp->write(e);
	  throw e;	  
	}
	// pure out of core: real factors only on rk=0
	if (retrieveGlobalRank()==0) {
	  if (dir < 0 || dir > gdom->getGrid().gdim-1) {
	    RVLException e;
	    e<<"Error: GridDerivOp constructor\n";
	    e<<"  direction index "<<dir<<" out of dimension range [0,"<<gdom->getGrid().gdim-1<<"\n";
	    throw e;
	  }
	  RPNT d;
	  get_d(d,gdom->getGrid());
	  fac.push_back(scale/d[dir]);
	}
	else {
	  fac.push_back(REAL_ZERO);
	}
      }
    }
    catch (RVLException e) {
      e<<"\ncalled from GridDerivOp constructor\n";
      throw e;
    }
  }

  GridDerivOp::GridDerivOp(GridDerivOp const & op)
    : dir(op.dir), fac(op.fac), dom(op.dom) {}

  void GridDerivOp::apply(Vector<float> const & x,
			  Vector<float> & y) const {
    try {
      Components<float> cx(x);
      Components<float> cy(y);
      for (int j=0;j<(int)cx.getSize();j++) {
	GridFwdDerivFO f(dir,fac[j]);
	MPISerialFunctionObject<float> mpif(f);
    	cy[j].eval(mpif,cx[j]);
      }
#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridDerivOp::apply\n";
      throw e;
    }
  }

  void GridDerivOp::applyAdj(Vector<float> const & x,
			     Vector<float> & y) const {
    try {
      Components<float> cx(x);
      Components<float> cy(y);
      for (int j=0;j<(int)cx.getSize();j++) {
	GridAdjDerivFO f(dir,fac[j]);
	MPISerialFunctionObject<float> mpif(f);
	cy[j].eval(mpif,cx[j]);
      }
#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridDerivOp::applyAdj\n";
      throw e;
    }
  }

  ostream & GridDerivOp::write(ostream & str) const {
    str<<"GridDerivOp: directional derivative, axis = "<<dir<<"\n";
    str<<"Domain:\n";
    dom.write(str);
    return str;
  }

  void GridFwdExtendFO::operator()(LocalDataContainer<float> & x,
				   LocalDataContainer<float> const & y) {
    try {
      // cerr<<"GridFwdExtendFO::operator() begin\n";
      ContentPackage< float, RARR > const & gy = 
	dynamic_cast<ContentPackage< float, RARR > const &>(y);

      ContentPackage< float, RARR > & gx = 
	dynamic_cast<ContentPackage< float, RARR > &>(x);

      RARR & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();	  

      // presumption: these are called only from GridExtendOp, which
      // does all sanity testing - could make them private data
      // exception: number of extended axes can be checked!!
      if (rax.ndim-ray.ndim != n_ext) {
	RVLException e;
	e<<"Error: GridExtendFwdFO::operator()\n"; 
	e<<"  dimension difference of input, output does not match\n";
	e<<"  number of extended axes\n";
	throw e;
      }

      IPNT s0, e0, i;
      ra_a_gse(&rax,s0,e0);

      // zero if internal
      if (!ext) {
	// sanity check - must have zero section in all internal extended axes
	for (int ii=ray.ndim;ii<rax.ndim;ii++) {
	  if (s0[ii] > 0 || e0[ii] < 0) {
	    RVLException e;
	    e<<"Error: GridExtendFwdFO::operator()\n"; 
	    e<<"  index range on axis "<<ii<<" of output extended array\n";
	    e<<"  = ["<<s0[ii]<<","<<e0[ii]<<"], does not contain zero as is\n";
	    e<<"  required for internal extended axes\n";
	    throw e;
	  }
	}
	ra_a_zero(&rax);
      }

      if (ray.ndim==1) {

#if RARR_MAX_NDIM > 1

	if (rax.ndim==2) {
	  if (ext) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s2[i[1]][i[0]] = ray._s0[i[0]];
	      }
	    }
	  }
	  else {
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
#pragma ivdep
	      rax._s2[0][i[0]] = fac*ray._s0[i[0]];
	    }
	  }
	}

#if RARR_MAX_NDIM > 2

	if (rax.ndim==3) {
	  if (ext) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  rax._s3[i[2]][i[1]][i[0]] = ray._s0[i[0]];
		}
	      }
	    }
	  }
	  else {
#pragma ivdep
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      rax._s3[0][0][i[0]] = fac*ray._s0[i[0]];
	    }
	  }
	}

#if RARR_MAX_NDIM > 3

	if (rax.ndim==4) {
	  if (ext) {
	    for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	      for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		    rax._s4[i[3]][i[2]][i[1]][i[0]] = ray._s0[i[0]];
		  }
		}
	      }
	    }
	  }
	  else {
#pragma ivdep
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      rax._s4[0][0][0][i[0]] = fac*ray._s0[i[0]];
	    }
	  }
	}

#if RARR_MAX_NDIM > 4

	if (rax.ndim==5) {
	  if (ext) {
	    for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	      for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
		for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		      rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = ray._s0[i[0]];
		    }
		  }
		}
	      }
	    }
	  }
	  else {
#pragma ivdep
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      rax._s5[0][0][0][0][i[0]] = fac*ray._s0[i[0]];
	    }
	  }
	}

#endif
#endif
#endif
#endif
      }

      else if (ray.ndim==2) {

#if RARR_MAX_NDIM > 2

	if (rax.ndim==3) {
	  if (ext) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  rax._s3[i[2]][i[1]][i[0]] = ray._s2[i[1]][i[0]];
		}
	      }
	    }
	  }
	  else {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s3[0][i[1]][i[0]] = fac*ray._s2[i[1]][i[0]];
	      }
	    }
	  }
	}

#if RARR_MAX_NDIM > 3

	if (rax.ndim==4) {
	  if (ext) {
	    for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	      for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		    rax._s4[i[3]][i[2]][i[1]][i[0]] = ray._s2[i[1]][i[0]];
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s4[0][0][i[1]][i[0]] = fac*ray._s2[i[1]][i[0]];
	      }
	    }
	  }
	}

#if RARR_MAX_NDIM > 4

	if (rax.ndim==5) {
	  if (ext) {
	    for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	      for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
		for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		      rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = ray._s2[i[1]][i[0]];
		    }
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		rax._s5[0][0][0][i[1]][i[0]] = fac*ray._s2[i[1]][i[0]];
	      }
	    }
	  }
	}

#endif
#endif
#endif

      }

      else if (ray.ndim==3) {

#if RARR_MAX_NDIM > 3

	if (rax.ndim==4) {
	  if (ext) {
	    for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	      for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		    rax._s4[i[3]][i[2]][i[1]][i[0]] = ray._s3[i[2]][i[1]][i[0]];
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  rax._s4[0][i[2]][i[1]][i[0]] = fac*ray._s3[i[2]][i[1]][i[0]];
		}
	      }
	    }
	  }
	}

#if RARR_MAX_NDIM > 4

	if (rax.ndim==5) {
	  if (ext) {
	    for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	      for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
		for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		      rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]] = ray._s3[i[2]][i[1]][i[0]];
		    }
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  rax._s5[0][0][i[2]][i[1]][i[0]] = fac*ray._s3[i[2]][i[1]][i[0]];
		}
	      }
	    }
	  }
	}

#endif
#endif

      }

      else {
	RVLException e;
	e<<"Error: GridFwdExtendFO::operator()\n";
	e<<"  input rarr dimension not 1, 2, or 3 - only permitted spatial dims\n";
	throw e;
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridFwdExtendFO::operator()\n";
      e<<"  input type error - not CP<float,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridFwdExtendFO::operator()\n";
      throw e;
    }
  }

  void GridAdjExtendFO::operator()(LocalDataContainer<float> & y,
				   LocalDataContainer<float> const & x) {
    try {

      ContentPackage< float, RARR > const & gx = 
	dynamic_cast<ContentPackage< float, RARR > const &>(x);

      ContentPackage< float, RARR > & gy = 
	dynamic_cast<ContentPackage< float, RARR > &>(y);

      RARR & ray = gy.getMetadata();
      RARR const & rax = gx.getMetadata();	  

      // presumption: these are called only from GridExtendOp, which
      // does all sanity testing - could make them private data
      // exception: number of extended axes can be checked!!
      if (rax.ndim-ray.ndim != n_ext) {
	RVLException e;
	e<<"Error: GridExtendAdjFO::operator()\n"; 
	e<<"  dimension difference of input, output does not match\n";
	e<<"  number of extended axes\n";
	throw e;
      }

      IPNT s0, e0, i;
      ra_a_gse(&rax,s0,e0);

      if (!ext) {
	// sanity check - must have zero section in all internal extended axes
	for (int ii=ray.ndim;ii<rax.ndim;ii++) {
	  if (s0[ii] > 0 || e0[ii] < 0) {
	    RVLException e;
	    e<<"Error: GridExtendAdjFO::operator()\n"; 
	    e<<"  index range on axis "<<ii<<" of output extended array\n";
	    e<<"  = ["<<s0[ii]<<","<<e0[ii]<<"], does not contain zero as is\n";
	    e<<"  required for internal extended axes\n";
	    throw e;
	  }
	}
      }

      // zero output anyway
      ra_a_zero(&ray);

      if (ray.ndim==1) {

#if RARR_MAX_NDIM > 1

	if (rax.ndim==2) {
	  if (ext) {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		ray._s0[i[0]] += fac*rax._s2[i[1]][i[0]];
	      }
	    }
	  }
	  else {
#pragma ivdep
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      ray._s0[i[0]] = rax._s2[0][i[0]];
	    }
	  }
	}

#if RARR_MAX_NDIM > 2

	if (rax.ndim==3) {
	  if (ext) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  ray._s0[i[0]]+= fac*rax._s3[i[2]][i[1]][i[0]];
		}
	      }
	    }
	  }
	  else {
#pragma ivdep
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      ray._s0[i[0]]=rax._s3[0][0][i[0]];
	    }
	  }
	}

#if RARR_MAX_NDIM > 3

	if (rax.ndim==4) {
	  if (ext) {
	    for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	      for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		    ray._s0[i[0]] += fac* rax._s4[i[3]][i[2]][i[1]][i[0]] ;
		  }
		}
	      }
	    }
	  }
	  else {
#pragma ivdep
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      ray._s0[i[0]]=rax._s4[0][0][0][i[0]];
	    }
	  }
	}

#if RARR_MAX_NDIM > 4

	if (rax.ndim==5) {
	  if (ext) {
	    for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	      for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
		for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		      ray._s0[i[0]] += fac* rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]];
		    }
		  }
		}
	      }
	    }
	  }
	  else {
#pragma ivdep
	    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	      ray._s0[i[0]]=rax._s5[0][0][0][0][i[0]];
	    }
	  }
	}

#endif
#endif
#endif
#endif
      }

      else if (ray.ndim==2) {

#if RARR_MAX_NDIM > 2

	if (rax.ndim==3) {
	  if (ext) {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  ray._s2[i[1]][i[0]] += fac* rax._s3[i[2]][i[1]][i[0]];
		}
	      }
	    }
	  }
	  else {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		ray._s2[i[1]][i[0]] = rax._s3[0][i[1]][i[0]];
	      }
	    }
	  }
	}

#if RARR_MAX_NDIM > 3

	if (rax.ndim==4) {
	  if (ext) {
	    for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	      for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		    ray._s2[i[1]][i[0]] += fac* rax._s4[i[3]][i[2]][i[1]][i[0]];
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		ray._s2[i[1]][i[0]] = rax._s4[0][0][i[1]][i[0]];
	      }
	    }
	  }
	}

#if RARR_MAX_NDIM > 4

	if (rax.ndim==5) {
	  if (ext) {
	    for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	      for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
		for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		      ray._s2[i[1]][i[0]] += fac* rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]];
		    }
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
	      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		ray._s2[i[1]][i[0]] = rax._s5[0][0][0][i[1]][i[0]];
	      }
	    }
	  }
	}

#endif
#endif
#endif

      }

      else if (ray.ndim==3) {

#if RARR_MAX_NDIM > 3

	if (rax.ndim==4) {
	  if (ext) {
	    for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
	      for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		  for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		    ray._s3[i[2]][i[1]][i[0]] += fac* rax._s4[i[3]][i[2]][i[1]][i[0]];
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  ray._s3[i[2]][i[1]][i[0]] = rax._s4[0][i[2]][i[1]][i[0]];
		}
	      }
	    }
	  }
	}

#if RARR_MAX_NDIM > 4

	if (rax.ndim==5) {
	  if (ext) {
	    for (i[4]=s0[4];i[4]<=e0[4];i[4]++) {
	      for (i[3]=s0[3];i[3]<=e0[3];i[3]++) {
		for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
		  for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		    for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		      ray._s3[i[2]][i[1]][i[0]] += fac* rax._s5[i[4]][i[3]][i[2]][i[1]][i[0]];
		    }
		  }
		}
	      }
	    }
	  }
	  else {
	    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
	      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
#pragma ivdep
		for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
		  ray._s3[i[2]][i[1]][i[0]] = rax._s5[0][0][i[2]][i[1]][i[0]];
		}
	      }
	    }
	  }
	}

#endif
#endif

      }

      else {
	RVLException e;
	e<<"Error: GridAdjExtendFO::operator()\n";
	e<<"  input rarr dimension not 1, 2, or 3 - only permitted spatial dims\n";
	throw e;
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridAdjExtendFO::operator()\n";
      e<<"  input type error - not CP<float,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridAdjExtendFO::operator()\n";
      throw e;
    }
  }

  
  GridExtendOp::GridExtendOp(Space<float> const & _dom,
			     Space<float> const & _rng)
    : dom(_dom), rng(_rng), n_ext(0), ext(0), fac(0) {
    try {

      // step 1: compatibility of product structures
      ProductSpace<float> const * pdom = NULL;
      ProductSpace<float> const * prng = NULL;
      pdom = dynamic_cast<ProductSpace<float> const *>(&dom);
      prng = dynamic_cast<ProductSpace<float> const *>(&rng);
      int n_dom=1;
      int n_rng=1;
      if (pdom) n_dom=pdom->getSize();
      if (prng) n_rng=prng->getSize();
      if (n_dom != n_rng) {
	RVLException e;
	e<<"Error: GridExtendOp constructor\n";
	e<<"  domain, range spaces have differing numbers of factors\n";
	e<<"  domain space:\n";
	dom.write(e);
	e<<"  range space:\n";
	rng.write(e);
	throw e;
      }

      // step 2: every factor of range is really extension of 
      // corresponding factor in domain, and domain factors are
      // all simple spatial grids with no extended axes
      Space<float> const * dsp = NULL;
      Space<float> const * rsp = NULL;
      for (int j=0; j<n_dom; j++) {
	if (pdom) dsp = &((*pdom)[j]);
	else dsp=&dom;
	if (prng) rsp = &((*prng)[j]); 
	else rsp=&rng;
	myGridSpace const & gdom = dynamic_cast<myGridSpace const &>(*dsp);    
	myGridSpace const & grng = dynamic_cast<myGridSpace const &>(*rsp);    
	if (retrieveGlobalRank()==0) {
	  // temporary limitation: only incore range
	  if (!(grng.isIncore())) {
	    RVLException e;
	    e<<"Error: GridExtendOp constructor\n";
	    e<<"  only incore range space allowed\n";
	    e<<"  this range space:\n";
	    rsp->write(e);
	    throw e;
	  }
	  // no extended axes in domain
	  if (gdom.getGrid().gdim != gdom.getGrid().dim) {
	    RVLException e;
	    e<<"Error: GridExtendOp constructor\n";
	    e<<"  domain factor "<<j<<" has extended axes - not permitted\n";
	    e<<"  this domain factor:\n";
	    dsp->write(e);
	    throw e;
	  }
	  // range is extension of domain
	  if (gdom.getGrid().dim != grng.getGrid().dim) {
	    RVLException e;
	    e<<"Error: GridExtendOp constructor\n";
	    e<<"  spatial dims (keyword dim) of domain, range factors "<<j<<"\n";
	    e<<"  differ - must be same\n";
	    e<<"  domain factor "<<j<<":\n";
	    dsp->write(e);
	    e<<"  range factor "<<j<<":\n";
	    rsp->write(e);
	    throw e;
	  }
	  // find spatial axes, must be in same order, and come before any extended
	  // also must be geometrically same
	  for (int i=0;i<gdom.getGrid().dim;i++) {
	    int idom=-1;
	    int irng=-1;
	    for (int k=0;k<gdom.getGrid().gdim;k++) 
	      if (gdom.getGrid().axes[k].id == i) idom=k;
	    for (int k=0;k<grng.getGrid().gdim;k++) 
	      if (grng.getGrid().axes[k].id == i) irng=k;
	    if (idom<0 || irng<0 || idom > gdom.getGrid().dim-1 || idom != irng ) {
	      RVLException e;
	      e<<"Error: GridExtendOp constructor\n";
	      if (idom<0 || irng<0)
		e<<"  failed to identify axis in domain or range with id="<<i<<"\n";
	      if (idom>gdom.getGrid().dim-1)
		e<<"  spatial grid axes must be ordered first - but axis id="<<i<<" has index "<<idom<<"\n";
	      if (idom != irng)
		e<<"  indices for id="<<i<<" differ: dom="<<idom<<" rng="<<irng<<"\n";
	      throw e;
	    }
	    if ((compare_axis(gdom.getGrid().axes[idom],grng.getGrid().axes[irng]))) {
	      RVLException e;
	      e<<"Error: GridExtendOp constructor\n";
	      e<<"  axes differ for spatial axis id = "<<i<<" dom = rng axis index = "<<irng<<"\n";
	      throw e;
	    }
	  }
	  // determine number of extended axes, whether external or internal, and 
	  // scale factor
	  float tmpfac = REAL_ONE;
	  // set tmp external flag by first external axis, if there is one
	  int tmp_n_ext = grng.getGrid().gdim - grng.getGrid().dim;
	  bool tmpext = false;
	  if (tmp_n_ext > 0) {
	    tmpext = (grng.getGrid().axes[grng.getGrid().dim].id < EXTINT);
	  }
	  for (int i=gdom.getGrid().dim;i<grng.getGrid().gdim;i++) {
	    if (grng.getGrid().axes[i].id > 99) {
	      if (tmpext) {
		RVLException e;
		e<<"Error: GridExtendOp constructor\n";
		e<<"  range space mixes internal, external extended axes\n";
		e<<"  not permitted in current design\n";
		//		grng.write(e);
		throw e;
	      }
	      float newtmpfac = REAL_ONE;
	      if (ProtectedDivision<float>(REAL_ONE,grng.getGrid().axes[i].d,newtmpfac)) {
		RVLException e;
		e<<"Error: GridExtendOp constructor\n";
		e<<"  zerodivide by cell vol in axis="<<i<<" id="
		 <<grng.getGrid().axes[i].id<<" of GridSpace:\n";
		//		grng.write(e);
		throw e;
	      }
	      tmpfac*=newtmpfac;
	    }
	    if (grng.getGrid().axes[i].id < 100) {
	      if (!tmpext) {
		RVLException e;
		e<<"Error: GridExtendOp constructor\n";
		e<<"  range space mixes internal, external extended axes\n";
		e<<"  not permitted in current design\n";
		//		grng.write(e);
		throw e;
	      }
	      tmpfac*=grng.getGrid().axes[i].d;
	    }
	  }
	  n_ext.push_back(tmp_n_ext);
	  ext.push_back(tmpext);
	  // TEMPORARY: do not do quadrature on internal extd axes
	  fac.push_back(tmpfac);
	  //	  fac.push_back(REAL_ONE);
	}
	else {
	  n_ext.push_back(0);
	  ext.push_back(false);
	  fac.push_back(REAL_ONE);
	}
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridExtendOp constructor\n";
      e<<"  at least one factor in input space is not a GridSpace\n";
      throw e;
    }
  }

  GridExtendOp::GridExtendOp(GridExtendOp const & op)
    : dom(op.dom), rng(op.rng), n_ext(op.n_ext), ext(op.ext), fac(op.fac) {}

  void GridExtendOp::apply(Vector<float> const & x,
			   Vector<float> & y) const {
    try {
      // note that LinearOp base class takes care of space
      // membership tests, hence sanity tests
      // here x is extended to y
      Components<float> cx(x);
      Components<float> cy(y);
      for (int j=0;j<(int)cx.getSize();j++) {
	if (n_ext[j] <= 0) {
	  cy[j].copy(cx[j]);
	}
	else {
	  GridFwdExtendFO f(n_ext[j], ext[j], fac[j]);
	  MPISerialFunctionObject<float> mpif(f);
	  cy[j].eval(mpif,cx[j]);
	}
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridExtendOp::apply\n";
      throw e;
    }
  }

  void GridExtendOp::applyAdj(Vector<float> const & x,
			      Vector<float> & y) const {
    try {
      Components<float> cx(x);
      Components<float> cy(y);
      for (int j=0;j<(int)cx.getSize();j++) {
	if (n_ext[j] <= 0) {
	  cy[j].copy(cx[j]);
	}
	else {
	  GridAdjExtendFO f(n_ext[j], ext[j], fac[j]);
	  MPISerialFunctionObject<float> mpif(f);
	  cy[j].eval(mpif,cx[j]);
	}
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridExtendOp::applyAdj\n";
      throw e;
    }
  }

  ostream & GridExtendOp::write(ostream & str) const {
    str<<"GridExtendOp: inject spatial grid into extended grid\n";
    str<<"Domain:\n";
    dom.write(str);
    str<<"Range:\n";
    rng.write(str);
    return str;
  }

  /*
    void HelmFO::operator()(LocalDataContainer<float> & x,
    LocalDataContainer<float> const & y){
    try{
    float *indata=NULL;
    float *outdata=NULL;
    float *work=NULL;
    integer f2c_n1;
    integer f2c_n2;
    integer lenwork;
    ContentPackage<float, RARR>  & gx =
    dynamic_cast<ContentPackage <float, RARR>  &> (x);
    ContentPackage<float, RARR> const & gy =
    dynamic_cast<ContentPackage <float, RARR> const &> (y);
        
    // precondition - metadata are same dimn
    RARR  & rax = gx.getMetadata();
    RARR const & ray = gy.getMetadata();
    int dimx; int dimy;
    int lendom;
    ra_ndim(&rax,&dimx);
    ra_ndim(&ray,&dimy);
    //cerr << "\n xdim=" << dimx << endl;
    //cerr << "\n ydim=" << dimy << endl;
    if (dimx != dimy) {
    RVLException e;
    e<<"Error: HelmFO::operator()\n";
    e<<"arguments have different dims:\n";
    e<<"dimx="<<dimx<<" dimy="<<dimy<<"\n";
    throw e;
    }
        
    // compute grid params
    IPNT gsx; IPNT gex;
    IPNT gsy; IPNT gey;
    IPNT s; IPNT e;
    ra_a_gse(&rax,gsx,gex);
    ra_a_gse(&ray,gsy,gey);
    //        cerr << "\n===========================\n";
    //        cerr << "\n gsx[0]=" << gsx[0] << endl;
    //        cerr << "\n gex[0]=" << gex[0] << endl;
    //        cerr << "\n gsx[1]=" << gsx[1] << endl;
    //        cerr << "\n gex[1]=" << gex[1] << endl;
    //        cerr << "\n===========================\n";
    //        cerr << "\n gsy[0]=" << gsy[0] << endl;
    //        cerr << "\n gey[0]=" << gey[0] << endl;
    //        cerr << "\n gsy[1]=" << gsy[1] << endl;
    //        cerr << "\n gey[1]=" << gey[1] << endl;
    //        cerr << "\n===========================\n";
    // calculate grid overlap
    for (int ii=0;ii<dimx;ii++)  {
    s[ii]=max(gsy[ii],gsx[ii]);
    e[ii]=min(gey[ii],gex[ii]);
    }
        
    f2c_n1 = n_arr[0];
    f2c_n2 = n_arr[1];
    lendom=f2c_n1*f2c_n2;
    float _scale1=scale1;
    float _scale2=scale2;
    float _power=power;
    float _datum=datum;
    integer iter=0;
        
    // initialize workspace
    lenwork = 6*n_arr[1]*n_arr[0]+3*iwave_max(n_arr[1],2*n_arr[0])+21;
    //        cerr << "\n lenwork=" << lenwork << endl;
    //        cerr << "\n length of data = " << get_datasize_grid(gdom) << endl;
    //        cerr << "\n n_arr[0] = " << n_arr[0] << endl;
    //        cerr << "\n n_arr[1] = " << n_arr[1] << endl;
    //cerr << "\n physical domain size=" << lendom << endl;
    //cerr << "\n retrieveGlobalRank()=" << retrieveGlobalRank() << endl;        
    if (!(work = (float *)malloc(lenwork*sizeof(float)))) {
    RVLException e;
    e<<"Error: HelmOp::apply - failed to allocate " << lenwork << " floats for work buffer\n";
    throw e;
    }
    // allocate data arrays
    if (!(indata = (float *)malloc(lendom*sizeof(float)))) {
    RVLException e;
    e<<"Error: HelmOp::apply - failed to allocate " << lendom << " floats for input data\n";
    throw e;
    }
    if (!(outdata = (float *)malloc(lendom*sizeof(float)))) {
    RVLException e;
    e<<"Error: HelmOp::apply - failed to allocate " << lendom << " floats for output data\n";
    throw e;
    }
    IPNT i;
    integer idx;
    #if RARR_MAX_NDIM > 0
    if (dimx==1) {
    for (i[0]=s[0];i[0]<=e[0];i[0]++) {
    indata[i[0]-s[0]]=ray._s1[i[0]];
    }
    helm_(DirichletSides,&f2c_n1,&f2c_n2,
    &(d_arr[0]),&(d_arr[1]),
    &(_scale1),&(_scale2),
    &_power,&_datum,
    indata,
    outdata,
    work,
    &lenwork,
    &iter);
    fprintf(stderr, "\n indata [100] = %f\n", indata[100]);
    fprintf(stderr, "\n outdata [100] = %f\n", outdata[100]);
    for (i[0]=s[0];i[0]<=e[0];i[0]++) {
    rax._s1[i[0]]=outdata[i[0]-s[0]];
    }
    }
    #endif
    #if RARR_MAX_NDIM > 1
    if (dimx==2) {
    for (i[1]=s[1];i[1]<=e[1];i[1]++) {
    #pragma ivdep
    for (i[0]=s[0];i[0]<=e[0];i[0]++) {
    idx = (i[1]-s[1])*n_arr[0] + i[0]-s[0];
    indata[idx]=ray._s2[i[1]][i[0]];
    }
    }
    helm_(DirichletSides,&f2c_n1,&f2c_n2,
    &(d_arr[0]),&(d_arr[1]),
    &(_scale1),&(_scale2),
    &_power,&_datum,
    indata,
    outdata,
    work,
    &lenwork,
    &iter);
    fprintf(stderr, "\n indata [100] = %f\n", indata[100]);
    fprintf(stderr, "\n outdata [100] = %f\n", outdata[100]);
    // copy data back
    for (i[1]=s[1];i[1]<=e[1];i[1]++) {
    #pragma ivdep
    for (i[0]=s[0];i[0]<=e[0];i[0]++) {
    idx = (i[1]-s[1])*n_arr[0] + i[0]-s[0];
    rax._s2[i[1]][i[0]]=outdata[idx];
    }
    }
    }
    #endif
    #if RARR_MAX_NDIM > 2
    if (dimx==3) {
    //cerr << "\n dim3=" << e[2] << endl;
    for (i[2]=s[2];i[2]<=e[2];i[2]++) {
    for (i[1]=s[1];i[1]<=e[1];i[1]++) {
    #pragma ivdep
    for (i[0]=s[0];i[0]<=e[0];i[0]++) {
    idx = (i[1]-s[1])*n_arr[0] + i[0]-s[0];
    indata[idx]=ray._s3[i[2]][i[1]][i[0]];
    }
    }
    helm_(DirichletSides,&f2c_n1,&f2c_n2,
    &(d_arr[0]),&(d_arr[1]),
    &(_scale1),&(_scale2),
    &_power,&_datum,
    indata,
    outdata,
    work,
    &lenwork,
    &iter);
    // copy data back
    for (i[1]=s[1];i[1]<=e[1];i[1]++) {
    #pragma ivdep
    for (i[0]=s[0];i[0]<=e[0];i[0]++) {
    idx = (i[1]-s[1])*n_arr[0] + i[0]-s[0];
    rax._s3[i[2]][i[1]][i[0]]=outdata[idx];
    }
    }    
    }
    fprintf(stderr, "\n indata [100] = %f\n", indata[10]);
    fprintf(stderr, "\n outdata [100] = %f\n", outdata[10]);
    }
    #endif
    if (dimx<1 || dimx>3) {
    RVLException e;
    e<<"Error: HelmFO::operator()\n";
    e<<"dim = "<<dimx<<" outside of admissible set {1, 2}\n";
    throw e;
    }
    }
    catch (bad_cast) {
    RVLException e;
    e<<"\nError: HelmFO::operator()\n";
    e<<"at least one arg is not ContentPackage<float,RARR>\n";
    throw e;
    }
    catch (RVLException & e) {
    e<<"\ncalled from HelmFO::operator()\n";
    throw e;
    }
        
    }
    

    void GridHelmOp::apply(const Vector<float> & x,
    Vector<float> & y) const {
    try {
    // extract components - fine even if only one!
    Components<float> cx(x);
    Components<float> cy(y);

    // detect product structure
    ProductSpace<float> const * pdom = NULL;
    pdom = dynamic_cast<ProductSpace<float> const *>(&dom);
    int n_fac=1;
    if (pdom) n_fac=pdom->getSize();
    Space<float> const * sp = NULL;
   
    // component loop
    for (int j=0; j<n_fac; j++) {
    if (pdom) sp = &((*pdom)[j]);
    else sp = &dom;

    // late tests
    myGridSpace const * gdom = dynamic_cast<myGridSpace const *>(sp);
    if (!gdom) {
    RVLException e;
    e<<"Error: GridHelmOp::apply\n";
    e<<"  factor "<<j<<" of input space is not a GridSpace\n";
    e<<"  description:\n";
    sp->write(e);
    throw e;	  
    }
    if (retrieveGlobalRank() == 0) {
    if (gdom->getGrid().dim != 2) {
    RVLException e;
    e<<"Error: GridHelmOp::apply\n";
    e<<"  current implementation is 2D only\n";
    throw e;
    }
    }

    IPNT n_arr;
    RPNT d_arr;
    if (retrieveGlobalRank() == 0) {
    get_d(d_arr,gdom->getGrid());
    get_n(n_arr,gdom->getGrid());
    }
    HelmFO fo(n_arr,d_arr,weights[0],weights[1],power,datum,DirichletSides);
    MPISerialFunctionObject<float> mpifo(fo);
    cy[j].eval(mpifo,cx[j]);    
    }
    }
    catch (RVLException & e) {
    e<<"\ncalled in GridHelmOp::apply\n";
    throw e;
    }
            
    }
        
    void GridHelmOp::applyAdj(const Vector<float> & x,
    Vector<float> & y) const {
    try {
    apply(x,y);
    }
    catch (RVLException & e) {
    e<<"\ncalled in GridHelmOp::applyAdj\n";
    throw e;
    }
    }
  */

  
  void GridOverSprayFO::operator()(LocalDataContainer<float> & x,
				   LocalDataContainer<float> const & y) {
    try {

      ContentPackage< float, RARR > const & gy =
	dynamic_cast<ContentPackage< float, RARR > const &>(y);
            
      ContentPackage< float, RARR > & gx =
	dynamic_cast<ContentPackage< float, RARR > &>(x);
            
      // precondition - metadata are same dimn
      RARR & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();
      int dimx; int dimy;
      ra_ndim(&rax,&dimx);
      ra_ndim(&ray,&dimy);
      if (dimx != dimy) {
	RVLException e;
	e<<"Error: GridWeightOverlapFO::operator()\n";
	e<<"arguments have different dims:\n";
	e<<"dimx="<<dimx<<" dimy="<<dimy<<"\n";
	throw e;
      }
      // compute grid params
      IPNT gsx; IPNT gex;
      IPNT gsy; IPNT gey;
      IPNT s; IPNT e;
      ra_a_gse(&rax,gsx,gex);
      ra_a_gse(&ray,gsy,gey);
      // calculate grid overlap
      for (int ii=0;ii<dimx;ii++)  {
	s[ii]=max(gsy[ii],gsx[ii]);
	e[ii]=min(gey[ii],gex[ii]);
      }

      IPNT i;      
#if RARR_MAX_NDIM > 0
      if (dimx==1) {
#pragma ivdep
	for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	  rax._s1[i[0]]=ray._s1[i[0]];
	}
	if (axis==0) {
	  for (i[0]=gsx[0];i[0]<s[0];i[0]++) {
	    rax._s1[i[0]]=rax._s1[s[0]];
	  }
	  for (i[0]=e[0]+1;i[0]<=gex[0];i[0]++) {
	    rax._s1[i[0]]=rax._s1[e[0]];
	  }
	}
      }
#endif
#if RARR_MAX_NDIM > 1
      if (dimx==2) {
	for (i[1]=s[1];i[1]<=e[1];i[1]++) {
#pragma ivdep
	  for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	    rax._s2[i[1]][i[0]]=ray._s2[i[1]][i[0]];
	  }
	}
	if (axis==0) {
	  for (i[1]=s[1];i[1]<=e[1];i[1]++) {
#pragma ivdep
	    for (i[0]=gsx[0];i[0]<s[0];i[0]++) {
	      rax._s2[i[1]][i[0]]=rax._s2[i[1]][s[0]];
	    }
#pragma ivdep	    
	    for (i[0]=e[0]+1;i[0]<=gex[0];i[0]++) {
	      rax._s2[i[1]][i[0]]=rax._s2[i[1]][e[0]];
	    }
	  }
	}
	if (axis==1) {
	  for (i[1]=gsx[1];i[1]<s[1];i[1]++) {
#pragma ivdep	    
	    for (i[0]=gsx[0];i[0]<=gex[0];i[0]++) {
	      rax._s2[i[1]][i[0]]=rax._s2[s[1]][i[0]];
	    }
	  }
	  for (i[1]=e[1]+1;i[1]<=gex[1];i[1]++) {
#pragma ivdep	    
	    for (i[0]=gsx[0];i[0]<=gex[0];i[0]++) {
	      rax._s2[i[1]][i[0]]=rax._s2[e[1]][i[0]];
	    }
	  }	  
	}
      }
#endif
#if RARR_MAX_NDIM > 2
      if (dimx==3) {
	for (i[2]=s[2];i[2]<=e[2];i[2]++) {
	  for (i[1]=s[1];i[1]<=e[1];i[1]++) {
#pragma ivdep
	    for (i[0]=s[0];i[0]<=e[0];i[0]++) {
	      rax._s3[i[2]][i[1]][i[0]]=ray._s3[i[2]][i[1]][i[0]];
	    }
	  }
	}
	if (axis==0) {
	  for (i[2]=s[2];i[2]<=e[2];i[2]++) {
	    for (i[1]=s[1];i[1]<=e[1];i[1]++) {
#pragma ivdep
	      for (i[0]=gsx[0];i[0]<s[0];i[0]++) {
		rax._s3[i[2]][i[1]][i[0]]=rax._s3[i[2]][i[1]][s[0]];
	      }
#pragma ivdep	    
	      for (i[0]=e[0]+1;i[0]<=gex[0];i[0]++) {
		rax._s3[i[2]][i[1]][i[0]]=rax._s3[i[2]][i[1]][e[0]];
	      }
	    }
	  }
	}
	if (axis==1) {
	  for (i[2]=s[2];i[2]<=e[2];i[2]++) {
	    for (i[1]=gsx[1];i[1]<s[1];i[1]++) {
#pragma ivdep	    
	      for (i[0]=gsx[0];i[0]<=gex[0];i[0]++) {
		rax._s3[i[2]][i[1]][i[0]]=rax._s3[i[2]][s[1]][i[0]];
	      }
	    }
	    for (i[1]=e[1]+1;i[1]<=gex[1];i[1]++) {
#pragma ivdep	    
	      for (i[0]=gsx[0];i[0]<=gex[0];i[0]++) {
		rax._s3[i[2]][i[1]][i[0]]=rax._s3[i[2]][e[1]][s[0]];
	      }
	    }	  
	  }
	}
	if (axis==2) {
	  for (i[2]=gsx[2];i[2]<s[2];i[2]++) {
	    for (i[1]=gsx[1];i[1]<=gex[1];i[1]++) {
#pragma ivdep	    
	      for (i[0]=gsx[0];i[0]<=gex[0];i[0]++) {
		rax._s3[i[2]][i[1]][i[0]]=rax._s3[s[2]][i[1]][i[0]];
	      }
	    }
	  }
	  for (i[2]=e[2]+1;i[2]<=gex[2];i[2]++) {	  
	    for (i[1]=gsx[1];i[1]<=gex[1];i[1]++) {
#pragma ivdep	    
	      for (i[0]=gsx[0];i[0]<=gex[0];i[0]++) {
		rax._s3[i[2]][i[1]][i[0]]=rax._s3[e[2]][i[1]][s[0]];
	      }
	    }	  
	  }
	}
      }
      
#endif
      if (dimx<1 || dimx>3) {
	RVLException e;
	e<<"Error: GridOverSprayFO::operator()\n";
	e<<"dim = "<<dimx<<" outside of admissible set {1, 2, 3}\n";
	throw e;
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"\nError: GridOverSprayFO::operator()\n";
      e<<"at least one arg is not ContentPackage<float,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridOverSprayFO::operator()\n";
      throw e;
    }
  }
        
  void GridScaleFO::operator()(LocalDataContainer<float> & x,
			       LocalDataContainer<float> const & y) {
    
    try {
#ifdef IWAVE_VERBOSE
      cerr<<"in GridScaleFO\n";
#endif

      size_t nx = x.getSize();
      size_t ny = y.getSize();
      
#ifdef IWAVE_VERBOSE
      cerr<<"nx="<<nx<<" ny="
	  <<ny<<endl;
#endif	
      int nh = nx/ny;	
      if (nx != nh*ny) {
	RVLException e;
	e<<"Error: GridScaleFO\n";
	e<<"  length of target not integer multiple of length source\n";
	throw e;
      }      
      
#ifdef IWAVE_VERBOSE	
      cerr<<" nh="<<nh<<endl;
#endif
      float * sc = (float *) usermalloc_(ny*sizeof(float));
      for (int j=0; j<ny; j++) {
	sc[j]=fac*powf((y.getData())[j],pwr);
      }
      size_t loc = 0;
      for (int i=0; i<nh; i++) {
	for (int j=0; j<ny; j++) {
	  (x.getData())[loc+j] *= sc[j];
	}
	loc += ny;
      }
      userfree_(sc);
    }
    
    catch(RVL::RVLException & e) {
      e<<"\ncalled from GridScaleFO::operator()\n";
      throw e;
    }
  }

  // fac is data members
  void GridFwdZDerivFO::operator()(LocalDataContainer<float> & x,
				LocalDataContainer<float> const & y) {
    
    try {
#ifdef IWAVE_VERBOSE
      cerr<<"in GridFwdZDerivFO\n";
#endif
      ContentPackage< float, RARR > const & gy =
	dynamic_cast<ContentPackage< float, RARR > const &>(y);
            
      ContentPackage< float, RARR > & gx =
	dynamic_cast<ContentPackage< float, RARR > &>(x);
            
      // precondition - metadata are same dimn
      RARR & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();
      int dimx; int dimy;
      ra_ndim(&rax,&dimx);
      ra_ndim(&ray,&dimy);
      if (dimx != dimy) {
	RVLException e;
	e<<"Error: GridFwdZDerivFO::operator()\n";
	e<<"arguments have different dims:\n";
	e<<"dimx="<<dimx<<" dimy="<<dimy<<"\n";
	throw e;
      }
      // compute grid params
      IPNT nxa; IPNT nya;
      ra_a_size(&rax,nxa);
      ra_a_size(&ray,nya);
      for (int i=0;i<dimx;i++) {
	if (nxa[i] != nya[i]) {
	  RVLException e;
	  e<<"Error GridZderivFO::operator()";
	  e<<"  axis "<<i<<" lengths of input, output differ\n";
	  throw e;
	}
      }
      
      size_t nz = nxa[0];
      size_t ny = 1;
      
      for (int i=1;i<dimx;i++) ny *= nxa[i]; 
      
#ifdef IWAVE_VERBOSE
      size_t nx = x.getSize();
      cerr<<"nx="<<nx<<" ny="
	  <<ny<<" nz="<<nz<<endl;
#endif
      
      size_t loc = 0;
      float fac2=0.5f*fac;
      for (int i=0; i<ny; i++) {
	x.getData()[loc]=fac*(y.getData()[loc+1]-y.getData()[loc]);
	x.getData()[loc+nz-1]=fac*
	  (y.getData()[loc+nz-1]-y.getData()[loc+nz-2]);
#pragma ivdep
	for (int j=1; j<nz-1; j++) {
	  (x.getData())[loc+j] =fac2*
	    (y.getData()[loc+j+1]-y.getData()[loc+j-1]);
	}
	loc += nz;
      }
    }
    
    catch(RVL::RVLException & e) {
      e<<"\ncalled from GridFwdZDerivFO::operator()\n";
      throw e;
    }
  }

  GridZDerivOp::GridZDerivOp(Space<float> const & _dom)
    : dom(_dom), fac(1.0f) {
    try {
      myGridSpace const & gdom = dynamic_cast<myGridSpace const &>(dom);
      if (retrieveGlobalRank()==0) {
	RPNT d;
	get_d(d,gdom.getGrid());
	fac /= d[0];
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridZDerivOp constructor\n";
      e<<"  input space is not a GridSpace\n";
      e<<"  description:\n";
      dom.write(e);
      throw e;
    }
    catch (RVLException e) {
      e<<"\ncalled from GridZDerivOp constructor\n";
      throw e;
    }
  }

  GridZDerivOp::GridZDerivOp(GridZDerivOp const & op)
    : fac(op.fac), dom(op.dom) {}

  void GridZDerivOp::apply(Vector<float> const & x,
			   Vector<float> & y) const {
    try {
      GridFwdZDerivFO f(fac);
      MPISerialFunctionObject<float> mpif(f);
      y.eval(mpif,x);

#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridZDerivOp::apply\n";
      throw e;
    }
  }

  void GridZDerivOp::applyAdj(Vector<float> const & x,
			     Vector<float> & y) const {
    try {
      RVLException e;
      e<<"Error GridZDerivOp::applyAdj\n";
      e<<"  not implemented yet\n";
#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridZDerivOp::applyAdj\n";
      throw e;
    }
  }

  ostream & GridZDerivOp::write(ostream & str) const {
    str<<"GridZDerivOp: directional derivative, axis = 0\n";
    str<<"Domain:\n";
    dom.write(str);
    return str;
  }
  
  void GridZTaperFO::operator()(LocalDataContainer<float> & x) {
    
    try {
#ifdef IWAVE_VERBOSE
      cerr<<"in GridZTaperFO\n";
#endif
            
      ContentPackage< float, RARR > & gx =
	dynamic_cast<ContentPackage< float, RARR > &>(x);
            
      // precondition - metadata are same dimn
      RARR const & rax = gx.getMetadata();
      int dimx; 
      ra_ndim(&rax,&dimx);

      // compute grid params
      IPNT nxa;
      ra_a_size(&rax,nxa);
      
      size_t nz = nxa[0];
      size_t ny = 1;
      
      for (int i=1;i<dimx;i++) ny *= nxa[i]; 
      
#ifdef IWAVE_VERBOSE
      size_t nx = x.getSize();
      cerr<<"nx="<<nx<<" ny="
	  <<ny<<" nz="<<nz<<endl;
#endif
      
      size_t loc = 0;
      for (int i=0; i<ny; i++) {
#pragma ivdep
	for (int j=0; j<itop; j++) (x.getData())[loc+j] = 0.0f;
	for (int j=itop; j<ibot; j++) (x.getData())[loc+j] *= ((float)(j-itop+1)/(float)(ibot-itop+1));
	loc += nz;
      }
    }
    
    catch(RVL::RVLException & e) {
      e<<"\ncalled from GridZTaperFO::operator()\n";
      throw e;
    }
  }  

  GridZTaperOp::GridZTaperOp(Space<float> const & _dom,
			     float top, float bot)
    : dom(_dom) {
    try {
      RPNT d;
      RVL::ProductSpace<float> const * pdom = NULL;
      if (pdom = dynamic_cast<RVL::ProductSpace<float> const *>(&dom)) {;
	myGridSpace const & gdom = dynamic_cast<myGridSpace const &>((*pdom)[0]);
	if (retrieveGlobalRank()==0) {
	  get_d(d,gdom.getGrid());
	}
      }
      else {
	myGridSpace const & gdom = dynamic_cast<myGridSpace const &>(dom);
	if (retrieveGlobalRank()==0) {
	  get_d(d,gdom.getGrid());
	}
      }
      itop = int(0.1 +(top/d[0]));
      ibot = int(0.1 +(bot/d[0]));	
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridZTaperOp constructor\n";
      e<<"  input space is not a GridSpace\n";
      e<<"  description:\n";
      dom.write(e);
      throw e;
    }
    catch (RVLException e) {
      e<<"\ncalled from GridZTaperOp constructor\n";
      throw e;
    }
  }

  GridZTaperOp::GridZTaperOp(GridZTaperOp const & op)
    : dom(op.dom), itop(op.itop), ibot(op.ibot) {}

  void GridZTaperOp::apply(Vector<float> const & x,
			   Vector<float> & y) const {
    try {
      y.copy(x);
      GridZTaperFO f(itop,ibot);
      MPISerialFunctionObject<float> mpif(f);
      y.eval(mpif);

#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridZTaperOp::apply\n";
      throw e;
    }
  }

  void GridZTaperOp::applyAdj(Vector<float> const & x,
			     Vector<float> & y) const {
    try {
      this->apply(x,y);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridZTaperOp::applyAdj\n";
      throw e;
    }
  }

  ostream & GridZTaperOp::write(ostream & str) const {
    str<<"GridZTaperOp: taper on 0 axis\n";
    str<<"Domain:\n";
    dom.write(str);
    str<<"top index = "<<itop<<" bottom index = "<<ibot<<"\n";
    return str;
  }
  
  void GridGardnerFO::operator()(LocalDataContainer<float> & x,
				 LocalDataContainer<float> const & y) {
    try {
      if (x.getSize() != y.getSize()) {
	RVL::RVLException e;
	e<<"Error: GridGardnerFO::operator()\n";
	e<<"  input different size from output\n";
	e<<"  input:\n";
	y.write(e);
	e<<"  output:\n";
	x.write(e);
	throw e;
      }
      for (size_t i = 0; i < x.getSize(); i++) 
	x.getData()[i] = ((y.getData()[i] < 1.6) ? 1.0 : 0.31*powf(fabs(1000.0*y.getData()[i]),0.25));
    }
    catch (RVL::RVLException & e) {
      e<<"\ncalled from GridGardnerFO::operator()\n";
      throw e;
    }
  }


  /**************************** coord scale **********************/

  /* helper function */
  void next_idx(int dimx, IPNT n, IPNT nxa) {
    bool flag = true;
    for (int k=0; k < dimx; k++) {
      if (flag) {
	n[k]++;
	if (n[k] == nxa[k]) n[k]=0;
	else flag=false;
      }
    }
  }
  
  void GridCoordScaleFO::operator()(LocalDataContainer<float> & x) {

    try {
#ifdef IWAVE_VERBOSE
      cerr<<"in GridCoordScaleFO\n";
#endif
            
      ContentPackage< float, RARR > & gx =
	dynamic_cast<ContentPackage< float, RARR > &>(x);
            
      // precondition - metadata are same dimn
      RARR const & rax = gx.getMetadata();
      int dimx; 
      ra_ndim(&rax,&dimx);

      size_t len;
      ra_a_datasize(&rax,&len);

      // compute grid params
      IPNT nxa;
      ra_a_size(&rax,nxa);
      
      IPNT n;
      IASN(n,IPNT_0);
      
      for (size_t i=0;i<len;i++) {
	float fac=0;
	for (int k=0;k<dimx;k++) 
	  fac += (o[k] + n[k]*d[k])*a[k];
	fac += b;
	x.getData()[i] *= fac;
	next_idx(dimx,n,nxa);
      }

    }
    
    catch(RVL::RVLException & e) {
      e<<"\ncalled from GridCoordScaleFO::operator()\n";
      throw e;
    }
  }  

  GridCoordScaleOp::GridCoordScaleOp(Space<float> const & _dom,
				     std::vector<float> _a,
				     float _b)
    : dom(_dom), a(_a), b(_b) {
    try {
      RVL::ProductSpace<float> const * pdom = NULL;
      if (pdom = dynamic_cast<RVL::ProductSpace<float> const *>(&dom)) {;
	myGridSpace const & gdom = dynamic_cast<myGridSpace const &>((*pdom)[0]);
	if (retrieveGlobalRank()==0) {
	  get_d(d,gdom.getGrid());
	  get_o(o,gdom.getGrid());	  
	}
      }
      else {
	myGridSpace const & gdom = dynamic_cast<myGridSpace const &>(dom);
	if (retrieveGlobalRank()==0) {
	  get_d(d,gdom.getGrid());
	  get_o(o,gdom.getGrid());	  	  
	}
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridCoordScaleOp constructor\n";
      e<<"  input space is not a GridSpace\n";
      e<<"  description:\n";
      dom.write(e);
      throw e;
    }
    catch (RVLException e) {
      e<<"\ncalled from GridCoordScaleOp constructor\n";
      throw e;
    }
  }

  void GridCoordScaleOp::apply(Vector<float> const & x,
			   Vector<float> & y) const {
    try {
      y.copy(x);
      GridCoordScaleFO f(a,b,o,d);
      MPISerialFunctionObject<float> mpif(f);
      y.eval(mpif);

#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridCoordScaleOp::apply\n";
      throw e;
    }
  }

  void GridCoordScaleOp::applyAdj(Vector<float> const & x,
			     Vector<float> & y) const {
    try {
      this->apply(x,y);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridCoordScaleOp::applyAdj\n";
      throw e;
    }
  }

  ostream & GridCoordScaleOp::write(ostream & str) const {
    str<<"GridCoordScaleOp:\n";
    str<<"Domain:\n";
    dom.write(str);
    return str;
  }

  
  void GridRadialScaleFO::operator()(LocalDataContainer<float> & x) {

    try {
#ifdef IWAVE_VERBOSE
      cerr<<"in GridRadialScaleFO\n";
#endif
            
      ContentPackage< float, RARR > & gx =
	dynamic_cast<ContentPackage< float, RARR > &>(x);
            
      // precondition - metadata are same dimn
      RARR const & rax = gx.getMetadata();
      int dimx; 
      ra_ndim(&rax,&dimx);

      size_t len;
      ra_a_datasize(&rax,&len);

      // compute grid params
      IPNT nxa;
      ra_a_size(&rax,nxa);
      
      IPNT n;
      IASN(n,IPNT_0);
      
      for (size_t i=0;i<len;i++) {
	float fac=0;
	for (int k=0;k<dimx;k++) 
	  fac += (o[k] + n[k]*d[k])*(o[k] + n[k]*d[k]);
	fac = alpha * sqrt(fac);
	x.getData()[i] *= fac;
	next_idx(dimx,n,nxa);
      }

    }
    
    catch(RVL::RVLException & e) {
      e<<"\ncalled from GridRadialScaleFO::operator()\n";
      throw e;
    }
  }  

  GridRadialScaleOp::GridRadialScaleOp(Space<float> const & _dom,
				       float _alpha)
    : dom(_dom), alpha(_alpha) {
    try {
      RVL::ProductSpace<float> const * pdom = NULL;
      if (pdom = dynamic_cast<RVL::ProductSpace<float> const *>(&dom)) {;
	myGridSpace const & gdom = dynamic_cast<myGridSpace const &>((*pdom)[0]);
	if (retrieveGlobalRank()==0) {
	  get_d(d,gdom.getGrid());
	  get_o(o,gdom.getGrid());	  
	}
      }
      else {
	myGridSpace const & gdom = dynamic_cast<myGridSpace const &>(dom);
	if (retrieveGlobalRank()==0) {
	  get_d(d,gdom.getGrid());
	  get_o(o,gdom.getGrid());	  	  
	}
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: GridRadialScaleOp constructor\n";
      e<<"  input space is not a GridSpace\n";
      e<<"  description:\n";
      dom.write(e);
      throw e;
    }
    catch (RVLException e) {
      e<<"\ncalled from GridRadialScaleOp constructor\n";
      throw e;
    }
  }

  void GridRadialScaleOp::apply(Vector<float> const & x,
			   Vector<float> & y) const {
    try {
      y.copy(x);
      GridRadialScaleFO f(alpha,o,d);
      MPISerialFunctionObject<float> mpif(f);
      y.eval(mpif);

#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridRadialScaleOp::apply\n";
      throw e;
    }
  }

  void GridRadialScaleOp::applyAdj(Vector<float> const & x,
			     Vector<float> & y) const {
    try {
      this->apply(x,y);
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridRadialScaleOp::applyAdj\n";
      throw e;
    }
  }

  ostream & GridRadialScaleOp::write(ostream & str) const {
    str<<"GridRadialScaleOp:\n";
    str<<"Domain:\n";
    dom.write(str);
    return str;
  }

  // moving avg - stride 1
  void ma1(const float * restrict x, float * restrict y, int n, int rad) {
    // start - fill in half-radii on ends
    // example: rad=1
    // y[0] = x[0] + x[1] + x[2];
    // y[n-1] = x[n-1] + x[n-2] + x[n-3];
    // y[1] = x[0] + x[1] + x[2] + x[3];
    // y[n-2] = x[n-1] + x[n-2] + x[n-3] + x[n-4];
    y[0]=x[0];
    y[n-1]=x[n-1];
    for (int i=1; i<rad+1; i++) {
      y[0] += x[i]; y[0] += x[n-1-i];
      y[n-1] += x[i-1]; y[n-1] += x[n-1-i];
    }
    for (int i=1; i<rad+1; i++) {
      y[i] = y[i-1] + x[i+rad] - x[i-rad-1+n];
      y[n-1-i] = y[n-i] + x[n-1-i-rad] - x[rad-i];
    }
    // expl: have y[0], y[1], y[2], y[n-3], y[n-2], y[n-1]
    // middle
    // expl: y[3] = y[2] + x[5] - x[0]
    //       y[n-4] = y[n-5] + x[n-2] - x[n-7]
#pragma ivdep
#pragma GCC ivdep      
    for (int i=rad; i<n-rad-1; i++) {
      y[i] = y[i-1] + x[i+rad] - x[i-rad-1];
    }
    // scale
#pragma ivdep
#pragma GCC ivdep      
    for (int i=0; i<n; i++) y[i] /= (2.0*rad+1.0);
  }

  // moving avg - arbitrary origin, stride
  void mas(const float * restrict x, float * restrict y, int n, int o, int d, int rad) {
    if (d==1) {
      ma1(x+o, y+o, n, rad);
    }
    else {

      y[o]=x[o];
      y[o+(n-1)*d]=x[o+(n-1)];

      for (int i=1; i<rad+1; i++) {
	y[o] += x[o+i*d]; y[o] += x[o+(n-1-i)*d];
	y[o+(n-1)*d] += x[o+(i-1)*d]; y[o+(n-1)*d] += x[o+(n-1-i)*d];
      }

      for (int i=1; i<rad+1; i++) {
	y[o+i*d] = y[o+(i-1)*d] + x[o+(i+rad)*d] - x[o+(i-rad-1+n)*d];
	y[o+(n-1-i)*d] = y[o+(n-i)*d] + x[o+(n-1-i-rad)*d] - x[o+(rad-i)*d];
      }

#pragma ivdep
#pragma GCC ivdep      
      for (int i=rad+1; i<n-rad-1; i++) {
	y[o+i*d] = y[o+(i-1)*d] + x[o+(i+rad)*d] - x[o+(i-rad-1)*d];
      }

#pragma ivdep
#pragma GCC ivdep
      for (int i=0; i<n; i++) y[o+i*d] /= (2.0*rad+1.0);
    }
  }
  
  void ma(const float * restrict x, float * restrict y, int dim, int * n, int * rad) {
    
    if (dim==1) {
      ma1(x,y,n[0],rad[0]);
    }
    else {

      size_t len=1;
      for (int i=0; i< dim; i++) len *= n[i];
      float * buf = (float *)usermalloc_(len*sizeof(float));
      
      if (dim==2) {

	for (int i1=0; i1<n[1]; i1++) {
	  mas(x, y, n[0], i1*n[0], 1, rad[0]);
	}
	
	memcpy(buf,y,len*sizeof(float));

	for (int i0=0; i0<n[0]; i0++) {
	  mas(buf, y, n[1], i0, n[0], rad[1]);
	}

	userfree_(buf);
      }
    
      else if (dim==3) {
	// average over dim 0
	for (int i2=0; i2<n[2]; i2++) {
	  for (int i1=0; i1<n[1]; i1++) {
	    mas(x, y, n[0], i1*n[0]+i2*n[0]*n[1], 1, rad[0]);
	  }
	}

	memcpy(buf,y,len*sizeof(float));
	
	// average over dim 1
	for (int i2=0; i2<n[2]; i2++) {      
	  for (int i0=0; i0<n[0]; i0++) {
	    mas(buf, y, n[1], i0+i2*n[0]*n[1], n[0], rad[1]);
	  }
	}

	memcpy(buf,y,len*sizeof(float));
	
	// average over dim 2
	for (int i2=0; i2<n[2]; i2++) {
	  for (int i1=0; i1<n[1]; i1++) {
	    mas(buf, y, n[2], i1*n[0]+i2*n[0]*n[1], n[0]*n[1], rad[2]);
	  }
	}

	userfree_(buf);
      }
      else {
	userfree_(buf);
	RVLException e;
	e<<"Error: ma [moving average, C function]\n";
	e<<"  dim="<<dim<<" not legit (1, 2, or 3)\n";
	throw e;
      }
    }
  }
  
  void GridMAFO::operator()(LocalDataContainer<float> & x,
			    LocalDataContainer<float> const & y) {
    try {

      // cerr<<"GridMAFO::operator() begin\n";
      ContentPackage< float, RARR > const & gy =
	dynamic_cast<ContentPackage< float, RARR > const &>(y);
      
      ContentPackage< float, RARR > & gx =
	dynamic_cast<ContentPackage< float, RARR > &>(x);
      
      // precondition - metadata are same dimn
      RARR & rax = gx.getMetadata();
      RARR const & ray = gy.getMetadata();
      int dimx; int dimy;
      ra_ndim(&rax,&dimx);
      ra_ndim(&ray,&dimy);
      if (dimx != dimy) {
	RVLException e;
	e<<"Error: GridMAFO::operator()\n";
	e<<"arguments have different dims:\n";
	e<<"dimx="<<dimx<<" dimy="<<dimy<<"\n";
	throw e;
      }

      for (int i=0; i<dimx; i++) {
	if (rad[i]<0) {
	  RVL::RVLException e;
	  e<<"Error: GridMAFO::operator()\n";
	  e<<"  averaging radius for axis "<<i<<" = "<<rad[i]
	   <<" must be non-negative int\n"; 
	  throw e;
	}
      }

      // compute grid params
      IPNT nx;
      IPNT ny;
      ra_a_size(&rax,nx);
      ra_a_size(&ray,ny);
      
      // check same lattice geom
      for (int i=0; i<dimx; i++) {
	if (nx[i] != ny[i]) {
	  RVLException e;
	  e<<"Error: GridMAFO::operator()\n";
	  e<<"arguments have different sizes in dim "<<i<<":\n";
	  e<<"nx[dim]="<<nx[i]<<" ny[dim]="<<ny[i]<<"\n";
	  throw e;
	}
      }

      // call raw function
      ma(y.getData(), x.getData(), dimx, nx, rad);

      // to repeat, need a buffer, since y is const
      if (rep > 1) {
	size_t len;
	ra_a_datasize(&rax,&len);
	float * buf = (float *)usermalloc_(len*sizeof(float));
	for (int i=1; i<rep; i++) {
	  memcpy(buf,x.getData(),len*sizeof(float));
	  ma(buf, x.getData(), dimx, nx, rad);
	}
	userfree_(buf);
      }

    }
    catch (bad_cast) {
      RVL::RVLException e;
      e<<"Error: GridMAFO::operator()\n";
      e<<"  Failed to cast LDC to GridLDC\n";
      throw e;
    }
    catch (RVL::RVLException & e) {
      e<<"\ncalled from GridMAFO::operator()\n";
      throw e;
    }
  }

  std::shared_ptr<LinearOp<float> >
  GridMAOpBuilder::build(std::shared_ptr<GridOpData> d) const {

    // declare lopM - null by default
    std::shared_ptr<RVL::LinearOp<float> > lopM;
    
    // if f is not initialized, initialize it, unless
    // MA defined by param choices is trivial
    if (!f) {
      
      IPNT rad;
      IASN(rad,IPNT_0);
      int rmax=0;
      for (int i=1; i < RARR_MAX_NDIM+1; i++) {
	std::stringstream foo;
	foo<<i;
	std::string barf = "rect";
	rad[i-1] = RVL::valparse<int>(*(d->pars),barf+foo.str(),0);
	if (rad[i-1]<0) {
	  RVL::RVLException e;
	  e<<"Error: GridMAOpBuilder\n";
	  e<<"  input rect"<<i<<" < 0\n";
	  throw e;
	}
	rmax += rad[i-1];
      }
         
      int rep = RVL::valparse<int>(*(d->pars),"repeat",1);
      if (rep < 1) {
	RVL::RVLException e;
	e<<"Error: GridMAOpBuilder\n";
	e<<"  input repeat"<<" < 1\n";
	throw e;
      }

      // if params define nontrivial MA and f not already
      // constructed, construct it
      if (rmax >0) f = make_shared<TSOpt::GridMAFO>(rep,rad);

    }

    // if f points to FO, build LOp - note that op owns dom,
    // this owns f, so no dangles
    
    if (f) {
      lopM = make_shared<RVL::LinearOpFO<float> >(d->dom,
						  d->dom,
						  *f, *f);
    }
    
    return lopM;

  }
  
  std::ostream & GridMAOpBuilder::write(std::ostream & str) const {
    if (f) {
      str<<"Moving Average builder: returns LinearOp based on FunctionObject\n";
      f->write(str);
    }
    else {
      str<<"Moving average builder: trivial or uninitialized, returns "
	 <<"null pointer\n";
    }
    return str;
  }

}
