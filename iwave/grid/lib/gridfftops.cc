#include "gridfftops.hh"

#undef IWAVE_VERBOSE

namespace TSOpt {

  void GridZFTScaleFO::operator()(RVL::LocalDataContainer<float> & x) {
    
    try {
#ifdef IWAVE_VERBOSE
      cerr<<"in GridZFTScaleFO\n";
#endif

      if (ab==0 && adj==1) {
	RVL::RVLException e;
	e<<"Error: GridZFTScaleFO::operator()\n";
	e<<"  adjoint not implemented yet for non-abs case\n";
	throw e;
      }
      
      RVL::ContentPackage< float, RARR > & gx =
	dynamic_cast<RVL::ContentPackage< float, RARR > &>(x);
            
      RARR & rax = gx.getMetadata();
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
      cerr<<"GridZFTScaleFO::operator() nx="<<nx<<" ny="
	  <<ny<<" nz="<<nz<<endl;
#endif

      size_t loc = 0;
      for (int i=0; i<ny; i++) {
	FTScale(&(x.getData()[loc]),nz,dz,ab,pow,band);
	loc += nz;
      }
    }
    
    catch(RVL::RVLException & e) {
      e<<"\ncalled from GridZFTScaleFO::operator()\n";
      throw e;
    }
  }
  
  GridZFTScaleOp::GridZFTScaleOp(RVL::Space<float> const & _dom,
				 float _pow,
				 int _ab,
				 float locut,
				 float lopas,
				 float hipas,
				 float hicut)
    : dom(_dom), pow(_pow), ab(_ab), band(new float[4]),
      dz(1.0f) {
    try {
      myGridSpace const & gdom = dynamic_cast<myGridSpace const &>(dom);
      band[0]=locut;
      band[1]=lopas;
      band[2]=hipas;
      band[3]=hicut;
      if (retrieveGlobalRank()==0) {
	RPNT d;
	get_d(d,gdom.getGrid());
	dz = d[0];
      }
    }
    catch (bad_cast) {
      RVL::RVLException e;
      e<<"Error: GridZFTScaleOp constructor\n";
      e<<"  input space is not a GridSpace\n";
      e<<"  description:\n";
      dom.write(e);
      throw e;
    }
    catch (RVL::RVLException e) {
      e<<"\ncalled from GridZFTScaleOp constructor\n";
      throw e;
    }
  }
  GridZFTScaleOp::GridZFTScaleOp(GridZFTScaleOp const & op)
    : dom(op.dom), pow(op.pow), ab(op.ab), band(new float[4]), dz(op.dz) {
    band[0]=op.band[0];
    band[1]=op.band[1];
    band[2]=op.band[2];
    band[3]=op.band[3];
  }

  void GridZFTScaleOp::apply(RVL::Vector<float> const & x,
			     RVL::Vector<float> & y) const {
    try {
      GridZFTScaleFO f(dz,pow,ab,band,0);
      RVL::MPISerialFunctionObject<float> mpif(f);
      y.copy(x);
      y.eval(mpif);

#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridZFTScaleOp::apply\n";
      throw e;
    }
  }

  void GridZFTScaleOp::applyAdj(RVL::Vector<float> const & x,
				RVL::Vector<float> & y) const {
    try {
      GridZFTScaleFO f(dz,pow,ab,band,1);
      RVL::MPISerialFunctionObject<float> mpif(f);
      y.copy(x);
      y.eval(mpif);

#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from GridZFTScaleOp::apply\n";
      throw e;
    }
  }

  ostream & GridZFTScaleOp::write(ostream & str) const {
    str<<"GridZFTScaleOp: directional derivative, axis = 0\n";
    str<<"Domain:\n";
    dom.write(str);
    return str;
  }

}
