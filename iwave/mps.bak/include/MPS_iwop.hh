// MPS_iwop.hh
// Author: Mario J. Bencomo
// last modified: 12/12/16

/**
 * \file MPS_iwop.hh
 * Implementation of an MPS version of IWaveLOVOp.
 *
 * Defines MPS_IWaveLOVOp via multichannel convolution with proper 
 * Green's functions computed via operator composition between 
 * MPS_to_RHS and IWaveLOVOp. These operators are initialized internally.
 *
 * The MPS_IWaveLOVOp constructor will need the following MPS_KEYS initialized
 * with their values supplied by user (see main constructor for more info):
 *   - MPS_KEYS.appx_ord
 *   - MPS_KEYS.MPS_file
 *   - MPS_KEYS.MPS_grid
 *   - MPS_KEYS.hdr_files
 *   - MPS_KEYS.delta_file
 *   - MPS_KEYS.G_files
 * Other MPS_KEYS do not need to be initialized.
 */

#ifndef __MPS_IWOP_HH_
#define __MPS_IWOP_HH_

#include "MPS_includes.hh"
#include "MPS_to_RHS.hh"
#include "MPS_conv.hh"

//#define VERBOSE_MPS_IWOP_HH

namespace TSOpt {

  using TSOpt::MPS_conv;
  using TSOpt::MPS_to_RHS;

  /**
   * This function zeroes out all traces of a given SEGY file, except trace idx. 
   * Useful for setting temporary internal deltas.
   *
   * @param [in] filename - filename of SEGY input file
   * @param [in] idx - index of trace to not be zeroed out
   * @param [in] Nsrc - number of sources
   */
  void zero_out_ex(string filename, int idx, int Nsrc);

  //==============================================================================//
  // MPS_IWaveLOVOp class
  //==============================================================================//
  /**
   * MPS version of IWaveLOVOp.
   * Forward map is computed as a multichannel convolution operator with Green's functions.
   * Green's functions are computed and stored during initializing. 
   * Recomputation of Green's functions occur whenever forward map is computed for different models.
   */
  template< class T_MPS_Space >
  class MPS_IWaveLOVOp: public LinOpValOp<ireal> {

  private:

    PARARRAY *pars;        /**< key-val pairs */
    IWaveInfo ic;          /**< IWaveInfo struct */
    mutable FILE *stream;  /**< stream */
    MPS_KEYS mks;          /**< MPS keys */

    StdProductSpace<ireal> dom; /**< domain space, product space of source and model parameters */
    IWaveSpace rng;             /**< data range space */
    T_MPS_Space MPS_dom;        /**< underlying MPS space */

    shared_ptr<T_MPS_Space> delta_sp_ptr;             /**< MPS spaces related to delta source parameter input*/
    vector< shared_ptr< Vector<ireal> > > delta_ptrs; /**< delta source paramter input*/

    StdProductSpace<ireal> int_x_sp;       /**< space of internal copy of model-source parameter */
    shared_ptr< Vector<ireal> > int_x_ptr; /**< internal copy of model-source parameter*/

    shared_ptr<IWaveSpace> G_sp_ptr;                       /**< space associated with Green's functions*/
    mutable vector< shared_ptr< Vector<ireal> > > G_ptrs;  /**< Green's functions */

    mutable shared_ptr<MPS_conv> conv_ptr;          /**< internal convolution operator */
    shared_ptr< LinCompLOVOp<ireal> > lovop_ptr;    /**< internal op from linear composition of IWaveLOVOp with MPS_to_RHS for output data*/
    shared_ptr< LinCompLOVOp<ireal> > G_lovop_ptr;  /**< internal op from linear composition of IWaveLOVOp with MPS_to_RHS for green's function*/

    //handy functions for setting data members
    void set_delta();
    void set_int_x();
    void set_G();
    void set_ops();

    /**
     * Handy function call for recomputing Green's functions.
     */
    void compute_G() const;
    //    void sanity_checks() const;

  protected:

    void apply0(const Vector<ireal> & x0,
		const Vector<ireal> & x1,
		Vector<ireal> & y) const;

    void applyAdj0(const Vector<ireal> & x0,
		   const Vector<ireal> & y,
		   Vector<ireal> & x1) const;

    void applyPartialDeriv0(const Vector<ireal> & x0,
			    const Vector<ireal> & x1,
			    const Vector<ireal> & dx0,
			    Vector<ireal> & dy) const;

    void applyAdjPartialDeriv0(const Vector<ireal> & x0,
			       const Vector<ireal> & x1,
			       const Vector<ireal> & dy,
			       Vector<ireal> & dx0) const;
      
    void applyPartialDeriv20(const Vector<ireal> & x0, 
			     const Vector<ireal> & x1,
			     const Vector<ireal> & dx00,
			     Vector<ireal> dy) const;

    void applyAdjPartialDeriv20(const Vector<ireal> & x0, 
				const Vector<ireal> & x1,
				const Vector<ireal> & dy,
				Vector<ireal> dx00) const;

    RVL::OperatorProductDomain<ireal> * clonePD() const{
      return new MPS_IWaveLOVOp(*this);
    }

  public:

    /**
     * Copy constructor.
     */
    MPS_IWaveLOVOp( MPS_IWaveLOVOp<T_MPS_Space> const &op );

    /**
     * Main constructor.
     *
     * MPS_IWaveLOVOp main constructor requires the following MPS_KEYS be initialized:
     *     - MPS_KEYS.appx_ord   - key for singular source approximation order
     *     - MPS_KEYS.MPS_file   - key for MPS filename
     *     - MPS_KEYS.grid_file  - key for RSF file containing spatial grid info
     *     - MPS_KEYS.RHS_files  - keys for RHS filenames required for MPS_to_RHS op
     *     - MPS_KEYS.hdr_files  - keys for data filenames 
     *     - MPS_KEYS.delta_file - key for SU file containing delta MPS 
     *     - MPS_KEYS.G_files    - keys for Green's function filenames
     *
     * @param [in] _mks - keys related to initialization of MPS spaces, see \ref MPS_KEYS
     * @param [in] _pars - input key-val pairs
     * @param [in] _stream
     */
    MPS_IWaveLOVOp( MPS_KEYS _mks,
		    PARARRAY const &_pars,
		    FILE *_stream);

    /**
     * Destructor.
     */
    ~MPS_IWaveLOVOp(){ ps_delete(&pars);}
    
    /**
     * Handy function for retrieving domain space as product space.
     */
    const ProductSpace<ireal> & getProductDomain() const{ return dom;}

    /**
     * Handy function for retrieving range space.
     */
    const Space<ireal> & getRange() const { return rng;}

    /**
     * Handy function for retrieving nonlinear part of domain space.
     */
    const IWaveSpace & getNonLinDomain() const;

    /**
     * Handy function for retrieving linear part of domain space.
     */
    const T_MPS_Space & getLinDomain() const;
    //PARARRAY & getPar(){}
    //PARARRAY const & getPar() const{}

    ostream &write(ostream &str) const;
      
  };


  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  ostream & MPS_IWaveLOVOp<T_MPS_Space>::write(ostream &str) const{
  //----------------------------------------------------------------------//
    str<<"MPS_IWaveLOVOp object\n";
    str<<"\n-Domain: \n";           dom.write(str); str<<"\n";
    str<<"\n-Range: \n";            rng.write(str); str<<"\n";
    str<<"\n-MPS domain: \n";       MPS_dom.write(str); str<<"\n";
    str<<"\n-delta vector: \n";     delta_ptrs[0]->write(str); str<<"\n";
    str<<"\n-internal x: \n";       int_x_ptr->write(str); str<<"\n";
    str<<"\n-Green's function: \n"; G_ptrs[0]->write(str); str<<"\n";
    str<<"\n-internal LOVOp: \n";   lovop_ptr->write(str); str<<"\n";
    str<<"\n-internal GLOVOp: \n";  G_lovop_ptr->write(str); str<<"\n";
    str<<"\n-convolution op: \n";  conv_ptr->write(str); str<<"\n";
    return str;
  }

  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  MPS_IWaveLOVOp<T_MPS_Space>::MPS_IWaveLOVOp( MPS_IWaveLOVOp<T_MPS_Space> const &op)
    : stream(op.stream),
      pars(NULL),
      mks(op.mks),
      rng(op.rng),
      dom(op.dom),
      MPS_dom(op.MPS_dom),
      delta_sp_ptr(op.delta_sp_ptr),
      delta_ptrs(op.delta_ptrs),
      int_x_sp(op.int_x_sp),
      int_x_ptr(op.int_x_ptr),
      G_sp_ptr(op.G_sp_ptr),
      G_ptrs(op.G_ptrs),
      conv_ptr(op.conv_ptr),
      lovop_ptr(op.lovop_ptr),
      G_lovop_ptr(op.G_lovop_ptr)
  {
  //----------------------------------------------------------------------//
    try{
      //copying pars
      pars = ps_new();
      ps_copy(&pars,*op.pars);
    }
    catch(RVLException &e){
      e << "ERROR from MPS_IWaveLOVOp copy constructor\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  MPS_IWaveLOVOp<T_MPS_Space>::MPS_IWaveLOVOp(MPS_KEYS _mks,
					      PARARRAY const &_pars,
					      FILE *_stream)
    : stream(_stream),
      pars(NULL),
      mks(_mks),
      rng(_pars,ic,false,2),
      dom(2),
      MPS_dom(_mks,_pars,false),
      delta_sp_ptr(),
      delta_ptrs(MPS_dom.get_dim()),
      int_x_sp(2),
      int_x_ptr(),
      G_sp_ptr(),
      G_ptrs(MPS_dom.get_dim()),
      conv_ptr(),
      lovop_ptr(),
      G_lovop_ptr()
  {
  //----------------------------------------------------------------------//
    try{

      //copying pars
      pars = ps_new();
      ps_copy(&pars,_pars);

      //setting domain
      IWaveSpace IW_dom(*pars,ic,true,1);
      dom.set(IW_dom,0);
      dom.set(MPS_dom,1);

      set_delta();
      set_int_x();
      set_G();
      set_ops();
      //sanity_checks();
    }
    catch(RVLException &e){
      e << "ERROR from MPS_IWaveLOVOp constructor\n";
      throw e;
    }
  }
  
  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  const IWaveSpace & MPS_IWaveLOVOp<T_MPS_Space>::getNonLinDomain() const{
  //----------------------------------------------------------------------//
    try{
      IWaveSpace const & iwsp = 
	dynamic_cast<IWaveSpace const&> (this->getProductDomain()[0]);
      return iwsp;
    }
    catch(RVLException & e){
      e << "ERROR from MPS_IWaveLOVOp::getNonLinDomain\n";
      throw e;
  }
    catch(bad_cast){
      RVLException e;
      e << "ERROR from MPS_IWaveLOVOp::getNonLinDomain\n"
      << " trouble type casting first component as IWaveSpace\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  const T_MPS_Space & MPS_IWaveLOVOp<T_MPS_Space>::getLinDomain() const{
  //----------------------------------------------------------------------//
    try{
      T_MPS_Space const & mps_sp = 
	dynamic_cast<T_MPS_Space const&> (this->getProductDomain()[1]);
      return mps_sp;
    }
    catch(RVLException & e){
      e << "ERROR from MPS_IWaveLOVOp::getLinDomain\n";
      throw e;
    }
    catch(bad_cast){
      RVLException e;
      e << "ERROR from MPS_IWaveLOVOp::getLinDomain\n"
	<< " trouble type casting second component as MPS_Space\n";
      throw e;
    }
  }
  
  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  void MPS_IWaveLOVOp<T_MPS_Space>::compute_G() const{
  //----------------------------------------------------------------------//
    try{
      
      for( int i=0; i<G_ptrs.size(); i++){
	Components<ireal> c_int_x(*int_x_ptr);
	c_int_x[1].copy(*delta_ptrs[i]);
	Operator<ireal>::export_apply(*G_lovop_ptr,
				      *int_x_ptr,
				      *G_ptrs[i]);
      }
    }
    catch(RVLException & e){
      e << "ERROR from MPS_IWaveLOVOp::compute_G\n";
      throw e;
    }
  }


  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  void MPS_IWaveLOVOp<T_MPS_Space>::set_delta() {
  //----------------------------------------------------------------------//
    try{
      
      PARARRAY *locpars;
      locpars = ps_new();
      ps_copy(&locpars,*pars);

      //setting MPS delta space
      string ref_file = valparse<string>(*pars,mks.delta_file);
      add_to_pars(*locpars,
		  mks.MPS_file,
		  ref_file);
      
      T_MPS_Space tmp_sp(mks,*locpars);
      shared_ptr< Space<ireal> > tmp_ptr = Space<ireal>::clonePtr(tmp_sp);
      delta_sp_ptr = dynamic_pointer_cast< T_MPS_Space >(tmp_ptr);

      //setting reference delta vector
      Vector<ireal> ref_delta(*delta_sp_ptr);
      AssignFilename af_ref(ref_file);
      ref_delta.eval(af_ref);

      //initializing deltas
      for( int i=0; i<delta_ptrs.size(); i++){
	delta_ptrs[i] = Vector<ireal>::newPtr(*delta_sp_ptr);
	delta_ptrs[i]->copy(ref_delta);
	
	//injecting filename into locpars under key mks.MPS_file
	AssignParams ap(*locpars,mks.MPS_file);
	delta_ptrs[i]->eval(ap);

#ifdef IWAVE_USE_MPI
	if(retrieveGlobalRank()==0){
#endif
	  zero_out_ex(valparse<string>(*locpars,mks.MPS_file),i,
		      MPS_dom.get_s_pos().size());
#ifdef IWAVE_USE_MPI
	}
	MPI_Barrier(retrieveGlobalComm());
#endif	
      }
      
      ps_delete(&locpars);
    }
    catch(RVLException &e){
      e << "ERROR from MPS_IWaveLOVOp::set_delta\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  void MPS_IWaveLOVOp<T_MPS_Space>::set_int_x() {
  //----------------------------------------------------------------------//
    try{
      //initializing int_sp
      int_x_sp.set(dom[0],0);
      int_x_sp.set(*delta_sp_ptr,1);

      //initializing int_x
      int_x_ptr = Vector<ireal>::newPtr(int_x_sp);
      int_x_ptr->zero();
    }
    catch(RVLException &e){
      e << "ERROR from MPS_IWaveLOVOp::set_int_x\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  void MPS_IWaveLOVOp<T_MPS_Space>::set_G() {
  //----------------------------------------------------------------------//
    try{

      PARARRAY *locpars;
      locpars = ps_new();
      ps_copy(&locpars,*pars);

      //setting IWaveSpace for G
      vector<string> G_files(mks.G_files.size());
      for( int i=0; i<G_files.size(); i++){
	G_files[i] = valparse<string>(*pars,mks.G_files[i]);
      }

      add_to_pars(*locpars,
		  mks.hdr_files,
		  G_files);
      
      IWaveSpace tmp_sp(*locpars,ic,false,2);
      shared_ptr< Space<ireal> > tmp_ptr = Space<ireal>::clonePtr(tmp_sp);
      G_sp_ptr = dynamic_pointer_cast< IWaveSpace >(tmp_ptr);

      //initializing G's
      for(int i=0; i<G_ptrs.size(); i++){
	G_ptrs[i] = Vector<ireal>::newPtr(*G_sp_ptr);
	G_ptrs[i]->zero();
      }
      ps_delete(&locpars);
    }
    catch(RVLException &e){
      e << "ERROR from MPS_IWaveLOVOp::set_G\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  void MPS_IWaveLOVOp<T_MPS_Space>::set_ops() {
  //----------------------------------------------------------------------//
    try{

      PARARRAY *locpars;
      locpars = ps_new();
      ps_copy(&locpars,*pars);
      
      ////////////////////////////
      // initializing lovop_ptr //
      ////////////////////////////
      //resetting RHS filenames
      for(int i=0; i<mks.RHS_files.size(); i++)
	add_to_pars(*locpars,mks.RHS_files[i],"empty");

      MPS_to_RHS m2r(mks,*locpars,MPS_dom);      

      //adding generated RHS files to locpars for IWaveOp initialization      
      add_to_pars(*locpars,
		  m2r.get_IWave_Range().getKeys(),
		  m2r.get_rng_filenames());

      IWaveLOVOp iwop(*locpars,stream);
      LinCompLOVOp<ireal> lovop(m2r,iwop);
      shared_ptr< Operator<ireal> > tmp_ptr = Operator<ireal>::clonePtr(lovop);
      lovop_ptr = dynamic_pointer_cast< LinCompLOVOp<ireal> >(tmp_ptr);

      //////////////////////////////
      // initializing G_lovop_ptr //
      //////////////////////////////
      //resetting RHS filenames on locpars
      for(int i=0; i<mks.RHS_files.size(); i++)
	add_to_pars(*locpars,mks.RHS_files[i],"empty");
      //adding MPS delta filename to locpars
      add_to_pars(*locpars,
		  mks.MPS_file,
		  delta_sp_ptr->get_coreFilename());

      MPS_to_RHS G_m2r(mks,*locpars,*delta_sp_ptr);

      //adding generated RHS files to locpars for IWaveOp initialization      
      add_to_pars(*locpars,
		  G_m2r.get_IWave_Range().getKeys(),
		  G_m2r.get_rng_filenames());

      //extracting filename corresponding to G_sp_ptr
      vector<string> G_filenames(G_sp_ptr->getSize());
      for( int i=0; i<G_filenames.size(); i++){
#ifdef IWAVE_USE_MPI
	MPISEGYSpace const* segy_sp =dynamic_cast<MPISEGYSpace const *>(&(*G_sp_ptr)[i]);
#else
	SEGYSpace const *segy_sp = dynamic_cast<SEGYSpace const *>(&(*G_sp_ptr)[i]);
#endif	
	G_filenames[i] = segy_sp->getPrototypeFilename();
      }
      
      add_to_pars(*locpars,
		  G_sp_ptr->getKeys(),
		  G_filenames );

      IWaveLOVOp G_iwop(*locpars,stream);
      LinCompLOVOp<ireal> G_lovop(G_m2r,G_iwop);
      shared_ptr< Operator<ireal> > G_tmp_ptr = Operator<ireal>::clonePtr(G_lovop);
      G_lovop_ptr = dynamic_pointer_cast< LinCompLOVOp<ireal> >(G_tmp_ptr);

      ps_delete(&locpars);

      /////////////////////////// 
      // initializing conv_ptr //
      ///////////////////////////
      MPS_conv tmp_conv(*pars,*G_ptrs[0],MPS_dom,rng);
      conv_ptr = MPS_conv::newPtr(tmp_conv);
      
    }
    catch(RVLException &e){
      e << "ERROR from MPS_IWaveLOVOp::set_ops\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  void MPS_IWaveLOVOp<T_MPS_Space>::apply0(const Vector<ireal> & x0,
					   const Vector<ireal> & x1,
					   Vector<ireal> & y) const{
  //----------------------------------------------------------------------//
    try{
#ifdef VERBOSE_MPS_IWOP_HH
      cerr << "\n<<<< Inside MPS_IWaveLOVOp::apply0 >>>>\n";
#endif

      SpaceTest(this->getNonLinDomain(),x0,"MPS_IWaveLOVOp::apply0 (nonlin dom)");
      SpaceTest(this->getLinDomain(),x1,"MPS_IWaveLOVOp::apply0 (lin dom)");
      SpaceTest(this->getRange(),y,"MPS_IWaveLOVOp::apply0 (rng)");
      
      y.zero();

      //checking int_x_ptr has been initialized
      if( !(int_x_ptr) ){
	RVLException e;
	e << "int_x_ptr is null!\n";
	throw e;
      }

      //recompute Greens functions if x0 is diff from internal x0
      Components<ireal> c_int_x(*int_x_ptr);
      Vector<ireal> tmp_x0(x0);
      tmp_x0.linComb(-1.0,c_int_x[0]);
      float diff_x0 = tmp_x0.norm()/x0.norm();

      if(diff_x0>1e-5){
	/*
	cerr << "Difference between x0 and int_x0 detected:\n"
	     << "|int_x0-x0|/|x0|="<<diff_x0<<"\n"
	     << "|x0|=" << x0.norm() << "\n"
	     << "|int_x0|=" << c_int_x[0].norm() << "\n";
	*/
	c_int_x[0].copy(x0);
	compute_G();
      }

      //computing output
      Vector<ireal> tmp_y(y);
      for(int i=0; i<G_ptrs.size(); i++){
	tmp_y.zero();
	conv_ptr->set_kernel(*G_ptrs[i]);
	conv_ptr->set_MPS_idx(i);
	conv_ptr->applyOp(x1,tmp_y);
	y.linComb(1.0,tmp_y);
      }
    }
    catch(RVLException &e){
      e << "ERROR from MPS_IWaveLOVOp::apply0\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  void MPS_IWaveLOVOp<T_MPS_Space>::applyAdj0(const Vector<ireal> & x0,
					      const Vector<ireal> & y,
					      Vector<ireal> & x1) const{
  //----------------------------------------------------------------------//
    try{
#ifdef VERBOSE_MPS_IWOP_HH
      cerr << "\n<<<< Inside MPS_IWaveLOVOp::applyAdj0 >>>>\n";
#endif

      SpaceTest(this->getNonLinDomain(),x0,"MPS_IWaveLOVOp::applyAdj0 (nonlin dom)");
      SpaceTest(this->getLinDomain(),x1,"MPS_IWaveLOVOp::applyAdj0 (lin dom)");
      SpaceTest(this->getRange(),y,"MPS_IWaveLOVOp::applyAdj0 (rng)");
      
      x1.zero();

      //checking int_x_ptr has been initialized
      if( !(int_x_ptr) ){
	RVLException e;
	e << "int_x_ptr is null!\n";
	throw e;
      }

      //recompute Greens functions if x0 is diff from internal x0
      Components<ireal> c_int_x(*int_x_ptr);
      Vector<ireal> tmp_x0(x0);
      tmp_x0.linComb(-1.0,c_int_x[0]);
      float diff_x0 = tmp_x0.norm()/x0.norm();

      if(diff_x0>1e-5){
	/*
      	cerr << "Difference between x0 and int_x0 detected:\n"
	     << "|int_x0-x0|/|x0|="<<diff_x0<<"\n"
	     << "|x0|=" << x0.norm() << "\n"
	     << "|int_x0|=" << c_int_x[0].norm() << "\n";
	*/

	c_int_x[0].copy(x0);
	compute_G();
      }
      
      //computing output
      Vector<ireal> tmp_x1(x1);

      for(int i=0; i<G_ptrs.size(); i++){
	tmp_x1.zero();
	conv_ptr->set_kernel(*G_ptrs[i]);
	conv_ptr->set_MPS_idx(i);	
	conv_ptr->applyAdjOp(y,tmp_x1);
	x1.linComb(1.0,tmp_x1);
      }
    }
    catch(RVLException &e){
      e << "ERROR from MPS_IWaveLOVOp::applyAdj0\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  void MPS_IWaveLOVOp<T_MPS_Space>::applyPartialDeriv0(const Vector<ireal> & x0,
						       const Vector<ireal> & x1,
						       const Vector<ireal> & dx0,
						       Vector<ireal> & dy) const{
  //----------------------------------------------------------------------//
    try{
      LinOpValOp<ireal>::export_applyPartialDeriv0(*lovop_ptr,
						   x0,
						   x1,
						   dx0,
						   dy);
    }
    catch(RVLException &e){
      e << "ERROR from MPS_IWaveLOVOp::applyPartialDeriv0\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  void MPS_IWaveLOVOp<T_MPS_Space>::applyAdjPartialDeriv0(const Vector<ireal> & x0,
							  const Vector<ireal> & x1,
							  const Vector<ireal> & dy,
							  Vector<ireal> & dx0) const{
  //----------------------------------------------------------------------//
    try{
      LinOpValOp<ireal>::export_applyAdjPartialDeriv0(*lovop_ptr,
						      x0,
						      x1,
						      dy,
						      dx0);
    }
    catch(RVLException &e){
      e << "ERROR from MPS_IWaveLOVOp::applyAdjPartialDeriv0\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  void MPS_IWaveLOVOp<T_MPS_Space>::applyPartialDeriv20(const Vector<ireal> & x0, 
							const Vector<ireal> & x1,
							const Vector<ireal> & dx00,
							Vector<ireal> dy) const {
  //----------------------------------------------------------------------//
    try{
      LinOpValOp<ireal>::export_applyPartialDeriv20(*lovop_ptr,
						    x0,
						    x1,
						    dx00,
						    dy);
    }
    catch(RVLException &e){
      e << "ERROR from MPS_IWaveLOVOp::applyPartialDeriv20\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  template<class T_MPS_Space>
  void MPS_IWaveLOVOp<T_MPS_Space>::applyAdjPartialDeriv20(const Vector<ireal> & x0, 
							   const Vector<ireal> & x1,
							   const Vector<ireal> & dy,
							   Vector<ireal> dx00) const {
  //----------------------------------------------------------------------//
    try{
      LinOpValOp<ireal>::export_applyAdjPartialDeriv20(*lovop_ptr,
						       x0,
						       x1,
						       dy,
						       dx00);
    }
    catch(RVLException &e){
      e << "ERROR from MPS_IWaveLOVOp::applyAdjPartialDeriv20\n";
      throw e;
    }
  }

}
#endif
