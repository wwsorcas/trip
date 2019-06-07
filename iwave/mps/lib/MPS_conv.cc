// MPS_conv.cc
// Author: Mario J. Bencomo
// last modified: 11/10/16

#include "MPS_conv.hh"

//#define VERBOSE_MPS_CONV_CC

namespace TSOpt{

  //----------------------------------------------------------------------// 
  void MPS_conv_kern( vector<Rtuple> s_pos, 
		      int MPS_dim, 
		      int MPS_idx,
		      string x_file,
		      vector<string> y_files,
		      vector<string> k_files ){
  //----------------------------------------------------------------------// 
    try{
      
#ifdef VERBOSE_MPS_CONV_CC
      cerr << "*************************************************\n"
           << "  Inside MPS_conv_kern\n"
           << "*************************************************\n\n";
      
      cerr << "x file = " << x_file << "\n"
	   << "y files = "; print_vec(y_files);
      cerr << "k files = "; print_vec(k_files);
      cerr << "MPS_dim = " << MPS_dim << "\n"
	   << "MPS_idx = " << MPS_idx << "\n";
#endif
      
      FILE *fp_x, *fp_y, *fp_k;
      segy tr_x, tr_y, tr_k;
      float * x_arr = NULL;
      float * y_arr = NULL;
      float * k_arr = NULL;
      int NTR_x, NTR_y, NTR_k;
      int NS_x, NS_y, NS_k;
      int NTR_src; //# traces per source
      int i_shift, i_shift_x, i_shift_y, i_shift_k;
      fpos_t pos_y;

      if( MPS_idx<0 || MPS_idx>=MPS_dim ){
	RVLException e;
	e << "MPS_idx="<< MPS_idx << " out of bounds!\n";
	throw e;
      }

      //getting size of input and setting up buffer
      fp_x = iwave_const_fopen(x_file.c_str(),"r",NULL,stderr);
      if( fp_x==NULL ){
	RVLException e;
	e << "Could not open input file "<< x_file <<"\n";
	throw e;
      }
      fgettr(fp_x,&tr_x);
      NTR_x = get_ntr(x_file); 
      NS_x = tr_x.ns; 
      iwave_fclose(fp_x); fp_x=NULL;

#ifdef VERBOSE_MPS_CONV_CC
      cerr << "NTR_x = " << NTR_x << "\n"
	   << "NS_x = " << NS_x << "\n";
#endif

      x_arr = new float[NTR_x*NS_x];
      if(x_arr==NULL){
	RVLException e;
	e << "Trouble allocating memory for input buffer array.\n";
	throw e;
      }
      
      //reading in input buffer
      fp_x = iwave_const_fopen(x_file.c_str(),"r",NULL,stderr);
      
      for( int i_tr=0; i_tr<NTR_x; i_tr++ ){
	fgettr(fp_x,&tr_x);
	for( int i_t=0; i_t<NS_x; i_t++ ){
	  x_arr[ i_t+i_tr*NS_x ] = tr_x.data[i_t];
	}
      }
      iwave_fclose(fp_x); fp_x=NULL;


      ///////////////////////////
      //looping over components//
      ///////////////////////////
      for(int i_c=0; i_c<y_files.size(); i_c++){
	
#ifdef VERBOSE_MPS_CONV_CC
	cerr << "\nInside components-loop: i_c="<<i_c<<"\n";
#endif
	
	//getting size of output and setting up buffer
	fp_y = iwave_const_fopen(y_files[i_c].c_str(),"r",NULL,stderr);
	if( fp_y==NULL ){
	  RVLException e;
	  e << "Could not open output file "<< y_files[i_c] <<"\n";
	  throw e;
	}
	fgettr(fp_y,&tr_y);
	NTR_y = get_ntr(y_files[i_c]);
	NS_y = tr_y.ns;
	iwave_fclose(fp_y); fp_y=NULL;

#ifdef VERBOSE_MPS_CONV_CC
	cerr << "NTR_y = " << NTR_y << "\n"
	     << "NS_y = "  << NS_y << "\n";
#endif

	y_arr = new float[NTR_y*NS_y];
	if(y_arr==NULL){
	  RVLException e;
	  e << "Trouble allocating memory for output buffer array.\n";
	  throw e;
	}

	//getting size of kernel and setting up buffer
	fp_k = iwave_const_fopen(k_files[i_c].c_str(),"r",NULL,stderr);
	if( fp_k==NULL ){
	  RVLException e;
	  e << "Could not open kernel file "<< k_files[i_c] <<"\n";
	  throw e;
	}
	fgettr(fp_k,&tr_k);
	NTR_k = get_ntr(k_files[i_c]);
	NS_k = tr_k.ns;
	iwave_fclose(fp_k); fp_k=NULL;

#ifdef VERBOSE_MPS_CONV_CC
	cerr << "NTR_k = " << NTR_k << "\n"
	     << "NS_k = " << NS_k << "\n";
#endif

	k_arr = new float[NTR_k*NS_k];
	if(k_arr==NULL){
	  RVLException e;
	  e << "Trouble allocating memory for kernel buffer array.\n";
	  throw e;
	}

	//reading in kernel buffer
	fp_k = iwave_const_fopen(k_files[i_c].c_str(),"r",NULL,stderr);
	
	for( int i_tr=0; i_tr<NTR_k; i_tr++ ){
	  fgettr(fp_k,&tr_k);
	  for( int i_t=0; i_t<NS_k; i_t++ ){
	    k_arr[ i_t+i_tr*NS_k ] = tr_k.data[i_t];
	  }
	}
	iwave_fclose(fp_k); fp_k=NULL;

	if( tr_x.dt!=tr_y.dt || tr_x.dt!=tr_k.dt ){
	  RVLException e;
	  e << "Trace dt's are not the same for input, output, or kernel for convolution!\n"
	    << "component idx="<< i_c <<"\n";
	  throw e;
	}
	
	//computing number of traces per source
	NTR_src = NTR_k/s_pos.size();

	////////////////////////
	//looping over sources//
	////////////////////////
	for( int i_s=0; i_s<s_pos.size(); i_s++){
	  
	  i_shift = i_s*MPS_dim + MPS_idx;
	  i_shift_x = i_shift*NS_x;

#ifdef VERBOSE_MPS_CONV_CC
	  cerr << "Inside source-loop: i_s="<<i_s<<"\n";
#endif

	  //////////////////////////
	  //looping over receivers//
	  //////////////////////////
	  for( int i_rc=0; i_rc<NTR_src; i_rc++ ){

	    i_shift = i_rc + NTR_src*i_s;
	    i_shift_y = i_shift*NS_y;
	    i_shift_k = i_shift*NS_k;
	    
	    //actual convolution
	    float dt   = tr_x.dt; 
	    int n_in   = NS_x;
	    int i0_in  = int( 1000 * tr_x.delrt/dt );
	    int n_out  = NS_y;
	    int i0_out = int( 1000 * tr_y.delrt/dt );
	    int n_ker  = NS_k;
	    int i0_ker = int( 1000 * tr_k.delrt/dt );
	    float scal = dt/1e6;
	    
#ifdef VERBOSE_MPS_CONV_CC
	    cerr << "    n_in   ="<< n_in   <<"\n"
		 << "    i0_in  ="<< i0_in  <<"\n"
		 << "    n_out  ="<< n_out  <<"\n"
		 << "    i0_out ="<< i0_out <<"\n"
		 << "    n_ker  ="<< n_ker  <<"\n"
		 << "    i0_ker ="<< i0_ker <<"\n\n"
		 << "    tr_x.delrt/dt="<<tr_x.delrt/dt<<"\n"
		 << "    tr_y.delrt/dt="<<tr_y.delrt/dt<<"\n"
		 << "    tr_k.delrt/dt="<<tr_k.delrt/dt<<"\n"
		 << "    scal = " << scal <<"\n\n"
		 << "    i_shift_x ="<<i_shift_x<<"\n"
		 << "    i_shift_y ="<<i_shift_y<<"\n"
		 << "    i_shift_k ="<<i_shift_k<<"\n";
#endif
	    	
	    int shift = i0_out - i0_in - i0_ker;
	    TSOpt::conv( shift, 
			 n_out, 
			 n_in, 
			 n_ker,
			 y_arr+i_shift_y,
			 x_arr+i_shift_x,
			 k_arr+i_shift_k,
			 scal );

	  }//receiver loop

	}//source loop

	//writing output buffer to file
	fp_y = iwave_const_fopen(y_files[i_c].c_str(),"r+",NULL,stderr);
	for( int i_tr=0; i_tr<NTR_y; i_tr++ ){
	  fgetpos(fp_y,&pos_y);
	  fgettr(fp_y,&tr_y);
	  for( int i_t=0; i_t<NS_y; i_t++ ){
	    tr_y.data[i_t] = y_arr[ i_t+i_tr*NS_y ];
	  }
	  fsetpos(fp_y,&pos_y);
	  fputtr(fp_y,&tr_y);
	}
	fflush(fp_y);
	iwave_fclose(fp_y); fp_y=NULL;

	delete[] y_arr;
	delete[] k_arr;

      }//component loop

      delete[] x_arr;

    }
    catch(RVLException &e){
      e << "ERROR from MPS_conv_kern!\n";
      throw e;
    }
  }


  //----------------------------------------------------------------------// 
  void MPS_corr_kern( vector<Rtuple> s_pos, 
		      int MPS_dim, 
		      int MPS_idx,
		      string x_file,
		      vector<string> y_files,
		      vector<string> k_files ){
  //----------------------------------------------------------------------// 
    try{
      
#ifdef VERBOSE_MPS_CONV_VV
      cerr << "*************************************************\n"
           << "  Inside MPS_corr_kern\n"
           << "*************************************************\n\n";
      
      cerr << "x file = " << x_file << "\n"
	   << "y files = "; print_vec(y_files);
      cerr << "k files = "; print_vec(k_files);
      cerr << "MPS_dim = " << MPS_dim << "\n"
	   << "MPS_idx = " << MPS_idx << "\n";
#endif

      FILE *fp_x, *fp_y, *fp_k;
      segy tr_x, tr_y, tr_k;
      float * x_arr = NULL;
      float * y_arr = NULL;
      float * k_arr = NULL;
      float * buff_arr = NULL;
      int NTR_x, NTR_y, NTR_k;
      int NS_x, NS_y, NS_k;
      int NTR_src; //# traces per source
      int i_shift, i_shift_x, i_shift_y, i_shift_k;
      fpos_t pos_x;

      if( MPS_idx<0 || MPS_idx>=MPS_dim ){
	RVLException e;
	e << "MPS_idx="<< MPS_idx << " out of bounds!\n";
	throw e;
      }

      //getting size of output and setting up buffer
      fp_x = iwave_const_fopen(x_file.c_str(),"r",NULL,stderr);
      if( fp_x==NULL ){
	RVLException e;
	e << "Could not open output file "<< x_file <<"\n";
	throw e;
      }
      fgettr(fp_x,&tr_x);
      NTR_x = get_ntr(x_file); 
      NS_x = tr_x.ns; 
      iwave_fclose(fp_x); fp_x=NULL;

#ifdef VERBOSE_MPS_CONV_CC
      cerr << "NTR_x = " << NTR_x << "\n"
	   << "NS_x = " << NS_x << "\n";
#endif

      x_arr = new float[NTR_x*NS_x];
      if(x_arr==NULL){
	RVLException e;
	e << "Trouble allocating memory for output buffer array.\n";
	throw e;
      }
      for(int i_tr=0; i_tr<NTR_x; i_tr++){
	for( int i_t=0; i_t<NS_x; i_t++){
	  x_arr[i_t+i_tr*NS_x] = 0.0;
	}
      }

      buff_arr = new float[NS_x];
      if(buff_arr==NULL){
	RVLException e;
	e << "Trouble allocating memory for output second buffer array.\n";
	throw e;
      }
      

      ///////////////////////////
      //looping over components//
      ///////////////////////////
      for(int i_c=0; i_c<y_files.size(); i_c++){
	
#ifdef VERBOSE_MPS_CONV_CC
	cerr << "Inside components-loop: i_c="<<i_c<<"\n";
#endif

	//getting size of input and setting up buffer
	fp_y = iwave_const_fopen(y_files[i_c].c_str(),"r",NULL,stderr);
	if( fp_y==NULL ){
	  RVLException e;
	  e << "Could not open input file "<< y_files[i_c] <<"\n";
	  throw e;
	}
	fgettr(fp_y,&tr_y);
	NTR_y = get_ntr(y_files[i_c]);
	NS_y = tr_y.ns;
	iwave_fclose(fp_y); fp_y=NULL;

#ifdef VERBOSE_MPS_CONV_CC
	cerr << "NTR_y = " << NTR_y << "\n"
	     << "NS_y = "  << NS_y << "\n";
#endif

	y_arr = new float[NTR_y*NS_y];
	if(y_arr==NULL){
	  RVLException e;
	  e << "Trouble allocating memory for input buffer array.\n";
	  throw e;
	}

	//reading in input buffer
	fp_y = iwave_const_fopen(y_files[i_c].c_str(),"r",NULL,stderr);
	for( int i_tr=0; i_tr<NTR_y; i_tr++ ){
	  fgettr(fp_y,&tr_y);
	  for( int i_t=0; i_t<NS_y; i_t++ ){
	    y_arr[ i_t+i_tr*NS_y ] = tr_y.data[i_t];
	  }
	}
	iwave_fclose(fp_y); fp_y=NULL;

	//getting size of kernel and setting up buffer
	fp_k = iwave_const_fopen(k_files[i_c].c_str(),"r",NULL,stderr);
	if( fp_k==NULL ){
	  RVLException e;
	  e << "Could not open kernel file "<< k_files[i_c] <<"\n";
	  throw e;
	}
	fgettr(fp_k,&tr_k);
	NTR_k = get_ntr(k_files[i_c]);
	NS_k = tr_k.ns;
	iwave_fclose(fp_k); fp_k=NULL;

#ifdef VERBOSE_MPS_CONV_CC
	cerr << "NTR_k = " << NTR_k << "\n"
	     << "NS_k = " << NS_k << "\n";
#endif

	k_arr = new float[NTR_k*NS_k];
	if(k_arr==NULL){
	  RVLException e;
	  e << "Trouble allocating memory for kernel buffer array.\n";
	  throw e;
	}

	//reading in kernel buffer
	fp_k = iwave_const_fopen(k_files[i_c].c_str(),"r",NULL,stderr);
	
	for( int i_tr=0; i_tr<NTR_k; i_tr++ ){
	  fgettr(fp_k,&tr_k);
	  for( int i_t=0; i_t<NS_k; i_t++ ){
	    k_arr[ i_t+i_tr*NS_k ] = tr_k.data[i_t];
	  }
	}
	iwave_fclose(fp_k); fp_k=NULL;

	//checking dt's
	if( tr_x.dt!=tr_y.dt || tr_x.dt!=tr_k.dt ){
	  RVLException e;
	  e << "Trace dt's are not the same for input, output, or kernel for convolution!\n"
	    << "component idx="<< i_c <<"\n";
	  throw e;
	}
	
	//computing number of traces per source
	NTR_src = NTR_k/s_pos.size();


	////////////////////////
	//looping over sources//
	////////////////////////
	for( int i_s=0; i_s<s_pos.size(); i_s++){

	  i_shift = i_s*MPS_dim + MPS_idx;
	  i_shift_x = i_shift*NS_x;	  

#ifdef VERBOSE_MPS_CONV_CC
	  cerr << "Inside source-loop: i_s="<<i_s<<"\n";
	  cerr << "i_shift_x="<<i_shift_x<<"\n";
#endif
	  
	  //////////////////////////
	  //looping over receivers//
	  //////////////////////////
	  for( int i_rc=0; i_rc<NTR_src; i_rc++ ){

	    i_shift = i_rc + NTR_src*i_s;
	    i_shift_y = i_shift*NS_y;
	    i_shift_k = i_shift*NS_k;

	    //actual cross-correlation
	    float dt  = tr_x.dt;
	    int n_in   = NS_y;
	    int i0_in  = int( 1000 * tr_y.delrt/dt );
	    int n_out  = NS_x;
	    int i0_out = int( 1000 * tr_x.delrt/dt );
	    int n_ker  = NS_k;
	    int i0_ker = int( 1000 * tr_k.delrt/dt );
	    float scal = dt/1e6;

#ifdef VERBOSE_MPS_CONV_CC
	    cerr << "    n_in   ="<< n_in   <<"\n"
		 << "    i0_in  ="<< i0_in  <<"\n"
		 << "    n_out  ="<< n_out  <<"\n"
		 << "    i0_out ="<< i0_out <<"\n"
		 << "    n_ker  ="<< n_ker  <<"\n"
		 << "    i0_ker ="<< i0_ker <<"\n"
		 << "    tr_y.delrt/dt="<<tr_y.delrt/dt<<"\n"
		 << "    tr_x.delrt/dt="<<tr_x.delrt/dt<<"\n"
		 << "    tr_k.delrt/dt="<<tr_k.delrt/dt<<"\n"
		 << "    scal = " << scal <<"\n";
#endif
	    
	    int shift = i0_in - i0_out - i0_ker;
	    TSOpt::corr( shift,
			 n_out,
			 n_in,
			 n_ker,
			 buff_arr,
			 y_arr+i_shift_y,
			 k_arr+i_shift_k,
			 scal );
			
	    for(int i_t=0; i_t<NS_x; i_t++){
	      x_arr[i_t+i_shift_x] += buff_arr[i_t];
	    }

	  }//receiver loop
	    
	}//source loop
	
	delete[] y_arr;
	delete[] k_arr;

      }//component loop
      
      //writing output buffer to file
      fp_x = iwave_const_fopen(x_file.c_str(),"r+",NULL,stderr);
      for( int i_s=0; i_s<s_pos.size(); i_s++ ){

	i_shift = i_s*MPS_dim + MPS_idx;
	i_shift_x = i_shift*NS_x;	

	if(i_shift!=0) fgettra(fp_x,&tr_x,i_shift-1);
	fgetpos(fp_x,&pos_x);
	fgettr(fp_x,&tr_x);

	for( int i_t=0; i_t<NS_x; i_t++ ){
	    tr_x.data[i_t] = x_arr[ i_t+i_shift_x ];
	}
	fsetpos(fp_x,&pos_x);
	fputtr(fp_x,&tr_x);
      }
      fflush(fp_x);
      iwave_fclose(fp_x); fp_x=NULL;
      
      delete[] x_arr;
      delete[] buff_arr;
    }
    catch(RVLException &e){
      e << "ERROR from MPS_conv_kern!\n";
      throw e;
    }
  }



  //////////////////////////////////////////////////////////////////////////
  // MPS_conv //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  
  //----------------------------------------------------------------------// 
  MPS_conv::MPS_conv( MPS_conv const & conv )
    : pars(NULL), 
      ker_ptr(conv.ker_ptr), 
      ker_keys(conv.ker_keys),
      dom(conv.dom),
      rng(conv.rng),
      MPS_idx(conv.MPS_idx) {
  //----------------------------------------------------------------------// 
    try{
      //deep copy of pars
      pars = ps_new();
      if(ps_copy(&pars,*conv.pars)){
	RVLException e;
	e << "failed to copy parameter table\n";
	throw e;
      }
    }
    catch(RVLException & e){
      e << "ERROR from MPS_conv copy constructor\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------// 
  MPS_conv::MPS_conv( PARARRAY const &_pars, 
		      Vector<ireal> &k,
		      MPS_Space const &d,
		      IWaveSpace const &r )
    : pars(NULL), 
      ker_ptr(), 
      ker_keys(0),
      dom(d),
      rng(r),
      MPS_idx(-1) {
  //----------------------------------------------------------------------// 
    try{
 
      //copying pars
      pars = ps_new();
      if(ps_copy(&pars,_pars)){
	RVLException e;
	e << "failed to copy parameter table\n";
	throw e;
      }

      //initializing and setting ker
      set_kernel(k);

      //running sanity check
      sanity_checks();

      //extracting kernel keys 
      IWaveSpace const* sp_ker 
	= dynamic_cast<IWaveSpace const*>(&(ker_ptr->getSpace()));
      for( int i=0; i<sp_ker->getSize(); i++){
	string key = sp_ker->getKeys()[i];
	key += "_ker";
	ker_keys.push_back(key);
      }
      
    }
    catch(RVLException & e){
      e << "ERROR from MPS_conv constructor\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------// 
  MPS_conv::~MPS_conv(){
  //----------------------------------------------------------------------// 
    ps_delete(&pars);
  }

  //----------------------------------------------------------------------// 
  void MPS_conv::set_kernel(Vector<ireal> const &k){
  //----------------------------------------------------------------------// 
    try{
      ker_ptr = Vector<ireal>::newPtr(k.getSpace());
      change_kernel(k);
    }
    catch(RVLException &e){
      e << "ERROR from MPS_conv::set_kernel\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void MPS_conv::change_kernel(Vector<ireal> const &k){
  //----------------------------------------------------------------------//
    try{
      SpaceTest(ker_ptr->getSpace(),k,"MPS_conv::change_kernel (kern sp)");
      ker_ptr->copy(k);
      //sanity_checks();
    }
    catch(RVLException &e){
      e << "ERROR from MPS_conv::change_kernel\n";
      throw e;
    }
  }


  //----------------------------------------------------------------------// 
  void MPS_conv::sanity_checks(){
  //----------------------------------------------------------------------// 
    try{
      
      //Checking that size of greens and rng spaces matches.
      IWaveSpace const* sp_ker 
	= dynamic_cast<IWaveSpace const*>(&(ker_ptr->getSpace()));
    
      if(sp_ker->getSize()!=rng.getSize()){
	RVLException e;
	e << "size of kernel space and rng space does not match!\n";
	throw e;
      }
      
      
      //Checking that MPS_Space has same sources as
      //green's kernel, and rng IWaveSpace.
      
      //number of sources from MPS_Space domain
      vector<Rtuple> s_pos_MPS = dom.get_s_pos();
      
      /////////////////////////////////
      //looping over subspaces of rng//
      /////////////////////////////////
      for( int i=0; i<rng.getSize(); i++){

	//extracting src pos for rng space i
#ifdef IWAVE_USE_MPI
	MPISEGYSpace const *segy_rng 
	  = dynamic_cast<MPISEGYSpace const *>(&(rng[i]));
#else
        SEGYSpace const *segy_rng
	  = dynamic_cast<SEGYSpace const *>(&(rng[i]));
#endif
      
	if(segy_rng==NULL){
	  RVLException e;
	  e << "trouble getting segy space from rng["<<i<<"] space\n";
	  throw e;
	}
	string file_rng = segy_rng->getPrototypeFilename();
	vector<Rtuple> s_pos_rng = ex_pos(file_rng,false);
	
	//extracting src pos for green sp i
#ifdef IWAVE_USE_MPI
	MPISEGYSpace const *segy_ker 
	  = dynamic_cast<MPISEGYSpace const *>(&((*sp_ker)[i]));
#else	
	SEGYSpace const *segy_ker 
	  = dynamic_cast<SEGYSpace const*>(&((*sp_ker)[i]));
#endif
	if(segy_ker==NULL){
	  RVLException e;
	  e << "trouble getting segy space from sp_ker["<<i<<"] space\n";
	  throw e;
	}
	string file_ker = segy_ker->getPrototypeFilename();
	vector<Rtuple> s_pos_ker = ex_pos(file_ker,false);
	
	//checking number of sources
	if( s_pos_MPS.size()!=s_pos_rng.size() ||
	    s_pos_MPS.size()!=s_pos_ker.size() ){
	  RVLException e;
	  e << "Number of sources does not match!\n"
	    << "Nsrc in MPS = "<< s_pos_MPS.size() << "\n"
	    << "Nsrc in rng = "<< s_pos_rng.size() << "\n"
	    << "Nsrc in ker = "<< s_pos_ker.size() << "\n";
	  throw e;
	}
	
	//checking src pos
	for( int j=0; j<s_pos_MPS.size(); j++){
	  if( s_pos_MPS[j].coor[0]!=s_pos_rng[j].coor[0] ||
	      s_pos_MPS[j].coor[1]!=s_pos_rng[j].coor[1] ||
	      s_pos_MPS[j].coor[2]!=s_pos_rng[j].coor[2] ){
	    RVLException e;
	    e << "Source j="<<j<<" does not match between domain and range spaces!\n"
	      << "dom->s_pos["<<j<<"].coor[0]="<<s_pos_MPS[j].coor[0]<<"\n"
	      << "dom->s_pos["<<j<<"].coor[1]="<<s_pos_MPS[j].coor[1]<<"\n"
	      << "dom->s_pos["<<j<<"].coor[2]="<<s_pos_MPS[j].coor[2]<<"\n";
	    throw e;
	  }
	  if( s_pos_MPS[j].coor[0]!=s_pos_ker[j].coor[0] ||
	      s_pos_MPS[j].coor[1]!=s_pos_ker[j].coor[1] ||
	      s_pos_MPS[j].coor[2]!=s_pos_ker[j].coor[2] ){
	    RVLException e;
	    e << "Source j="<<j<<" does not match between domain and kernel spaces!\n";
	    throw e;
	  }
	}
	
	//checking that number of traces in kernel and range space matches
	int NTR_rng = get_ntr(file_rng);
	int NTR_ker = get_ntr(file_ker);
	
	if(NTR_rng!=NTR_ker){
	  RVLException e;
	  e << "Number of traces does not match between kernel and rng space["<<i<<"]\n"
	    << "file_rng="<<file_rng<<", NTR_rng="<<NTR_rng<<"\n"
	    << "file_ker="<<file_ker<<", NTR_ker="<<NTR_ker<<"\n";
	  throw e;
	}
      }
    } 
    catch(bad_cast){
      RVLException e;
      e << "bad cast!\n";
      throw e;
    }
    catch(RVLException &e){
      e << "ERROR from MPS_conv::sanity_checks()\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void MPS_conv::param_set( Vector<ireal> const &x,
			    Vector<ireal> const &y ) const{
  //----------------------------------------------------------------------//  
    try{
      
      AssignParams ap_x(*pars,dom.get_MPS_keys().MPS_file,stderr);
      x.eval(ap_x);

      Components<ireal> cy(y);
      for (size_t i=0; i<rng.getKeys().size(); i++) {
	AssignParams ap_y(*pars,rng.getKeys()[i],stderr);
	cy[i].eval(ap_y);
      }

      Components<ireal> ck(*ker_ptr);
      for (size_t i=0; i<ker_keys.size(); i++) {
	AssignParams ap_k(*pars,ker_keys[i],stderr);
	ck[i].eval(ap_k);
      }

    }
    catch(RVLException &e){
      e << "ERROR from MPS_conv::param_set\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  ostream & MPS_conv::write(ostream & str) const {
  //----------------------------------------------------------------------//
    str<<"MPS_conv object\n";
    str<<"domain:\n";
    dom.write(str);
    str<<"range:\n";
    rng.write(str);
    str<<"kernel:\n";
    ker_ptr->write(str);

    return str; 
  }  

  //----------------------------------------------------------------------//
  void MPS_conv::apply( const Vector<ireal> &x,
			Vector<ireal> &y ) const{
  //----------------------------------------------------------------------//
    try{
     
      if( !(ker_ptr) ){
	RVLException e;
	e << "Empty kernel vector!\n";
	throw e;
      }

      SpaceTest(this->get_MPS_Domain(),x,"MPS_conv::apply (dom)");
      SpaceTest(this->get_IWave_Range(),y,"MPS_conv::apply (rng)");

      param_set(x,y);

      // zero output 
      y.zero();

      //extracting file names
      string x_file = valparse<string>(*pars,dom.get_MPS_keys().MPS_file);
      vector<string> y_files(0);
      vector<string> k_files(0);
      for(int i=0; i<ker_keys.size(); i++){
	y_files.push_back(valparse<string>(*pars,rng.getKeys()[i]));
	k_files.push_back(valparse<string>(*pars,ker_keys[i]));
      }

      //convolution kernel
#ifdef IWAVE_USE_MPI
      MPS_conv_kern_MPI( dom.get_s_pos(),
			 dom.get_dim(),
			 MPS_idx,
			 x_file,
			 y_files,
			 k_files );
#else
      MPS_conv_kern( dom.get_s_pos(),
		     dom.get_dim(),
		     MPS_idx,
		     x_file,
		     y_files,
		     k_files );
#endif

    }
    catch (RVLException & e) {
      e<<"ERROR from MPS_conv::apply\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void MPS_conv::applyAdj( const Vector<ireal> &y,
			   Vector<ireal> &x ) const {
  //----------------------------------------------------------------------//
    try{
     
      if( !(ker_ptr) ){
	RVLException e;
	e << "Empty kernel vector!\n";
	throw e;
      }

      SpaceTest(this->get_MPS_Domain(),x,"MPS_conv::applyAdj (dom)");
      SpaceTest(this->get_IWave_Range(),y,"MPS_conv::applyAdj (rng)");
      param_set(x,y);

      // zero output 
      x.zero();

      //extracting file names
      string x_file = valparse<string>(*pars,dom.get_MPS_keys().MPS_file);
      vector<string> y_files(0);
      vector<string> k_files(0);
      for(int i=0; i<ker_keys.size(); i++){
	y_files.push_back(valparse<string>(*pars,rng.getKeys()[i]));
	k_files.push_back(valparse<string>(*pars,ker_keys[i]));
      }

      //cross-correlation kernel
#ifdef IWAVE_USE_MPI
      MPS_corr_kern_MPI( dom.get_s_pos(),
			 dom.get_dim(),
			 MPS_idx,
			 x_file,
			 y_files,
			 k_files );
#else
      MPS_corr_kern( dom.get_s_pos(),
		     dom.get_dim(),
		     MPS_idx,
		     x_file,
		     y_files,
		     k_files );
#endif

    }
    catch (RVLException & e) {
      e<<"ERROR from MPS_conv::applyAdj\n";
      throw e;
    }

  }



}
