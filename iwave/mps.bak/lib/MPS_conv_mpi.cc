// MPS_conv_mpi.cc
// Author: Mario J. Bencomo
// last modified: 11/10/16

#include "MPS_conv.hh"

//#define VERBOSE_MPS_CONV_MPI_CC
#ifdef IWAVE_USE_MPI

namespace TSOpt{

  //----------------------------------------------------------------------// 
  void MPS_conv_kern_MPI( vector<Rtuple> s_pos, 
			  int MPS_dim, 
			  int MPS_idx,
			  string x_file,
			  vector<string> y_files,
			  vector<string> k_files ){
  //----------------------------------------------------------------------// 
    try{
      
#ifdef VERBOSE_MPS_CONV_MPI_CC
      if( retrieveGlobalRank()==0 ){
	cerr << "*************************************************\n"
	     << "  Inside MPS_conv_kern_MPI\n"
	     << "*************************************************\n\n";
	
	cerr << "x file = " << x_file << "\n"
	     << "y files = "; print_vec(y_files);
	cerr << "k files = "; print_vec(k_files);
	cerr << "MPS_dim = " << MPS_dim << "\n"
	     << "MPS_idx = " << MPS_idx << "\n";
      }
      MPI_Barrier(retrieveGlobalComm());
#endif
      
      FILE *fp_x, *fp_y, *fp_k;
      segy tr_x, tr_y, tr_k;
      
      //MPI variables
      int RANK = retrieveGlobalRank();
      int rank = retrieveRank();
      int SIZE = retrieveGlobalSize();
      int size = retrieveSize();
      int N_gr = retrieveNumGroups();
      int grID = retrieveGroupID();
      MPI_Comm COMM = retrieveGlobalComm();
      MPI_Comm comm = retrieveComm();
      
      int N_src = s_pos.size(); //total number of sources
      int N_src_loc; //number of sources at current process
      
      int NS_x,      NS_y,      NS_k; //total number of time sample points per trace
      int NTR_x,     NTR_y,     NTR_k; //total number of traces
      int NTR_x_loc, NTR_y_loc, NTR_k_loc; //total number of traces at current process
      int NTPS; //number of traces per source in kernel and y 

      float dt_x,    dt_y,    dt_k;
      float delrt_x, delrt_y, delrt_k;

      //global arrays
      float * x_arr = NULL; 
      float * y_arr = NULL;
      float * k_arr = NULL;
      //local arrays
      float * x_arr_loc = NULL;
      float * y_arr_loc = NULL;
      float * k_arr_loc = NULL;
      //other arrays
      int *NSPG  = NULL; //array containing number of sources per group, for each group
      int *DATA_SZ_k  = NULL; //array containing size of data related to k_arr to be sent via scatterv
      int *DATA_OFF_k = NULL; //array of offsets for data related to k_arr to be sent via scatterv
      int *DATA_SZ_y  = NULL; //array containing size of data related to y_arr to be received via gatherv
      int *DATA_OFF_y = NULL; //array of offsets for data related to y_arr to be received via gatherv
      fpos_t pos_y;

      
      if( MPS_idx<0 || MPS_idx>=MPS_dim ){
	RVLException e;
	e << "MPS_idx="<< MPS_idx << " out of bounds!\n";
	throw e;
      }
            
      //setting up number of sources per group (NSPG) array
      NSPG = new int[N_gr];
      for( int g=0; g<N_gr; g++ ){
	NSPG[g] = ( g<(N_src%N_gr) )?
	  ( N_src/N_gr + 1 ):
	  ( N_src/N_gr );
      }
      N_src_loc = NSPG[grID];
      
      
      //////////////////////
      // setting up x_arr //
      //////////////////////
      NTR_x = get_ntr(x_file); 

      if(MPS_dim!=(NTR_x/N_src)){
	RVLException e;
	e << "MPS_dim="<<MPS_dim<<" does not match number of traces per source,  NTR_x/N_src="
	  << NTR_x/N_src<<"\n";
	throw e;
      }

      if( RANK==0 ){

	//getting size
	fp_x = iwave_const_fopen(x_file.c_str(),"r",NULL,stderr);
	if( fp_x==NULL ){
	  RVLException e;
	  e << "Could not open input file "<< x_file <<"\n";
	  throw e;
	}
	fgettr(fp_x,&tr_x);
	iwave_fclose(fp_x); fp_x=NULL;

	NS_x = tr_x.ns; 
	dt_x = tr_x.dt;
	delrt_x = tr_x.delrt;
	
	x_arr = new float[NS_x*N_src];
	if(x_arr==NULL){
	  RVLException e;
	  e << "Trouble allocating memory for input buffer array.\n";
	  throw e;
	}
      
	//reading in input buffer
	int i_tr;
	fp_x = iwave_const_fopen(x_file.c_str(),"r",NULL,stderr);     
	for( int i_src=0; i_src<N_src; i_src++ ){
	  i_tr = i_src*MPS_dim + MPS_idx;
	  fgettra(fp_x,&tr_x,i_tr);

	  for( int i_t=0; i_t<NS_x; i_t++ ){
	    x_arr[ i_t+i_src*NS_x ] = tr_x.data[i_t];
	  }
	}
	iwave_fclose(fp_x); fp_x=NULL;

      }//end RANK==0 if

      //broadcasting 
      MPI_Bcast( &NS_x, 1, MPI_INT, 0, COMM );
      MPI_Bcast( &dt_x, 1, MPI_FLOAT, 0, COMM );
      MPI_Bcast( &delrt_x, 1, MPI_FLOAT, 0, COMM );

#ifdef VERBOSE_MPS_CONV_MPI_CC
      for( int r=0; r<SIZE; r++){
	if(RANK==r){
	  cerr <<"rank="<<r<<", After setting up x_arr.\n"
	       <<"  NTR_x   = "<<NTR_x<<"\n"
	       <<"  NS_x    = "<<NS_x<<"\n"
	       <<"  dt_x    = "<<dt_x<<"\n"
	       <<"  delrt_x = "<<delrt_x<<"\n";
	}
	MPI_Barrier(COMM);
      }
#endif

      ////////////////////////
      // distributing x_arr //
      ////////////////////////
      NTR_x_loc = N_src_loc;
      x_arr_loc = new float[NS_x*NTR_x_loc];
      
      //distributing x_arr to group root processes
      int tag = 0;
      if( RANK==0 ){

	int dest;
	int offset_g=0; //trace offset due to group ID
	int offset = 0; //total offset

	//looping over groups
	for( int g=0; g<N_gr; g++ ){
	  dest = g*size;
	  
	  if(dest==0){
	    for( int i=0; i<NS_x*NTR_x_loc; i++ )
		x_arr_loc[i] = x_arr[i];
	  }
	  else{	      	
	    offset = offset_g*NS_x;
	    MPI_Send( x_arr+offset,
		      NS_x*NSPG[g],
		      MPI_FLOAT,
		      dest,
		      tag,
		      COMM );
	  }
	  offset_g += NSPG[g];
	}
      }
      else if( rank==0 ){
	MPI_Recv( x_arr_loc,
		  NS_x*NTR_x_loc,
		  MPI_FLOAT,
		  0,
		  tag,
		  COMM,
		  MPI_STATUS_IGNORE );
      }

#ifdef VERBOSE_MPS_CONV_MPI_CC
      for( int r=0; r<SIZE; r++){
	if(RANK==r){
	  cerr <<"rank="<<r<<", After distributing x_arr to group roots\n"
	       <<"  NTR_x_loc="<<NTR_x_loc<<"\n";
	}
	MPI_Barrier(COMM);
      }
#endif

      //redistributing from group root processes to the rest
      MPI_Bcast( x_arr_loc, NTR_x_loc*NS_x, MPI_FLOAT, 0, comm );
      

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      // looping over RHS components //
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      for(int i_c=0; i_c<y_files.size(); i_c++){
	
#ifdef VERBOSE_MPS_CONV_MPI_CC
	for( int r=0; r<SIZE; r++){
	  if(RANK==r){
	    cerr << "rank="<<r<<", Inside components-loop: i_c="<<i_c<<"\n";	
	  }
	  MPI_Barrier(COMM);
	}
#endif
      
	//////////////////////
	// setting up k_arr //
	//////////////////////
	NTR_k = get_ntr(k_files[i_c]);
	NTPS = NTR_k/N_src;

	if( RANK==0 ){
	  
	  //getting size
	  fp_k = iwave_const_fopen(k_files[i_c].c_str(),"r",NULL,stderr);
	  if( fp_k==NULL ){
	    RVLException e;
	    e << "Could not open kernel file "<< k_files[i_c] <<"\n";
	    throw e;
	  }
	  fgettr(fp_k,&tr_k);
	  iwave_fclose(fp_k); fp_k=NULL;
	  
	  NS_k = tr_k.ns;
	  dt_k = tr_k.dt;
	  delrt_k = tr_k.delrt;
	  
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
	
	}//end RANK==0 if

	//broadcasting
	MPI_Bcast( &NS_k, 1, MPI_INT, 0, COMM );
	MPI_Bcast( &dt_k, 1, MPI_FLOAT, 0, COMM );
	MPI_Bcast( &delrt_k, 1, MPI_FLOAT, 0, COMM );
	
#ifdef VERBOSE_MPS_CONV_MPI_CC
	for( int r=0; r<SIZE; r++){
	  if(RANK==r){
	    cerr <<"rank="<<r<<", After setting up k_arr.\n"
		 <<"  NTPS    = "<<NTPS<<"\n"
		 <<"  NTR_k   = "<<NTR_k<<"\n"
		 <<"  NS_k    = "<<NS_k<<"\n"
		 <<"  dt_k    = "<<dt_k<<"\n"
		 <<"  delrt_k = "<<delrt_k<<"\n";
	  }
	  MPI_Barrier(COMM);
	}
#endif
      	
	////////////////////////
	// distributing k_arr //
	////////////////////////
	NTR_k_loc = N_src_loc*NTPS;
	k_arr_loc = new float[NS_k*NTR_k_loc];

	DATA_SZ_k = new int[SIZE];
	for( int g=0; g<N_gr; g++ ){
	  for( int r=0; r<size; r++ ){
	    int tmp = NSPG[g]*NTPS; //number of traces for group g	    
	    DATA_SZ_k[ r + g*size ] = ( r<(tmp%size) )?//number of traces for rank k in group g
	      ( tmp/size + 1 ):
	      ( tmp/size );
	    DATA_SZ_k[ r + g*size ] *= NS_k;
	  }
	}
	DATA_OFF_k = new int[SIZE];
	for( int r=0; r<SIZE; r++ ){
	  DATA_OFF_k[r]=0;
	  for( int k=0; k<r; k++ ){
	    DATA_OFF_k[r] += DATA_SZ_k[k];
	  }
	}

	MPI_Scatterv( k_arr, DATA_SZ_k, DATA_OFF_k, MPI_FLOAT,
		      k_arr_loc, NS_k*NTR_k_loc, MPI_FLOAT,
		      0, COMM );

#ifdef VERBOSE_MPS_CONV_MPI_CC
	for( int r=0; r<SIZE; r++){
	  if(RANK==r){
	    cerr <<"rank="<<r<<", After distributing (scatterv) k_arr.\n"
		 <<"  NTR_k_loc = "<<NTR_k_loc<<"\n"
		 <<"  recvcount = NS_k*NTR_k_loc = "<<NS_k*NTR_k_loc<<"\n";
	  }
	  MPI_Barrier(COMM);
	}
	
	if(RANK==0){
	  cerr <<"printing out DATA_SZ_k:\n";
	  for(int i=0; i<SIZE; i++)
	    cerr <<"  DATA_SZ_k["<<i<<"]="<<DATA_SZ_k[i]<<"\n";
       
	  cerr <<"printing out DATA_OFF_k:\n";
	  for(int i=0; i<SIZE; i++)
	    cerr <<"  DATA_OFF_k["<<i<<"]="<<DATA_OFF_k[i]<<"\n";	 
	}
	MPI_Barrier(COMM);
#endif
	
	
	////////////////////////////////////
	// setting up y_arr and y_arr_loc //
	////////////////////////////////////
	NTR_y = get_ntr(y_files[i_c]);

	if( RANK==0 ){

	  //getting size
	  fp_y = iwave_const_fopen(y_files[i_c].c_str(),"r",NULL,stderr);
	  if( fp_y==NULL ){
	    RVLException e;
	    e << "Could not open output file "<< y_files[i_c] <<"\n";
	    throw e;
	  }
	  fgettr(fp_y,&tr_y);
	  iwave_fclose(fp_y); fp_y=NULL;
	  
	  NS_y = tr_y.ns;
	  dt_y = tr_y.dt;
	  delrt_y = tr_y.delrt;
	  
	  y_arr = new float[NTR_y*NS_y];
	  if(y_arr==NULL){
	    RVLException e;
	    e << "Trouble allocating memory for output buffer array.\n";
	    throw e;
	  }

	}//end RANK==0 if
      
	//broadcasting 
	MPI_Bcast( &NS_y, 1, MPI_INT, 0, COMM );
	MPI_Bcast( &dt_y, 1, MPI_FLOAT, 0, COMM );
	MPI_Bcast( &delrt_y, 1, MPI_FLOAT, 0, COMM );
	
	NTR_y_loc = N_src_loc * NTPS;
	y_arr_loc = new float[NS_y*NTR_y_loc];


#ifdef VERBOSE_MPS_CONV_MPI_CC
	for( int r=0; r<SIZE; r++){
	  if(RANK==r){
	    cerr <<"rank="<<r<<", After setting up y_arr and y_arr_loc.\n"
		 <<"  NTR_y     = "<<NTR_y<<"\n"
		 <<"  NTR_y_loc = "<<NTR_y_loc<<"\n"
		 <<"  NS_y      = "<<NS_y<<"\n"
		 <<"  dt_y      = "<<dt_y<<"\n"
		 <<"  delrt_y   = "<<delrt_y<<"\n";
	  }
	  MPI_Barrier(COMM);
	}
#endif
	
	
	//some sanity checks
	if( dt_x!=dt_y || dt_x!=dt_k ){
	  RVLException e;
	  e << "Trace dt's are not the same for input, output, or kernel for convolution!\n"
	    << "component idx="<< i_c <<"\n";
	  throw e;
	}	  
	if( NTR_k!=NTR_y ){
	  RVLException e;
	  e << "NTR is not the same between kernel and output SUs!\n";
	  throw e;
	}
      
       
	//////////////////////////
	// args for convolution //
	//////////////////////////

	float dt   = dt_x; 
	int n_in   = NS_x;
	int i0_in  = int( 1000 * tr_x.delrt/dt );
	int n_out  = NS_y;
	int i0_out = int( 1000 * tr_y.delrt/dt );
	int n_ker  = NS_k;
	int i0_ker = int( 1000 * tr_k.delrt/dt );
	float scal = dt/1e6;
	    
#ifdef VERBOSE_MPS_CONV_MPI_CC
	if( RANK==0 ){
	  cerr << "\nConvolution args:\n"
	       << "    n_in   ="<< n_in   <<"\n"
	       << "    i0_in  ="<< i0_in  <<"\n"
	       << "    n_out  ="<< n_out  <<"\n"
	       << "    i0_out ="<< i0_out <<"\n"
	       << "    n_ker  ="<< n_ker  <<"\n"
	       << "    i0_ker ="<< i0_ker <<"\n"
	       << "    tr_x.delrt/dt="<<tr_x.delrt/dt<<"\n"
	       << "    tr_y.delrt/dt="<<tr_y.delrt/dt<<"\n"
	       << "    tr_k.delrt/dt="<<tr_k.delrt/dt<<"\n"
	       << "    scal = " << scal <<"\n";
	}
	MPI_Barrier(COMM);
#endif
      
	//~~~~~~~~~~~~~~~~~~~~~~//
	// looping over traces  //
	//~~~~~~~~~~~~~~~~~~~~~~//
	for( int i_tr=0; i_tr<NTR_k_loc; i_tr++ ){

	  //current source index, local to group
	  int i_src_loc = ( i_tr + DATA_OFF_k[RANK]/NS_k )/NTPS;
	  for( int g=0; g<grID; g++ ) i_src_loc -= NSPG[g];	  

	  //actual convolution
	  int shift = i0_out - i0_in - i0_ker;
	  TSOpt::conv( shift, 
		       n_out, 
		       n_in, 
		       n_ker,
		       y_arr_loc+i_tr*NS_y,
		       x_arr_loc+i_src_loc*NS_x,
		       k_arr_loc+i_tr*NS_k,
		       scal );
	  
	}//trace loop
	MPI_Barrier(COMM);

	
	////////////////////
	// gathering data //
	////////////////////
	DATA_SZ_y = new int[SIZE];
	for( int g=0; g<N_gr; g++ ){
	  for( int r=0; r<size; r++ ){
	    int tmp = NSPG[g]*NTPS; //number of traces for group g	    
	    DATA_SZ_y[ r + g*size ] = ( r<(tmp%size) )?//number of traces for rank k in group g
	      ( tmp/size + 1 ):
	      ( tmp/size );
	    DATA_SZ_y[ r + g*size ] *= NS_y;
	  }
	}
	DATA_OFF_y = new int[SIZE];
	for( int r=0; r<SIZE; r++ ){
	  DATA_OFF_y[r]=0;
	  for( int k=0; k<r; k++ ){
	    DATA_OFF_y[r] += DATA_SZ_y[k];
	  }
	}

	MPI_Gatherv( y_arr_loc, NS_y*NTR_y_loc, MPI_FLOAT,
		     y_arr, DATA_SZ_y, DATA_OFF_y, MPI_FLOAT,
		     0, COMM );

#ifdef VERBOSE_MPS_CONV_MPI_CC
	for( int r=0; r<SIZE; r++ ){
	  if(RANK==r){
	    cerr <<"rank="<<r<<", Post gatherv call\n"
		 <<"  sendcount = NS_y*NTR_y_loc ="<<NS_y*NTR_y_loc<<"\n";
	  }
	  MPI_Barrier(COMM);
	}

	if(RANK==0){
	  cerr <<"printing out DATA_SZ_y:\n";
	  for(int i=0; i<SIZE; i++)
	    cerr <<"  DATA_SZ_y["<<i<<"]="<<DATA_SZ_y[i]<<"\n";
       
	  cerr <<"printing out DATA_OFF_y:\n";
	  for(int i=0; i<SIZE; i++)
	    cerr <<"  DATA_OFF_y["<<i<<"]="<<DATA_OFF_y[i]<<"\n";	 
	}
	MPI_Barrier(COMM);
#endif

	//////////////////////
	// writing data out //
	//////////////////////
	if( RANK==0 ){

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
	}

#ifdef VERBOSE_MPS_CONV_MPI_CC
	MPI_Barrier(COMM);
	for( int r=0; r<SIZE; r++ ){
	  if(RANK==r){
	    cerr <<"rank="<<r<<", Post writing data out\n";
	  }
	  MPI_Barrier(COMM);
	}
#endif
    
	if(y_arr!=NULL) delete[] y_arr;
	if(k_arr!=NULL) delete[] k_arr;	
	if(y_arr_loc!=NULL) delete[] y_arr_loc;
	if(k_arr_loc!=NULL) delete[] k_arr_loc;
	if(DATA_SZ_k!=NULL) delete[] DATA_SZ_k;
	if(DATA_OFF_k!=NULL) delete[] DATA_OFF_k;
	if(DATA_SZ_y!=NULL) delete[] DATA_SZ_y;
	if(DATA_OFF_y!=NULL) delete[] DATA_OFF_y;

      }//component loop
      MPI_Barrier(COMM);
      

      if(x_arr!=NULL) delete[] x_arr;
      if(x_arr_loc!=NULL) delete[] x_arr_loc;
      if(NSPG!=NULL) delete[] NSPG;

#ifdef VERBOSE_MPS_CONV_MPI_CC
      for( int r=0; r<SIZE; r++ ){
	if(RANK==r){
	  cerr <<"rank="<<r<<", On my way out of kernel!\n";
	}
	MPI_Barrier(COMM);
      }
#endif

    }
    catch(RVLException &e){
      e << "ERROR from MPS_conv_kern!\n";
      throw e;
    }
  }


  //----------------------------------------------------------------------// 
  void MPS_corr_kern_MPI( vector<Rtuple> s_pos, 
			  int MPS_dim, 
			  int MPS_idx,
			  string x_file,
			  vector<string> y_files,
			  vector<string> k_files ){
  //----------------------------------------------------------------------// 
    try{

#ifdef VERBOSE_MPS_CONV_MPI_CC
      MPI_Barrier(retrieveGlobalComm());
      if( retrieveGlobalRank()==0 ){
	cerr << "*************************************************\n"
	     << "  Inside MPS_corr_kern_MPI\n"
	     << "*************************************************\n\n";
	
	cerr << "x file = " << x_file << "\n"
	     << "y files = "; print_vec(y_files);
	cerr << "k files = "; print_vec(k_files);
	cerr << "MPS_dim = " << MPS_dim << "\n"
	     << "MPS_idx = " << MPS_idx << "\n";
      }
      MPI_Barrier(retrieveGlobalComm());
#endif
      
      FILE *fp_x, *fp_y, *fp_k;
      segy tr_x, tr_y, tr_k;
      
      //MPI variables
      int RANK = retrieveGlobalRank();
      int rank = retrieveRank();
      int SIZE = retrieveGlobalSize(); 
      int size = retrieveSize();       
      int N_gr = retrieveNumGroups();  
      int grID = retrieveGroupID(); 
      MPI_Comm COMM = retrieveGlobalComm();
      MPI_Comm comm = retrieveComm();
      
      int N_src = s_pos.size(); //total number of sources
      int N_src_loc; //number of sources at current process
      
      int NS_x,      NS_y,      NS_k; //total number of time sample points per trace
      int NTR_x,     NTR_y,     NTR_k; //total number of traces
      int NTR_x_loc, NTR_y_loc, NTR_k_loc; //total number of traces at current process
      int NTPS; //number of traces per source

      float dt_x,    dt_y,    dt_k;
      float delrt_x, delrt_y, delrt_k;

      //global arrays
      float * x_arr = NULL; 
      float * y_arr = NULL;
      float * k_arr = NULL;     
      //local arrays
      float * x_arr_loc = NULL;
      float * y_arr_loc = NULL;
      float * k_arr_loc = NULL;
      float * x_arr_red = NULL; //buffer for group root reduction
      float * x_buff = NULL; //single trace buffer
      //other arrays
      int *NSPG = NULL; //array containing number of sources per group, for each group
      int *DATA_SZ_k  = NULL; //array containing size of data related to y_arr to be sent via scatterv
      int *DATA_OFF_k = NULL; //array of offsets for data to be sent via scatterv
      int *DATA_SZ_y   = NULL; //array containing size of data related to kernel to be sent via scatterv
      int *DATA_OFF_y  = NULL; //array of offsets for data to be sent via scatterv
      fpos_t pos_x;

      if( MPS_idx<0 || MPS_idx>=MPS_dim ){
	RVLException e;
	e << "MPS_idx="<< MPS_idx << " out of bounds!\n";
	throw e;
      }

      //setting up number of sources per group (NSPG) array
      NSPG = new int[N_gr];
      for( int g=0; g<N_gr; g++ ){
	NSPG[g] = ( g<(N_src%N_gr) )?
	  ( N_src/N_gr + 1 ):
	  ( N_src/N_gr );
      }
      N_src_loc = NSPG[grID];


      //////////////////////
      // setting up x_arr //
      //////////////////////
      NTR_x = get_ntr(x_file); 

      if(MPS_dim!=(NTR_x/N_src)){
	RVLException e;
	e << "MPS_dim="<<MPS_dim<<" does not match number of traces per source,  NTR_x/N_src="
	  << NTR_x/N_src<<"\n";
	throw e;
      }

      if( RANK==0 ){

	//getting size
	fp_x = iwave_const_fopen(x_file.c_str(),"r",NULL,stderr);
	if( fp_x==NULL ){
	  RVLException e;
	  e << "Could not open input file "<< x_file <<"\n";
	  throw e;
	}
	fgettr(fp_x,&tr_x);
	iwave_fclose(fp_x); fp_x=NULL;

	NS_x = tr_x.ns; 
	dt_x = tr_x.dt;
	delrt_x = tr_x.delrt;
	
	x_arr = new float[NS_x*N_src];
	if(x_arr==NULL){
	  RVLException e;
	  e << "Trouble allocating memory for input buffer array.\n";
	  throw e;
	}
      
	iwave_fclose(fp_x); fp_x=NULL;

      }//end RANK==0 if

      //broadcasting 
      MPI_Bcast( &NS_x, 1, MPI_INT, 0, COMM );
      MPI_Bcast( &dt_x, 1, MPI_FLOAT, 0, COMM );
      MPI_Bcast( &delrt_x, 1, MPI_FLOAT, 0, COMM );
      
#ifdef VERBOSE_MPS_CONV_MPI_CC
      for( int r=0; r<SIZE; r++){
	if(RANK==r){
	  cerr <<"rank="<<r<<", After setting up x_arr.\n"
	       <<"  NTR_x   = "<<NTR_x<<"\n"
	       <<"  NS_x    = "<<NS_x<<"\n"
	       <<"  dt_x    = "<<dt_x<<"\n"
	       <<"  delrt_x = "<<delrt_x<<"\n";
	}
	MPI_Barrier(COMM);
      }
#endif


      /////////////////////////////////////////////////
      // setting up x_arr_loc, x_buff, and x_arr_red //
      /////////////////////////////////////////////////
      NTR_x_loc = N_src_loc;
      x_arr_loc = new float[NS_x*NTR_x_loc];
      for( int i=0; i<NS_x*NTR_x_loc; i++)
	x_arr_loc[i] = 0.0f;

      x_buff = new float[NS_x];

      if( rank==0 )
	x_arr_red = new float[NS_x*NTR_x_loc];



      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      // looping over RHS components //
      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      for(int i_c=0; i_c<y_files.size(); i_c++){
	
#ifdef VERBOSE_MPS_CONV_MPI_CC
	for( int r=0; r<SIZE; r++){
	  if(RANK==r){
	    cerr << "rank="<<r<<", Inside components-loop: i_c="<<i_c<<"\n";	
	  }
	  MPI_Barrier(COMM);
	}
#endif
      
	//////////////////////
	// setting up k_arr //
	//////////////////////
	NTR_k = get_ntr(k_files[i_c]);
	NTPS = NTR_k/N_src;

	if( RANK==0 ){
	  
	  //getting size
	  fp_k = iwave_const_fopen(k_files[i_c].c_str(),"r",NULL,stderr);
	  if( fp_k==NULL ){
	    RVLException e;
	    e << "Could not open kernel file "<< k_files[i_c] <<"\n";
	    throw e;
	  }
	  fgettr(fp_k,&tr_k);
	  iwave_fclose(fp_k); fp_k=NULL;
	  
	  NS_k = tr_k.ns;
	  dt_k = tr_k.dt;
	  delrt_k = tr_k.delrt;
	  
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
	
	}//end RANK==0 if

	//broadcasting
	MPI_Bcast( &NS_k, 1, MPI_INT, 0, COMM );
	MPI_Bcast( &dt_k, 1, MPI_FLOAT, 0, COMM );
	MPI_Bcast( &delrt_k, 1, MPI_FLOAT, 0, COMM );

#ifdef VERBOSE_MPS_CONV_MPI_CC
	for( int r=0; r<SIZE; r++){
	  if(RANK==r){
	    cerr <<"rank="<<r<<", After setting up k_arr.\n"
		 <<"  NTPS    = "<<NTPS<<"\n"
		 <<"  NTR_k   = "<<NTR_k<<"\n"
		 <<"  NS_k    = "<<NS_k<<"\n"
		 <<"  dt_k    = "<<dt_k<<"\n"
		 <<"  delrt_k = "<<delrt_k<<"\n";
	  }
	  MPI_Barrier(COMM);
	}
#endif

	////////////////////////
	// distributing k_arr //
	////////////////////////
	NTR_k_loc = N_src_loc * NTPS;
	k_arr_loc = new float[NS_k*NTR_k_loc];

	DATA_SZ_k = new int[SIZE];
	for( int g=0; g<N_gr; g++ ){
	  for( int r=0; r<size; r++ ){
	    int tmp = NSPG[g]*NTPS;
	    DATA_SZ_k[ r + g*size ] = ( r<(tmp%size) )?
	      ( tmp/size + 1 ):
	      ( tmp/size );
	    DATA_SZ_k[ r + g*size ] *= NS_k;
	  }
	}
	DATA_OFF_k = new int[SIZE];
	for( int r=0; r<SIZE; r++ ){
	  DATA_OFF_k[r]=0;
	  for( int k=0; k<r; k++ ){
	    DATA_OFF_k[r] += DATA_SZ_k[k];
	  }
	}

	//cerr << "----> Scatter 1\n";
	MPI_Scatterv( k_arr, DATA_SZ_k, DATA_OFF_k, MPI_FLOAT,
		      k_arr_loc, NS_k*NTR_k_loc, MPI_FLOAT,
		      0, COMM );
      
#ifdef VERBOSE_MPS_CONV_MPI_CC
	for( int r=0; r<SIZE; r++){
	  if(RANK==r){
	    cerr <<"rank="<<r<<", After distributing (scatterv) k_arr.\n"
		 <<"  NTR_k_loc = "<<NTR_k_loc<<"\n"
		 <<"  recvcount = NS_k*NTR_k_loc = "<<NS_k*NTR_k_loc<<"\n";
	  }
	  MPI_Barrier(COMM);
	}
	
	if(RANK==0){
	  cerr <<"printing out DATA_SZ_k:\n";
	  for(int i=0; i<SIZE; i++)
	    cerr <<"  DATA_SZ_k["<<i<<"]="<<DATA_SZ_k[i]<<"\n";
       
	  cerr <<"printing out DATA_OFF_k:\n";
	  for(int i=0; i<SIZE; i++)
	    cerr <<"  DATA_OFF_k["<<i<<"]="<<DATA_OFF_k[i]<<"\n";	 
	}
	MPI_Barrier(COMM);
#endif


	//////////////////////
	// setting up y_arr //
	//////////////////////
	NTR_y = get_ntr(y_files[i_c]);

	if( RANK==0 ){

	  //getting size
	  fp_y = iwave_const_fopen(y_files[i_c].c_str(),"r",NULL,stderr);
	  if( fp_y==NULL ){
	    RVLException e;
	    e << "Could not open output file "<< y_files[i_c] <<"\n";
	    throw e;
	  }
	  fgettr(fp_y,&tr_y);
	  iwave_fclose(fp_y); fp_y=NULL;
	  
	  NS_y = tr_y.ns;
	  dt_y = tr_y.dt;
	  delrt_y = tr_y.delrt;
	  
	  y_arr = new float[NTR_y*NS_y];
	  if(y_arr==NULL){
	    RVLException e;
	    e << "Trouble allocating memory for output buffer array.\n";
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

	}//end RANK==0 if
      
	//broadcasting 
	MPI_Bcast( &NS_y, 1, MPI_INT, 0, COMM );
	MPI_Bcast( &dt_y, 1, MPI_FLOAT, 0, COMM );
	MPI_Bcast( &delrt_y, 1, MPI_FLOAT, 0, COMM );

#ifdef VERBOSE_MPS_CONV_MPI_CC
	for( int r=0; r<SIZE; r++){
	  if(RANK==r){
	    cerr <<"rank="<<r<<", After setting up y_arr.\n"
		 <<"  NTR_y   = "<<NTR_y<<"\n"
		 <<"  NS_y    = "<<NS_y<<"\n"
		 <<"  dt_y    = "<<dt_y<<"\n"
		 <<"  delrt_y = "<<delrt_y<<"\n";
	  }
	  MPI_Barrier(COMM);
	}
#endif

	////////////////////////
	// distributing y_arr //
	////////////////////////
	NTR_y_loc = N_src_loc * NTPS;
	y_arr_loc = new float[NS_y*NTR_y_loc];


	DATA_SZ_y = new int[SIZE];
	for( int g=0; g<N_gr; g++ ){
	  for( int r=0; r<size; r++ ){
	    int tmp = NSPG[g]*NTPS; //number of traces for group g	    
	    DATA_SZ_y[ r + g*size ] = ( r<(tmp%size) )?//number of traces for rank k in group g
	      ( tmp/size + 1 ):
	      ( tmp/size );
	    DATA_SZ_y[ r + g*size ] *= NS_y;
	  }
	}
	DATA_OFF_y = new int[SIZE];
	for( int r=0; r<SIZE; r++ ){
	  DATA_OFF_y[r]=0;
	  for( int k=0; k<r; k++ ){
	    DATA_OFF_y[r] += DATA_SZ_y[k];
	  }
	}

	//cerr <<"----> Scatter 2\n";
	MPI_Scatterv( y_arr, DATA_SZ_y, DATA_OFF_y, MPI_FLOAT,
		      y_arr_loc, NS_y*NTR_y_loc, MPI_FLOAT,
		      0, COMM );
      
#ifdef VERBOSE_MPS_CONV_MPI_CC
	for( int r=0; r<SIZE; r++ ){
	  if(RANK==r){
	    cerr <<"rank="<<r<<", Post scatterv call\n"
		 <<"  recvcount = NS_y*NTR_y_loc ="<<NS_y*NTR_y_loc<<"\n";
	  }
	  MPI_Barrier(COMM);
	}

	if(RANK==0){
	  cerr <<"printing out DATA_SZ_y:\n";
	  for(int i=0; i<SIZE; i++)
	    cerr <<"  DATA_SZ_y["<<i<<"]="<<DATA_SZ_y[i]<<"\n";
       
	  cerr <<"printing out DATA_OFF_y:\n";
	  for(int i=0; i<SIZE; i++)
	    cerr <<"  DATA_OFF_y["<<i<<"]="<<DATA_OFF_y[i]<<"\n";	 
	}
	MPI_Barrier(COMM);
#endif

	//checking dt's
	if( dt_x!=dt_y || dt_x!=dt_k ){
	  RVLException e;
	  e << "Trace dt's are not the same for input, output, or kernel for convolution!\n"
	    << "component idx="<< i_c <<"\n";
	  throw e;
	}
	if( NTR_k!=NTR_y ){
	  RVLException e;
	  e << "NTR is not the same between kernel and output SUs!\n";
	  throw e;
	}


	//////////////////////////
	// args for correlation //
	//////////////////////////
	float dt  = dt_x;
	int n_in   = NS_y;
	int i0_in  = int( 1000 * tr_y.delrt/dt );
	int n_out  = NS_x;
	int i0_out = int( 1000 * tr_x.delrt/dt );
	int n_ker  = NS_k;
	int i0_ker = int( 1000 * tr_k.delrt/dt );
	float scal = dt/1e6;

#ifdef VERBOSE_MPS_CONV_MPI_CC
	if( RANK==0 ){
	  cerr << "\nCorrelation args:\n"
	       << "    n_in   ="<< n_in   <<"\n"
	       << "    i0_in  ="<< i0_in  <<"\n"
	       << "    n_out  ="<< n_out  <<"\n"
	       << "    i0_out ="<< i0_out <<"\n"
	       << "    n_ker  ="<< n_ker  <<"\n"
	       << "    i0_ker ="<< i0_ker <<"\n"
	       << "    tr_y.delrt/dt="<<tr_y.delrt/dt<<"\n"
	       << "    tr_x.delrt/dt="<<tr_x.delrt/dt<<"\n"
	       << "    tr_k.delrt/dt="<<tr_k.delrt/dt<<"\n"
	       << "    scal = " << scal <<"\n";
	}
	MPI_Barrier(COMM);
#endif


	//~~~~~~~~~~~~~~~~~~~~~~//
	// looping over traces  //
	//~~~~~~~~~~~~~~~~~~~~~~//
	for( int i_tr=0; i_tr<NTR_k_loc; i_tr++ ){

	  //current source index, local to group
	  int i_src_loc = ( i_tr + DATA_OFF_k[RANK]/NS_k )/NTPS;
	  for( int g=0; g<grID; g++ ) i_src_loc -= NSPG[g];	  
	  
	  //actual convolution
	  int shift = i0_in - i0_out - i0_ker;
	  TSOpt::corr( shift,
		       n_out,
		       n_in,
		       n_ker,
		       x_buff,
		       y_arr_loc+i_tr*NS_y,
		       k_arr_loc+i_tr*NS_k,
		       scal );
	  
	  //accumulating x_arr_loc
	  for(int i=0; i<NS_x; i++){
	    x_arr_loc[i+i_src_loc*NS_x] += x_buff[i];
	  }

	}//trace loop
	MPI_Barrier(COMM);
	
	if(y_arr!=NULL) delete[] y_arr;
	if(k_arr!=NULL) delete[] k_arr;	
	if(y_arr_loc!=NULL) delete[] y_arr_loc;
	if(k_arr_loc!=NULL) delete[] k_arr_loc;
	if(DATA_SZ_k!=NULL) delete[] DATA_SZ_k;
	if(DATA_OFF_k!=NULL) delete[] DATA_OFF_k;
	if(DATA_SZ_y!=NULL) delete[] DATA_SZ_y;
	if(DATA_OFF_y!=NULL) delete[] DATA_OFF_y;

      }//component loop
      MPI_Barrier(COMM);

      
      //////////////////////////////
      // reduction to group roots //
      //////////////////////////////
      MPI_Reduce( x_arr_loc,
		  x_arr_red,
		  N_src_loc*NS_x,
		  MPI_FLOAT,
		  MPI_SUM,
		  0,
		  comm );
      
      //////////////////////////////////////
      // sending reduction to global root //
      //////////////////////////////////////
      int tag=0;
      if( rank==0 && RANK!=0 ){
	MPI_Send( x_arr_red,
		  N_src_loc*NS_x,
		  MPI_FLOAT,
		  0,
		  tag,
		  COMM );
      }
      if( RANK==0 ){

	int offset=0;
	int from;
	for(int g=0; g<N_gr; g++){

	  from = g*size;
	  if(g==0){
	    for(int i=0; i<N_src_loc*NS_x; i++){
	      x_arr[i] = x_arr_red[i];
	    }
	  }
	  else{
	    MPI_Recv( x_arr+offset,
		      NSPG[g]*NS_x,
		      MPI_FLOAT,
		      from,
		      tag,
		      COMM,
		      MPI_STATUS_IGNORE );
	  }	 
	  offset+=NSPG[g]*NS_x;
	}
      }
      MPI_Barrier(COMM);


      ///////////////////////////////////
      // writing output buffer to file //
      ///////////////////////////////////
      if( RANK==0 ){

	/*	
	cerr << "printing out x_arr\n";
	for(int i=0; i<NS_x*NTR_x; i++){
	  cerr << "  x_arr["<<i<<"]="<<x_arr[i]<<"\n";
	}
	*/

	int offset;
	int tr_shift;

	fp_x = iwave_const_fopen(x_file.c_str(),"r+",NULL,stderr);
	for( int i_src=0; i_src<N_src; i_src++ ){
	  
	  tr_shift = i_src*MPS_dim + MPS_idx;
	  offset = i_src*NS_x;	

	  /*	  
	  cerr << "i_src = "<< i_src<<"\n"
	       << "  tr_shift="<<tr_shift<<"\n"
	       << "  offset="<<offset<<"\n";
	  */

	  if(tr_shift!=0) fgettra(fp_x,&tr_x,tr_shift-1);
	  fgetpos(fp_x,&pos_x);
	  fgettr(fp_x,&tr_x);
	  
	  for( int i_t=0; i_t<NS_x; i_t++ ){
	    tr_x.data[i_t] = x_arr[ i_t+offset ];
	  }
	  fsetpos(fp_x,&pos_x);
	  fputtr(fp_x,&tr_x);
	}
	fflush(fp_x);
	iwave_fclose(fp_x); fp_x=NULL;
      }


      if(x_arr!=NULL) delete[] x_arr;
      if(x_arr_loc!=NULL) delete[] x_arr_loc;
      if(x_arr_red!=NULL) delete[] x_arr_red;
      if(x_buff!=NULL) delete[] x_buff;
      if(NSPG!=NULL) delete[] NSPG;

#ifdef VERBOSE_MPS_CONV_MPI_CC
      for( int r=0; r<SIZE; r++ ){
	if(RANK==r){
	  cerr <<"rank="<<r<<", On my way out of kernel!\n";
	}
	MPI_Barrier(COMM);
      }
#endif
    }
    catch(RVLException &e){
      e << "ERROR from MPS_conv_kern!\n";
      throw e;
    }
  }



}

#endif
