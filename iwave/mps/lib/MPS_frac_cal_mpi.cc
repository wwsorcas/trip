// MPS_frac_cal_mpi.cc
// Author: Mario J. Bencomo
// last modified: 11/10/16

#include "MPS_frac_cal.hh"

#ifdef IWAVE_USE_MPI
//#define VERBOSE_MJB

namespace TSOpt{

  //----------------------------------------------------------------------// 
  void MPS_frac_cal_kern_MPI( vector<float> order,
			      vector<float> scalars,
			      int N_src,
			      string in_file,
			      string out_file,
			      bool adj,
			      bool inv ){
  //----------------------------------------------------------------------// 
    try{
      
#ifdef VERBOSE_MJB
      if( retrieveGlobalRank()==0 ){
	cerr << "*************************************************\n"
	     << "  Inside MPS_frac_cal_kern_MPI\n"
	     << "*************************************************\n\n"
	     << "in  file = " << in_file << "\n"
	     << "out file = " << out_file << "\n"
	     << "adj flag = " << adj << "\n"
	     << "inv flag = " << inv << "\n"
	     << "order:";
	for(int i=0; i<order.size(); i++) cerr << order[i] << ",";
	cerr <<"\n";
	cerr <<"scalars:";
	for(int i=0; i<order.size(); i++) cerr << order[i] << ",";
	cerr <<"\n";
      }
      MPI_Barrier(retrieveGlobalComm());
#endif

      FILE *fp_x, *fp_y;
      segy tr;
      int NTR;
      int MPS_dim;
      int NS;
      float dt;
      float q;
      float scal;
      float *x_arr = NULL;
      float *y_arr = NULL;
      float *x_arr_loc = NULL;
      float *y_arr_loc = NULL;

      //MPI variables
      int RANK = retrieveGlobalRank();
      //int rank = retrieveRank();
      int SIZE = retrieveGlobalSize();
      //int size = retrieveSize();
      int N_gr = retrieveNumGroups();
      //int grID = retrieveGroupID();
      MPI_Comm COMM = retrieveGlobalComm();
      //MPI_Comm comm = retrieveComm();

      int N_src_loc;
      int *snd_sz = NULL;
      int *offsets = NULL;
      fpos_t pos;

      //setting NTR and MPS_dim
      NTR = get_ntr(in_file);
      MPS_dim = order.size();
      if( NTR!=(N_src*MPS_dim) ){
	RVLException e;
	e << "Number of traces is wrong for input file: "<< in_file<<"\n"
	  << "Expected NTR = "<< N_src*MPS_dim <<"\n"
	  << "Actual NTR = "<< NTR <<"\n";
	throw e;
      }

      //////////////////////
      // Setting up x_arr //
      //////////////////////
      if( RANK==0 ){
	//getting NS and dt
	fp_x = iwave_const_fopen(in_file.c_str(),"r",NULL,stderr);
	if( fp_x==NULL ){
	  RVLException e;
	  e << "Could not open input file "<< in_file<<".\n";
	  throw e;
	}      
	fgettr(fp_x,&tr);
	iwave_fclose(fp_x); fp_x=NULL;
	NS = tr.ns;
	dt = tr.dt*1e-6; //converting into units [sec]
      
	//allocating x_arr
	x_arr = new float[NTR*NS];
	fp_x = iwave_const_fopen(in_file.c_str(),"r",NULL,stderr);
	for( int i=0; i<NTR; i++){
	  fgettr(fp_x,&tr);
	  for( int i_t=0; i_t<NS; i_t++ ){
	    x_arr[i_t+i*NS] = tr.data[i_t]; 
	  }
	}
      }

      //broadcasting 
      MPI_Bcast( &NS, 1, MPI_INT, 0, COMM );
      MPI_Bcast( &dt, 1, MPI_FLOAT, 0, COMM );


      ////////////////////////
      // distributing x_arr //
      ////////////////////////
      N_src_loc = ( RANK<(N_src%N_gr) )?
	(N_src/N_gr+1):
	(N_src/N_gr);
      snd_sz = new int[SIZE];
      for( int i=0; i<SIZE; i++ ){
	snd_sz[i]=0;
	if( i<N_gr ){
	  snd_sz[i] = (i<(N_src%N_gr))?
	    (N_src/N_gr+1):
	    (N_src/N_gr);
	  snd_sz[i] *= NS*MPS_dim;
	}
      }
      offsets = new int[SIZE];
      for( int i=0; i<SIZE; i++ ){
	offsets[i]=0;
	if( i<N_gr ){
	  for( int j=0; j<i; j++ ){
	    offsets[i] += snd_sz[j]; 
	  }
	}
      }
      
#ifdef VERBOSE_MJB
      for( int r=0; r<SIZE; r++ ){
	if(RANK==r){
	  cerr <<"rank="<<r<<"\n";
	  for(int i=0; i<SIZE; i++){
	    cerr <<"  snd_sz["<<i<<"]="<<snd_sz[i]<<" "
		 <<"  offsets["<<i<<"]="<<offsets[i]<<"\n";
	  }
	}
	MPI_Barrier(COMM);
      }
#endif

      x_arr_loc = new float[N_src_loc*MPS_dim*NS];//snd_sz[RANK]];

      MPI_Scatterv( x_arr,snd_sz,offsets,MPI_FLOAT,
		    x_arr_loc,snd_sz[RANK],MPI_FLOAT,
		    0, COMM );

      /////////////////////////////////
      // Setting y_arr and y_arr_loc //
      /////////////////////////////////
      if( RANK==0 ){
	y_arr = new float[NTR*NS];
      }
      y_arr_loc = new float[N_src_loc*MPS_dim*NS];//snd_sz[RANK]];
      

#ifdef VERBOSE_MJB
      for( int r=0; r<SIZE; r++ ){
	if(RANK==r){
	  cerr <<"rank="<<r<<", N_src_loc="<< N_src_loc <<"!\n";
	}
	MPI_Barrier(COMM);
      }
#endif


      if( RANK<N_gr ){
	int offset;
	////////////////////////
	//looping over sources//
	////////////////////////
	for(int i_src=0; i_src<N_src_loc; i_src++){
	  
	  //////////////////////////
	  //looping over MPS basis//
	  //////////////////////////
	  for(int i_b=0; i_b<MPS_dim; i_b++){
	    
	    q = order[i_b];
	    offset = (i_src+i_b)*NS;
	    
	    q *= (inv)?(-1):(1);
	    scal = (inv)?(1.0f/scalars[i_b]):(scalars[i_b]); 

	    if(q>=0){
	      frac_deriv( q, scal, NS,
			  x_arr_loc+offset,
			  y_arr_loc+offset, 
			  dt, adj );
	    }
	    else{
	      q *= -1;
	      frac_integ( q, scal, NS,
			  x_arr_loc+offset,
			  y_arr_loc+offset,
			  dt, adj );
	    }
	    
	  }//MPS basis loop	
	}//sources loop
      }

      ////////////////////
      // gathering data //
      ////////////////////
      MPI_Gatherv( y_arr_loc, snd_sz[RANK], MPI_FLOAT,
		   y_arr, snd_sz, offsets, MPI_FLOAT,
		   0, COMM );
      
#ifdef VERBOSE_MJB
      for( int r=0; r<SIZE; r++ ){
	if(RANK==r){
	  cerr <<"rank="<<r<<", Post gatherv call\n";
	}
	MPI_Barrier(COMM);
      }
#endif
      
      /////////////////
      // Writing out //
      /////////////////
      if( RANK==0 ){

	fp_y= iwave_const_fopen(out_file.c_str(),"r+",NULL,stderr);     
	for( int i=0; i<NTR; i++ ){
	  fgetpos(fp_y,&pos);
	  fgettr(fp_y,&tr);
	  for( int i_t=0; i_t<NS; i_t++ ){
	    tr.data[i_t] = y_arr[ i_t+i*NS ];
	  }
	  fsetpos(fp_y,&pos);
	  fputtr(fp_y,&tr);
	}
	fflush(fp_y);
	iwave_fclose(fp_y); fp_y=NULL;
      }

#ifdef VERBOSE_MJB
      for( int r=0; r<SIZE; r++ ){
	if(RANK==r){
	  cerr <<"rank="<<r<<", Post writing data out\n";
	}
	MPI_Barrier(COMM);
      }
#endif

      if(x_arr!=NULL) delete[] x_arr;
      if(y_arr!=NULL) delete[] y_arr;      
      if(x_arr_loc!=NULL) delete[] x_arr_loc;
      if(y_arr_loc!=NULL) delete[] y_arr_loc;
      if(snd_sz!=NULL)  delete[] snd_sz;
      if(offsets!=NULL) delete[] offsets;

#ifdef VERBOSE_MJB
      for( int r=0; r<SIZE; r++ ){
	if(RANK==r){
	  cerr <<"rank="<<r<<", On my way out of kernel!\n";
	}
	MPI_Barrier(COMM);
      }
#endif

    }
    catch(RVLException &e){
      e << "ERROR from MPS_frac_cal_kern_MPI!\n";
      throw e;
    }
  }


}//end TSOpt

#endif
