// MPS_iwop.cc
// Author: Mario J. Bencomo
// last modified: 11/10/16

#include "MPS_iwop.hh"

//#define VERBOSE_MJB

namespace TSOpt {
  
  //----------------------------------------------------------------------//
  void zero_out_ex(string filename, 
		   int idx,
		   int Nsrc){
  //----------------------------------------------------------------------//
    try{
#ifdef IWAVE_USE_MPI
      if(retrieveGlobalRank()!=0){
	RVLException e;
	e << "ERROR! zero_out_ex() should only be called from root process.\n";
	throw e;
      }
#endif  

#ifdef VERBOSE_MJB    
      cerr << ">>>> Inside zero_out_ex <<<<\n"
	   << "  filename="<<filename<<"\n"
	   << "  idx="<<idx<<"\n"
	   << "  Nsrc="<<Nsrc<<"\n";
#endif

      FILE *fp = NULL;
      segy tr;
      fpos_t pos;
      //fpos_t pos_d;

      //getting total number of traces
      int ntr = get_ntr_serial(filename);
      int ntr_nsrc = ntr/Nsrc;
      if(idx>=(ntr_nsrc)){
	RVLException e;
	e << "Index idx="<<idx<<" is greater than ntr/Nsrc="<<ntr/Nsrc<<"\n"
	  << "ntr="<<ntr <<", Nsrc="<<Nsrc<<"\n";
	throw e;
      }
      if(ntr%Nsrc!=0){
	RVLException e;
	e << "Number of traces is not divisible by number of sources!\n";
	throw e;
      }

#ifdef VERBOSE_MJB
      cerr << "ntr="<<ntr<<"\n"
	   << "ntr/Nsrc="<<ntr_nsrc<<"\n";
#endif

      //opening file 
      if( !(fp=iwave_const_fopen(filename.c_str(),"r+",NULL,stderr)) ){
	RVLException e;
	e << "failed to open file: "<< filename <<"\n";
	throw e;
      }
      
      /*
      //making a trace of zeroes
      fgetpos(fp,&pos);
      fgettr(fp,&tr);
      for( int it=0; it<tr.ns; it++)
	tr.data[it]=0;
      
      //getting pos_d
      fgetpos(fp,&pos_d);
      pos_d = pos_d - pos;
      */

      //zeroing out all traces except trace idx
      for( int i=0; i<ntr; i++ ){

	if( idx==(i%ntr_nsrc) ) {
	  //pos = pos + pos_d;
	  continue;
	}
	
	if(i!=0) fgettra(fp,&tr,i-1);
	fgetpos(fp,&pos);
	fgettr(fp,&tr);

	//zeroing out trace
	for( int i_t=0; i_t<tr.ns; i_t++){
	  tr.data[i_t] = 0.0;
	}
	
	//writting out trace
	fsetpos(fp,&pos);
	fputtr(fp,&tr);
      }

      iwave_fclose(fp);
    }
    catch(RVLException &e){
      e << "ERROR from zero_out_ex\n";
      throw e;
    }
  }
}
