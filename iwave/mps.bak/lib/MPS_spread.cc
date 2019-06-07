// MPS_spread.cc
// Author: Mario J. Bencomo
// last modified: 02/20/17

#include "MPS_spread.hh"

//#define VERBOSE_MJB

namespace TSOpt{

  //----------------------------------------------------------------------// 
  void MPS_spread_kern( int N_src, 
			int MPS_dim, 
			string x_file,
			string y_file ){
  //----------------------------------------------------------------------// 
    try{
      
#ifdef VERBOSE_MJB
      cerr << "*************************************************\n"
           << "  Inside MPS_spread_kern\n"
           << "*************************************************\n"
	   << "x file  = " << x_file << "\n"
	   << "y file  = " << y_file << "\n"
	   << "MPS_dim = " << MPS_dim << "\n"
	   << "N_src   = " << N_src << "\n";
#endif
      
      FILE *fp_x, *fp_y;
      segy tr_x, tr_y;
      fpos_t pos_y; 

      //opening input and output files
      fp_x = iwave_const_fopen(x_file.c_str(),"r",NULL,stderr);
      if( fp_x==NULL ){
	RVLException e;
	e << "Could not open input file "<< x_file <<"\n";
	throw e;
      }

      fp_y = iwave_const_fopen(y_file.c_str(),"r+",NULL,stderr);
      if( fp_x==NULL ){
	RVLException e;
	e << "Could not open output file "<< y_file <<"\n";
	throw e;
      }

      //checking time axes are the same
      fgettr(fp_x,&tr_x);
      fgettr(fp_y,&tr_y);

      if( tr_x.ns!=tr_y.ns || 
	  tr_x.dt!=tr_y.dt ||
	  tr_x.delrt!=tr_y.delrt ){
	RVLException e;
	e << "Time axes do not match!\n"
	  << "tr_x.ns="<<tr_x.ns<<", tr_x.dt="<<tr_x.dt<<", tr_x.delrt="<<tr_x.delrt<<"\n"
	  << "tr_y.ns="<<tr_y.ns<<", tr_y.dt="<<tr_y.dt<<", tr_y.delrt="<<tr_y.delrt<<"\n";
	throw e;
      }

      fseeko(fp_x,0L,SEEK_SET);
      fseeko(fp_y,0L,SEEK_SET);

      //applying MPS spread operation
      for( int i_src=0; i_src<N_src; i_src++ ){
	fseeko(fp_x,0L,SEEK_SET);
	for( int i_b=0; i_b<MPS_dim; i_b++ ){
	  fgettr(fp_x,&tr_x);
	  fgetpos(fp_y,&pos_y);
	  fgettr(fp_y,&tr_y);
	  
	  for( int i_t=0; i_t<tr_x.ns; i_t++ ){
	    tr_y.data[i_t] = tr_x.data[i_t];
	  }
	  fsetpos(fp_y,&pos_y);
	  fputtr(fp_y,&tr_y);
	}
      }
      fflush(fp_y);
      iwave_fclose(fp_y);
      iwave_fclose(fp_x);

#ifdef VERBOSE_MJB
      cerr << "Leaving MPS_spread_kern\n";
#endif

    }
    catch(RVLException &e){
      e << "ERROR from MPS_spread_kern!\n";
      throw e;
    }
  }

  
  //----------------------------------------------------------------------// 
  void MPS_stack_kern( int N_src, 
		       int MPS_dim, 
		       string x_file,
		       string y_file ){
  //----------------------------------------------------------------------// 
    try{
      
#ifdef VERBOSE_MJB
      cerr << "*************************************************\n"
           << "  Inside MPS_stack_kern\n"
           << "*************************************************\n"
	   << "x file  = " << x_file << "\n"
	   << "y file  = " << y_file << "\n"
	   << "MPS_dim = " << MPS_dim << "\n"
	   << "N_src   = " << N_src << "\n";
#endif
      
      FILE *fp_x, *fp_y;
      segy tr_x, tr_y;
      fpos_t pos_x; 
      float *buff;

      //opening files
      fp_x = iwave_const_fopen(x_file.c_str(),"r+",NULL,stderr);
      if( fp_x==NULL ){
	RVLException e;
	e << "Could not open output file "<< x_file <<"\n";
	throw e;
      }

      fp_y = iwave_const_fopen(y_file.c_str(),"r",NULL,stderr);
      if( fp_x==NULL ){
	RVLException e;
	e << "Could not open input file "<< y_file <<"\n";
	throw e;
      }

      //checking time axes are the same
      fgettr(fp_x,&tr_x);
      fgettr(fp_y,&tr_y);

      if( tr_x.ns!=tr_y.ns || 
	  tr_x.dt!=tr_y.dt ||
	  tr_x.delrt!=tr_y.delrt ){
	RVLException e;
	e << "Time axes do not match!\n"
	  << "tr_x.ns="<<tr_x.ns<<", tr_x.dt="<<tr_x.dt<<", tr_x.delrt="<<tr_x.delrt<<"\n"
	  << "tr_y.ns="<<tr_y.ns<<", tr_y.dt="<<tr_y.dt<<", tr_y.delrt="<<tr_y.delrt<<"\n";
	throw e;
      }

      fseeko(fp_x,0L,SEEK_SET);
      fseeko(fp_y,0L,SEEK_SET);

      //allocating buffer
      buff = new float[tr_x.ns*MPS_dim];
      if(buff==NULL){
	RVLException e;
	e << "Trouble allocating memory for buffer array.\n";
	throw e;
      }
      for( int i_b=0; i_b<MPS_dim; i_b++ ){
	for( int i_t=0; i_t<tr_x.ns; i_t++ ){
	  buff[i_t+i_b*tr_x.ns] = 0.f;
	}
      }
      
      //applying stacking MPS operation
      for( int i_src=0; i_src<N_src; i_src++ ){	

	for( int i_b=0; i_b<MPS_dim; i_b++ ){
	  fgettr(fp_y,&tr_y);

	  for( int i_t=0; i_t<tr_y.ns; i_t++ ){
	    buff[i_t+i_b*tr_y.ns] += tr_y.data[i_t];
	  }
	}
      }
      
      //writing out
      for( int i_b=0; i_b<MPS_dim; i_b++ ){
	fgetpos(fp_x,&pos_x);
	fgettr(fp_x,&tr_x);
	for( int i_t=0; i_t<tr_x.ns; i_t++ ){
	  tr_x.data[i_t] = buff[ i_t+i_b*tr_x.ns ];
	}
	fsetpos(fp_x,&pos_x);
	fputtr(fp_x,&tr_x);	  
      }
      
      fflush(fp_x);
      iwave_fclose(fp_y);
      iwave_fclose(fp_x);

      delete[] buff;

#ifdef VERBOSE_MJB
      cerr << "Leaving MPS_stack_kern\n";
#endif

    }
    catch(RVLException &e){
      e << "ERROR from MPS_stack_kern!\n";
      throw e;
    }
  }


  //////////////////////////////////////////////////////////////////////////
  // MPS_spread ////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  
  //----------------------------------------------------------------------//
  void MPS_spread::param_set( Vector<ireal> const &x,
			      Vector<ireal> const &y ) const{
  //----------------------------------------------------------------------//  
    try{
      
      AssignParams ap_x(*pars,"spread_x",stderr);
      x.eval(ap_x);

      AssignParams ap_y(*pars,"spread_y",stderr);
      y.eval(ap_y);

    }
    catch(RVLException &e){
      e << "ERROR from MPS_spread::param_set\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  ostream & MPS_spread::write(ostream & str) const {
  //----------------------------------------------------------------------//
    str<<"MPS_spread object\n";
    str<<"domain:\n";
    dom.write(str);
    str<<"range:\n";
    rng.write(str);
    return str; 
  }  

  //----------------------------------------------------------------------//
  void MPS_spread::apply( const Vector<ireal> &x,
			  Vector<ireal> &y ) const{
  //----------------------------------------------------------------------//
    try{

      SpaceTest(this->getMPSDomain(),x,"MPS_spread::apply (dom)");
      SpaceTest(this->getMPSRange(),y,"MPS_spread::apply (rng)");

      param_set(x,y);

      //zero output 
      y.zero();

      //extracting file names
      string x_file = valparse<string>(*pars,"spread_x");
      string y_file = valparse<string>(*pars,"spread_y");

      //stack kernel
      MPS_spread_kern( rng.get_s_pos().size(),
		       dom.get_dim(),
		       x_file,
		       y_file );

    }
    catch (RVLException & e) {
      e<<"ERROR from MPS_spread::apply\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void MPS_spread::applyAdj( const Vector<ireal> &y,
			     Vector<ireal> &x ) const {
  //----------------------------------------------------------------------//
    try{

      SpaceTest(this->getMPSDomain(),x,"MPS_spread::applyAdj (dom)");
      SpaceTest(this->getMPSRange(),y,"MPS_spread::applyAdj (rng)");

      param_set(x,y);

      //zero output 
      x.zero();

      //extracting file names
      string x_file = valparse<string>(*pars,"spread_x");
      string y_file = valparse<string>(*pars,"spread_y");

      //stack kernel
      MPS_stack_kern( rng.get_s_pos().size(),
		      dom.get_dim(),
		      x_file,
		      y_file );
    }
    catch (RVLException & e) {
      e<<"ERROR from MPS_spread::applyAdj\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void MPS_spread::sanity_checks() const {
  //----------------------------------------------------------------------//
    try{

      if(dom.get_s_pos().size()!=1){
	RVLException e;
	e << "number of sources in domain space is not equal to 1!\n";
	throw e;
      }
      
      vector<MPS_base> dom_basis = dom.get_MPS_basis();
      vector<MPS_base> rng_basis = rng.get_MPS_basis();
      
      int rng_dim = rng.get_dim();
      int dom_dim = dom.get_dim();

      if(rng_dim!=dom_dim){
	RVLException e;
	e << "range and domain MPS spaces are not of the same dimension!\n"
	  << "rng_dim="<<rng_dim <<", dom_dim="<<dom_dim<<"\n";      
	throw e;
      }

      for( int i_b=0; i_b<dom_dim; i_b++ ){
	if(rng_basis[i_b].NT!=dom_basis[i_b].NT){
	  RVLException e;
	  e << "range and domain MPS base i_b"<<i_b<<" do not have the same NT!\n"
	    << "rng NT="<<rng_basis[i_b].NT<<", dom NT="<<dom_basis[i_b].NT<<"\n";
	  throw e;
	}
	
	for( int i_t=0; i_t<dom_basis[i_b].NT; i_t++ ){
	  for( int d=0; d<3; d++ ){
	    if(rng_basis[i_b].derv[i_t].coor[d]!=dom_basis[i_b].derv[i_t].coor[d]){
	      RVLException e;
	      e << "range and domain MPS base i_b"<<i_b<<", term i_t="<<i_t<<" do not have the same derv!\n"
		<< "rng derv["<<d<<"]="<<rng_basis[i_b].derv[i_t].coor[d] <<", "
		<< "dom derv["<<d<<"]="<<dom_basis[i_b].derv[i_t].coor[d] <<"\n";
	      throw e;	      
	    }
	  }
	}

	if( dom_basis[i_b].RHS_flags.size()!= rng_basis[i_b].RHS_flags.size() ){
	  RVLException e;
	  e << "range and domain MPS base i_b"<<i_b<<" do not have the same number of RHS flags!\n"
	    << "rng NT="<<rng_basis[i_b].RHS_flags.size()<<", dom NT="<<dom_basis[i_b].RHS_flags.size()<<"\n";
	  throw e;	  
	}

	for( int i_r=0; i_r<dom_basis[i_b].RHS_flags.size(); i_r++ ){
	  if(dom_basis[i_b].RHS_flags[i_r]!=rng_basis[i_b].RHS_flags[i_r]){
	    RVLException e;
	    e << "range and domain MPS base i_b"<<i_b<<" RHS flag i_r="<<i_r<<" are not equal!\n"
	      << "rng RHS_flag="<<rng_basis[i_b].RHS_flags[i_r]<<", dom RHS_flag="<<dom_basis[i_b].RHS_flags[i_r]<<"\n";
	    throw e;	  	
	  }  
	}
      }
      
    }
    catch (RVLException & e) {
      e<<"ERROR from MPS_spread::sanity_checks\n";
      throw e;
    }
  }
}
