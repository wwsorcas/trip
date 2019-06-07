// MPS_Space_Examples.cc
// Author: Mario J. Bencomo
// last modified: 12/12/16

#include "MPS_Space_Examples.hh"

namespace TSOpt {
  

  //////////////////////////////////////////////////////////////////////////
  // CanScal_MPS_Space /////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////  

  //----------------------------------------------------------------------//
  CanScal_MPS_Space::CanScal_MPS_Space(	MPS_KEYS _mks,
					PARARRAY const& pars,
					bool make )
    : MPS_Space(_mks,pars,make) {
  //----------------------------------------------------------------------//
    try{
      
      //setting order
      order = valparse<int>(pars,"MPS_ord");

      //setting type
      std::stringstream ss;
      ss << "Canonical scalar MPS_Space in " << s_dim << "D of order "<< order;
      type = ss.str();

      //setting MPS_basis
      build_MPS_basis();
      	
      //case where SEGY core is made from scratch
      if( make ){
	gen_SEGY_core(pars);
      }
      sanity_checks();
    }
    catch(RVLException &e){
      e << "ERROR from CanScal_MPS_Space constructor\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void CanScal_MPS_Space::build_MPS_basis(){
  //----------------------------------------------------------------------//
    try{
      
      //initializing monopole basis
      if(order>=0){
	
	//monopole term delta(x)
	MPS_base b; 
	b.NT = 1;

	b.RHS_flags.resize(1);
	b.RHS_flags[0]=1;

	b.derv.resize(b.NT);
	for(int i=0; i<s_dim; i++)
	  b.derv[0].coor[i]=0;

	MPS_basis.push_back(b);
      }

      //initializing dipole basis
      if(order>=1){
		
	for(int i=0; i<s_dim; i++){
	  //dipole term ∂x_i delta(x)
	  MPS_base b; 
	  b.NT = 1;

	  b.RHS_flags.resize(1);
	  b.RHS_flags[0]=1;

	  b.derv.resize(b.NT);
	  for(int j=0; j<s_dim; j++)
	    b.derv[0].coor[j]=0;
	  b.derv[0].coor[i]++;

	  MPS_basis.push_back(b);
	}
      }

      //initializing quadrupole basis
      if(order>=2){

	MPS_base b;
	b.NT = 1;
	
	b.RHS_flags.resize(1);
	b.RHS_flags[0]=1;
	b.derv.resize(b.NT);
	
	if(s_dim==2){
	  //quadrupole term ∂x_1 ∂x_1 delta(x)
	  b.derv[0].coor[0]=2;
	  b.derv[0].coor[1]=0;
	  MPS_basis.push_back(b);

	  //quadrupole term ∂x_2 ∂x_2 delta(x)
	  b.derv[0].coor[0]=0;
	  b.derv[0].coor[1]=2;
	  MPS_basis.push_back(b);

	  //quadrupole term ∂x_1 ∂x_2 delta(x)
	  b.derv[0].coor[0]=1;
	  b.derv[0].coor[1]=1;
	  MPS_basis.push_back(b);
	}

	if(s_dim==3){

	  //quadrupole term ∂x_1 ∂x_1 delta(x)
	  b.derv[0].coor[0]=2;
	  b.derv[0].coor[1]=0;
	  b.derv[0].coor[2]=0;
	  MPS_basis.push_back(b);

	  //quadrupole term ∂x_2 ∂x_2 delta(x)
	  b.derv[0].coor[0]=0;
	  b.derv[0].coor[1]=2;	  
	  b.derv[0].coor[2]=0;
	  MPS_basis.push_back(b);

	  //quadrupole term ∂x_3 ∂x_3 delta(x)
	  b.derv[0].coor[0]=0;
	  b.derv[0].coor[1]=0;	  
	  b.derv[0].coor[2]=2;
	  MPS_basis.push_back(b);

	  //quadrupole term ∂x_1 ∂x_2 delta(x)
	  b.derv[0].coor[0]=1;
	  b.derv[0].coor[1]=1;
	  b.derv[0].coor[2]=0;
	  MPS_basis.push_back(b);

	  //quadrupole term ∂x_1 ∂x_3 delta(x)
	  b.derv[0].coor[0]=1;
	  b.derv[0].coor[1]=0;
	  b.derv[0].coor[2]=1;
	  MPS_basis.push_back(b);

	  //quadrupole term ∂x_2 ∂x_3 delta(x)
	  b.derv[0].coor[0]=0;
	  b.derv[0].coor[1]=1;	  
	  b.derv[0].coor[2]=1;
	  MPS_basis.push_back(b);
	}
	/*
	for(int i=0; i<s_dim; i++){
	  for(int j=0; j<s_dim; j++){
	  
	    //dipole term ∂x_i ∂x_j delta(x)
	    MPS_base b; 
	    b.NT = 1;

	    b.RHS_flags.resize(1);
	    b.RHS_flags[0]=1;
	    
	    b.derv.resize(b.NT);
	    for(int k=0; k<s_dim; k++)
	      b.derv[0].coor[k]=0;
	    
	    b.derv[0].coor[i]++;
	    b.derv[0].coor[j]++;
	    
	    MPS_basis.push_back(b);
	  }
	}
	*/
      }

      if(order>2){
	RVLException e;
	e << "MPS order > 2, currently not available!\n";
	throw e;
      }

    }
    catch(RVLException &e){
      e << "ERROR from CanScal_MPS_Space::build_MPS_basis()!\n";
      throw e;
    }
  }


  //////////////////////////////////////////////////////////////////////////
  // CanVec_MPS_Space /////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////  

  //----------------------------------------------------------------------//
  CanVec_MPS_Space::CanVec_MPS_Space( MPS_KEYS _mks,
				      PARARRAY const& pars,
				      bool make )
    : MPS_Space(_mks,pars,make) {
  //----------------------------------------------------------------------//
    try{
      
      //setting order
      order = valparse<int>(pars,"MPS_ord");
      
      //setting type
      std::stringstream ss;
      ss << "Canonical vector MPS_Space in " << s_dim << "D of order " << order;
      type = ss.str();

      //setting MPS_basis
      build_MPS_basis();
      	
      //case where SEGY core is made from scratch
      if( make ){
	gen_SEGY_core(pars);
      }
      sanity_checks();
    }
    catch(RVLException &e){
      e << "ERROR from CanVec_MPS_Space constructor type\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void CanVec_MPS_Space::build_MPS_basis(){
  //----------------------------------------------------------------------//
    try{

      //initializing monopole basis
      if(order>=0){
	
	//looping over canonical direction e_i
	for( int i=0; i<s_dim; i++){

	  //monopole term e_i delta(x)
	  MPS_base b; 
	  b.NT = 1;

	  b.RHS_flags.resize(s_dim);
	  b.derv.resize(b.NT);
	  
	  for(int j=0; j<s_dim; j++){
 	    b.RHS_flags[j]=0;
	    b.derv[0].coor[j]=0;
	  }
	  b.RHS_flags[i]=1;
	  
	  MPS_basis.push_back(b);
	}
      }

      //initializing dipole basis
      if(order>=1){

	//looping over canonical directions e_i
	for(int i=0; i<s_dim; i++){
	
	  //looping over partial derivatives ∂x_j
	  for(int j=0; j<s_dim; j++){

	    //dipole term e_i ∂x_j delta(x)
	    MPS_base b; 
	    b.NT = 1;
	   
	    b.RHS_flags.resize(s_dim);
	    b.derv.resize(b.NT);

	    for(int k=0; k<s_dim; k++){
	      b.derv[0].coor[k]=0;
	      b.RHS_flags[k]=0;
	    }
	    b.RHS_flags[i]=1;
	    b.derv[0].coor[j]=1;

	    MPS_basis.push_back(b);
	  }
	}
      }

      if(order>1){
	RVLException e;
	e << "MPS oder > 1, currently not available!\n";
	throw e;
      }

    }
    catch(RVLException &e){
      e << "ERROR from CanVec_MPS_Space::build_MPS_basis()!\n";
      throw e;
    }
  }


  //////////////////////////////////////////////////////////////////////////
  // ExVec_MPS_Space ///////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////  

  //----------------------------------------------------------------------//
  ExVec_MPS_Space::ExVec_MPS_Space( MPS_KEYS _mks, 
				    PARARRAY const& pars,
				    bool make )
    : MPS_Space(_mks,pars,make) {
  //----------------------------------------------------------------------//
    try{


      if(s_dim!=2){
	RVLException e;
	e << "input spatial dimension is not 2!\n";
	throw e;
      }

      //setting type
      type = "Example vector MPS_Space in 2D";

      //setting MPS_basis
      build_MPS_basis();
      	
      //case where SEGY core is made from scratch
      if( make ){
	gen_SEGY_core(pars);
      }
      sanity_checks();
    }
    catch(RVLException &e){
      e << "ERROR from ExVec_MPS_Space constructor\n";
      throw e;
    }
  }
      
  //----------------------------------------------------------------------//
  void ExVec_MPS_Space::build_MPS_basis(){
  //----------------------------------------------------------------------//
    try{

      //b1 term = [1,0] delta(x_0,x_1)
      MPS_base b1; 
      b1.NT = 1;

      b1.RHS_flags.resize(2);
      b1.RHS_flags[0]=1;
      b1.RHS_flags[1]=0;

      b1.derv.resize(b1.NT);
      b1.derv[0].coor[0]=0;
      b1.derv[0].coor[1]=0;

      MPS_basis.push_back(b1);
      
      //b2 term = [0,1] delta(x_0,x_1)
      MPS_base b2;
      b2.NT = 1;

      b2.RHS_flags.resize(2);
      b2.RHS_flags[0]=0;
      b2.RHS_flags[1]=1;

      b2.derv.resize(b2.NT);
      b2.derv[0].coor[0]=0;
      b2.derv[0].coor[1]=0;
      
      MPS_basis.push_back(b2);

      //b3 term = [1,1] { d/dx_0 + d/dx_1 } delta
      MPS_base b3;
      b3.NT = 2;

      b3.RHS_flags.resize(2);
      b3.RHS_flags[0]=1;
      b3.RHS_flags[1]=1;

      b3.derv.resize(b3.NT);
      b3.derv[0].coor[0]=1;
      b3.derv[0].coor[1]=0;
      b3.derv[1].coor[0]=0;
      b3.derv[1].coor[1]=1;

      MPS_basis.push_back(b3);

    }
    catch(RVLException &e){
      e << "ERROR from ExVec_MPS_Space::build_MPS_basis()!\n";
      throw e;
    }
  }


  //////////////////////////////////////////////////////////////////////////
  // Scal_MPS_Space ////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////  

  //----------------------------------------------------------------------//
  Scal_MPS_Space::Scal_MPS_Space( MPS_KEYS _mks,
				  PARARRAY const& pars,
				  bool make )
    : MPS_Space(_mks,pars,make) {
  //----------------------------------------------------------------------//
    try{
      
      //setting order
      order.coor[0] = valparse<int>(pars,"MPS_ord_0");
      if(s_dim>1) order.coor[1] = valparse<int>(pars,"MPS_ord_1");
      if(s_dim>2) order.coor[2] = valparse<int>(pars,"MPS_ord_2");

      //settting type
      std::stringstream ss;
      ss << "Single scalar MPS_Space in " << s_dim << "D of order ";
      for(int i=0; i<s_dim; i++) ss << order.coor[i] << ",";
      type = ss.str();

      //setting MPS_basis
      build_MPS_basis();
      	
      //case where SEGY core is made from scratch
      if( make ){
	gen_SEGY_core(pars);
      }
      sanity_checks();
    }
    catch(RVLException &e){
      e << "ERROR from Scal_MPS_Space constructor\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void Scal_MPS_Space::build_MPS_basis(){
  //----------------------------------------------------------------------//
    try{

      MPS_base b;
      b.NT = 1;
 
      b.RHS_flags.resize(1);
      b.RHS_flags[0]=1;
      
      b.derv.resize(b.NT);
      for( int i=0; i<s_dim; i++ ){
	b.derv[0].coor[i] = order.coor[i];
      }

      MPS_basis.push_back(b);

    }
    catch(RVLException &e){
      e << "ERROR from Scal_MPS_Space::build_MPS_basis()!\n";
      throw e;
    }
  }

  
}//end TSOpt
