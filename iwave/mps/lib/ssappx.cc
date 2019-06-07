// ssappx.cc
// Author: Mario J. Bencomo
// last modified: 01/26/16

#include "ssappx.hh"

//#define VERBOSE_MJB
//#define MINUS_TO_POW(p) ((p)%2)?(-1):(1)

namespace TSOpt {

  //----------------------------------------------------------------------// 
  int factorial(int n){
  //----------------------------------------------------------------------// 
    try{
      if(n<0){
	RVLException e;
	e << "n is negative!\n";
	throw e;
      }
      if(n>30){
	RVLException e;
	e << "n is greater than 30!\n";
	throw e;
      }

      int result = 1;
      for (int i=1; i<=n; i++)
	result = result * i;
      return result;
    }
    catch(RVLException &e){
      e << " ERROR from factorial()\n";
      throw e;
    }
  }

  //-----------------------------------------------------------------------------------
  void midx_update( vector< vector<int> > &MIDX, 
		    vector<int> const &midx, 
		    int N, int i, int ell ){
  //-----------------------------------------------------------------------------------
    try{
      int sz = midx.size();
      
      if( i>=0 ){

	int diff = -midx[i];
	if( i==(sz-1) )
	  diff += N;	  
	else
	  diff += midx[i+1];

	for( int k=1; k<diff; k++ ){
	  if( (midx[i]+k)!=ell ){
	    vector<int> new_midx = midx;
	    new_midx[i] += k;
	    MIDX.push_back(new_midx);
	    midx_update( MIDX, new_midx, N, i-1, ell);
	  }
	}
      }
    }
    catch(RVL::RVLException &e){
      e << "ERROR from midx_update!\n";
      throw e;
    }
  }


  //-----------------------------------------------------------------------------------
  vector< vector<int> > midx_set( int q, int s, int ell ){
  //-----------------------------------------------------------------------------------
    try{

      int N = q+s;
      int sz = q-1;
      vector< vector<int> > MIDX;
      
      //starting midx
      vector<int> midx(sz);
      for(int i=0; i<sz; i++){
	if( i<ell )
	  midx[i] = i;
	else
	  midx[i] = i+1;
      }
      MIDX.push_back(midx);
      midx_update(MIDX,midx,N,sz-1,ell);

      return MIDX;
    }
    catch(RVL::RVLException &e){
      e << "ERROR from midx_set!\n";
      throw e;
    }
  }


  //-----------------------------------------------------------------------------------
  float poly_ell( float x, int ell, float h, int s, int q){
  //-----------------------------------------------------------------------------------
    try{

      int N = q+s;
      
      //computing numerator
      float tmp = 1.f;

      if( q>1 ){
	
	vector< vector<int> > MIDX = midx_set(q,s,ell);
	float tmp_sum = 0.f;

	//looping over multi-indexes
	for( int i=0; i<MIDX.size(); i++){

	  float tmp_prd = 1.f;

	  //looping over multi-index elements
	  for( int j=0; j<MIDX[i].size(); j++){
	    tmp_prd *= h*(MIDX[i][j]-ell) + x;
	  }
	  tmp_sum += tmp_prd;
	}
	tmp = pow(-1,q-1) * tmp_sum;
      }
      else if(q!=1){
	RVL::RVLException e;
	e << "q="<<q<<" must be >=1\n";
	throw e;
      }

      float numer = factorial(s) * pow(-1,s) * tmp;

      //computing denominator
      float denom = 1.f;

      for( int m=0; m<N; m++ ){
	if( m!= ell ){
	  denom *= ell-m;
	}
      }
      denom *= pow(h,N-1);
      
      return numer/denom;
    }
    catch(RVL::RVLException &e){
      e << "ERROR from poly_ell!\n";
      throw e;
    }
  }

  //-----------------------------------------------------------------------------------
  float ssappx( float x, float h, int s, int q){
  //-----------------------------------------------------------------------------------
    try{

#ifdef VERBOSE_MJB
      cerr << "        Inside ssappx, with \n"
	   << "          x="<<x<<"\n"
	   << "          h="<<h<<"\n" 
	   << "          s="<<s<<"\n"
	   << "          q="<<q<<"\n";
#endif
      if(h<=0){
	RVL::RVLException e;
	e << "h="<<h<<" must be >0\n";
	throw e;
      }

      if(s<0){
	RVL::RVLException e;
	e << "s="<<h<<" must be >=0\n";
	throw e;	
      }

      if(q<=0){
	RVL::RVLException e;
	e << "q="<<q<<" must be >0\n";
	throw e;
      }

      int N = q+s;
      float a1 = -N*h*0.5;
      int ell=-1;
      float y=0.f;
      
      //cheking what interval x is in
      for( int j=0; j<N; j++){
	float aj = a1 + j*h;
	float aj1 = aj + h;
	if( (aj<=x) && (x<aj1) ){
	  ell = j;
	  break;
	} 
      } 

#ifdef VERBOSE_MJB
      cerr << "        x is in ell="<<ell<<" interval.\n";
#endif
      //calling ell-th polynomial
      if( ell>=0 ){
	y = poly_ell(x,ell,h,s,q);
      }

#ifdef VERBOSE_MJB
      cerr <<"        ouput="<<y<<"\n";
#endif
      
      return y;
    }
    catch(RVL::RVLException &e){
      e << "ERROR from ssappx!\n";
      throw e;
    }
  }
}
