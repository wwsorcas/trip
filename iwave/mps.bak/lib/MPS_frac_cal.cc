// MPS_frac_cal.cc
// Author: Mario J. Bencomo
// last modified: 04/10/18

#include "MPS_frac_cal.hh"

//#define VERBOSE_MPS_FRAC_CAL_CC
#define MY_MAX(x,y) ((x)>(y))?(x):(y)
#define MY_MIN(x,y) -(MY_MAX(-(x),-(y)))
#define MINUS_TO_POW(p) ((p)%2)?(-1):(1)

namespace TSOpt{


  //-----------------------------------------------------------------------------------
  int trunc_fact(int m,int s){
  //-----------------------------------------------------------------------------------
    if( m<0 || s<0 ){
      RVL::RVLException e;
      e << "ERROR from trunc_fact: s="<<s<<" or m="<<m<<" is negative!\n";
      throw e;
    }
    if( s>m ){
      RVL::RVLException e;
      e << "ERROR from trunc_fact: s="<<s<<" is greater than m="<<m<<"!\n";
      throw e;
    }
    int prod = 1;
    for(int i=m; i>(m-s); i--) prod*=i;
    return prod;
  }


  //----------------------------------------------------------------------// 
  int nchoosek(int n, int k){
  //----------------------------------------------------------------------// 
    try{
      if(n<0||k<0){
	RVLException e;
	e << "n and/or k is/are negative!\n"
	  << "k="<< k <<", n="<< n <<"\n";
	throw e;
      }
      if( (k>n) && (n!=0)){
	RVLException e;
	e << "k is greater than n!\n"
	  << "k=" <<k<<", n="<<n<<"\n";
	throw e;
      }
      if( n>30 ){
	RVLException e;
	e << "n is greater that 30!\n";
	throw e;
      }

      if( (k==0) || (k==n)) return 1;
      return nchoosek(n-1, k-1)+nchoosek(n-1,k);
    }
    catch(RVLException &e){
      e << " ERROR from nchoosek()\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------// 
  void frac_deriv( float q, 
		   float scal,
		   int N, 
		   const float *x, 
		   float *y, 
		   float dt, 
		   bool adj ){
  //----------------------------------------------------------------------//     
    try{

      if( q<0 || N<=0 ){
	RVLException e;
	e << "Derivative order (q="<<q<<") or size of arrays (N="<<N<<") "
	  << "must be greater than 0!\n";
	throw e;	
      }

      bool is_int = (floor(q)-q)==0;
      double sign;
      double arg;
      int int_arg;
      float tmp;
      
      ////////////
      // 0 case //
      ////////////
      if( q==0 ){
	for(int i=0; i<N; i++) y[i]=x[i];
      }

      else{
	
	//////////////////
	// Forward case //
	//////////////////
	if(!adj){
	  
	  /* Integer case */
	  if(is_int){
	  
	    for(int i=0; i<N; i++){
	      y[i]=0.0;
	      
	      for(int k=(MY_MAX(i-q,0)); k<=i; k++){
		y[i] += (MINUS_TO_POW((i-k))) * nchoosek(q,(i-k)) * x[k];
	      }
	      y[i] *= pow(dt,-q)*scal;
	    }
	  }

	  /* Fractional case */
	  else{
	    
	    for(int i=0; i<N; i++){	  
	      y[i]=0.0;
	      
	      for(int k=0; k<=i; k++){
		sign = 1.0;
		arg = q-i+k+1;
		int_arg = (int)(floor(arg));
		
		if( arg<0 && int_arg%2 ) sign = -1.0;
		
		tmp = lgamma((double)(q+1)) - lgamma((double)(arg)) - lgamma((double)(i-k+1));
		tmp = exp(tmp) * sign;
		y[i] += (MINUS_TO_POW((i-k))) * tmp * x[k];
	      }
	      y[i] *= pow(dt,-q)*scal;
	    }
	  }
	}

	//////////////////
	// Adjoint case //
	//////////////////
	else{
	  
	  /* Integer case */
	  if(is_int){

	    for(int i=0; i<N; i++){
	      y[i]=0.0;
	      
	      for(int k=i; k<=(MY_MIN(N-1,i+q)); k++){
		y[i] += (MINUS_TO_POW(k-i)) * nchoosek(q,(k-i)) * x[k];
	      }
	      y[i] *= pow(dt,-q)*scal;
	    }
	  }
	  
	  /* Fractional case */
	  else{
	    for(int i=0; i<N; i++){	      
	      y[i]=0.0;

	      for( int k=i; k<=N-1; k++){
		double sign = 1.0;
		double arg = q-k+i+1;
		int int_arg = (int)(floor(arg));

		if( arg<0 && int_arg%2 ) sign = -1.0;
		
		tmp = lgamma((double)(q+1)) - lgamma((double)(arg)) - lgamma((double)(k-i+1));
		tmp = exp(tmp) * sign;
		y[i] += (MINUS_TO_POW((i-k))) * tmp * x[k];
	      }
	      y[i] *= pow(dt,-q)*scal;
	    } 
	  }
	}
      }
    }
    catch(RVLException &e){
      e << " ERROR from frac_deriv!\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------// 
  void frac_integ( float q, 
		   float scal,
		   int N, 
		   const float *x, 
		   float *y,
		   float dt,
		   bool adj ){
  //----------------------------------------------------------------------// 
    try{

      if( q<=0 || N<=0 ){
	RVLException e;
	e << "Integral order (q="<<q<<") or size of arrays (N="<<N<<") "
	  << "must be greater than 0!\n";
	throw e;	
      }

      bool is_int = (floor(q)-q)==0;
      float tmp;
      
      //////////////////
      // Forward case //
      //////////////////
      if(!adj){
	
	/* Integer case */
	if(is_int){
	  
	  for(int i=0; i<N; i++){
	    y[i]=0.0;
	    
	    for(int k=0; k<=i; k++){
	      y[i] += trunc_fact(q+i-k-1,q-1) * (1.0/factorial(q-1)) * x[k];
	    }
	    y[i] *= pow(dt,q)*scal;
	  }
	}
	
	/* Fractional case */
	else{
	  
	  for(int i=0; i<N; i++){
	    y[i] = 0.0;
	    
	    for(int k=0; k<=i; k++){
	      tmp = lgamma((double)(q+i-k)) - lgamma((double)(i-k+1)) - lgamma(q);
	      tmp = exp(tmp);
	      y[i] +=  tmp * x[k];
	    }
	    y[i] *= pow(dt,q)*scal;
	  }
	}
      }
      
      //////////////////
      // Adjoint case //
      //////////////////
      else{
	
	/* Integer case */
	if(is_int){

	  for(int i=0; i<N; i++){
	    y[i]=0.0;

	    for(int k=i; k<N; k++){
	      y[i] += trunc_fact(q+k-i-1,q-1) * (1.0/factorial(q-1)) * x[k];
	    }
	    y[i] *= pow(dt,q)*scal;
	  }
	}
	
	/* Fractional case */
	else{
	  
	  for(int i=0; i<N; i++){
	    y[i]=0.0;
		
	    for(int k=i; k<N; k++){
	      tmp = lgamma((double)(q+k-i)) - lgamma((double)(k-i+1)) - lgamma(q);
	      tmp = exp(tmp);
	      y[i] +=  tmp * x[k];
	    }
	    y[i] *= pow(dt,q)*scal;
	  }
	}	       
      }
    }
    catch(RVLException &e){
      e << " ERROR from frac_integ!\n";
      throw e;
    }
  }  

  //----------------------------------------------------------------------// 
  void MPS_frac_cal_kern( vector<float> const order,
			  vector<float> const scalars,
			  int N_src,
			  string in_file,
			  string out_file,
			  bool adj,
			  bool inv ){
  //----------------------------------------------------------------------// 
    try{
      
#ifdef VERBOSE_MPS_FRAC_CAL_CC
      cerr << "*************************************************\n"
           << "  Inside MPS_frac_cal_kern\n"
           << "*************************************************\n\n"
	   << "in  file = " << in_file << "\n"
	   << "out file = " << out_file << "\n"
	   << "adj flag = " << adj << "\n"
	   << "inv flag = " << inv << "\n"
	   << "order:";
      for(int i=0; i<order.size(); i++)	cerr << order[i] << ",";
      cerr <<"\n";
      cerr << "scalars:";
      for(int i=0; i<scalars.size(); i++) cerr << scalars[i] << ",";
      cerr << "\n";
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
      
      //setting x_arr and y_arr
      x_arr = new float[NTR*NS];
      fp_x = iwave_const_fopen(in_file.c_str(),"r",NULL,stderr);
      for( int i=0; i<NTR; i++){
	fgettr(fp_x,&tr);
	for( int i_t=0; i_t<NS; i_t++ ){
	  x_arr[i_t+i*NS] = tr.data[i_t]; 
	}
      }
      y_arr = new float[NTR*NS];

      
      int offset;
      ////////////////////////
      //looping over sources//
      ////////////////////////
      for(int i_src=0; i_src<N_src; i_src++){

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
			x_arr+offset,
			y_arr+offset,
			dt, adj );
	  }
	  else{
	    q *= -1;
	    frac_integ( q, scal, NS,
			x_arr+offset,
			y_arr+offset,
			dt, adj );
	  }

	}//MPS basis loop	
      }//sources loop
      
      //writing out  
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

      if(x_arr!=NULL) delete[] x_arr;
      if(y_arr!=NULL) delete[] y_arr;
    }
    catch(RVLException &e){
      e << "ERROR from MPS_frac_cal_kern!\n";
      throw e;
    }
  }


  //////////////////////////////////////////////////////////////////////////
  // MPS_frac_cal //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  
  //----------------------------------------------------------------------// 
  MPS_frac_cal::MPS_frac_cal( MPS_frac_cal const & op )
    : pars(NULL), 
      mps_sp(op.mps_sp),
      scalars(op.scalars),
      orders(op.orders),
      in_key(op.in_key),
      out_key(op.out_key){
  //----------------------------------------------------------------------// 
    try{
      //deep copy of pars
      pars = ps_new();
      if(ps_copy(&pars,*op.pars)){
	RVLException e;
	e << "failed to copy parameter table\n";
	throw e;
      }
    }
    catch(RVLException & e){
      e << "ERROR from MPS_frac_cal copy constructor\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------// 
  MPS_frac_cal::MPS_frac_cal( PARARRAY const &_pars, 
			      MPS_Space const &_sp,
			      float c,
			      float order_0,
			      float order_d,
			      string _in_key,
			      string _out_key)
    : pars(NULL), 
      mps_sp(_sp),
      in_key(_in_key),
      out_key(_out_key){
  //----------------------------------------------------------------------// 
    try{
 
      if(c==0){
	RVLException e;
	e << "c=0?!\n";
	throw e;
      }

      //copying pars
      pars = ps_new();
      if(ps_copy(&pars,_pars)){
	RVLException e;
	e << "failed to copy parameter table\n";
	throw e;
      }
      
      int derv_ord;
      int tmp;
      float q, scal;
      vector<MPS_base> MPS_basis = mps_sp.get_MPS_basis();

      //loop over MPS basis
      for(int i=0; i<mps_sp.get_dim(); i++){
	derv_ord=0;
	//loop over basis terms
	for(int j=0; j<MPS_basis[i].NT; j++){
	  tmp=0;
	  //loop over spatial coors
	  for(int d=0; d<mps_sp.get_s_dim(); d++){
	     tmp += MPS_basis[i].derv[j].coor[d]; 
	  }
	  derv_ord = MY_MAX(derv_ord,tmp);
	}

	q = order_0 + derv_ord*order_d;
	orders.push_back(q);

	//initializing scalar in [m/s], assuming c in [km/s]
	if(derv_ord==0)
	  scal = 1.0f;
	else{
	  scal = c*1000;
	  for(int d=1; d<derv_ord; d++) scal*=scal;
	  scal = 1.0f/scal;
	  //scal = (inv)?(scal):(1.0f/scal);
	}
	scalars.push_back(scal);
      }

      //running sanity check
      sanity_checks();
      
    }
    catch(RVLException & e){
      e << "ERROR from MPS_frac_cal constructor\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------// 
  MPS_frac_cal::~MPS_frac_cal(){
  //----------------------------------------------------------------------// 
    ps_delete(&pars);
  }

  //----------------------------------------------------------------------// 
  void MPS_frac_cal::sanity_checks(){
  //----------------------------------------------------------------------// 
    try{
      
    } 
    catch(bad_cast){
      RVLException e;
      e << "bad cast!\n";
      throw e;
    }
    catch(RVLException &e){
      e << "ERROR from MPS_frac_cal::sanity_checks()\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void MPS_frac_cal::param_set( Vector<ireal> const &x,
				Vector<ireal> const &y ) const{
  //----------------------------------------------------------------------//  
    try{
  
      AssignParams ap_x(*pars,in_key,stderr);
      x.eval(ap_x);

      AssignParams ap_y(*pars,out_key,stderr);
      y.eval(ap_y);
  
    }
    catch(RVLException &e){
      e << "ERROR from MPS_frac_cal::param_set\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  ostream & MPS_frac_cal::write(ostream & str) const {
  //----------------------------------------------------------------------//
    str<<"MPS_frac_cal object\n";
    str<<"order:";
    for(int i=0; i<orders.size(); i++)
      str << orders[i] << ",";
    str<<"\nscalars:\n";
    for( int ib=0; ib<scalars.size(); ib++)
      str << scalars[ib] <<",";
    str<<"\ndomain and range:\n";
    mps_sp.write(str);

    return str; 
  }  

  //----------------------------------------------------------------------//
  void MPS_frac_cal::apply( const Vector<ireal> &x,
			    Vector<ireal> &y ) const{
  //----------------------------------------------------------------------//
    try{
     
      SpaceTest(this->get_MPS_Domain(),x,"MPS_frac_cal::apply (dom)");
      SpaceTest(this->get_MPS_Range(),y,"MPS_frac_cal::apply (rng)");

      //setting filenames in param
      param_set(x,y);

      // zero output 
      y.zero();

      //extracting file names
      string x_file = valparse<string>(*pars,in_key);
      string y_file = valparse<string>(*pars,out_key);
      
      bool adj = false;
      bool inv = false;

      //diffint kernel
#ifdef IWAVE_USE_MPI
      MPS_frac_cal_kern_MPI( orders,
			     scalars,
			     mps_sp.get_s_pos().size(),
			     x_file,
			     y_file,
			     adj,
			     inv );
#else
      MPS_frac_cal_kern( orders,
			 scalars,
			 mps_sp.get_s_pos().size(),
			 x_file,
			 y_file,
			 adj,
			 inv );
#endif

    }
    catch (RVLException & e) {
      e<<"ERROR from MPS_frac_cal::apply\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void MPS_frac_cal::applyAdj( const Vector<ireal> &y,
			   Vector<ireal> &x ) const {
  //----------------------------------------------------------------------//
    try{
     
      SpaceTest(this->get_MPS_Domain(),x,"MPS_conv::applyAdj (rng)");
      SpaceTest(this->get_MPS_Range(),y,"MPS_conv::applyAdj (dom)");

      //setting filenames in param
      param_set(y,x);

      // zero output 
      x.zero();

      //extracting file names
      string y_file = valparse<string>(*pars,in_key);
      string x_file = valparse<string>(*pars,out_key);

      bool adj = true;
      bool inv = false;

      //diffint adj kernel
#ifdef IWAVE_USE_MPI
      MPS_frac_cal_kern_MPI( orders,
			     scalars,
			     mps_sp.get_s_pos().size(),
			     y_file,
			     x_file,
			     adj,
			     inv );
#else
      MPS_frac_cal_kern( orders,
			 scalars,
			 mps_sp.get_s_pos().size(),
			 y_file,
			 x_file,
			 adj,
			 inv );
#endif
    }
    catch (RVLException & e) {
      e<<"ERROR from MPS_frac_cal::applyAdj\n";
      throw e;
    }
  }


  //----------------------------------------------------------------------//
  void MPS_frac_cal::applyInv( const Vector<ireal> &y,
			       Vector<ireal> &x ) const{
  //----------------------------------------------------------------------//
    try{
     
      SpaceTest(this->get_MPS_Domain(),x,"MPS_frac_cal::apply (dom)");
      SpaceTest(this->get_MPS_Range(),y,"MPS_frac_cal::apply (rng)");

      //setting filenames in param
      param_set(y,x);

      // zero output 
      x.zero();

      //extracting file names
      string y_file = valparse<string>(*pars,in_key);
      string x_file = valparse<string>(*pars,out_key);

      bool adj = false;
      bool inv = true;

      //diffint kernel
#ifdef IWAVE_USE_MPI
      MPS_frac_cal_kern_MPI( orders,
			     scalars,
			     mps_sp.get_s_pos().size(),
			     y_file,
			     x_file,
			     adj,
			     inv );
#else
      MPS_frac_cal_kern( orders,
			 scalars,
			 mps_sp.get_s_pos().size(),
			 y_file,
			 x_file,
			 adj,
			 inv );
#endif

    }
    catch (RVLException & e) {
      e<<"ERROR from MPS_frac_cal::applyInv\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void MPS_frac_cal::applyInvAdj( const Vector<ireal> &x,
				  Vector<ireal> &y ) const {
  //----------------------------------------------------------------------//
    try{
     
      SpaceTest(this->get_MPS_Domain(),x,"MPS_conv::applyAdj (rng)");
      SpaceTest(this->get_MPS_Range(),y,"MPS_conv::applyAdj (dom)");

      //setting filenames in param
      param_set(x,y);

      // zero output 
      y.zero();

      //extracting file names
      string x_file = valparse<string>(*pars,in_key);
      string y_file = valparse<string>(*pars,out_key);

      bool adj = true;
      bool inv = true;

      //diffint adj kernel
#ifdef IWAVE_USE_MPI
      MPS_frac_cal_kern_MPI( orders,
			     scalars,
			     mps_sp.get_s_pos().size(),
			     x_file,
			     y_file,
			     adj,
			     inv );
#else
      MPS_frac_cal_kern( orders,
			 scalars,
			 mps_sp.get_s_pos().size(),
			 x_file,
			 y_file,
			 adj,
			 inv );
#endif
    }
    catch (RVLException & e) {
      e<<"ERROR from MPS_frac_cal::applyInvAdj\n";
      throw e;
    }
  }


}//end TSOpt
