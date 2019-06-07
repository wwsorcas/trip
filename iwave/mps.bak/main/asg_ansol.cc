// asg_ansol.cc
// Author: Mario J. Bencomo
// last modification: 11/10/16

#include "asg_defn.hh"
#include "asg.hh"
#include "MPS_includes.hh"
#include "MPS_Space.hh"

using RVL::valparse;
using RVL::RVLException;
using TSOpt::Rtuple;
using TSOpt::Ituple;
using TSOpt::ex_pos;

const char *sdoc[] = {
  " ============================================================================",
  " asg_ansol.x ",
  " ============================================================================",
  " Authors: Mario J. Bencomo",
  " ",
  " Analytical solution for acoustics with variable density in first order form,",
  " for homogenous unbounded medium.",
  " ",
  " PDE system: ",
  "  dp/dt + bulk{ div v } = f(t,x)",
  "  dv/dt + buoy grad p = 0",
  " where p=pressure, v=velocity, and ",
  "   2-D, f(x,t) = f_0(t) (d/dx_0)^{s_0} (d/dx_1)^{s_1} delta(x-x^*)",
  "   3-D, f(x,t) = f_0(t) (d/dx_0)^{s_0} (d/dx_1)^{s_1} (d/dx_2)^{s_2} delta(x-x^*)",
  " Output is p(x_r,t) for given reciever locations specified by input SU file. ",
  " ",
  " Typical parameter list. May be copied, edited, and used for input: either",
  " include parameters on command line (for example in Flow), or place",
  " in file <foo> and include \"par=<foo>\" on command line. Any parameter",
  " included explicitly in command line overrides parameter with same key",
  " in par file.",
  " ",
  "  Invoke single threaded execution by ",
  " \"asg_ansol.x [parameters] [standalone install]\"",
  " ",
  " --------------------------- begin parameters ---------------------------",
  "       pressure = <path>   output SU file for writing in pressure field, must exist",
  "         source = <path>   input SU file containing source time function",
  "           bulk = 1.0      bulk modulus in GPa",
  "           buoy = 1.0      buoyancy in cm^3/g",
  "      dimension = 3        spatial dimension",
  "      MPS_ord_0 = 0        MPS order in z-direction, s_0",
  "      MPS_ord_1 = 0        MPS order in x-direction, s_1",
  "      MPS_ord_2 = 0        MPS order in y-direction, s_2",
  " ---------------------------end parameters ------------------------------",
  NULL };



IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"bulkmod",  0, true,  1 },
  {"buoyancy", 1, true,  1 },
  {"source_p", 2, true,  2 },
  {"data_p",   2, false, 2 },
  {"",         0, false, 0 }
};


float lininterp(int nt, float dt, float t0, float t, float *f) {
  float s = t/dt;
  int it = int(s);
  if (it < 0 || it > nt-2) return 0.0f;
  else return (s-float(it))*f[it]+(float(it)+1.0f-s)*f[it+1];
}

/*

*/

//----------------------------------------------------------------------//
void exact( int nt_out,     
	    int nt_src, 
	    float t0_out,   
	    float t0_src,
	    float dt, 
	    RPNT x, 
	    float c,
	    IPNT derv,
	    int dim,
	    float * f,
	    float * p ){
//----------------------------------------------------------------------//
  try{

    int d_ord = 0;
    for(int d=0; d<dim; d++) d_ord += derv[d];

    if( (dim!=2) && (dim!=3) ){
      RVLException e;
      e << "Spatial dimension ("<<dim<<") must be 2 or 3.\n";
      throw e;
    }

    if( dim==2 ){
      if( d_ord>1 ){
	RVLException e;
	e << "Derivative order must be <= 1 for dim="<<dim<<"\n"
	  << "Currently, d=";
	for(int d=0; d<dim; d++) e << derv[d]<<",";
	e << "\n";
	throw e;
      }      
    } 
    if( dim==3 ){
      if( d_ord>2 ){
	RVLException e;
	e << "Derivative order must be <= 2 for dim="<<dim<<"\n"
	  << "Currently, d=";
	for(int d=0; d<dim; d++) e << derv[d]<<",";
	e << "\n";
	throw e;
      }
    }
    
    float * df=NULL;
    float * ddf=NULL;
    float * dddf=NULL;

    //compute df
    df = new float[nt_src];
    df[0]=(f[1]-f[0])/dt;
    df[nt_src-1]=(f[nt_src-1]-f[nt_src-2])/dt;
    for (int i=1;i<nt_src-1;i++)
      df[i]=(f[i+1]-f[i-1])/(2.0f*dt);
    
    if(d_ord>0){
      //compute ddf
      ddf = new float[nt_src];
      ddf[0]=(df[1]-df[0])/dt;
      ddf[nt_src-1]=(df[nt_src-1]-df[nt_src-2])/dt;
      for (int i=1;i<nt_src-1;i++)
	ddf[i]=(df[i+1]-df[i-1])/(2.0f*dt);
    }
    if(d_ord>1){
      //compute dddf
      dddf = new float[nt_src];
      dddf[0]=(ddf[1]-ddf[0])/dt;
      dddf[nt_src-1]=(ddf[nt_src-1]-ddf[nt_src-2])/dt;
      for (int i=1;i<nt_src-1;i++)
	dddf[i]=(ddf[i+1]-ddf[i-1])/(2.0f*dt);
    }
   
    // r = |x|
    float r = 0.0f;
    for(int d=0; d<dim; d++) r += x[d]*x[d];
    r = sqrt(r);

    // gamma_i = x_i/r
    // gamma_ij = x_i*x_j/r^2
    float gamma = 1.0f;
    for(int d=0;d<dim;d++){
      for(int k=0;k<derv[d];k++)
	gamma *= x[d]/r;
    }

    // delta_ij
    float delta = 0.0f;
    bool tmp = false;
    for(int d=0; d<dim; d++) tmp = tmp || derv[d]==2;
    if(tmp) delta = 1.0f;

    float xoc = r/c; //traveltime
    float T_src = (nt_src-1)*dt + t0_src; //final source time
    float pi = 4.0f*atan(1.0);

    /////////////
    // 2D case //
    /////////////
    if( dim==2 ){
      //time loop
      for( int i=0; i<nt_out; i++ ){
	
	p[i] = 0.0f;
	float tau = i*dt - t0_out - xoc;

	if( tau>=0 ){
	  float lim = sqrt(tau);
	  int i_lim = int(lim/dt);
	  float a = 0.0f;

	  //integral over sigma
	  for( int k=0; k<=i_lim; k++ ){

	    float sigma2 = k*dt; 
	    sigma2*=sigma2;
	    float Omega_inv = sqrt( sigma2 + 2.0f*r/c );
	    float arg = tau - sigma2;
	    
	    if(d_ord==0){
	      a += lininterp(nt_src,dt,t0_src,arg,df)/Omega_inv;
	    }
	    else if(d_ord==1){
	      a += -(gamma/c) * lininterp(nt_src,dt,t0_src,arg,ddf)/Omega_inv 
		-(gamma/c) * lininterp(nt_src,dt,t0_src,arg,df)/(Omega_inv*Omega_inv*Omega_inv);
	    }
	  }
	  p[i] = dt*a/(pi*c*c);
	}//end if
      }
      
    }
    /////////////
    // 3D case //
    /////////////
    if( dim==3 ){
      //time loop
      for( int i=0; i<nt_out; i++ ) {
	
	p[i]=0.0f;
	float tau = i*dt - t0_out - xoc;
      
	// only nonzero for t-r/c in (t0_src,T_src)
	if( tau>=t0_src && tau<=T_src ) {
	  int j = int((tau+t0_src)/dt);
	  
	  if(d_ord==0){
	    p[i] = df[j];
	  }
	  else if(d_ord==1){
	    p[i] = -ddf[j]*gamma/c - df[j]*gamma/r;
	  }
	  else if(d_ord==2){
	    p[i] = dddf[j]*gamma/(c*c) 
	      - ddf[j]*(delta-3.0f*gamma)/(c*r)
	      - dddf[j]*(delta-3.0f*gamma)/(r*r);
	  }
	  p[i] /= 4.0f*pi*c*c*r;
	}
      }
    }

    // cleanup
    delete[] df;
    delete[] ddf;
    delete[] dddf;
    
  }
  catch(RVLException &e){
    e << "ERROR from exact!\n";
    throw e;
  }
}



int xargc;
char **xargv;

//----------------------------------------------------------------------//
int main(int argc, char ** argv) {
//----------------------------------------------------------------------//
  try {

#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);        
#endif

    PARARRAY *pars=NULL;
    FILE *stream=NULL;
    char *str=NULL;
    Value val;
    FILE *fp_src=NULL;
    FILE *fp_trc=NULL;
    float bulk, buoy;
    string src, trc;
    segy tr_src, tr_trc;
    int nt_src, nt_trc;
    float dt_src, dt_trc, dt;
    float delrt_src, delrt_trc;
    float *pbuf=NULL;
    float *gbuf=NULL;
    fpos_t pos_trc;
    RPNT x;
    IPNT derv;
    int freesurf;
    float surf_z;

    str = (char *)usermalloc_(128*sizeof(char));    
    xargc=argc; xargv=argv;
    TSOpt::IWaveEnvironment(argc, argv, 0, &pars, &stream);

#ifdef IWAVE_USE_MPI
    if(retrieveGlobalRank()==0){
#endif

      if(argc<2){
	pagedoc();
	exit(0);
      }

      cerr << "\n////////////////////////////////////////////////////////////\n"
	   << "Running asg_ansol.x\n";

#ifdef IWAVE_USE_MPI
    }
#endif


    //filenames
    src = valparse<string>(*pars,"source");
    trc = valparse<string>(*pars,"pressure");
    
    //spatial dimension
    int dim = valparse<int>(*pars,"dimension");

    //derivative information
    derv[0] = valparse<int>(*pars,"MPS_ord_0");
    derv[1] = valparse<int>(*pars,"MPS_ord_1");
    if(dim==3)
      derv[2] = valparse<int>(*pars,"MPS_ord_2");

    //location of free surface
    freesurf = valparse<int>(*pars,"freesurf",0);
    surf_z = valparse<float>(*pars,"surf_z",0.0);


    //extract nt, dt, delrt from source file
    if (!(fp_src=iwave_const_fopen(src.c_str(),"r",NULL,stderr))) {
      RVLException e;
      e<<"failed to open source trace file = "<<src<<"\n";
      throw e;
    }

    if (!(fvgettr(fp_src,&tr_src))) {
      RVLException e;
      e<<"failed to read first trace from source file = "<<src<<"\n";
      throw e;
    }
    iwave_fclose(fp_src);    
    
    strcpy(str,"ns");
    gethdval(&tr_src,str,&val);
    nt_src = vtoi(hdtype(str),val);

    strcpy(str,"dt");
    gethdval(&tr_src,str,&val);
    dt_src = vtof(hdtype(str),val);

    strcpy(str,"delrt");
    gethdval(&tr_src,str,&val);
    delrt_src = vtof(hdtype(str),val); //in ms


    //extract nt, dt, delrt from output file 
    if (!(fp_trc=iwave_const_fopen(trc.c_str(),"r",NULL,stderr))) {
      RVLException e;
      e<<"failed to open pressure trace file = "<<trc<<"\n";
      throw e;
    }

    if (!(fvgettr(fp_trc,&tr_trc))) {
      RVLException e;
      e<<"failed to read first trace from pressure file = "<<trc<<"\n";
      throw e;
    }
    iwave_fclose(fp_trc);

    strcpy(str,"ns");
    gethdval(&tr_trc,str,&val);
    nt_trc = vtoi(hdtype(str),val);

    strcpy(str,"dt");
    gethdval(&tr_trc,str,&val);
    dt_trc = vtof(hdtype(str),val); //in us

    strcpy(str,"delrt");
    gethdval(&tr_trc,str,&val);
    delrt_trc = vtof(hdtype(str),val); //in us
    
    //sanity check
    if(dt_src!=dt_trc){
      RVLException e;
      e << "dt for source and pressure files differ!\n"
	<< "dt_src="<<dt_src<<", dt_trc="<<dt_trc<<"\n";
      throw e;
    }

    //convert dt to ms 
    dt = (0.001)*dt_src;

    //computing sound speed
    bulk = valparse<float>(*pars,"bulk");
    buoy = valparse<float>(*pars,"buoy");
    float c = sqrt(bulk*buoy); 

    //extracting receiver positions and source position
    vector<Rtuple> g_pos = ex_pos(trc,true);
    vector<Rtuple> s_pos = ex_pos(trc,false);
   
    //data buffer for output pressure trace (and ghost)
    pbuf = (float *)usermalloc_(g_pos.size()*nt_trc*sizeof(float));    
    gbuf = (float *)usermalloc_(g_pos.size()*nt_trc*sizeof(float));    

    /////////////////////////
    // loop over receivers //
    /////////////////////////
    for( int i=0; i<g_pos.size(); i++){

      //computing x = x_r-x_s
      for(int d=0; d<dim; d++)
	x[d] = g_pos[i].coor[d] - s_pos[0].coor[d];

      //computing pressure trace
      exact( nt_trc, 
	     nt_src, 
	     delrt_trc,
	     delrt_src,
	     dt, x, c, derv, dim,
	     tr_src.data, 
	     pbuf+i*nt_trc );

      //computing ghost trace
      if( freesurf){

	//some sanity checks
	if( ( surf_z>s_pos[0].coor[0] && surf_z<g_pos[i].coor[0] ) || 
	    ( surf_z<s_pos[0].coor[0] && surf_z>g_pos[i].coor[0] ) ){
	  RVLException e;
	  e << "Surface depth ("<<surf_z<<") cannot be in between "
	    << "source depth ("<<s_pos[0].coor[0]<<") and receiver depth ("<<g_pos[i].coor[0]<<")\n";
	  throw e;
	}
	x[0] -= 2*(surf_z-s_pos[0].coor[0]);

	exact( nt_trc,
	       nt_src,
	       delrt_trc,
	       delrt_src,
	       dt, x, c, derv, dim,
	       tr_src.data,
	       gbuf+i*nt_trc );
      }
    }

    //copying buffer to trace and writing out
    if (!(fp_trc=iwave_const_fopen(trc.c_str(),"r+",NULL,stderr))) {
      RVLException e;
      e<<"failed to open pressure trace file = "<<trc<<" for writing\n";
      throw e;
    }

    for( int i_tr=0; i_tr<g_pos.size(); i_tr++ ){
      fgetpos(fp_trc,&pos_trc);
      fvgettr(fp_trc,&tr_trc);
      
      for( int i_t=0; i_t<nt_trc; i_t++ ){
	tr_trc.data[i_t] = pbuf[ i_t+i_tr*nt_trc ];
	if(freesurf){
	  tr_trc.data[i_t] -= gbuf[ i_t+i_tr*nt_trc ];
	}
      }
      fsetpos(fp_trc,&pos_trc);
      fputtr(fp_trc,&tr_trc);
    }
    fflush(fp_trc);
    iwave_fclose(fp_trc);
    
    //clean up
    userfree_(str);
    userfree_(pbuf);
    ps_delete(&pars);
    iwave_fdestroy();

#ifdef IWAVE_USE_MPI
    if(retrieveGlobalRank()==0){
#endif

      cerr << "\nFinishing asg_ansol.x\n"
	   << "////////////////////////////////////////////////////////////\n";

#ifdef IWAVE_USE_MPI
    }
    MPI_Finalize();
#endif

  }
  catch (RVL::RVLException & e) {
    e << "Exiting with error!\n";
#ifdef IWAVE_USE_MPI
    if( retrieveGlobalRank()==0 ){
#endif
      e.write(cerr);
#ifdef IWAVE_USE_MPI
    }
    MPI_Barrier(retrieveGlobalComm());
    MPI_Abort(retrieveGlobalComm(),0);
#endif
    exit(1);
  }
}
