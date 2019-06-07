// ansol_wave.cc
// Author: Mario J. Bencomo
// last modification: 11/30/16

//#include "asg_defn.hh"
//#include "asg.hh"
#include "MPS_includes.hh"
#include "MPS_Space.hh"
#include "my_waveform.hh"

using RVL::valparse;
using RVL::RVLException;
using TSOpt::Rtuple;
using TSOpt::Ituple;
using TSOpt::ex_pos;

const char *sdoc[] = {
  " ============================================================================",
  " ansol_wave.x ",
  " ============================================================================",
  " Authors: Mario J. Bencomo",
  " ",
  " Analytical solution for acoustic systems, either 1st or 2nd order form, in a", 
  " homogeneous unbounded domain, ony 3-D.",
  " ",
  " 2nd order PDE system: ",
  "     d^2p/dt^2 - c^2*nabla^2(p) = f(x,t)",
  " ",
  " 1st order PDE system: ",
  "     dp/dt - kappa*div(v) = f(x,t)",
  "     dv/dt - buoy*grad(p) = 0",
  " ",
  " where ",
  "     c = velocity, ",
  " kappa = bulk modulus, ",
  "  buoy = buoyancy, "
  "     p = pressure field, ",
  "     v = velocity field, ",
  "     f = source term ",
  "       = w(t) (d/dx_0)^{s_0} (d/dx_1)^{s_1} (d/dx_2)^{s_2} delta(x-x^*)",
  "",
  " Source waveform w(t) is assumed to be a scaled gaussian or some derivative thereof."
  " Output is p(x_r,t) for given reciever locations x_r specified by input SU file. ",
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
  "       pressure = <path>   Output SU file for writing in pressure field, ",
  "                           must exist.",
  "",
  "        src_type = 0       Type of source wavelet,",
  "                           (0=gaussian,1=dgaussian,2=ricker,3=dddgaussian)",
  "       src_fpeak =         Source wavelet peak frequency, in Hertz",
  "         src_off = 0       Source wavelet time offset/center, in seconds."
  "        src_scal = 1       Scaling factor of source wavelet."
  "",
  "         PDE_ord = 2       Order of PDE system,",
  "                           1=first order system, this adds an extra time ",
  "                             derivative on the source wavelet in the ",
  "                             analytical solutions",
  "                           2=second order system",
  "",
  "              c = 1.0      velocity in km/s",
  "            s_0 = 0        MPS order in z-direction",
  "            s_1 = 0        MPS order in x-direction",
  "            s_2 = 0        MPS order in y-direction",
  " ---------------------------end parameters ------------------------------",
  NULL };


//----------------------------------------------------------------------//
void ansol( vector<float> t,   
	    RPNT  x, 
	    float c,
	    IPNT  s,
	    int   PDE_ord,
	    int   src_type,
	    float src_fpeak,
	    float src_off,
	    float src_scal,
	    float *p ){
//----------------------------------------------------------------------//
  try{
   
    //cerr << "Inside ansol():\n";

    //s_ord=|s|
    int s_ord=0;
    for( int d=0; d<3; d++ )
      s_ord += s[d];
    cerr << "Inside ansol, |s|="<<s_ord<<"\n";
    
    // r = |x|
    float r = 0.0f;
    for( int d=0; d<3; d++ ) 
      r += x[d]*x[d];
    r = sqrt(r);

    // gamma = x_i/r for dipole
    // gamma = x_i*x_j/r^2 for quadrapole
    float gamma = 1.0f;
    for( int d=0; d<3; d++ ){
      for( int k=0; k<s[d]; k++ )
	gamma *= x[d]/r;
    }

    //kronocker delta
    int k_delta = 0;
    if(s_ord==2){
      for( int d=0; d<3; d++ ) 
	if(s[d]==2) k_delta=1;
    }

    float tau;
    float tmp = 4*PI*c*c*r;
    int type = 2 - PDE_ord + src_type;

    if( (type+s_ord)>3 || (type+s_ord)<0 ){
      RVLException e;
      e << "Total waveform type must be >=0 and <=3! type="<<type<<", |s|="<<s_ord<<"\n";
      throw e;
    }

    //monopole case
    if( s_ord==0 ){
      for( int i=0; i<t.size(); i++ ) {
	tau   = t[i] - r/c - src_off;
        p[i]  = src_scal*my_waveform(type,tau,src_fpeak);
	p[i] /= tmp;
      }
    }
    //dipole case
    else if( s_ord==1 ){
      for( int i=0; i<t.size(); i++ ) {
	tau   = t[i] - r/c - src_off;
	p[i]  = gamma/c * src_scal*my_waveform(type+1,tau,src_fpeak);
	p[i] += gamma/r * src_scal*my_waveform(type,tau,src_fpeak);
	p[i] /= -tmp;

      }
    }
    //quadrapole case
    else if( s_ord==2 ){
      for( int i=0; i<t.size(); i++ ) {
      	tau   = t[i] - r/c - src_off;
	p[i]  = gamma/(c*c) * src_scal*my_waveform(type+2,tau,src_fpeak);
	p[i] += -(k_delta-3*gamma)/(c*r) * src_scal*my_waveform(type+1,tau,src_fpeak);
	p[i] += -(k_delta-3*gamma)/(r*r) * src_scal*my_waveform(type,tau,src_fpeak);
	p[i] /= tmp;
      }
    }
  }
  catch(RVLException &e){
    e << "ERROR from ansol!\n";
    throw e;
  }
}



int xargc;
char **xargv;

//----------------------------------------------------------------------//
int main(int argc, char ** argv) {
//----------------------------------------------------------------------//
  try {

    FILE  *fp=NULL;   //output file pointer
    string filename;  //output filename 

    int    src_type;  //source wavelet type
    float  src_fpeak; //source wavelet peak frequency
    float  src_off;   //source wavelet offset/center
    float  src_scal;  //source wavelet scaling factor
    float  c;         //velocity
    int    PDE_ord;   //asg key
    IPNT   s;         //multipole derivative multi-index

    segy   trc;       //output segy trace
    int    nt;        //number of time samples in output trace
    float  dt;        //time step size in output trace
    float  t0;        //initial time
    vector<float> t;  //time axis
    fpos_t pos;       //position of trace
    RPNT   x;         //R3 vector given by x_r-x_s

    char  *str;
    Value  val;

    str = (char *)usermalloc_(128*sizeof(char));    

    //reading in parameters
    PARARRAY *par = ps_new();
    if (ps_createargs(par,argc-1,argv+1)) {
      fprintf(stderr,
	      "ERROR. could not process command line\n");
      exit(1);
    }


    requestdoc(0);

    cerr << "\n////////////////////////////////////////////////////////////\n"
	 << "Running ansol_wave.x\n";

    //output filename
    filename = valparse<string>(*par,"pressure");

    //source info
    src_type  = valparse<int>(*par,"src_type",0);
    src_fpeak = valparse<float>(*par,"src_fpeak");
    src_off   = valparse<float>(*par,"src_off",0.0); //already in seconds
    src_scal  = valparse<float>(*par,"src_scal",1.0);

    //PDE order
    PDE_ord = valparse<int>(*par,"PDE_ord",2);

    //medium velocity
    c = valparse<float>(*par,"c",1.0);
    c *= 1000; //converting from km/s to m/s

    //derivative information
    s[0] = valparse<int>(*par,"s_0",0);
    s[1] = valparse<int>(*par,"s_1",0);
    s[2] = valparse<int>(*par,"s_2",0);

    /*
    cerr <<"Printing out parameters:\n"
         << "  pressure = "<<filename<<"\n"
	 << "  src_type = "<<src_type<<"\n"
	 << "   src_off = "<<src_off<<"\n"
	 << "  src_scal = "<<src_scal<<"\n"
	 << "   PDE_ord = "<<PDE_ord<<"\n"
	 << "         c = "<<c<<"\n"
	 << "         s = ["<<s[0]<<","<<s[1]<<","<<s[2]<<"]\n";
    */

    //extract nt, dt, t0 from output file 
    if (!(fp=iwave_const_fopen(filename.c_str(),"r",NULL,stderr))) {
      RVLException e;
      e<<"failed to open pressure trace file = "<< filename <<"\n";
      throw e;
    }

    if (!(fvgettr(fp,&trc))) {
      RVLException e;
      e<<"failed to read first trace from pressure file = "<< filename <<"\n";
      throw e;
    }
    iwave_fclose(fp);

    strcpy(str,"ns");
    gethdval(&trc,str,&val);
    nt = vtoi(hdtype(str),val);

    strcpy(str,"dt");
    gethdval(&trc,str,&val);
    dt = vtof(hdtype(str),val); //in us
    dt *= 1e-6; //converting to seconds

    strcpy(str,"delrt");
    gethdval(&trc,str,&val);
    t0 = vtof(hdtype(str),val); //in us
    t0 *= 1e-6; //converting to seconds

    //setting time axis
    t.resize(nt);
    for( int i=0; i<nt; i++ )
      t[i] = t0 + i*dt;

    //extracting receiver positions and source position
    vector<Rtuple> g_pos = ex_pos(filename,true);
    vector<Rtuple> s_pos = ex_pos(filename,false);


    if (!(fp=iwave_const_fopen(filename.c_str(),"r+",NULL,stderr))) {
      RVLException e;
      e<<"failed to open pressure trace file = "<< filename <<" for writing\n";
      throw e;
    }    

    /////////////////////////
    // loop over receivers //
    /////////////////////////
    for( int i=0; i<g_pos.size(); i++){

      fgetpos(fp,&pos);
      fvgettr(fp,&trc);
      
      //computing x = x_r-x_s
      for(int d=0; d<3; d++)
	x[d] = g_pos[i].coor[d] - s_pos[0].coor[d];
      
      //computing pressure trace
      ansol( t, x, c, s,
	     PDE_ord,
	     src_type,
	     src_fpeak,
	     src_off,
	     src_scal,
	     trc.data );

      fsetpos(fp,&pos);
      fputtr(fp,&trc);      

    }
    fflush(fp);
    iwave_fclose(fp);
    
    //clean up
    userfree_(str);
    ps_delete(&par);
    iwave_fdestroy();

    cerr << "\nFinishing ansol_wave.x\n"
         << "////////////////////////////////////////////////////////////\n";

  }
  catch (RVL::RVLException & e) {
    e << "Exiting with error!\n";
    e.write(cerr);
    exit(1);
  }
}
