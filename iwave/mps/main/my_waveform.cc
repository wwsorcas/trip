#include <my_waveform.hh>
#include <except.hh>
#include <su.h>
#include <segy.h>
#include <header.h>
#include <parser.h>
//#include <math.h>

const char * sdoc[] = { 
  "Usage: my_wlt.x CWPROOT= wavelet= dt=",
  " ",
  "Required parameters:",
  "  CWPROOT  = path to SU root directory",
  "  output   = filename of output wavelet",
  "  type     = waveform type ",
  "             (0=gauss,1=dgauss,2=ricker,3=dddgauss,-1=cinf)",
  "  nt       = number of time grid",
  "  dt       = time increment",
  "  cit      = time index of center",
  "  scal     = scaling factor",
  "  width    = width of cinf support",
  "  fpeak    = peak freq for gauss,dgauss,ricker",
NULL};  

using RVL::RVLException;
using RVL::valparse;
using std::string;

int xargc;
char **xargv;


//--------------------------------------------------------------------
int main(int argc, char ** argv) {
//--------------------------------------------------------------------
  try {

    xargc=argc; xargv=argv;
    requestdoc(1);

    PARARRAY * par = ps_new();
    if ( ps_createargs(par, argc - 1, argv + 1) ) {
      printf("Error parsing input data. ABORT.\n");
      exit(1);
    }
    //ps_printall(*par,stderr);

    string cwp    = valparse<string>(*par,"CWPROOT");
    string output = valparse<string>(*par,"output");
    int    type   = valparse<int>(*par,"type");
    int    nt     = valparse<int>(*par,"nt");  
    float  dt     = valparse<float>(*par,"dt"); //units [s]

    //scaling factor
    float scal = valparse<float>(*par,"scal",1.0);    

    //center time index,
    // required for type=0,1,2,3,4
    int cit = valparse<int>(*par,"cit",0);

    //peak frequency,
    // required for type=0,1,2,3,6
    int fpeak = valparse<int>(*par,"fpeak",0); //units [Hz]
    if( fpeak==0 && 
	(type==0 || type==1 || type==2 || type==3 || type==6)){
      RVLException e;
      e << "fpeak is required for type="<< type <<"\n";
      throw;
    }
    
    //width of support for type=-1
    //width of ramp for type=4
    float width = valparse<float>(*par,"width",0.1); //units [s]
    // int it0 = cit - (int)(width*0.5/dt)-1;
    // int itF = cit + (int)(width*0.5/dt)+1;

    //bottom and top of ramp,
    // required for type=4
    float bot = valparse<float>(*par,"ramp_bot",0.0f);
    float top = valparse<float>(*par,"ramp_top",1.0f);

    //creating su file    
    string sunull = cwp + "/bin/sunull";
    stringstream cmd;
    cmd << sunull << " ntr=1 nt="<< nt <<" dt="<< dt
	<< " >" << output; 
    system(cmd.str().c_str());

    //reading in trace
    FILE* fp;
    segy tr;
   
    fp = fopen(output.c_str(),"r"); 
    fvgettr(fp,&tr);
    fclose(fp);

    fp = fopen(output.c_str(),"w");

    //writing out traces
    float t;
    float val;

    for( int it=0; it<nt; it++ ){
      
      t = (it-cit)*dt;

      if(type==-1)
	val=my_cinf(t,width);
      else if(type==0)
	val=my_gauss(t,fpeak);  
      else if(type==1)
	val=my_dgauss(t,fpeak);
      else if(type==2)
	val=my_ricker(t,fpeak);
      else if(type==3)
	val=my_dddgauss(t,fpeak);
      else if(type==4)
	val=my_ramp(t,width,bot,top);
      else if(type==5)
	val=my_cinf_ramp(t,width,bot,top);
      else if(type==6)
	val=my_wpack(t,width,fpeak);
      else{
	RVLException e;
	e << "type must be between -1 and 5!\n";
	throw e;
      } 
      tr.data[it]=val*scal;
    }

    
    /*
    if(type==-1){
      for(int it=0; it<nt; it++) 
	tr.data[it]=my_cinf(dt,it,it0,itF)*scal;
    }
    else
    if(type==0){
      for(int it=0; it<nt; it++) 
	tr.data[it]=my_gauss(dt,it,cit,fpeak)*scal;
    } 
    else
    if(type==1){
      for(int it=0; it<nt; it++) 
	tr.data[it]=my_dgauss(dt,it,cit,fpeak)*scal;
    }
    else
    if(type==2){
      for(int it=0; it<nt; it++) 
	tr.data[it]=my_ricker(dt,it,cit,fpeak)*scal;
    }
    else
    if(type==3){
      for(int it=0; it<nt; it++) 
	tr.data[it]=my_dddgauss(dt,it,cit,fpeak)*scal;
    }
    else{
      RVLException e;
      e << "type must be between -1 and 3!\n";
      throw e;
    }
    */

    fvputtr(fp,&tr);
    fclose(fp);
  }
  catch (RVLException & e) {
    e << "Error caught in main for my_waveform.cc\n";
    e.write(cerr);
    exit(1);
  }
}
