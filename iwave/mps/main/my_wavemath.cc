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
  "  in1 = left input SU file",
  "  in2 = right input SU file",
  "  out = output SU file",
  "  op  = operator flag (pointwise operations),",
  "       options=add,sub,mul,div"
  "  a=1  left constant in add and sub operations",
  "  b=1  right constant in add and sub operations",
  "",
  "Available operations (pointwise):",
  "  add: out[i] = a*in1[i] + b*in2[i]",
  "  sub: out[i] = a*in1[i] - b*in2[i]",
  "  mul: out[i] = in1[i] * in2[i]",
  "  div: out[i] = in1[i] / in2[i]",
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

    string in1    = valparse<string>(*par,"in1");
    string in2    = valparse<string>(*par,"in2");
    string out    = valparse<string>(*par,"out");
    string op     = valparse<string>(*par,"op");
    float  a      = valparse<float>(*par,"a",1.0f);
    float  b      = valparse<float>(*par,"b",1.0f);

    //reading in in1
    FILE* fp_in1=NULL;
    segy tr_in1;
    
    if(!(fp_in1 = fopen(in1.c_str(),"r"))){
      RVLException e;
      e << "Could not open in1="<<in1<<"\n";
      throw e;
    }
    fvgettr(fp_in1,&tr_in1);

    int   nt_in1    = tr_in1.ns;
    float dt_in1    = tr_in1.dt;
    float delrt_in1 = tr_in1.delrt;

    //reading in in2
    FILE* fp_in2=NULL;
    segy tr_in2;
    if(!(fp_in2 = fopen(in2.c_str(),"r"))){
      RVLException e;
      e << "Could not open in2="<<in2<<"\n";
      throw e;      
    }
    fvgettr(fp_in2,&tr_in2);
    
    int   nt_in2    = tr_in2.ns;
    float dt_in2    = tr_in2.dt;
    float delrt_in2 = tr_in2.delrt;

    //checking for compatability
    if( nt_in1!=nt_in2 || 
	dt_in1!=dt_in2 || 
	delrt_in1!=delrt_in2 ){
      RVLException e;
      e << "Input traces are not compatible!\n"
	<< "in1="<<in1<<"\n"
	<< " nt="<<nt_in1<<"\n"
	<< " dt="<<dt_in1<<"\n"
	<< " delrt="<<delrt_in1<<"\n"
	<< "in2="<<in2<<"\n"
	<< " nt="<<nt_in2<<"\n"
	<< " dt="<<dt_in2<<"\n"
	<< " delrt="<<delrt_in2<<"\n";
      throw e;
    }

    //creating output su file
    stringstream cmd;
    cmd << "/bin/cp "<<in1<<" "<<out;
    system(cmd.str().c_str());

    //Writing output
    FILE* fp_out=NULL;
    if(!(fp_out = fopen(out.c_str(),"w"))){
      RVLException e;
      e << "Could not open out="<<out<<"\n";
      throw e;     
    }

    do{
      for( int it=0; it<nt_in1; it++ ){
	if(op.compare("add")==0){
	  tr_in1.data[it] = tr_in1.data[it]*a + tr_in2.data[it]*b;
	}
	else if(op.compare("sub")==0){
	  tr_in1.data[it] = tr_in1.data[it]*a - tr_in2.data[it]*b;
	}
	else if(op.compare("mul")==0){
	  tr_in1.data[it] *= tr_in2.data[it];
	}
	else if(op.compare("div")==0){
	  tr_in1.data[it] /= tr_in2.data[it];
	}
      }
      fvputtr(fp_out,&tr_in1);      
    }
    while(fvgettr(fp_in1,&tr_in1)&&
	  fvgettr(fp_in2,&tr_in2));

    fclose(fp_in1);
    fclose(fp_in2);
    fclose(fp_out);
  }
  catch (RVLException & e) {
    e << "Error caught in main for my_wavemath.cc\n";
    e.write(cerr);
    exit(1);
  }
}
