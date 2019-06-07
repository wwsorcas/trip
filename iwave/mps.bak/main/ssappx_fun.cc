#include <except.hh>
#include <su.h>
#include <segy.h>
#include <header.h>
#include <parser.h>

#include "MPS_includes.hh"
#include "ssappx.hh"


const char * sdoc[] = { 
  "Usage: ssappx_fun.x h= a_ord= d_ord= ox= dx= nx= filename=",
NULL};  

using TSOpt::ssappx;
using RVL::valparse;

int xargc;
char **xargv;

//-----------------------------------------------------------------------------------
int main(int argc, char ** argv) {
//-----------------------------------------------------------------------------------
  try {
 
    xargc=argc; xargv=argv;
    requestdoc(1);

    PARARRAY * par = ps_new();
    if ( ps_createargs(par, argc - 1, argv + 1) ) {
      printf("Error parsing input data. ABORT.\n");
      exit(1);
    }
    //ps_printall(*par,stderr);

    //parameters    
    float h = valparse<float>(*par,"h");
    int a_ord = valparse<int>(*par,"a_ord");
    string filename = valparse<string>(*par,"filename");

    int s0 = valparse<int>(*par,"s0",0);
    int s1 = valparse<int>(*par,"s1",0);

    float o_x0 = valparse<float>(*par,"o_x0",0);
    float d_x0 = valparse<float>(*par,"d_x0",1);
    float f_x0 = valparse<float>(*par,"f_x0",0);
    int n_x0 = (int)((f_x0-o_x0)/d_x0)+1;

    float o_x1 = valparse<float>(*par,"o_x1",0);
    float d_x1 = valparse<float>(*par,"d_x1",1);
    float f_x1 = valparse<float>(*par,"f_x1",0);
    int n_x1 = (int)((f_x1-o_x1)/d_x1)+1;

    ofstream ofs;
    ofs.open (filename, ios::out);

    if(h<=0){
      RVLException e;
      e << "grid size h="<<h<<" has to be >0\n";
      throw e;
    }
    if(a_ord<1){
      RVLException e;
      e << "approximation order a_ord="<<a_ord<<" has to be >=1\n";
      throw e;
    }
    if(s0<0||s1<0){
      RVLException e;
      e << "derivative order s=("<<s0<<","<<s1<<" has to be >=0\n";
      throw e;
    }


    float x0, x1;
    float W;

    for( int i0=0; i0<n_x0; i0++){
      
      for( int i1=0; i1<n_x1; i1++){
	x0 = o_x0+ i0*d_x0;	
	W = ssappx( x0, h, s0, a_ord );

	x1 = o_x1+ i1*d_x1;	
	W *= ssappx( x1, h, s1, a_ord );

	ofs << W<<" ";
	//cerr << "x="<<x<<", y="<<y<<"\n";
      }
      ofs << "\n";
    }

    ofs.close();

  }
  catch (RVLException & e) {
    e << "ERROR from ssappx_fun.x!\n";
    e.write(cerr);
    exit(1);
  }
}
