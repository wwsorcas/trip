#include <except.hh>
#include <su.h>
#include <segy.h>
#include <header.h>
#include <parser.h>
#include <math.h>

const char * sdoc[] = { 
  "Usage: my_diffint.x ",
  " ",
  "Required parameters:",
  "  output = filename of output wavelet",
  "  order  = order of diffint operation",
NULL};  

using RVL::RVLException;
using RVL::valparse;
using std::string;

int xargc;
char **xargv;


//--------------------------------------------------------------------
float diffint(float dt, int i, float *data, int order){
//--------------------------------------------------------------------
  float ans = 0;
  if(order==-1){
    for(int j=0; j<=i; j++){
      ans += data[j];
    }
    ans *= dt;
  }
  else if(order==1){
    ans = data[i];
    if(i>0) ans += -data[i-1];
    ans /= dt;
    /*
    ans = 3*data[i];
    if(i>0) ans += -4*data[i-1];
    if(i>1) ans += data[i-2];
    ans /= 2*dt;
    */
  }
  else if(order==2){
    
    ans = data[i];
    if(i>0) ans += -2*data[i-1];
    if(i>1) ans += data[i-2];
    ans /= dt*dt;
    /*
    ans = 9*data[i];
    if(i>0) ans += -24*data[i-1];
    if(i>1) ans += 22*data[i-2];
    if(i>2) ans += -8*data[i-3];
    if(i>3) ans += data[i-4];
    ans /= 4*dt*dt;
    */
  }

  return ans;
}

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

    string input  = valparse<string>(*par,"input");
    string output = valparse<string>(*par,"output");
    int    order  = valparse<int>(*par,"order");

#ifdef VERBOSE_MJB
    cerr << "Inside my_diffint.x\n"
	 << "  input="<<input<<"\n"
	 << "  output="<<output<<"\n"
	 << "  order="<<order<<"\n";
#endif

    FILE *fp_in, *fp_out;
    segy tr;
    float dt;
    int nt;
    float *buff = NULL;
    
    fp_in  = fopen(input.c_str(),"r"); 
    fp_out = fopen(output.c_str(),"w"); 

    fgettr(fp_in,&tr);
    dt = tr.dt*(1.e-6); // convert units, from [mus] to [s]
    nt = tr.ns;

    buff = new float[nt];

    do{
      
      for( int i=0; i<nt; i++ )
	buff[i] = diffint(dt,i,tr.data,order);
      for( int i=0; i<nt; i++ )
	tr.data[i] = buff[i];

      fputtr(fp_out,&tr);
    }
    while(fgettr(fp_in,&tr));
      
    delete[] buff;
    fclose(fp_in);
    fclose(fp_out);

  }
  catch (RVLException & e) {
    e << "Error caught in main for my_diffint.cc\n";
    e.write(cerr);
    exit(1);
  }
}
