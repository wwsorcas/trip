#include <except.hh>
#include <su.h>
#include <segy.h>
#include <header.h>
#include <parser.h>

#define VERBOSE

const char * sdoc[] = {
  "This executable concatenates several *.su files into a single *.su file",
  "with a consistent time grid for all.",
  " ",
  "Sample usage: tgs_n_cat.x ",
  "       CWPROOT=$CWPROOT ",
  "       out= output.su ",
  "       N=2 ",
  "       wlt_0=input1.su ",
  "       wlt_1=input2.su ",
  " ",
  "Required parameters:",
  "  CWPROOT = path to SU root directory",
  "  out     = output filename",
  "  N       = number of input wavelets",
  "  wlt_0   = filename for 1st input source wavelet",
  "  wlt_1   = filename for 2nd input source wavelet",
  "  ... ",
NULL};  

using RVL::RVLException;
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

    //Reading in parameters from command line
    std::string cwp = valparse<std::string>(*par,"CWPROOT");
    
    int N = valparse<int>(*par,"N");
    std::string out = valparse<string>(*par,"out");

    std::string *wlts = new std::string[N];
    float *ot = new float[N];
    float *dt = new float[N];
    float *et = new float[N];
    int   *nt = new int[N];
    
    //largest dt, smallest ot, and largest et
    float max_dt = -99999;
    float min_ot =  99999;
    float max_et = -99999;
    
    for(int n=0; n<N; n++){
      std::stringstream ss;
      ss << "wlt_" << n ;
      std::string wlt_i = ss.str();
      wlts[n] = valparse<std::string>(*par,wlt_i);

      FILE* fp;
      segy tr;

      fp = fopen(wlts[n].c_str(),"r"); 
      fvgettr(fp,&tr);

      ot[n] = tr.delrt*(1.e-3); //s
      dt[n] = tr.dt*(1.e-6); //convert from mus to s
      nt[n] = tr.ns;
      et[n] = ot[n] + nt[n]*dt[n]; //s
      
      if(dt[n]>max_dt) max_dt = dt[n];
      if(ot[n]<min_ot) min_ot = ot[n];
      if(et[n]>max_et) max_et = et[n];

      fclose(fp);
    }
    int tot_nt = (int)((max_et-min_ot)/max_dt);
    
    //set up paths to SU commands
    std::string sushw   = cwp + "/bin/sushw";
    std::string suresamp= cwp + "/bin/suresamp";

    //resampling wavelets
    for(int n=0; n<N; n++){

      std::stringstream cmd;
      cmd << suresamp <<" <"<< wlts[n] <<" >> tmp.su"
	  <<" tmin="<< min_ot <<" dt="<< max_dt <<" nt="<< tot_nt;
      
#ifdef VERBOSE
      cerr << cmd.str() << endl;
#endif
    
      system(cmd.str().c_str());
    }

    //numbering traces for final output su
    std::stringstream cmd;
    cmd << sushw << " <tmp.su key=tracr a=1 b=1 | "
	<< sushw << " key=tracl a=1 b=1 >"<< out 
	<< " && /bin/rm tmp.su";

#ifdef VERBOSE
      cerr << cmd.str() << endl;
#endif

    system(cmd.str().c_str());

    delete [] wlts;
    delete [] ot;
    delete [] nt;
    delete [] et;
    delete [] dt;
  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
