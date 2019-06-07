#include <except.hh>
#include <su.h>
#include <segy.h>
#include <header.h>
#include <parser.h>

#define VERBOSE

const char * sdoc[] = { 
  "Usage: tg_setup.x ",
  "       CWPROOT=$CWPROOT "
  "       N=2 ",
  "       wlt_0=input1.su ",
  "       wlt_1=input2.su ",
  " ",
  "Required parameters:",
  "  CWPROOT = path to SU root directory",
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
    std::string *wlts    = new std::string[N];
    std::string *wlts_cp = new std::string[N];
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
      ss << "wlt_" << n;
      std::string wlt_i = ss.str();
      wlts[n] = valparse<std::string>(*par,wlt_i);
      wlts_cp[n] = "cp_" + wlts[n];

      std::stringstream move;
      move << "/bin/mv "<< wlts[n] <<" "<< wlts_cp[n];
      system(move.str().c_str());

#ifdef VERBOSE
      move << move.str() << endl;
#endif

      FILE* fp;
      segy tr;

      fp = fopen(wlts_cp[n].c_str(),"r"); 
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
    std::string suplane = cwp + "/bin/suplane";
    std::string sushw   = cwp + "/bin/sushw";
    std::string suconv  = cwp + "/bin/suconv";
    std::string sugain  = cwp + "/bin/sugain";
    std::string suresamp= cwp + "/bin/suresamp";

    //resampling wavelets
    for(int n=0; n<N; n++){
      std::stringstream cmd;
      cmd << suresamp <<" <"<< wlts_cp[n] <<" >"<< wlts[n];
      cmd << " tmin="<< min_ot <<" dt="<< max_dt <<" nt="<< tot_nt;
  
#ifdef VERBOSE
      cerr << cmd.str() << endl;
#endif
      
      system(cmd.str().c_str());
    
      //clean up wavelets_re
      std::stringstream clean;
      clean << "/bin/rm "<< wlts_cp[n];
      system(clean.str().c_str());
    }

    delete [] wlts;
    delete [] wlts_cp;
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
