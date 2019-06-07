#include <except.hh>
#include <su.h>
#include <segy.h>
#include <header.h>
#include <parser.h>

#define VERBOSE

const char * sdoc[] = {
  "This executable concatenates replicas of itself with itself.",
  " ",
  "Sample usage: replicate.x ",
  "       CWPROOT=$CWPROOT ",
  "       in=input.su ",
  "       out=output.su ",
  "       N=2 ",
  " ",
  "Required parameters:",
  "  CWPROOT = path to SU root directory",
  "  in      = input filename",
  "  out     = output filename",
  "  N       = number replicas",
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
    
    cerr << "Inside replicate:\n";
    ps_printall(*par,stderr);

    //Reading in parameters from command line
    std::string cwp = valparse<std::string>(*par,"CWPROOT");
    
    int N = valparse<int>(*par,"N");
    std::string in  = valparse<string>(*par,"in");
    std::string out = valparse<string>(*par,"out");

    //set up paths to SU commands
    std::string sushw = cwp + "/bin/sushw";
    std::string suresamp = cwp + "/bin/suresamp";

    //concatenating
    for(int n=0; n<N; n++){
      std::stringstream cmd;
      cmd << suresamp <<" <"<< in <<" >> tmp.su";

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
  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
