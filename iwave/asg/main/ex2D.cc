#include <su.h>
#include <segy.h>
#include <header.h>
#include <parser.h>

extern void exact(int nt, float dt, float r,
		  float bulk, float buoy,
		  float * f,
		  float * p);

int xargc;
char **xargv;

const char * sdoc[] = {
  "ex2D.x: computes one trace of a 2D point source acoustic field in homogeneous",
  "fluid. Field equations are ",
  " D_t p + k (D_x v_x + D_z v_z) = w(t) delta(x) delta(z)",
  " D_t v_x + b (D_x p) = 0",
  " D_t v_z + b (D_z p) = 0",
  "Usage:",
  " ex2D.x bulk=<float> buoy=<float> distance=<float> source=<string> pressure=<string>",
  "where",
  " bulk = bulk modulus",
  " buoy = buoyancy = 1/density",
  " distance = source-receiver distance",
  " source = w(t) = input source wavelet",
  " pressure = output pressure trace at given distance from source point",
  " ",
  "Notes:",
  " 1. pressure trace is created, same nt and dt as source",
  " 2. scaling as above, i.e. source is NOT scaled by bulk modulus.",
  "This choice of scale is consistent with the asg modeling code sim.x,",
  "that is, if the wavelet trace is input as source_p to sim.x, should",
  "obtain a close approximation to the output of this utility.",
  NULL};
  
int main(int argc, char ** argv) {

  try {

    xargc=argc; xargv=argv;
    requestdoc(1);

    PARARRAY * pars = ps_new();
    if (ps_createargs(pars,argc-1,&(argv[1]))) {
      RVL::RVLException e;
      e<<"ERROR: poisson from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++) 
	e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }

    // fish constants out of command line
    float bulk = RVL::valparse<float>(*pars,"bulk");
    float buoy = RVL::valparse<float>(*pars,"buoy");
    float distance = RVL::valparse<float>(*pars,"distance");

    // filenames for source and pressure trace
    std::string src = RVL::valparse<std::string>(*pars,"source");
    std::string trc = RVL::valparse<std::string>(*pars,"pressure");

    // open source file
    FILE * fpsrc = NULL;
    if (!(fpsrc=fopen(src.c_str(),"r"))) {
      RVL::RVLException e;
      e<<"ERROR: poisson from fopen\n";
      e<<"  failed to open source trace file = "<<src<<"\n";
      throw e;
    }

    // open output trace file
    FILE * fptrc = NULL;
    if (!(fptrc=fopen(trc.c_str(),"w"))) {
      RVL::RVLException e;
      e<<"ERROR: poisson from fopen\n";
      e<<"  failed to open pressure trace file = "<<trc<<"\n";
      throw e;
    }

    // extract nt, dt, data samples from src
    segy tr;
    if (!(fgettr(fpsrc,&tr))) {
      RVL::RVLException e;
      e<<"ERROR: poisson from fgettr\n";
      e<<"  failed to read first trace from source file = "<<src<<"\n";
      throw e;
    }
    fclose(fpsrc);

    char * str = (char *)usermalloc_(128*sizeof(char));
    Value val;
    strcpy(str,"ns");
    gethdval(&tr,str,&val);
    int nt = vtoi(hdtype(str),val);

    //cerr<<"1\n";

    strcpy(str,"dt");
    gethdval(&tr,str,&val);
    float dt = vtof(hdtype(str),val);

    //cerr<<"2\n";
    
    // convert dt to ms
    dt *= 0.001;
    
    // data buffer for output pressure trace
    float * pbuf = (float *)usermalloc_(nt*sizeof(float));

    //cerr<<"3\n";    
    // compute trace, copy to trace data buffer
    exact(nt, dt, distance, bulk, buoy,
	  tr.data, pbuf);
    // scale by 1/bulk to remove assumed scaling in function
    for (int i=0; i<nt; i++) (tr.data)[i] = pbuf[i]/bulk;
    //cerr<<"4\n";
    // write trace to file
    fputtr(fptrc,&tr);

    //cerr<<"5\n";
    
    // clean up
    fclose(fptrc);
    userfree_(str);
    userfree_(pbuf);

    exit(0);
  }
  catch (RVL::RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
