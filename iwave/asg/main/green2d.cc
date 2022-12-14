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
  "green2d.x: computes one trace of a 2D point source acoustic green's function in homogeneous",
  "fluid. Field equations are ",
  " D_t p + k (D_x v_x + D_z v_z) = H(t) delta(x) delta(z)",
  " D_t v_x + b (D_x p) = 0",
  " D_t v_z + b (D_z p) = 0",
  "Usage:",
  " green2d.x bulk=<float> buoy=<float> out=<string>",
  "where",
  " bulk = bulk modulus",
  " buoy = buoyancy = 1/density",
  " out  = output trace file",
  " ",
  "Notes:",
  " source-receiver distance computed from trace headers gelev, selev,",
  " gx, sx, times from nt, dt, all read from output file. Therefore output",
  " file must contain correct trace geometry before this command is invoked.",
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

    // filename
    std::string out = RVL::valparse<std::string>(*pars,"out");

    // open source file
    FILE * fpout = NULL;
    if (!(fpout=fopen(out.c_str(),"r+"))) {
      RVL::RVLException e;
      e<<"ERROR: green2d from fopen\n";
      e<<"  failed to open output trace file = "<<out<<"\n";
      throw e;
    }

    // extract nt, dt, data samples from src
    segy tr;
    if (!(fgettr(fpout,&tr))) {
      RVL::RVLException e;
      e<<"ERROR: green2d from fgettr\n";
      e<<"  failed to read first trace from output file = "<<out<<"\n";
      throw e;
    }

    char * str = (char *)usermalloc_(128*sizeof(char));
    Value val;
    strcpy(str,"ns");
    gethdval(&tr,str,&val);
    int nt = vtoi(hdtype(str),val);

    strcpy(str,"dt");
    gethdval(&tr,str,&val);
    float dt = vtof(hdtype(str),val);

    strcpy(str,"gelev");
    gethdval(&tr,str,&val);
    float gz = vtof(hdtype(str),val);

    strcpy(str,"selev");
    gethdval(&tr,str,&val);
    float sz = vtof(hdtype(str),val);

    strcpy(str,"gx");
    gethdval(&tr,str,&val);
    float gx = vtof(hdtype(str),val);

    strcpy(str,"sx");
    gethdval(&tr,str,&val);
    float sx = vtof(hdtype(str),val);
    
    // compugte distance
    float r = sqrt((gx-sx)*(gx-sx) + (gz-sz)*(gz-sz));
    
    // convert dt to ms
    dt *= 0.001;
    
    float pi=4.0*atan(1.0);
    float fac=1.0/(2.0*bulk*buoy*pi);
    float c = sqrt(bulk*buoy);
    int ia = int(r/(c*dt))+1;
    
    //cerr<<"3\n";    
    // compute trace, copy to trace data buffer
    for (int i=0; i<nt; i++) {
      if (i<ia) tr.data[i]=0.0f;
      else {
	float t = i*dt;
	float s = iwave_max(dt*dt,t-r/c);
	tr.data[i] = fac/sqrt(s*(t+r/c));
      }
    }
    // scale by 1/bulk to remove assumed scaling in function
    //    for (int i=0; i<nt; i++) (tr.data)[i] = pbuf[i]/bulk;
    //cerr<<"4\n";
    // write trace to file
    fseeko(fpout,0L,SEEK_SET);
    fputtr(fpout,&tr);

    //cerr<<"5\n";
    
    // clean up
    fclose(fpout);
    userfree_(str);

    exit(0);
  }
  catch (RVL::RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
