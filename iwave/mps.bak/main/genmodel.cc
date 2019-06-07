/* generates simple 2D-3D models  
 * Code adapted from standardmodel.cc
 * Mario J. Bencomo
*/

#include <par.h>
#include <parser.h>
#include <math.h>

/* define number of model types */
#define NMODEL 5

#define C_MAX 4
#define C_MIN 2
#define RHO   1 //2.3

#define o1_def 0
#define o2_def 0
#define o3_def 0
#define d1_def 5
#define d2_def 5
#define d3_def 5
#define n1_def 101
#define n2_def 101
#define n3_def 1

#define relsigmax1_def .2
#define relsigmax2_def .2

#define cycles_def 2

#ifndef PI
#define PI 3.14159265359
#endif

/*********************** self documentation **********************/
const char *sdoc[] = {
    "                                                                ",
    " GENMODEL - generates simple 2-D and 3-D models.",
    "                                                                ",
    " This module builds some standard velocity and density cubes    ",
    " useful for comparing analytic solutions to finite difference   ",
    " variable-density acoustic simulators.                          ",
    "                                                                ",
    " Units:                                                         ",
    " densities: g/cm^3.                                             ",
    " velocities: m/ms = km/s.                                       ",
    " bulk moduli: GPa                                               ",
    " buoyancies: cm^3/g                                             ",
    "                                                                ",
    " The output format is native binary floats.                     ",
    "                                                                ",
    " Example run:   genmodel model=4 choose=2 > vp.bin         ",
    "          or:   genmodel model=4 choose=2 hfile=vp.rsf     ",
    "                                                                ",
    " Optional parameters:                                           ",
    "   model=1                    = choice of standard model        ",
    "                              =0; homogeneous                   ",
    "                              =1; linear velocity in depth, with",
    "                                  const density                 ",
    "                              =2; negative gaussian lense in    ",
    "                                  velocity, with const density  ",
    "                              =3; Sinusodial velocity in depth, ",
    "                                  with const density            ",
    "                              =4; half-plane model              ",
    "                                                                ",
    "   choose=1                   =0; output vp                     ",
    "                              =1; output bulk modulus           ",
    "                              =2; output density                ",
    "                              =3; output buoyancy               ",
    "                              =4; output vp*vp                  ",
    "                                                                ",
    "   hfile=\"\"                 =<hfile>: output in rsf/sep format",
    "                                  consisting of header file     ",
    "                                  hfile containing grid info,   ",
    "                                  and binary data file. If null ",
    "                                  string, output to stdout. If  ",
    "                                  non-null hfile given, path to ",
    "                                  data file is of form          ",
    "                                                                ",
    "                                  DATAPATH<hfile>@              ",
    "                                                                ",
    "                                  with DATAPATH read from env   ",
    "                                  (must end with / if set).     ",
    "                                                                ",
    "                                  Typical example: hfile=vp.rsf,",
    "                                  DATAPATH=/var/tmp/, so data   ",
    "                                  written to /var/tmp/vp.rsf@.  ",
    " ",
    "                                                                ",
    " Fast Dimension (depth)                                         ",
    "   o1=0.                      starting coordinate on grid       ",
    "   d1=5.0                     increment on grid in meters       ",
    "   e1=500.0                   end coordinate on grid            ",
    "   n1=(int)(1.5+(e1-o1)/d1)   number of gridpoints              ",
    "                                                                ",
    " Middle Dimension (x)                                           ",
    "   o2=0.                      starting coordinate on grid       ",
    "   d2=5.0                     increment on grid in meters       ",
    "   e2=500.0                   end coordinate on grid            ",
    "   n2=(int)(1.5+(e2-o2)/d2)   number of gridpoints              ",
    "                                                                ",
    " Slow Dimension (y)                                             ",
    "   o3=0.                      starting coordinate on grid       ",
    "   d3=5.0                     increment on grid in meters       ",
    "   e3=500.0                   end coordinate on grid            ",
    "   n3=(int)(1.5+(e3-o3)/d3)   number of gridpoints              ",
    "                                                                ",
    NULL};

//----------------------------------------------------------------------------//
// model 0: Homogeneous model
static inline float homogeneous(int choose, float * modelpars)
{
  float c_min = modelpars[0];
  float c_max = modelpars[1];
  float rho = modelpars[2];

  float vel = (c_max + c_min)*0.5;

  if (choose==1) return vel*vel*rho;
  if (choose==2) return rho;
  if (choose==3) return 1./rho;
  if (choose==4) return vel*vel;
  return vel;
}
//----------------------------------------------------------------------------//
// model 1: linear depth velocity model with constant density
static inline float lineardepth(float x1, int choose, float * modelpars)
{
  float c_min = modelpars[0];
  float c_max = modelpars[1];
  float rho = modelpars[2];
  float o1 = modelpars[3];
  float e1 = modelpars[4];
  
  float vel = (e1-x1)/(e1-o1)*c_min + (x1-o1)/(e1-o1)*c_max;
  
  if (choose==1) return vel*vel*rho;
  if (choose==2) return rho;
  if (choose==3) return 1./rho;
  if (choose==4) return vel*vel;
  
  return vel;
}
//----------------------------------------------------------------------------//
// model 2: negative gaussian lense in velocity, with constant density
static inline float gaussianlense2D(float x1, float x2, int choose, float * modelpars)
{
  float c_min = modelpars[0];
  float c_max = modelpars[1];
  float rho = modelpars[2];
  float o1 = modelpars[3];
  float e1 = modelpars[4];
  float o2 = modelpars[5];
  float e2 = modelpars[6];
  float cx1 = modelpars[7];
  float cx2 = modelpars[8];
  float relsigmax1 = modelpars[9];
  float relsigmax2 = modelpars[10];
  
  float sigmax1 = relsigmax1*(e1-o1);
  float sigmax2 = relsigmax2*(e2-o2);  
  float gaussx1 = exp( -(x1-cx1)*(x1-cx1)/(2.*sigmax1*sigmax1) );
  float gaussx2 = exp( -(x2-cx2)*(x2-cx2)/(2.*sigmax2*sigmax2) );
  float vel = c_max - (c_max-c_min)*gaussx1*gaussx2;
  
  if (choose==1) return vel*vel*rho;
  if (choose==2) return rho;
  if (choose==3) return 1./rho;
  if (choose==4) return vel*vel;
  return vel;
}
//----------------------------------------------------------------------------//
// model 3: Sinusodial-depth velocity with constant density
static inline float Sindepth(float x1, int choose, float * modelpars)
{

  float c_min = modelpars[0];
  float c_max = modelpars[1];
  float rho = modelpars[2];
  float o1 = modelpars[3];
  float e1 = modelpars[4];
  float cycles = modelpars[5];
  
  float wnum = cycles*2*PI;
  float vel = (c_max-c_min)*.5*sin(wnum*(x1-e1)/(e1-o1)) + (c_max+c_min)*.5;
  
  if (choose==1) return vel*vel*rho;
  if (choose==2) return rho;
  if (choose==3) return 1./rho;
  if (choose==4) return vel*vel;
  return vel;
}

//----------------------------------------------------------------------------//
// model 3: Sinusodial-depth velocity with constant density
static inline float halfspace(float x1, int choose, float * modelpars)
{

  float c_top = modelpars[0];
  float c_bot = modelpars[1];
  float rho_top = modelpars[2];
  float rho_bot = modelpars[3];
  float depth = modelpars[4];

  float vel;
  float rho;
  
  if( (depth-x1)/depth > 1.e-8 ){
    //  if(x1<=depth){
    vel = c_top;
    rho = rho_top;
  }
  else{
    vel = c_bot;
    rho = rho_bot;
  }
  
  if (choose==1) return vel*vel*rho;
  if (choose==2) return rho;
  if (choose==3) return 1./rho;
  if (choose==4) return vel*vel;
  return vel;
}

/******************************************************************************/
float get_zvalue(float x1, float x2, float x3,
                 int model, int choose, float * modelpars)
{
    float v;
	
    switch (model){
    case 0: v = homogeneous(choose,modelpars);           break;
    case 1: v = lineardepth(x1,choose,modelpars);        break;
    case 2: v = gaussianlense2D(x1,x2,choose,modelpars); break;
    case 3: v = Sindepth(x1,choose,modelpars);           break;
    case 4: v = halfspace(x1,choose,modelpars);          break;
    default: v = homogeneous(choose,modelpars);
    }
    return v;
}
/******************** end self doc **************/

int writemodel(
               int choose,
               int model,      /* choice of model                              */
               int n1,         /* number of samples in depth (fast dimension)  */
               int n2,         /* number of samples in x (middle dimension)    */
               int n3,         /* number of samples in y (slow dimenssion)     */
               float o1,       /* start depth                                  */
               float o2,       /* start x                                      */
               float o3,       /* start y                                      */
               float e1,       /* end depth                                    */
               float e2,       /* end x                                        */
               float e3,       /* end y                                        */
               float d1,       /* depth increment                              */
               float d2,       /* x increment                                  */
               float d3,       /* y increment                                  */
               FILE * fp,      /* output stream                                */
               float * omv     /* other model values                           */
	       )
{
  int j1, j2, j3;
  float x1, x2, x3;
  float *v   = NULL;
  float * modelpars = NULL;
  
  if (model==0) {
    modelpars = (float*)emalloc(3*sizeof(float));
    modelpars[0] = omv[0]; //c_min
    modelpars[1] = omv[1]; //c_max
    modelpars[2] = omv[2]; //rho
  }
  
  if (model==1) {
    modelpars = (float*)emalloc(5*sizeof(float));
    modelpars[0] = omv[0]; //c_min
    modelpars[1] = omv[1]; //c_max
    modelpars[2] = omv[2]; //rho
    modelpars[3] = o1;
    modelpars[4] = e1;
  }
  if (model==2) {
    modelpars = (float*)emalloc(11*sizeof(float));
    modelpars[0] = omv[0]; //c_min
    modelpars[1] = omv[1]; //c_max
    modelpars[2] = omv[2]; //rho
    modelpars[3] = o1;
    modelpars[4] = e1;
    modelpars[5] = o2;
    modelpars[6] = e2;
    modelpars[7] = omv[3]; //x1 center of gauss
    modelpars[8] = omv[4]; //x2 center of gauss
    modelpars[9] = omv[5]; //relative sigmax1
    modelpars[10] = omv[6]; //relative sigmax2
  }
  if (model==3) {
    modelpars = (float*)emalloc(6*sizeof(float));
    modelpars[0] = omv[0]; //c_min
    modelpars[1] = omv[1]; //c_max
    modelpars[2] = omv[2]; //rho
    modelpars[3] = o1;
    modelpars[4] = e1;
    modelpars[5] = omv[3];
  }
  if (model==4) {
    modelpars = (float*)emalloc(5*sizeof(float));
    modelpars[0] = omv[0]; //c_top
    modelpars[1] = omv[1]; //c_bot
    modelpars[2] = omv[2]; //rho_top
    modelpars[3] = omv[3]; //rho_bot
    modelpars[4] = omv[4]; //depth
  }

  v   = (float *)emalloc(n1*sizeof(float));
  
  for(j3 = 0; j3 < n3; j3++){
    x3 = o3 + j3 * d3;
    for(j2 = 0; j2 < n2; j2++){
      x2 = o2 + j2 * d2;
      for(j1 = 0; j1 < n1; j1++){
	x1 = o1 + j1 * d1;
	
	v[j1] = get_zvalue(x1, x2, x3, model, choose, modelpars);
      }
      if (fwrite(v, sizeof(float), n1, fp) != n1){
	fprintf(stderr, "write error\n");
	free(v);
	return 1;
      }
    }
  }
  free(v);
  free(modelpars);
  return 0;
}

/******************************************************************************/
int main(int argc, char **argv) {
    
  float o1,o2,o3,d1,d2,d3,e1,e2,e3;
  int n1,n2,n3,model,choose;
  float * omv; //other model values
  
  /* WWS */
  char * fname;
  char * dname;
  FILE * fp;
  /* end WWS */
  
  //int i,j;
  
  char * cwdpath;
  //  char * pathptr;
    
  /******************
   * get parameters
   ******************/
  
  PARARRAY * par = ps_new();
  
  if (ps_createargs(par,argc-1,argv+1)) {
    fprintf(stderr,
	    "ERROR. could not process command line\n");
    exit(1);
  }
  
  /***********************
   * end get parameters
   ***********************/
  
  xargc=argc; xargv=argv;
  requestdoc(0);
  
  if (!(ps_ffint(*par,"model",&model))) {
    if (model<0 || model >= NMODEL) {
      fprintf(stderr,"Error: genmodel.x\n");
      fprintf(stderr,"model index %d not defined: must lie in range [1,%d]\n",model,NMODEL);
      exit(1);
    }
  }
  else {
    fprintf(stdout,"Warning: genmodel.x\n");
    fprintf(stdout,"no model index given, so using model=1\n");
    model=1;
  }
  if (!(ps_ffint(*par,"choose",&choose))) {
    if (choose<0 || choose > 4) {
      fprintf(stderr,"Error: genmodel.x\n");
      fprintf(stderr,"choose index must be 0 (velocity), 1 (bulkmod), or 2 (density), or 3 (buoyancy), or 4 (vel2) \n");
            exit(1);
    }
  }
  else {
    fprintf(stdout,"Warning: genmodel.x\n");
    fprintf(stdout,"no choose index given, so using choose=1 (density)\n");
    choose=1;
  }
  
  if(ps_fffloat(*par,"o1",&o1)) {o1=o1_def; fprintf(stdout, "WARNING: o1 is set to the default value!\n"); }
  if(ps_fffloat(*par,"o2",&o2)) {o2=o2_def; fprintf(stdout, "WARNING: o2 is set to the default value!\n"); }
  if(ps_fffloat(*par,"o3",&o3)) {o3=o3_def; fprintf(stdout, "WARNING: o3 is set to the default value!\n"); }
  
  if(ps_fffloat(*par,"d1",&d1)) {d1=d1_def; fprintf(stdout, "WARNING: d1 is set to the default value!\n"); }
  if(ps_fffloat(*par,"d2",&d2)) {d2=d2_def; fprintf(stdout, "WARNING: d2 is set to the default value!\n"); }
  if(ps_fffloat(*par,"d3",&d3)) {d3=d3_def; fprintf(stdout, "WARNING: d3 is set to the default value!\n"); }
  
  
  if(ps_ffint(*par,"n1",&n1)) {n1 = n1_def; fprintf(stdout, "WARNING: n1 is set to the default value!\n"); }
  if(ps_ffint(*par,"n2",&n2)) {n2 = n2_def; fprintf(stdout, "WARNING: n2 is set to the default value!\n"); }
  if(ps_ffint(*par,"n3",&n3)) {n3 = n3_def; fprintf(stdout, "WARNING: n3 is set to the default value!\n"); }
  
  e1 = o1 + d1 * (n1 - 1);
  e2 = o2 + d2 * (n2 - 1);
  e3 = o3 + d3 * (n3 - 1);

  if (model==0||model==1){
    omv = (float*)emalloc(3*sizeof(float));
  }
  if (model==2) {
    omv = (float*)emalloc(7*sizeof(float));
    if(ps_fffloat(*par,"cx1",omv+3)) {
      omv[3]=(e1-o1)*.5; 
      fprintf(stdout, "WARNING: cx1 is set to the default value!\n"); 
    }
    if(ps_fffloat(*par,"cx2",omv+4)) {
      omv[4]=(e2-o2)*.5; 
      fprintf(stdout, "WARNING: cx2 is set to the default value!\n"); 
    }
    if(ps_fffloat(*par,"relsigmax1",omv+5)) {
      omv[5]=relsigmax1_def; 
      fprintf(stdout, "WARNING: relsigmax1 is set to the default value!\n"); 
    }
    if(ps_fffloat(*par,"relsigmax2",omv+6)) {
      omv[6]=relsigmax2_def; 
      fprintf(stdout, "WARNING: relsigmax2 is set to the default value!\n"); 
    }
  }
  if (model==3) {
    omv = (float*) emalloc(4*sizeof(float));
    if(ps_fffloat(*par,"cycles",omv+3)) {
      omv[3]=cycles_def; 
      fprintf(stdout, "WARNING: cycles is set to the default value!\n"); 
    }
  }
  if (model>=0&&model<=3){
    if(ps_fffloat(*par,"c_min",omv)) {
      omv[0]=C_MIN; 
      fprintf(stdout,"WARNING: c_min is set to the default value!\n"); 
    }
    if(ps_fffloat(*par,"c_max",omv+1)) {
      omv[1]=C_MAX; 
      fprintf(stdout,"WARNING: c_max is set to the default value!\n"); 
    }
    if(ps_fffloat(*par,"rho",omv+2)) {
      omv[2]=RHO; 
      fprintf(stdout,"WARNING: rho is set to the default value!\n"); 
    }
  }

  if (model==4){
    omv = (float*) emalloc(5*sizeof(float));
    if(ps_fffloat(*par,"c_top",omv)) {
      omv[0]=C_MIN; 
      fprintf(stdout, "WARNING: c_top is set to the default value!\n"); 
    }
    if(ps_fffloat(*par,"c_bot",omv+1)) {
      omv[1]=C_MAX; 
      fprintf(stdout, "WARNING: c_bot is set to the default value!\n"); 
    }
    if(ps_fffloat(*par,"rho_top",omv+2)) {
      omv[2]=RHO; 
      fprintf(stdout, "WARNING: rho_top is set to the default value!\n"); 
    }
    if(ps_fffloat(*par,"rho_bot",omv+3)) {
      omv[3]=RHO; 
      fprintf(stdout, "WARNING: rho_bot is set to the default value!\n"); 
    }
    if(ps_fffloat(*par,"depth",omv+4)) {
      omv[4]=(e1+o1)*0.5; 
      fprintf(stdout, "WARNING: depth is set to the default value!\n"); 
    }
  }
  
  fprintf(stdout," o1=%f e1=%f d1=%f n1=%d\n",o1,e1,d1,n1);
  fprintf(stdout," o2=%f e2=%f d2=%f n2=%d\n",o2,e2,d2,n2);
  fprintf(stdout," o3=%f e3=%f d3=%f n3=%d\n",o3,e3,d3,n3);
  //fprintf(stdout," c_min=%f c_max=%f rho=%f\n",omv[0],omv[1],omv[2]);
  
  /* WWS */
  if (!(ps_ffcstring(*par,"hfile",&fname))) {
    
    /* DATAPATH deprecated - WWS 08.01.12*/
    /* DATAPATH revived - WWS 28.05.12 */
    /*    if (getenv("DATAPATH")) {
    //      dname=malloc(strlen(getenv("DATAPATH"))+strlen(fname)+2);
    //      strcpy(dname,getenv("DATAPATH"));
    // DATAPATH made preferentially a parameter WWS 08.28.12
    // else env variable else */
    
    cwdpath = NULL;
    /* if you can't get it from the parfile */
    if (ps_flcstring(*par,"datapath",&cwdpath)) {
      /* try to get it from the environment */
      cwdpath = (char *)malloc(128*sizeof(char));
      memset(cwdpath,'\0',128);
      
      //pathptr = getenv("DATAPATH");
      //if (pathptr) strcpy(cwdpath,pathptr);
      /* otherwise set to cwd */
      //else strcpy(cwdpath,".");
      strcpy(cwdpath,".");
    }

    dname=(char *)malloc(sizeof(char)*(strlen(cwdpath)+strlen(fname)+2));
    strcpy(dname,cwdpath);
    
    if (cwdpath[strlen(cwdpath)-1] != '/') strcat(dname,"/");
    strcat(dname,fname);
    strcat(dname,"@");
    
    fprintf(stdout,"writing header file %s\n",fname);
    if (!(fp=fopen(fname,"w"))) {
      fprintf(stderr,"Error: genmodel\n");
      fprintf(stderr,"failed to open new header file %s\n",fname);
      exit(1);
    }
    fprintf(fp,"n1=%d d1=%e o1=%e\n",
	    n1,d1,o1);
    fprintf(fp,"n2=%d d2=%e o2=%e\n",
	    n2,d2,o2);
    fprintf(fp,"n3=%d d3=%e o3=%e\n",
	    n3,d3,o3);
    switch (choose){
    case 0: fprintf(fp,"data_type=velocity\n"); break;
    case 1: fprintf(fp,"data_type=bulkmod\n"); break;
    case 2: fprintf(fp,"data_type=density\n"); break;
    case 3: fprintf(fp,"data_type=buoyancy\n"); break;
    case 4: fprintf(fp,"data_type=velocity squared\n"); break;
    default: fprintf(stderr,"choose failed\n"); exit(1);
    }
    
    fprintf(fp,"data_format=native_float\n");
    fprintf(fp,"in=%s\n",dname);
    fclose(fp);
    
    fprintf(stdout,"writing data to %s\n",dname);
    if (!(fp=fopen(dname,"w"))) {
      fprintf(stderr,"Error: genmodel\n");
      fprintf(stderr,"failed to open new data file %s\n",dname);
      exit(1);
    }
  }
  
  else fp=stdout;
  
  ps_delete(&par);
  
  writemodel(choose,model,n1,n2,n3,o1,o2,o3,e1,e2,e3,d1,d2,d3,fp,omv);
  
  free(omv);
  exit(0);
}
