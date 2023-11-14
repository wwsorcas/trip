#include <su.h>
#include <segy.h>
#include <header.h>
#include <cubic.h>
#include <parser.h>
#include <utility.hh>

const char * sdoc[] = { 
  "Usage: gain.x in= scale= out= spower= tpower=",
  "",
  "Purpose: apply scaling of input traces by trace-dep factor read from ",
  "segy file with trace length 1, raised to specified power, and also by a",
  "specified power of t",
  "",
  "NOTE: if it exists on call output file is overwritten, else it is created.",
  "",
  "Required parameters:",
  "  in [string]  = filename for input",
  "  sfile [string] = filename for rms values (data[0] in each trace)",
  "  out [string] = filename for output",
  "  spower [int] = 0 power of scale factor",
  "  tpower [int] = 1 power of t (must be non-negative)",
  "",
  "output samples = input samples scaled by the tpower-th power of t, and",
  "optionally by a power of a factor read from sfile.",
  "",
  "If spower=0, then sfile is an inactive argument (and can be left out), ",
  "and further scaling is not applied. If spower != 0, then the file ",
  "named by sfile must exist and have at least the number of traces in the ",
  "input file. The first data sample in each trace should be positive, and ",
  "output samples in corresponding trace are also scaled by the",
  "spower-th power of that sfile sample. So",
  "",
  "output(t) = input(t) * (t**tpower) * (sfile.data[0]**spower)",
  "",
  "NOTE: this function is required for an implementation of AWI.",
  NULL};


int main(int argc, char ** argv) {
  
  /* ARGUMENTS */
  char * in;       /* input file name */
  char * sfile;    /* scale data file name */
  char * out;      /* output file name */
  int spower;      /* power of scale factor */
  int tpower;      /* power of t */

  /* INTERNAL VARIABLES */
  PARARRAY * par;  /* param array */
  FILE * fpin;     /* input file pointer */
  FILE * fpsc;     /* scale file */
  FILE * fpout;    /* output file pointer */
  segy tr;         /* trace workspace */
  segy tr1;        /* trace workspace */
  Value val;       /* header word workspace */
  float dt;        /* time step (ms)*/
  float delrt;     /* time of first sample */
  int nt=0;        /* number of samples per trace */
  int i;           /* sample counter */
  int idelrt;      /* index of first sample on real line */
  int j;           /* trace counter */
  float fac;       /* scale factor */
  float * tp;      /* time scale trace */

  xargc=argc; xargv=argv;
  requestdoc(1);

  /* extract input parameters */
  par=ps_new();
  if ( ps_createargs(par, argc - 1, argv + 1) ) {
    fprintf(stderr,"Error parsing input data. ABORT.\n");
    exit(1);
  }
             
  if (ps_flcstring(*par,"in",&in)) {
    fprintf(stderr,"error reading input file name ABORT.\n");
    exit(1);
  }

  // default: spower=0 (no scale factor)
  if (ps_flint(*par,"spower",&spower)) {
    spower = 0;
  }

  // default: tpower=1 (scale by t)
  if (ps_flint(*par,"tpower",&tpower)) {
    tpower = 1;
  }
  if (tpower < 0) {
    fprintf(stderr,"error tpower < 0 ABORT.\n");
    exit(1);
  }
  
  if (spower != 0) {
    if (ps_flcstring(*par,"sfile",&sfile)) {
      fprintf(stderr,"error reading scale file name ABORT.\n");
      exit(1);
    }
  }
  
  if (ps_flcstring(*par,"out",&out)) {
    fprintf(stderr,"error reading output file name. ABORT.\n");
    exit(1);
  }

  /* open data files */
  if (!(fpin=fopen(in,"r"))) {
    fprintf(stderr,"failed to open input file = %s for read. ABORT.\n",in);
    exit(1);
  }

  if (spower) {
    if (!(fpsc=fopen(sfile,"r"))) {
      fprintf(stderr,"failed to open input file = %s for read. ABORT.\n",sfile);
      exit(1);
    }
  }
  
  if (!(fpout=fopen(out,"w"))) {
    fprintf(stderr,"failed to open output file = %s for write. ABORT.\n",out);
    exit(1);
  }

  //fprintf(stderr,"read loop\n");

  // initialize workspace
  tp = NULL;

  j=0;

  char * str = (char *)usermalloc_(128*sizeof(char));

  while (fgettr(fpin,&tr)) {

    if (spower) {
      if (!(fgettr(fpsc,&tr1))) {
	fprintf(stderr,"failed to read scale trace %d \n",j);
	exit(1);
      }
    }
    
    /* read time geometry */
    strcpy(str,"ns");
    gethdval(&tr,str,&val);
    nt = vtoi(hdtype(str),val);

    strcpy(str,"dt");
    gethdval(&tr,str,&val);
    dt = 0.001*vtof(hdtype(str),val);

    strcpy(str,"delrt");
    gethdval(&tr,str,&val);
    delrt = vtof(hdtype(str),val);
    
    // initialize time scale trace on first read
    if ((tp == NULL) && (tpower != 0)) {
      tp = (float *)usermalloc_(nt*sizeof(float));
      idelrt = (int)((delrt/dt) + 0.1f);
      for (i=0;i<nt;i++) {
	tp[i] = powf((float)(i+idelrt)*dt,tpower);
      }
    }
      
    // scale factor = (first sample of scale trace ** spower) * (dt ** tpower)
    if (spower > 0) {
      fac = powf(tr1.data[0],spower);
    }
    else if (spower < 0) {
      if (RVL::ProtectedDivision<float>(1.0,
					powf(tr1.data[0], abs(spower)),
					fac)) {
	fac=0.0f;
      }
    }
    else fac = 1.0;

    if (tpower == 0) {
      if (spower != 0) {
	for (i=0;i<nt;i++) {
	  tr.data[i] *= fac;
	}
      }
    }
    else {
      for (i=0;i<nt;i++) {
	tr.data[i] *= tp[i]*fac;
      }
    }
    /* write to output */
    fputtr(fpout,&tr);

    j++;
  }

  userfree_(str);
  if (tpower) userfree_(tp);
  fclose(fpin);
  if (spower) fclose(fpsc);
  fclose(fpout);
  ps_delete(&par);

  exit(0);
}
  
