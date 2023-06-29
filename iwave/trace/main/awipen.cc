#include <su.h>
#include <segy.h>
#include <header.h>
#include <cubic.h>
#include <parser.h>
#include <utility.hh>

const char * sdoc[] = { 
  "Usage: rms.x in= out=",
  "",
  "Purpose: apply awi penalty op, using trace rms in segy file with ",
  "traces of length 1 to scale the mult-by-t operator",
  "",
  "NOTE: if it exists on call output file is overwritten, else it is created.",
  "",
  "Required parameters:",
  "  in [string]  = filename for input",
  "  rms [string] = filename for rms values (data[0] in each trace)",
  "  out [string] = filename for output (traces with nt=1)",
  "  alpha [float] = penalty parameter",
  "  precond [int] = preconditioning flag",
  "",
  "If precond=0, then rms is an inactive argument (and can be left out), ",
  "and awi preconditioning is not applied. If precond != 0, then the file ",
  "named by rms must exist and have at least the number of traces in the ",
  "input file. The first data word in each trace should be positive, and ",
  "is used as the rms reciprocal scaling in awi.",
  NULL};


int main(int argc, char ** argv) {
  
  /* ARGUMENTS */
  char * in;       /* input file name */
  char * rms;      /* rms data file name */
  char * out;      /* output file name */
  float alpha;     /* penalty param */
  int precond;     /* precondition (1) or not (0) */

  /* INTERNAL VARIABLES */
  PARARRAY * par;  /* param array */
  FILE * fpin;     /* input file pointer */
  FILE * fprms;    /* rms file */
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

  xargc=argc; xargv=argv;
  requestdoc(1);

  /* extract input parameters */
  par=ps_new();
  if ( ps_createargs(par, argc - 1, argv + 1) ) {
    printf("Error parsing input data. ABORT.\n");
    exit(1);
  }
             
  if (ps_flcstring(*par,"in",&in)) {
    printf("error reading input file name ABORT.\n");
    exit(1);
  }

    // default: use awi preconditioning
  if (ps_flint(*par,"precond",&precond)) {
    precond = 1;
  }

  if (precond) {
    if (ps_flcstring(*par,"rms",&rms)) {
      printf("error reading input file name ABORT.\n");
      exit(1);
    }
  }
  
  if (ps_flcstring(*par,"out",&out)) {
    printf("error reading output file name. ABORT.\n");
    exit(1);
  }

  if (ps_flfloat(*par,"alpha",&alpha)) {
    printf("error reading output file name. ABORT.\n");
    exit(1);
  }


  /* open data files */
  if (!(fpin=fopen(in,"r"))) {
    printf("failed to open input file = %s for read. ABORT.\n",in);
    exit(1);
  }

  if (precond) {
    if (!(fprms=fopen(rms,"r"))) {
      printf("failed to open input file = %s for read. ABORT.\n",rms);
      exit(1);
    }
  }
  
  if (!(fpout=fopen(out,"w"))) {
    printf("failed to open output file = %s for write. ABORT.\n",out);
    exit(1);
  }

  /* read loop */

  j=0;

  char * str = (char *)malloc(128*sizeof(char));

  while (fgettr(fpin,&tr)) {

    if (precond) {
      if (!(fgettr(fprms,&tr1))) {
	printf("failed to read rms trace %d \n",j);
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
      
    // if rms underflows, set scale = 0
    if (precond) {
      if (RVL::ProtectedDivision<float>(alpha*dt,tr1.data[0],fac)) fac=0.0f;
    }
    else fac=alpha*dt;

    idelrt = (int)((delrt/dt) + 0.1f);

    //cout<<"nt = "<<nt<<" dt = "<<dt<<" fac = "<<fac<<endl;
    /* trace rms */
    for (i=0;i<nt;i++) {
      tr.data[i] *= (i+idelrt)*fac;
    }

    /* write to output */
    fputtr(fpout,&tr);

    j++;
  }

  free(str);
  fclose(fpin);
  fclose(fprms);
  fclose(fpout);
  ps_delete(&par);

  exit(0);
}
  
