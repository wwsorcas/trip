#include <su.h>
#include <segy.h>
#include <header.h>
#include <cubic.h>
#include <parser.h>

const char * sdoc[] = { 
  "Usage: rms.x in= out=",
  "",
  "Purpose: output trace rms in segy file with traces of length 1",
  "",
  "NOTE: if it exists on call output file is overwritten, else it is created.",
  "",
  "Required parameters:",
  "  in [string]  = filename for input",
  "  out [string]  = filename for output (traces with nt=1)",
  "",
  NULL};


int main(int argc, char ** argv) {
  
  /* ARGUMENTS */
  char * in;      /* input file name */
  char * out;      /* output file name */


  /* INTERNAL VARIABLES */
  PARARRAY * par;  /* param array */
  FILE * fpin;     /* input file pointer */
  FILE * fpout;    /* output file pointer */
  segy tr;         /* trace workspace */
  segy tr1;         /* trace workspace */
  Value val;       /* header word workspace */
  float dt=0.0f;   /* time step (ms)*/
  int nt=0;        /* number of samples per trace */
  int i;           /* sample counter */

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

  if (ps_flcstring(*par,"out",&out)) {
    printf("error reading output file name. ABORT.\n");
    exit(1);
  }

  /* open data files */
  if (!(fpin=fopen(in,"r"))) {
    printf("failed to open input file = %s for read. ABORT.\n",in);
    exit(1);
  }

  if (!(fpout=fopen(out,"w"))) {
    printf("failed to open output file = %s for write. ABORT.\n",out);
    exit(1);
  }

  /* read loop */

  char * str = (char *)malloc(128*sizeof(char));

  while (fgettr(fpin,&tr)) {

    /* read time geometry, reset delrt */
    strcpy(str,"ns");
    gethdval(&tr,str,&val);
    nt = vtoi(hdtype(str),val);
    val.u = 1;
    puthdval(&tr1,str,&val);

    strcpy(str,"dt");
    gethdval(&tr,str,&val);
    puthdval(&tr1,str,&val);
    dt = 0.001*vtof(hdtype(str),val);
      
    /* trace rms */
    float tmp = 0;
    for (i=0;i<nt;i++) {
      tmp += tr.data[i]*tr.data[i];
    }
    tr1.data[0] = sqrt(dt*tmp);

    /* write to output */
    fputtr(fpout,&tr1);
  }

  free(str);
  fclose(fpin);
  fclose(fpout);
  ps_delete(&par);

  exit(0);
}
  
