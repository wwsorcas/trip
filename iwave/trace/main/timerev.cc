#include <su.h>
#include <segy.h>
#include <header.h>
#include <cubic.h>
#include <parser.h>

#define DT_TOL 0.001

const char * sdoc[] = { 
  "Usage: timerev.x in= out=",
  "",
  "Purpose: time-reverse the traces in the input file, write to specified",
  "output file. Unlike suflip, correctly sets delrt, therefore geometry of",
  "time axis.",
  "",
  "NOTE: if it exists on call output file is overwritten, else it is created."
  "",
  "Required parameters:",
  "  in [string]  = filename for input",
  "  out [string]  = filename for output (time reversal of input)",
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
  Value val;       /* header word workspace */
  float dt=0.0f;   /* time step (ms)*/
  int nt=0;        /* number of samples per trace */
  float delrt;     /* time of first sample in OUTPUT */
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

  printf("******************************************\n");
  printf("*    Time reversal of SEGY data sets      \n");
  printf("*                                         \n");
  printf("* Input  = %s\n",in);
  printf("* Output = %s\n",out);
  printf("******************************************\n");
  printf("\n");

  /* read loop */

  while (fgettr(fpin,&tr)) {

    char * str = (char *)malloc(128*sizeof(char));

    /* read time geometry, reset delrt */
    strcpy(str,"ns");
    gethdval(&tr,str,&val);
    nt = vtoi(hdtype(str),val);

    strcpy(str,"dt");
    gethdval(&tr,str,&val);
    dt = 0.001*vtof(hdtype(str),val);

    strcpy(str,"delrt");
    gethdval(&tr,str,&val);
    delrt = vtof(hdtype(str),val);
    
    /* reset delrt */

    delrt = -(delrt + (nt-1)*dt);
    val.h=(int)(delrt);
    puthdval(&tr,str,&val);

    free(str);
    
    /* reverse samples */
    float tmp;
    for (i=0;i<(nt+1)/2;i++) {
      tmp = tr.data[i];
      tr.data[i] = tr.data[nt-i];
      tr.data[nt-i] = tmp;
    }

    /* write to output */
    fputtr(fpout,&tr);
  }

  fclose(fpin);
  fclose(fpout);
  ps_delete(&par);

  exit(0);
}
  
