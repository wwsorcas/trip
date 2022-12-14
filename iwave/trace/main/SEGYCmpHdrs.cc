#include <su.h>
#include <segy.h>
#include <header.h>
#include <cubic.h>
#include <parser.h>

const char * sdoc[] = { 
  "Usage: SEGYCmpHdrs.x in1= in2=",
  "",
  "Compares two SU files for compatibility (defining the same SEGYSpace).",
  "First checks that files have same lengths, then that trace headers are",
  "the same byte-wise."
  " ",
  "Returns 0 on success.",
  "",
  NULL};

int main(int argc, char ** argv) {
  
  /* ARGUMENTS */
  char * in1;      /* input 1 file name */
  char * in2;      /* input 2 file name */

  /* INTERNAL VARIABLES */
  PARARRAY * par;  /* param array */
  FILE * fp1;      /* input 1 file pointer */
  FILE * fp2;      /* input 2 file pointer */
  segy tr1;        /* input 1 trace workspace */
  segy tr2;        /* input 2 trace workspace */
  int nb1;         /* number of bytes read on trace 1 */
  int nb2;         /* number of bytes read on trace 2 */
  int itr;         /* trace counter */

  xargc=argc; xargv=argv;
  requestdoc(1);

  /* extract input parameters */
  par=ps_new();
  if ( ps_createargs(par, argc - 1, argv + 1) ) {
    printf("Error parsing input data. ABORT.\n");
    exit(1);
  }
             
  if (ps_flcstring(*par,"in1",&in1)) {
    printf("COMP reading in1. ABORT.\n");
    exit(1);
  }

  if (ps_flcstring(*par,"in2",&in2)) {
    printf("COMP reading in2. ABORT.\n");
    exit(1);
  }

  /* open data files */
  if (!(fp1=fopen(in1,"r"))) {
    printf("COMP: failed to open 1st input file = %s. ABORT.\n",in1);
    exit(1);
  }

  if (!(fp2=fopen(in2,"r"))) {
    printf("COMP: failed to open 2nd input file = %s. ABORT.\n",in2);
    exit(1);
  }

  /* first check that file lengths are same */
  if ((fseeko(fp1,0L,SEEK_END)) || (fseeko(fp2,0L,SEEK_END))) {
    printf("COMP: seek-to-end failed - files comprmised\n");
    exit(1);
  }
  if (ftello(fp1) != ftello(fp2)) {
    printf("COMP: files have different lengths\n");
    exit(1);
  }
  if ((fseeko(fp1,0L,SEEK_SET)) || (fseeko(fp2,0L,SEEK_SET))) {
    printf("COMP: seek-to-start failed - files comprmised\n");
    exit(1);
  }  

  itr=0;
  /* read loop */
  while ((nb1=fgettr(fp1,&tr1))*(nb2=fgettr(fp2,&tr2))) {
    // traces read - compare headers byte-wise using
    // code borrowed from sucmp.c
    if (memcmp( &tr1, &tr2, HDRBYTES)) {
      printf("COMP: Files %s & %s differ in headers at trace %d\n",
	     in1, in2, itr);
      exit(1);
    }
    itr++;
  } 
  if (((nb1 == 0) && (nb2 != 0)) ||
      ((nb1 != 0) && (nb2 == 0))) {
    // not same number of traces - in principle can't happen
    printf("COMP: number of traces differs\n");
    exit(1);
  }
  fclose(fp1);
  fclose(fp2);
  iwave_fdestroy();
  ps_delete(&par);

  exit(0);
}
  
