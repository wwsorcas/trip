//#include "cstd.h"

#define ABSORB
//#define BARF

/*void asg_pstep2d(float ** restrict bulk, 
		 float ** restrict p0, float ** restrict p1,
		 float ** restrict v0, float ** restrict v1,
		 float ** restrict ep, float ** restrict epp,
		 float * restrict sdiv,
		 int * gsc, int * gec, 
		 int * lbc, int * rbc,
		 int maxoff, float ** restrict c) {
*/
/* signature for tapenade*/
void asg_pstep2d(float ** bulk, 
		 float ** p0, float ** p1,
		 float ** v0, float ** v1,
		 float ** ep, float ** epp,
		 float * sdiv,
		 int * gsc, int * gec, 
		 int * lbc, int * rbc,
		 int maxoff, float ** c) {
  

  int i0, i1;
  int ioff;


#ifdef BARF
  fprintf(stderr,"gsc[0]=%d gec[0]=%d\n",gsc[0],gec[0]);
  for (i0=gsc[0]; i0<=gec[0]; i0++) 
    fprintf(stderr,"i0=%d ep0=%e epp0=%e\n",i0,ep[0][i0],epp[0][i0]);
  fprintf(stderr,"gsc[1]=%d gec[1]=%d\n",gsc[1],gec[1]);
  for (i1=gsc[1]; i1<=gec[1]; i1++)
    fprintf(stderr,"i1=%d ep1=%e epp1=%e\n",i1,ep[1][i1],epp[1][i1]);
#endif

  for (i1=gsc[1]; i1 <= gec[1]; i1++) {
    for (i0=gsc[0]; i0 <= gec[0]; i0++) {
      sdiv[i0] = 
	(c[0][0]*(v0[i1][i0]-v0[i1][i0-1]) +
	 c[1][0]*(v1[i1][i0]-v1[i1-1][i0]));
    }
    for (ioff = 1; ioff<maxoff; ioff++) {
      for (i0=gsc[0]; i0 <= gec[0]; i0++) {
	sdiv[i0] += 
	  (c[0][ioff]*(v0[i1][i0+ioff]-v0[i1][i0-ioff-1]) +
	   c[1][ioff]*(v1[i1+ioff][i0]-v1[i1-ioff-1][i0]));
      }
    }
    for (i0=gsc[0]; i0 <= gec[0]; i0++) {
#ifdef ABSORB
      //      fprintf(stderr,"i1=%d i0=%d bulk=%e\n",i1,i0,bulk[i1][i0]);
      p0[i1][i0] = (p0[i1][i0]*ep[0][i0] - bulk[i1][i0]*sdiv[i0])*epp[0][i0];
      p1[i1][i0] = (p1[i1][i0]*ep[1][i1] - bulk[i1][i0]*sdiv[i0])*epp[1][i1];
#else
      p0[i1][i0] = p0[i1][i0] - bulk[i1][i0]*sdiv[i0];
      p1[i1][i0] = p0[i1][i0];
#endif
    }
  }

  /* boundary conditions - p is odd about index just before/after comp domain */
  if (lbc[0]) {
    for (i1=gsc[1];i1<=gec[1];i1++) {
      p0[i1][gsc[0]-1]=0.0f;
      for (ioff=1;ioff<maxoff;ioff++) {
	p0[i1][gsc[0]-ioff-1]=-p0[i1][gsc[0]+ioff-1];
      }
    }
  }
  if (rbc[0]) {
    for (i1=gsc[1];i1<=gec[1];i1++) {
      p0[i1][gec[0]+1]=0.0f;
      for (ioff=1;ioff<maxoff;ioff++) {
	p0[i1][gec[0]+ioff+1]=-p0[i1][gec[0]-ioff+1];
      }
    }
  }
  if (lbc[1]) {
    for (i0=gsc[0];i0<=gec[0];i0++) {
      p1[gsc[1]-1][i0]=0.0f;
      for (ioff=1;ioff<maxoff;ioff++) {
	p1[gsc[1]-ioff-1][i0]=-p1[gsc[1]+ioff-1][i0];
      }
    }
  }
  if (rbc[1]) {
    for (i0=gsc[0];i0<=gec[0];i0++) {
      p1[gec[1]+1][i0]=0.0f;
      for (ioff=1;ioff<maxoff;ioff++) {
	p1[gec[1]+ioff+1][i0]=-p1[gec[1]-ioff+1][i0];
      }
    }
  }
}
