#define SUPERVERBOSE

/*        Generated by TAPENADE     (INRIA, Ecuador team)
	  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
*/

/* hand clean-up 2016 Oct 21 WWS */
//#include "GlobalDeclarations_b.c"

/*
  Differentiation of asg_pstep2d in reverse (adjoint) mode:
  gradient     of useful results: **p0 **p1
  with respect to varying inputs: **v0 **v1 **bulk **p0 **p1
  RW status of diff variables: **v0:out **v1:out **bulk:out **p0:in-out
  **p1:in-out
  Plus diff mem management of: v0:in *v0:in v1:in *v1:in bulk:in
  *bulk:in p0:in *p0:in p1:in *p1:in sdiv:in

  void asg_pstep2d(float ** restrict bulk, 
  float ** restrict p0, float ** restrict p1,
  float ** restrict v0, float ** restrict v1,
  float ** restrict ep, float ** restrict epp,
  float * restrict sdiv,
  int * gsc, int * gec, 
  int * lbc, int * rbc,
  int maxoff, float ** restrict c) {

  signature for tapenade*/
#include "cstd.h"

void asg_pstep2d_b(float ** restrict bulk,
		   float ** restrict bulkb,
		   float ** restrict p0,
		   float ** restrict p0b,
		   float ** restrict p1,
		   float ** restrict p1b,
		   float ** restrict v0,
		   float ** restrict v0b,
		   float ** restrict v1,
		   float ** restrict v1b, 
		   float ** restrict ep,
		   float ** restrict epp,
		   float *sdiv, float *sdivb,
		   int *gsc, int *gec,
		   int *lbc, int *rbc,
		   int maxoff, float ** restrict c) {
  int i0, i1;
  int ioff;
  float tempb;
  float tempb0;
  float tempb1;
  float tempb2;
  float tmpb;
  float tmpb0;
  float tmpb1;
  float tmpb2;

#ifdef SUPERVERBOSE
  fprintf(stderr,"asg_pstep2d_b: l2 norms squared\n");
  float jnk;
  jnk=0.0;
  for (i1 = gsc[1]; i1 < gec[1]+1; ++i1) {
    for (i0 = gsc[0]; i0 < gec[0]+1; ++i0) {
      jnk += bulkb[i1][i0]*bulkb[i1][i0];
    }
  }
  fprintf(stderr,"  bulkb=%e",jnk);
  jnk=0.0;
  for (i1 = gsc[1]; i1 < gec[1]+1; ++i1) {
    for (i0 = gsc[0]; i0 < gec[0]+1; ++i0) {
      jnk += p0b[i1][i0]*p0b[i1][i0];
    }
  }
  fprintf(stderr,"  p0b=%e",jnk);
  fprintf(stderr,"\n");
#endif
  
  for (i1 = gsc[1]; i1 < gec[1]+1; ++i1) {

  }
  /* boundary conditions - p is odd about index just before/after comp domain 
   */
  /*
    if (lbc[0])
    pushcontrol1b(0);
    else
    pushcontrol1b(1);
    if (rbc[0])
    pushcontrol1b(0);
    else
    pushcontrol1b(1);
    if (lbc[1])
    pushcontrol1b(0);
    else
    pushcontrol1b(1);
  */
  if (rbc[1])
    for (i0 = gec[0]; i0 > gsc[0]-1; --i0) {
      for (ioff = maxoff-1; ioff > 0; --ioff) {
	tmpb2 = p1b[gec[1] + ioff + 1][i0];
	p1b[gec[1] + ioff + 1][i0] = 0.0;
	p1b[gec[1] - ioff + 1][i0] = p1b[gec[1] - ioff + 1][i0] - 
	  tmpb2;
      }
      p1b[gec[1] + 1][i0] = 0.0;
    }
  //    popcontrol1b(&branch);
  //    if (branch == 0)
  if (lbc[1]) 
    for (i0 = gec[0]; i0 > gsc[0]-1; --i0) {
      for (ioff = maxoff-1; ioff > 0; --ioff) {
	tmpb1 = p1b[gsc[1] - ioff - 1][i0];
	p1b[gsc[1] - ioff - 1][i0] = 0.0;
	p1b[gsc[1] + ioff - 1][i0] = p1b[gsc[1] + ioff - 1][i0] - 
	  tmpb1;
      }
      p1b[gsc[1] - 1][i0] = 0.0;
    }
  //    popcontrol1b(&branch);
  //    if (branch == 0)
  if (rbc[0])
    for (i1 = gec[1]; i1 > gsc[1]-1; --i1) {
      for (ioff = maxoff-1; ioff > 0; --ioff) {
	tmpb0 = p0b[i1][gec[0] + ioff + 1];
	p0b[i1][gec[0] + ioff + 1] = 0.0;
	p0b[i1][gec[0] - ioff + 1] = p0b[i1][gec[0] - ioff + 1] - 
	  tmpb0;
      }
      p0b[i1][gec[0] + 1] = 0.0;
    }
  //    popcontrol1b(&branch);
  //    if (branch == 0)
  if (lbc[0])
    for (i1 = gec[1]; i1 > gsc[1]-1; --i1) {
      for (ioff = maxoff-1; ioff > 0; --ioff) {
	tmpb = p0b[i1][gsc[0] - ioff - 1];
	p0b[i1][gsc[0] - ioff - 1] = 0.0;
	p0b[i1][gsc[0] + ioff - 1] = p0b[i1][gsc[0] + ioff - 1] - tmpb
	  ;
      }
      p0b[i1][gsc[0] - 1] = 0.0;
    }
  /*
  **v0b = 0.0;
  **v1b = 0.0;
  **bulkb = 0.0;
  *sdivb = 0.0;
  */
  for (i1 = gec[1]; i1 > gsc[1]-1; --i1) {

    for (i0 = gsc[0]; i0 < gec[0]+1; ++i0) {
      sdiv[i0] =
	c[0][0]*(v0[i1][i0]-v0[i1][i0-1]) +
	c[1][0]*(v1[i1][i0]-v1[i1-1][i0]);
    }
    for (ioff = 1; ioff < maxoff; ++ioff) {
      for (i0 = gsc[0]; i0 < gec[0]+1; ++i0) {
	sdiv[i0] = sdiv[i0] +
	  (c[0][ioff]*(v0[i1][i0+ioff]-v0[i1][i0-ioff-1]) +
	   c[1][ioff]*(v1[i1+ioff][i0]-v1[i1-ioff-1][i0]));
      }
    }
    for (i0 = gec[0]; i0 > gsc[0]-1; --i0) {
      tempb1 = epp[1][i1]*p1b[i1][i0];
      bulkb[i1][i0] = bulkb[i1][i0] - sdiv[i0]*tempb1;
      sdivb[i0] = sdivb[i0] - bulk[i1][i0]*tempb1;
      p1b[i1][i0] = ep[1][i1]*tempb1;
      tempb2 = epp[0][i0]*p0b[i1][i0];
      bulkb[i1][i0] = bulkb[i1][i0] - sdiv[i0]*tempb2;
      sdivb[i0] = sdivb[i0] - bulk[i1][i0]*tempb2;
      p0b[i1][i0] = ep[0][i0]*tempb2;
    }
    for (ioff = maxoff-1; ioff > 0; --ioff) {
      for (i0 = gec[0]; i0 > gsc[0]-1; --i0) {
	tempb = c[0][ioff]*sdivb[i0];
	tempb0 = c[1][ioff]*sdivb[i0];
	v0b[i1][i0 + ioff] = v0b[i1][i0 + ioff] + tempb;
	v0b[i1][i0 - ioff - 1] = v0b[i1][i0 - ioff - 1] - tempb;
	v1b[i1 + ioff][i0] = v1b[i1 + ioff][i0] + tempb0;
	v1b[i1 - ioff - 1][i0] = v1b[i1 - ioff - 1][i0] - tempb0;
      }
    }
    for (i0 = gec[0]; i0 > gsc[0]-1; --i0) {
      v0b[i1][i0] = v0b[i1][i0] + c[0][0]*sdivb[i0];
      v0b[i1][i0 - 1] = v0b[i1][i0 - 1] - c[0][0]*sdivb[i0];
      v1b[i1][i0] = v1b[i1][i0] + c[1][0]*sdivb[i0];
      v1b[i1 - 1][i0] = v1b[i1 - 1][i0] - c[1][0]*sdivb[i0];
      sdivb[i0] = 0.0;
    }
  }
}
