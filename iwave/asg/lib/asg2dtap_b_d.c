/*        Generated by TAPENADE     (INRIA, Ecuador team)
	  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
*/
#include "cstd.h"

/*
  Differentiation of asg2dtap_b in forward (tangent) mode:
  variations   of useful results: **v1b **p0b **p0 **p1 **p1b
  **v0b
  with respect to varying inputs: **v0 **v1 ***bulk **v1b **p0b
  **p0 **p1 **p1b **v0b
  RW status of diff variables: **v0:in **v1:in ***bulk:in **v1b:in-out
  **p0b:in-out **p0:in-out **p1:in-out **p1b:in-out
  **v0b:in-out
  Plus diff mem management of: gradp0b:in v0:in *v0:in v1:in
  *v1:in sdivb0:in *sdivb0:in bulk:in *bulk:in **bulk:in
  v1b:in *v1b:in p0b:in *p0b:in p0:in *p0:in p1:in
  *p1:in gradp1b:in sdiv:in *sdiv:in p1b:in *p1b:in
  sdivbb:in *sdivbb:in sdivb:in *sdivb:in v0b:in
  *v0b:in


  Differentiation of asg2dtap in reverse (adjoint) mode:
  gradient     of useful results: **v0 **v1 **p0 **p1
  with respect to varying inputs: **v0 **v1 ***bulk **buoy **p0
  **p1
  RW status of diff variables: **v0:in-out **v1:in-out ***bulk:out
  **buoy:out **p0:in-out **p1:in-out
  Plus diff mem management of: v0:in *v0:in v1:in *v1:in bulk:in
  *bulk:in **bulk:in buoy:in *buoy:in p0:in *p0:in
  p1:in *p1:in sdiv:in *sdiv:in sdivb:in *sdivb:in
  gradp0:in gradp1:in

  Generated by TAPENADE     (INRIA, Ecuador team)
  Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
*/
void asg2dtap_b_d(float *** restrict bulk, float *** restrict bulkd, float ***
		  restrict bulkb, float ** restrict buoy, float ** restrict buoyb, float
		  ** restrict p0, float ** restrict p0d, float ** restrict p0b, float **
		  restrict p0bd, float ** restrict p1, float ** restrict p1d, float ** 
		  restrict p1b, float ** restrict p1bd, float ** restrict v0, float ** 
		  restrict v0d, float ** restrict v0b, float ** restrict v0bd, float ** 
		  restrict v1, float ** restrict v1d, float ** restrict v1b, float ** 
		  restrict v1bd, float ** restrict ep, float ** restrict epp, float ** 
		  restrict ev, float ** restrict evp, float ** restrict sdiv, float ** 
		  restrict sdivd, float ** restrict sdivb0, float ** restrict sdivb0d, 
		  float ** restrict sdivb, float ** restrict sdivbd, float ** restrict 
		  sdivbb, float ** restrict sdivbbd, float * restrict gradp0, float * 
		  restrict gradp0b, float * restrict gradp0bd, float * restrict gradp1, 
		  float * restrict gradp1b, float * restrict gradp1bd, int *gsc_p, int *
		  gec_p, int *gsc_v0, int *gec_v0, int *gsc_v1, int *gec_v1, int *lbc, 
		  int *rbc, int ihmax, int maxoff, float ** restrict c, FILE *stream) {
  int i0, i1;
  int ioff;

  {
    int ih;
    int sh;
    int eh;
    float tempb;
    float tempbd;
    float tempb0;
    float tempb0d;
    float tempb1;
    float tempb1d;
    float tempb2;
    float tempb2d;
    float tmp;
    float tmpd;
    float tmpb;
    float tmpbd;
    float tmp0;
    float tmp0d;
    float tmpb0;
    float tmpb0d;
    float tmp1;
    float tmp1d;
    float tmpb1;
    float tmpb1d;
    float tmp2;
    float tmp2d;
    float tmpb2;
    float tmpb2d;
    float tempb3;
    float tempb3d;
    float tempb4;
    float tempb4d;
    float tempb5;
    float tempb6;
    float tempb6d;
    float tempb7;
    float tempb7d;
    float tempb8;
    float tmpb3;
    float tmpb3d;
    float tmpb4;
    float tmpb4d;
    float tmpb5;
    float tmpb5d;
    float tmpb6;
    float tmpb6d;

    //      **sdivd = 0.0;
    // fprintf(stderr,"1\n");
    for (i1 = gsc_p[1]; i1 < gec_p[1]+1; ++i1) {
      //      // fprintf(stderr,"1a i1=%d gsc_p[1]=%d gec_p[1]=%d\n",i1,gsc_p[1],gec_p[1]);
      /*
	memset_d((char *)sdiv[i1 - gsc_p[1]], (char *)sdivd[i1 - gsc_p[1]], 
	0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
	memset_d0((char *)sdivb[i1 - gsc_p[1]], (char *)sdivbd[i1 - gsc_p[1]
	], 0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
      */
      memset((char *)sdiv[i1 - gsc_p[1]],   0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
      memset((char *)sdivbd[i1 - gsc_p[1]], 0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
      memset((char *)sdivd[i1 - gsc_p[1]],  0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
      memset((char *)sdivb[i1 - gsc_p[1]],  0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
      //      // fprintf(stderr,"1b\n");      
      for (ioff = 0; ioff < maxoff; ++ioff)
#pragma ivdep
	for (i0 = gsc_p[0]; i0 < gec_p[0]+1; ++i0) {
	  sdivd[i1 - gsc_p[1]][i0 - gsc_p[0]] = sdivd[i1 - gsc_p[1]][
								     i0 - gsc_p[0]] + c[0][ioff]*(v0d[i1][i0+ioff]-v0d[i1][i0
															   -ioff-1]) + c[1][ioff]*(v1d[i1+ioff][i0]-v1d[i1-ioff-1][
																						   i0]);
	  sdiv[i1 - gsc_p[1]][i0 - gsc_p[0]] += c[0][ioff]*(v0[i1][i0+
								   ioff]-v0[i1][i0-ioff-1]) + c[1][ioff]*(v1[i1+ioff][i0]-v1[i1
															     -ioff-1][i0]);
	}
    }
    //      **sdivbd = 0.0;
    // fprintf(stderr,"2\n");    
    for (ih = -ihmax; ih < ihmax+1; ++ih) {
      //      sh = ((ih <= 0) ?  gsc_p[1]-ih : gsc_p[1]);
      //      eh = ((ih >= 0) ?  gec_p[1]-ih : gec_p[1]);
      if (ih <= 0) {
	sh = gsc_p[1] - ih;
	eh = gec_p[1];
      } else {
	sh = gsc_p[1];
	eh = gec_p[1] - ih;
      }
      for (i1 = sh; i1 < eh+1; ++i1)
#pragma ivdep
	for (i0 = gsc_p[0]; i0 < gec_p[0]+1; ++i0) {
	  sdivbd[i1 - gsc_p[1]][i0 - gsc_p[0]] = sdivbd[i1 - gsc_p[1]]
	    [i0 - gsc_p[0]] + bulkd[ih+ihmax][i1][i0]*sdiv[i1-gsc_p[
								    1]+ih][i0-gsc_p[0]] + bulk[ih+ihmax][i1][i0]*sdivd[i1-
														       gsc_p[1]+ih][i0-gsc_p[0]];
	  sdivb[i1 - gsc_p[1]][i0 - gsc_p[0]] += bulk[ih+ihmax][i1][i0
								    ]*sdiv[i1-gsc_p[1]+ih][i0-gsc_p[0]];
	}
    }
    // fprintf(stderr,"3\n");
    for (i1 = gsc_p[1]; i1 < gec_p[1]+1; ++i1)
#pragma ivdep
      for (i0 = gsc_p[0]; i0 < gec_p[0]+1; ++i0) {
	p0d[i1][i0] = epp[0][i0]*(ep[0][i0]*p0d[i1][i0]-sdivbd[i1-gsc_p[
									1]][i0-gsc_p[0]]);
	p0[i1][i0] = (ep[0][i0]*p0[i1][i0]-sdivb[i1-gsc_p[1]][i0-gsc_p[0
								       ]])*epp[0][i0];
	p1d[i1][i0] = epp[1][i1]*(ep[1][i1]*p1d[i1][i0]-sdivbd[i1-gsc_p[
									1]][i0-gsc_p[0]]);
	p1[i1][i0] = (ep[1][i1]*p1[i1][i0]-sdivb[i1-gsc_p[1]][i0-gsc_p[0
								       ]])*epp[1][i1];
      }
    // fprintf(stderr,"4\n");
    if (lbc[0])
      for (i1 = gsc_p[1]; i1 < gec_p[1]+1; ++i1) {
	p0d[i1][gsc_p[0] - 1] = 0.0;
	p0[i1][gsc_p[0] - 1] = 0.0f;
	for (ioff = 1; ioff < maxoff; ++ioff) {
	  tmpd = -p0d[i1][gsc_p[0]+ioff-1];
	  tmp = -p0[i1][gsc_p[0]+ioff-1];
	  p0d[i1][gsc_p[0] - ioff - 1] = tmpd;
	  p0[i1][gsc_p[0] - ioff - 1] = tmp;
	}
      }
    if (rbc[0])
      for (i1 = gsc_p[1]; i1 < gec_p[1]+1; ++i1) {
	p0d[i1][gec_p[0] + 1] = 0.0;
	p0[i1][gec_p[0] + 1] = 0.0f;
	for (ioff = 1; ioff < maxoff; ++ioff) {
	  tmp0d = -p0d[i1][gec_p[0]-ioff+1];
	  tmp0 = -p0[i1][gec_p[0]-ioff+1];
	  p0d[i1][gec_p[0] + ioff + 1] = tmp0d;
	  p0[i1][gec_p[0] + ioff + 1] = tmp0;
	}
      }
    if (lbc[1])
      for (ioff = 0; ioff < maxoff; ++ioff)
#pragma ivdep
	for (i0 = gsc_p[0]; i0 < gec_p[0]+1; ++i0) {
	  tmp1d = -p1d[gsc_p[1]+ioff-1][i0];
	  tmp1 = -p1[gsc_p[1]+ioff-1][i0];
	  p1d[gsc_p[1] - ioff - 1][i0] = (ioff < 1 ?  0.0f : tmp1d);
	  p1[gsc_p[1] - ioff - 1][i0] = (ioff < 1 ? 0.0f : tmp1);
	}
    if (rbc[1])
      for (ioff = 1; ioff < maxoff; ++ioff)
#pragma ivdep
	for (i0 = gsc_p[0]; i0 < gec_p[0]+1; ++i0) {
	  tmp2d = -p1d[gec_p[1]-ioff+1][i0];
	  tmp2 = -p1[gec_p[1]-ioff+1][i0];
	  p1d[gec_p[1] + ioff + 1][i0] = (ioff < 1 ?  0.0f : tmp2d);
	  p1[gec_p[1] + ioff + 1][i0] = (ioff < 1 ? 0.0f : tmp2);
	}
    // fprintf(stderr,"6\n");
    if (rbc[1])
      for (ioff = maxoff-1; ioff > 0; --ioff)
#pragma ivdep
	for (i0 = gec_v1[0]; i0 > gsc_v1[0]-1; --i0) {
	  tmpb6d = v1bd[gec_v1[1] + ioff][i0];
	  tmpb6 = v1b[gec_v1[1] + ioff][i0];
	  v1bd[gec_v1[1] + ioff][i0] = 0.0;
	  v1b[gec_v1[1] + ioff][i0] = 0.0;
	  v1bd[gec_v1[1] - ioff + 1][i0] = v1bd[gec_v1[1] - ioff + 1][
								      i0] + tmpb6d;
	  v1b[gec_v1[1] - ioff + 1][i0] = v1b[gec_v1[1] - ioff + 1][i0
								    ] + tmpb6;
	}
    if (lbc[1])
      for (ioff = maxoff-1; ioff > 0; --ioff)
#pragma ivdep
	for (i0 = gec_v1[0]; i0 > gsc_v1[0]-1; --i0) {
	  tmpb5d = v1bd[gsc_v1[1] - ioff][i0];
	  tmpb5 = v1b[gsc_v1[1] - ioff][i0];
	  v1bd[gsc_v1[1] - ioff][i0] = 0.0;
	  v1b[gsc_v1[1] - ioff][i0] = 0.0;
	  v1bd[gsc_v1[1] + ioff - 1][i0] = v1bd[gsc_v1[1] + ioff - 1][
								      i0] + tmpb5d;
	  v1b[gsc_v1[1] + ioff - 1][i0] = v1b[gsc_v1[1] + ioff - 1][i0
								    ] + tmpb5;
	}
    // fprintf(stderr,"8\n");
    if (rbc[0])
      for (ioff = maxoff-1; ioff > 0; --ioff)
	for (i1 = gec_v0[1]; i1 > gsc_v0[1]-1; --i1) {
	  tmpb4d = v0bd[i1][gec_v0[0] + ioff];
	  tmpb4 = v0b[i1][gec_v0[0] + ioff];
	  v0bd[i1][gec_v0[0] + ioff] = 0.0;
	  v0b[i1][gec_v0[0] + ioff] = 0.0;
	  v0bd[i1][gec_v0[0] - ioff + 1] = v0bd[i1][gec_v0[0] - ioff +
						    1] + tmpb4d;
	  v0b[i1][gec_v0[0] - ioff + 1] = v0b[i1][gec_v0[0] - ioff + 1
						  ] + tmpb4;
	}
    // fprintf(stderr,"9\n");
    if (lbc[0]) {
      for (ioff = maxoff-1; ioff > 0; --ioff)
	for (i1 = gec_v0[1]; i1 > gsc_v0[1]-1; --i1) {
	  tmpb3d = v0bd[i1][gsc_v0[0] - ioff];
	  tmpb3 = v0b[i1][gsc_v0[0] - ioff];
	  v0bd[i1][gsc_v0[0] - ioff] = 0.0;
	  v0b[i1][gsc_v0[0] - ioff] = 0.0;
	  v0bd[i1][gsc_v0[0] + ioff - 1] = v0bd[i1][gsc_v0[0] + ioff -
						    1] + tmpb3d;
	  v0b[i1][gsc_v0[0] + ioff - 1] = v0b[i1][gsc_v0[0] + ioff - 1
						  ] + tmpb3;
	}
    }

    // fprintf(stderr,"10\n");
    for (i1 = gec_v1[1]; i1 > gsc_v1[1]-1; --i1) {
      memset((char *)gradp1, 0, (gec_v1[0]-gsc_v1[0]+1)*sizeof(float));
      /*
	memset_d((char *)gradp1b, (char *)gradp1bd, 0, (gec_v1[0]-gsc_v1[0]+
	1)*sizeof(float));
      */
      memset((char *)gradp1b,  0, (gec_v1[0]-gsc_v1[0]+1)*sizeof(float));
      memset((char *)gradp1bd, 0, (gec_v1[0]-gsc_v1[0]+1)*sizeof(float));
      for (ioff = 0; ioff < maxoff; ++ioff)
#pragma ivdep
	for (i0 = gsc_v1[0]; i0 < gec_v1[0]+1; ++i0)
	  gradp1[i0 - gsc_v1[0]] = gradp1[i0 - gsc_v1[0]] + c[1][ioff]
	    *(p1[i1+ioff+1][i0]-p1[i1-ioff][i0]);
      //#pragma ivdep
      for (i0 = gec_v1[0]; i0 > gsc_v1[0]-1; --i0) {
	tempb7d = evp[1][i1]*v1bd[i1][i0];
	tempb7 = evp[1][i1]*v1b[i1][i0];
	tempb8 = -(0.5f*gradp1[i0-gsc_v1[0]]*tempb7);
	buoyb[i1][i0] = buoyb[i1][i0] + tempb8;
	buoyb[i1 + 1][i0] = buoyb[i1 + 1][i0] + tempb8;
	gradp1bd[i0 - gsc_v1[0]] = gradp1bd[i0 - gsc_v1[0]] - 0.5f*(buoy
								    [i1][i0]+buoy[i1+1][i0])*tempb7d;
	gradp1b[i0 - gsc_v1[0]] = gradp1b[i0 - gsc_v1[0]] - 0.5f*(buoy[
								       i1][i0]+buoy[i1+1][i0])*tempb7;
	v1bd[i1][i0] = ev[1][i1]*tempb7d;
	v1b[i1][i0] = ev[1][i1]*tempb7;
      }
      for (ioff = maxoff-1; ioff > -1; --ioff)
#pragma ivdep
	for (i0 = gec_v1[0]; i0 > gsc_v1[0]-1; --i0) {
	  tempb6d = c[1][ioff]*gradp1bd[i0-gsc_v1[0]];
	  tempb6 = c[1][ioff]*gradp1b[i0-gsc_v1[0]];
	  p1bd[i1 + ioff + 1][i0] = p1bd[i1 + ioff + 1][i0] + tempb6d;
	  p1b[i1 + ioff + 1][i0] = p1b[i1 + ioff + 1][i0] + tempb6;
	  p1bd[i1 - ioff][i0] = p1bd[i1 - ioff][i0] - tempb6d;
	  p1b[i1 - ioff][i0] = p1b[i1 - ioff][i0] - tempb6;
	}
    }
    //      *gradp0bd = 0.0;
    // fprintf(stderr,"12\n");
    for (i1 = gec_v0[1]; i1 > gsc_v0[1]-1; --i1) {
      memset((char *)gradp0, 0, (gec_v0[0]-gsc_v0[0]+1)*sizeof(float));
      /*
	memset_d((char *)gradp0b, (char *)gradp0bd, 0, (gec_v0[0]-gsc_v0[0]+
	1)*sizeof(float));
      */
      memset((char *)gradp0b,  0, (gec_v0[0]-gsc_v0[0]+1)*sizeof(float));
      memset((char *)gradp0bd, 0, (gec_v0[0]-gsc_v0[0]+1)*sizeof(float));
      for (ioff = 0; ioff < maxoff; ++ioff)
#pragma ivdep
	for (i0 = gsc_v0[0]; i0 < gec_v0[0]+1; ++i0)
	  gradp0[i0 - gsc_v0[0]] = gradp0[i0 - gsc_v0[0]] + c[0][ioff]
	    *(p0[i1][i0+ioff+1]-p0[i1][i0-ioff]);
#pragma ivdep
      for (i0 = gec_v0[0]; i0 > gsc_v0[0]-1; --i0) {
	tempb4d = evp[0][i0]*v0bd[i1][i0];
	tempb4 = evp[0][i0]*v0b[i1][i0];
	tempb5 = -(0.5f*gradp0[i0-gsc_v0[0]]*tempb4);
	buoyb[i1][i0] = buoyb[i1][i0] + tempb5;
	buoyb[i1][i0 + 1] = buoyb[i1][i0 + 1] + tempb5;
	gradp0bd[i0 - gsc_v0[0]] = gradp0bd[i0 - gsc_v0[0]] - 0.5f*(buoy
								    [i1][i0]+buoy[i1][i0+1])*tempb4d;
	gradp0b[i0 - gsc_v0[0]] += -0.5f*(buoy[i1][i0]+buoy[i1][i0+1])*
	  tempb4;
	v0bd[i1][i0] = ev[0][i0]*tempb4d;
	v0b[i1][i0] = ev[0][i0]*tempb4;
      }
      for (ioff = maxoff-1; ioff > -1; --ioff)
#pragma ivdep
	for (i0 = gec_v0[0]; i0 > gsc_v0[0]-1; --i0) {
	  tempb3d = c[0][ioff]*gradp0bd[i0-gsc_v0[0]];
	  tempb3 = c[0][ioff]*gradp0b[i0-gsc_v0[0]];
	  p0bd[i1][i0 + ioff + 1] = p0bd[i1][i0 + ioff + 1] + tempb3d;
	  p0b[i1][i0 + ioff + 1] = p0b[i1][i0 + ioff + 1] + tempb3;
	  p0bd[i1][i0 - ioff] = p0bd[i1][i0 - ioff] - tempb3d;
	  p0b[i1][i0 - ioff] = p0b[i1][i0 - ioff] - tempb3;
	}
    }
    // fprintf(stderr,"13\n");
    if (rbc[1])
      for (ioff = maxoff-1; ioff > 0; --ioff) {
#pragma ivdep
	for (i0 = gec_p[0]; i0 > gsc_p[0]-1; --i0) {
	  tmpb2d = p1bd[gec_p[1] + ioff + 1][i0];
	  tmpb2 = p1b[gec_p[1] + ioff + 1][i0];
	  p1bd[gec_p[1] + ioff + 1][i0] = 0.0;
	  p1b[gec_p[1] + ioff + 1][i0] = 0.0;
	  p1bd[gec_p[1] - ioff + 1][i0] = (ioff < 1 ? 0.0f : p1bd[gec_p[1]
							      - ioff + 1][i0] - tmpb2d);
	  p1b[gec_p[1] - ioff + 1][i0] = (ioff < 1 ? 0.0f : p1b[gec_p[
								      1] - ioff + 1][i0] - tmpb2);
	}
      }
#pragma ivdep
    for (i0 = gec_p[0]; i0 > gsc_p[0]-1; --i0) {      
      p1bd[gec_p[1] + 1][i0] = 0.0;
      p1b[gec_p[1] + 1][i0] = 0.0;
    }
      
    // fprintf(stderr,"14\n");
    if (lbc[1])
      for (ioff = maxoff-1; ioff > 0; --ioff) {
#pragma ivdep	
	for (i0 = gec_p[0]; i0 > gsc_p[0]-1; --i0) {
	  tmpb1d = p1bd[gsc_p[1] - ioff - 1][i0];
	  tmpb1 = p1b[gsc_p[1] - ioff - 1][i0];
	  p1bd[gsc_p[1] - ioff - 1][i0] = 0.0;
	  p1b[gsc_p[1] - ioff - 1][i0] = 0.0;
	  p1bd[gsc_p[1] + ioff - 1][i0] = p1bd[gsc_p[1] + ioff - 1][i0
								    ] - tmpb1d;
	  p1b[gsc_p[1] + ioff - 1][i0] = p1b[gsc_p[1] + ioff - 1][i0] 
	    - tmpb1;
	}
      }
#pragma ivdep	
    for (i0 = gec_p[0]; i0 > gsc_p[0]-1; --i0) {	  
      p1bd[gsc_p[1] - 1][i0] = 0.0;
      p1b[gsc_p[1] - 1][i0] = 0.0;
    }
    // fprintf(stderr,"15\n");
    if (rbc[0])
      for (i1 = gec_p[1]; i1 > gsc_p[1]-1; --i1) {
	for (ioff = maxoff-1; ioff > 0; --ioff) {
	  tmpb0d = p0bd[i1][gec_p[0] + ioff + 1];
	  tmpb0 = p0b[i1][gec_p[0] + ioff + 1];
	  p0bd[i1][gec_p[0] + ioff + 1] = 0.0;
	  p0b[i1][gec_p[0] + ioff + 1] = 0.0;
	  p0bd[i1][gec_p[0] - ioff + 1] = p0bd[i1][gec_p[0] - ioff + 1
						   ] - tmpb0d;
	  p0b[i1][gec_p[0] - ioff + 1] = p0b[i1][gec_p[0] - ioff + 1] 
	    - tmpb0;
	}
	p0bd[i1][gec_p[0] + 1] = 0.0;
	p0b[i1][gec_p[0] + 1] = 0.0;
      }
    // fprintf(stderr,"16\n");
    if (lbc[0]) {
      for (i1 = gec_p[1]; i1 > gsc_p[1]-1; --i1) {
	for (ioff = maxoff-1; ioff > 0; --ioff) {
	  tmpbd = p0bd[i1][gsc_p[0] - ioff - 1];
	  tmpb = p0b[i1][gsc_p[0] - ioff - 1];
	  p0bd[i1][gsc_p[0] - ioff - 1] = 0.0;
	  p0b[i1][gsc_p[0] - ioff - 1] = 0.0;
	  p0bd[i1][gsc_p[0] + ioff - 1] = p0bd[i1][gsc_p[0] + ioff - 1
						   ] - tmpbd;
	  p0b[i1][gsc_p[0] + ioff - 1] = p0b[i1][gsc_p[0] + ioff - 1] 
	    - tmpb;
	}
	p0bd[i1][gsc_p[0] - 1] = 0.0;
	p0b[i1][gsc_p[0] - 1] = 0.0;
      }
      //          **sdivbbd = 0.0;
    } //else
      //          **sdivbbd = 0.0;
      // fprintf(stderr,"17\n");
    for (i1 = gec_p[1]; i1 > gsc_p[1]-1; --i1) {
      /*
	memset_d((char *)sdivbb[i1 - gsc_p[1]], (char *)sdivbbd[i1 - gsc_p[1
	]], 0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
      */
      /*
	memset_d0((char *)sdivb0[i1 - gsc_p[1]], (char *)sdivb0d[i1 - gsc_p[
	1]], 0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
      */
      memset((char *)sdivbb[i1 - gsc_p[1]], 0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
      memset((char *)sdivbbd[i1 - gsc_p[1]], 0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
      memset((char *)sdivb0[i1 - gsc_p[1]],  0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
      memset((char *)sdivb0d[i1 - gsc_p[1]], 0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
      // fprintf(stderr,"18\n");
#pragma ivdep
      for (i0 = gec_p[0]; i0 > gsc_p[0]-1; --i0) {
	tempb1d = epp[1][i1]*p1bd[i1][i0];
	tempb1 = epp[1][i1]*p1b[i1][i0];
	sdivbbd[i1 - gsc_p[1]][i0 - gsc_p[0]] = sdivbbd[i1 - gsc_p[1]][
								       i0 - gsc_p[0]] - tempb1d;
	sdivbb[i1 - gsc_p[1]][i0 - gsc_p[0]] = sdivbb[i1 - gsc_p[1]][i0 
								     - gsc_p[0]] - tempb1;
	p1bd[i1][i0] = ep[1][i1]*tempb1d;
	p1b[i1][i0] = ep[1][i1]*tempb1;
	tempb2d = epp[0][i0]*p0bd[i1][i0];
	tempb2 = epp[0][i0]*p0b[i1][i0];
	sdivbbd[i1 - gsc_p[1]][i0 - gsc_p[0]] = sdivbbd[i1 - gsc_p[1]][
								       i0 - gsc_p[0]] - tempb2d;
	sdivbb[i1 - gsc_p[1]][i0 - gsc_p[0]] = sdivbb[i1 - gsc_p[1]][i0 
								     - gsc_p[0]] - tempb2;
	p0bd[i1][i0] = ep[0][i0]*tempb2d;
	p0b[i1][i0] = ep[0][i0]*tempb2;
      }
    }
    //      **sdivb0d = 0.0;
    // fprintf(stderr,"19\n");
    for (ih = ihmax; ih > -1-ihmax; --ih) {
      if (ih <= 0) {
	sh = gsc_p[1] - ih;
	eh = gec_p[1];
      } else {
	sh = gsc_p[1];
	eh = gec_p[1] - ih;
      }
      for (i1 = eh; i1 > sh-1; --i1)
#pragma ivdep
	for (i0 = gec_p[0]; i0 > gsc_p[0]-1; --i0) {
	  bulkb[ih + ihmax][i1][i0] +=
	    sdiv[i1-gsc_p[1]+ih][i0-gsc_p[0]]*sdivbb[i1-gsc_p[1]][i0-gsc_p[0]];
	  sdivb0d[i1 - gsc_p[1] + ih][i0 - gsc_p[0]] +=
	    bulkd[ih+ihmax][i1][i0]*sdivbb[i1-gsc_p[1]][i0-gsc_p[0]] +
	    bulk[ih+ihmax][i1][i0]*sdivbbd[i1-gsc_p[1]][i0-gsc_p[0]];
	  sdivb0[i1 - gsc_p[1] + ih][i0 - gsc_p[0]] +=
	    bulk[ih+ihmax][i1][i0]*sdivbb[i1-gsc_p[1]][i0-gsc_p[0]];
	}
    }
    for (i1 = gec_p[1]; i1 > gsc_p[1]-1; --i1)
      for (ioff = maxoff-1; ioff > -1; --ioff)
#pragma ivdep
	for (i0 = gec_p[0]; i0 > gsc_p[0]-1; --i0) {
	  tempbd = c[0][ioff]*sdivb0d[i1-gsc_p[1]][i0-gsc_p[0]];
	  tempb = c[0][ioff]*sdivb0[i1-gsc_p[1]][i0-gsc_p[0]];
	  tempb0d = c[1][ioff]*sdivb0d[i1-gsc_p[1]][i0-gsc_p[0]];
	  tempb0 = c[1][ioff]*sdivb0[i1-gsc_p[1]][i0-gsc_p[0]];
	  v0bd[i1][i0 + ioff] = v0bd[i1][i0 + ioff] + tempbd;
	  v0b[i1][i0 + ioff] = v0b[i1][i0 + ioff] + tempb;
	  v0bd[i1][i0 - ioff - 1] = v0bd[i1][i0 - ioff - 1] - tempbd;
	  v0b[i1][i0 - ioff - 1] = v0b[i1][i0 - ioff - 1] - tempb;
	  v1bd[i1 + ioff][i0] = v1bd[i1 + ioff][i0] + tempb0d;
	  v1b[i1 + ioff][i0] = v1b[i1 + ioff][i0] + tempb0;
	  v1bd[i1 - ioff - 1][i0] = v1bd[i1 - ioff - 1][i0] - tempb0d;
	  v1b[i1 - ioff - 1][i0] = v1b[i1 - ioff - 1][i0] - tempb0;
	}
  }
}
