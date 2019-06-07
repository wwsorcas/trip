#include <stdio.h>

/* signature for tapenade */
void asg2d(float ** bulk, float ** buoy,
	   float ** p0, float ** p1,
	   float ** v0, float ** v1,
	   float ** ep, float ** epp,
	   float ** ev, float ** evp,
	   float * sdiv,
	   float ** gradp,
	   int * gsc_p, int * gec_p, 
	   int * gsc_v0, int * gec_v0,
	   int * gsc_v1, int * gec_v1,
	   int * lbc, int * rbc,
	   int maxoff, float ** c,
	   FILE * stream) {

  int i0, i1;
  int ioff;

  for (i1=gsc_p[1]; i1 <= gec_p[1]; i1++) {
    for (i0=gsc_p[0]; i0 <= gec_p[0]; i0++) {
      sdiv[i0] = 
	(c[0][0]*(v0[i1][i0]-v0[i1][i0-1]) +
	 c[1][0]*(v1[i1][i0]-v1[i1-1][i0]));
    }
    for (ioff = 1; ioff<maxoff; ioff++) {
      for (i0=gsc_p[0]; i0 <= gec_p[0]; i0++) {
	sdiv[i0] += 
	  (c[0][ioff]*(v0[i1][i0+ioff]-v0[i1][i0-ioff-1]) +
	   c[1][ioff]*(v1[i1+ioff][i0]-v1[i1-ioff-1][i0]));
      }
    }
    for (i0=gsc_p[0]; i0 <= gec_p[0]; i0++) {
      p0[i1][i0] = (p0[i1][i0]*ep[0][i0] - bulk[i1][i0]*sdiv[i0])*epp[0][i0];
      p1[i1][i0] = (p1[i1][i0]*ep[1][i1] - bulk[i1][i0]*sdiv[i0])*epp[1][i1];
    }
  }

  if (lbc[0]) {
    for (i1=gsc_p[1];i1<=gec_p[1];i1++) {
      p0[i1][gsc_p[0]-1]=0.0f;
      for (ioff=1;ioff<maxoff;ioff++) {
	p0[i1][gsc_p[0]-ioff-1]=-p0[i1][gsc_p[0]+ioff-1];
      }
    }
  }
  if (rbc[0]) {
    for (i1=gsc_p[1];i1<=gec_p[1];i1++) {
      p0[i1][gec_p[0]+1]=0.0f;
      for (ioff=1;ioff<maxoff;ioff++) {
	p0[i1][gec_p[0]+ioff+1]=-p0[i1][gec_p[0]-ioff+1];
      }
    }
  }
  if (lbc[1]) {
    for (i0=gsc_p[0];i0<=gec_p[0];i0++) {
      p1[gsc_p[1]-1][i0]=0.0f;
      for (ioff=1;ioff<maxoff;ioff++) {
	p1[gsc_p[1]-ioff-1][i0]=-p1[gsc_p[1]+ioff-1][i0];
      }
    }
  }
  if (rbc[1]) {
    for (i0=gsc_p[0];i0<=gec_p[0];i0++) {
      p1[gec_p[1]+1][i0]=0.0f;
      for (ioff=1;ioff<maxoff;ioff++) {
	p1[gec_p[1]+ioff+1][i0]=-p1[gec_p[1]-ioff+1][i0];
      }
    }
  }

  for (i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++ ) {
    
    for (i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ ) gradp[0][i0]=0;
    
    for (ioff=0; ioff<maxoff; ioff++) {
      for (i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ )
	gradp[0][i0] += c[0][ioff]*(p0[i1][i0+ioff+1]-p0[i1][i0-ioff]);
    }
    
    for (i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ ) {
      v0[i1][i0] =
	evp[0][i0]*
	(v0[i1][i0]*ev[0][i0] - 0.5*(buoy[i1][i0]+buoy[i1][i0+1])*gradp[0][i0]);
    }
  }
  
  for (i1=gsc_v1[1]; i1 <= gec_v1[1]; i1++ ) {
    for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) gradp[1][i0]=0;
    for (ioff=0; ioff<maxoff; ioff++) {
      for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) {
	gradp[1][i0] += c[1][ioff]*(p1[i1+ioff+1][i0]-p1[i1-ioff][i0]);
      }
    }
    for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) {
      v1[i1][i0] =
	evp[1][i1]*
	(v1[i1][i0]*ev[1][i1] - 0.5*(buoy[i1][i0]+buoy[i1+1][i0])*gradp[1][i0]);
    }
  }
    
  if (lbc[0]) {
    for (i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++) {
      for (ioff=1; ioff<maxoff;ioff++) {
	v0[i1][gsc_v0[0]-ioff]=v0[i1][gsc_v0[0]+ioff-1];
      }
    }
  }
  if (rbc[0]) {
    for (i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++) {
      for (ioff=1; ioff<maxoff;ioff++) {
	v0[i1][gec_v0[0]+ioff]=v0[i1][gec_v0[0]-ioff+1];
      }
    }
  }
  if (lbc[1]) {
    for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++) {
      for (ioff=1; ioff<maxoff;ioff++) {
	v1[gsc_v1[1]-ioff][i0]=v1[gsc_v1[1]+ioff-1][i0];
      }
    }
  }
  if (rbc[1]) {
    for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++) {
      for (ioff=1; ioff<maxoff;ioff++) {
	v1[gec_v1[1]+ioff][i0]=v1[gec_v1[1]-ioff+1][i0];
      }
    }
  }
    
}
