#include "cstd.h"

#define ABSORB
//#define LOCAL_VERBOSE
#undef LOCAL_VERBOSE

void asg_vstep2d(float ** restrict buoy,
                 float ** restrict p0,float ** restrict p1,
                 float ** restrict v0,float ** restrict v1,
                 float ** restrict ev, float ** restrict evp,
                 float ** restrict gradp,
                 int * gsc_v0, int * gec_v0,
                 int * gsc_v1, int * gec_v1,
                 int * lbc, int * rbc,
                 int maxoff,float ** restrict c,
                 FILE * stream) {
  /*  
    // signature for tapenade

     void asg_vstep2d(float **buoy,
     float **p0,float **p1,
     float **v0,float **v1,
     float **ev, float **evp,
     float **gradp,
     int * gsc_v0, int * gec_v0,
     int * gsc_v1, int * gec_v1,
     int * lbc, int * rbc,
     int maxoff,float **c) {
  */
    int i0, i1;
    int ioff;
    
#ifdef LOCAL_VERBOSE
    fprintf(stream,"asg_vstep2d v0\n");
    fprintf(stream,"idim=0 gsc_v0=%d gec_v0=%d\n",gsc_v0[0],gec_v0[0]);
    fprintf(stream,"idim=1 gsc_v0=%d gec_v0=%d\n",gsc_v0[1],gec_v0[1]);
    fprintf(stream,"maxoff=%d\n",maxoff);
    fflush(stream);
#endif
    
    for (i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++ ) {
        
      for (i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ ) gradp[0][i0]=0;
      
      for (ioff=0; ioff<maxoff; ioff++) {
	for (i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ )
	  gradp[0][i0] += c[0][ioff]*(p0[i1][i0+ioff+1]-p0[i1][i0-ioff]);
      }
      
      for (i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ ) {
	v0[i1][i0] =
#ifdef ABSORB
	  evp[0][i0]*
	  (v0[i1][i0]*ev[0][i0] - 0.5*(buoy[i1][i0]+buoy[i1][i0+1])*gradp[0][i0]);
#else
          v0[i1][i0] - 0.5*(buoy[i1][i0]+buoy[i1][i0+1])*gradp[0][i0];
#endif
      }
    }
    
#ifdef LOCAL_VERBOSE
    fprintf(stream,"asg_vstep2d v1\n");
    fprintf(stream,"idim=0 gsc_v1=%d gec_v1=%d\n",gsc_v1[0],gec_v1[0]);
    fprintf(stream,"idim=1 gsc_v1=%d gec_v1=%d\n",gsc_v1[1],gec_v1[1]);
    fprintf(stream,"maxoff=%d\n",maxoff);
    fflush(stream);
#endif

    for (i1=gsc_v1[1]; i1 <= gec_v1[1]; i1++ ) {
      for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) gradp[1][i0]=0;
      for (ioff=0; ioff<maxoff; ioff++) {
	for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) {
	  gradp[1][i0] += c[1][ioff]*(p1[i1+ioff+1][i0]-p1[i1-ioff][i0]);
        }
      }
      for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) {
	v1[i1][i0] =
#ifdef ABSORB
	  evp[1][i1]*
	  (v1[i1][i0]*ev[1][i1] - 0.5*(buoy[i1][i0]+buoy[i1+1][i0])*gradp[1][i0]);
#else
	  v1[i1][i0] - 0.5*(buoy[i1][i0]+buoy[i1+1][i0])*gradp[1][i0];
#endif
      }
    }
#ifdef LOCAL_VERBOSE
    fprintf(stream,"local p at depth 2:\n");
    fprintf(stream,"idim=0 gsc_v0=%d gsc_v1=%d p0[gsc_v0][2]=%g p1[gsc_v1][2]=%g p1[gsc_v1+1][2]=%g\n",gsc_v0[1],gsc_v1[1],p0[gsc_v0[1]][2],p1[gsc_v1[1]][2],p1[gsc_v1[1]+1][2]);
    fprintf(stream,"idim=1 gec_v0=%d gec_v1=%d p0[gec_v0][2]=%g p1[gec_v1][2]=%g p1[gec_v1+1][2]=%g\n",gec_v0[1],gec_v1[1],p0[gec_v0[1]][2],p1[gec_v1[1]][2],p1[gec_v1[1]+1][2]);
    fprintf(stream,"local v0 at depth 2:\n");
    fprintf(stream,"idim=0 gsc_v0=%d v0[gsc_v0][2]=%g v0[gsc_v0+1][2]=%g\n",gsc_v0[1],v0[gsc_v0[1]][2],v0[gsc_v0[1]+1][2]);
    fprintf(stream,"idim=0 gec_v0=%d v0[gec_v0][2]=%g v0[gec_v0-1][2]=%g\n",gec_v0[1],v0[gec_v0[1]][2],v0[gec_v0[1]-1][2]);
    fprintf(stream,"local v1 at depth 2:\n");
    fprintf(stream,"idim=0 gsc_v1=%d v1[gsc_v1][2]=%g v1[gsc_v1+1][2]=%g\n",gsc_v1[1],v1[gsc_v1[1]][2],v1[gsc_v1[1]+1][2]);
    fprintf(stream,"idim=0 gec_v1=%d v1[gec_v1][2]=%g v1[gec_v1-1][2]=%g\n",gec_v1[1],v1[gec_v1[1]][2],v1[gec_v1[1]-1][2]);
#endif
    
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
    
#ifdef LOCAL_VERBOSE
    fprintf(stream,"asg_vstep2d end\n");
    fflush(stream);
#endif
}

