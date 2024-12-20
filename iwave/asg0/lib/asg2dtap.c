#include "cstd.h"

void asg2dtap(float *** restrict  bulk, float ** restrict  buoy,
	      float **  restrict  p0,   float ** restrict  p1,
	      float **  restrict  v0,   float ** restrict  v1,
	      float **  restrict  ep,   float ** restrict  epp,
	      float **  restrict  ev,   float ** restrict  evp,
	      float **  restrict  sdiv, float ** restrict  sdivb,
	      float *   restrict  gradp0,
	      float *   restrict  gradp1,
	      int * gsc_p, int * gec_p, 
	      int * gsc_v0, int * gec_v0,
	      int * gsc_v1, int * gec_v1,
	      int * lbc, int * rbc,
	      float dh, int ihmax,
	      int maxoff, float ** restrict  c,
	      FILE * stream) {
  
  int i0, i1;
  int ioff;
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
    int ih;
    int sh;
    int eh;
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
    //    fprintf(stderr,"1\n");
    for (i1=gsc_p[1]; i1<=gec_p[1]; i1++) {
      //      fprintf(stderr,"1a i1=%d gsc_p[1]=%d gec_p[1]=%d\n",i1,gsc_p[1],gec_p[1]);
      memset((char *)(sdiv[i1-gsc_p[1]]),  0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
      memset((char *)(sdivb[i1-gsc_p[1]]), 0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
      /*
      for (i0=gsc_p[0]; i0 <= gec_p[0]; i0++) {
	sdiv[i1-gsc_p[1]][i0-gsc_p[0]] = 0.0f;
	sdivb[i1-gsc_p[1]][i0-gsc_p[0]] = 0.0f;	
      }
      */
      //      fprintf(stderr,"1b\n");      
      for (ioff = 0; ioff<maxoff; ioff++) {
#pragma ivdep
	for (i0=gsc_p[0]; i0 <= gec_p[0]; i0++) {
	  sdiv[i1-gsc_p[1]][i0-gsc_p[0]] +=
	    c[0][ioff]*(v0[i1][i0+ioff]-v0[i1][i0-ioff-1]) +
	    c[1][ioff]*(v1[i1+ioff][i0]-v1[i1-ioff-1][i0]);
	}
      }
    }
    //    fprintf(stderr,"2\n");    
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
    // wws mod 2017.10.29: modify indexing of bulk to reflect fact that
    // origin is not correctly accounted for
    for (ih=-ihmax; ih<=ihmax; ih++) {
      //      sh = ((ih <= 0) ?  gsc_p[1]-ih : gsc_p[1]);
      //      eh = ((ih >= 0) ?  gec_p[1]-ih : gec_p[1]);
      if (ih <=0) {
	sh = gsc_p[1]; eh = gec_p[1] + 2*ih;
      }
      else {
	sh = gsc_p[1] + 2*ih; eh = gec_p[1];
      }
      // 2017.11.11 WWS modify indexing to correctly mirror action of bulkmod
      // in midpoint-offset parametrization
      for (i1=sh; i1<=eh; i1++) {    
#pragma ivdep
	for (i0=gsc_p[0]; i0 <= gec_p[0]; i0++) {      
	  sdivb[i1-gsc_p[1]][i0-gsc_p[0]] += dh * bulk[ih+ihmax][i1-ih][i0] * sdiv[i1-gsc_p[1]-2*ih][i0-gsc_p[0]];
	} 
      }
    }
    //    fprintf(stderr,"3\n");
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
    for (i1=gsc_p[1]; i1<=gec_p[1]; i1++) {      
#pragma ivdep
      for (i0=gsc_p[0]; i0 <= gec_p[0]; i0++) {      
	p0[i1][i0] = (ep[0][i0]*p0[i1][i0] - sdivb[i1-gsc_p[1]][i0-gsc_p[0]])*epp[0][i0];
	p1[i1][i0] = (ep[1][i1]*p1[i1][i0] - sdivb[i1-gsc_p[1]][i0-gsc_p[0]])*epp[1][i1]; 
      }
    }
    
    if (lbc[0]) {
#ifdef _OPENMP
#pragma omp for schedule(static,1)
#endif
      for (i1=gsc_p[1];i1<=gec_p[1];i1++) {
	p0[i1][gsc_p[0]-1]=0.0f;
	for (ioff=1;ioff<maxoff;ioff++) {
	  p0[i1][gsc_p[0]-ioff-1]=-p0[i1][gsc_p[0]+ioff-1];
	}
      }
    }

    if (rbc[0]) {
#ifdef _OPENMP
#pragma omp for schedule(static,1)
#endif
      for (i1=gsc_p[1];i1<=gec_p[1];i1++) {
	p0[i1][gec_p[0]+1]=0.0f;
	for (ioff=1;ioff<maxoff;ioff++) {
	  p0[i1][gec_p[0]+ioff+1]=-p0[i1][gec_p[0]-ioff+1];
	}
      }
    }

    if (lbc[1]) {
#ifdef _OPENMP
#pragma omp for schedule(static,1)
#endif
#pragma ivdep
      for (i0=gsc_p[0];i0<=gec_p[0];i0++) {
	p1[gsc_p[1]-1][i0]=0.0f;
	for (ioff=1;ioff<maxoff;ioff++) {
	  p1[gsc_p[1]-ioff-1][i0]=-p1[gsc_p[1]+ioff-1][i0];
	}
      }
    }

    if (rbc[1]) {
#ifdef _OPENMP
#pragma omp for schedule(static,1)
#endif
#pragma ivdep
      for (i0=gsc_p[0];i0<=gec_p[0];i0++) {
	p1[gec_p[1]+1][i0]=0.0f;
	for (ioff=1;ioff<maxoff;ioff++) {
	  p1[gec_p[1]+ioff+1][i0]=-p1[gec_p[1]-ioff+1][i0];
	}
      }
    }

#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
    for (i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++ ) {
      memset((char *)gradp0, 0, (gec_v0[0]-gsc_v0[0]+1)*sizeof(float));
      /*
      for (i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ ) {
	gradp0[i0-gsc_v0[0]]=0.0f;
      }
      */
      for (ioff=0; ioff<maxoff; ioff++) {
#pragma ivdep	
	for (i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ ) {
	  gradp0[i0-gsc_v0[0]] += c[0][ioff]*(p0[i1][i0+ioff+1]-p0[i1][i0-ioff]);
	}
      }
      for (i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ ) {      
	v0[i1][i0] = evp[0][i0]*
	  (ev[0][i0]*v0[i1][i0] -
	   0.5f*(buoy[i1][i0]+buoy[i1][i0+1])*gradp0[i0-gsc_v0[0]]);
      }
    }
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif      
    for (i1=gsc_v1[1]; i1 <= gec_v1[1]; i1++ ) {
      memset((char *)gradp1, 0, (gec_v1[0]-gsc_v1[0]+1)*sizeof(float));
      /*
      for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) {
	gradp1[i0-gsc_v1[0]]=0.0f;
      }
      */
      for (ioff=0; ioff<maxoff; ioff++) {
#pragma ivdep
	for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) {
	  gradp1[i0-gsc_v1[0]] +=  c[1][ioff]*(p1[i1+ioff+1][i0]-p1[i1-ioff][i0]);
	}
      }
      for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) {      
	v1[i1][i0] = evp[1][i1]*
	  (ev[1][i1]* v1[i1][i0] - 0.5f*(buoy[i1][i0]+buoy[i1+1][i0])*gradp1[i0-gsc_v1[0]]);

      }
    }
    
    if (lbc[0]) {
#ifdef _OPENMP
#pragma omp for schedule(static,1)
#endif
      for (ioff=1; ioff<maxoff;ioff++) {
#pragma ivdep
	for (i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++) {
	  v0[i1][gsc_v0[0]-ioff]=v0[i1][gsc_v0[0]+ioff-1];
	}
      }
    }
    
    if (rbc[0]) {
#ifdef _OPENMP
#pragma omp for schedule(static,1)
#endif
      for (ioff=1; ioff<maxoff;ioff++) {
#pragma ivdep      
	for (i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++) {
	  v0[i1][gec_v0[0]+ioff]=v0[i1][gec_v0[0]-ioff+1];
	}
      }
    }
    
    if (lbc[1]) {
#ifdef _OPENMP
#pragma omp for schedule(static,1)
#endif
      for (ioff=1; ioff<maxoff;ioff++) {
#pragma ivdep
	for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++) {
	  v1[gsc_v1[1]-ioff][i0]=v1[gsc_v1[1]+ioff-1][i0];
	}
      }
    }
    
    if (rbc[1]) {
#ifdef _OPENMP
#pragma omp for schedule(static,1)
#endif
      for (ioff=1; ioff<maxoff;ioff++) {
#pragma ivdep
	for (i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++) {
	  v1[gec_v1[1]+ioff][i0]=v1[gec_v1[1]-ioff+1][i0];
	}
      }
    }
    
  }

}
