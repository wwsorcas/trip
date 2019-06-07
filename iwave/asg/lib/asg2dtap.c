#include "cstd.h"

#define MY_CHUNK_SIZE 16

#define my_min(a, b) ((a) < (b) ? (a) : (b))
#define my_max(a, b) ((a) > (b) ? (a) : (b))

//#undef IHFIRST
#define IHFIRST

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

  // offset loop: accumulate
  //	    dh * bulk[ih+ihmax][i1-ih][i0] * sdiv[i1-gsc_p[1]-2*ih][i0-gsc_p[0]];  
  // two ways to handle ih vs i1 loops: if you go ih first, then
#ifdef IHFIRST
  int * sh = (int *)malloc((2*ihmax+1)*sizeof(int));
  int * eh = (int *)malloc((2*ihmax+1)*sizeof(int));  
  for (int ih=-ihmax; ih<=ihmax; ih++) {
    //      sh = ((ih <= 0) ?  gsc_p[1]-ih : gsc_p[1]);
    //      eh = ((ih >= 0) ?  gec_p[1]-ih : gec_p[1]);
    if (ih <=0) {
      sh[ih+ihmax] = gsc_p[1]; eh[ih+ihmax] = gec_p[1] + 2*ih;
    }
    else {
      sh[ih+ihmax] = gsc_p[1] + 2*ih; eh[ih+ihmax] = gec_p[1];
    }
  }
#else
  // if you go i1 first, then ih must satisfy
  // -ihmax<=ih<=ihmax
  // gsc_p[1] <= i1-ih <= gec_p[1]
  // 0 <= i1-gsc_p[1]-2*ih <= gec_p[1]-gsc_p[1]+1
  // from the second we get i1-gec_p[1] <= ih <= i1-gsc_p[1]
  // from the third we get i1-gec_p[1]-1 <= 2*ih <= i1-gsc_p[1]
  // since i1-gsc_p[1] >=0 and i1-gec_p[1]-1 < 0, only the second condition pertains
  // so
  int * sh = (int *)malloc((gec_p[1]-gsc_p[1]+1)*sizeof(float));
  int * eh = (int *)malloc((gec_p[1]-gsc_p[1]+1)*sizeof(float));
  for (int i1=gsc_p[1]; i1<=gec_p[1]; i1++) {
    //    sh[i1]=my_max(-ihmax,(i1-gec_p[1]-1)/2);
    sh[i1-gsc_p[1]]=my_max(-ihmax,(i1-gec_p[1])/2);
    eh[i1-gsc_p[1]]=my_min(ihmax,(i1-gsc_p[1])/2);
  }
  // check: if ihmax=0, sh=eh=0
  // if i1 is in middle of range, then i1-gec_p << 0, i1-gsc_p >> 0, so sh=-ihmax, eh=ihmax
#endif

  // begin parallel region
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
    // fprintf(stderr,"1\n");
#ifdef _OPENMP
#pragma omp for schedule(static, MY_CHUNK_SIZE)
#endif
    for (int i1=gsc_p[1]; i1<=gec_p[1]; i1++) {
      memset((char *)(sdiv[i1-gsc_p[1]]),  0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
      memset((char *)(sdivb[i1-gsc_p[1]]), 0, (gec_p[0]-gsc_p[0]+1)*sizeof(float));
    }
    // fprintf(stderr,"2\n");
#ifdef _OPENMP
    //#pragma omp for schedule(static, MY_CHUNK_SIZE) collapse(2)
#pragma omp for schedule(static, MY_CHUNK_SIZE)
#endif
    for (int i1=gsc_p[1]; i1<=gec_p[1]; i1++) {    
      for (int ioff = 0; ioff<maxoff; ioff++) {
#pragma ivdep
#pragma GCC ivdep
	for (int i0=gsc_p[0]; i0 <= gec_p[0]; i0++) {
	  sdiv[i1-gsc_p[1]][i0-gsc_p[0]] +=
	    c[0][ioff]*(v0[i1][i0+ioff]-v0[i1][i0-ioff-1]) +
	    c[1][ioff]*(v1[i1+ioff][i0]-v1[i1-ioff-1][i0]);
	  //	  fprintf(stream,"i0=%d i1=%d ioff=%d v0=%g v1=%g\n",
	  //		  i0,i1,ioff,v0[i1][i0+ioff],v1[i1+ioff][i0]);
	}
      }
    }

#ifdef IHFIRST

  } // end par region 
  

  for (int ih=-ihmax; ih<=ihmax; ih++) {

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif

    {     // begin parallel region
    
#ifdef _OPENMP
#pragma omp for schedule(static, MY_CHUNK_SIZE)
#endif
      for (int i1=sh[ih+ihmax]; i1<=eh[ih+ihmax]; i1++) {
#pragma ivdep
#pragma GCC ivdep
	for (int i0=gsc_p[0]; i0 <= gec_p[0]; i0++) {      
	  sdivb[i1-gsc_p[1]][i0-gsc_p[0]] +=
	    dh * bulk[ih+ihmax][i1-ih][i0] * sdiv[i1-gsc_p[1]-2*ih][i0-gsc_p[0]];
	  //	  fprintf(stream,"i0=%d i1=%d ih=%d bulk[ih+ihmax][i1-ih][i0]=%g\n",
	  //		  i0,i1,bulk[ih+ihmax][i1-ih][i0]);
	} 
      }
    } // end par region

  } // end ih loop

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif

  {     // begin parallel region
  
#else
    
    // fprintf(stderr,"3\n");
#ifdef _OPENMP
#pragma omp for schedule(static, MY_CHUNK_SIZE)
#endif
    for (int i1=gsc_p[1]; i1<=gec_p[1]; i1++) {
      for (int ih=sh[i1-gsc_p[1]]; ih<=eh[i1-gsc_p[1]]; ih++) {
#pragma ivdep
#pragma GCC ivdep
	for (int i0=gsc_p[0]; i0 <= gec_p[0]; i0++) {      
	  sdivb[i1-gsc_p[1]][i0-gsc_p[0]] +=
	    dh * bulk[ih+ihmax][i1-ih][i0] * sdiv[i1-gsc_p[1]-2*ih][i0-gsc_p[0]];
	} 
      }
    }
#endif

    int gscp1=gsc_p[1];
    int gscp0=gsc_p[0];
    // fprintf(stderr,"4\n");
    //    float pmax=0.0;
#ifdef _OPENMP
#pragma omp for schedule(static, MY_CHUNK_SIZE)
#endif
    for (int i1=gsc_p[1]; i1<=gec_p[1]; i1++) {      
#pragma ivdep
#pragma GCC ivdep
      for (int i0=gsc_p[0]; i0 <= gec_p[0]; i0++) {
	p0[i1][i0] = (ep[0][i0]*p0[i1][i0] - sdivb[i1-gscp1][i0-gscp0])*epp[0][i0];
	p1[i1][i0] = (ep[1][i1]*p1[i1][i0] - sdivb[i1-gscp1][i0-gscp0])*epp[1][i1];
	//	pmax += fabs(p0[i1][i0]);
	//	fprintf(stream,"i0=%d i1=%d p0=%g p1=%g\n",i0,i1,p0[i1][i0],p1[i1][i0]);
      }
    }
    //    fprintf(stream,"pmax=%g\n",pmax);

  }

    // fprintf(stderr,"5-8\n");
  if (lbc[0]) {
#ifdef _OPENMP
#pragma omp for schedule(static,MY_CHUNK_SIZE)
#endif
    for (int i1=gsc_p[1];i1<=gec_p[1];i1++) {
      p0[i1][gsc_p[0]-1]=0.0f;
      for (int ioff=1;ioff<maxoff;ioff++) {
	p0[i1][gsc_p[0]-ioff-1]=-p0[i1][gsc_p[0]+ioff-1];
      }
    }
  }

  if (rbc[0]) {
#ifdef _OPENMP
#pragma omp for schedule(static,MY_CHUNK_SIZE)
#endif
      for (int i1=gsc_p[1];i1<=gec_p[1];i1++) {
	p0[i1][gec_p[0]+1]=0.0f;
	for (int ioff=1;ioff<maxoff;ioff++) {
	  p0[i1][gec_p[0]+ioff+1]=-p0[i1][gec_p[0]-ioff+1];
	}
      }
    }

  if (lbc[1]) {
#pragma ivdep
#pragma GCC ivdep
    for (int i0=gsc_p[0];i0<=gec_p[0];i0++) {
      p1[gsc_p[1]-1][i0]=0.0f;
    }
    //#ifdef _OPENMP
    //#pragma omp for schedule(static,1)
    //#endif
    for (int ioff=1;ioff<maxoff;ioff++) {
#pragma ivdep
#pragma GCC ivdep
      for (int i0=gsc_p[0];i0<=gec_p[0];i0++) {
	p1[gsc_p[1]-ioff-1][i0]=-p1[gsc_p[1]+ioff-1][i0];
      }
    }
  }

  if (rbc[1]) {
#pragma ivdep
#pragma GCC ivdep
    for (int i0=gsc_p[0];i0<=gec_p[0];i0++) {
      p1[gec_p[1]+1][i0]=0.0f;
    }
    //#ifdef _OPENMP
    //#pragma omp for schedule(static,1)
    //#endif
    for (int ioff=1;ioff<maxoff;ioff++) {
#pragma ivdep
#pragma GCC ivdep
      for (int i0=gsc_p[0];i0<=gec_p[0];i0++) {
	p1[gec_p[1]+ioff+1][i0]=-p1[gec_p[1]-ioff+1][i0];
      }
    }
  }
  
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
    
  // fprintf(stderr,"9\n");
#ifdef _OPENMP
#pragma omp for schedule(static, MY_CHUNK_SIZE)
#endif
    for (int i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++ ) {
      
      memset((char *)gradp0, 0, (gec_v0[0]-gsc_v0[0]+1)*sizeof(float));

      for (int ioff=0; ioff<maxoff; ioff++) {
#pragma ivdep	
#pragma GCC ivdep
	for (int i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ ) {
	  gradp0[i0-gsc_v0[0]] += c[0][ioff]*(p0[i1][i0+ioff+1]-p0[i1][i0-ioff]);
	}
      }

#pragma ivdep	
#pragma GCC ivdep
      for (int i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ ) {      
	v0[i1][i0] = evp[0][i0]*
	  (ev[0][i0]*v0[i1][i0] -
	   0.5f*(buoy[i1][i0]+buoy[i1][i0+1])*gradp0[i0-gsc_v0[0]]);
      }
    }

    // fprintf(stderr,"12\n");
#ifdef _OPENMP
#pragma omp for schedule(static, MY_CHUNK_SIZE)
#endif      
    for (int i1=gsc_v1[1]; i1 <= gec_v1[1]; i1++ ) {
      memset((char *)gradp1, 0, (gec_v1[0]-gsc_v1[0]+1)*sizeof(float));

      for (int ioff=0; ioff<maxoff; ioff++) {
#pragma ivdep
#pragma GCC ivdep
	for (int i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) {
	  gradp1[i0-gsc_v1[0]] +=  c[1][ioff]*(p1[i1+ioff+1][i0]-p1[i1-ioff][i0]);
	}
      }

    // fprintf(stderr,"14\n");
#pragma ivdep
#pragma GCC ivdep
      for (int i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) {      
	v1[i1][i0] = evp[1][i1]*
	  (ev[1][i1]* v1[i1][i0] - 0.5f*(buoy[i1][i0]+buoy[i1+1][i0])*gradp1[i0-gsc_v1[0]]);

      }
    }

  } // end parallel region

    // fprintf(stderr,"15-18\n");
  if (lbc[0]) {
#ifdef _OPENMP
#pragma omp for schedule(static,MY_CHUNK_SIZE)
#endif
    for (int i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++) {
      for (int ioff=1; ioff<maxoff;ioff++) {
	v0[i1][gsc_v0[0]-ioff]=v0[i1][gsc_v0[0]+ioff-1];
      }
    }
  }

  if (rbc[0]) {
#ifdef _OPENMP
#pragma omp for schedule(static,MY_CHUNK_SIZE)
#endif
    for (int i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++) {
      for (int ioff=1; ioff<maxoff;ioff++) {
	v0[i1][gec_v0[0]+ioff]=v0[i1][gec_v0[0]-ioff+1];
      }
    }
  }
    
  if (lbc[1]) {
    //#ifdef _OPENMP
    //#pragma omp for schedule(static,1)
    //#endif
    for (int ioff=1; ioff<maxoff;ioff++) {
#pragma ivdep
#pragma GCC ivdep
      for (int i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++) {
	v1[gsc_v1[1]-ioff][i0]=v1[gsc_v1[1]+ioff-1][i0];
      }
    }
  }
    
  if (rbc[1]) {
    //#ifdef _OPENMP
    //#pragma omp for schedule(static,1)
    //#endif
    for (int ioff=1; ioff<maxoff;ioff++) {
#pragma ivdep
#pragma GCC ivdep
      for (int i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++) {
	v1[gec_v1[1]+ioff][i0]=v1[gec_v1[1]-ioff+1][i0];
      }
    }
  }
  //  fprintf(stderr,"end fd\n");

}
