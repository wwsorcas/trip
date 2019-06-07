#include "cstd.h"

void acd_2d(float ** uc, 
	    float ** up, 
	    float ** csq, 
	    int * s, 
	    int * e, 
	    float ** c, 
	    int k,
	    int * lbc,
	    int * rbc,
	    float * lap) {
  
  int i0, i1, ioff;
  int s0=s[0];
  int e0=e[0];

  //  fprintf(stderr,"acd_2d\n");

  float c0 = c[0][0]+c[1][0];

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  //  begin parallel region
  {
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif

    for (i1=s[1]; i1<=e[1]; i1++) {
      // offset independent part of stencil
#pragma ivdep
      for (i0=s0; i0<=e0; i0++) {    
	up[i1][i0] = 2.0*uc[i1][i0] - up[i1][i0];
	lap[i0] = c0*uc[i1][i0];
      }
      // offset dependent part of stencil
      for (ioff=1;ioff<=k;ioff++) {
#pragma ivdep
	for (i0=s0; i0<=e0; i0++) {
	  lap[i0] +=
	    c[0][ioff]*(uc[i1][i0+ioff] + uc[i1][i0-ioff]) +
	    c[1][ioff]*(uc[i1+ioff][i0] + uc[i1-ioff][i0]);
	}
      }
#pragma ivdep      
      for (i0=s0; i0<=e0; i0++) {    
	up[i1][i0] += csq[i1][i0]*lap[i0];
      }
    }

    /* boundary conditions - note that uc[-1][i]=0 etc. */
    if (lbc[1]) {
      
#ifdef _OPENMP
#pragma omp for schedule(static, 1)
#endif
      
      for (ioff=0;ioff<k-1;ioff++) {
#pragma ivdep
	for (i0=s0;i0<=e0;i0++) {
	  up[s[1]-ioff-2][i0]=-up[s[1]+ioff][i0];
	}
      }
    }
    if (rbc[1]) {

#ifdef _OPENMP
#pragma omp for schedule(static, 1)
#endif
    
      for (ioff=0;ioff<k-1;ioff++) {
#pragma ivdep
	for (i0=s0;i0<=e0;i0++) {
	  up[e[1]+2+ioff][i0]=-up[e[1]-ioff][i0];
	}
      }
    }

    if (lbc[0]) {

#ifdef _OPENMP
#pragma omp for schedule(static, 1)
#endif
  
      for (i1=s[1];i1<=e[1];i1++) {
	for (ioff=0;ioff<k-1;ioff++) {
	  up[i1][s[0]-ioff-2]=-up[i1][s[0]+ioff];
	}
      }
    }

    if (rbc[0]) {

#ifdef _OPENMP
#pragma omp for schedule(static, 1)
#endif  

      for (i1=s[1];i1<=e[1];i1++) {
	for (ioff=0;ioff<k-1;ioff++) {
	  up[i1][e[0]+ioff+2]=-up[i1][e[0]-ioff];
	}
      }
    }
    
  } // end parallel region
    
}


