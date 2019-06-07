#include "cstd.h"

void acde_2d_4_b(float **uc, float **ucb, float **up, float **upb,
		 float **csq, float ***csqb,
		 int *s, int *e, int *z,
		 float c0, float *c1, float *c2, int *lbc, int *rbc,
		 float *dx, float *lap) {
  
  int i0, i1, i2;
  int s0 = s[0];
  int e0 = e[0];
  
  /* boundary conditions - note that uc[-1][i]=0 etc. */
  if (rbc[0]) {
    for (i1 = e[1]; i1 >= s[1]; i1--) {
      upb[i1][e[0]]  -= upb[i1][e[0]+2];
      upb[i1][e[0]+2] = 0.0;
    }
  }
  if (lbc[0]) {
    for (i1 = e[1]; i1 >= s[1]; i1--) {
      upb[i1][s[0]]  -= upb[i1][s[0]-2];
      upb[i1][s[0]-2] = 0.0;
    }
  }
  if (rbc[1]) {
#pragma ivdep
    for (i0 = e0; i0 >= s0; i0--) {
      upb[e[1]][i0]  -= upb[e[1]+2][i0];
      upb[e[1]+2][i0] = 0.0;
    }
  }
  if (lbc[1]) {
#pragma ivdep
    for (i0 = e0; i0 >= s0; i0--) {
      upb[s[1]][i0]  -= upb[s[1]-2][i0];
      upb[s[1]-2][i0] = 0.0;
    }
  }
  
  /* calculate the Laplacian operator as a function of x */
  for (i1 = e[1]; i1 >= s[1]; i1--) {
#pragma ivdep
    for (i0 = e0; i0 >= s0; i0--) {
      lap[i0] = c0*uc[i1][i0] + 
	        c1[0]*(uc[i1][i0+1] + uc[i1][i0-1]) + 
	        c1[1]*(uc[i1+1][i0] + uc[i1-1][i0]) + 
	        c2[0]*(uc[i1][i0+2] + uc[i1][i0-2]) + 
	        c2[1]*(uc[i1+2][i0] + uc[i1-2][i0]);
    } /* end loop over x */
    
    /* calculate the extended csqb operator */
    for (i2 = s[2]-z[2]; i2 <= e[2]-z[2]; i2++) {
      int locs0 = (s[0] > s[0] + 2*i2 ? s[0] : s[0] + 2*i2);
      int loce0 = (e[0] < e[0] + 2*i2 ? e[0] : e[0] + 2*i2);
      int loci2 = i2 + z[2];
      //      memset(csqb[loci2][i1], 0.f, (e0-s0+1)*sizeof(float));
#pragma ivdep          
      for (i0 = loce0; i0 >= locs0; i0--) {
	csqb[loci2][i1][i0-i2] += lap[i0-2*i2] * upb[i1][i0];
      }  /* end loop over x */
    } /* end loop over h */
    
    /* extrapolate the previous time step */
#pragma ivdep
    for (i0 = e0; i0 >= s0; i0--) {            
      float csqupb = csq[i1][i0] * upb[i1][i0];
      ucb[i1][i0+2] += c2[0] * csqupb;
      ucb[i1][i0+1] += c1[0] * csqupb;
      ucb[i1+2][i0] += c2[1] * csqupb;
      ucb[i1+1][i0] += c1[1] * csqupb;
      ucb[i1][i0]   += c0 * csqupb + 2.0*upb[i1][i0];
      ucb[i1][i0-1] += c1[0] * csqupb;
      ucb[i1][i0-2] += c2[0] * csqupb;
      ucb[i1-1][i0] += c1[1] * csqupb;
      ucb[i1-2][i0] += c2[1] * csqupb;

      upb[i1][i0] = -upb[i1][i0];
    } /* end loop over x */
  } /* end loop over z */

}
