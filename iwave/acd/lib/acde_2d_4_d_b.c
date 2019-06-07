#include "cstd.h"

void acde_2d_4_d_b(float **uc, float **ucb, float **ucd, float **ucdb, 
		   float **up, float **upb, float **upd, float **updb, 
		   float **csq, float **csqb, float ***csqd, 
		   int *s, int *e, int *z,
		   float c0, float *c1, float *c2, int *lbc, int *rbc,
		   float *dx, float *lapb) {
  
  int i0, i1, i2;
  int s0 = s[0];
  int e0 = e[0];
  
  /* boundary conditions - note that uc[-1][i]=0 etc. */
  if (rbc[0]) {
    for (i1 = e[1]; i1 >= s[1]; i1--) {
      upb[i1][e[0]]   -= upb[i1][e[0]+2];
      upb[i1][e[0]+2]  = 0.0;
      updb[i1][e[0]]  -= updb[i1][e[0]+2];
      updb[i1][e[0]+2] = 0.0;
    }
  }
  if (lbc[0]) {
    for (i1 = e[1]; i1 >= s[1]; i1--) {
      upb[i1][s[0]]   -= upb[i1][s[0]-2];
      upb[i1][s[0]-2]  = 0.0;
      updb[i1][s[0]]  -= updb[i1][s[0]-2];
      updb[i1][s[0]-2] = 0.0;
    }
  }
  if (rbc[1]) {
#pragma ivdep
    for (i0 = e0; i0 >= s0; i0--) {
      upb[e[1]][i0]   -= upb[e[1]+2][i0];
      upb[e[1]+2][i0]  = 0.0;
      updb[e[1]][i0]  -= updb[e[1]+2][i0];
      updb[e[1]+2][i0] = 0.0;
    }
  }
  if (lbc[1]) {
#pragma ivdep
    for (i0 = e0; i0 >= s0; i0--) {
      upb[s[1]][i0]   -= upb[s[1]-2][i0];
      upb[s[1]-2][i0]  = 0.0;
      updb[s[1]][i0]  -= updb[s[1]-2][i0];
      updb[s[1]-2][i0] = 0.0;
    }
  }

  /* calculate the Laplacian operators */
  for (i1 = e[1]; i1 >= s[1]; i1--) {
#pragma ivdep
    for (i0 = e0; i0 >= s0; i0--) {
      float lap  = c0*uc[i1][i0] + 
	           c1[0]*(uc[i1][i0+1] + uc[i1][i0-1]) + 
	           c1[1]*(uc[i1+1][i0] + uc[i1-1][i0]) + 
	           c2[0]*(uc[i1][i0+2] + uc[i1][i0-2]) + 
	           c2[1]*(uc[i1+2][i0] + uc[i1-2][i0]);
      float lapd = c0*ucd[i1][i0] + 
	           c1[0]*(ucd[i1][i0+1] + ucd[i1][i0-1]) + 
	           c1[1]*(ucd[i1+1][i0] + ucd[i1-1][i0]) + 
	           c2[0]*(ucd[i1][i0+2] + ucd[i1][i0-2]) + 
	           c2[1]*(ucd[i1+2][i0] + ucd[i1-2][i0]);
      csqb[i1][i0] += lap * upb[i1][i0] + lapd * updb[i1][i0];

      lapb[i0] = csq[i1][i0] * upb[i1][i0];
    } /* end loop over x */
    
    /* apply the extended csqd operator */
    for (i2 = s[2]-z[2]; i2 <= e[2]-z[2]; i2++) {
      int locs0 = (s[0] > s[0] + 2*i2 ? s[0] : s[0] + 2*i2);
      int loce0 = (e[0] < e[0] - 2*i2 ? e[0] : e[0] - 2*i2);
      int loci2 = i2 + z[2];
#pragma ivdep
      for (i0 = loce0; i0 >= locs0; i0--) {
	lapb[i0-2*i2] += csqd[loci2][i1][i0-i2] * updb[i1][i0];
      } /* end loop over x */
    } /* end loop over h */
    
    /* extrapolate the previous time step */
#pragma ivdep
    for (i0 = e0; i0 >= s0; i0--) {
      ucb[i1][i0+2] += c2[0] * lapb[i0];
      ucb[i1][i0+1] += c1[0] * lapb[i0];
      ucb[i1+2][i0] += c2[1] * lapb[i0];
      ucb[i1+1][i0] += c1[1] * lapb[i0];
      ucb[i1][i0]   += c0*lapb[i0] + 2.0*upb[i1][i0];
      ucb[i1][i0-1] += c1[0] * lapb[i0];
      ucb[i1][i0-2] += c2[0] * lapb[i0];
      ucb[i1-1][i0] += c1[1] * lapb[i0];
      ucb[i1-2][i0] += c2[1] * lapb[i0];
      
      float lapdb = csq[i1][i0] * updb[i1][i0];
      ucdb[i1][i0+2] += c2[0] * lapdb;
      ucdb[i1][i0+1] += c1[0] * lapdb;
      ucdb[i1+2][i0] += c2[1] * lapdb;
      ucdb[i1+1][i0] += c1[1] * lapdb;
      ucdb[i1][i0]   += c0 * lapdb + 2.0*updb[i1][i0];
      ucdb[i1][i0-1] += c1[0] * lapdb;
      ucdb[i1][i0-2] += c2[0] * lapdb;
      ucdb[i1-1][i0] += c1[1] * lapdb;
      ucdb[i1-2][i0] += c2[1] * lapdb;

      upb[i1][i0]  = -upb[i1][i0];
      updb[i1][i0] = -updb[i1][i0];
    } /* end loop over x */
  } /* end loop over z */
  
}
