#include "cstd.h"

void acde_2d_4_d_d(float **uc,  float **ucd0,  float **ucd,  float **ucdd,
		   float **up,  float **upd0,  float **upd,  float **updd,
		   float **csq, float **csqd0, float ***csqd, 
		   int *s, int *e, int *z,
		   float c0, float *c1, float *c2, int *lbc, int *rbc,
		   float *dx, float *lap, float *lapd0, float *csqdlap, float *csqdlapd0) {

  int i0, i1, i2;
  int s0 = s[0];
  int e0 = e[0];
  float lapd, lapdd;
  
  /* calculate the Laplacian operator as a function of x */
  for (i1 = s[1]; i1 <= e[1]; i1++) {
#pragma ivdep
    for (i0 = s0; i0 <= e0; i0++) {
      lapd0[i0] = c0*ucd0[i1][i0] + 
	              c1[0]*(ucd0[i1][i0+1] + ucd0[i1][i0-1]) + 
	              c1[1]*(ucd0[i1+1][i0] + ucd0[i1-1][i0]) + 
    		      c2[0]*(ucd0[i1][i0+2] + ucd0[i1][i0-2]) + 
    		      c2[1]*(ucd0[i1+2][i0] + ucd0[i1-2][i0]);
      lap[i0] =   c0*uc[i1][i0] + 
    		      c1[0]*(uc[i1][i0+1] + uc[i1][i0-1]) + 
    		      c1[1]*(uc[i1+1][i0] + uc[i1-1][i0]) + 
    		      c2[0]*(uc[i1][i0+2] + uc[i1][i0-2]) + 
    		      c2[1]*(uc[i1+2][i0] + uc[i1-2][i0]);
      csqdlap[i0]   = 0.f;
      csqdlapd0[i0] = 0.f;	      
    } /* end loop over x */
    
    /* apply the extended csqd operator on the laplacians */
    for (i2 = s[2]-z[2]; i2 <= e[2]-z[2]; i2++) {
      int locs0 = (s[0] > s[0] + 2*i2 ? s[0] : s[0] + 2*i2);
      int loce0 = (e[0] < e[0] - 2*i2 ? e[0] : e[0] - 2*i2);
      int loci2 = i2 + z[2];
#pragma ivdep
      for (i0 = locs0; i0 <= loce0; i0++) {
	csqdlapd0[i0] += csqd[loci2][i1][i0-i2] * lapd0[i0-2*i2];
	csqdlap[i0]   += csqd[loci2][i1][i0-i2] * lap[i0-2*i2];
      } /* end loop over x */
    } /* end loop over h */
    
#pragma ivdep
    for (i0 = s0; i0 <= e0; i0++) {
      /* extrapolate the next time step */
      lapdd = c0*ucdd[i1][i0] + 
	          c1[0]*(ucdd[i1][i0+1] + ucdd[i1][i0-1]) + 
	          c1[1]*(ucdd[i1+1][i0] + ucdd[i1-1][i0]) + 
	          c2[0]*(ucdd[i1][i0+2] + ucdd[i1][i0-2]) + 
    		  c2[1]*(ucdd[i1+2][i0] + ucdd[i1-2][i0]);
      lapd =  c0*ucd[i1][i0] + 
    		  c1[0]*(ucd[i1][i0+1] + ucd[i1][i0-1]) + 
    		  c1[1]*(ucd[i1+1][i0] + ucd[i1-1][i0]) + 
    		  c2[0]*(ucd[i1][i0+2] + ucd[i1][i0-2]) + 
    		  c2[1]*(ucd[i1+2][i0] + ucd[i1-2][i0]);
      updd[i1][i0] = 2.0*ucdd[i1][i0] - updd[i1][i0] + csqdlapd0[i0] + 
    		  csqd0[i1][i0]*lapd + csq[i1][i0]*lapdd;
      upd[i1][i0]  = 2.0*ucd[i1][i0]  - upd[i1][i0]  + csqdlap[i0] + 
    		  csq[i1][i0]*lapd;
      upd0[i1][i0] = 2.0*ucd0[i1][i0] - upd0[i1][i0] + csqd0[i1][i0]*lap[i0] + 
    		  csq[i1][i0]*lapd0[i0];
      up[i1][i0]   = 2.0*uc[i1][i0]   - up[i1][i0]   + csq[i1][i0]*lap[i0];
    } /* end loop over x */
  } /* end loop over z */

    /* boundary conditions - note that uc[-1][i]=0 etc. */
  if (lbc[1])
#pragma ivdep
    for (i0 = s0; i0 <= e0; i0++) {
      updd[s[1]-2][i0] = -updd[s[1]][i0];
      upd[s[1]-2][i0]  = -upd[s[1]][i0];
      upd0[s[1]-2][i0] = -upd0[s[1]][i0];
      up[s[1]-2][i0]   = -up[s[1]][i0];
    }
  if (rbc[1])
#pragma ivdep
    for (i0 = s0; i0 <= e0; i0++) {
      updd[e[1]+2][i0] = -updd[e[1]][i0];
      upd[e[1]+2][i0]  = -upd[e[1]][i0];
      upd0[e[1]+2][i0] = -upd0[e[1]][i0];
      up[e[1]+2][i0]   = -up[e[1]][i0];
    }
  if (lbc[0])
    for (i1 = s[1]; i1 <= e[1]; i1++) {
      updd[i1][s[0]-2] = -updd[i1][s[0]];
      upd[i1][s[0]-2]  = -upd[i1][s[0]];
      upd0[i1][s[0]-2] = -upd0[i1][s[0]];
      up[i1][s[0]-2]   = -up[i1][s[0]];
    }
  if (rbc[0])
    for (i1 = s[1]; i1 <= e[1]; i1++) {
      updd[i1][e[0]+2] = -updd[i1][e[0]];
      upd[i1][e[0]+2]  = -upd[i1][e[0]];
      upd0[i1][e[0]+2] = -upd0[i1][e[0]];
      up[i1][e[0]+2]   = -up[i1][e[0]];
    }
  
}
