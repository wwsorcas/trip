#include "cstd.h"

void acde_2d_4(float **uc, 
	       float **up, 
	       float ***csq, 
	       int *s, 
	       int *e, 
	       int *z,
	       float c0, 
	       float *c1,
	       float *c2,
	       int *lbc,
	       int *rbc,
	       float *dx,
	       float *lap,
	       float *csqlap) {
  
  int i0, i1, i2;
  int s0 = s[0];
  int e0 = e[0];

  //  fprintf(stderr,"s0=%d s1=%d s2=%d e0=%d e1=%d e2=%d\n",s[0],s[1],s[2],e[0],e[1],e[2]);
  /* calculate the Laplacian operator as a function of x */
  for (i1 = s[1]; i1 <= e[1]; i1++) {
#pragma ivdep
    for (i0 = s0; i0 <= e0; i0++) {
      lap[i0] = c0*uc[i1][i0] +
                c1[0]*(uc[i1][i0+1] + uc[i1][i0-1]) +
                c1[1]*(uc[i1+1][i0] + uc[i1-1][i0]) +
                c2[0]*(uc[i1][i0+2] + uc[i1][i0-2]) +
                c2[1]*(uc[i1+2][i0] + uc[i1-2][i0]);
      csqlap[i0] = 0.0f;
    } /* end loop over x */
    
    //    fprintf(stderr,"s[2]=%d e[2]=%d\n",s[2],e[2]);
    /* apply the extended csq operator on the laplacian */
    
    for (i2 = s[2]-z[2]; i2 <= e[2]-z[2]; i2++) {
      int locs0=s[0]>s[0]+2*i2?s[0]:s[0]+2*i2;
      int loce0=e[0]<e[0]-2*i2?e[0]:e[0]-2*i2;
      int loci2=i2+z[2];
#pragma ivdep
      for (i0 = locs0; i0 <= loce0; i0++) {
	csqlap[i0] += csq[loci2][i1][i0-i2] * lap[i0-2*i2];
      } // end loop over x 
    } // end loop over h       

    /* extrapolate the next time step */
    /* NOTE: this version assumes csq has already been scaled by dx[2]*/
    for (i0=s0; i0<=e0; i0++) {
      up[i1][i0] = 2.0*uc[i1][i0] - up[i1][i0] + csqlap[i0];//dx[2]*csqlap[i0];
    } /* end loop over x */
  } /* end loop over z */
  
  /* boundary conditions - note that uc[-1][i]=0 etc. */
  if (lbc[1]) {
#pragma ivdep
    for (i0 = s0; i0 <= e0; i0++) {
      up[s[1]-2][i0] = -up[s[1]][i0];
    }
  }
  if (rbc[1]) {
#pragma ivdep
    for (i0 = s0; i0 <= e0; i0++) {
      up[e[1]+2][i0] = -up[e[1]][i0];
    }
  }
  if (lbc[0]) {
    for (i1 = s[1]; i1 <= e[1]; i1++) {
      up[i1][s[0]-2] = -up[i1][s[0]];
    }
  }
  if (rbc[0]) {
    for (i1 = s[1]; i1 <= e[1]; i1++) {
      up[i1][e[0]+2] = -up[i1][e[0]];
    }
  }
}

 

