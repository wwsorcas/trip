#include "cstd.h"

//#undef REORDER
#define REORDER

#define my_min(a, b) ((a) < (b) ? (a) : (b))
#define my_max(a, b) ((a) > (b) ? (a) : (b))


/* staggered grid derivative in dimension idx
 */
void sggrad(int * gsx, int * gex,
	    int * gsy0, int * gey0,
	    int * csy0, int * cey0,
	    int * gsy1, int * gey1,
	    int * csy1, int * cey1,
	    int m, int dim, 
	    float * x,
	    float * a0,
	    float * a1,
	    float * y0,
	    float * y1)
{

  // compute offsets, lengths
  int * oxy0 = (int *)malloc(dim*sizeof(int));
  int * oxy1 = (int *)malloc(dim*sizeof(int));
  int * nx = (int *)malloc(dim*sizeof(int));
  int * ny0 = (int *)malloc(dim*sizeof(int));
  int * ny1 = (int *)malloc(dim*sizeof(int));
  for (int i=0;i<dim; i++) {
    oxy0[i]=gsy0[i]-gsx[i];
    oxy1[i]=gsy1[i]-gsx[i];
    nx[i]=gex[i]-gsx[i]+1;
    ny0[i]=gey0[i]-gsy0[i]+1;
    ny1[i]=gey1[i]-gsy1[i]+1;
  }

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
#ifdef _OPENMP    
#pragma omp for schedule(static, OMP_CHUNK_SIZE)    
#endif
    for (int i1=csy0[1]-gsy0[1]; i1<cey0[1]-gsy0[1]+1; i1++) {
      for (int j=0; j<m; j++) {
#pragma ivdep
#pragma GCC ivdep
	for (int i0=csy0[0]-gsy0[0]; i0< cey0[0]-gsy0[0]+1; i0++) {
	  y0[i0 + i1*ny0[0]] +=
	    a0[j]*(x[i0+oxy0[0]+1+j + (i1+oxy0[1])*nx[0]] - x[i0+oxy0[0] -j + (i1+oxy0[1])*nx[0]]);	     
	}
      }
    }
#ifdef _OPENMP    
#pragma omp for schedule(static, OMP_CHUNK_SIZE)    
#endif
    for (int i1=csy1[1]-gsy1[1]; i1<cey1[1]-gsy1[1]+1; i1++) {
      for (int j=0; j<m; j++) {
#pragma ivdep
#pragma GCC ivdep
	for (int i0=csy1[0]-gsy1[0]; i0< cey1[0]-gsy1[0]+1; i0++) {
	  y1[i0 + i1*ny1[0]] +=
	    a1[j]*(x[i0+oxy1[0] + (i1+oxy1[1]+1+j)*nx[0]] - x[i0+oxy1[0] + (i1+oxy1[1]-j)*nx[0]]);	     
	}
      }
    } 
  }
  free(oxy0);
  free(oxy1);
  free(nx);
  free(ny0);
  free(ny1);

}

int sgdiv(int * gsx, int * gex,
	  int * gsy0, int * gey0,
	  int * csy0, int * cey0,
	  int * gsy1, int * gey1,
	  int * csy1, int * cey1,
	  int m, int dim, 
	  float * x,
	  float * a0,
	  float * a1,
	  float * y0,
	  float * y1)
{

  // only one error return, if idx out-of-range
  int ierr=0;
  
  // compute offsets, lengths
  int * oxy0 = (int *)malloc(dim*sizeof(int));
  int * oxy1 = (int *)malloc(dim*sizeof(int));
  int * nx = (int *)malloc(dim*sizeof(int));
  int * ny0 = (int *)malloc(dim*sizeof(int));
  int * ny1 = (int *)malloc(dim*sizeof(int));
  for (int i=0;i<dim; i++) {
    oxy0[i]=gsy0[i]-gsx[i];
    oxy1[i]=gsy1[i]-gsx[i];
    nx[i]=gex[i]-gsx[i]+1;
    ny0[i]=gey0[i]-gsy0[i]+1;
    ny1[i]=gey1[i]-gsy1[i]+1;
  }

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
#ifdef _OPENMP    
#pragma omp for schedule(static, OMP_CHUNK_SIZE)    
#endif
      
    for (int i1=csy0[1]-gsy0[1]; i1<cey0[1]-gsy0[1]+1; i1++) {
      
#ifndef REORDER // naive adjoint loop

      for (int j=0; j<m; j++) {
#pragma ivdep
#pragma GCC ivdep
	for (int i0=csy0[0]-gsy0[0]; i0< cey0[0]-gsy0[0]+1; i0++) {
	  x[i0+oxy0[0]+j+1 + (i1+oxy0[1])*nx[0]] += a0[j]*y0[i0 + i1*ny0[0]];
	  x[i0+oxy0[0]-j   + (i1+oxy0[1])*nx[0]] -= a0[j]*y0[i0 + i1*ny0[0]];
	}
      }
      
      /* basic adjoint: same loop and
	 x[i0+oxy[0]+1+j + (i1+oxy[1])*nx[0]] += a[j]*y[i0 + i1*ny[0]]
	 x[i0+oxy[0] -j + (i1+oxy[1])*nx[0]]  -= a[j]*y[i0 + i1*ny[0]]
	 i0+j+1=k0 in first eqn i0=k0-j-1
	 limits: k0 >= csy[0]-gsy[0]+j+1 k0+oxy[0]>=0,
	 k0 <  cey[0]-gsy[0]+1+J+1 k0+oxy[0]<nx[0]
	 x[k0+oxy[0] + (i1+oxy[1])*nx[0]] += a[j]*y[k0-j-1 + i1*ny[0]]
	 i0-j=k0 in 2nd eqn i0=k0+j
	 limits: k0 >= csy[0]-gsy[0]-j k0+oxy[0]>=0
	 k0 <  cey[0]-gsy[0]+1-j k0+oxy[0]<nx[0]
	 x[k0+oxy[0] + (i1+oxy[1])*nx[0]] -= a[j]*y[k0+j  + i1*ny[0]]
      */
	
#else // reordered loop
      
      for (int j=0; j<m; j++) {
#pragma ivdep
#pragma GCC ivdep
	for (int k0=my_max(csy0[0]-gsy0[0]+j+1,-oxy0[0]);
	     k0<my_min(cey0[0]-gsy0[0]+1+j+1,nx[0]-oxy0[0]); k0++) {
	  x[k0+oxy0[0] + (i1+oxy0[1])*nx[0]] +=	a0[j]*y0[k0-j-1 + i1*ny0[0]];
	}
#pragma ivdep
#pragma GCC ivdep      
	for (int k0=my_max(csy0[0]-gsy0[0]-j,-oxy0[0]);
	     k0<my_min(cey0[0]-gsy0[0]+1-j,nx[0]-oxy0[0]); k0++) {
	  x[k0+oxy0[0] + (i1+oxy0[1])*nx[0]] -=
	    a0[j]*y0[k0+j  + i1*ny0[0]];
	}		   
      }
	
#endif
    }

#ifndef REORDER // naive adjoint loop
      
#ifdef _OPENMP    
#pragma omp for schedule(static, OMP_CHUNK_SIZE)    
#endif
    for (int i1=csy1[1]-gsy1[1]; i1<cey1[1]-gsy1[1]+1; i1++) {

      for (int j=0; j<m; j++) {
#pragma ivdep
#pragma GCC ivdep
	for (int i0=csy1[0]-gsy1[0]; i0< cey1[0]-gsy1[0]+1; i0++) {
	  x[i0+oxy1[0] + (i1+oxy1[1]+1+j)*nx[0]] += a1[j]*y1[i0 + i1*ny1[0]]; 
	  x[i0+oxy1[0] + (i1+oxy1[1]  -j)*nx[0]] -= a1[j]*y1[i0 + i1*ny1[0]]; 	     
	}
      }

    }
  }
  
#else
  }

  for (int j=0; j<m; j++) {

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
    {
    
#ifdef _OPENMP    
#pragma omp for schedule(static, OMP_CHUNK_SIZE)    
#endif
    //    for (int i1=csy1[1]-gsy1[1]; i1<cey1[1]-gsy1[1]+1; i1++) {
      for (int k1=my_max(csy1[1]-gsy1[1]+j+1,-oxy1[1]);
	   k1<my_min(cey1[1]-gsy1[1]+1+j+1,nx[1]-oxy1[1]); k1++) {
#pragma ivdep
#pragma GCC ivdep
	for (int i0=csy1[0]-gsy1[0]; i0< cey1[0]-gsy1[0]+1; i0++) {
	  x[i0+oxy1[0] + (k1+oxy1[1])*nx[0]] += a1[j]*y1[i0 + (k1-j-1)*ny1[0]]; 
	}
      }
#ifdef _OPENMP    
#pragma omp for schedule(static, OMP_CHUNK_SIZE)    
#endif      
      for (int k1=my_max(csy1[1]-gsy1[1]-j,-oxy1[1]);
	     k1<my_min(cey1[1]-gsy1[1]+1-j,nx[1]-oxy1[1]); k1++) {
#pragma ivdep
#pragma GCC ivdep      
	for (int i0=csy1[0]-gsy1[0]; i0< cey1[0]-gsy1[0]+1; i0++) {
	  x[i0+oxy0[0] + (k1+oxy1[1])*nx[0]] -= a1[j]*y1[i0+(k1+j)*ny1[0]];
	}		   
      }
    }
  }
#endif
	
  free(oxy0);
  free(oxy1);  
  free(nx);
  free(ny0);
  free(ny1);
  return ierr; 
}
