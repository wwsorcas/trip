#include "parser.h"
#include "except.hh"

#define MY_CHUNK_SIZE 16

/* y = x conv a */
void scale(int n0, int n1,
	 float * x,
	 float * y) {
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
#ifdef _OPENMP
#pragma omp for schedule(static, MY_CHUNK_SIZE)
#endif  
    for (int i1=0; i1<n1; i1++) 
#pragma ivdep
#pragma GCC ivdep
      for (int i0=0; i0<n0; i0++)
	x[i1*n0+i0]*=y[i0];
  }
}
      
int main(int argc, char ** argv) {

  PARARRAY * par = ps_new();
  ps_createargs(par,argc-1,&(argv[1]));
  int num_threads = 1;
  ps_flint(*par,"num_threads",&num_threads);
  
#ifdef _OPENMP
  omp_set_num_threads(num_threads);
  printf("Number of OMP threads in use = %d\n",omp_get_max_threads());
#endif

  int n0 = 2000;
  int n1 = 300;
  int nt = 30000;
  
  float * x = (float *)malloc(n0*n1*sizeof(float));
  float * y = (float *)malloc(n0*sizeof(float));
  
  srand(getpid());
  for (int i=0; i<n0*n1; i++) { x[i] = (-0.5+ rand()/(RAND_MAX+1.0)); }
  for (int i=0; i<n0; i++) { y[i] = (-0.5+ rand()/(RAND_MAX+1.0)); }

  for (int k=0; k<nt; k++) {
    x[0]= (-0.5+ rand()/(RAND_MAX+1.0));
    scale(n0,n1,x,y);
  }
}
	 
	 
