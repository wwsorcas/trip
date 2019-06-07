#include "parser.h"

//--- omp ---
#define OMP_CHUNK_SIZE 32
//#define OMP_NUM_THREADS 8

//#undef REORDER
#define REORDER

#define N0 1000
#define N1 1000
#define NT 10000
#define M 8

/* y = x conv a */
void grad(int n00, int n01,
	  int n10, int n11,
	  int m, 
	  float ** x0,
	  float ** x1,
	  float * a0,
	  float * a1,
	  float ** y0,
	  float ** y1) {
  //  fprintf(stderr,"grad\n");
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
#ifdef _OPENMP    
#pragma omp for schedule(static, OMP_CHUNK_SIZE)    
#endif
    for (int i1=0; i1<n01; i1++) {
#pragma ivdep
#pragma GCC ivdep
      for (int i0=m; i0<n00-m; i0++) {
	y0[i1][i0] = a0[0]*(x0[i1][i0]-x0[i1][i0-1]);
      }
      for (int j=1; j<m; j++) {
#pragma ivdep
#pragma GCC ivdep
	for (int i0=m; i0<n00-m; i0++) {
	  y0[i1][i0] += a0[j]*(x0[i1][i0+j]-x0[i1][i0-j-1]);
	}
      }
    }
#ifdef _OPENMP    
#pragma omp for schedule(static, OMP_CHUNK_SIZE)    
#endif
    for (int i1=m; i1<n11-m; i1++) {
#pragma ivdep
#pragma GCC ivdep
      for (int i0=0; i0<n10; i0++) {
	y1[i1][i0] = a1[0]*(x1[i1][i0]-x1[i1-1][i0]);
      }
      for (int j=1; j<m; j++) {
#pragma ivdep
#pragma GCC ivdep
	for (int i0=0; i0<n10; i0++) {
	  y1[i1][i0] += a1[j]*(x1[i1+j][i0]-x1[i1-j-1][i0]);
	}
      }
    }    
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

  int n00 = N0+M;
  int n01 = N0;
  int n10 = N1;
  int n11 = N1+M;
  int m = M;

  float * a0 = (float *)malloc(m*sizeof(float));
  float * a1 = (float *)malloc(m*sizeof(float));  
  
  float * x0 = (float *)malloc(n00*n01*sizeof(float));
  float ** xx0 = (float **)malloc(n01*sizeof(float*));
  for (int i=0;i<n01;i++) xx0[i]=&(x0[i*n00]);

  float * x1 = (float *)malloc(n10*n11*sizeof(float));
  float ** xx1 = (float **)malloc(n11*sizeof(float*));
  for (int i=0;i<n11;i++) xx1[i]=&(x1[i*n10]);

  float * y0 = (float *)malloc(n00*n01*sizeof(float));
  float ** yy0 = (float **)malloc(n01*sizeof(float*));
  for (int i=0;i<n01;i++) yy0[i]=&(y0[i*n00]);

  float * y1 = (float *)malloc(n10*n11*sizeof(float));
  float ** yy1 = (float **)malloc(n11*sizeof(float*));
  for (int i=0;i<n11;i++) yy1[i]=&(y1[i*n10]);

  srand(getpid());
  for (int j=0; j<m; j++) { a0[j] = (-0.5+ rand()/(RAND_MAX+1.0)); }
  for (int j=0; j<m; j++) { a1[j] = (-0.5+ rand()/(RAND_MAX+1.0)); }

  for (int i=0; i<n00*n01; i++) {
    x0[i] = (-0.5+ rand()/(RAND_MAX+1.0));
    y0[i] = 0.0f;
  }
  for (int i=0; i<n10*n11; i++) {
    x1[i] = (-0.5+ rand()/(RAND_MAX+1.0));
    y1[i] = 0.0f;
  }

  for (int k=0; k<NT; k++) {
    x0[0] = (-0.5+ rand()/(RAND_MAX+1.0));
    x1[0] = (-0.5+ rand()/(RAND_MAX+1.0));     
    grad(n00,n01,n10,n11,m,xx0,xx1,a0,a1,yy0,yy1);
  }
}
	 
	 
