#include "parser.h"
#include "except.hh"

//--- omp ---
#define OMP_CHUNK_SIZE 32
//#define OMP_NUM_THREADS 8

//#undef REORDER
#define REORDER

#define N0 1000
#define N1 1000
#define NT 10000
#define M 8

extern "C" int sggrad(int * gsx, int * gex,
		      int * gsy0, int * gey0,
		      int * csy0, int * cey0,
		      int * gsy1, int * gey1,
		      int * csy1, int * cey1,
		      int m, int dim, 
		      float * x,
		      float * a0,
		      float * a1,
		      float * y0,
		      float * y1);

extern "C" int sgdiv(int * gsx, int * gex,
		      int * gsy0, int * gey0,
		      int * csy0, int * cey0,
		      int * gsy1, int * gey1,
		      int * csy1, int * cey1,
		      int m, int dim, 
		      float * x,
		      float * a0,
		      float * a1,
		      float * y0,
		      float * y1);

/* dot product */
float dot(int n,
	  float * x,
	  float * y) {
  float p = 0.0f;
#pragma ivdep
#pragma GCC ivdep            
  for (int i=0; i<n; i++) {
    p += x[i]*y[i];
  }
  return p;
}

/* dimension check for staggered grid partial derivatives

   for active index (idx),

   y[i]=sum_j=0,m-1 a[j](x[i+j+1]-x[i-j])

   max i+j+1 = i+m, min = i+1
   max i-j   = i,   min = i-m+1

   so in active dim, x index range must go from <= gsy-m+1 to gey+m

   if y range is 0 to n-1 in active range, then x range should be -m+1 to n+m-1, and x length should be
   n+m-1-(-m+1)+1 = n+2m-1

   in other dimensions index ranges should be same, or y's should be subset of x's

   offset test for start: gsy-gsx >= m-1

   relation between input, output indices: in each dimension, oxy = iy-ix
   so in loop over iy, set ix=iy-oxy

   gs, ge = global start, end indices
   os, oe = desired output start indices - must define subrange of global
 */
void gradcheck(int * gsx, int * gex,
	       int * gsy, int * gey,
	       int * csy, int * cey,
	       int m, int idx, int dim) 
{
  // compute offsets, lengths
  IPNT oxy;
  for (int i=0;i<dim; i++) {
    oxy[i]=gsy[i]-gsx[i];
  }
  if (oxy[idx]<m-1) {
    RVL::RVLException e;
    e<<"Error: oxy["<<idx<<"] < m="<<m<<"\n";
    throw e;
  }
  for (int i=0;i<dim;i++) {
    if ((csy[i]<gsy[i]) || (cey[i]>gey[i])) {
      RVL::RVLException e;
      e<<"Error: dim "<<i<<" output range [cs,ce]=["<<csy[i]<<","<<cey[i]<<"] not subrange\n";
      e<<"  of global range [gs,ge]=["<<gsy[i]<<","<<gey[i]<<"]\n";
      throw e;
    }      
    if (i==idx) {
      if ((gsx[i]>gsy[i]-m+1) || (gex[i]<gey[i]+m)) {
        RVL::RVLException e;
        e<<"Error: in active dim = "<<idx<<"gsx="<<gsx[i]<<" gsy[i]-m+1="<<gsy[i]-m+1<<"\n";
        e<<"Error: in active dim = "<<idx<<"gex="<<gex[i]<<" gey[i]+m  ="<<gey[i]+m<<"\n";
        throw e;
      }
    }
    else {
      if ((gsx[i]>gsy[i] || gex[i]<gey[i])) {
        RVL::RVLException e;
        e<<"Error: in inactive dim = "<<idx<<"gsx="<<gsx[i]<<" gsy[i]="<<gsy[i]<<"\n";
        e<<"Error: in inactive dim = "<<idx<<"gex="<<gex[i]<<" gey[i]  ="<<gey[i]<<"\n";
        throw e;
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

  int n0 = N0;
  int n1 = N1;
  int m  = M;

  float * a0 = (float *)malloc(m*sizeof(float));
  float * a1 = (float *)malloc(m*sizeof(float));  
  
  float * y00 = (float *)malloc(n0*n1*sizeof(float));
  float * y01 = (float *)malloc(n0*n1*sizeof(float));  
  float * y10 = (float *)malloc(n0*n1*sizeof(float));
  float * y11 = (float *)malloc(n0*n1*sizeof(float));

  IPNT gsy0; IPNT gey0; IPNT csy0; IPNT cey0;
  gsy0[0]=0; gey0[0]=n0-1;
  gsy0[1]=0; gey0[1]=n1-1;
  csy0[0]=gsy0[0];
  cey0[0]=gey0[0];
  csy0[1]=gsy0[1];
  cey0[1]=gey0[1];

  IPNT gsy1; IPNT gey1; IPNT csy1; IPNT cey1;
  gsy1[0]=0; gey1[0]=n0-1;
  gsy1[1]=0; gey1[1]=n1-1;
  csy1[0]=gsy1[0];
  cey1[0]=gey1[0];
  csy1[1]=gsy1[1];
  cey1[1]=gey1[1];

  // fprintf(stderr,"allocate x\n");
  float * x0 = (float *)malloc((n0+2*m-1)*(n1+2*m-1)*sizeof(float));
  float * x1 = (float *)malloc((n0+2*m-1)*(n1+2*m-1)*sizeof(float));
  IPNT gsx; IPNT gex;
  gsx[0]=-m+1; gex[0]=n0+m-1;
  gsx[1]=-m+1; gex[1]=n1+m-1;
  
  srand(getpid());
  //srand(19490615);
  for (int j=0; j<m; j++) { a0[j] = (-0.5+ rand()/(RAND_MAX+1.0)); }
  for (int j=0; j<m; j++) { a1[j] = (-0.5+ rand()/(RAND_MAX+1.0)); }

  // fprintf(stderr,"initialize x\n");
  for (int i=1; i<(n0+2*m-1)*(n1+2*m-1); i++) {
    x0[i] = (-0.5+ rand()/(RAND_MAX+1.0));
    x1[i] = 0.0f;
  }

  for (int i=0;i<n0*n1; i++) {
    y00[i]=0.0f;
    y01[i]=0.0f;    
    y10[i]=(-0.5+ rand()/(RAND_MAX+1.0));
    y11[i]=(-0.5+ rand()/(RAND_MAX+1.0));
    //    y11[i]=0.0f;
  }

  int idx=0;
  int dim=2;
  
  try {
    idx=0;
    gradcheck(gsx, gex, gsy0, gey0, csy0, cey0, m, idx, dim);
    idx=1;
    gradcheck(gsx, gex, gsy1, gey1, csy1, cey1, m, idx, dim);      

    // first test: dot product of 0 axis deriv, adjoint
    cerr<<"dot test\n";
    sggrad(gsx, gex, gsy0, gey0, csy0, cey0, gsy1, gey1, csy1, cey1, m, dim, x0, a0, a1, y00, y01);
    sgdiv(gsx, gex, gsy0, gey0, csy0, cey0, gsy1, gey1, csy1, cey1, m, dim, x1, a0, a1, y10, y11);
    cerr<<"y1^T grad x0 = "<<dot(n1*n0,y10,y00)+dot(n1*n0,y11,y01)<<endl;
    cerr<<"x0^T div y1  = "<<dot((n1+2*m-1)*(n0+2*m-1),x1,x0)<<endl;

    cerr<<"grad loop\n";
    for (int k=0;k<NT;k++) {
      x0[0]=(-0.5+ rand()/(RAND_MAX+1.0));
      sggrad(gsx, gex, gsy0, gey0, csy0, cey0, gsy1, gey1, csy1, cey1, m, dim, x0, a0, a1, y00, y01);
    }

    cerr<<"div loop\n";
    for (int k=0;k<NT;k++) {    
      y10[0]=(-0.5+ rand()/(RAND_MAX+1.0));      
      sgdiv(gsx, gex, gsy0, gey0, csy0, cey0, gsy1, gey1, csy1, cey1, m, dim, x1, a0, a1, y10, y11);    
    }
    
  }
  catch (RVL::RVLException e) {
    e.write(cerr);
    exit(1);
  }
}
	 
	 
