#include "parser.h"
#include "except.hh"

//--- omp ---
#define OMP_CHUNK_SIZE 32
//#define OMP_NUM_THREADS 8

#undef DIVONLY
//#define DIVONLY

//#undef PUPDONLY
#define PUPDONLY

extern "C" void asgdiv(float *  restrict  v0,   float *  restrict  v1,
		       float *  restrict sdiv,
		       int * gsc_p,  int * gec_p, 
		       int * gsa_v0, int * gea_v0,
		       int * gsa_v1, int * gea_v1,
		       int maxoff, float ** restrict  c);

extern "C" void asgpupd(float *  restrict  bulk,
			float *  restrict  p0,   float *  restrict  p1,
			float ** restrict ep,    float ** restrict  epp,
			float *  restrict sdiv,
			int * gsa_b,  int * gea_b,
			int * gsa_p0, int * gea_p0,
			int * gsa_p1, int * gea_p1,
			int * gsc_p,  int * gec_p);

  extern "C" void asg2dtaptest(float *  restrict  bulk, float *  restrict  buoy,
			     float *  restrict  p0,   float *  restrict  p1,
			     float *  restrict  v0,   float *  restrict  v1,
			     float ** restrict ep,    float ** restrict  epp,
			     float ** restrict ev,    float ** restrict  evp,
			     float *  restrict sdiv,  float * restrict  sdivb,
			     float *  restrict gradp0,
			     float *  restrict gradp1,
			     int * gsa_b,  int * gea_b,
			     int * gsa_p0, int * gea_p0,
			     int * gsa_p1, int * gea_p1,
			     int * gsc_p,  int * gec_p, 
			     int * gsa_v0, int * gea_v0,
			     int * gsc_v0, int * gec_v0,
			     int * gsa_v1, int * gea_v1,
			     int * gsc_v1, int * gec_v1,		  
			     int * lbc, int * rbc,
			     float dh, int ihmax,
			     int maxoff, float ** restrict  c);
			     //		  int maxoff, float ** restrict  c,
			     //		  FILE * stream);

int main(int argc, char ** argv) {

  PARARRAY * par = ps_new();
  ps_createargs(par,argc-1,&(argv[1]));
  int num_threads = 1;
  ps_flint(*par,"num_threads",&num_threads);
  
#ifdef _OPENMP
  omp_set_num_threads(num_threads);
  printf("Number of OMP threads in use = %d\n",omp_get_max_threads());
#endif

  // params
  int ihmax=0;
  float dh=12.0;
  int maxoff=4;
  int nt=100000;
  
  // bulk, buoy
  IPNT gsa_b; IPNT gea_b;
  gsa_b[0]=0; gea_b[0]=290;
  gsa_b[1]=0; gea_b[1]=766;
  int nb=(gea_b[0]-gsa_b[0]+1)*(gea_b[1]-gsa_b[1]+1);
  float * bulk = (float *)malloc((2*ihmax+1)*nb*sizeof(float));
  float * buoy = (float *)malloc(nb*sizeof(float));
  for (int i=0; i<ihmax*nb; i++) bulk[i]=0.0;
  if (ihmax > 0) {
    for (int i=ihmax*nb; i<(ihmax+1)*nb; i++) bulk[i]=4.0/dh;
  }
  else {
    for (int i=ihmax*nb; i<(ihmax+1)*nb; i++) bulk[i]=4.0;    
  }
  for (int i=(ihmax+1)*nb; i<(2*ihmax+1)*nb; i++) bulk[i]=0.0;
  for (int i=0; i<nb; i++) buoy[i]=1.0f;

  // random initialization of acoustic fields
  srand(getpid());
  
  // p0
  IPNT gsa_p0; IPNT gea_p0;
  gsa_p0[0]=-3; gea_p0[0]=293;
  gsa_p0[1]=1; gea_p0[1]=765;
  int np0=(gea_p0[0]-gsa_p0[0]+1)*(gea_p0[1]-gsa_p0[1]+1);
  float * p0 = (float *)malloc(np0*sizeof(float));
  
  // p1
  IPNT gsa_p1; IPNT gea_p1;
  gsa_p1[0]=1; gea_p1[0]=289;
  gsa_p1[1]=-3; gea_p1[1]=769;
  int np1=(gea_p1[0]-gsa_p1[0]+1)*(gea_p1[1]-gsa_p1[1]+1);
  float * p1 = (float *)malloc(np1*sizeof(float));
  
  // common comp grid
  IPNT gsc_p; IPNT gec_p;
  gsc_p[0]=1; gec_p[0]=289;
  gsc_p[1]=1; gec_p[1]=765;

  // v0
  IPNT gsa_v0; IPNT gea_v0;
  IPNT gsc_v0; IPNT gec_v0;
  gsa_v0[0]=-3; gea_v0[0]=292;
  gsa_v0[1]=1; gea_v0[1]=765;
  gsc_v0[0]=0; gec_v0[0]=289;
  gsc_v0[1]=1; gec_v0[1]=765;
  int nv0=(gea_v0[0]-gsa_v0[0]+1)*(gea_v0[1]-gsa_v0[1]+1);
  float * v0 = (float *)malloc(nv0*sizeof(float));
  
  // v1
  IPNT gsa_v1; IPNT gea_v1;  
  IPNT gsc_v1; IPNT gec_v1;  
  gsa_v1[0]=1; gea_v1[0]=289;
  gsa_v1[1]=-3; gea_v1[1]=768;
  gsc_v1[0]=1; gec_v1[0]=289;
  gsc_v1[1]=0; gec_v1[1]=765;
  int nv1=(gea_v1[0]-gsa_v1[0]+1)*(gea_v1[1]-gsa_v1[1]+1);
  float * v1 = (float *)malloc(nv1*sizeof(float));
  
  // ep, epp == 1
  // ep[0][gsc_p[0],gec_p[0]]
  // ep[1][gsc_p[1],gec_p[1]]
  float ** ep = (float **)malloc(2*sizeof(float *));
  ep[0]=(float *)malloc((gec_p[0]-gsc_p[0]+1)*sizeof(float));
  ep[1]=(float *)malloc((gec_p[1]-gsc_p[1]+1)*sizeof(float));
  float ** epp = (float **)malloc(2*sizeof(float *));
  epp[0]=(float *)malloc((gec_p[0]-gsc_p[0]+1)*sizeof(float));
  epp[1]=(float *)malloc((gec_p[1]-gsc_p[1]+1)*sizeof(float));
  
  // ev, evp
  // ev[0][gsc_v0[0],gec_v0[0]]
  // ev[1][gsc_v1[1],gec_v0[1]]
  float ** ev = (float **)malloc(2*sizeof(float *));
  ev[0]=(float *)malloc((gec_v0[0]-gsc_v0[0]+1)*sizeof(float));
  ev[1]=(float *)malloc((gec_v1[1]-gsc_v1[1]+1)*sizeof(float));
  float ** evp = (float **)malloc(2*sizeof(float *));
  evp[0]=(float *)malloc((gec_v0[0]-gsc_v0[0]+1)*sizeof(float));
  evp[1]=(float *)malloc((gec_v1[1]-gsc_v1[1]+1)*sizeof(float));

  // sdiv, sdivb [gsc_p[0], gec_p[0]] x [gsc_p[1], gec_p[1]]
  float * sdiv = (float *)malloc((gec_p[0]-gsc_p[0]+1)*(gec_p[1]-gsc_p[1]+1)*sizeof(float));
  float * sdivb = (float *)malloc((gec_p[0]-gsc_p[0]+1)*(gec_p[1]-gsc_p[1]+1)*sizeof(float));

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif  
  for (int i=0; i<np0; i++) p0[i]=-0.5+ rand()/(RAND_MAX+1.0);
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
  for (int i=0; i<np1; i++) p1[i]=-0.5+ rand()/(RAND_MAX+1.0);
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
  for (int i=0; i<nv0; i++) v0[i]=-0.5+ rand()/(RAND_MAX+1.0);
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
  for (int i=0; i<nv1; i++) v1[i]=-0.5+ rand()/(RAND_MAX+1.0);
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
  for (int i0=gsc_p[0]; i0<=gec_p[0]; i0++) ep[0][i0-gsc_p[0]]=1.0f;
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
  for (int i1=gsc_p[1]; i1<=gec_p[1]; i1++) ep[1][i1-gsc_p[1]]=1.0f;
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
  for (int i0=gsc_p[0]; i0<=gec_p[0]; i0++) epp[0][i0-gsc_p[0]]=1.0f;
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
  for (int i1=gsc_p[1]; i1<=gec_p[1]; i1++) epp[1][i1-gsc_p[1]]=1.0f;
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
  for (int i0=gsc_v0[0]; i0<=gec_v0[0];i0++) ev[0][i0-gsc_v0[0]]=1.0f;
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
  for (int i1=gsc_v1[1]; i1<=gec_v0[0];i1++) ev[1][i1-gsc_v1[1]]=1.0f;
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
  for (int i0=gsc_v0[0]; i0<=gec_v0[0];i0++) evp[0][i0-gsc_v0[0]]=1.0f;
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
  for (int i1=gsc_v1[1]; i1<=gec_v0[0];i1++) evp[1][i1-gsc_v1[1]]=1.0f;
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
  for (int i=0; i<(gec_p[0]-gsc_p[0]+1)*(gec_p[1]-gsc_p[1]+1); i++) sdiv[i]=-0.5+ rand()/(RAND_MAX+1.0);
  
  // gradp0 [gsc_v0[0],gec_v0[0]]
  float * gradp0 = (float *)malloc((gec_v0[0]-gsc_v0[0]+1)*sizeof(float));
  // no initialization needed
  
  // gradp1 [gsc_v1[1],gec_v1[1]]
  float * gradp1 = (float *)malloc((gec_v1[1]-gsc_v1[1]+1)*sizeof(float));
  // no initialization needed
  
  // lbc, rbc = all 1's
  int * lbc = (int *)malloc(2*sizeof(int));
  int * rbc = (int *)malloc(2*sizeof(int));  
  lbc[0]=1; lbc[1]=1;
  rbc[0]=1; rbc[1]=1;

  // c 2 x maxoff
  float ** c = (float **)malloc(2*sizeof(float));
  c[0] = (float *)malloc(maxoff*sizeof(float));
  c[1] = (float *)malloc(maxoff*sizeof(float));
  c[0][0]=0.070492;
  c[0][1]=-0.00469947;
  c[0][2]=0.000563936;
  c[0][3]=-4.11032e-05;
  c[1][0]=0.070492;
  c[1][1]=-0.00469947;
  c[1][2]=0.000563936;
  c[1][3]=-4.11032e-05;

#ifdef DIVONLY
  for (int k=0; k<nt; k++) {
    v0[((gea_v0[1]-gsa_v0[1]+1)*(gea_v0[0]-gsa_v0[0]+1)/2 +(gea_v0[0]-gsa_v0[0]+1)/2)]
      +=-0.5+ rand()/(RAND_MAX+1.0);
    asgdiv(v0, v1,
	   sdiv,
	   gsc_p, gec_p, 
	   gsa_v0, gea_v0,
	   gsa_v1, gea_v1,
	   maxoff, c);
  }
#endif
#ifdef PUPDONLY
  for (int k=0;k<nt;k++) {
    sdiv[0]+=-0.5+ rand()/(RAND_MAX+1.0);
    asgpupd(bulk,
	    p0, p1,
	    ep, epp,
	    sdiv,
	    gsa_b, gea_b,
	    gsa_p0, gea_p0,
	    gsa_p1, gea_p1,
	    gsc_p,  gec_p);
  }
#endif
      
#ifdef TIMESTEP
  // time loop
  for (int k=0;k<nt;k++) {
    //    fprintf(stderr,"time step k = %d\n",k);
    p0[((gea_p0[1]-gsa_p0[1]+1)*(gea_p0[0]-gsa_p0[0]+1)/2 +(gea_p0[0]-gsa_p0[0]+1)/2)]
      +=-0.5+ rand()/(RAND_MAX+1.0);
    asg2dtaptest(bulk,buoy,
		 p0, p1,
		 v0, v1,
		 ep, epp,
		 ev, evp,
		 sdiv, sdivb,
		 gradp0,
		 gradp1,
		 gsa_b,  gea_b,
		 gsa_p0, gea_p0,
		 gsa_p1, gea_p1,
		 gsc_p,  gec_p, 
		 gsa_v0, gea_v0,
		 gsc_v0, gec_v0,
		 gsa_v1, gea_v1,
		 gsc_v1, gec_v1,		  
		 lbc, rbc,
		 dh, ihmax,
		 maxoff,  c);
  }
#endif
  exit(0);
}
	 
	 
