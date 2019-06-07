#include "cstd.h"

/* version with explicit index arithmetic
 */

#define OMP_CHUNK_SIZE 32

#define my_min(a, b) ((a) < (b) ? (a) : (b))
#define my_max(a, b) ((a) > (b) ? (a) : (b))

void asg2dtaptest(float *  restrict  bulk, float *  restrict  buoy,
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
		  int maxoff, float ** restrict  c) {		  
		  //		  int maxoff, float ** restrict  c,
		  //		  FILE * stream) {

  int * sh = (int *)malloc((2*ihmax+1)*sizeof(int));
  int * eh = (int *)malloc((2*ihmax+1)*sizeof(int));  
  for (int ih=-ihmax; ih<=ihmax; ih++) {
    //      sh = ((ih <= 0) ?  gsc_p[1]-ih : gsc_p[1]);
    //      eh = ((ih >= 0) ?  gec_p[1]-ih : gec_p[1]);
    if (ih <=0) {
      sh[ih+ihmax] = gsc_p[1]; eh[ih+ihmax] = gec_p[1] + 2*ih;
    }
    else {
      sh[ih+ihmax] = gsc_p[1] + 2*ih; eh[ih+ihmax] = gec_p[1];
    }
  }
  
  int gscp1=gsc_p[1];
  int gscp0=gsc_p[0];
	
  // fprintf(stderr,"1\n");
#pragma ivdep
#pragma GCC ivdep
  for (int i=0;i<(gec_p[1]-gsc_p[1]+1)*(gec_p[0]-gsc_p[0]+1);i++) {
    sdiv[i]=0.0f;
  }
  
  //  memset((char *)sdiv,  0, (gec_p[1]-gsc_p[1]+1)*(gec_p[0]-gsc_p[0]+1)*sizeof(float));
  
  #ifdef _OPENMP
  #pragma omp parallel default(shared)
  #endif
  {
    // fprintf(stderr,"2\n");
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
    for (int i1=gsc_p[1]; i1<=gec_p[1]; i1++) {
      for (int ioff = 0; ioff<maxoff; ioff++) {
#pragma ivdep
#pragma GCC ivdep
	for (int i0=gsc_p[0]; i0 <= gec_p[0]; i0++) {
	  sdiv[(i1-gsc_p[1])*(gec_p[0]-gsc_p[0]+1) + i0-gsc_p[0]] +=
	    c[0][ioff]*(v0[(i1-gsa_v0[1])*(gea_v0[0]-gsa_v0[0]+1) + i0+ioff  -gsa_v0[0]] -
			v0[(i1-gsa_v0[1])*(gea_v0[0]-gsa_v0[0]+1) + i0-ioff-1-gsa_v0[0]]) +
	    c[1][ioff]*(v1[(i1+ioff-gsa_v1[1])*(gea_v1[0]-gsa_v1[0]+1) +   i0-gsa_v1[0]] -
			v1[(i1-ioff-1-gsa_v1[1])*(gea_v1[0]-gsa_v1[0]+1) + i0-gsa_v1[0]]);
	}
      }
#pragma ivdep
#pragma GCC ivdep
      for (int i0=gsc_p[0]; i0 <= gec_p[0]; i0++) {
	p0[(i1-gsa_p0[1])*(gea_p0[0]-gsa_p0[0]+1) + i0-gsa_p0[0]] =
	  (ep[0][i0-gscp0]*p0[(i1-gsa_p0[1])*(gea_p0[0]-gsa_p0[0]+1) + i0-gsa_p0[0]] -
	   bulk[(i1-gsa_b[1])*(gea_b[0]-gsa_b[0]+1) + i0-gsa_b[0]]*
	   sdiv[(i1-gscp1)*(gec_p[0]-gsc_p[0]+1)+i0-gscp0])*epp[0][i0-gscp0];
	p1[(i1-gsa_p1[1])*(gea_p1[0]-gsa_p1[0]+1) + i0-gsa_p1[0]] =
	  (ep[1][i1-gscp0]*p1[(i1-gsa_p1[1])*(gea_p1[0]-gsa_p1[0]+1) + i0-gsa_p1[0]] -
	   bulk[(i1-gsa_b[1])*(gea_b[0]-gsa_b[0]+1) + i0-gsa_b[0]]*
	   sdiv[(i1-gscp1)*(gec_p[0]-gsc_p[0]+1)+i0-gscp0])*epp[1][i1-gscp1]; 
      }      
    }

  } // end par region 

  // separate logic for physical sim

  if (ihmax != 0) {
    
    for (int ih=-ihmax; ih<=ihmax; ih++) {
      
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
      
      {     // begin parallel region
	
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
	for (int i1=sh[ih+ihmax]; i1<=eh[ih+ihmax]; i1++) {
#pragma ivdep
#pragma GCC ivdep
	  for (int i0=gsc_p[0]; i0 <= gec_p[0]; i0++) {      
	    sdivb[(i1-gsc_p[1])*(gec_p[0]-gsc_p[0]+1) + i0-gsc_p[0]] +=
	      dh * bulk[(ih+ihmax)*(gea_b[1]-gsa_b[1]+1)*(gea_b[0]-gsa_b[0]+1) +(i1-ih-gsa_b[1])*(gea_b[0]-gsa_b[0]+1) + i0-gsa_b[0]] * sdiv[(i1-gsc_p[1]-2*ih)*(gec_p[0]-gsc_p[0]+1) +i0-gsc_p[0]];
	  } 
	}
      } // end par region
      
    } // end ih loop
    
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
    
    {     // begin parallel region
      
      // fprintf(stderr,"4\n");
#ifdef _OPENMP
#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif
      for (int i1=gsc_p[1]; i1<=gec_p[1]; i1++) {      
#pragma ivdep
#pragma GCC ivdep
	for (int i0=gsc_p[0]; i0 <= gec_p[0]; i0++) {
	  p0[(i1-gsa_p0[1])*(gea_p0[0]-gsa_p0[0]+1) + i0-gsa_p0[0]] =
	    (ep[0][i0-gscp0]*p0[(i1-gsa_p0[1])*(gea_p0[0]-gsa_p0[0]+1) + i0-gsa_p0[0]] -
	     sdivb[(i1-gscp1)*(gec_p[0]-gsc_p[0]+1)+i0-gscp0])*epp[0][i0-gscp0];
	  p1[(i1-gsa_p1[1])*(gea_p1[0]-gsa_p1[0]+1) + i0-gsa_p1[0]] =
	    (ep[1][i1-gscp0]*p1[(i1-gsa_p1[1])*(gea_p1[0]-gsa_p1[0]+1) + i0-gsa_p1[0]] -
	     sdivb[(i1-gscp1)*(gec_p[0]-gsc_p[0]+1)+i0-gscp0])*epp[1][i1-gscp1]; 
	}
      }
    
    }
  } // end ihmax>0 branch

  /*
    // fprintf(stderr,"5-8\n");
  if (lbc[0]) {
    //#ifdef _OPENMP
    //#pragma omp for schedule(static,OMP_CHUNK_SIZE)
    //#endif
    for (int i1=gsc_p[1];i1<=gec_p[1];i1++) {
      p0[i1][gsc_p[0]-1]=0.0f;
      for (int ioff=1;ioff<maxoff;ioff++) {
	p0[i1][gsc_p[0]-ioff-1]=-p0[i1][gsc_p[0]+ioff-1];
      }
    }
  }

  if (rbc[0]) {
    //#ifdef _OPENMP
    //#pragma omp for schedule(static,OMP_CHUNK_SIZE)
    //#endif
      for (int i1=gsc_p[1];i1<=gec_p[1];i1++) {
	p0[i1][gec_p[0]+1]=0.0f;
	for (int ioff=1;ioff<maxoff;ioff++) {
	  p0[i1][gec_p[0]+ioff+1]=-p0[i1][gec_p[0]-ioff+1];
	}
      }
    }

  if (lbc[1]) {
#pragma ivdep
#pragma GCC ivdep
    for (int i0=gsc_p[0];i0<=gec_p[0];i0++) {
      p1[gsc_p[1]-1][i0]=0.0f;
    }
    //#ifdef _OPENMP
    //#pragma omp for schedule(static,1)
    //#endif
    for (int ioff=1;ioff<maxoff;ioff++) {
#pragma ivdep
#pragma GCC ivdep
      for (int i0=gsc_p[0];i0<=gec_p[0];i0++) {
	p1[gsc_p[1]-ioff-1][i0]=-p1[gsc_p[1]+ioff-1][i0];
      }
    }
  }

  if (rbc[1]) {
#pragma ivdep
#pragma GCC ivdep
    for (int i0=gsc_p[0];i0<=gec_p[0];i0++) {
      p1[gec_p[1]+1][i0]=0.0f;
    }
    //#ifdef _OPENMP
    //#pragma omp for schedule(static,1)
    //#endif
    for (int ioff=1;ioff<maxoff;ioff++) {
#pragma ivdep
#pragma GCC ivdep
      for (int i0=gsc_p[0];i0<=gec_p[0];i0++) {
	p1[gec_p[1]+ioff+1][i0]=-p1[gec_p[1]-ioff+1][i0];
      }
    }
  }
  */

  /*
  //#ifdef _OPENMP
  //#pragma omp parallel default(shared)
  //#endif
  {
    
  // fprintf(stderr,"9\n");
    //#ifdef _OPENMP
    //#pragma omp for schedule(static, OMP_CHUNK_SIZE)
    //#endif
    for (int i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++ ) {

      memset((char *)gradp0, 0, (gec_v0[0]-gsc_v0[0]+1)*sizeof(float));

      for (int ioff=0; ioff<maxoff; ioff++) {
#pragma ivdep	
#pragma GCC ivdep
	for (int i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ ) {
	  gradp0[i0-gsc_v0[0]] += c[0][ioff]*
	    (p0[(i1-gsa_p0[1])*(gea_p0[0]-gsa_p0[0]+1) + i0+ioff+1-gsa_p0[0]]-
	     p0[(i1-gsa_p0[1])*(gea_p0[0]-gsa_p0[0]+1) + i0-ioff  -gsa_p0[0]]);
	}
      }
    
#pragma ivdep	
#pragma GCC ivdep
      for (int i0=gsc_v0[0]; i0 <= gec_v0[0]; i0++ ) {      
	v0[(i1-gsa_v0[1])*(gea_v0[0]-gsa_v0[0]+1) + i0-gsa_v0[0]] = evp[0][i0-gsc_v0[0]]*
	  (ev[0][i0-gsc_v0[0]]*v0[(i1-gsa_v0[1])*(gea_v0[0]-gsa_v0[0]+1) + i0-gsa_v0[0]] -
	   0.5f*(buoy[(i1-gsa_b[1])*(gea_b[0]-gsa_b[0]+1) + i0  -gsa_b[0]]+
		 buoy[(i1-gsa_b[1])*(gea_b[0]-gsa_b[0]+1) + i0+1-gsa_b[0]])*
		 gradp0[i0-gsc_v0[0]]);
      }
    }
    
    // fprintf(stderr,"12\n");
    //#ifdef _OPENMP
    //#pragma omp for schedule(static, OMP_CHUNK_SIZE)
    //#endif      
    for (int i1=gsc_v1[1]; i1 <= gec_v1[1]; i1++ ) {
      memset((char *)gradp1, 0, (gec_v1[0]-gsc_v1[0]+1)*sizeof(float));

      for (int ioff=0; ioff<maxoff; ioff++) {
#pragma ivdep
#pragma GCC ivdep
	for (int i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) {
	  gradp1[i0-gsc_v0[0]] += c[1][ioff]*
	    (p1[(i1+ioff+1-gsa_p1[1])*(gea_p1[0]-gsa_p1[0]+1) + i0-gsa_p1[0]]-
	     p1[(i1-ioff  -gsa_p1[1])*(gea_p1[0]-gsa_p1[0]+1) + i0-gsa_p1[0]]);
	}	
      }
#pragma ivdep
#pragma GCC ivdep
      for (int i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++ ) {      
	v1[(i1-gsa_v1[1])*(gea_v1[0]-gsa_v1[0]+1) + i0-gsa_v1[0]] = evp[1][i1-gsc_v1[1]]*
	  (ev[1][i1-gsc_v1[1]]*v1[(i1-gsa_v1[1])*(gea_v1[0]-gsa_v1[0]+1) + i0-gsa_v1[0]] -
	   0.5f*(buoy[(i1-gsa_b[1])*(gea_b[0]-gsa_b[0]+1) + i0  -gsa_b[0]]+
		 buoy[(i1+1-gsa_b[1])*(gea_b[0]-gsa_b[0]+1) + i0-gsa_b[0]])*
		 gradp1[i0-gsc_v1[0]]);
      }
    }

  } // end parallel region
  */
/*
    // fprintf(stderr,"15-18\n");
  if (lbc[0]) {
      //#ifdef _OPENMP
  //#pragma omp for schedule(static,OMP_CHUNK_SIZE) collapse(2)
  //#endif
    for (int i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++) {
      for (int ioff=1; ioff<maxoff;ioff++) {
	v0[i1][gsc_v0[0]-ioff]=v0[i1][gsc_v0[0]+ioff-1];
      }
    }
  }

  if (rbc[0]) {
    //#ifdef _OPENMP
    //#pragma omp for schedule(static,OMP_CHUNK_SIZE) collapse(2)
    //#endif
    for (int i1=gsc_v0[1]; i1 <= gec_v0[1]; i1++) {
      for (int ioff=1; ioff<maxoff;ioff++) {
	v0[i1][gec_v0[0]+ioff]=v0[i1][gec_v0[0]-ioff+1];
      }
    }
  }
    
  if (lbc[1]) {
    //#ifdef _OPENMP
    //#pragma omp for schedule(static,1)
    //#endif
    for (int ioff=1; ioff<maxoff;ioff++) {
#pragma ivdep
#pragma GCC ivdep
      for (int i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++) {
	v1[gsc_v1[1]-ioff][i0]=v1[gsc_v1[1]+ioff-1][i0];
      }
    }
  }
    
  if (rbc[1]) {
    //#ifdef _OPENMP
    //#pragma omp for schedule(static,1)
    //#endif
    for (int ioff=1; ioff<maxoff;ioff++) {
#pragma ivdep
#pragma GCC ivdep
      for (int i0=gsc_v1[0]; i0 <= gec_v1[0]; i0++) {
	v1[gec_v1[1]+ioff][i0]=v1[gec_v1[1]-ioff+1][i0];
      }
    }
  }
*/
  //  fprintf(stderr,"end fd\n");

}
