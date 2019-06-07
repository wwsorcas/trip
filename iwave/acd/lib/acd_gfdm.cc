#include "acd_gfdm.h"

#undef OLD
#undef DOLD

/*--- time step functions ---------------------------------------------------*/
extern "C" {
void acd_2d_d(float **uc, float **ucd, float **up, float **upd, float **csq, 
              float **csqd, int *s, int *e, float **c, int k, int *lbc, int *rbc,
  	      float *lap, float *lapd);

void acd_2d_b(float **uc, float **ucb, float **up, float **upb, float **csq, 
	      float **csqb, int *s, int *e, float **c, int k, int *lbc, int *rbc,
	      float *lap, float *lapb);

#ifdef DOLD

void acde_2d_d(float **uc, float **ucd, float **up, float **upd, float **csq,
	       float ***csqd, int *s, int *e, int *z, float ** c, int k,
	       int *lbc, int *rbc, float *lap, float *csqdlap);

void acde_2d_b(float **uc, float **ucd, float **up, float **upd, float **csq,
               float ***csqd, int *s, int *e, int *z, float ** c, int k,
	       int *lbc, int *rbc, float *lap);

#else

void acde_2d_d(float **uc, float **ucd, float **up, float **upd, float ***csq,
	       float ***csqd, int *s, int *e, int *z, float ** c, int k, int pbg,
	       int *lbc, int *rbc, float *lap, float *lapd, float *csqdlap,
	       float *csqlapd);

void acde_2d_b(float **uc, float **ucd, float **up, float **upd, float ***csq,
               float ***csqd, int *s, int *e, int *z, float ** c, int k, int pbg,
	       int *lbc, int *rbc, float *lap, float * lapb);

#endif

void acd_2d_2_d(float **uc, float **ucd, float **up, float **upd, float **csq,
		float **csqd, int *s, int *e, float c0, float *c1);

void acd_2d_2_b(float **uc, float **ucb, float **up, float **upb, float **csq,
		float **csqb, int *s, int *e, float c0, float *c1);

void acd_2d_4_d(float **uc, float **ucd, float **up, float **upd, float **csq,
		float **csqd, int *s, int *e, float c0, float *c1, float *c2, 
		int *lbc, int *rbc);

void acd_2d_4_b(float **uc, float **ucb, float **up, float **upb, float **csq,
		float **csqb, int *s, int *e, float c0, float *c1, float *c2, 
		int *lbc, int *rbc);

void acd_2d_8_d(float **uc, float **ucd, float **up, float **upd, float **csq,
		float **csqd, int *s, int *e, float c0, float *c1, float *c2, 
		float *c3, float *c4, int *lbc, int *rbc);

void acd_2d_8_b(float **uc, float **ucb, float **up, float **upb, float **csq,
		float **csqb, int *s, int *e, float c0, float *c1, float *c2, 
		float *c3, float *c4, int *lbc, int *rbc);

void acde_2d_4_d(float **uc, float **ucd, float **up, float **upd, float **csq,
		float ***csqd, int *s, int *e, int *z, float c0, float *c1, float *c2, 
		int *lbc, int *rbc, float *dx, float *lap, float *csqdlap);

void acde_2d_4_b(float **uc, float **ucb, float **up, float **upb, float **csq,
		float ***csqb, int *s, int *e, int *z, float c0, float *c1, float *c2, 
		int *lbc, int *rbc, float *dx, float *lap);

}


int acd_tsfm(RDOM * p, RDOM * r, int ia, void * fdpars) {

  /* swap pointer workspace */
  ireal tmp;
  IPNT i;

  /* pointers for 2D case */
  register ireal ** restrict uc2;
  register ireal ** restrict up2;
  register ireal ** restrict csq2;
  register ireal ** restrict uc2d;
  register ireal ** restrict up2d;
  register ireal ** restrict csq2d;
  /* pointers for 3D case */
  //  register ireal *** restrict uc3;
  //  register ireal *** restrict up3;
  register ireal *** restrict csq3;
  //  register ireal *** restrict uc3d;
  //  register ireal *** restrict up3d;
  register ireal *** restrict csq3d;
  
  int ndim;                       /* problem dmn */
  IPNT s, s0, cs0;                /* loop starts */ 
  IPNT e, e0, ce0;                /* loop ends */

  /* acd struct */
  ACD_TS_PARS * acdpars = (ACD_TS_PARS *)fdpars;
  /* extract dimn info */
  ra_ndim (&(r->_s[D_UC]),&ndim);
  ra_gse  (&(r->_s[D_UC]),s,e);
  ra_a_gse(&(p->_s[D_CSQ]),cs0,ce0);
  ra_a_gse(&(r->_s[D_UC]),s0,e0);

  /* pass extended axes */
  for (int i=0;i<acdpars->next;i++) {
    s[i+acdpars->ndim]=cs0[i+acdpars->ndim];
    e[i+acdpars->ndim]=ce0[i+acdpars->ndim];
  }

  
  if (ndim == 2) {
    
    /* 2D computational arrays */
    uc2    = (r->_s)[D_UC ]._s2;
    up2    = (r->_s)[D_UP ]._s2;
    csq2   = (r->_s)[D_CSQ]._s2;
    csq3   = (r->_s)[D_CSQ]._s3;
    uc2d   = (p->_s)[D_UC ]._s2;
    up2d   = (p->_s)[D_UP ]._s2;
    csq2d  = (p->_s)[D_CSQ]._s2;
    if (acdpars->next == 1) {
      csq3d  = (p->_s)[D_CSQ]._s3;
      csq2   = ((r->_s)[D_CSQ]._s3)[(acdpars->z)[2]];
    }

#ifndef OLD
    // extended modeling
    if (acdpars->next == 1) {

#ifdef DOLD
      acde_2d_d(uc2, uc2d,
		up2, up2d,
		csq2, csq3d,
		s, e,
		acdpars->z,
		acdpars->c,
		acdpars->k,
		acdpars->lbc, 
		acdpars->rbc,
		acdpars->lap,
		acdpars->lap1);
#else
      acde_2d_d(uc2, uc2d,
		up2, up2d,
		csq3, csq3d,
		s, e,
		acdpars->z,
		acdpars->c,
		acdpars->k,
		acdpars->pbg,
		acdpars->lbc, 
		acdpars->rbc,
		acdpars->lap,
		acdpars->lap1,
		acdpars->lap2,
		acdpars->lap3);
#endif
    }
    // physical modeling
    else {
      acd_2d_d(uc2, uc2d,
	       up2, up2d,
	       csq2, csq2d,
	       s, e,     
	       acdpars->c, 
	       acdpars->k,
	       acdpars->lbc,
	       acdpars->rbc,
               acdpars->lap,
               acdpars->lap1);
    }

#else //end NEW branch

    /* 2nd order case */
    if (acdpars->k == 1) {
      acd_2d_2_d(uc2, uc2d,
		 up2, up2d,
		 csq2, csq2d, 
		 s, e, 
		 acdpars->c0, 
		 acdpars->c1);
    }
    /* 4th order case */
    else if (acdpars->k == 2) {
      if (acdpars->next == 1) {
	acde_2d_4_d(uc2, uc2d,
		    up2, up2d,
		    csq2, csq3d,
		    s, e, acdpars->z,
		    acdpars->c0, 
		    acdpars->c1, acdpars->c2,
		    acdpars->lbc, acdpars->rbc,
		    acdpars->dx, acdpars->lap, acdpars->lap1);
      } else {
	acd_2d_4_d(uc2, uc2d,
		   up2, up2d,
		   csq2, csq2d,
		   s, e, 
		   acdpars->c0, 
		   acdpars->c1, acdpars->c2,
		   acdpars->lbc, acdpars->rbc);
      }
    }
    /* 8th order case */
    else if (acdpars->k == 4) {
      acd_2d_8_d(uc2, uc2d,
		 up2, up2d,
		 csq2, csq2d,
		 s, e, 
		 acdpars->c0, 
		 acdpars->c1, acdpars->c2,
		 acdpars->c3, acdpars->c4,
		 acdpars->lbc, acdpars->rbc);
    }
    else {
      fprintf(stderr,"ERROR: acd_step\n");
      fprintf(stderr,"called with half-order != 1, 2, 4\n");
      return E_BADINPUT;
    }

#endif //end OLD branch

    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	tmp=uc2d[i[1]][i[0]];
	uc2d[i[1]][i[0]]=up2d[i[1]][i[0]];
	up2d[i[1]][i[0]]=tmp;
          tmp=uc2[i[1]][i[0]];
          uc2[i[1]][i[0]]=up2[i[1]][i[0]];
          up2[i[1]][i[0]]=tmp;
      }
    }
  }
  else if (ndim == 3) {
    
      /*
    uc3    = (r->_s)[D_UC ]._s3;
    up3    = (r->_s)[D_UP ]._s3;
    csq3   = (r->_s)[D_CSQ]._s3;
    uc3d   = (p->_s)[D_UC ]._s3;
    up3d   = (p->_s)[D_UP ]._s3;
    csq3d  = (p->_s)[D_CSQ]._s3;
      */    

    // 2nd order case 
    if (acdpars->k == 1) {
      /*
      acd_3d_2(uc3, up3, csq3, 
	       s, e, 
	       acdpars->c0, 
	       acdpars->c1);
      */
    }
    // 4th order case
    else if (acdpars->k == 2) {
      /*
      acd_3d_4(uc3, up3, csq3,
	       s, e, 
	       acdpars->c0, 
	       acdpars->c1, acdpars->c2,
	       acdpars->lbc, acdpars->rbc);
      */
    }
    // 8th order case
    else if (acdpars->k == 4) {
      /*
      acd_3d_8(uc3, up3, csq3,
	       s, e, 
	       acdpars->c0, 
	       acdpars->c1, acdpars->c2,
	       acdpars->c3, acdpars->c4,
	       acdpars->lbc, acdpars->rbc);
      */
    }
    else {
      fprintf(stderr,"ERROR: acd_step\n");
      fprintf(stderr,"called with half-order != 1, 2, 4\n");
      return E_BADINPUT;
    }

    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	  tmp=((p->_s)[D_UC]._s3)[i[2]][i[1]][i[0]];
	  ((p->_s)[D_UC]._s3)[i[2]][i[1]][i[0]]=((p->_s)[D_UP]._s3)[i[2]][i[1]][i[0]];
	  ((p->_s)[D_UP]._s3)[i[2]][i[1]][i[0]]=tmp;
	}
      }
    }
  }
  else {
    fprintf(stderr,"ERROR: acd_step\n");
    fprintf(stderr,"called with space dim != 2 or 3\n");
    return E_BADINPUT;
  }
  
  return 0;
}

int acd_tsam(RDOM * p, RDOM * r, int ia, void * fdpars) {

  /* swap pointer workspace */
  ireal tmp;
  IPNT i;
  IPNT n;

  /* pointers for 2D case */
  register ireal ** restrict uc2;
  register ireal ** restrict up2;
  register ireal ** restrict csq2;
  register ireal ** restrict uc2d;
  register ireal ** restrict up2d;
  register ireal ** restrict csq2d;
  /* pointers for 3D case */
  //  register ireal *** restrict uc3;
  //  register ireal *** restrict up3;
  register ireal *** restrict csq3;
  //  register ireal *** restrict uc3d;
  //  register ireal *** restrict up3d;
  register ireal *** restrict csq3d;

  int ndim;                       /* problem dmn */
  IPNT s, s0, cs0;                /* loop starts */ 
  IPNT e, e0, ce0;                /* loop ends */

  //////////
  //  cerr<<"acg_gfdm"<<endl;

  /* acd struct */
  ACD_TS_PARS * acdpars = (ACD_TS_PARS *)fdpars;
  /* extract dimn info */
  ra_ndim  (&(r->_s[D_UC]),&ndim);
  ra_gse   (&(r->_s[D_UC]),s,e);
  ra_a_gse (&(p->_s[D_CSQ]),cs0,ce0);
  ra_a_gse (&(r->_s[D_UC]),s0,e0);
  ra_a_size(&(r->_s[D_UC]),n);
  
  /* pass extended axes */
  for (int i=0;i<acdpars->next;i++) {
    s[i+acdpars->ndim]=cs0[i+acdpars->ndim];
    e[i+acdpars->ndim]=ce0[i+acdpars->ndim];
  }
  

  if (ndim == 2) {
	  
    /* 2D computational arrays */
    uc2    = (r->_s)[D_UC ]._s2;
    up2    = (r->_s)[D_UP ]._s2;
    csq2   = (r->_s)[D_CSQ]._s2;
    csq3   = (r->_s)[D_CSQ]._s3;
    uc2d   = (p->_s)[D_UC ]._s2;
    up2d   = (p->_s)[D_UP ]._s2;
    csq2d  = (p->_s)[D_CSQ]._s2;
    if (acdpars->next == 1) {
      csq3d  = (p->_s)[D_CSQ]._s3;
      csq2   = ((r->_s)[D_CSQ]._s3)[(acdpars->z)[2]];
    }

    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	tmp=uc2d[i[1]][i[0]];
	uc2d[i[1]][i[0]]=up2d[i[1]][i[0]];
	up2d[i[1]][i[0]]=tmp;
      }
    }

#ifndef OLD
    // extended imaging
    if (acdpars->next == 1) {
#ifdef DOLD
      acde_2d_b(uc2, uc2d,
		up2, up2d,
		csq2, csq3d,
		s, e,
		acdpars->z,
		acdpars->c,
		acdpars->k,
		acdpars->lbc, 
		acdpars->rbc,
		acdpars->lap);
#else
      acde_2d_b(uc2, uc2d,
		up2, up2d,
		csq3, csq3d,
		s, e,
		acdpars->z,
		acdpars->c,
		acdpars->k,
		acdpars->pbg,
		acdpars->lbc, 
		acdpars->rbc,
		acdpars->lap,
		acdpars->lap1);
#endif
    }
    // physical imaging
    else {
      acd_2d_b(uc2, uc2d,
	       up2, up2d,
	       csq2, csq2d,
	       s, e,     
	       acdpars->c, 
	       acdpars->k,
	       acdpars->lbc,
	       acdpars->rbc,
      	       acdpars->lap,
      	       acdpars->lap1);
    }

#else // end NEW branch
    
    /* 2nd order case */ 
    if (acdpars->k == 1) {
      acd_2d_2_b(uc2, uc2d,
		 up2, up2d,
		 csq2, csq2d, 
		 s, e, 
		 acdpars->c0, 
		 acdpars->c1);
    }
    /* 4th order case */
    else if (acdpars->k == 2) {
      if (acdpars->next == 1) {
	acde_2d_4_b(uc2, uc2d,
		    up2, up2d,
		    csq2, csq3d,
		    s, e, acdpars->z,
		    acdpars->c0, 
		    acdpars->c1, acdpars->c2,
		    acdpars->lbc, acdpars->rbc,
		    acdpars->dx, acdpars->lap);
      } else {
	acd_2d_4_b(uc2, uc2d,
		   up2, up2d,
		   csq2, csq2d,
		   s, e, 
		   acdpars->c0, 
		   acdpars->c1, acdpars->c2,
		   acdpars->lbc, acdpars->rbc);
      }
    }
    /* 8th order case */
    else if (acdpars->k == 4) {
      acd_2d_8_b(uc2, uc2d,
		 up2, up2d,
		 csq2, csq2d,
		 s, e, 
		 acdpars->c0, 
		 acdpars->c1, acdpars->c2,
		 acdpars->c3, acdpars->c4,
		 acdpars->lbc, acdpars->rbc);
    }
    else {
      fprintf(stderr,"ERROR: acd_step\n");
      fprintf(stderr,"called with half-order != 1, 2, 4\n");
      return E_BADINPUT;
    }
#endif //end OLD branch
  }
  else if (ndim == 3) {

    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	  tmp=((p->_s)[D_UC]._s3)[i[2]][i[1]][i[0]];
	  ((p->_s)[D_UC]._s3)[i[2]][i[1]][i[0]]=((p->_s)[D_UP]._s3)[i[2]][i[1]][i[0]];
	  ((p->_s)[D_UP]._s3)[i[2]][i[1]][i[0]]=tmp;
	}
      }
    }

/*
    uc3    = (r->_s)[D_UC ]._s3;
    up3    = (r->_s)[D_UP ]._s3;
    csq3   = (r->_s)[D_CSQ]._s3;
    uc3d   = (p->_s)[D_UC ]._s3;
    up3d   = (p->_s)[D_UP ]._s3;
    csq3d  = (p->_s)[D_CSQ]._s3;
*/
  
    // 2nd order case 
    if (acdpars->k == 1) {
      /*
      acd_3d_2(uc3, up3, csq3, 
	       s, e, 
	       acdpars->c0, 
	       acdpars->c1);
      */
    }
    // 4th order case
    else if (acdpars->k == 2) {
      /*
      acd_3d_4(uc3, up3, csq3,
	       s, e, 
	       acdpars->c0, 
	       acdpars->c1, acdpars->c2,
	       acdpars->lbc, acdpars->rbc);
      */
    }
    // 8th order case
    else if (acdpars->k == 4) {
      /*
      acd_3d_8(uc3, up3, csq3,
	       s, e, 
	       acdpars->c0, 
	       acdpars->c1, acdpars->c2,
	       acdpars->c3, acdpars->c4,
	       acdpars->lbc, acdpars->rbc);
      */
    } 
    else {
      fprintf(stderr,"ERROR: acd_step\n");
      fprintf(stderr,"called with half-order != 1, 2, 4\n");
      return E_BADINPUT;
    }

  }
  else {
    fprintf(stderr,"ERROR: acd_step\n");
    fprintf(stderr,"called with space dim != 2 or 3\n");
    return E_BADINPUT;
  }
  
  return 0;
}
