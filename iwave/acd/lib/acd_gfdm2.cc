#include "acd_gfdm2.h"

#undef OLD

/*--- time step functions ---------------------------------------------------*/
extern "C" {
    void acd_2d_d_d(float **uc, float **ucd0, float **ucd, float **ucdd, 
                    float **up, float **upd0, float **upd, float **updd, 
                    float **csq, float **csqd0, float **csqd, 
                    int *s, int *e, float **c, int k, int *lbc, int *rbc, 
                    float *lap, float *lapd0, float *lapd, float *lapdd);

    void acd_2d_d_b(float **uc, float **ucb, float **ucd, float **ucdb, 
                    float **up, float **upb, float **upd, float **updb, 
                    float **csq, float **csqb,float **csqd, 
                    int *s, int *e, float **c, int k, int *lbc, int *rbc,
		    float *lap, float *lapb, float *lapd, float *lapdb);

    void acd_2d_2_d_d(float **uc, float **ucd0, float **ucd, float **ucdd,
                      float **up, float **upd0, float **upd, float **updd,
                      float **csq, float **csqd0, float **csqd,
                      int *s, int *e, float c0, float *c1);

    void acd_2d_2_d_b(float **uc, float **ucb, float **ucd, float **ucdb,
                      float **up, float **upb, float **upd, float **updb,
                      float **csq, float **csqb,float **csqd,
                      int *s, int *e, float c0, float *c1);

    void acd_2d_4_d_d(float **uc, float **ucd0, float **ucd, float **ucdd,
                      float **up, float **upd0, float **upd, float **updd,
                      float **csq, float **csqd0, float **csqd,
                      int *s, int *e, float c0, float *c1, float *c2,
                      int *lbc, int *rbc);

    void acd_2d_4_d_b(float **uc, float **ucb, float **ucd, float **ucdb,
                      float **up, float **upb, float **upd, float **updb,
                      float **csq, float **csqb, float **csqd,
                      int *s, int *e, float c0, float *c1, float *c2,
                      int *lbc, int *rbc);
    
    void acd_2d_8_d_d(float **uc, float **ucd0, float **ucd, float **ucdd,
                      float **up, float **upd0, float **upd, float **updd,
                      float **csq, float **csqd0, float **csqd,
                      int *s, int *e, float c0, float *c1, float *c2,
                      float *c3, float *c4, int *lbc, int *rbc);

    void acd_2d_8_d_b(float **uc, float **ucb, float **ucd, float **ucdb,
                      float **up, float **upb, float **upd, float **updb,
                      float **csq, float **csqb, float **csqd,
                      int *s, int *e, float c0, float *c1, float *c2,
                      float *c3, float *c4, int *lbc, int *rbc);

    void acde_2d_4_d_d(float **uc,  float **ucd0,  float **ucd,  float **ucdd,
    		   	   	   float **up,  float **upd0,  float **upd,  float **updd,
    		   	   	   float **csq, float **csqd0, float ***csqd, 
    		   	   	   int *s, int *e, int *z, float c0, float *c1, float *c2, 
    		   	   	   int *lbc, int *rbc, float *dx, float *lap, float *lapd0,
    		   	   	  float *csqdlap, float *csqdlapd0);
    
    void acde_2d_4_d_b(float **uc, float **ucb, float **ucd, float **ucdb, 
    				   float **up, float **upb, float **upd, float **updb, 
    				   float **csq, float **csqb, float ***csqd, 
    				   int *s, int *e, int *z, float c0, float *c1, float *c2,
    				   int *lbc, int *rbc, float *dx, float *lapb);
    
}

int acd_tsfm2(RDOM * dd, RDOM *d0, RDOM * d, RDOM * r, int ia, void * fdpars) {

  /* swap pointer workspace */
  ireal tmp;
  IPNT i;

  /* pointers for 2D case */
  ireal ** restrict uc2;
  ireal ** restrict up2;
  ireal ** restrict csq2;
  ireal ** restrict uc2d;
  ireal ** restrict up2d;
  ireal ** restrict csq2d;
    
  ireal ** restrict uc2d0;
  ireal ** restrict up2d0;
  ireal ** restrict csq2d0;
  ireal ** restrict uc2dd;
  ireal ** restrict up2dd;

  /* pointers for 3D case */
  //  ireal *** restrict uc3;
  //  ireal *** restrict up3;
  //  ireal *** restrict csq3;
  //  ireal *** restrict uc3d;
  //  ireal *** restrict up3d;
  ireal *** restrict csq3d;
 
  int ndim;                       /* problem dmn */
  IPNT s, s0, cs0;                /* loop starts */
  IPNT e, e0, ce0;                /* loop ends */

  /* acd struct */
  ACD_TS_PARS * acdpars = (ACD_TS_PARS *)fdpars;
  /* extract dimn info */
  ra_ndim (&(r->_s[D_UC]),&ndim);
  ra_gse  (&(r->_s[D_UC]),s,e);
  ra_a_gse(&(d->_s[D_CSQ]),cs0,ce0);
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
    uc2d   = (d->_s)[D_UC ]._s2;
    up2d   = (d->_s)[D_UP ]._s2;
    csq2d  = (d->_s)[D_CSQ]._s2;
    csq3d  = (d->_s)[D_CSQ]._s3;

    uc2d0  = (d0->_s)[D_UC ]._s2;
    up2d0  = (d0->_s)[D_UP ]._s2;
    csq2d0 = (d0->_s)[D_CSQ]._s2;
    uc2dd  = (dd->_s)[D_UC ]._s2;
    up2dd  = (dd->_s)[D_UP ]._s2;

#ifndef OLD
    // extended modeling
    if (acdpars->next == 1) {
      acde_2d_4_d_d(uc2, uc2d0,  uc2d, uc2dd,
                    up2, up2d0, up2d,  up2dd,
                    csq2, csq2d0, csq3d,
                    s, e, acdpars->z,
                    acdpars->c0, acdpars->c1, acdpars->c2,
                    acdpars->lbc, acdpars->rbc,
                    acdpars->dx, &acdpars->lap[0], &acdpars->lap[acdpars->n[0]],
                    &acdpars->lap1[0], &acdpars->lap1[acdpars->n[0]]); 
    }
    // physical modeling
    else {
      acd_2d_d_d(   uc2, uc2d0,  uc2d, uc2dd,
                    up2, up2d0, up2d,  up2dd,
                    csq2, csq2d0, csq2d,
                    s, e, acdpars->c, acdpars->k,
                    acdpars->lbc, acdpars->rbc,
                    &acdpars->lap[0], &acdpars->lap[acdpars->n[0]],
                    &acdpars->lap1[0], &acdpars->lap1[acdpars->n[0]]);

    }

#else //end NEW branch

    /* 2nd order case */
    if (acdpars->k == 1) {
        acd_2d_2_d_d(uc2, uc2d0,  uc2d, uc2dd,
                     up2, up2d0, up2d,  up2dd,
                     csq2, csq2d0, csq2d,
                     s, e,
                     acdpars->c0, acdpars->c1);
    }
    /* 4th order case */
    else if (acdpars->k == 2) {
    	if (acdpars->next == 1) {
            acde_2d_4_d_d(uc2, uc2d0,  uc2d, uc2dd,
                    up2, up2d0, up2d,  up2dd,
                    csq2, csq2d0, csq3d,
                    s, e, acdpars->z,
                    acdpars->c0, acdpars->c1, acdpars->c2,
                    acdpars->lbc, acdpars->rbc,
                    acdpars->dx, &acdpars->lap[0], &acdpars->lap[acdpars->n[0]],
                    &acdpars->lap1[0], &acdpars->lap1[acdpars->n[0]]);
    	} else {
            acd_2d_4_d_d(uc2, uc2d0,  uc2d, uc2dd,
                    up2, up2d0, up2d,  up2dd,
                    csq2, csq2d0, csq2d,
                    s, e,
                    acdpars->c0, acdpars->c1, acdpars->c2,
                    acdpars->lbc, acdpars->rbc);
    	}
    }
    /* 8th order case */
    else if (acdpars->k == 4) {
        acd_2d_8_d_d(uc2, uc2d0,  uc2d, uc2dd,
                     up2, up2d0, up2d,  up2dd,
                    csq2, csq2d0, csq2d,
                    s, e, 
                    acdpars->c0, acdpars->c1, acdpars->c2,
                    acdpars->c3, acdpars->c4,
                    acdpars->lbc, acdpars->rbc);
    }
    else {
      fprintf(stderr,"ERROR: acd_step\n");
      fprintf(stderr,"called with half-order != 1, 2, 4\n");
      return E_BADINPUT;
    }

#endif //end OLD branch

 /* should add these lines when do timestepping */
 for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
          tmp=uc2dd[i[1]][i[0]];
          uc2dd[i[1]][i[0]]=up2dd[i[1]][i[0]];
          up2dd[i[1]][i[0]]=tmp;
          
          tmp=uc2d0[i[1]][i[0]];
          uc2d0[i[1]][i[0]]=up2d0[i[1]][i[0]];
          up2d0[i[1]][i[0]]=tmp;
          
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

/*    
 // should add these lines when do timestepping
 for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	  tmp=((p->_s)[D_UC]._s3)[i[2]][i[1]][i[0]];
	  ((p->_s)[D_UC]._s3)[i[2]][i[1]][i[0]]=((p->_s)[D_UP]._s3)[i[2]][i[1]][i[0]];
	  ((p->_s)[D_UP]._s3)[i[2]][i[1]][i[0]]=tmp;
	}
      }
    }
*/  }
  else {
    fprintf(stderr,"ERROR: acd_step\n");
    fprintf(stderr,"called with space dim != 2 or 3\n");
    return E_BADINPUT;
  }
  
  return 0;
}

int acd_tsam2(RDOM * db, RDOM *b, RDOM * d, RDOM * r, int ia, void * fdpars) {

  /* swap pointer workspace */
  ireal tmp;
  IPNT i;
  IPNT n;

  /* pointers for 2D case */
  ireal ** restrict uc2;
  ireal ** restrict up2;
  ireal ** restrict csq2;
  ireal ** restrict uc2d;
  ireal ** restrict up2d;
  ireal ** restrict csq2d;
    
  ireal ** restrict uc2b;
  ireal ** restrict up2b;
  ireal ** restrict csq2b;
  ireal ** restrict uc2db;
  ireal ** restrict up2db;

  /* pointers for 3D case */
  //  ireal *** restrict uc3;
  //  ireal *** restrict up3;
  //  ireal *** restrict csq3;
  //  ireal *** restrict uc3d;
  //  ireal *** restrict up3d;
  ireal *** restrict csq3d;

  int ndim;                       /* problem dmn */
  IPNT s, s0, cs0;                /* loop starts */
  IPNT e, e0, ce0;                /* loop ends */
  
  /* acd struct */
  ACD_TS_PARS * acdpars = (ACD_TS_PARS *)fdpars;
  /* extract dimn info */
  ra_ndim  (&(r->_s[D_UC]),&ndim);
  ra_gse   (&(r->_s[D_UC]),s,e);
  ra_a_gse (&(d->_s[D_CSQ]),cs0,ce0);
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
    uc2d   = (d->_s)[D_UC ]._s2;
    up2d   = (d->_s)[D_UP ]._s2;
    csq2d  = (d->_s)[D_CSQ]._s2;
    csq3d  = (d->_s)[D_CSQ]._s3;

    uc2b  = (b->_s)[D_UC ]._s2;
    up2b  = (b->_s)[D_UP ]._s2;
    csq2b = (b->_s)[D_CSQ]._s2;
    uc2db  = (db->_s)[D_UC ]._s2;
    up2db  = (db->_s)[D_UP ]._s2;
    
    for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
      for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
          tmp=uc2b[i[1]][i[0]];
          uc2b[i[1]][i[0]]=up2b[i[1]][i[0]];
          up2b[i[1]][i[0]]=tmp;
          
          tmp=uc2db[i[1]][i[0]];
          uc2db[i[1]][i[0]]=up2db[i[1]][i[0]];
          up2db[i[1]][i[0]]=tmp;
      }
    }
    
#ifndef OLD
    // extended imaging
    if (acdpars->next == 1) {
      acde_2d_4_d_b(uc2,  uc2b,  uc2d, uc2db,
                    up2,  up2b,  up2d, up2db,
                    csq2, csq2b, csq3d,
                    s, e, acdpars->z,
                    acdpars->c0,
		    acdpars->c1, acdpars->c2,
                    acdpars->lbc, acdpars->rbc,
		    acdpars->dx, acdpars->lap); 

    }
    // physical imaging
    else {
      acd_2d_d_b(uc2,  uc2b,  uc2d, uc2db,
		 up2,  up2b,  up2d, up2db,
		 csq2, csq2b, csq2d,
	         s, e, acdpars->c, acdpars->k, 
                 acdpars->lbc, acdpars->rbc,
		 &acdpars->lap[0], &acdpars->lap[acdpars->n[0]],
		 &acdpars->lap1[0], &acdpars->lap1[acdpars->n[0]]);
    }

#else // end NEW branch

    /* 2nd order case */ 
    if (acdpars->k == 1) {
      acd_2d_2_d_b(uc2,  uc2b,  uc2d, uc2db,
                   up2,  up2b,  up2d, up2db,
                   csq2, csq2b, csq2d,
                  s, e,
                  acdpars->c0, acdpars->c1);
    }
    /* 4th order case */
    else if (acdpars->k == 2) {
      if (acdpars->next == 1) {
        acde_2d_4_d_b(uc2,  uc2b,  uc2d, uc2db,
                     up2,  up2b,  up2d, up2db,
                     csq2, csq2b, csq3d,
                     s, e, acdpars->z,
                     acdpars->c0,
                     acdpars->c1, acdpars->c2,
                     acdpars->lbc, acdpars->rbc,
                     acdpars->dx, acdpars->lap);
      } else {
        acd_2d_4_d_b(uc2,  uc2b,  uc2d, uc2db,
                     up2,  up2b,  up2d, up2db,
                     csq2, csq2b, csq2d,
                     s, e,
                     acdpars->c0,
                     acdpars->c1, acdpars->c2,
                     acdpars->lbc, acdpars->rbc);
      }
    }
    /* 8th order case */
    else if (acdpars->k == 4) {
        acd_2d_8_d_b(uc2,  uc2b,  uc2d, uc2db,
                     up2,  up2b,  up2d, up2db,
                    csq2, csq2b, csq2d, 
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

/*    for (i[2]=s0[2];i[2]<=e0[2];i[2]++) {
      for (i[1]=s0[1];i[1]<=e0[1];i[1]++) {
	for (i[0]=s0[0];i[0]<=e0[0];i[0]++) {
	  tmp=((p->_s)[D_UC]._s3)[i[2]][i[1]][i[0]];
	  ((p->_s)[D_UC]._s3)[i[2]][i[1]][i[0]]=((p->_s)[D_UP]._s3)[i[2]][i[1]][i[0]];
	  ((p->_s)[D_UP]._s3)[i[2]][i[1]][i[0]]=tmp;
	}
      }
    }
*/
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
