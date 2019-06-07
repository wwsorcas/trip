#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "contextAD.h"

double nextRandom() {
  dbad_currentSeed += dbad_seed ;
  if (dbad_currentSeed>1.0) dbad_currentSeed-=1.0 ;
  return dbad_currentSeed ;
}

void context_tgt_init(double epsilon, double seed) {
  dbad_mode = 1 ;
  dbad_ddeps = epsilon ;
  dbad_seed = seed ;
  char* phase = getenv("DBAD_PHASE") ;
  if (phase==NULL) {
    printf("Please set DBAD_PHASE environment variable to 1 (perturbed) or 2 (tangent)\n") ;
    exit(0) ;
  } else if (strcmp(phase,"2")==0) {
    printf("Tangent code,  seed=%7.1e\n", seed) ;
    printf("=============================================\n") ;
    dbad_phase = 2 ;
    dbad_currentSeed = 0.0 ;
  } else if (strcmp(phase,"1")==0) {
    printf("Perturbed run, seed=%7.1e, epsilon=%7.1e\n", seed, epsilon) ;
    printf("=============================================\n") ;
    dbad_phase = 1 ;
    dbad_currentSeed = 0.0 ;
  } else if (strcmp(phase,"99")==0) {
    printf("INTERNAL INTERFACE TESTS, seed=%7.1e, epsilon=%7.1e\n", seed, epsilon) ;
    printf("=============================================\n") ;
    dbad_phase = 99 ;
  } else {
    printf("DBAD_PHASE environment variable must be set to 1 or 2\n") ;
    exit(0) ;
  }
}
/** Version of context_tgt_init called from Fortran */
void context_tgt_init_(double *epsilon, double *seed) {
  context_tgt_init(*epsilon, *seed) ;
}

void context_tgt_initreal8(char* varname, double *indep, double *indepd) {
  *indepd = nextRandom() ;
  if (dbad_phase==1)
    *indep = (*indep)+dbad_ddeps*(*indepd) ;
  else if (dbad_phase==99)
    printf("initreal8_ of %s: %24.16e //%24.16e\n", varname, *indep, *indepd) ;
}
/** Version of context_tgt_initreal8 called from Fortran */
void context_tgt_initreal8_(char* varname, double *indep, double *indepd) {
  context_tgt_initreal8(varname, indep, indepd) ;
}

void context_tgt_initreal8array(char* varname, double *indep, double *indepd, int length) {
  int i ;
  for (i=0 ; i<length ; ++i) {
    indepd[i] = nextRandom() ;
  }
  if (dbad_phase==1) {
    for (i=0 ; i<length ; ++i) {
      indep[i] = indep[i]+dbad_ddeps*indepd[i] ;
    }
  } else if (dbad_phase==99) {
    printf("initreal8array_ of %s, length=%i:\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e //%24.16e",i,indep[i],indepd[i]) ;
    printf("\n") ;
  }
}
/** Version of context_tgt_initreal8array called from Fortran */
void context_tgt_initreal8array_(char* varname, double *indep, double *indepd, int *length) {
  context_tgt_initreal8array(varname, indep, indepd, *length) ;
}

void context_tgt_concludestart() {
  dbad_currentSeed= 0.0 ;
  dbad_condensed_val = 0.0 ;
  dbad_condensed_tgt = 0.0 ;
}
/** Version of context_tgt_concludestart called from Fortran */
void context_tgt_concludestart_() {
  context_tgt_concludestart() ;
}

void context_tgt_concludereal8(char* varname, double dep, double depd) {
  double depb = nextRandom() ;
  dbad_condensed_val += depb*(dep) ;
  if (dbad_phase==2)
    dbad_condensed_tgt += depb*(depd) ;
  else if (dbad_phase==99)
    printf("concludereal8_ %24.16e //%24.16e //%24.16e\n", depb, dep, depd) ;
}
/** Version of context_tgt_concludereal8 called from Fortran */
void context_tgt_concludereal8_(char* varname, double *dep, double *depd) {
  context_tgt_concludereal8(varname, *dep, *depd) ;
}

void context_tgt_concludereal8array(char* varname, double *dep, double *depd, int length) {
  int i ;
  double depb ;
  for (i=0 ; i<length ; ++i) {
    depb = nextRandom() ;
    dbad_condensed_val += depb*dep[i] ;
    if (dbad_phase==2) {
       dbad_condensed_tgt += depb*depd[i] ;
    }
  }
}
/** Version of context_tgt_concludereal8array called from Fortran */
void context_tgt_concludereal8array_(char* varname, double *dep, double *depd, int *length) {
  context_tgt_concludereal8array(varname, dep, depd, *length) ;
}

void context_tgt_conclude() {
  if (dbad_phase==2) {
    printf("[seed:%7.1e] Condensed result : %24.16e\n", dbad_seed, dbad_condensed_val) ;
    printf("[seed:%7.1e] Condensed tangent: %24.16e\n", dbad_seed, dbad_condensed_tgt) ;
  } else if (dbad_phase==1) {
    printf("[seed:%7.1e] Condensed perturbed result : %24.16e (epsilon:%7.1e)\n",
           dbad_seed, dbad_condensed_val, dbad_ddeps) ;
  }
}
/** Version of context_tgt_conclude called from Fortran */
void context_tgt_conclude_() {
  context_tgt_conclude() ;
}

void context_adj_init(double seed) {
  dbad_mode = 0 ;
  dbad_seed = seed ;
  char* phase = getenv("DBAD_PHASE") ;
  if (phase==NULL) {
    dbad_phase = 0 ;
  } else if (strcmp(phase,"99")==0) {
    dbad_phase = 99 ;
    printf("INTERNAL INTERFACE TESTS, seed=%7.1e\n", seed) ;
  } else {
    dbad_phase = 0 ;
  }
  printf("Adjoint code,  seed=%7.1e\n", seed) ;
  printf("===================================\n") ;
  dbad_currentSeed = 0.0 ;
}
/** Version of context_adj_init called from Fortran */
void context_adj_init_(double *seed) {
  context_adj_init(*seed) ;
}

void context_adj_initreal8(char* varname, double *dep, double *depb) {
  *depb = nextRandom() ;
  if (dbad_phase==99)
    printf("initreal8_ %24.16e\n", *depb) ;
}
/** Version of context_adj_initreal8 called from Fortran */
void context_adj_initreal8_(char* varname, double *dep, double *depb) {
  context_adj_initreal8(varname, dep, depb) ;
}

void context_adj_initreal8array(char* varname, double *dep, double *depb, int length) {
  int i ;
  for (i=0 ; i<length ; ++i) {
    depb[i] = nextRandom() ;
  }
  if (dbad_phase==99) {
    printf("initreal8array_ length=%i\n", length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e",i,depb[i]) ;
    printf("\n") ;
  }
}
/** Version of context_adj_initreal8array called from Fortran */
void context_adj_initreal8array_(char* varname, double *dep, double *depb, int *length) {
  context_adj_initreal8array(varname, dep, depb, *length) ;
}

void context_adj_concludestart() {
  dbad_currentSeed= 0.0 ;
  dbad_condensed_adj = 0.0 ;
}
/** Version of context_adj_concludestart called from Fortran */
void context_adj_concludestart_() {
  context_adj_concludestart() ;
}

void context_adj_concludereal8(char* varname, double dep, double depb) {
  double depd = nextRandom() ;
  dbad_condensed_adj += depd*depb ;
  if (dbad_phase==99)
    printf("concludereal8_ %24.16e //%24.16e\n", depb, depd) ;
}
/** Version of context_adj_concludereal8 called from Fortran */
void context_adj_concludereal8_(char* varname, double *dep, double *depb) {
  context_adj_concludereal8(varname, *dep, *depb) ;
}

void context_adj_concludereal8array(char* varname, double *dep, double *depb, int length) {
  int i ;
  double depd ;
  for (i=0 ; i<length ; ++i) {
    depd = nextRandom() ;
    dbad_condensed_adj += depd*depb[i] ;
  }
}
/** Version of context_adj_concludereal8array called from Fortran */
void context_adj_concludereal8array_(char* varname, double *dep, double *depb, int *length) {
  context_adj_concludereal8array(varname, dep, depb, *length) ;
}

void context_adj_conclude() {
  printf("[seed:%7.1e] Condensed adjoint: %24.16e\n", dbad_seed, dbad_condensed_adj) ;
}
/** Version of context_adj_conclude called from Fortran */
void context_adj_conclude_() {
  context_adj_conclude() ;
}
