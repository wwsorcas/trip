#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include "debugAD.h"

extern void pushreal8(double x) ;
extern void popreal8(double *x) ;
extern void pushinteger4(int x)  ;
extern void popinteger4(int *x) ;

/* The "call stack" used by debugging to
 * keep track of the position in the call tree */
typedef struct _DBAD_CallStackElem {
  char *funcname ;
  int   deltadepth ;
  int   code ;
  struct _DBAD_CallStackElem *context ;
} DBAD_CallStackElem ;

static DBAD_CallStackElem dbad_topContext ;
static DBAD_CallStackElem *dbad_callStack ;
static int dbad_calltracedepth = 1 ;

static int dbad_mode, dbad_phase, dbad_nberrors ;
static int dbad_trace = 0 ;
static int dbad_nocommunication = 0 ;
static FILE *dbad_file ;
static double dbad_almostzero, dbad_errormax, dbad_ddeps = 1.e-6 ;
static double dbad_seed = 0.137 ;
static double dbad_currentSeed = 0.0 ;
static double dbad_condensed_dd, dbad_condensed_tgt, dbad_condensed_adj ;
static double dbad_refsum, dbad_nextrefsum ;
static double dbad_coeff = 1.0 ;
static double dbad_sum = 0.0 ;

void dbad_pushCallFrame(char* unitname, int deltadepth, int forcetraced) {
  DBAD_CallStackElem *newCallLevel = (DBAD_CallStackElem*)malloc(sizeof(DBAD_CallStackElem)) ;
  newCallLevel->funcname = (char*)malloc(100) ;
  sprintf(newCallLevel->funcname, "%s", unitname) ;
  newCallLevel->deltadepth = (dbad_calltracedepth>0?1-deltadepth:0) ;
  dbad_calltracedepth -= newCallLevel->deltadepth ;
  // forcing mechanism:
  if (forcetraced>0 && forcetraced>dbad_calltracedepth) {
    newCallLevel->deltadepth -= (forcetraced-dbad_calltracedepth) ;
    dbad_calltracedepth = forcetraced ;
  }
  newCallLevel->code = 0 ;
  newCallLevel->context = dbad_callStack ;
  dbad_callStack = newCallLevel ;
}

void dbad_popCallFrame() {
  dbad_calltracedepth += dbad_callStack->deltadepth ;
  DBAD_CallStackElem *newCallLevel = dbad_callStack->context ;
  free(dbad_callStack->funcname) ;
  free(dbad_callStack) ;
  dbad_callStack = newCallLevel ;
}

int dbad_debughere(int forcetraced) {
  return (dbad_calltracedepth>0 || forcetraced) ;
}

int dbad_debugabove() {
  return (dbad_calltracedepth+dbad_callStack->deltadepth)>0 ;
}

int dbad_callstacksize() {
  DBAD_CallStackElem *incallstack = dbad_callStack ;
  int depth = 1 ;
  while (incallstack) {++depth ; incallstack=incallstack->context ;}
  return depth ;
}

void dbad_resetCondensors() {
  dbad_currentSeed = 0.0 ;
  dbad_condensed_dd = 0.0 ;
  dbad_condensed_tgt = 0.0 ;
  dbad_condensed_adj = 0.0 ;
}

double dbad_nextRandom() {
  dbad_currentSeed += dbad_seed ;
  if (dbad_currentSeed>1.0) dbad_currentSeed-=1.0 ;
  return dbad_currentSeed ;
}

void dbad_ddcheckvarname(char* varname) {
  char ddvarname[40] ;
  int ret = fscanf(dbad_file, " %s\n", ddvarname) ;
  if (strcmp(varname, ddvarname)!=0) {
    printf("Control mismatch, expecting a variable named \"%s\", got \"%s\"\n",varname,ddvarname) ;
    exit(0) ;
  }
}

//TODO: detect NaNs and if so, set *areNaNs :
void dbad_ddpicktwo8(double var, double *varwr, FILE *dbad_file, double *ddvar, int *areNaNs) {
  FILE *wrfile ;
  int ret ;
  wrfile = fopen("ddwrfile", "w") ;
  fprintf(wrfile, "%24.16e", var) ;
  freopen("ddwrfile", "r", wrfile) ;
  ret = fscanf(wrfile," %lf",varwr) ;
  fclose(wrfile) ;
  ret = fscanf(dbad_file," %lf\n", ddvar) ;
}

void dbad_display_location(char *placename) {
  int i ;
  char* enclosproc = dbad_callStack->funcname ;
  for (i=dbad_callstacksize() ; i>0 ; --i) printf("  ") ;
  printf(" AT:%s OF %s\n", placename, enclosproc) ;
}

//###################### DEBUG OF TANGENT ##############################

void debug_tgt_init(double epsilon, double ezero, double errmax, double seed) {
  dbad_mode = 1 ;
  dbad_ddeps = epsilon ;
  dbad_almostzero = ezero ;
  dbad_errormax = errmax ;
  dbad_seed = seed ;
  dbad_topContext.funcname = "TESTED CODE\0" ;
  dbad_topContext.deltadepth = 0 ;
  dbad_topContext.code = 0 ;
  dbad_calltracedepth = 1 ;
  dbad_callStack = &dbad_topContext ;
  char* phase = getenv("DBAD_PHASE") ;
  if (phase==NULL) {
    printf("Please set DBAD_PHASE environment variable to 1 (perturbed), 2 (tangent), or 0 (no debug)\n") ;
    exit(0) ;
  } else if (strcmp(phase,"0")==0) {
    dbad_phase = 0 ;
  } else if (strcmp(phase,"1")==0) {
    dbad_phase = 1 ;
  } else if (strcmp(phase,"2")==0) {
    dbad_phase = 2 ;
  } else if (strcmp(phase,"-1")==0) {
    dbad_phase = 1 ;
    dbad_trace = 1 ;
  } else if (strcmp(phase,"-2")==0) {
    dbad_phase = 2 ;
    dbad_trace = 1 ;
  } else {
    printf("DBAD_PHASE environment variable must be set to 1 or 2 (or 0 for no debug)\n") ;
    exit(0) ;
  }
  if (dbad_trace) {
    printf("INTERNAL TESTS, epsilon=%7.1e, zero=%7.1e, errmax=%4.1f%\n", epsilon, ezero, errmax) ;
    printf("===========================================================\n") ;
  }
  if (dbad_phase==1) {
    printf("Starting TGT test, phase one, epsilon=%7.1e, zero=%7.1e, errmax=%4.1f%, seed=%7.1e\n",
           epsilon, ezero, errmax, seed) ;
    printf("===========================================================\n") ;
    dbad_resetCondensors() ;
    mkfifo("/tmp/DBAD_fifo", S_IWUSR | S_IRUSR | S_IRGRP | S_IROTH | S_IRWXU | S_IRWXO) ;
    dbad_file = fopen("/tmp/DBAD_fifo", "a") ;
    if (dbad_file==NULL) {
      char errbuf[20] ;
      strerror_r(errno, errbuf, 20) ;
      printf("FIFO ERROR %i: %s  OR  %s\n",errno,strerror(errno),errbuf) ;
      exit(0) ;
    }
  } else if (dbad_phase==2) {
    printf("Starting TGT test, phase two, epsilon=%7.1e, zero=%7.1e, errmax=%4.1f%, seed=%7.1e\n",
           epsilon, ezero, errmax, seed) ;
    printf("===========================================================\n") ;
    dbad_resetCondensors() ;
    dbad_file = fopen("/tmp/DBAD_fifo", "r") ;
  }
}

void debug_tgt_call(char* unitname, int deltadepth, int forcetraced) {
  if (dbad_trace) {
    printf("call_ %s deltadepth:%i forcetraced:%i\n", unitname, deltadepth, forcetraced) ;
  }
  if (dbad_phase!=0) {
    dbad_pushCallFrame(unitname, deltadepth, forcetraced) ;
  }
}

void debug_tgt_exit() {
  if (dbad_trace) {
    printf("exit_\n") ;
  }
  if (dbad_phase!=0) {
    dbad_popCallFrame() ;
  }
  dbad_resetCondensors() ;
}

int debug_tgt_here(char* placename, int forcetraced) {
  if (dbad_trace) {
    printf("here_?? %s forcetraced:%i\n", placename, forcetraced) ;
    printf("    will return %i\n", (dbad_phase==0?0:dbad_debughere(forcetraced))) ;
  }
  if (dbad_phase==0) return 0 ;
  return dbad_debughere(forcetraced) ;
}

void debug_tgt_initreal8(char* varname, double *indep, double *indepd) {
  *indepd = dbad_nextRandom() ;
  if (dbad_phase==1)
    *indep = (*indep)+dbad_ddeps*(*indepd) ;
  if (dbad_trace)
    printf("initreal8_ of %s: %24.16e //%24.16e\n", varname, *indep, *indepd) ;
}

void debug_tgt_initreal8array(char* varname, double *indep, double *indepd, int length) {
  int i ;
  for (i=0 ; i<length ; ++i) {
    indepd[i] = dbad_nextRandom() ;
  }
  if (dbad_phase==1) {
    for (i=0 ; i<length ; ++i) {
      indep[i] = indep[i]+dbad_ddeps*indepd[i] ;
    }
  }
  if (dbad_trace) {
    printf("initreal8array_ of %s, length=%i:\n", varname, length) ;
    for (i=0 ; i<length ; ++i)
      printf("    %i:%24.16e //%24.16e",i,indep[i],indepd[i]) ;
    printf("\n") ;
  }
}

void debug_tgt_concludereal8(char* varname, double dep, double depd) {
  double depb = dbad_nextRandom() ;
  if (dbad_trace) {
    printf("concludereal8_ %s %24.16e //%24.16e\n", varname, dep, depd) ;
  }
  if (dbad_phase==1) {
    fprintf(dbad_file, "%s\n", varname) ;
    fprintf(dbad_file, "%24.16e\n", dep) ;
  } else if (dbad_phase==2) {
    double ddvar, dd, diff, varwr, maxabs, absvard, absdd, absvardmdd ;
    int areNaNs = 0;
    dbad_ddcheckvarname(varname) ;
    dbad_ddpicktwo8(dep, &varwr, dbad_file, &ddvar, &areNaNs) ;
    if (!areNaNs) {
      dd = (ddvar-varwr)/dbad_ddeps ;
      absvard = (depd>=0.0?depd:-depd) ;
      absdd = (dd>=0.0?dd:-dd) ;
      maxabs = (absvard>absdd?absvard:absdd) ;
      absvardmdd = depd-dd ;
      if (absvardmdd<0.0) absvardmdd=-absvardmdd ;
      diff = (absvardmdd*100.0)/ maxabs ;
      if (diff>dbad_errormax) {
        printf("Result variable %s sum: %24.16e (ad)%5.2f% DIFF WITH (dd)%24.16e\n",
               varname,depd,diff,dd) ;
      }
      dbad_condensed_dd += dd*depb;
      dbad_condensed_tgt += depd*depb;
    } else {
      printf("Variable %s has NaNs\n", varname) ;
    }
  }
}

void debug_tgt_concludereal8array(char* varname, double *dep, double *depd, int length) {
  if (dbad_trace) {
    printf("concludereal8array_ %s[%i]:", varname, length) ;
  }
  if (dbad_phase==1) {
    fprintf(dbad_file, "%s\n", varname) ;
  } else if (dbad_phase==2) {
    dbad_ddcheckvarname(varname) ;
  }
  int i ;
  double depb ;
  for (i=0 ; i<length ; ++i) {
    depb = dbad_nextRandom() ;
    if (dbad_trace) {
      printf("    %i:%24.16e //%24.16e",i,dep[i],depd[i]) ;
    }
    if (dbad_phase==1) {
      fprintf(dbad_file, "%24.16e\n", dep[i]) ;
    } else if (dbad_phase==2) {
      double ddvar, dd, diff, varwr, maxabs, absvard, absdd, absvardmdd ;
      int areNaNs = 0;
      dbad_ddpicktwo8(dep[i], &varwr, dbad_file, &ddvar, &areNaNs) ;
      if (!areNaNs) {
        dd = (ddvar-varwr)/dbad_ddeps ;
        absvard = (depd[i]>=0.0?depd[i]:-depd[i]) ;
        absdd = (dd>=0.0?dd:-dd) ;
        maxabs = (absvard>absdd?absvard:absdd) ;
        absvardmdd = depd[i]-dd ;
        if (absvardmdd<0.0) absvardmdd=-absvardmdd ;
        diff = (absvardmdd*100.0)/ maxabs ;
        if (diff>dbad_errormax) {
          printf("Result variable %s[%i] sum: %24.16e (ad)%5.2f% DIFF WITH (dd)%24.16e\n",
                 varname,i,depd[i],diff,dd) ;
        }
        dbad_condensed_dd += dd*depb;
        dbad_condensed_tgt += depd[i]*depb;
      } else {
        printf("Variable %s[%i] has NaNs\n", varname, i) ;
      }
    }
  }
  if (dbad_trace) {
    printf("\n") ;
  }
}

void debug_tgt_conclude() {
  printf("===========================================================\n") ;
  if (dbad_trace) {
    printf("    condensed result %24.16e //%24.16e",dbad_condensed_tgt,dbad_condensed_dd) ;
  }
  if (dbad_phase==2) {
    double abstgt, absdd, maxabs, abserror, diff ;
    abstgt = (dbad_condensed_tgt>=0.0?dbad_condensed_tgt:-dbad_condensed_tgt) ;
    absdd = (dbad_condensed_dd>=0.0?dbad_condensed_dd:-dbad_condensed_dd) ;
    maxabs = (abstgt>absdd?abstgt:absdd) ;
    abserror = dbad_condensed_tgt-dbad_condensed_dd ;
    if (abserror<0.0) abserror=-abserror ;
    diff = (abserror*100.0)/ maxabs ;
    printf("[seed:%7.1e] Condensed tangent : %24.16e (ad)%5.2f% DIFF WITH (dd)%24.16e\n",
           dbad_seed,dbad_condensed_tgt,diff,dbad_condensed_dd) ;
  }
}

void debug_tgt_real8(char *varname, double var, double vard) {
  if (dbad_phase==1) {
    fprintf(dbad_file, "%s\n", varname) ;
    fprintf(dbad_file, "%24.16e\n", var) ;
  } else if (dbad_phase==2) {
    double ddvar, dd, diff, varwr, maxabs, absvard, absdd, absvardmdd ;
    int areNaNs = 0 ;
    dbad_ddcheckvarname(varname) ;
    dbad_ddpicktwo8(var, &varwr, dbad_file, &ddvar, &areNaNs) ;
    if (!areNaNs || dbad_trace) {
      dd = (ddvar-varwr)/dbad_ddeps ;
      absvard = (vard>=0.0?vard:-vard) ;
      absdd = (dd>=0.0?dd:-dd) ;
      if (absvard>dbad_almostzero || absdd>dbad_almostzero || dbad_trace) {
        maxabs = (absvard>absdd?absvard:absdd) ;
        absvardmdd = vard-dd ;
        if (absvardmdd<0.0) absvardmdd=-absvardmdd ;
        diff = (absvardmdd*100.0)/maxabs ;
        if (dbad_trace) {
          printf("%s v-eps:%24.16e v-loc:%24.16e ->(dd)%24.16e=?=(ad)%24.16e %5.1f%\n",
                 varname,ddvar,varwr,dd,vard,diff) ;
        } else if (diff>dbad_errormax) {
          printf("%s:  %24.16e (ad) %5.1f% DIFF WITH (dd) %24.16e\n",
                 varname, vard, diff, dd) ;
        }
      }
    }
  }
}

void debug_tgt_real8array(char *varname, double* var, double* vard, int length) {
  int i ;
  if (dbad_phase==1) {
    fprintf(dbad_file, "%s\n", varname) ;
    for (i=0 ; i<length ; ++i) {
      fprintf(dbad_file, "%24.16e\n", var[i]) ;
    }
  } else if (dbad_phase==2) {
    double ddvar, dd, diff, varwr, maxabs, absvard, absdd, absvardmdd ;
    int areNaNs = 0 ;
    double vardbuf[10], ddbuf[10] ;
    double varepsbuf[10], varlocbuf[10] ;
    int indexbuf[10] ;
    dbad_ddcheckvarname(varname) ;
    int ibuf = 0;
    int printedheader = 0 ;
    for (i=0 ; i<length ; ++i) {
      dbad_ddpicktwo8(var[i],&varwr,dbad_file,&ddvar,&areNaNs) ;
      if (!areNaNs || dbad_trace) {
        dd = (ddvar-varwr)/dbad_ddeps ;
        absvard = (vard[i]>=0.0?vard[i]:-vard[i]) ;
        absdd = (dd>=0.0?dd:-dd) ;
        if (absvard>dbad_almostzero || absdd>dbad_almostzero || dbad_trace) {
          maxabs = (absvard>absdd?absvard:absdd) ;
          absvardmdd = vard[i]-dd ;
          if (absvardmdd<0.0) absvardmdd=-absvardmdd ;
          diff = (absvardmdd*100.0)/maxabs ;
          if (dbad_trace) {
            vardbuf[ibuf] = vard[i] ;
            varepsbuf[ibuf] = ddvar ;
            varlocbuf[ibuf] = varwr ;
            ddbuf[ibuf] = dd ;
            indexbuf[ibuf] = i ;
            ++ibuf ;
          } else if (diff>dbad_errormax) {
            vardbuf[ibuf] = vard[i] ;
            ddbuf[ibuf] = dd ;
            indexbuf[ibuf] = i ;
            ++ibuf ;
          }
        }
      }
      if (ibuf>=10 || (i==length-1 && ibuf>0)) {
        int j ;
        if (!printedheader) {
          printf("%s:\n", varname) ;
          printedheader = 1 ;
        }
        printf("    ") ;
        for (j=0 ; j<ibuf ; ++j)
          printf(" %4i->%11.4e", indexbuf[j], vardbuf[j]) ;
        printf("\n    ") ;
        for (j=0 ; j<ibuf ; ++j)
          printf("  (eps)%12.5e", varepsbuf[j]) ;
        printf("\n    ") ;
        for (j=0 ; j<ibuf ; ++j)
          printf("  (loc)%12.5e", varlocbuf[j]) ;
        printf("\n    ") ;
        for (j=0 ; j<ibuf ; ++j)
          printf("  (dd:)%11.4e", ddbuf[j]) ;
        printf("\n") ;
        ibuf = 0 ;
      }
    }
  }
}

void debug_tgt_display(char *placename) {
  if (dbad_trace) {
    printf("display_ %s\n", placename) ;
  }
  if (dbad_phase==2) {
    dbad_display_location(placename) ;
  }
}

//############## DEBUG OF ADJOINT, FIRST SWEEP: ADJOINT RUN ################

void debug_bwd_init(double ezero, double errmax, double seed) {
  dbad_mode = -1 ;
  dbad_phase = 1 ;
  dbad_almostzero = ezero ;
  dbad_errormax = errmax ;
  dbad_seed = seed ;
  dbad_topContext.funcname =  "TESTED CODE\0" ;
  dbad_topContext.deltadepth = 0 ;
  dbad_topContext.code = 0 ;
  dbad_calltracedepth = 1 ;
  dbad_callStack = &dbad_topContext ;
  char* phase = getenv("DBAD_PHASE") ;
  if (phase==NULL) {
    printf("Please set DBAD_PHASE environment variable to 0 (no debug), 1 (sendToTgt), -1 (plusTraces)\n") ;
    dbad_phase = 1 ;
  } else if (strcmp(phase,"0")==0) {
    dbad_phase = 1 ;
    dbad_nocommunication = 1 ;
  } else if (strcmp(phase,"1")==0) {
    dbad_phase = 1 ;
  } else if (strcmp(phase,"-1")==0) {
    dbad_phase = 1 ;
    dbad_trace = 1 ;
  } else {
    printf("DBAD_PHASE environment variable must be set to 1 or -1 (or 0 for no debug)\n") ;
    exit(0) ;
  }
  printf("Starting ADJ test, phase one (bwd), zero=%7.1e, errmax=%4.1f%, seed=%7.1e\n", ezero, errmax, seed) ;
  printf("===========================================================\n") ;
  if (dbad_nocommunication) {
    dbad_file = NULL ;
    printf("FIFO COMMUNICATION TURNED OFF !\n") ;
  } else {
    mkfifo("/tmp/DBAD_fifo", S_IWUSR | S_IRUSR | S_IRGRP | S_IROTH | S_IRWXU | S_IRWXO) ;
    dbad_file = fopen("/tmp/DBAD_fifo", "a") ;
    if (dbad_file==NULL) {
      char errbuf[20] ;
      strerror_r(errno, errbuf, 20) ;
      printf("FIFO ERROR %i: %s  OR  %s\n",errno,strerror(errno),errbuf) ;
      exit(0) ;
    }
  }
  dbad_resetCondensors() ;
}

void debug_bwd_call(char *funcname, int deltadepth) {
  dbad_pushCallFrame(funcname, deltadepth, 0) ;
}

void debug_bwd_exit() {
  if (dbad_debugabove()) {
    if (dbad_nocommunication) {
      printf("debugAD would send (%i %s)\n", (dbad_debughere(0)?2:-2), dbad_callStack->funcname) ;
    } else {
      fprintf(dbad_file, "%i\n", (dbad_debughere(0)?2:-2)) ;
      fprintf(dbad_file, "%s\n", dbad_callStack->funcname) ;
    }
  }
  dbad_popCallFrame() ;
}

int debug_bwd_here(char* placename) {
  dbad_resetCondensors() ;
  return debug_tgt_here(placename, 0) ;
}

//############## DEBUG OF ADJOINT, SECOND SWEEP: TANGENT RUN ################

void debug_fwd_init(double ezero, double errmax, double seed) {
  dbad_mode = -1 ;
  dbad_phase = 2 ;
  dbad_almostzero = ezero ;
  dbad_errormax = errmax ;
  dbad_seed = seed ;
  dbad_topContext.funcname = "TESTED CODE\0" ;
  dbad_topContext.deltadepth = 0 ;
  dbad_topContext.code = 0 ;
  dbad_calltracedepth = 1 ;
  dbad_callStack = &dbad_topContext ;
  char* phase = getenv("DBAD_PHASE") ;
  if (phase==NULL) {
    printf("Please set DBAD_PHASE environment variable to 0 (no debug), 2 (readFromAdj), -2 (plusTraces)\n") ;
    dbad_phase = 2 ;
  } else if (strcmp(phase,"0")==0) {
    dbad_phase = 2 ;
    dbad_nocommunication = 1 ;
  } else if (strcmp(phase,"2")==0) {
    dbad_phase = 2 ;
  } else if (strcmp(phase,"-2")==0) {
    dbad_phase = 2 ;
    dbad_trace = 1 ;
  } else {
    printf("DBAD_PHASE environment variable must be set to 2 or -2 (or 0 for no debug)\n") ;
    exit(0) ;
  }
  dbad_nberrors = 0 ;
  printf("Starting ADJ test, phase two (fwd), zero=%7.1e, errmax=%4.1f%, seed=%7.1e\n", ezero, errmax, seed) ;
  printf("===========================================================\n") ;
  if (dbad_nocommunication) {
    dbad_file = NULL ;
    printf("FIFO COMMUNICATION TURNED OFF !\n") ;
  } else {
    dbad_file = fopen("/tmp/DBAD_fifo", "r") ;
    dbad_resetCondensors() ;
    /* Convention on the meaning of labels:
       -1 -> a debug point, skipped
        0 -> a debug point, traced but no associated value.
        1 -> a debug point, traced, with an associated value.
       -2 -> a call, skipped
        2 -> a call, traced
    */
    int ret=0 ;
    int label ;
    char placename[40] ;
    double value ;
    while (1) {
      ret = fscanf(dbad_file, "%i\n", &label) ;
      if (ret!=1) break ;
      ret = fscanf(dbad_file, "%s\n", placename) ;
      if (label==1) {
        ret = fscanf(dbad_file, "%lf\n", &value) ;
        pushreal8(value) ;
      }
      pushcharacterarray(placename, 40) ;
      pushinteger4(label) ;
    }
  }
}

void debug_fwd_call(char *funcname) {
  int label ;
  char funcnamefrom[40] ;
  char funcnamehere[40] ;
  // In the special debug^2 case, on the 2nd phase (tangent) of the debugAdj, with DBAD_PHASE=0,
  // push the call frame but do essentially nothing!
  if (dbad_debughere(0) && !(dbad_nocommunication && dbad_phase==2)) {
    popinteger4(&label) ;
    if (label!=2 && label!=-2) {
      printf("Control mismatch, expecting a trace (-2or2) from %s bwd call exit, got %i\n",funcname,label) ;
      exit(0) ;
    }
    popcharacterarray(funcnamefrom, 40) ;
    sprintf(funcnamehere,"%s",funcname) ;
    if (strcmp(funcnamefrom,funcnamehere)!=0) {
      printf("Control mismatch, expecting a call to %s, got %s\n",funcnamehere,funcnamefrom) ;
      exit(0) ;
    }
    dbad_pushCallFrame(funcname, 0, 0) ;
    if (label==2) { // then the call is traced:
      dbad_callStack->deltadepth += (dbad_calltracedepth-1) ;
      dbad_calltracedepth = 1 ;
    } else { // then label==-2: the call is not traced:
      dbad_callStack->deltadepth += dbad_calltracedepth ;
      dbad_calltracedepth = 0 ;
    }
  } else {
    dbad_pushCallFrame(funcname, 0, 0) ;
  }
}

void debug_fwd_exit() {
  dbad_popCallFrame() ;
}

int debug_fwd_here(char* placename) {
  // In the special debug^2 case, on the 2nd phase (tangent) of the debugAdj, with DBAD_PHASE=0,
  // never go into the derivative manipulation body, except to st the inputs at the very "start"
  // and to print the result at the very "end".
  if (dbad_nocommunication && dbad_phase==2) {
    if (strcmp(placename,"end")==0 || strcmp(placename,"start")==0)
      return 1 ;
    else
      return 0 ;
  } else {
    if (dbad_debughere(0)) {
      int label ;
      char placenamefrom[40] ;
      char placenamehere[40] ;
      dbad_resetCondensors() ;
      popinteger4(&label) ;
      if (label!=1 && label!=-1 && label!=0) {
        printf("Control mismatch, expecting a trace (-1or0or1) from place %s, got %i\n",placename,label) ;
        exit(0) ;
      }
      popcharacterarray(placenamefrom, 40) ;
      sprintf(placenamehere,"%s",placename) ;
      if (strcmp(placenamefrom,placenamehere)!=0) {
        printf("Control mismatch, expecting place %s, got %s\n",placenamehere,placenamefrom) ;
        exit(0) ;
      }
      if (label==1) {
        popreal8(&dbad_nextrefsum) ;
      }
      return label!=-1 ;
    } else {
      return 0 ;
    }
  }
}

//############## DEBUG OF ADJOINT, FOR BOTH SWEEPS: ################

void debug_adj_rwreal8(double *vard) {
  double varb = dbad_nextRandom() ;
  dbad_condensed_adj += varb*(*vard) ;
  *vard = varb ;
}

/** Although at present this routine doesn't modify its argument,
 * we still expect a reference, just in case, and for consistency
 * with debug_adj_wreal8() */
void debug_adj_rreal8(double *vard) {
  double varb = dbad_nextRandom() ;
  dbad_condensed_adj += varb*(*vard) ;
}

void debug_adj_wreal8(double *vard) {
  *vard = dbad_nextRandom() ;
}

void debug_adj_rwreal8array(double *vard, int length) {
  int i ;
  if (vard)
    for (i=0 ; i<length ; ++i)
      debug_adj_rwreal8(&(vard[i])) ;
}

void debug_adj_rreal8array(double *vard, int length) {
  int i ;
  if (vard)
    for (i=0 ; i<length ; ++i)
      debug_adj_rreal8(&(vard[i])) ;
}

void debug_adj_wreal8array(double *vard, int length) {
  int i ;
  if (vard)
    for (i=0 ; i<length ; ++i)
      debug_adj_wreal8(&(vard[i])) ;
}

void debug_adj_rwdisplay(char *placename, int indent) {
  debug_adj_rdisplay(placename, indent) ;
  if (dbad_phase==2)
    dbad_refsum = dbad_nextrefsum ;
}

void debug_adj_rdisplay(char *placename, int indent) {
  if (dbad_phase==1) {
    if (dbad_nocommunication) {
      printf("debugAD would send (1 %s %24.16e)\n", placename, dbad_condensed_adj) ;
    } else {
      fprintf(dbad_file, "1\n") ;
      fprintf(dbad_file, "%s\n", placename) ;
      fprintf(dbad_file, "%24.16e\n", dbad_condensed_adj) ;
    }
  } else if (dbad_phase==2) {
    // In the special debug^2 case, on the 2nd phase (tangent) of the debugAdj, with DBAD_PHASE=0,
    // debug_adj_wdisplay is called ony on the "end" location. Print the tangent result:
    if (dbad_nocommunication) {
      printf("Condensed tangent result is %24.16e\n", dbad_condensed_adj) ;
    } else {
      double absref = (dbad_refsum>=0.0?dbad_refsum:-dbad_refsum) ;
      double absadj = (dbad_condensed_adj>=0.0?dbad_condensed_adj:-dbad_condensed_adj) ;
      double absdiff = dbad_refsum - dbad_condensed_adj ;
      if (absdiff<0.0) absdiff = -absdiff ;
      double reldiff = (absdiff*200.0)/(absref+absadj) ;
      if (reldiff>dbad_errormax) {
        printf("                         %5.1f% DIFFERENCE!!  tgt:%24.16e  adj:%24.16e\n",
               reldiff, dbad_condensed_adj, dbad_refsum) ;
        ++dbad_nberrors ;
      } else if (strcmp(placename,"end")==0 && dbad_nberrors==0) {
        // When we are at end and no errors were found, always show the compared values
        printf("                         difference is just %7.3f% between tgt:%24.16e and adj:%24.16e\n",
               reldiff, dbad_condensed_adj, dbad_refsum) ;
      }
      if (indent==0) dbad_display_location(placename) ;
    }
  }
  dbad_resetCondensors() ;
}

void debug_adj_wdisplay(char *placename, int indent) {
  if (dbad_phase==1) {
    if (dbad_nocommunication) {
      printf("debugAD would send (0 %s)\n", placename) ;
    } else {
      fprintf(dbad_file, "0\n") ;
      fprintf(dbad_file, "%s\n", placename) ;
    }
  } else if (dbad_phase==2) {
      if (indent==0) dbad_display_location(placename) ;
      dbad_refsum = dbad_nextrefsum ;
  }
  dbad_resetCondensors() ;
}

void debug_adj_skip(char *placename) {
  if (dbad_phase==1 && dbad_debughere(0)) {
    if (dbad_nocommunication) {
      printf("debugAD would send (-1 %s)\n", placename) ;
    } else {
      fprintf(dbad_file, "-1\n") ;
      fprintf(dbad_file, "%s\n", placename) ;
    }
  }
}

void debug_adj_conclude() {
  if (dbad_phase==2) {
    if (!dbad_nocommunication) {
      // In the special debug^2 case, on the 2nd phase (tangent) of the debugAdj, with DBAD_PHASE=0,
      // don't claim that any testing has been done!! but show the expected condensed tangent:
      printf("End of ADJ test, %i error(s) found.\n", dbad_nberrors) ;
      printf("===========================================================\n") ;
    }
  }
}

/* void debug_adj_show() { */
/*   printf("Present sum %24.16e, current seed is %f (%f)\n", dbad_condensed_adj, dbad_currentSeed, dbad_seed) ; */
/* } */

//############## INTERFACE PROCEDURES CALLED FROM FORTRAN ################

void debug_tgt_init_(double *epsilon, double *ezero, double *errmax, double *seed) {
  debug_tgt_init(*epsilon, *ezero, *errmax, *seed) ;
}

void debug_tgt_call_(char* unitname, int *deltadepth, int *forcetraced) {
  debug_tgt_call(unitname, *deltadepth, *forcetraced) ;
}

void debug_tgt_exit_() {
  debug_tgt_exit() ;
}

int debug_tgt_here_(char* placename, int *forcetraced) {
  return debug_tgt_here(placename, *forcetraced) ;
}

void debug_tgt_initreal8_(char* varname, double *indep, double *indepd) {
  debug_tgt_initreal8(varname, indep, indepd) ;
}

void debug_tgt_initreal8array_(char* varname, double *indep, double *indepd, int *length) {
  debug_tgt_initreal8array(varname, indep, indepd, *length) ;
}

void debug_tgt_concludereal8_(char* varname, double *dep, double *depd) {
  debug_tgt_concludereal8(varname, *dep, *depd) ;
}

void debug_tgt_concludereal8array_(char* varname, double *dep, double *depd, int *length) {
  debug_tgt_concludereal8array(varname, dep, depd, *length) ;
}

void debug_tgt_conclude_() {
  debug_tgt_conclude() ;
}

void debug_tgt_real8_(char *varname, double *var, double *vard) {
  debug_tgt_real8(varname, *var, *vard) ;
}

void debug_tgt_real8array_(char *varname, double* var, double* vard, int *length) {
  debug_tgt_real8array(varname, var, vard, *length) ;
}

void debug_tgt_display_(char *placename) {
  debug_tgt_display(placename) ;
}

void debug_bwd_init_(double *ezero, double *errmax, double *seed) {
  debug_bwd_init(*ezero, *errmax, *seed) ;
}

void debug_bwd_call_(char *funcname, int *deltadepth) {
  debug_bwd_call(funcname, *deltadepth) ;
}

void debug_bwd_exit_() {
  debug_bwd_exit() ;
}

int debug_bwd_here_(char* placename) {
  return debug_bwd_here(placename) ;
}

void debug_fwd_init_(double *ezero, double *errmax, double *seed) {
  debug_fwd_init(*ezero, *errmax, *seed) ;
}

void debug_fwd_call_(char *funcname) {
  debug_fwd_call(funcname) ;
}

void debug_fwd_exit_() {
  debug_fwd_exit() ;
}

int debug_fwd_here_(char* placename) {
  return debug_fwd_here(placename) ;
}

void debug_adj_rwreal8_(double *vard) {
  debug_adj_rwreal8(vard) ;
}

void debug_adj_rreal8_(double *vard) {
  debug_adj_rreal8(vard) ;
}

void debug_adj_wreal8_(double *vard) {
  debug_adj_wreal8(vard) ;
}

void debug_adj_rwreal8array_(double *vard, int *length) {
  debug_adj_rwreal8array(vard, *length) ;
}

void debug_adj_rreal8array_(double *vard, int *length) {
  debug_adj_rreal8array(vard, *length) ;
}

void debug_adj_wreal8array_(double *vard, int *length) {
  debug_adj_wreal8array(vard, *length) ;
}

void debug_adj_rwdisplay_(char *placename, int *indent) {
  debug_adj_rwdisplay(placename, *indent) ;
}

void debug_adj_rdisplay_(char *placename, int *indent) {
  debug_adj_rdisplay(placename, *indent) ;
}

void debug_adj_wdisplay_(char *placename, int *indent) {
  debug_adj_wdisplay(placename, *indent) ;
}

void debug_adj_skip_(char* placename) {
  debug_adj_skip(placename) ;
}

void debug_adj_conclude_() {
  debug_adj_conclude() ;
}
