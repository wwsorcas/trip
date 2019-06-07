#ifndef CONTEXTAD_INCLUDED
#define CONTEXTAD_INCLUDED

static int dbad_mode, dbad_phase ;
static double dbad_ddeps = 1.e-6 ;
static double dbad_seed = 0.137 ;
static double dbad_currentSeed = 0.0 ;
static double dbad_condensed_val, dbad_condensed_tgt, dbad_condensed_adj ;

void context_tgt_init(double epsilon, double seed) ;
void context_tgt_initreal8(char* varname, double *indep, double *indepd) ;
void context_tgt_initreal8array(char* varname, double *indep, double *indepd, int length) ;
void context_tgt_concludestart() ;
void context_tgt_concludereal8(char* varname, double dep, double depd) ;
void context_tgt_concludereal8array(char* varname, double *dep, double *depd, int length) ;
void context_tgt_conclude() ;
void context_adj_init(double seed) ;
void context_adj_initreal8(char* varname, double *dep, double *depb) ;
void context_adj_initreal8array(char* varname, double *dep, double *depb, int length) ;
void context_adj_concludestart() ;
void context_adj_concludereal8(char* varname, double dep, double depb) ;
void context_adj_concludereal8array(char* varname, double *dep, double *depb, int length) ;
void context_adj_conclude() ;

#endif
