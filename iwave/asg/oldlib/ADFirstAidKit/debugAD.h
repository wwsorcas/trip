#ifndef DEBUGAD_INCLUDED
#define DEBUGAD_INCLUDED

void debug_tgt_init(double epsilon, double ezero, double errmax, double seed) ;
void debug_tgt_call(char *unitname, int traced, int forcetraced) ;
void debug_tgt_exit() ;
int debug_tgt_here(char* placename, int forcetraced) ;

void debug_tgt_initreal8(char* varname, double *indep, double *indepd) ;
void debug_tgt_initreal8array(char* varname, double *indep, double *indepd, int length) ;
void debug_tgt_concludereal8(char *varname, double dep, double depd) ;
void debug_tgt_concludereal8array(char *varname, double *tvar, double *tvard, int length) ;
void debug_tgt_conclude() ;
void debug_tgt_real8(char *varname, double var, double vard) ;
void debug_tgt_real8array(char *varname, double *var, double *vard, int length) ;
void debug_tgt_display(char *placename) ;

void debug_bwd_init(double ezero, double errmax, double seed) ;
void debug_bwd_call(char *funcname, int deltadepth) ;
void debug_bwd_exit() ;
int debug_bwd_here(char* placename) ;

void debug_fwd_init(double ezero, double errmax, double seed) ;
void debug_fwd_call(char *funcname) ;
void debug_fwd_exit() ;
int debug_fwd_here(char* placename) ;

void debug_adj_rwreal8(double *vard) ;
void debug_adj_rreal8(double *vard) ;
void debug_adj_wreal8(double *vard) ;
void debug_adj_rwreal8array(double *vard, int length) ;
void debug_adj_rreal8array(double *vard, int length) ;
void debug_adj_wreal8array(double *vard, int length) ;
void debug_adj_rwdisplay(char *placename, int indent) ;
void debug_adj_rdisplay(char *placename, int indent) ;
void debug_adj_wdisplay(char *placename, int indent) ;
void debug_adj_skip(char *placename) ;
void debug_adj_conclude() ;

#endif
