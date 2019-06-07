
#ifndef ADBUFFER_LOADED
#define ADBUFFER_LOADED 1

extern void printtraffic() ;

extern void pushinteger4(int x)  ;
extern void lookinteger4(int *x) ;
extern void popinteger4(int *x) ;

extern void pushcontrol1b(int cc) ;
extern void lookcontrol1b(int *cc) ;
extern void popcontrol1b(int *cc) ;
extern void pushcontrol2b(int cc) ;
extern void lookcontrol2b(int *cc) ;
extern void popcontrol2b(int *cc) ;
extern void pushcontrol3b(int cc) ;
extern void lookcontrol3b(int *cc) ;
extern void popcontrol3b(int *cc) ;
extern void pushcontrol4b(int cc) ;
extern void lookcontrol4b(int *cc) ;
extern void popcontrol4b(int *cc) ;
extern void pushcontrol5b(int cc) ;
extern void lookcontrol5b(int *cc) ;
extern void popcontrol5b(int *cc) ;

extern void pushreal4(float x) ;
extern void lookreal4(float *x) ;
extern void popreal4(float *x) ;

extern void pushreal8(double x) ;
extern void lookreal8(double *x) ;
extern void popreal8(double *x) ;

extern void pushcharacter(char x) ;
extern void lookcharacter(char *x) ;
extern void popcharacter(char *x) ;

extern void pushpointer4(void *x) ;
extern void lookpointer4(void **x) ;
extern void poppointer4(void **x) ;

extern void pushpointer8(void *x) ;
extern void lookpointer8(void **x) ;
extern void poppointer8(void **x) ;

/** Very complete display of the current size in bytes of
 * the global C stack followed by the auxiliary stacks. */
extern void printallbuffers() ;

extern void printbuffertop() ;
extern void showallstacks() ;

#endif
