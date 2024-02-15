#include "parser.h"
#ifdef IWAVE_USE_MPI
#include "mpisegypp.hh"
#else
#include "segypp.hh"
#endif
#include "segyops.hh"

#ifdef IWAVE_USE_MPI
using TSOpt::MPISEGYSpace;
typedef TSOpt::MPISEGYSpace gsp;
#else
using TSOpt::SEGYSpace;
typedef TSOpt::SEGYSpace gsp;
#endif

int xargc;
char **xargv;

const char * sdoc[] = {
  "SEGYConv: convolve or cross-correlate (adjoint-convolve) a kernel gather ( = collection",
  "of traces) with an input gather to produce an output gather. The input and output gathers",
  "should have the same geometry as the kernel, although the command will run to completion if",
  "both have the same time step and AT LEAST as many traces as the kernel gather."
  " ",
  "Can be parallelized over blocks of traces by looping over min and max trace indices."
  " ",
  "required inputs:",
  "in       = name of su file containing inpur traces:",
  "out      = name of su file containing output traces",
  "ker      = name of file containing signature trace", 
  "mode     = fwd -   convolve",
  "           adj -   adjoint-convolve (cross-correlate)",
  "min      = 0 - min tracl",
  "max      = std::numeric_limits<int>::max() - max tracl",
  NULL};

int main(int argc, char ** argv) {

  /* ARGUMENTS */
  char * inp;      /* input file name */
  char * out;      /* output file name */
  char * ker;      /* kernel file name */
  int adj;         /* 0 = fwd, 1 = adj */
  int mnx;         /* min tracl - default = 0 */
  int mxx;         /* max tracl - default = max int */

  /* INTERNAL VARIABLES */
  PARARRAY * par;  /* param array */
  FILE * fpinp;    /* input file pointer */
  FILE * fpout;    /* output file pointer */
  FILE * fpker;    /* kernel file pointer */

  xargc=argc; xargv=argv;
  requestdoc(1);

  try {

#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);   
    storeGlobalComm(MPI_COMM_WORLD);
#endif

    PARARRAY * pars = ps_new();
    
    if (ps_createargs(pars,argc-1,&(argv[1]))) {
      RVL::RVLException e;
      e<<"ERROR: SEGYCGDecon from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++) 
	e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }

    /* extract input parameters */
    par=ps_new();

    if ( ps_createargs(par, argc - 1, argv + 1) ) {
      printf("Error parsing input data. ABORT.\n");
      exit(1);
    }

    if (ps_flcstring(*par,"in",&inp)) {
      printf("error reading input file name ABORT.\n");
      exit(1);
    }

    if (ps_flcstring(*par,"out",&out)) {
      printf("error reading output file name. ABORT.\n");
      exit(1);
    }

    if (ps_flcstring(*par,"ker",&ker)) {
      printf("error reading output file name. ABORT.\n");
      exit(1);
    }

    adj = RVL::valparse<int>(*par,"adj",0);
    mnx = RVL::valparse<int>(*par,"min",0);
    mxx = RVL::valparse<int>(*par,"max",std::numeric_limits<int>::max());

    /*
    cerr<<"adj="<<adj<<endl;
    cerr<<"inp="<<inp<<endl;
    cerr<<"out="<<out<<endl;
    cerr<<"ker="<<ker<<endl;
    cerr<<"mnx="<<mnx<<endl;
    cerr<<"mxx="<<mxx<<endl;
    */
    
    /* open data files */
    if (!(fpinp=fopen(inp,"r"))) {
      printf("SEGYConv: failed to open input file = %s for read. ABORT.\n",inp);
      exit(1);
    }

    if (!(fpout=fopen(out,"r+"))) {
      printf("SEGYConv: failed to open output file = %s for read. ABORT.\n",out);
      exit(1);
    }    

    if (!(fpker=fopen(ker,"r"))) {
      printf("SEGYConv: failed to open kernel file = %s for read. ABORT.\n",ker);
      exit(1);
    }

    // trace counter
    int nker=0;

    // trace workspace
    segy trinp;
    segy trout;
    segy trker;

    // file offsets
    off_t outoff = 0L;
    off_t inpoff = 0L;
    off_t keroff = 0L;

    // header workspace
    Value val;
    float dt;
    int nsinp;
    int nsout;
    int nsker;
    string dtstr="dt";
    string nsstr="ns";
    string rtstr="delrt";
    float tmp;

    // time offsets
    int inbeg;
    int outbeg;
    int kerbeg;
    
    // read first traces
    if (0 == fgettr(fpker,&trker)) {
      RVL::RVLException e;
      e<<"ERROR: SEGYConv\n";
      e<<"  fwd mode\n";
      e<<"  failed to read kernel trace 0\n";
      throw e;
    }
    if (0 == fgettr(fpout,&trout)) {
      RVL::RVLException e;
      e<<"ERROR: SEGYConv\n";
      e<<"  fwd mode\n";
      e<<"  failed to read output trace 0\n";
      throw e;
    }
    if (0 == fgettr(fpinp,&trinp)) {
      RVL::RVLException e;
      e<<"ERROR: SEGYConv\n";
      e<<"  fwd mode\n";
      e<<"  failed to read input trace 0\n";
      throw e;
    }
    
    // read trace headers, determine time, trace offsets
    gethdval(&trout,  (char*)(dtstr.c_str()), &val);
    dt=vtof(hdtype((char*)(dtstr.c_str())), val);
    dt*=0.001;
    // cerr<<"conv: dt="<<dt<<endl;
    
    gethdval(&trinp,  (char*)(nsstr.c_str()), &val);
    nsinp=vtoi(hdtype(nsstr.c_str()), val);
    gethdval(&trinp,  (char*)(rtstr.c_str()), &val);
    tmp=vtof(hdtype(rtstr.c_str()), val);
    if (tmp<0.0f) inbeg=int(-0.01 + tmp/dt);
    else inbeg=int(0.01 + tmp/dt);
    
    gethdval(&trout,  (char*)(nsstr.c_str()), &val);
    nsout=vtoi(hdtype(nsstr.c_str()), val);
    gethdval(&trout,  (char*)(rtstr.c_str()), &val);
    tmp=vtof(hdtype(rtstr.c_str()), val);
    if (tmp<0.0f) outbeg=int(-0.01 + tmp/dt);
    else outbeg=int(0.01 + tmp/dt);
    
    gethdval(&trker,  (char*)(nsstr.c_str()), &val);
    nsker=vtoi(hdtype(nsstr.c_str()), val);
    gethdval(&trker,  (char*)(rtstr.c_str()), &val);
    tmp=vtof(hdtype(rtstr.c_str()), val);
    if (tmp<0.0f) kerbeg=int(-0.01 + tmp/dt);
    else kerbeg=int(0.01 + tmp/dt);
    
    // compute offsets to min tracl - do this on first trace read
    inpoff = mnx*(nsinp*sizeof(float) + HDRBYTES);
    outoff = mnx*(nsout*sizeof(float) + HDRBYTES);
    keroff = mnx*(nsker*sizeof(float) + HDRBYTES);

    // seek to initial trace
    if (fseeko(fpinp,inpoff,SEEK_SET)) {
      RVL::RVLException e;
      e<<"Error: SEGYConv\n";
      e<<"  adj mode\n";	  
      e<<"reset failed on file "<<inp<<"\n";
      throw e;
    }
    if (fseeko(fpout,outoff,SEEK_SET)) {
      RVL::RVLException e;
      e<<"Error: SEGYConv\n";
      e<<"  adj mode\n";	  
      e<<"reset failed on file "<<out<<"\n";
      throw e;
    }
    if (fseeko(fpker,keroff,SEEK_SET)) {
      RVL::RVLException e;
      e<<"Error: SEGYConv\n";
      e<<"  adj mode\n";	  
      e<<"reset failed on file "<<out<<"\n";
      throw e;
    }

    nker = mnx;
      
    while (fgettr(fpker,&trker) && (nker < mxx)) {
      nker++;
	
      if (0 == fgettr(fpinp,&trinp)) {
	RVL::RVLException e;
	e<<"ERROR: SEGYConv\n";
	e<<"  fwd mode\n";
	e<<"  failed to read output trace "<<nker<<"\n";
	throw e;
      }

      // read output trace, then back up, to get header right
      off_t tmpoutoff = ftello(fpout);
      if (0 == fgettr(fpout,&trout)) {
	RVL::RVLException e;
	e<<"ERROR: SEGYConv\n";
	e<<"  fwd mode\n";
	e<<"  failed to read output trace "<<nker<<"\n";
	throw e;
      }
      if (fseeko(fpout,tmpoutoff,SEEK_SET)) {
	RVL::RVLException e;
	e<<"Error: SEGYConv\n";
	e<<"  adj mode\n";	  
	e<<"reset failed on file "<<out<<"\n";
	throw e;
      }
	
      // cerr<<"SEGYConv nout="<<nout<<" nker="<<nker<<" ninp="<<ninp<<endl;
       
      // cerr<<"SEGYConv inbeg="<<inbeg<<" outbeg="<<outbeg<<" kerbeg="<<kerbeg<<"\n";	
      // printf("SEGYConv ninp = %d nout = %d nker = %d outoff=%d\n", ninp, nout, nker, outoff);
      // printf("SEGYConv nsinp = %d nsout = %d nsker = %d\n",nsinp, nsout, nsker);
      
      if (adj==0) {
	int ishift=outbeg-inbeg-kerbeg;
	// cerr<<"SEGYConvolve: ishift="<<ishift<<endl;
	TSOpt::conv(ishift,nsout,nsinp,nsker,
		    trout.data,trinp.data,trker.data,dt);
      }
      else {
	int ishift=inbeg-kerbeg-outbeg;
	// cerr<<"SEGYConvolve: ishift="<<ishift<<endl;
	TSOpt::corr(ishift,nsout,nsinp,nsker,
		    trout.data,trinp.data,trker.data,dt);
      }
      fputtr(fpout,&trout);      
    }

    fclose(fpinp);
    fclose(fpout);
    fclose(fpker);
    ps_delete(&par);
    //iwave_fdestroy();
    
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
    exit(0);
  }
    
  catch (RVL::RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
  
}


