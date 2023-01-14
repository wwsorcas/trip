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
  "of traces) with an input gather to produce an output gather. For convolution, the input",
  "gather is constrained to contain either a single trace or (at least) the same number of",
  "traces as the kernel. The output gather should have the same geometry as the kernel.",
  "For adjoint convolution, the roles of, and constraints on, input and output are",
  "interchanged.",
  " ",
  "Intended use: kernel is Green's function of a causal point source simulation. For ",
  "convolution, input is either point source wavelet (single trace), or a wavelet for",
  "each receiver. For adjoint convolution, single trace output is stacked.",
  " ",
  "required inputs:",
  "in       = name of su file containing inpur traces:",
  "out      = name of su file containing output traces",
  "ker      = name of file containing signature trace", 
  "mode     = fwd -   convolve",
  "           adj -   adjoint-convolve (cross-correlate)",
  NULL};

int main(int argc, char ** argv) {

  /* ARGUMENTS */
  char * inp;      /* input file name */
  char * out;      /* output file name */
  char * ker;      /* kernel file name */
  int adj;         /* 0 = fwd, 1 = adj */

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

    adj=0;
    ps_flint(*par,"adj",&adj);

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

    // single trace flag
    int single=1;

    int ninp=0;
    int nout=0;
    int nker=0;

    segy trinp;
    segy trout;
    segy trker;

    off_t outoff = 0L;

    Value val;
    float dt;
    int nsinp;
    int nsout;
    int nsker;
    
    int inbeg;
    int outbeg;
    int kerbeg;
    
    string dtstr="dt";
    string nsstr="ns";
    string rtstr="delrt";
    float tmpinp;
    float tmpout;
    float tmpker;
    float tmpacc;
    float tmp;
    
    // workspace for accumulation of adjoint in single trace case
    float * buf = NULL;

    // case single: either (adj=0 and ninp=1) or (adj=1 and nout=1)
    // determine if file has 1 or > 1 traces
    // for input file, fwd case, if 1 trace then don't meess
    // with it - will be neither read nor written further,
    // data is in segy

    // in single trace case, at the end of this block trinp (adj=0)
    // or trout (adj=1) contains correct headers for first (and only)
    // input/output trace
    
    if (adj==0) {
      while (fgettr(fpinp,&trinp) && ninp<3) ninp++;
      if (ninp==0) {
	RVL::RVLException e;
	e<<"ERROR: SEGYConv\n";
	e<<"  fwd mode\n";
	e<<"  failed to read first input trace\n";
	throw e;
      }
      if (ninp==1) single=1;
      else {
	single=0;
	// reset count
	ninp=0;
      }
      // reset input file since it will need to be re-read
      if (fseeko(fpinp,0L,SEEK_SET)) {
	RVL::RVLException e;
	e<<"Error: SEGYConv\n";
	e<<"  fwd mode\n";	  
	e<<"reset failed on file "<<inp<<"\n";
	throw e;
      }
      
    }
    else {
      while (fgettr(fpout,&trout) && nout<3) nout++;
      if (nout==0) {
	RVL::RVLException e;
	e<<"ERROR: SEGYConv\n";
	e<<"  adj mode\n";
	e<<"  failed to read first output trace\n";
	throw e;
      }
      if (nout==1) single=1;
      else {
	single=0;
	// reset count
	nout=0;
      }

      // reset output file in any case since it will need to be written
      if (fseeko(fpout,0L,SEEK_SET)) {
	RVL::RVLException e;
	e<<"Error: SEGYConv\n";
	e<<"  adj mode\n";	  
	e<<"reset failed on file "<<out<<"\n";
	throw e;
      }
      
    }

    if (single==0) {
      
      while (fgettr(fpker,&trker)) {
	nker++;
	
	// reading block
	outoff = ftello(fpout);
	if (fgettr(fpout,&trout)) nout++;
	else {
	  RVL::RVLException e;
	  e<<"ERROR: SEGYConv\n";
	  e<<"  fwd mode\n";
	  e<<"  failed to read output trace "<<nker<<"\n";
	  throw e;
	}
	if (fseeko(fpout,outoff,SEEK_SET)) {
	  RVL::RVLException e;
	  e<<"Error: SEGYConv\n";
	  e<<"  adj mode\n";	  
	  e<<"reset failed on file "<<out<<"\n";
	  throw e;
	}

	// cerr<<"SEGYConv outoff="<<outoff<<endl;
	
	if (fgettr(fpinp,&trinp)) ninp++;
	else {
	  RVL::RVLException e;
	  e<<"ERROR: SEGYConv\n";
	  e<<"  fwd mode\n";
	  e<<"  failed to read output trace "<<nker<<"\n";
	  throw e;
	}
	
	// cerr<<"SEGYConv nout="<<nout<<" nker="<<nker<<" ninp="<<ninp<<endl;
	
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

    }
    // single trace case

    else {

      while (fgettr(fpker,&trker)) {
	nker++;
	
	// reading block
	
	if (adj==0) {
	  // read output trace, reset for write - input trace already read  
	  outoff = ftello(fpout);
	  if (fgettr(fpout,&trout)) nout++;
	  else {
	    RVL::RVLException e;
	    e<<"ERROR: SEGYConv\n";
	    e<<"  fwd mode\n";
	    e<<"  failed to read output trace "<<nker<<"\n";
	    throw e;
	  }
	  if (fseeko(fpout,outoff,SEEK_SET)) {
	    RVL::RVLException e;
	    e<<"Error: SEGYConv\n";
	    e<<"  adj mode\n";	  
	    e<<"reset failed on file "<<out<<"\n";
	    throw e;
	  }
	  // cerr<<"SEGYConv outoff="<<outoff<<endl;
	}
	else {
          // read input trace - output trace already read
	  if (fgettr(fpinp,&trinp)) ninp++;
	  else {
	    RVL::RVLException e;
	    e<<"ERROR: SEGYConv\n";
	    e<<"  fwd mode\n";
	    e<<"  failed to read input trace "<<nker<<"\n";
	    throw e;
	  }

	}

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
	
	// cerr<<"SEGYConv inbeg="<<inbeg<<" outbeg="<<outbeg<<" kerbeg="<<kerbeg<<"\n";	
	// cerr<<"SEGYConv nout="<<nout<<" nker="<<nker<<" ninp="<<ninp<<endl;
	// printf("SEGYConv nsinp = %d nsout = %d nsker = %d\n",nsinp, nsout, nsker);      
	
	if (adj==0) {
	  int ishift=outbeg-inbeg-kerbeg;
	  // cerr<<"SEGYConvolve: ishift="<<ishift<<endl;
	  TSOpt::conv(ishift,nsout,nsinp,nsker,
		      trout.data,trinp.data,trker.data,dt);
	  // cerr<<"SEGYConv adj="<<adj<<" write output trace\n";	  
	  fputtr(fpout,&trout);	  
	}
	else {
	  // accumulation buffer
	  if (!buf) {
	    // cerr<<"SEGYConv adj allocate accum buffer\n";
	    buf = (float *)usermalloc_(nsout*sizeof(float));
	  }
	  int ishift=inbeg-kerbeg-outbeg;
	  // cerr<<"SEGYConv ishift="<<ishift<<endl;
	  // cerr<<"SEGYConv save current trace accumulation to buffer\n";
	  //	  memcpy(buf,trout.data,nsout*sizeof(float));
	  for (int i=0; i<nsout; i++) {
	    if (nker==1) trout.data[i]=0.0f;
	    buf[i]=trout.data[i];
	  }
	  // overwrite output trace with current cross-corr
	  TSOpt::corr(ishift,nsout,nsinp,nsker,
		      buf,trinp.data,trker.data,dt);
	  // add accumulated trace back in
	  tmpacc=0.0f;
	  tmpinp=0.0f;
	  tmpker=0.0f;
	  for (int i=0; i<nsout; i++) {
	    trout.data[i] += buf[i];
	    tmpacc+=buf[i]*buf[i];
	    tmpout+=trout.data[i]*trout.data[i];
	  }
	  for (int i=0; i<nsinp; i++) {
	    tmpinp+=trinp.data[i]*trinp.data[i];
	  }
	  for (int i=0; i<nsker; i++) {
	    tmpker+=trker.data[i]*trker.data[i];
	  }
	  // cerr<<"SEGYConv single adj sq acc = "<<tmpacc<<endl;
	  // cerr<<"SEGYConv single adj sq inp = "<<tmpinp<<endl;
	  // cerr<<"SEGYConv single adj sq ker = "<<tmpker<<endl;
	  // cerr<<"SEGYConv single adj sq out = "<<tmpout<<endl;
	}
	  
      }
      
      if (adj) {
	// cerr<<"SEGYConv adj="<<adj<<" write output trace\n";
	fputtr(fpout,&trout);
      }
    }
      
    fclose(fpinp);
    fclose(fpout);
    fclose(fpker);
    if (buf) {
      userfree_(buf);
      buf=NULL;
    }
    ps_delete(&par);
    iwave_fdestroy();
    
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


