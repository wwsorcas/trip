#include "ssappx.hh"

int xargc;
char **xargv;

int main(int argc, char ** argv) {

  try {

#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
#endif

#ifdef IWAVE_USE_MPI    
    if( retrieveGlobalRank()==0 ){
#endif
      cerr << "\n------------------------------------------------------------\n" 
	   << "Running test kernel for ssappx\n"
	   << "------------------------------------------------------------\n";
      //ps_printall(*pars,stderr);
#ifdef IWAVE_USE_MPI
    }
#endif


    int s, q, N;
    vector< vector< vector<int> > > MIDX;

    //case 1
    s=0;
    q=2;
    N = q+s;

    cerr << "CASE 1: s="<<s<<", q="<<2<<"\n";    
    MIDX.resize(N);
    for( int ell=0; ell<N; ell++){
      cerr<< "ell="<<ell<<"\n";

      MIDX[ell] = TSOpt::midx_set(q,s,ell);

      cerr << "  MIDX for ell="<<ell<<"\n";
      for( int i=0; i<MIDX[ell].size(); i++ ){
	cerr << "    (";
	for( int j=0; j<MIDX[ell][i].size(); j++ ){
	  cerr << MIDX[ell][i][j] <<",";
	}
	cerr << ")\n";
      }
    }
    MIDX.clear();
    cerr << "\n";

    //case 2
    s=0;
    q=4;
    N = q+s;

    cerr << "CASE 2: s="<<s<<", q="<<q<<"\n";

    MIDX.resize(N);
    for( int ell=0; ell<N; ell++){
      MIDX[ell] = TSOpt::midx_set(q,s,ell);

      cerr << "  MIDX for ell="<<ell<<"\n";
      for( int i=0; i<MIDX[ell].size(); i++ ){
	cerr << "    (";
	for( int j=0; j<MIDX[ell][i].size(); j++ ){
	  cerr << MIDX[ell][i][j] <<",";
	}
	cerr << ")\n";
      }
    }
    MIDX.clear();
    cerr << "\n";

    //case 3
    s=1;
    q=2;
    N = q+s;

    cerr << "CASE 3: s="<<s<<", q="<<2<<"\n";    

    MIDX.resize(N);
    for( int ell=0; ell<N; ell++){
      MIDX[ell] = TSOpt::midx_set(q,s,ell);

      cerr << "  MIDX for ell="<<ell<<"\n";
      for( int i=0; i<MIDX[ell].size(); i++ ){
	cerr << "    (";
	for( int j=0; j<MIDX[ell][i].size(); j++ ){
	  cerr << MIDX[ell][i][j] <<",";
	}
	cerr << ")\n";
      }
    }
    MIDX.clear();
    cerr << "\n";

    //case 4
    s=1;
    q=4;
    N = q+s;

    cerr << "CASE 4: s="<<s<<", q="<<q<<"\n";

    MIDX.resize(N);
    for( int ell=0; ell<N; ell++){
      MIDX[ell] = TSOpt::midx_set(q,s,ell);

      cerr << "  MIDX for ell="<<ell<<"\n";
      for( int i=0; i<MIDX[ell].size(); i++ ){
	cerr << "    (";
	for( int j=0; j<MIDX[ell][i].size(); j++ ){
	  cerr << MIDX[ell][i][j] <<",";
	}
	cerr << ")\n";
      }
    }
    MIDX.clear();
    cerr << "\n";

    //case 5
    s=2;
    q=2;
    N = q+s;

    cerr << "CASE 5: s="<<s<<", q="<<q<<"\n";

    MIDX.resize(N);
    for( int ell=0; ell<N; ell++){
      MIDX[ell] = TSOpt::midx_set(q,s,ell);

      cerr << "  MIDX for ell="<<ell<<"\n";
      for( int i=0; i<MIDX[ell].size(); i++ ){
	cerr << "    (";
	for( int j=0; j<MIDX[ell][i].size(); j++ ){
	  cerr << MIDX[ell][i][j] <<",";
	}
	cerr << ")\n";
      }
    }
    MIDX.clear();
    cerr << "\n";

    //case 6
    s=2;
    q=4;
    N = q+s;

    cerr << "CASE 6: s="<<s<<", q="<<q<<"\n";

    MIDX.resize(N);
    for( int ell=0; ell<N; ell++){
      MIDX[ell] = TSOpt::midx_set(q,s,ell);

      cerr << "  MIDX for ell="<<ell<<"\n";
      for( int i=0; i<MIDX[ell].size(); i++ ){
	cerr << "    (";
	for( int j=0; j<MIDX[ell][i].size(); j++ ){
	  cerr << MIDX[ell][i][j] <<",";
	}
	cerr << ")\n";
      }
    }
    MIDX.clear();
    cerr << "\n";


#ifdef IWAVE_USE_MPI   
    if( retrieveGlobalRank()==0 ){
#endif        
      cerr << "------------------------------------------------------------\n" 
	   << "Exiting test kernel for assappx\n"
	   << "------------------------------------------------------------\n\n";
#ifdef IWAVE_USE_MPI    
    }
    MPI_Finalize();
#endif

  }
  catch(RVLException &e){
    e << "Exiting with error!\n";
#ifdef IWAVE_USE_MPI
    if( retrieveGlobalRank()==0 ){
#endif
      e.write(cerr);
#ifdef IWAVE_USE_MPI
    }   
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
}
