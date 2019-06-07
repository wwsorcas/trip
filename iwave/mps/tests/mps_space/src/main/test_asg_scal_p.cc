#include "asg_defn.hh"
#include "asg.hh"
#include "test_kernel.hh"

IOKEY IWaveInfo::iwave_iokeys[]
= {
  {BULKMOD,    0, true,  1 },
  {BUOYANCY,   1, true,  1 },
  {SOURCE_P,   2, true,  2 },
  {DATA_P,     2, false, 2 },
  {"",         0, false, 0 }
};


int xargc;
char **xargv;

int main(int argc, char ** argv) {

  try {

#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
#endif

    PARARRAY * pars = NULL;
    FILE * stream = NULL;
    TSOpt::IWaveEnvironment(argc, argv, 0, &pars, &stream);

#ifdef IWAVE_USE_MPI
    if( retrieveGlobalRank()==0 ){
#endif
      cerr << "\n------------------------------------------------------------\n" 
	   << "Running test_asg_scal_p.x\n"
	   << "------------------------------------------------------------\n";
      //ps_printall(*pars,stderr);
#ifdef IWAVE_USE_MPI
    }
#endif

    //consistency check of driver, data, and src type
    string driver_type = valparse<string>(*pars,"driver");
    string data_type   = valparse<string>(*pars,"data_type");
    string src_type    = valparse<string>(*pars,"src_type");

    if( driver_type.compare("asg")!=0    ||
	src_type.compare("pressure")!=0  || 
	data_type.compare("pressure")!=0 ){
      RVLException e;
      e << "Driver/src/data type is unexpected!\n"
	<< "  driver_type="<< driver_type <<"\n"
	<< "  src_type="   << src_type    <<"\n"
	<< "  data_type="  << data_type   <<"\n";
      throw e;
    }

    //running test kernel
    test_kernel<Scal_MPS_Space>(*pars,*stream);

    //clean up
    ps_delete(&pars);
    iwave_fdestroy();

#ifdef IWAVE_USE_MPI    
    if( retrieveGlobalRank()==0 ){
#endif            
      cerr << "------------------------------------------------------------\n" 
	   << "Exiting test_asg_scal_p.x\n"
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
