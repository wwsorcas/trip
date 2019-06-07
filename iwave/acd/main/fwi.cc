#include "iwopt.hh"
#include <par.h>

#include "acd_defn.hh"
#include "acdfwi_selfdoc.h"

IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csq",    0, 1, 1 },
  {"data",   1, 0, 1 },
  {"source", 1, 1, 0 },
  {"",       0, 0, 0 }
};


int main(int argc, char ** argv) {

  try {

#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
#endif

    if (retrieveGlobalRank()==0 && argc<2) {
      pagedoc();
      exit(0);
    }

    TSOpt::GridMAOpBuilder pb;

    TSOpt::IWaveOpApply(argc,argv,TSOpt::StraightOLS,pb);

#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
  }
  catch (RVL::RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
  
}
