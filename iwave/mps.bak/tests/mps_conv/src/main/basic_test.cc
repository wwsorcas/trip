#include "asg_defn.hh"
#include "asg.hh"
#include "MPS_conv.hh"

using RVL::valparse;

IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"bulkmod",    0, true,  1 },
  {"buoyancy",   1, true,  1 },
  {"source_p",   2, true,  2 },
  {"data_p",     2, false, 2 },
  {"",           0, false, 0 }
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
    cerr << "\n////////////////////////////////////////////////////////////\n"
	 << "Basic test of convolution kernels:\n"
	 << "////////////////////////////////////////////////////////////\n";
#ifdef IWAVE_USE_MPI
    }
#endif

    TSOpt::MPS_KEYS mks;
    mks.appx_ord  = "appx_ord";
    mks.grid_file = "bulkmod";
    mks.MPS_file  = "MPS_file";
    mks.RHS_files.push_back("source_p");
    mks.hdr_files.push_back("data_p");

    ps_printall(*pars,stderr);

    //testing convolution kernel
    cerr << "\n[] Testing convolution kernel\n\n";

    string x_file = valparse<string>(*pars,"x_file");
    vector<string> y_files(1);
    y_files[0] = valparse<string>(*pars,"y_files");
    vector<string> k_files(1);
    k_files[0] = valparse<string>(*pars,"k_files");
    vector<TSOpt::Rtuple> s_pos(1);
    s_pos[0].coor[0] = 0;
    s_pos[0].coor[1] = 0;
    s_pos[0].coor[2] = 0;

    int MPS_dim = 1;
    int MPS_idx = 0;


    MPS_conv_kern( s_pos,
		   MPS_dim,
		   MPS_idx,
		   x_file,
		   y_files,
		   k_files );

    //testing convolution kernel
    cerr << "\n[] Testing correlation kernel\n\n";

    string z_file = valparse<string>(*pars,"z_file");

    MPS_corr_kern( s_pos,
		   MPS_dim,
		   MPS_idx,
		   z_file,
		   y_files,
		   k_files );

    ps_delete(&pars);
    iwave_fdestroy();

#ifdef IWAVE_USE_MPI    
    MPI_Finalize();
#endif
    
  }
  catch(RVLException &e){
    e << "Exiting with error!\n";
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
}
