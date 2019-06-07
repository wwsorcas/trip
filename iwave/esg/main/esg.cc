#include "iwsim.hh"
#include <par.h>
#include "esg_defn.hh"
#include "esg_selfdoc.h"


enum{
    D_GRID  = 0,
    D_BUOY  = 1,
    D_C11   = 2,
    D_C22   = 3,
    D_C33   = 4,
    D_C23   = 5,
    D_C13   = 6,
    D_C12   = 7,
    D_C44   = 8,
    D_C55   = 9,
    D_C66   = 10,
    D_VX    = 11,
    D_VY    = 12,
    D_VZ    = 13,
    D_SXX   = 14,
    D_SYY   = 15,
    D_SZZ   = 16,
    D_SYZ   = 17,
    D_SXZ   = 18,
    D_SXY   = 19,
    D_TRACE = 20
};

IOKEY IWaveInfo::iwave_iokeys[]
={
    {"grid_info",       D_GRID, true, false},
    {"source_vx",       D_VX,  true,   false},
    {"source_vy",       D_VY,  true,   false},
    {"source_vz",       D_VZ,  true,   false},
    {"source_sxx",       D_SXX,  true,   false},
    {"source_syy",       D_SYY,  true,   false},
    {"source_szz",       D_SZZ,  true,   false},
    {"buoyancy",        D_BUOY,  true,   true},
    {"c11",       D_C11,  true,   true},
    {"c22",       D_C22,  true,   true},
    {"c33",       D_C33,  true,   true},
    {"c23",       D_C23,  true,   true},
    {"c13",       D_C13,  true,   true},
    {"c12",       D_C12,  true,   true},
    {"c44",       D_C44,  true,   true},
    {"c55",       D_C55,  true,   true},
    {"c66",       D_C66,  true,   true},
    {"data_vx",         D_VX,  false,  true},
    {"data_vy",         D_VY,  false,  true},
    {"data_vz",         D_VZ,  false,  true},
    {"data_sxx",        D_SXX, false,  true},
    {"data_syy",        D_SYY, false,  true},
    {"data_szz",        D_SZZ, false,  true},
    {"data_syz",        D_SYZ, false,  true},
    {"data_sxz",        D_SXZ, false,  true},
    {"data_sxy",        D_SXY, false,  true},
    {"movie_vx",        D_VX,  false,  false},
    {"movie_vy",        D_VY,  false,  false},
    {"movie_vz",        D_VZ,  false,  false},
    {"data_trace",        D_TRACE,  false,  true},
    {"",                0,  false,  false}
    
};

int main(int argc, char ** argv) {
    
    try {
        
#ifdef IWAVE_USE_MPI
        MPI_Init(&argc,&argv);
#endif
        int rk=0;
#ifdef IWAVE_USE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD,&rk);
#endif
        
        if (rk==0 && argc<2) {
            pagedoc();
            exit(0);
        }
        
        TSOpt::IWaveApply(argc,argv);
        
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

