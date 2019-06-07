#include "cstd.h"
#include "utils.h"

//  Update 2D velocities when iv==0
//  OUTPUT: vx, vz
void esg_vstep2d(ireal ** restrict buoyancy,
                 ireal ** restrict sxx,
                 ireal ** restrict szz,
                 ireal ** restrict sxz,
                 IPNT gsc_vx, IPNT gec_vx,
                 IPNT gsc_vz, IPNT gec_vz,
                 IPNT lbc, IPNT rbc,
                 int radius, ireal ** restrict coeff,
                 ireal ** restrict vx,
                 ireal ** restrict vz
                 ){
    
    //  Update vx in the computational domain
    for(int i1 = gsc_vx[1]; i1 <= gec_vx[1]; i1++)
        for(int ir = 0; ir < radius; ir++)
            for(int i0 = gsc_vx[0]; i0 <= gec_vx[0]; i0++)
                vx[i1][i0] += buoyancy[i1][i0]*(coeff[1][ir]*(sxx[i1+ir][i0]-sxx[i1-ir-1][i0])
                                                +coeff[0][ir]*(sxz[i1][i0+ir]-sxz[i1][i0-ir-1]));
    
    //  Update vz in the computational domain
    for(int i1 = gsc_vz[1]; i1 <= gec_vz[1]; i1++)
        for(int ir = 0; ir < radius; ir++){
            for(int i0 = gsc_vz[0]; i0 <= gec_vz[0]; i0++)
                vz[i1][i0] += 0.25*(buoyancy[i1][i0]+buoyancy[i1][i0+1]+buoyancy[i1+1][i0]+buoyancy[i1+1][i0+1])
                *(coeff[1][ir]*(sxz[i1+ir+1][i0]-sxz[i1-ir][i0])
                  +coeff[0][ir]*(szz[i1][i0+ir+1]-szz[i1][i0-ir]));
        }
    
    //  Update ghost velocities
    //  Case 1: If domain contains z[start]
    if(lbc[0]){
        int i0 = gsc_vx[0]-1;
        for(int ir = 1; ir < radius; ir++)
            for(int i1 = gsc_vx[1]; i1 <= gec_vx[1]; i1++)
                vx[i1][i0-ir] = -vx[i1][i0+ir];
        
        i0 = gsc_vz[0];
        for(int ir = 1; ir < radius; ir++)
            for(int i1 = gsc_vz[1]; i1 <= gec_vz[1]; i1++)
                vz[i1][i0-ir] = -vz[i1][i0+ir-1];
    }
    
    //  Case 2: If domain contains z[end]
    if(rbc[0]){
        int i0 = gec_vx[0]+1;
        for(int ir = 1; ir < radius; ir++)
            for(int i1 = gsc_vx[1]; i1 <= gec_vx[1]; i1++)
                vx[i1][i0+ir] = -vx[i1][i0-ir];
        
        i0 = gec_vz[0];
        for(int ir = 1; ir < radius; ir++)
            for(int i1 = gsc_vz[1]; i1 <= gec_vz[1]; i1++)
                vz[i1][i0+ir] = -vz[i1][i0-ir+1];
    }
    
    //   Case 3: If domain contains x[start]
    if(lbc[1]){
        int i1 = gsc_vx[1]-1;
        for(int ir = 1; ir < radius; ir++)
            for(int i0 = gsc_vx[0]; i0 <= gec_vx[0]; i0++)
                vx[i1-ir][i0] = -vx[i1+ir][i0];
        
        i1 = gsc_vz[1];
        for(int ir = 1; ir < radius; ir++)
            for(int i0 = gsc_vz[0]; i0 <= gec_vz[0]; i0++)
                vz[i1-ir][i0] = -vz[i1+ir-1][i0];
    }
    
    //   Case 4: If domain contains x[end]
    if(rbc[1]){
        int i1 = gec_vx[1]+1;
        for(int ir = 1; ir < radius; ir++)
            for(int i0 = gsc_vx[0]; i0 <= gec_vx[0]; i0++)
                vx[i1+ir][i0] = -vx[i1-ir][i0];
        
        i1 = gec_vz[1];
        for(int ir = 1; ir < radius; ir++)
            for(int i0 = gsc_vz[0]; i0 <= gec_vz[0]; i0++)
                vz[i1+ir][i0] = -vz[i1-ir+1][i0];
    }
}


