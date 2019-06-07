#include "cstd.h"
#include "utils.h"

//  update 3D velocity when iv==0
//  OUTPUT: vx,vy,vz
void esg_vstep3d(ireal*** restrict buoyancy,
                 ireal*** restrict sxx,
                 ireal*** restrict syy,
                 ireal*** restrict szz,
                 ireal*** restrict syz,
                 ireal*** restrict sxz,
                 ireal*** restrict sxy,
                 IPNT gsc_vx, IPNT gec_vx,
                 IPNT gsc_vy, IPNT gec_vy,
                 IPNT gsc_vz, IPNT gec_vz,
                 IPNT lbc, IPNT rbc,
                 int radius, ireal** restrict coeff,
                 ireal*** restrict vx,
                 ireal*** restrict vy,
                 ireal*** restrict vz
                 ){
    
    //  Update vz
    for(int i2 = gsc_vz[2]; i2 <= gec_vz[2]; i2++)
        for(int i1 = gsc_vz[1]; i1 <= gec_vz[1]; i1++)
            for(int ir = 0; ir < radius; ir++)
                for(int i0 = gsc_vz[0]; i0 <= gec_vz[0]; i0++)
                    vz[i2][i1][i0] += 0.25*(buoyancy[i2][i1][i0]
                                            +buoyancy[i2][i1][i0+1]
                                            +buoyancy[i2][i1+1][i0]
                                            +buoyancy[i2][i1+1][i0+1])
                    *(coeff[1][ir]*(sxz[i2][i1+ir+1][i0]-sxz[i2][i1-ir][i0])
                      +coeff[0][ir]*(szz[i2][i1][i0+ir+1]-szz[i2][i1][i0-ir])
                      +coeff[2][ir]*(syz[i2+ir][i1][i0]-syz[i2-ir-1][i1][i0]));
    
    //  Update vx
    for(int i2 = gsc_vx[2]; i2 <= gec_vx[2]; i2++)
        for(int i1 = gsc_vx[1]; i1 <= gec_vx[1]; i1++)
            for(int ir = 0; ir < radius; ir++)
                for(int i0 = gsc_vx[0]; i0 <= gec_vx[0]; i0++)
                    vx[i2][i1][i0] += buoyancy[i2][i1][i0]*(coeff[1][ir]*(sxx[i2][i1+ir][i0]-sxx[i2][i1-ir-1][i0])
                                                            +coeff[2][ir]*(sxy[i2+ir][i1][i0]-sxy[i2-ir-1][i1][i0])
                                                            +coeff[0][ir]*(sxz[i2][i1][i0+ir]-sxz[i2][i1][i0-ir-1]));
    
    //  Update vy
    for(int i2 = gsc_vy[2]; i2 <= gec_vy[2]; i2++)
        for(int i1 = gsc_vy[1]; i1 <= gec_vy[1]; i1++)
            for(int ir = 0; ir < radius; ir++)
                for(int i0 = gsc_vy[0]; i0 <= gec_vy[0]; i0++)
                    vy[i2][i1][i0] += 0.25*(buoyancy[i2][i1][i0]
                                            +buoyancy[i2][i1+1][i0]
                                            +buoyancy[i2+1][i1][i0]
                                            +buoyancy[i2+1][i1+1][i0])
                    *(coeff[1][ir]*(sxy[i2][i1+ir+1][i0]-sxy[i2][i1-ir][i0])
                      +coeff[0][ir]*(syz[i2][i1][i0+ir]-syz[i2][i1][i0-ir-1])
                      +coeff[2][ir]*(syy[i2+ir+1][i1][i0]-syy[i2-ir][i1][i0]));
    
    
    //  Update ghost velocities
    //  Case 1: If domain contains z[start]
    if(lbc[0]){
        int i0 = gsc_vz[0];
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_vz[2]; i2 <= gec_vz[2]; i2++)
                for(int i1 = gsc_vz[1]; i1 <= gec_vz[1]; i1++)
                    vz[i2][i1][i0-ir] = -vz[i2][i1][i0+ir-1];
        
        i0 = gsc_vx[0]-1;
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_vx[2]; i2 <= gec_vx[2]; i2++)
                for(int i1 = gsc_vx[1]; i1 <= gec_vx[1]; i1++)
                    vx[i2][i1][i0-ir] = -vx[i2][i1][i0+ir];
        
        i0 = gsc_vy[0]-1;
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_vy[2]; i2 <= gec_vy[2]; i2++)
                for(int i1 = gsc_vy[1]; i1 <= gec_vy[1]; i1++)
                    vy[i2][i1][i0-ir] = -vy[i2][i1][i0+ir];
        
    }
    
    //  Case 2: If domain contains z[end]
    if(rbc[0]){
        int i0 = gec_vz[0];
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_vz[2]; i2 <= gec_vz[2]; i2++)
                for(int i1 = gsc_vz[1]; i1 <= gec_vz[1]; i1++)
                    vz[i2][i1][i0+ir] = -vz[i2][i1][i0-ir+1];
        
        i0 = gec_vx[0]+1;
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_vx[2]; i2 <= gec_vx[2]; i2++)
                for(int i1 = gsc_vx[1]; i1 <= gec_vx[1]; i1++)
                    vx[i2][i1][i0+ir] = -vx[i2][i1][i0-ir];
        
        i0 = gec_vy[0]+1;
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_vy[2]; i2 <= gec_vy[2]; i2++)
                for(int i1 = gsc_vy[1]; i1 <= gec_vy[1]; i1++)
                    vy[i2][i1][i0+ir] = -vy[i2][i1][i0-ir];
        
    }
    
    //   Case 3: If domain contains x[start]
    if(lbc[1]){
        int i1 = gsc_vz[1];
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_vz[2]; i2 <= gec_vz[2]; i2++)
                for(int i0 = gsc_vz[0]; i0 <= gec_vz[0]; i0++)
                    vz[i2][i1-ir][i0] = - vz[i2][i1+ir-1][i0];
        
        i1 = gsc_vx[1]-1;
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_vx[2]; i2 <= gec_vx[2]; i2++)
                for(int i0 = gsc_vx[0]; i0 <= gec_vx[0]; i0++)
                    vx[i2][i1-ir][i0] = - vx[i2][i1+ir][i0];
        
        i1 = gsc_vy[1];
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_vy[2]; i2 <= gec_vy[2]; i2++)
                for(int i0 = gsc_vy[0]; i0 <= gec_vy[0]; i0++)
                    vy[i2][i1-ir][i0] = - vy[i2][i1+ir-1][i0];
    }
    
    //   Case 4: If domain contains x[end]
    if(rbc[1]){
        int i1 = gec_vz[1];
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_vz[2]; i2 <= gec_vz[2]; i2++)
                for(int i0 = gsc_vz[0]; i0 <= gec_vz[0]; i0++)
                    vz[i2][i1+ir][i0] = - vz[i2][i1-ir+1][i0];
        
        i1 = gec_vx[1]+1;
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_vx[2]; i2 <= gec_vx[2]; i2++)
                for(int i0 = gsc_vx[0]; i0 <= gec_vx[0]; i0++)
                    vx[i2][i1+ir][i0] = - vx[i2][i1-ir][i0];
        
        i1 = gec_vy[1];
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_vy[2]; i2 <= gec_vy[2]; i2++)
                for(int i0 = gsc_vy[0]; i0 <= gec_vy[0]; i0++)
                    vy[i2][i1+ir][i0] = - vy[i2][i1-ir+1][i0];
    }
    
    //   Case 5: If domain contains y[start]
    if(lbc[2]){
        int i2 = gsc_vz[2]-1;
        for(int ir = 1; ir < radius; ir++)
            for(int i1 = gsc_vz[1]; i1 <= gec_vz[1]; i1++)
                for(int i0 = gsc_vz[0]; i0 <= gec_vz[0]; i0++)
                    vz[i2-ir][i1][i0] = -vz[i2+ir][i1][i0];
        
        i2 = gsc_vx[2]-1;
        for(int ir = 1; ir < radius; ir++)
            for(int i1 = gsc_vx[1]; i1 <= gec_vx[1]; i1++)
                for(int i0 = gsc_vx[0]; i0 <= gec_vx[0]; i0++)
                    vx[i2-ir][i1][i0] = -vx[i2+ir][i1][i0];
        
        i2 = gsc_vy[2];
        for(int ir = 1; ir < radius; ir++)
            for(int i1 = gsc_vy[1]; i1 <= gec_vy[1]; i1++)
                for(int i0 = gsc_vy[0]; i0 <= gec_vy[0]; i0++)
                    vy[i2-ir][i1][i0] = -vy[i2+ir-1][i1][i0];
    }
    
    //   Case 6: If domain contains y[end]
    if(rbc[2]){
        int i2 = gec_vz[2]+1;
        for(int ir = 1; ir < radius; ir++)
            for(int i1 = gsc_vz[1]; i1 <= gec_vz[1]; i1++)
                for(int i0 = gsc_vz[0]; i0 <= gec_vz[0]; i0++)
                    vz[i2+ir][i1][i0] = -vz[i2-ir][i1][i0];
        
        i2 = gec_vx[2]+1;
        for(int ir = 1; ir < radius; ir++)
            for(int i1 = gsc_vx[1]; i1 <= gec_vx[1]; i1++)
                for(int i0 = gsc_vx[0]; i0 <= gec_vx[0]; i0++)
                    vx[i2+ir][i1][i0] = -vx[i2-ir][i1][i0];
        
        i2 = gec_vy[2];
        for(int ir = 1; ir < radius; ir++)
            for(int i1 = gsc_vy[1]; i1 <= gec_vy[1]; i1++)
                for(int i0 = gsc_vy[0]; i0 <= gec_vy[0]; i0++)
                    vy[i2+ir][i1][i0] = -vy[i2-ir+1][i1][i0];
    }
}

