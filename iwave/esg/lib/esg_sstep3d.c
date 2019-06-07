#include "cstd.h"
#include "utils.h"

//  Update 3D stresses when iv==1
//  OUTPUT: sxx, syy, szz, syz, sxz, sxy
void esg_sstep3d(ireal *** restrict c11,
                 ireal *** restrict c22,
                 ireal *** restrict c33,
                 ireal *** restrict c23,
                 ireal *** restrict c13,
                 ireal *** restrict c12,
                 ireal *** restrict c44,
                 ireal *** restrict c55,
                 ireal *** restrict c66,
                 ireal *** restrict vx,
                 ireal *** restrict vy,
                 ireal *** restrict vz,
                 IPNT gsc_szz, IPNT gec_szz,
                 IPNT gsc_syz, IPNT gec_syz,
                 IPNT gsc_sxz, IPNT gec_sxz,
                 IPNT gsc_sxy, IPNT gec_sxy,
                 IPNT lbc, IPNT rbc,
                 int radius, ireal ** restrict coeff,
                 ireal *** restrict sxx,
                 ireal *** restrict syy,
                 ireal *** restrict szz,
                 ireal *** restrict syz,
                 ireal *** restrict sxz,
                 ireal *** restrict sxy,
                 int eflag
                 ){
    
    //  Update sxx, szz
    for(int i2 = gsc_szz[2]; i2 <= gec_szz[2]; i2++)
        for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++)
            for(int ir = 0; ir < radius; ir++){
                for(int i0 = gsc_szz[0]; i0 <= gec_szz[0]; i0++){
                    //  Update szz
                    szz[i2][i1][i0] += 0.5*(c13[i2][i1][i0]+c13[i2][i1+1][i0])
                    *coeff[1][ir]*(vx[i2][i1+ir+1][i0]-vx[i2][i1-ir][i0])
                    +0.5*(c33[i2][i1][i0]+c33[i2][i1+1][i0])
                    *coeff[0][ir]*(vz[i2][i1][i0+ir]-vz[i2][i1][i0-ir-1])
                    +0.5*(c23[i2][i1][i0]+c23[i2][i1+1][i0])
                    *coeff[2][ir]*(vy[i2+ir][i1][i0]-vy[i2-ir-1][i1][i0]);
                    //  Update sxx
                    sxx[i2][i1][i0] += 0.5*(c11[i2][i1][i0]+c11[i2][i1+1][i0])
                    *coeff[1][ir]*(vx[i2][i1+ir+1][i0]-vx[i2][i1-ir][i0])
                    +0.5*(c13[i2][i1][i0]+c13[i2][i1+1][i0])
                    *coeff[0][ir]*(vz[i2][i1][i0+ir]-vz[i2][i1][i0-ir-1])
                    +0.5*(c12[i2][i1][i0]+c12[i2][i1+1][i0])
                    *coeff[2][ir]*(vy[i2+ir][i1][i0]-vy[i2-ir-1][i1][i0]);
                    //  Update syy
                    syy[i2][i1][i0] += 0.5*(c12[i2][i1][i0]+c12[i2][i1+1][i0])
                    *coeff[1][ir]*(vx[i2][i1+ir+1][i0]-vx[i2][i1-ir][i0])
                    +0.5*(c23[i2][i1][i0]+c23[i2][i1+1][i0])
                    *coeff[0][ir]*(vz[i2][i1][i0+ir]-vz[i2][i1][i0-ir-1])
                    +0.5*(c22[i2][i1][i0]+c22[i2][i1+1][i0])
                    *coeff[2][ir]*(vy[i2+ir][i1][i0]-vy[i2-ir-1][i1][i0]);
                }
            }
    
    //  Update syz
    for(int i2 = gsc_syz[2]; i2 <= gec_syz[2]; i2++)
        for(int i1 = gsc_syz[1]; i1 <= gec_syz[1]; i1++)
            for(int ir = 0; ir < radius; ir++)
                for(int i0 = gsc_syz[0]; i0 <= gec_syz[0]; i0++)
                    syz[i2][i1][i0] += 0.125*(c44[i2][i1][i0]+c44[i2][i1][i0+1]
                                              +c44[i2][i1+1][i0]+c44[i2][i1+1][i0+1]
                                              +c44[i2+1][i1][i0]+c44[i2+1][i1][i0+1]
                                              +c44[i2+1][i1+1][i0]+c44[i2+1][i1+1][i0+1])
                    *(coeff[0][ir]*(vy[i2][i1][i0+ir+1]-vy[i2][i1][i0-ir])
                      +coeff[2][ir]*(vz[i2+ir+1][i1][i0]-vz[i2-ir][i1][i0]));
    
    //  Update sxz
    for(int i2 = gsc_sxz[2]; i2 <= gec_sxz[2]; i2++)
        for(int i1 = gsc_sxz[1]; i1 <= gec_sxz[1]; i1++)
            for(int ir = 0; ir < radius; ir++)
                for(int i0 = gsc_sxz[0]; i0 <= gec_sxz[0]; i0++)
                    sxz[i2][i1][i0] += 0.5*(c55[i2][i1][i0]+c55[i2][i1][i0+1])
                    *(coeff[0][ir]*(vx[i2][i1][i0+ir+1]-vx[i2][i1][i0-ir])
                      +coeff[1][ir]*(vz[i2][i1+ir][i0]-vz[i2][i1-ir-1][i0]));
    
    //  Update sxy
    for(int i2 = gsc_sxy[2]; i2 <= gec_sxy[2]; i2++)
        for(int i1 = gsc_sxy[1]; i1 <= gec_sxy[1]; i1++)
            for(int ir = 0; ir < radius; ir++)
                for(int i0 = gsc_sxy[0]; i0 <= gec_sxy[0]; i0++)
                    sxy[i2][i1][i0] += 0.5*(c66[i2][i1][i0]+c66[i2+1][i1][i0])
                    *(coeff[2][ir]*(vx[i2+ir+1][i1][i0]-vx[i2-ir][i1][i0])
                      +coeff[1][ir]*(vy[i2][i1+ir][i0]-vy[i2][i1-ir-1][i0]));
    
    
    //  Update ghost stresses
    //  Case 1: If domain contains z[start]
    if(lbc[0]){
        int i0 = gsc_szz[0]-1;
        
        if(eflag==1)
            for(int i2 = gsc_szz[2]; i2 <= gec_szz[2]; i2++)
                for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++)
                    for(int ir = 0; ir < radius; ir++){
                        sxx[i2][i1][i0] += (c13[i2][i1][i0]+c13[i2][i1+1][i0])
                        *coeff[0][ir]*vz[i2][i1][i0+ir];
                        syy[i2][i1][i0] += (c23[i2][i1][i0]+c23[i2][i1+1][i0])
                        *coeff[0][ir]*vz[i2][i1][i0+ir];
                    }
        
        
        for(int i2 = gsc_szz[2]; i2 <= gec_szz[2]; i2++)
            for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++){
                for(int ir = 0; ir < radius; ir++)
                    szz[i2][i1][i0] += (c33[i2][i1][i0]+c33[i2][i1+1][i0])
                    *coeff[0][ir]*vz[i2][i1][i0+ir];
                for(int ir = 1; ir < radius; ir++)
                    szz[i2][i1][i0-ir] = szz[i2][i1][i0+ir];
            }
        
        
        i0 = gsc_syz[0];
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_syz[2]; i2 <= gec_syz[2]; i2++)
                for(int i1 = gsc_syz[1]; i1 <= gec_syz[1]; i1++)
                    syz[i2][i1][i0-ir] = syz[i2][i1][i0+ir-1];
        
        i0 = gsc_sxz[0];
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_sxz[2]; i2 <= gec_sxz[2]; i2++)
                for(int i1 = gsc_sxz[1]; i1 <= gec_sxz[1]; i1++)
                    sxz[i2][i1][i0-ir] = sxz[i2][i1][i0+ir-1];
    }
    
    //  If domain contains z[end]
    if(rbc[0]){
        int i0 = gec_szz[0]+1;
        
        if(eflag==1)
            for(int i2 = gsc_szz[2]; i2 <= gec_szz[2]; i2++)
                for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++)
                    for(int ir = 0; ir < radius; ir++){
                        sxx[i2][i1][i0] += -(c13[i2][i1][i0]+c13[i2][i1+1][i0])
                        *coeff[0][ir]*vz[i2][i1][i0-ir-1];
                        syy[i2][i1][i0] += -(c23[i2][i1][i0]+c23[i2][i1+1][i0])
                        *coeff[0][ir]*vz[i2][i1][i0-ir-1];
                    }
        
        for(int i2 = gsc_szz[2]; i2 <= gec_szz[2]; i2++)
            for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++){
                for(int ir = 0; ir < radius; ir++)
                    szz[i2][i1][i0] += -(c33[i2][i1][i0]+c33[i2][i1+1][i0])
                    *coeff[0][ir]*vz[i2][i1][i0-ir-1];
                for(int ir = 1; ir < radius; ir++)
                    szz[i2][i1][i0+ir] = szz[i2][i1][i0-ir];
            }
        
        i0 = gec_syz[0];
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_syz[2]; i2 <= gec_syz[2]; i2++)
                for(int i1 = gsc_syz[1]; i1 <= gec_syz[1]; i1++)
                    syz[i2][i1][i0+ir] = syz[i2][i1][i0-ir+1];
        
        i0 = gec_sxz[0];
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_sxz[2]; i2 <= gec_sxz[2]; i2++)
                for(int i1 = gsc_sxz[1]; i1 <= gec_sxz[1]; i1++)
                    sxz[i2][i1][i0+ir] = sxz[i2][i1][i0-ir+1];
    }
    
    //  If domain contains x[start]
    if(lbc[1]){
        int i1 = gsc_szz[1];
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_szz[2]; i2 <= gec_szz[2]; i2++)
                for(int i0 = gsc_szz[0]; i0 <= gec_szz[0]; i0++)
                    sxx[i2][i1-ir][i0] = sxx[i2][i1+ir-1][i0];
        
        i1 = gsc_sxz[1]-1;
        for(int i2 = gsc_sxz[2]; i2 <= gec_sxz[2]; i2++)
            for(int i0 = gsc_sxz[0]; i0 <= gec_sxz[0]; i0++){
                for(int ir = 0; ir < radius; ir++)
                    sxz[i2][i1][i0] += (c55[i2][i1][i0]+c55[i2][i1][i0+1])
                    *coeff[1][ir]*vz[i2][i1+ir][i0];
                for(int ir = 1; ir < radius; ir++)
                    sxz[i2][i1-ir][i0] = sxz[i2][i1+ir][i0];
            }
        
        i1 = gsc_sxy[1]-1;
        for(int i2 = gsc_sxy[2]; i2 <= gec_sxy[2]; i2++)
            for(int i0 = gsc_sxy[0]; i0 <= gec_sxy[0]; i0++){
                for(int ir = 0; ir < radius; ir++)
                    sxy[i2][i1][i0] += (c66[i2][i1][i0]+c66[i2+1][i1][i0])
                    *coeff[1][ir]*vy[i2][i1+ir][i0];
                for(int ir = 1; ir < radius; ir++)
                    sxy[i2][i1-ir][i0] = sxy[i2][i1+ir][i0];
            }
    }
    
    //  If domain contains x[end]
    if(rbc[1]){
        int i1 = gec_szz[1];
        for(int ir = 1; ir < radius; ir++)
            for(int i2 = gsc_szz[2]; i2 <= gec_szz[2]; i2++)
                for(int i0 = gsc_szz[0]; i0 <= gec_szz[0]; i0++)
                    sxx[i2][i1+ir][i0] = sxx[i2][i1-ir+1][i0];
        
        i1 = gec_sxz[1]+1;
        for(int i2 = gsc_sxz[2]; i2 <= gec_sxz[2]; i2++)
            for(int i0 = gsc_sxz[0]; i0 <= gec_sxz[0]; i0++){
                for(int ir = 0; ir < radius; ir++)
                    sxz[i2][i1][i0] += -(c55[i2][i1][i0]+c55[i2][i1][i0+1])
                    *coeff[1][ir]*vz[i2][i1-ir-1][i0];
                for(int ir = 1; ir < radius; ir++)
                    sxz[i2][i1+ir][i0] = sxz[i2][i1-ir][i0];
            }
        
        i1 = gec_sxy[1]+1;
        for(int i2 = gsc_sxy[2]; i2 <= gec_sxy[2]; i2++)
            for(int i0 = gsc_sxy[0]; i0 <= gec_sxy[0]; i0++){
                for(int ir = 0; ir < radius; ir++)
                    sxy[i2][i1][i0] += -(c66[i2][i1][i0]+c66[i2+1][i1][i0])
                    *coeff[1][ir]*vy[i2][i1-ir-1][i0];
                for(int ir = 1; ir < radius; ir++)
                    sxy[i2][i1+ir][i0] = sxy[i2][i1-ir][i0];
            }
    }
    
    
    //  If domain contains y[start]
    if(lbc[2]){
        int i2 = gsc_szz[2]-1;
        
        if(eflag==1)
            for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++)
                for(int i0 = gsc_szz[0]; i0 <= gec_szz[0]; i0++)
                    for(int ir = 0; ir < radius; ir++){
                        szz[i2][i1][i0] += (c23[i2][i1][i0]+c23[i2][i1+1][i0])
                        *coeff[2][ir]*vy[i2+ir][i1][i0];
                        sxx[i2][i1][i0] += (c12[i2][i1][i0]+c12[i2][i1+1][i0])
                        *coeff[2][ir]*vy[i2+ir][i1][i0];
                    }
        
        
        for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++)
            for(int i0 = gsc_szz[0]; i0 <= gec_szz[0]; i0++){
                for(int ir = 0; ir < radius; ir++)
                    syy[i2][i1][i0] += (c22[i2][i1][i0]+c22[i2][i1+1][i0])
                    *coeff[2][ir]*vy[i2+ir][i1][i0];
                for(int ir = 1; ir < radius; ir++)
                    syy[i2-ir][i1][i0] = syy[i2+ir][i1][i0];
            }
        
        i2 = gsc_syz[2];
        for(int ir = 1; ir < radius; ir++)
            for(int i1 = gsc_syz[1]; i1 <= gec_syz[1]; i1++)
                for(int i0 = gsc_syz[0]; i0 <= gec_syz[0]; i0++)
                    syz[i2-ir][i1][i0] = syz[i2+ir-1][i1][i0];
        
        i2 = gsc_sxy[2];
        for(int ir = 1; ir < radius; ir++)
            for(int i1 = gsc_sxy[1]; i1 <= gec_sxy[1]; i1++)
                for(int i0 = gsc_sxy[0]; i0 <= gec_sxy[0]; i0++)
                    sxy[i2-ir][i1][i0] = sxy[i2+ir-1][i1][i0];
    }
    
    
    //  If domain contains y[end]
    if(rbc[2]){
        int i2 = gec_szz[2]+1;
        
        if(eflag==1)
            for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++)
                for(int i0 = gsc_szz[0]; i0 <= gec_szz[0]; i0++)
                    for(int ir = 0; ir < radius; ir++){
                        szz[i2][i1][i0] += -(c23[i2][i1][i0]+c23[i2][i1+1][i0])
                        *coeff[2][ir]*vy[i2-ir-1][i1][i0];
                        sxx[i2][i1][i0] += -(c12[i2][i1][i0]+c12[i2][i1+1][i0])
                        *coeff[2][ir]*vy[i2-ir-1][i1][i0];
                    }
        
        for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++)
            for(int i0 = gsc_szz[0]; i0 <= gec_szz[0]; i0++){
                for(int ir = 0; ir < radius; ir++)
                    syy[i2][i1][i0] += -(c22[i2][i1][i0]+c22[i2][i1+1][i0])
                    *coeff[2][ir]*vy[i2-ir-1][i1][i0];
                for(int ir = 1; ir < radius; ir++)
                    syy[i2+ir][i1][i0] = syy[i2-ir][i1][i0];
            }
        
        i2 = gec_syz[2];
        for(int ir = 1; ir < radius; ir++)
            for(int i1 = gsc_syz[1]; i1 <= gec_syz[1]; i1++)
                for(int i0 = gsc_syz[0]; i0 <= gec_syz[0]; i0++)
                    syz[i2+ir][i1][i0] = syz[i2-ir+1][i1][i0];
        
        
        i2 = gec_sxy[2];
        for(int ir = 1; ir < radius; ir++)
            for(int i1 = gsc_sxy[1]; i1 <= gec_sxy[1]; i1++)
                for(int i0 = gsc_sxy[0]; i0 <= gec_sxy[0]; i0++)
                    sxy[i2+ir][i1][i0] = sxy[i2-ir+1][i1][i0];
    }
}

