#include "cstd.h"
#include "utils.h"
#include <float.h> //FLT_MIN
#include <stdlib.h> //malloc

//  Update 2D stresses when iv==1
//  OUTPUT: szz, sxx, sxz
void esg_sstep2d(ireal ** restrict c11,
                 ireal ** restrict c33,
                 ireal ** restrict c13,
                 ireal ** restrict c55,
                 ireal ** restrict vx,
                 ireal ** restrict vz,
                 IPNT gsc_szz, IPNT gec_szz,
                 IPNT gsc_sxz, IPNT gec_sxz,
                 IPNT lbc, IPNT rbc,
                 int radius, ireal ** restrict coeff,
                 ireal ** restrict sxx,
                 ireal ** restrict szz,
                 ireal ** restrict sxz,
                 int eflag,
                 int effective_media
                 ){
    //int effective_media = 1;
    
    //=============================================
    //  Approach1: using effective media parameters
    if(effective_media)
    {
        //  Update szz, sxx
        for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++)
            for(int ir = 0; ir < radius; ir++){
                for(int i0 = gsc_szz[0]; i0 <= gec_szz[0]; i0++){
                    //  Update szz
                    szz[i1][i0] += 2*c13[i1][i0]*c13[i1+1][i0]/(c13[i1][i0]+c13[i1+1][i0])
                    *coeff[1][ir]*(vx[i1+ir+1][i0]-vx[i1-ir][i0])
                    +2*c33[i1][i0]*c33[i1+1][i0]/(c33[i1][i0]+c33[i1+1][i0])
                    *coeff[0][ir]*(vz[i1][i0+ir]-vz[i1][i0-ir-1]);
                    //  Update sxx
                    sxx[i1][i0] += 2*c11[i1][i0]*c11[i1+1][i0]/(c11[i1][i0]+c11[i1+1][i0])
                    *coeff[1][ir]*(vx[i1+ir+1][i0]-vx[i1-ir][i0])
                    +2*c13[i1][i0]*c13[i1+1][i0]/(c13[i1][i0]+c13[i1+1][i0])
                    *coeff[0][ir]*(vz[i1][i0+ir]-vz[i1][i0-ir-1]);
                }
            }
        
        //  Update sxz
        ireal *c55_tmp = (ireal*)malloc((gec_sxz[0]-gsc_sxz[0]+1)*sizeof(ireal));
        for(int i1 = gsc_sxz[1]; i1 <= gec_sxz[1]; i1++)
            for(int ir = 0; ir < radius; ir++)
            {
                //precompute effective media parameters
                for(int i0 = gsc_sxz[0]; i0 <= gec_sxz[0]; i0++)
                {
                    if(c55[i1][i0]<FLT_MIN || c55[i1][i0+1]<FLT_MIN)
                        c55_tmp[i0-gsc_sxz[0]] = 0;
                    else
                        c55_tmp[i0-gsc_sxz[0]] = 2*c55[i1][i0]*c55[i1][i0+1]/(c55[i1][i0]+c55[i1][i0+1]);
                }
                for(int i0 = gsc_sxz[0]; i0 <= gec_sxz[0]; i0++)
                    sxz[i1][i0] += c55_tmp[i0-gsc_sxz[0]]
                    *(coeff[0][ir]*(vx[i1][i0+ir+1]-vx[i1][i0-ir])
                      +coeff[1][ir]*(vz[i1+ir][i0]-vz[i1-ir-1][i0]));
            }
        free(c55_tmp);
        
        //  Update ghost stresses
        //  Case 1: If domain contains z[start]
        if(lbc[0]){
            int i0 = gsc_szz[0]-1;
            
            if(eflag==1)
                for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++)
                    for(int ir = 0; ir < radius; ir++)
                        sxx[i1][i0] += 4*c13[i1][i0]*c13[i1+1][i0]/(c13[i1][i0]+c13[i1+1][i0])
                        *coeff[0][ir]*vz[i1][i0+ir];
            
            for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++){
                for(int ir = 0; ir < radius; ir++)
                    szz[i1][i0] += 4*c33[i1][i0]*c33[i1+1][i0]/(c33[i1][i0]+c33[i1+1][i0])
                    *coeff[0][ir]*vz[i1][i0+ir];
                for(int ir = 0; ir < radius-1; ir++)
                    szz[i1][i0-ir-1] = szz[i1][i0+ir+1];
            }
            
            i0 = gsc_sxz[0];
            for(int ir = 1; ir < radius; ir++)
                for(int i1 = gsc_sxz[1]; i1 <= gec_sxz[1]; i1++)
                    sxz[i1][i0-ir] = sxz[i1][i0+ir-1];
        }
        
        
        //  Case 2: If domain contains z[end]
        if(rbc[0]){
            int i0 = gec_szz[0]+1;
            
            if(eflag==1)
                for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++)
                    for(int ir = 0; ir < radius; ir++)
                        sxx[i1][i0] += -4*c13[i1][i0]*c13[i1+1][i0]/(c13[i1][i0]+c13[i1+1][i0])
                        *coeff[0][ir]*vz[i1][i0-ir-1];
            
            for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++){
                for(int ir = 0; ir < radius; ir++){
                    szz[i1][i0] += -4*c33[i1][i0]*c33[i1+1][i0]/(c33[i1][i0]+c33[i1+1][i0])
                    *coeff[0][ir]*vz[i1][i0-ir-1];
                    szz[i1][i0+ir+1] = szz[i1][i0-ir-1];
                }
                for(int ir = 0; ir < radius-1; ir++)
                    szz[i1][i0+ir+1] = szz[i1][i0-ir-1];
            }
            
            i0 = gec_sxz[0];
            for(int ir = 1; ir < radius; ir++)
                for(int i1 = gsc_sxz[1]; i1 <= gec_sxz[1]; i1++)
                    sxz[i1][i0+ir] = sxz[i1][i0-ir+1];
        }
        
        
        //   Case 3: If domain contains x[start]
        if(lbc[1]){
            int i1 = gsc_szz[1];
            for(int ir = 1; ir < radius; ir++)
                for(int i0 = gsc_szz[0]; i0 <= gec_szz[0]; i0++)
                    sxx[i1-ir][i0] = sxx[i1+ir-1][i0];
            
            
            i1 = gsc_sxz[1]-1;
            for(int i0 = gsc_sxz[0]; i0 <= gec_sxz[0]; i0++){
                
                //precompute effective media parameters
                ireal c55_tmp;
                if(c55[i1][i0] < FLT_MIN || c55[i1][i0+1] < FLT_MIN)
                    c55_tmp = 0;
                else
                    c55_tmp = 4*c55[i1][i0]*c55[i1][i0+1]/(c55[i1][i0]+c55[i1][i0+1]);
                
                for(int ir = 0; ir < radius; ir++)
                    sxz[i1][i0] += c55_tmp*coeff[1][ir]*vz[i1+ir][i0];
                for(int ir = 0; ir < radius-1; ir++)
                    sxz[i1-ir-1][i0] = sxz[i1+ir+1][i0];
            }
        }
        
        
        //   Case 4: If domain contains x[end]
        if(rbc[1]){
            int i1 = gec_szz[1];
            for(int ir = 1; ir < radius; ir++)
                for(int i0 = gsc_szz[0]; i0 <= gec_szz[0]; i0++)
                    sxx[i1+ir][i0] = sxx[i1-ir+1][i0];
            
            
            i1 = gec_sxz[1]+1;
            for(int i0 = gsc_sxz[0]; i0 <= gec_sxz[0]; i0++){
                
                //precompute effective media parameters
                ireal c55_tmp;
                if(c55[i1][i0] < FLT_MIN || c55[i1][i0+1] < FLT_MIN)
                    c55_tmp = 0;
                else
                    c55_tmp = 4*c55[i1][i0]*c55[i1][i0+1]/(c55[i1][i0]+c55[i1][i0+1]);
                
                for(int ir = 0; ir < radius; ir++)
                    sxz[i1][i0] += -c55_tmp*coeff[1][ir]*vz[i1-ir-1][i0];
                for(int ir = 0; ir < radius-1; ir++)
                    sxz[i1+ir+1][i0] = sxz[i1-ir-1][i0];
            }
        }
        
    }else{
        //=========================================================
        //  Approach2: using intuitive interpolation for parameters
        //  Update szz, sxx
        for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++)
            for(int ir = 0; ir < radius; ir++){
                for(int i0 = gsc_szz[0]; i0 <= gec_szz[0]; i0++){
                    //  Update szz
                    szz[i1][i0] += 0.5*(c13[i1][i0]+c13[i1+1][i0])
                    *coeff[1][ir]*(vx[i1+ir+1][i0]-vx[i1-ir][i0])
                    +0.5*(c33[i1][i0]+c33[i1+1][i0])
                    *coeff[0][ir]*(vz[i1][i0+ir]-vz[i1][i0-ir-1]);
                    //  Update sxx
                    sxx[i1][i0] += 0.5*(c11[i1][i0]+c11[i1+1][i0])
                    *coeff[1][ir]*(vx[i1+ir+1][i0]-vx[i1-ir][i0])
                    +0.5*(c13[i1][i0]+c13[i1+1][i0])
                    *coeff[0][ir]*(vz[i1][i0+ir]-vz[i1][i0-ir-1]);
                }
            }
        
        //  Update sxz
        for(int i1 = gsc_sxz[1]; i1 <= gec_sxz[1]; i1++)
            for(int ir = 0; ir < radius; ir++)
                for(int i0 = gsc_sxz[0]; i0 <= gec_sxz[0]; i0++)
                    sxz[i1][i0] += 0.5*(c55[i1][i0]+c55[i1][i0+1])
                    *(coeff[0][ir]*(vx[i1][i0+ir+1]-vx[i1][i0-ir])
                      +coeff[1][ir]*(vz[i1+ir][i0]-vz[i1-ir-1][i0]));
        
        //  Update ghost stresses
        //  Case 1: If domain contains z[start]
        if(lbc[0]){
            int i0 = gsc_szz[0]-1;
            
            if(eflag==1)
                for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++)
                    for(int ir = 0; ir < radius; ir++)
                        sxx[i1][i0] += (c13[i1][i0]+c13[i1+1][i0])
                        *coeff[0][ir]*vz[i1][i0+ir];
            
            for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++){
                for(int ir = 0; ir < radius; ir++)
                    szz[i1][i0] += (c33[i1][i0]+c33[i1+1][i0])
                    *coeff[0][ir]*vz[i1][i0+ir];
                for(int ir = 0; ir < radius-1; ir++)
                    szz[i1][i0-ir-1] = szz[i1][i0+ir+1];
            }
            
            i0 = gsc_sxz[0];
            for(int ir = 1; ir < radius; ir++)
                for(int i1 = gsc_sxz[1]; i1 <= gec_sxz[1]; i1++)
                    sxz[i1][i0-ir] = sxz[i1][i0+ir-1];
        }
        
        
        //  Case 2: If domain contains z[end]
        if(rbc[0]){
            int i0 = gec_szz[0]+1;
            
            if(eflag==1)
                for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++)
                    for(int ir = 0; ir < radius; ir++)
                        sxx[i1][i0] += -(c13[i1][i0]+c13[i1+1][i0])
                        *coeff[0][ir]*vz[i1][i0-ir-1];
            
            for(int i1 = gsc_szz[1]; i1 <= gec_szz[1]; i1++){
                for(int ir = 0; ir < radius; ir++){
                    szz[i1][i0] += -(c33[i1][i0]+c33[i1+1][i0])
                    *coeff[0][ir]*vz[i1][i0-ir-1];
                    szz[i1][i0+ir+1] = szz[i1][i0-ir-1];
                }
                for(int ir = 0; ir < radius-1; ir++)
                    szz[i1][i0+ir+1] = szz[i1][i0-ir-1];
            }
            
            i0 = gec_sxz[0];
            for(int ir = 1; ir < radius; ir++)
                for(int i1 = gsc_sxz[1]; i1 <= gec_sxz[1]; i1++)
                    sxz[i1][i0+ir] = sxz[i1][i0-ir+1];
        }
        
        
        //   Case 3: If domain contains x[start]
        if(lbc[1]){
            int i1 = gsc_szz[1];
            for(int ir = 1; ir < radius; ir++)
                for(int i0 = gsc_szz[0]; i0 <= gec_szz[0]; i0++)
                    sxx[i1-ir][i0] = sxx[i1+ir-1][i0];
            
            
            i1 = gsc_sxz[1]-1;
            for(int i0 = gsc_sxz[0]; i0 <= gec_sxz[0]; i0++){
                for(int ir = 0; ir < radius; ir++)
                    sxz[i1][i0] += (c55[i1][i0]+c55[i1][i0+1])*coeff[1][ir]*vz[i1+ir][i0];
                for(int ir = 0; ir < radius-1; ir++)
                    sxz[i1-ir-1][i0] = sxz[i1+ir+1][i0];
            }
        }
        
        
        //   Case 4: If domain contains x[end]
        if(rbc[1]){
            int i1 = gec_szz[1];
            for(int ir = 1; ir < radius; ir++)
                for(int i0 = gsc_szz[0]; i0 <= gec_szz[0]; i0++)
                    sxx[i1+ir][i0] = sxx[i1-ir+1][i0];
            
            
            i1 = gec_sxz[1]+1;
            for(int i0 = gsc_sxz[0]; i0 <= gec_sxz[0]; i0++){
                for(int ir = 0; ir < radius; ir++)
                    sxz[i1][i0] += -(c55[i1][i0]+c55[i1][i0+1])*coeff[1][ir]*vz[i1-ir-1][i0];
                for(int ir = 0; ir < radius-1; ir++)
                    sxz[i1+ir+1][i0] = sxz[i1-ir-1][i0];
            }
        }
    }
}

