#include "esg.hh"
#include "utils.h"
#include "rarray.h"

using RVL::valparse;
using RVL::RVLException;

enum{
    D_SIZE  = 21,
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

//  max FD half-order
#define MAXK 12


extern ireal * sgcoeffs(int k);
extern "C" void esg_sstep2d(ireal** restrict c11,
                            ireal** restrict c33,
                            ireal** restrict c13,
                            ireal** restrict c55,
                            ireal** restrict vx,
                            ireal** restrict vz,
                            IPNT gsc_szz, IPNT gec_szz,
                            IPNT gsc_sxz, IPNT gec_sxz,
                            IPNT lbc, IPNT rbc,
                            int radius,
                            ireal** restrict coeff,
                            ireal** restrict sxx,
                            ireal** restrict szz,
                            ireal** restrict sxz,
                            int eflag,
                            int effective_media
                            );

extern "C" void esg_sstep3d(ireal*** restrict c11,
                            ireal*** restrict c22,
                            ireal*** restrict c33,
                            ireal*** restrict c23,
                            ireal*** restrict c13,
                            ireal*** restrict c12,
                            ireal*** restrict c44,
                            ireal*** restrict c55,
                            ireal*** restrict c66,
                            ireal*** restrict vx,
                            ireal*** restrict vy,
                            ireal*** restrict vz,
                            IPNT gsc_szz, IPNT gec_szz,
                            IPNT gsc_syz, IPNT gec_syz,
                            IPNT gsc_sxz, IPNT gec_sxz,
                            IPNT gsc_sxy, IPNT gec_sxy,
                            IPNT lbc, IPNT rbc,
                            int radius, ireal** restrict coeff,
                            ireal*** restrict sxx,
                            ireal*** restrict syy,
                            ireal*** restrict szz,
                            ireal*** restrict syz,
                            ireal*** restrict sxz,
                            ireal*** restrict sxy,
                            int eflag
                            );

extern "C" void esg_vstep2d(ireal** restrict buoyancy,
                            ireal** restrict sxx,
                            ireal** restrict szz,
                            ireal** restrict sxz,
                            IPNT gsc_vx, IPNT gec_vx,
                            IPNT gsc_vz, IPNT gec_vz,
                            IPNT lbc, IPNT rbc,
                            int radius,
                            ireal** restrict coeff,
                            ireal** restrict vx,
                            ireal** restrict vz
                            );

extern "C" void esg_vstep3d(ireal*** restrict buoyancy,
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
                            );


//----------FUNCTION: initialize model
//----------INPUT:    pars
//----------OUTPUT:   stream, model
int esg_modelinit(PARARRAY pars,
                  FILE *stream,
                  IMODEL & model){
    
    //    RVLException e_init;
    //    e_init<<"\n\033[1;35m----> ERROR: init()\033[0m\n";
    //    throw e_init;
    
    try{
        int err=0;
        
        //  Allocate esgpars
        ESG_TS_PARS *esgpars=(ESG_TS_PARS*)usermalloc_(sizeof(ESG_TS_PARS));
        if(esgpars==NULL){
            err=E_ALLOC;
            fprintf(stream,"\nERROR:  esg_modelinit: failed to allocate esgpars\n");
            return err;
        }
        
        //  ==============================
        //  Set esgpars
        //  ==============================
        
        //  Set esgpars->cp,cs,k
        esgpars->cp = valparse<ireal>(pars,"cp",0);
        esgpars->cs = valparse<ireal>(pars,"cs",0);
        esgpars->k = valparse<int>(pars,"order",1);
#ifdef IWAVE_VERBOSE
        fprintf(stream, "\nNOTE: initializing ESG with half-order = %d\n", esgpars->k);
#endif
        if(esgpars->k<1 || esgpars->k>MAXK){
            RVLException e;
            e<<"\nERROR:   esgmodelinit: stencil radius is out of range\n";
            e<<"           Range [1,"<<MAXK<<"], stencil radius"<<esgpars->k<<"\n";
            throw e;
        }
        
        //  Get MPI cartesian coordinates
        IPNT cdims;
        IPNT crank;
#ifdef IWAVE_USE_MPI
        IPNT cpers;
#endif
        IASN(cdims, IPNT_1);
        IASN(crank, IPNT_0);
#ifdef IWAVE_USE_MPI
        MPI_Comm cm = retrieveComm();
        if(MPI_Cart_get(cm, RARR_MAX_NDIM, cdims, cpers, crank)!=MPI_SUCCESS){
            RVLException e;
            e<<"\nERROR:   esg_modelinit: cannot get MPI cartesian coordinates\n";
            throw e;
        }
        MPI_Bcast((void*)(&(model.g.dim)),1,MPI_INT,0,cm);
#endif
        //  Set esgpars->ndim,dt
        esgpars->ndim = model.g.dim;
        esgpars->dt = model.tsind.dt;
        
        
        //  Set boundary flags: esgpars->lbc,rbc
        IASN(esgpars->lbc, IPNT_0);
        IASN(esgpars->rbc, IPNT_0);
        for(int i=0; i<esgpars->ndim; i++){
            if(crank[i]==0) esgpars->lbc[i]=1;
            if(crank[i]==cdims[i]-1) esgpars->rbc[i]=1;
        }
        
        //  Set esgpars->coeffs
        esgpars->coeffs = (float**)usermalloc_(RARR_MAX_NDIM*sizeof(float*));
        
        for(int idim=0; idim<RARR_MAX_NDIM; idim++)
            (esgpars->coeffs)[idim]=NULL;
        
        RPNT dxs;
        get_d(dxs, model.g);
        for(int idim=0; idim<esgpars->ndim; idim++){
            
            if(dxs[idim]<=0.0){
                RVLException e;
                e<<"\nERROR:   esg_modelinit: negative grid step\n";
                e<<"           dim="<<idim<<", step="<<dxs[idim]<<"\n";
                throw e;
            }
            
            ireal lam = esgpars->dt/dxs[idim];
            (esgpars->coeffs)[idim] = sgcoeffs(esgpars->k);
            for (int is=0; is<esgpars->k; is++)
                (esgpars->coeffs)[idim][is] *= lam;
            
        }
        
        //  Set esgpars->eflag
        esgpars->eflag = valparse<int>(pars,"eflag",0);
        
        //  Set esgpars->effective_media
        esgpars->effective_media = valparse<int>(pars,"effective_media",1);
        
        //  Set esgpars->CellSize
        esgpars->CellSize = 1;
        for(int idim=0; idim < esgpars->ndim; idim++)
            esgpars->CellSize *= model.g.axes[idim].d;
        
        //  Set esgpars->h; assume h is the same on every axis
        esgpars->h = model.g.axes[0].d;
        
        //  ==============================
        //  Set model parameters
        //  ==============================
        model.active.resize(D_SIZE);
        
        model.active[D_GRID]="gridinfo";
        model.active[D_BUOY]="buoyancy";
        model.active[D_C11]="c11";
        model.active[D_C33]="c33";
        model.active[D_C13]="c13";
        model.active[D_C55]="c55";
        model.active[D_VX]="vx";
        model.active[D_VZ]="vz";
        model.active[D_SXX]="sxx";
        model.active[D_SZZ]="szz";
        model.active[D_SXZ]="sxz";
        model.active[D_TRACE]="trace";
        if(esgpars->ndim == 3){
            model.active[D_C22]="c22";
            model.active[D_C23]="c23";
            model.active[D_C12]="c12";
            model.active[D_C44]="c44";
            model.active[D_C66]="c66";
            model.active[D_VY]="vy";
            model.active[D_SYY]="syy";
            model.active[D_SYZ]="syz";
            model.active[D_SXY]="sxy";
        }
        
        //  Currently esg has no PML
        for(int idim=0; idim<esgpars->ndim; idim++){
            model.nls[idim]=0;
            model.nrs[idim]=0;
        }
        
        model.specs=(void*)esgpars;
        
        return 0;
    }
    catch (RVLException &e){
        e<<"\nERROR:    esg_modelinit\n";
        throw e;
    }
}







//----------FUNCTION: free model->specs
//----------INPUT:    fdpars
void esg_modeldest(void** fdpars){
    
    ESG_TS_PARS* esgpars = (ESG_TS_PARS*)(*fdpars);
    
    for (int idim=0; idim<RARR_MAX_NDIM; idim++)
        if ((esgpars->coeffs)[idim])
            userfree_((esgpars->coeffs)[idim]);
    
    userfree_(*fdpars);
}









//----------FUNCTION: compute time step: rhs, dt
//----------INPUT:    pars, g
//----------OUTPUT:   stream, rhs, dt
int esg_timegrid(PARARRAY* pars,
                 FILE* stream,
                 grid const &g,
                 ireal &dt,
                 ireal &rhs){
    
    try{
        ireal cp = valparse<ireal>(*pars,"cp",0);
        ireal cfl  = valparse<ireal>(*pars,"cfl",0.5);
        
        //  compute min grid step
        ireal a = g.axes[0].d;
        for(int idim=1; idim<g.dim; idim++)
            a = iwave_min(a, g.axes[idim].d);
        
        //  check cfl & a
        if ( (a < REAL_EPS) || (cfl < REAL_EPS) ) {
            fprintf(stream,"\nERROR:  esg_timegrid: either min dx=%e or cfl=%e "
                    " too small\n", a, cfl);
            return E_BADINPUT;
        }
        
        //  compute dt
        dt = a*cfl/(cp*sqrt((float)(g.dim)));
        
        //  read dt from pars
        ireal dttmp = dt;
        dt = valparse<ireal>(*pars,"dt",dttmp);
        if(dt<=0){
            fprintf(stream,"\nERROR:  esg_timegrid: dt=%f\n", dt);
            return E_BADINPUT;
        }
        
        rhs = dt;
        
        return 0;
    }
    catch (RVLException & e) {
        e<<"\nERROR:    esg_timegrid\n";
        throw e;
    }
    
}









//----------FUNCTION: simulate dom at internal time step iv
//----------INPUT:    dom, iv, fdpars, fwd
//----------OUTPUT:   dom
void esg_timestep(std::vector<RDOM *> dom,
                  bool fwd,
                  int iv,
                  void* fdpars){
    
    //  Check dom.size(), fwd, ndim, iv
    //  Currently assume dom.size()==1, fwd==true, ndim=2, iv=0,1
    if(dom.size() > 1){
        RVLException e;
        e<<"\nERROR:    esg_timestep(): dom.size()="<<dom.size()<<">1.\n";
        throw e;
    }
    if (fwd != true){
        RVLException e;
        e<<"\nERROR:    esg_timestep(): fwd="<<fwd<<">\n";
        throw e;
    }
    //  Get ndim
    int ndim;
    ra_ndim(&(dom[0]->_s[D_BUOY]),&ndim);
    
    if(iv<0 || iv>1){
        RVLException e;
        e<<"\nERROR:    esg_timestep(): iv="<<iv<<">\n";
        throw e;
    }
    
    
    ESG_TS_PARS *esgpars = (ESG_TS_PARS*) fdpars;
    
    //  Set start&end indices of simulation arrays
    IPNT gsc_vx, gec_vx;
    IPNT gsc_vy, gec_vy;
    IPNT gsc_vz, gec_vz;
    IPNT gsc_szz, gec_szz;
    IPNT gsc_syz, gec_syz;
    IPNT gsc_sxz, gec_sxz;
    IPNT gsc_sxy, gec_sxy;
    IPNT gsc_trace, gec_trace;
    
    IASN(gsc_vx,IPNT_0);
    IASN(gec_vx,IPNT_0);
    rd_gse(dom[0], D_VX, gsc_vx, gec_vx);
    IASN(gsc_vz,IPNT_0);
    IASN(gec_vz,IPNT_0);
    rd_gse(dom[0], D_VZ, gsc_vz, gec_vz);
    IASN(gsc_szz,IPNT_0);
    IASN(gec_szz,IPNT_0);
    rd_gse(dom[0], D_SXX, gsc_szz, gec_szz);
    IASN(gsc_sxz,IPNT_0);
    IASN(gec_sxz,IPNT_0);
    rd_gse(dom[0], D_SXZ, gsc_sxz, gec_sxz);
    IASN(gsc_trace,IPNT_0);
    IASN(gec_trace,IPNT_0);
    rd_gse(dom[0], D_TRACE, gsc_trace, gec_trace);
    
    if(ndim == 3){
        IASN(gsc_vy,IPNT_0);
        IASN(gec_vy,IPNT_0);
        rd_gse(dom[0], D_VY, gsc_vy, gec_vy);
        
        IASN(gsc_syz,IPNT_0);
        IASN(gec_syz,IPNT_0);
        rd_gse(dom[0], D_SYZ, gsc_syz, gec_syz);
        
        IASN(gsc_sxy,IPNT_0);
        IASN(gec_sxy,IPNT_0);
        rd_gse(dom[0], D_SXY, gsc_sxy, gec_sxy);
    }
    
    if(ndim==2){
        //  Update velocity (iv==0) or stresses (iv==1)
        register ireal ** restrict buoyancy = (dom[0]->_s)[D_BUOY]._s2;
        register ireal ** restrict c11 = (dom[0]->_s)[D_C11]._s2;
        register ireal ** restrict c33 = (dom[0]->_s)[D_C33]._s2;
        register ireal ** restrict c13 = (dom[0]->_s)[D_C13]._s2;
        register ireal ** restrict c55 = (dom[0]->_s)[D_C55]._s2;
        register ireal ** restrict vx = (dom[0]->_s)[D_VX]._s2;
        register ireal ** restrict vz = (dom[0]->_s)[D_VZ]._s2;
        register ireal ** restrict sxx = (dom[0]->_s)[D_SXX]._s2;
        register ireal ** restrict szz = (dom[0]->_s)[D_SZZ]._s2;
        register ireal ** restrict sxz = (dom[0]->_s)[D_SXZ]._s2;
        register ireal ** restrict trace = (dom[0]->_s)[D_TRACE]._s2;
        
        //  Update velocity vx, vz
        try{
            //  Currently velocities are at time level [n-1/2]
            //            stresses are at time step [n]
            
            
            if(iv==0){
                
                //  create vx_tmp, vz_tmp to record vx_tmp, vz_tmp at time level [n-1/2] when eflag==1
                RARR vx_tmp;
                RARR vz_tmp;
                ra_setnull(&vx_tmp);
                ra_setnull(&vz_tmp);
                
                if(esgpars->eflag==1){
                    if (ra_create(&vx_tmp, gsc_vx, gec_vx)) {
                        RVLException e;
                        e<<"\nERROR: esg_timestep: create tmp velocity arrays.\n";
                        throw e;
                    }
                    if (ra_create(&vz_tmp, gsc_vz, gec_vz)) {
                        RVLException e;
                        e<<"\nERROR: esg_timestep: create tmp velocity arrays.\n";
                        throw e;
                    }
                    ra_zero(&vx_tmp);
                    ra_zero(&vz_tmp);
                    if (ra_copy(&vx_tmp, &((dom[0]->_s)[D_VX]))){
                        RVLException e;
                        e<<"\nERROR: esg_timestep: copy velocity to tmp velocity arrays.\n";
                        throw e;
                    }
                    if (ra_copy(&vz_tmp, &((dom[0]->_s)[D_VZ]))){
                        RVLException e;
                        e<<"\nERROR: esg_timestep: copy velocity to tmp velocity arrays.\n";
                        throw e;
                    }
                }
                
                
                esg_vstep2d(
                            buoyancy, sxx, szz, sxz,
                            gsc_vx, gec_vx, gsc_vz, gec_vz,
                            esgpars->lbc, esgpars->rbc, esgpars->k, esgpars->coeffs,
                            vx, vz);
                
                if(esgpars->eflag==1){
                    
                    ireal energy=0;
                    ireal CellSize = esgpars->CellSize;
                    
                    //  add vx
                    for(int i1=gsc_vx[1]; i1<=gec_vx[1]; i1++)
                        for(int i0=gsc_vx[0]; i0<=gec_vx[0]; i0++)
                            energy += 0.5*CellSize/buoyancy[i1][i0]*(vx_tmp._s2)[i1][i0]*vx[i1][i0];
                    ra_destroy(&vx_tmp);
                    
                    //  add vz
                    for(int i1=gsc_vz[1]; i1<=gec_vz[1]; i1++)
                        for(int i0=gsc_vz[0]; i0<=gec_vz[0]; i0++)
                            energy += 0.5*CellSize*4/(buoyancy[i1][i0]+buoyancy[i1][i0+1]+buoyancy[i1+1][i0]+buoyancy[i1+1][i0+1])
                            *(vz_tmp._s2)[i1][i0]*vz[i1][i0];
                    ra_destroy(&vz_tmp);
                    
                    //  add txx,tzz
                    for(int i1=gsc_szz[1]; i1<=gec_szz[1]; i1++){
                        for(int i0=gsc_szz[0]; i0<=gec_szz[0]; i0++){
                            
                            ireal c11_tmp, c33_tmp, c13_tmp;
                            if(esgpars->effective_media){
                                c11_tmp = 2*c11[i1][i0]*c11[i1+1][i0]/(c11[i1][i0]+c11[i1+1][i0]);
                                c33_tmp = 2*c33[i1][i0]*c33[i1+1][i0]/(c33[i1][i0]+c33[i1+1][i0]);
                                c13_tmp = 2*c13[i1][i0]*c13[i1+1][i0]/(c13[i1][i0]+c13[i1+1][i0]);
                            }else{
                                c11_tmp = (c11[i1][i0]+c11[i1+1][i0])/2;
                                c33_tmp = (c33[i1][i0]+c33[i1+1][i0])/2;
                                c13_tmp = (c13[i1][i0]+c13[i1+1][i0])/2;
                            }
                            
                            ireal base = c11_tmp*c33_tmp - c13_tmp*c13_tmp;
                            energy += 0.5*CellSize*(c33_tmp*sxx[i1][i0]*sxx[i1][i0]
                                                    -2*c13_tmp*sxx[i1][i0]*szz[i1][i0]
                                                    +c11_tmp*szz[i1][i0]*szz[i1][i0])/base;
                        }
                    }
                    
                    if(esgpars->lbc[0]){
                        int i0 = gsc_szz[0]-1;
                        for(int i1=gsc_szz[1]; i1<=gec_szz[1]; i1++){
                            
                            ireal c11_tmp, c33_tmp, c13_tmp;
                            if(esgpars->effective_media){
                                c11_tmp = 2*c11[i1][i0]*c11[i1+1][i0]/(c11[i1][i0]+c11[i1+1][i0]);
                                c33_tmp = 2*c33[i1][i0]*c33[i1+1][i0]/(c33[i1][i0]+c33[i1+1][i0]);
                                c13_tmp = 2*c13[i1][i0]*c13[i1+1][i0]/(c13[i1][i0]+c13[i1+1][i0]);
                            }else{
                                c11_tmp = (c11[i1][i0]+c11[i1+1][i0])/2;
                                c33_tmp = (c33[i1][i0]+c33[i1+1][i0])/2;
                                c13_tmp = (c13[i1][i0]+c13[i1+1][i0])/2;
                            }
                            
                            ireal base = c11_tmp*c33_tmp - c13_tmp*c13_tmp;
                            energy += 0.25*CellSize*(c33_tmp*sxx[i1][i0]*sxx[i1][i0]
                                                     -2*c13_tmp*sxx[i1][i0]*szz[i1][i0]
                                                     +c11_tmp*szz[i1][i0]*szz[i1][i0])/base;
                        }
                    }
                    
                    if(esgpars->rbc[0]){
                        int i0 = gec_szz[0]+1;
                        for(int i1=gsc_szz[1]; i1<=gec_szz[1]; i1++){
                            
                            ireal c11_tmp, c33_tmp, c13_tmp;
                            if(esgpars->effective_media){
                                c11_tmp = 2*c11[i1][i0]*c11[i1+1][i0]/(c11[i1][i0]+c11[i1+1][i0]);
                                c33_tmp = 2*c33[i1][i0]*c33[i1+1][i0]/(c33[i1][i0]+c33[i1+1][i0]);
                                c13_tmp = 2*c13[i1][i0]*c13[i1+1][i0]/(c13[i1][i0]+c13[i1+1][i0]);
                            }else{
                                c11_tmp = (c11[i1][i0]+c11[i1+1][i0])/2;
                                c33_tmp = (c33[i1][i0]+c33[i1+1][i0])/2;
                                c13_tmp = (c13[i1][i0]+c13[i1+1][i0])/2;
                            }
                            
                            ireal base = c11_tmp*c33_tmp - c13_tmp*c13_tmp;
                            energy += 0.25*CellSize*(c33_tmp*sxx[i1][i0]*sxx[i1][i0]
                                                     -2*c13_tmp*sxx[i1][i0]*szz[i1][i0]
                                                     +c11_tmp*szz[i1][i0]*szz[i1][i0])/base;
                        }
                    }
                    
                    //  add txz
                    ireal *c55_tmp = (ireal*)malloc((gec_sxz[0]-gsc_sxz[0]+1)*sizeof(ireal));
                    for(int i1=gsc_sxz[1]; i1<=gec_sxz[1]; i1++)
                    {
                        for(int i0=gsc_sxz[0]; i0<=gec_sxz[0]; i0++)
                        {
                            if(esgpars->effective_media)
                            {
                                if(c55[i1][i0] < FLT_MIN || c55[i1][i0+1] < FLT_MIN)
                                    c55_tmp[i0-gsc_sxz[0]] = 0;
                                else
                                    c55_tmp[i0-gsc_sxz[0]] = 2*c55[i1][i0]*c55[i1][i0+1]/(c55[i1][i0]+c55[i1][i0+1]);
                            }
                            else
                                c55_tmp[i0-gsc_sxz[0]] = 0.5*(c55[i1][i0]+c55[i1][i0+1]);
                        }
                        
                        for(int i0=gsc_sxz[0]; i0<=gec_sxz[0]; i0++)
                            energy += 0.5*CellSize*sxz[i1][i0]*sxz[i1][i0]/c55_tmp[i0-gsc_sxz[0]];
                    }
                    
                    if(esgpars->lbc[1]){
                        int i1 = gsc_sxz[1]-1;
                        
                        for(int i0=gsc_sxz[0]; i0<=gec_sxz[0]; i0++)
                        {
                            if(esgpars->effective_media)
                            {
                                if(c55[i1][i0] < FLT_MIN || c55[i1][i0+1] < FLT_MIN)
                                    c55_tmp[i0-gsc_sxz[0]] = 0;
                                else
                                    c55_tmp[i0-gsc_sxz[0]] = 2*c55[i1][i0]*c55[i1][i0+1]/(c55[i1][i0]+c55[i1][i0+1]);
                            }
                            else
                                c55_tmp[i0-gsc_sxz[0]] = 0.5*(c55[i1][i0]+c55[i1][i0+1]);
                        }
                        
                        for(int i0=gsc_sxz[0]; i0<=gec_sxz[0]; i0++)
                            energy += 0.25*CellSize*sxz[i1][i0]*sxz[i1][i0]/c55_tmp[i0-gsc_sxz[0]];
                    }
                    
                    if(esgpars->rbc[1]){
                        int i1 = gec_sxz[1]+1;
                        
                        for(int i0=gsc_sxz[0]; i0<=gec_sxz[0]; i0++)
                        {
                            if(esgpars->effective_media)
                            {
                                if(c55[i1][i0] < FLT_MIN || c55[i1][i0+1] < FLT_MIN)
                                    c55_tmp[i0-gsc_sxz[0]] = 0;
                                else
                                    c55_tmp[i0-gsc_sxz[0]] = 2*c55[i1][i0]*c55[i1][i0+1]/(c55[i1][i0]+c55[i1][i0+1]);
                            }
                            else
                                c55_tmp[i0-gsc_sxz[0]] = 0.5*(c55[i1][i0]+c55[i1][i0+1]);
                        }
                        
                        for(int i0=gsc_sxz[0]; i0<=gec_sxz[0]; i0++)
                            energy += 0.25*CellSize*sxz[i1][i0]*sxz[i1][i0]/c55_tmp[i0-gsc_sxz[0]];
                    }
                    free(c55_tmp);
                    
                    //  record energy
                    float energy_sum=energy;
#ifdef IWAVE_USE_MPI
                    MPI_Comm cm = retrieveComm();
                    MPI_Reduce((void*)(&(energy)),(void*)(&(energy_sum)),1,MPI_FLOAT,MPI_SUM,0,cm);
#endif
                    trace[gsc_trace[1]][gsc_trace[0]] = energy_sum;
                }
            }
            else if(iv==1)
                esg_sstep2d(
                            c11, c33, c13, c55, vx, vz,
                            gsc_szz, gec_szz, gsc_sxz, gec_sxz,
                            esgpars->lbc, esgpars->rbc, esgpars->k, esgpars->coeffs,
                            sxx, szz, sxz, esgpars->eflag, esgpars->effective_media);
        }
        catch (RVLException & e){
            e<<"\nERROR: esg_timestep: simulation step.\n";
            throw e;
        }
    }else{  //  ndim == 3
        
        register ireal*** restrict buoyancy = (dom[0]->_s)[D_BUOY]._s3;
        register ireal*** restrict c11 = (dom[0]->_s)[D_C11]._s3;
        register ireal*** restrict c22 = (dom[0]->_s)[D_C22]._s3;
        register ireal*** restrict c33 = (dom[0]->_s)[D_C33]._s3;
        register ireal*** restrict c23 = (dom[0]->_s)[D_C23]._s3;
        register ireal*** restrict c13 = (dom[0]->_s)[D_C13]._s3;
        register ireal*** restrict c12 = (dom[0]->_s)[D_C12]._s3;
        register ireal*** restrict c44 = (dom[0]->_s)[D_C44]._s3;
        register ireal*** restrict c55 = (dom[0]->_s)[D_C55]._s3;
        register ireal*** restrict c66 = (dom[0]->_s)[D_C66]._s3;
        register ireal*** restrict vx = (dom[0]->_s)[D_VX]._s3;
        register ireal*** restrict vy = (dom[0]->_s)[D_VY]._s3;
        register ireal*** restrict vz = (dom[0]->_s)[D_VZ]._s3;
        register ireal*** restrict sxx = (dom[0]->_s)[D_SXX]._s3;
        register ireal*** restrict syy = (dom[0]->_s)[D_SYY]._s3;
        register ireal*** restrict szz = (dom[0]->_s)[D_SZZ]._s3;
        register ireal*** restrict syz = (dom[0]->_s)[D_SYZ]._s3;
        register ireal*** restrict sxz = (dom[0]->_s)[D_SXZ]._s3;
        register ireal*** restrict sxy = (dom[0]->_s)[D_SXY]._s3;
        register ireal*** restrict trace = (dom[0]->_s)[D_TRACE]._s3;
        
        
        try{
            if(iv == 0){
                
                RARR vz_tmp, vx_tmp, vy_tmp;
                
                ra_setnull(&vz_tmp);
                ra_setnull(&vx_tmp);
                ra_setnull(&vy_tmp);
                
                
                if(esgpars->eflag==1){
                    if (ra_create(&vx_tmp, gsc_vx, gec_vx)) {
                        RVLException e;
                        e<<"\nERROR: esg_timestep: create tmp velocity arrays.\n";
                        throw e;
                    }
                    if (ra_create(&vz_tmp, gsc_vz, gec_vz)) {
                        RVLException e;
                        e<<"\nERROR: esg_timestep: create tmp velocity arrays.\n";
                        throw e;
                    }
                    if (ra_create(&vy_tmp, gsc_vy, gec_vy)) {
                        RVLException e;
                        e<<"\nERROR: esg_timestep: create tmp velocity arrays.\n";
                        throw e;
                    }
                    
                    ra_zero(&vx_tmp);
                    ra_zero(&vz_tmp);
                    ra_zero(&vy_tmp);
                    
                    
                    if (ra_copy(&vz_tmp, &((dom[0]->_s)[D_VZ]))){
                        RVLException e;
                        e<<"\nERROR: esg_timestep: copy velocity to tmp velocity arrays.\n";
                        throw e;
                    }
                    if (ra_copy(&vx_tmp, &((dom[0]->_s)[D_VX]))){
                        RVLException e;
                        e<<"\nERROR: esg_timestep: copy velocity to tmp velocity arrays.\n";
                        throw e;
                    }
                    if (ra_copy(&vy_tmp, &((dom[0]->_s)[D_VY]))){
                        RVLException e;
                        e<<"\nERROR: esg_timestep: copy velocity to tmp velocity arrays.\n";
                        throw e;
                    }
                }
                
                
                // updating velocities
                esg_vstep3d(
                            buoyancy, sxx, syy, szz, syz, sxz, sxy,
                            gsc_vx, gec_vx, gsc_vy, gec_vy, gsc_vz, gec_vz,
                            esgpars->lbc, esgpars->rbc, esgpars->k, esgpars->coeffs,
                            vx, vy, vz);
                
                
                if(esgpars->eflag==1){
                    
                    ireal energy=0;
                    ireal CellSize=esgpars->CellSize;
                    
                    //  add vz
                    for(int i2=gsc_vz[2]; i2<=gec_vz[2]; i2++)
                        for(int i1=gsc_vz[1]; i1<=gec_vz[1]; i1++)
                            for(int i0=gsc_vz[0]; i0<=gec_vz[0]; i0++)
                                energy += 0.5*CellSize*4/(buoyancy[i2][i1][i0]+buoyancy[i2][i1][i0+1]+buoyancy[i2][i1+1][i0]+buoyancy[i2][i1+1][i0+1])
                                *(vz_tmp._s3)[i2][i1][i0]*vz[i2][i1][i0];
                    ra_destroy(&vz_tmp);
                    
                    //  add vx
                    for(int i2=gsc_vx[2]; i2<=gec_vx[2]; i2++)
                        for(int i1=gsc_vx[1]; i1<=gec_vx[1]; i1++)
                            for(int i0=gsc_vx[0]; i0<=gec_vx[0]; i0++)
                                energy += 0.5*CellSize/buoyancy[i2][i1][i0]*(vx_tmp._s3)[i2][i1][i0]*vx[i2][i1][i0];
                    ra_destroy(&vx_tmp);
                    
                    //   add vy
                    for(int i2=gsc_vy[2]; i2<=gec_vy[2]; i2++)
                        for(int i1=gsc_vy[1]; i1<=gec_vy[1]; i1++)
                            for(int i0=gsc_vy[0]; i0<=gec_vy[0]; i0++)
                                energy += 0.5*CellSize*4/(buoyancy[i2][i1][i0]+buoyancy[i2][i1+1][i0]+buoyancy[i2+1][i1][i0]+buoyancy[i2+1][i1+1][i0])
                                *(vy_tmp._s3)[i2][i1][i0]*vy[i2][i1][i0];
                    ra_destroy(&vy_tmp);
                    
                    
                    //   add sxx,syy,szz
                    for(int i2=gsc_szz[2]; i2<=gec_szz[2]; i2++)
                        for(int i1=gsc_szz[1]; i1<=gec_szz[1]; i1++)
                            for(int i0=gsc_szz[0]; i0<=gec_szz[0]; i0++){
                                ireal c11_tmp = (c11[i2][i1][i0] + c11[i2][i1+1][i0])/2;
                                ireal c22_tmp = (c22[i2][i1][i0] + c22[i2][i1+1][i0])/2;
                                ireal c33_tmp = (c33[i2][i1][i0] + c33[i2][i1+1][i0])/2;
                                ireal c23_tmp = (c23[i2][i1][i0] + c23[i2][i1+1][i0])/2;
                                ireal c13_tmp = (c13[i2][i1][i0] + c13[i2][i1+1][i0])/2;
                                ireal c12_tmp = (c12[i2][i1][i0] + c12[i2][i1+1][i0])/2;
                                ireal base = c11_tmp*c23_tmp*c23_tmp + c22_tmp*c13_tmp*c13_tmp + c33_tmp*c12_tmp*c12_tmp
                                - 2*c12_tmp*c13_tmp*c23_tmp - c11_tmp*c22_tmp*c33_tmp;
                                ireal c11_inv = c23_tmp*c23_tmp - c22_tmp*c33_tmp;
                                ireal c22_inv = c13_tmp*c13_tmp - c11_tmp*c33_tmp;
                                ireal c33_inv = c12_tmp*c12_tmp - c11_tmp*c22_tmp;
                                ireal c12_inv = c12_tmp*c33_tmp - c13_tmp*c23_tmp;
                                ireal c13_inv = c13_tmp*c22_tmp - c12_tmp*c23_tmp;
                                ireal c23_inv = c23_tmp*c11_tmp - c12_tmp*c13_tmp;
                                
                                energy += 0.5*CellSize*(c11_inv*sxx[i2][i1][i0]*sxx[i2][i1][i0]+
                                                        c22_inv*syy[i2][i1][i0]*syy[i2][i1][i0]+
                                                        c33_inv*szz[i2][i1][i0]*szz[i2][i1][i0]+
                                                        2*c23_inv*syy[i2][i1][i0]*szz[i2][i1][i0]+
                                                        2*c13_inv*sxx[i2][i1][i0]*szz[i2][i1][i0]+
                                                        2*c12_inv*sxx[i2][i1][i0]*syy[i2][i1][i0])/base;
                                
                            }
                    
                    if(esgpars->lbc[0]){
                        int i0 = gsc_szz[0]-1;
                        for(int i2=gsc_szz[2]; i2<=gec_szz[2]; i2++)
                            for(int i1=gsc_szz[1]; i1<=gec_szz[1]; i1++){
                                ireal c11_tmp = (c11[i2][i1][i0] + c11[i2][i1+1][i0])/2;
                                ireal c22_tmp = (c22[i2][i1][i0] + c22[i2][i1+1][i0])/2;
                                ireal c33_tmp = (c33[i2][i1][i0] + c33[i2][i1+1][i0])/2;
                                ireal c23_tmp = (c23[i2][i1][i0] + c23[i2][i1+1][i0])/2;
                                ireal c13_tmp = (c13[i2][i1][i0] + c13[i2][i1+1][i0])/2;
                                ireal c12_tmp = (c12[i2][i1][i0] + c12[i2][i1+1][i0])/2;
                                ireal base = c11_tmp*c23_tmp*c23_tmp + c22_tmp*c13_tmp*c13_tmp + c33_tmp*c12_tmp*c12_tmp
                                - 2*c12_tmp*c13_tmp*c23_tmp - c11_tmp*c22_tmp*c33_tmp;
                                ireal c11_inv = c23_tmp*c23_tmp - c22_tmp*c33_tmp;
                                ireal c22_inv = c13_tmp*c13_tmp - c11_tmp*c33_tmp;
                                ireal c33_inv = c12_tmp*c12_tmp - c11_tmp*c22_tmp;
                                ireal c12_inv = c12_tmp*c33_tmp - c13_tmp*c23_tmp;
                                ireal c13_inv = c13_tmp*c22_tmp - c12_tmp*c23_tmp;
                                ireal c23_inv = c23_tmp*c11_tmp - c12_tmp*c13_tmp;
                                
                                energy += 0.25*CellSize*(c11_inv*sxx[i2][i1][i0]*sxx[i2][i1][i0]+
                                                         c22_inv*syy[i2][i1][i0]*syy[i2][i1][i0]+
                                                         c33_inv*szz[i2][i1][i0]*szz[i2][i1][i0]+
                                                         2*c23_inv*syy[i2][i1][i0]*szz[i2][i1][i0]+
                                                         2*c13_inv*sxx[i2][i1][i0]*szz[i2][i1][i0]+
                                                         2*c12_inv*sxx[i2][i1][i0]*syy[i2][i1][i0])/base;
                                
                            }
                    }
                    
                    if(esgpars->rbc[0]){
                        int i0 = gec_szz[0]+1;
                        for(int i2=gsc_szz[2]; i2<=gec_szz[2]; i2++)
                            for(int i1=gsc_szz[1]; i1<=gec_szz[1]; i1++){
                                ireal c11_tmp = (c11[i2][i1][i0] + c11[i2][i1+1][i0])/2;
                                ireal c22_tmp = (c22[i2][i1][i0] + c22[i2][i1+1][i0])/2;
                                ireal c33_tmp = (c33[i2][i1][i0] + c33[i2][i1+1][i0])/2;
                                ireal c23_tmp = (c23[i2][i1][i0] + c23[i2][i1+1][i0])/2;
                                ireal c13_tmp = (c13[i2][i1][i0] + c13[i2][i1+1][i0])/2;
                                ireal c12_tmp = (c12[i2][i1][i0] + c12[i2][i1+1][i0])/2;
                                ireal base = c11_tmp*c23_tmp*c23_tmp + c22_tmp*c13_tmp*c13_tmp + c33_tmp*c12_tmp*c12_tmp
                                - 2*c12_tmp*c13_tmp*c23_tmp - c11_tmp*c22_tmp*c33_tmp;
                                ireal c11_inv = c23_tmp*c23_tmp - c22_tmp*c33_tmp;
                                ireal c22_inv = c13_tmp*c13_tmp - c11_tmp*c33_tmp;
                                ireal c33_inv = c12_tmp*c12_tmp - c11_tmp*c22_tmp;
                                ireal c12_inv = c12_tmp*c33_tmp - c13_tmp*c23_tmp;
                                ireal c13_inv = c13_tmp*c22_tmp - c12_tmp*c23_tmp;
                                ireal c23_inv = c23_tmp*c11_tmp - c12_tmp*c13_tmp;
                                
                                energy += 0.25*CellSize*(c11_inv*sxx[i2][i1][i0]*sxx[i2][i1][i0]+
                                                         c22_inv*syy[i2][i1][i0]*syy[i2][i1][i0]+
                                                         c33_inv*szz[i2][i1][i0]*szz[i2][i1][i0]+
                                                         2*c23_inv*syy[i2][i1][i0]*szz[i2][i1][i0]+
                                                         2*c13_inv*sxx[i2][i1][i0]*szz[i2][i1][i0]+
                                                         2*c12_inv*sxx[i2][i1][i0]*syy[i2][i1][i0])/base;
                                
                            }
                    }
                    
                    if(esgpars->lbc[2]){
                        int i2 = gsc_szz[2]-1;
                        for(int i1 = gsc_szz[1]; i1<=gec_szz[1]; i1++)
                            for(int i0 = gsc_szz[0]; i0<=gec_szz[0]; i0++){
                                ireal c11_tmp = (c11[i2][i1][i0] + c11[i2][i1+1][i0])/2;
                                ireal c22_tmp = (c22[i2][i1][i0] + c22[i2][i1+1][i0])/2;
                                ireal c33_tmp = (c33[i2][i1][i0] + c33[i2][i1+1][i0])/2;
                                ireal c23_tmp = (c23[i2][i1][i0] + c23[i2][i1+1][i0])/2;
                                ireal c13_tmp = (c13[i2][i1][i0] + c13[i2][i1+1][i0])/2;
                                ireal c12_tmp = (c12[i2][i1][i0] + c12[i2][i1+1][i0])/2;
                                ireal base = c11_tmp*c23_tmp*c23_tmp + c22_tmp*c13_tmp*c13_tmp + c33_tmp*c12_tmp*c12_tmp
                                - 2*c12_tmp*c13_tmp*c23_tmp - c11_tmp*c22_tmp*c33_tmp;
                                ireal c11_inv = c23_tmp*c23_tmp - c22_tmp*c33_tmp;
                                ireal c22_inv = c13_tmp*c13_tmp - c11_tmp*c33_tmp;
                                ireal c33_inv = c12_tmp*c12_tmp - c11_tmp*c22_tmp;
                                ireal c12_inv = c12_tmp*c33_tmp - c13_tmp*c23_tmp;
                                ireal c13_inv = c13_tmp*c22_tmp - c12_tmp*c23_tmp;
                                ireal c23_inv = c23_tmp*c11_tmp - c12_tmp*c13_tmp;
                                
                                energy += 0.25*CellSize*(c11_inv*sxx[i2][i1][i0]*sxx[i2][i1][i0]+
                                                         c22_inv*syy[i2][i1][i0]*syy[i2][i1][i0]+
                                                         c33_inv*szz[i2][i1][i0]*szz[i2][i1][i0]+
                                                         2*c23_inv*syy[i2][i1][i0]*szz[i2][i1][i0]+
                                                         2*c13_inv*sxx[i2][i1][i0]*szz[i2][i1][i0]+
                                                         2*c12_inv*sxx[i2][i1][i0]*syy[i2][i1][i0])/base;
                                
                            }
                    }
                    
                    if(esgpars->rbc[2]){
                        int i2 = gec_szz[2]+1;
                        for(int i1 = gsc_szz[1]; i1<=gec_szz[1]; i1++)
                            for(int i0 = gsc_szz[0]; i0<=gec_szz[0]; i0++){
                                ireal c11_tmp = (c11[i2][i1][i0] + c11[i2][i1+1][i0])/2;
                                ireal c22_tmp = (c22[i2][i1][i0] + c22[i2][i1+1][i0])/2;
                                ireal c33_tmp = (c33[i2][i1][i0] + c33[i2][i1+1][i0])/2;
                                ireal c23_tmp = (c23[i2][i1][i0] + c23[i2][i1+1][i0])/2;
                                ireal c13_tmp = (c13[i2][i1][i0] + c13[i2][i1+1][i0])/2;
                                ireal c12_tmp = (c12[i2][i1][i0] + c12[i2][i1+1][i0])/2;
                                ireal base = c11_tmp*c23_tmp*c23_tmp + c22_tmp*c13_tmp*c13_tmp + c33_tmp*c12_tmp*c12_tmp
                                - 2*c12_tmp*c13_tmp*c23_tmp - c11_tmp*c22_tmp*c33_tmp;
                                ireal c11_inv = c23_tmp*c23_tmp - c22_tmp*c33_tmp;
                                ireal c22_inv = c13_tmp*c13_tmp - c11_tmp*c33_tmp;
                                ireal c33_inv = c12_tmp*c12_tmp - c11_tmp*c22_tmp;
                                ireal c12_inv = c12_tmp*c33_tmp - c13_tmp*c23_tmp;
                                ireal c13_inv = c13_tmp*c22_tmp - c12_tmp*c23_tmp;
                                ireal c23_inv = c23_tmp*c11_tmp - c12_tmp*c13_tmp;
                                
                                energy += 0.25*CellSize*(c11_inv*sxx[i2][i1][i0]*sxx[i2][i1][i0]+
                                                         c22_inv*syy[i2][i1][i0]*syy[i2][i1][i0]+
                                                         c33_inv*szz[i2][i1][i0]*szz[i2][i1][i0]+
                                                         2*c23_inv*syy[i2][i1][i0]*szz[i2][i1][i0]+
                                                         2*c13_inv*sxx[i2][i1][i0]*szz[i2][i1][i0]+
                                                         2*c12_inv*sxx[i2][i1][i0]*syy[i2][i1][i0])/base;
                                
                            }
                    }
                    
                    // add syz
                    for(int i2=gsc_syz[2]; i2<=gec_syz[2]; i2++)
                        for(int i1=gsc_syz[1]; i1<=gec_syz[1]; i1++)
                            for(int i0=gsc_syz[0]; i0<=gec_syz[0]; i0++)
                                energy+= CellSize*4/(c44[i2][i1][i0] + c44[i2][i1][i0+1]
                                                     + c44[i2][i1+1][i0] + c44[i2][i1+1][i0+1]
                                                     + c44[i2+1][i1][i0] + c44[i2+1][i1][i0+1]
                                                     + c44[i2+1][i1+1][i0] + c44[i2+1][i1+1][i0+1])
                                *syz[i2][i1][i0]*syz[i2][i1][i0];
                    
                    // add sxz, sxy
                    for(int i2=gsc_sxz[2]; i2<=gec_sxz[2]; i2++)
                        for(int i1=gsc_sxz[1]; i1<=gec_sxz[1]; i1++)
                            for(int i0=gsc_sxz[0]; i0<=gec_sxz[0]; i0++)
                                energy+= CellSize/(c55[i2][i1][i0] + c55[i2][i1][i0+1])
                                *sxz[i2][i1][i0]*sxz[i2][i1][i0];
                    
                    for(int i2=gsc_sxy[2]; i2<=gec_sxy[2]; i2++)
                        for(int i1=gsc_sxy[1]; i1<=gec_sxy[1]; i1++)
                            for(int i0=gsc_sxy[0]; i0<=gec_sxy[0]; i0++)
                                energy+= CellSize/(c66[i2][i1][i0]+c66[i2+1][i1][i0])
                                *sxy[i2][i1][i0]*sxy[i2][i1][i0];
                    
                    if(esgpars->lbc[1]){
                        
                        int i1 = gsc_sxz[1]-1;
                        
                        for(int i2=gsc_sxz[2]; i2<=gec_sxz[2]; i2++)
                            for(int i0=gsc_sxz[0]; i0<=gec_sxz[0]; i0++)
                                energy+= 0.5*CellSize/(c55[i2][i1][i0] + c55[i2][i1][i0+1])
                                *sxz[i2][i1][i0]*sxz[i2][i1][i0];
                        
                        for(int i2=gsc_sxy[2]; i2<=gec_sxy[2]; i2++)
                            for(int i0=gsc_sxy[0]; i0<=gec_sxy[0]; i0++)
                                energy+= 0.5*CellSize/(c66[i2][i1][i0]+c66[i2+1][i1][i0])
                                *sxy[i2][i1][i0]*sxy[i2][i1][i0];
                    }
                    
                    if(esgpars->rbc[1]){
                        
                        int i1 = gec_sxz[1]+1;
                        
                        for(int i2=gsc_sxz[2]; i2<=gec_sxz[2]; i2++)
                            for(int i0=gsc_sxz[0]; i0<=gec_sxz[0]; i0++)
                                energy+= 0.5*CellSize/(c55[i2][i1][i0] + c55[i2][i1][i0+1])
                                *sxz[i2][i1][i0]*sxz[i2][i1][i0];
                        
                        for(int i2=gsc_sxy[2]; i2<=gec_sxy[2]; i2++)
                            for(int i0=gsc_sxy[0]; i0<=gec_sxy[0]; i0++)
                                energy+= 0.5*CellSize/(c66[i2][i1][i0]+c66[i2+1][i1][i0])
                                *sxy[i2][i1][i0]*sxy[i2][i1][i0];
                    }
                    
                    
                    
                    float energy_sum=energy;
#ifdef IWAVE_USE_MPI
                    MPI_Comm cm = retrieveComm();
                    MPI_Reduce((void*)(&(energy)),(void*)(&(energy_sum)),1,MPI_FLOAT,MPI_SUM,0,cm);
#endif
                    
                    trace[gsc_trace[2]][gsc_trace[1]][gsc_trace[0]] = energy_sum;
                    
                } //End: esgpars->flag==1
                
                
            }else{  //  iv==1
                esg_sstep3d(
                            c11, c22, c33, c23, c13, c12, c44, c55, c66, vx, vy, vz,
                            gsc_szz, gec_szz, gsc_syz, gec_syz, gsc_sxz, gec_sxz, gsc_sxy, gec_sxy,
                            esgpars->lbc, esgpars->rbc, esgpars->k, esgpars->coeffs,
                            sxx, syy, szz, syz, sxz, sxy, esgpars->eflag);
            }
        }catch(RVLException &e){
            e<<"\nERROR: esg_timestep: simulation step.\n";
            throw e;
        }
    }
    
}









//----------FUNCTION: compute stencil mask which is used for domain decomposition
//----------INPUT:    specs, ndim, gtype(not used in this function)
//----------OUTPUT:   sten, stream
int esg_create_sten(void * specs,
                    FILE * stream,
                    int ndim,
                    IPNT gtype[RDOM_MAX_NARR],
                    STENCIL * sten) {
    
    ESG_TS_PARS* esgpars = (ESG_TS_PARS*)specs;
    
    //  ===> Set number of masks (computed from below)
    int nmask = 27;
    if(esgpars->ndim == 3)
        nmask = 97;
    
    //  ===> Initialize sten
    sten_setnull(sten);
    int err;
    if((err = sten_create(sten, nmask))){
        fprintf(stream, "\nERROR:   esg_create_sten: failed to allocate stencil");
        return err;
    }
    
    //  ===> Set sten
    int len = 2*esgpars->k;
    int imask=-1;
    STENCIL_MASK mask;
    
    
    //  1.  sxz,sxx,sxy -> vx   +2/3
    int D_SZXY_VX[] = {D_SXZ, D_SXX, D_SXY};
    for(int idim = 0; idim < esgpars->ndim; idim++){
        imask++;
        if ((err=mask_create(&mask, D_SZXY_VX[idim], D_VX, len))) {
            fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
            sten_destroy(sten);
            return err;
        }
        for (int i=0; i<len; i++) {
            IPNT offs;
            IASN(offs,IPNT_0);
            offs[idim]=i-esgpars->k;
            if ((err=mask_set(&mask,i,offs))) {
                fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                sten_destroy(sten);
                return err;
            }
        }
        if ((err=sten_set(sten,imask, &mask))) {
            fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
            sten_destroy(sten);
            return err;
        }
    }
    
    
    //  2.  szz,sxz,syz -> vz   +2/3
    int D_SZXY_VZ[] = {D_SZZ, D_SXZ, D_SYZ};
    for(int idim = 0; idim < esgpars->ndim; idim++){
        imask++;
        if ((err=mask_create(&mask, D_SZXY_VZ[idim], D_VZ, len))) {
            fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
            sten_destroy(sten);
            return err;
        }
        for (int i=0; i<len; i++) {
            IPNT offs;
            IASN(offs,IPNT_0);
            offs[idim]=i-esgpars->k+1-idim/2;
            if ((err=mask_set(&mask,i,offs))) {
                fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                sten_destroy(sten);
                return err;
            }
        }
        if ((err=sten_set(sten,imask, &mask))) {
            fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
            sten_destroy(sten);
            return err;
        }
    }
    
    
    //  3.  buoy -> vz  +2/3
    //  along x,z direction
    for(int idim = 0; idim < 2; idim++){
        imask++;
        if ((err=mask_create(&mask, D_BUOY, D_VZ, 2))) {
            fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
            sten_destroy(sten);
            return err;
        }
        for (int i=0; i<2; i++) {
            IPNT offs;
            IASN(offs,IPNT_0);
            offs[idim]=i;
            if ((err=mask_set(&mask,i,offs))) {
                fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                sten_destroy(sten);
                return err;
            }
        }
        if ((err=sten_set(sten,imask, &mask))) {
            fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
            sten_destroy(sten);
            return err;
        }
    }
    
    //  along y direction
    if(esgpars->ndim == 3){
        imask++;
        if ((err=mask_create(&mask, D_BUOY, D_VZ, 3))) {
            fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
            sten_destroy(sten);
            return err;
        }
        for (int i=0; i<3; i++) {
            IPNT offs;
            IASN(offs,IPNT_0);
            offs[2]=i-1;
            if ((err=mask_set(&mask,i,offs))) {
                fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                sten_destroy(sten);
                return err;
            }
        }
        if ((err=sten_set(sten,imask, &mask))) {
            fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
            sten_destroy(sten);
            return err;
        }
    }
    
    
    //  4.  syz, sxy, syy, buoy -> vy   +0/6
    if(esgpars->ndim == 3){
        
        //  4.1 syz, sxy, syy -> vy
        int D_SZXY_VY[] = {D_SYZ, D_SXY, D_SYY};
        for(int idim = 0; idim < 3; idim++){
            imask++;
            if ((err=mask_create(&mask, D_SZXY_VY[idim], D_VY, len))) {
                fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
                sten_destroy(sten);
                return err;
            }
            for (int i=0; i<len; i++) {
                IPNT offs;
                IASN(offs,IPNT_0);
                offs[idim]=i-esgpars->k+(idim!=0);
                if ((err=mask_set(&mask,i,offs))) {
                    fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                    sten_destroy(sten);
                    return err;
                }
            }
            if ((err=sten_set(sten,imask, &mask))) {
                fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
                sten_destroy(sten);
                return err;
            }
        }
        
        //  4.2 buoy -> vy
        //  Along x,y direction
        for(int idim = 1; idim < 3; idim++){
            imask++;
            if ((err=mask_create(&mask, D_BUOY, D_VY, 2))) {
                fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
                sten_destroy(sten);
                return err;
            }
            for (int i=0; i<2; i++) {
                IPNT offs;
                IASN(offs,IPNT_0);
                offs[idim]=i;
                if ((err=mask_set(&mask,i,offs))) {
                    fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                    sten_destroy(sten);
                    return err;
                }
            }
            if ((err=sten_set(sten,imask, &mask))) {
                fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
                sten_destroy(sten);
                return err;
            }
        }
        //  Along z direction
        imask++;
        if ((err=mask_create(&mask, D_BUOY, D_VY, 3))) {
            fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
            sten_destroy(sten);
            return err;
        }
        for (int i=0; i<3; i++) {
            IPNT offs;
            IASN(offs,IPNT_0);
            offs[0]=i-1;
            if ((err=mask_set(&mask,i,offs))) {
                fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                sten_destroy(sten);
                return err;
            }
        }
        if ((err=sten_set(sten,imask, &mask))) {
            fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
            sten_destroy(sten);
            return err;
        }
    }
    
    //  5.  Create a fake sxx,syy,szz->vz,vy mask so they get extended to the boundary,
    //  so the program can compute the energy of the solution on the domain +1/4
    if(esgpars->eflag==1){
        //  5.1 sxx -> vz, vy   +1/2
        int D_SXX_VZY[] = {D_VZ, D_VY};
        for(int idim = 0; idim < esgpars->ndim-1; idim++){
            imask++;
            if ((err=mask_create(&mask, D_SXX, D_SXX_VZY[idim], 2))) {
                fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
                sten_destroy(sten);
                return err;
            }
            
            IPNT offs;
            IASN(offs,IPNT_0);
            offs[idim*2]=0;
            if ((err=mask_set(&mask,0,offs))) {
                fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                sten_destroy(sten);
                return err;
            }
            IASN(offs,IPNT_0);
            offs[idim*2]=1;
            if ((err=mask_set(&mask,1,offs))) {
                fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                sten_destroy(sten);
                return err;
            }
            
            if ((err=sten_set(sten,imask, &mask))) {
                fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
                sten_destroy(sten);
                return err;
            }
        }
        //  5.1 syy -> vz, szz -> vy    +0/2
        if(esgpars->ndim == 3){
            int D_SYZ2[] = {D_SYY, D_SZZ};
            int D_VZY[] = {D_VZ, D_VY};
            for(int idim = 0; idim < 2; idim++){
                imask++;
                if ((err=mask_create(&mask, D_SYZ2[idim], D_VZY[idim], 2))) {
                    fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
                    sten_destroy(sten);
                    return err;
                }
                
                IPNT offs;
                IASN(offs,IPNT_0);
                offs[idim*2]=0;
                if ((err=mask_set(&mask,0,offs))) {
                    fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                    sten_destroy(sten);
                    return err;
                }
                IASN(offs,IPNT_0);
                offs[idim*2]=1;
                if ((err=mask_set(&mask,1,offs))) {
                    fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                    sten_destroy(sten);
                    return err;
                }
                
                if ((err=sten_set(sten,imask, &mask))) {
                    fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
                    sten_destroy(sten);
                    return err;
                }
            }
        }
    }
    
    
    //  6.  vz, vx, vy -> sxx, syy, szz +4/9
    int D_VZXY[] = {D_VZ, D_VX, D_VY};
    int D_SZXY[] = {D_SZZ, D_SXX, D_SYY};
    for(int idim = 0; idim < esgpars->ndim; idim++)
        for(int j = 0; j < esgpars->ndim; j++){
            imask++;
            if ((err=mask_create(&mask, D_VZXY[idim], D_SZXY[j], len))) {
                fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
                sten_destroy(sten);
                return err;
            }
            for (int i=0; i<len; i++) {
                IPNT offs;
                IASN(offs,IPNT_0);
                offs[idim]=i-esgpars->k+(idim==1);
                if ((err=mask_set(&mask,i,offs))) {
                    fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                    sten_destroy(sten);
                    return err;
                }
            }
            if ((err=sten_set(sten,imask, &mask))) {
                fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
                sten_destroy(sten);
                return err;
            }
        }
    
    
    //  7.  c11, c33, c13, c22, c23, c12 -> sxx, syy, szz   +12/54
    int D_C[] = {D_C11, D_C33, D_C13, D_C22, D_C23, D_C12};
    int D_SZXY2[] = {D_SZZ, D_SXX, D_SYY};
    for(int i1 = 0; i1 < (esgpars->ndim-1)*3; i1++)
        for(int i2 = 0; i2 < esgpars->ndim; i2++){
            //  Along x direction   +6/18
            imask++;
            if ((err=mask_create(&mask, D_C[i1], D_SZXY2[i2], 2))) {
                fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
                sten_destroy(sten);
                return err;
            }
            for (int i=0; i<2; i++) {
                IPNT offs;
                IASN(offs,IPNT_0);
                offs[1]=i;
                if ((err=mask_set(&mask,i,offs))) {
                    fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                    sten_destroy(sten);
                    return err;
                }
            }
            if ((err=sten_set(sten,imask, &mask))) {
                fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
                sten_destroy(sten);
                return err;
            }
            //  Along z and y direction +6/36
            for(int idim = 0; idim < esgpars->ndim-1; idim++){
                imask++;
                if ((err=mask_create(&mask, D_C[i1], D_SZXY2[i2], 3))) {
                    fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
                    sten_destroy(sten);
                    return err;
                }
                for (int i=0; i<3; i++) {
                    IPNT offs;
                    IASN(offs,IPNT_0);
                    offs[idim*2]=i-1;
                    if ((err=mask_set(&mask,i,offs))) {
                        fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                        sten_destroy(sten);
                        return err;
                    }
                }
                if ((err=sten_set(sten,imask, &mask))) {
                    fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
                    sten_destroy(sten);
                    return err;
                }
            }
        }
    
    
    //  8.  c44, vy, vz -> syz  +0/5
    if(esgpars->ndim == 3){
        for(int idim = 0; idim < 3; idim++){
            imask++;
            if ((err=mask_create(&mask, D_C44, D_SYZ, 2))) {
                fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
                sten_destroy(sten);
                return err;
            }
            for (int i=0; i<2; i++) {
                IPNT offs;
                IASN(offs,IPNT_0);
                offs[idim]=i;
                if ((err=mask_set(&mask,i,offs))) {
                    fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                    sten_destroy(sten);
                    return err;
                }
            }
            if ((err=sten_set(sten,imask, &mask))) {
                fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
                sten_destroy(sten);
                return err;
            }
        }
        
        
        int D_VYZ[] = {D_VY, D_VZ};
        for(int idim = 0; idim < 2; idim++){
            imask++;
            if ((err=mask_create(&mask, D_VYZ[idim], D_SYZ, len))) {
                fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
                sten_destroy(sten);
                return err;
            }
            for (int i=0; i<len; i++) {
                IPNT offs;
                IASN(offs,IPNT_0);
                offs[idim*2]=i-esgpars->k+1;
                if ((err=mask_set(&mask,i,offs))) {
                    fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                    sten_destroy(sten);
                    return err;
                }
            }
            if ((err=sten_set(sten,imask, &mask))) {
                fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
                sten_destroy(sten);
                return err;
            }
        }
    }
    
    
    //  9.  c55, vx, vz -> sxz  +4/5
    imask++;
    if ((err=mask_create(&mask, D_C55, D_SXZ, 2))) {
        fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
        sten_destroy(sten);
        return err;
    }
    for (int i=0; i<2; i++) {
        IPNT offs;
        IASN(offs,IPNT_0);
        offs[0]=i;
        if ((err=mask_set(&mask,i,offs))) {
            fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
            sten_destroy(sten);
            return err;
        }
    }
    if ((err=sten_set(sten,imask, &mask))) {
        fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
        sten_destroy(sten);
        return err;
    }
    
    //  extend c55 to the boundary
    for(int idim = 1; idim < esgpars->ndim; idim++){
        imask++;
        if ((err=mask_create(&mask, D_C55, D_SXZ, 3))) {
            fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
            sten_destroy(sten);
            return err;
        }
        for (int i=0; i<3; i++) {
            IPNT offs;
            IASN(offs,IPNT_0);
            offs[idim]=i-1;
            if ((err=mask_set(&mask,i,offs))) {
                fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                sten_destroy(sten);
                return err;
            }
        }
        if ((err=sten_set(sten,imask, &mask))) {
            fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
            sten_destroy(sten);
            return err;
        }
    }
    
    int D_VXZ[] = {D_VX, D_VZ};
    for(int idim = 0; idim < 2; idim++){
        imask++;
        if ((err=mask_create(&mask, D_VXZ[idim], D_SXZ, len))) {
            fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
            sten_destroy(sten);
            return err;
        }
        for (int i=0; i<len; i++) {
            IPNT offs;
            IASN(offs,IPNT_0);
            offs[idim]=i-esgpars->k+(idim==0);
            if ((err=mask_set(&mask,i,offs))) {
                fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                sten_destroy(sten);
                return err;
            }
        }
        if ((err=sten_set(sten,imask, &mask))) {
            fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
            sten_destroy(sten);
            return err;
        }
    }
    
    
    //  10.  c66, vx, vy -> sxy +0/5
    
    if(esgpars->ndim == 3){
        imask++;
        if ((err=mask_create(&mask, D_C66, D_SXY, 2))) {
            fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
            sten_destroy(sten);
            return err;
        }
        for (int i=0; i<2; i++) {
            IPNT offs;
            IASN(offs,IPNT_0);
            offs[2]=i;
            if ((err=mask_set(&mask,i,offs))) {
                fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                sten_destroy(sten);
                return err;
            }
        }
        if ((err=sten_set(sten,imask, &mask))) {
            fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
            sten_destroy(sten);
            return err;
        }
        
        //  extend c66 to the boundary
        for(int idim = 0; idim < 2; idim++){
            imask++;
            if ((err=mask_create(&mask, D_C66, D_SXY, 3))) {
                fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
                sten_destroy(sten);
                return err;
            }
            for (int i=0; i<3; i++) {
                IPNT offs;
                IASN(offs,IPNT_0);
                offs[idim]=i-1;
                if ((err=mask_set(&mask,i,offs))) {
                    fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                    sten_destroy(sten);
                    return err;
                }
            }
            if ((err=sten_set(sten,imask, &mask))) {
                fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
                sten_destroy(sten);
                return err;
            }
        }
        
        int D_VYX[] = {D_VY, D_VX};
        for(int idim = 0; idim < 2; idim++){
            imask++;
            if ((err=mask_create(&mask, D_VYX[idim], D_SXY, len))) {
                fprintf(stream,"\nERROR:    esg_create_sten: mask_create\n");
                sten_destroy(sten);
                return err;
            }
            for (int i=0; i<len; i++) {
                IPNT offs;
                IASN(offs,IPNT_0);
                offs[idim+1]=i-esgpars->k+(idim==1);
                if ((err=mask_set(&mask,i,offs))) {
                    fprintf(stream,"\nERROR:    esg_create_sten from mask_set\n");
                    sten_destroy(sten);
                    return err;
                }
            }
            if ((err=sten_set(sten,imask, &mask))) {
                fprintf(stream,"\nERROR:    esg_create_sten from sten_set\n");
                sten_destroy(sten);
                return err;
            }
        }
    }
    
    return 0;
}




//----------FUNCTION:
//----------INPUT:    dom, specs
//----------OUTPUT:   stream
void esg_check(RDOM * dom,
               void * specs,
               FILE * stream){
    
    ESG_TS_PARS* esgpars = (ESG_TS_PARS*)specs;
    
    int ndim;
    ra_ndim(&(dom->_s[D_BUOY]),&ndim);
    
    
    if(ndim < 2 || ndim > 3){
        RVLException e;
        e<<"\nError:    esg_check: can't handle ndim = "<<ndim<<".\n";
        throw e;
    }
    
}

