/* 
exchange.c
Igor Terentyev.
*/

//#define EXC_VERBOSE

/*============================================================================*/

#include "exchange.h"

/*----------------------------------------------------------------------------*/

int ex_compute(int iarr, STENCIL *sten, RDOM *dom,
               IPNT a_gs, IPNT a_ge, IPNT r_gs[], IPNT r_ge[],
	       int frcvempty[], FILE * stream)
{
    int ndim, idim, nmask, imask, nrcv, ircv, n, i;
    IPNT rcvpnt, ind, gsr, ger, gs, ge; /* receive area number, mask element */
    STENCIL_MASK mask;  /* current mask */
    int fempty, fstasn; /* empty computational array, first assignment flags */
    int frempty;        /* empty recomputed array in the mask */
    int fake[RARR_MAX_3NDIM - 1];
    int bdrs[RARR_MAX_NDIM][4]; /* boundaries */

    /* dump input */
#ifdef EXC_VERBOSE
    fprintf(stream,"ex_compute: iarr=%d\n",iarr);
    fflush(stream);
#endif
    /* prepare ---------------------------------------------------------------*/
    if ( rd_ndim(dom, iarr, &ndim) ) return E_BADINPUT; /* get number of dims */
    if ( gen_3n1(ndim, &nrcv) ) return E_BADINPUT; /* number of receives */
    if ( rd_empty(dom, iarr, &fempty) ) return E_INTERNAL; /* set empty flag */
    if ( !fempty ) if ( rd_gse(dom, iarr, a_gs, a_ge) ) return E_INTERNAL; /* set alloc */
    if ( frcvempty == NULL ) frcvempty = fake;

#ifdef EXC_VERBOSE
    fprintf(stream,"ex_compute: initial a_gs, a_ge\n");
    for (int idim=0;idim<ndim; idim++) 
      fprintf(stream,"  idim=%d a_gs=%d a_ge=%d\n",idim,a_gs[idim],a_ge[idim]);
#endif
    /* make receives empty ---------------------------------------------------*/
    for ( ircv = 0; ircv < nrcv; ++ircv )     /* loop through receives   */
    {
        if ( gen_i2pnt(ndim, ircv, rcvpnt) ) return E_INTERNAL;
        for ( idim = 0; idim < ndim; ++idim ) /* loop through dimensions */
        {
            if ( rcvpnt[idim] == 0 )                    /* center */
            {
                r_gs[ircv][idim] = a_gs[idim];
                r_ge[ircv][idim] = a_ge[idim]; 
            }
            else if ( rcvpnt[idim] < 0 )                /* left  */
            {
                r_ge[ircv][idim] = a_gs[idim];
                r_gs[ircv][idim] = r_ge[ircv][idim] + 1; 
            }
            else /*if ( rcvpnt[idim] > 0 )*/            /* right */
            {
                r_gs[ircv][idim] = a_ge[idim];
                r_ge[ircv][idim] = r_gs[ircv][idim] - 1; 
            }
        }
        frcvempty[ircv] = 1;
    }
    
    /* loop over masks -------------------------------------------------------*/
    fstasn = 1; /* first assignment (used if empty computational array) */
    nmask = sten->nmask;
    
    for ( imask = 0; imask < nmask; ++imask )
    {
        if ( sten_get(sten, imask, &mask) ) return E_INTERNAL;
        if ( mask.ip != iarr ) continue;  /* mask for another array */
        /* get recomputed array global coordinates, skip if empty */
        if ( rd_empty(dom, mask.ir, &frempty) ) return E_BADINPUT;
        if ( frempty ) continue;
        if ( rd_gse(dom, mask.ir, gsr, ger) ) return E_INTERNAL;
#ifdef EXC_VERBOSE
	fprintf(stream,"ex_compute: imask=%d\n",imask);
	fprintf(stream,"  index of this array = mask.ip=%d\n",mask.ip);
	fprintf(stream,"  index of influenced array = mask.ir=%d\n",mask.ir);
	fprintf(stream,"  gs, gr of influenced array:\n");
	for ( idim = 0; idim < ndim; ++idim )	
	  fprintf(stream,"  idim=%d gsr=%d ger=%d\n",idim,gsr[idim],ger[idim]);
#endif         
        
        /* loop over mask elements */
        n = mask.n;
        for ( i = 0; i < n; ++i )
        {
            /* get mask element */
            if ( mask_get(&mask, i, ind) ) return E_INTERNAL;

            /* compute participating array coordinates */
            for ( idim = 0; idim < ndim; ++idim )
            {
                gs[idim] = gsr[idim] + ind[idim];
                ge[idim] = ger[idim] + ind[idim];		
            }
#ifdef EXC_VERBOSE
	    fprintf(stream,"ex_compute: gs, ge from array %d and mask element %d\n",mask.ir,i);
	    for (int idim=0;idim<ndim; idim++) 
	      fprintf(stream,"  idim=%d ind=%d gs=%d ge=%d\n",idim,ind[idim],gs[idim],ge[idim]);
#endif            
            /* empty computational array case */
            if ( fempty )
            {
                if ( fstasn ) /* first assignment */
                {
                    IASN(a_gs, gs);
                    IASN(a_ge, ge);
                    fstasn = 0;
                }
                else
                {
                    for ( idim = 0; idim < ndim; ++idim )
                    {
		      a_gs[idim] = iwave_min(a_gs[idim], gs[idim]);
		      a_ge[idim] = iwave_max(a_ge[idim], ge[idim]);
                    }
                }
#ifdef EXC_VERBOSE
		fprintf(stream,"ex_compute: empty, new computational array bounds\n");
		for (int idim=0;idim<ndim; idim++) 
		  fprintf(stream,"  idim=%d a_gs=%d a_ge=%d\n",idim,a_gs[idim],a_ge[idim]);
#endif            
            }
            /* non-empty computational array case */
            else
            {
                /* compute boundaries - 6 cases per dimension */
                for ( idim = 0; idim < ndim; ++idim )
                {
                    if      ( (gs[idim] <  a_gs[idim]) && (ge[idim] >  a_ge[idim]) ) /* [X[X[X[ */
                    {
#ifdef EXC_VERBOSE
		      fprintf(stream,"  case 1\n");
#endif            		    

		      bdrs[idim][0] = gs[idim];
		      bdrs[idim][1] = a_gs[idim];
		      bdrs[idim][2] = a_ge[idim] + 1;
		      bdrs[idim][3] = ge[idim] + 1;
                    }
                    else if ( (gs[idim] <  a_gs[idim]) && (ge[idim] >= a_gs[idim]) ) /* [X[X[ [ */
                    {
#ifdef EXC_VERBOSE
		      fprintf(stream,"  case 2\n");
#endif            		    
                        bdrs[idim][0] = gs[idim];
                        bdrs[idim][1] = a_gs[idim];
                        bdrs[idim][2] = bdrs[idim][3] = ge[idim] + 1;
                    }
                    else if ( (gs[idim] <  a_gs[idim]) && (ge[idim] <  a_gs[idim]) ) /* [X[ [ [ */
                    {
#ifdef EXC_VERBOSE
		      fprintf(stream,"  case 3\n");
#endif            		    
                        bdrs[idim][0] = gs[idim];
                        bdrs[idim][1] = bdrs[idim][2] = bdrs[idim][3] = ge[idim] + 1;
                    }
                    else if ( (gs[idim] <= a_ge[idim]) && (ge[idim] >  a_ge[idim]) ) /* [ [X[X[ */
                    {
#ifdef EXC_VERBOSE
		      fprintf(stream,"  case 4\n");
#endif            		    
                        bdrs[idim][0] = bdrs[idim][1] = gs[idim];
                        bdrs[idim][2] = a_ge[idim] + 1;
                        bdrs[idim][3] = ge[idim] + 1;
                    }
                    else if ( (gs[idim] <= a_ge[idim]) && (ge[idim] >= a_gs[idim]) ) /* [ [X[ [ */
                    {
#ifdef EXC_VERBOSE
		      fprintf(stream,"  case 5\n");
#endif            		    
                        bdrs[idim][0] = bdrs[idim][1] = gs[idim];
                        bdrs[idim][2] = bdrs[idim][3] = ge[idim] + 1;
                    }
                    else if ( (gs[idim] >  a_ge[idim]) && (ge[idim] >  a_ge[idim]) ) /* [ [ [X[ */
                    {
#ifdef EXC_VERBOSE
		      fprintf(stream,"  case 6\n");
#endif            		    
                        bdrs[idim][0] = bdrs[idim][1] = bdrs[idim][2] = gs[idim];
                        bdrs[idim][3] = ge[idim] + 1;
                    }
                    else return E_INTERNAL;
                }
                /* loop through receives */
                for ( ircv = 0; ircv < nrcv; ++ircv )
                {
                    if ( gen_i2pnt(ndim, ircv, rcvpnt) ) return E_INTERNAL;
                    
                    /* get corresponding part of the participating array */ 
                    for ( idim = 0; idim < ndim; ++idim )
                    {
                        gs[idim] = bdrs[idim][1+rcvpnt[idim]];
                        ge[idim] = bdrs[idim][2+rcvpnt[idim]] - 1;
#ifdef EXC_VERBOSE
			fprintf(stream,"ex_compute: gs, ge from bdrs, ircv=%d idim=%d gs=%d ge=%d\n",ircv,idim,gs[idim],ge[idim]);
#endif
			if ( ge[idim] < gs[idim] ) break; /* empty part */
                    }
                    if ( idim < ndim ) continue;/* nothing to do - empty part */
                    
                    /* so, non-empty part of the participating array */
                    if ( frcvempty[ircv] ) /* receive yet empty */
                    {
                        IASN(r_gs[ircv], gs);
                        IASN(r_ge[ircv], ge);
                        frcvempty[ircv] = 0;
#ifdef EXC_VERBOSE
			fprintf(stream,"ex_compute: ircv=%d, empty\n",ircv);
			for (int idim=0;idim<ndim; idim++) 
			  fprintf(stream,"  idim=%d gs=%d ge=%d r_gs=%d r_ge=%d\n",idim,gs[idim],ge[idim],r_gs[ircv][idim],r_ge[ircv][idim]);
#endif
                    }
                    else
                    {
                        for ( idim = 0; idim < ndim; ++idim )
                        {
			  r_gs[ircv][idim] = iwave_min(r_gs[ircv][idim], gs[idim]);
			  r_ge[ircv][idim] = iwave_max(r_ge[ircv][idim], ge[idim]);
#ifdef EXC_VERBOSE
			  fprintf(stream,"ex_compute: ircv=%d, non-empty\n",ircv);
			  for (int idim=0;idim<ndim; idim++) 
			    fprintf(stream,"  idim=%d gs=%d ge=%d r_gs=%d r_ge=%d\n",idim,gs[idim],ge[idim],r_gs[ircv][idim],r_ge[ircv][idim]);
#endif

                        }
                    }
                }/* receive areas */
            }/* if non-empty */
        }/* elements */
    }/* masks */
    
    /* compute allocted cube in case of non-empty computational area ---------*/
    if ( !fempty )
    {
        for ( ircv = 0; ircv < nrcv; ++ircv )
        {
            for ( idim = 0; idim < ndim; ++idim )
            {
                if ( r_gs[ircv][idim] < a_gs[idim] ) a_gs[idim] = r_gs[ircv][idim];
                if ( r_ge[ircv][idim] > a_ge[idim] ) a_ge[idim] = r_ge[ircv][idim];
            }
        }    
    }

    return 0;
}
