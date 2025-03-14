/*        Generated by TAPENADE     (INRIA, Ecuador team)
    Tapenade 3.10 (r5717) - 30 Jul 2015 16:03
*/
/*        Generated by TAPENADE     (INRIA, Ecuador team)
    Tapenade 3.10 (r5717) - 30 Jul 2015 16:03
*/
#include "cstd.h"

/*
  Differentiation of acd_2d_d in reverse (adjoint) mode:
   gradient     of useful results: **ucd **upd **uc **up
   with respect to varying inputs: **ucd **upd **csq **uc **up
   RW status of diff variables: **ucd:incr **upd:in-out **csq:out
                **uc:incr **up:in-out
   Plus diff mem management of: ucd:in *ucd:in lap:in upd:in *upd:in
                lapd:in csq:in *csq:in uc:in *uc:in up:in *up:in

*/
void acd_2d_d_b(float **uc, float **ucb, float **ucd, float **ucdb, float **up
        , float **upb, float **upd, float **updb, float **csq, float **csqb, 
        float **csqd, int *s, int *e, float **c, int k, int *lbc, int *rbc, 
        float *lap, float *lapb, float *lapd, float *lapdb) {

    int i0, i1, ioff;
    int s0 = s[0];
    int e0 = e[0];
    float c0 = c[0][0] + c[1][0];

    float tempb2;
    float tempb1;
    float tempb0;
    float tempb;

    /* boundary conditions - note that uc[-1][i]=0 etc. */
    if (rbc[0])
        for (i1 = e[1]; i1 > s[1]-1; --i1)
            for (ioff = k-2; ioff > -1; --ioff) {
                upb[i1][e[0] - ioff] = upb[i1][e[0] - ioff] - upb[i1][e[0] + ioff + 2];
                upb[i1][e[0] + ioff + 2] = 0.0;
                updb[i1][e[0] - ioff] = updb[i1][e[0] - ioff] - updb[i1][e[0] + ioff + 2];
                updb[i1][e[0] + ioff + 2] = 0.0;
            }

    if (lbc[0])
        for (i1 = e[1]; i1 > s[1]-1; --i1)
            for (ioff = k-2; ioff > -1; --ioff) {
                upb[i1][s[0] + ioff] = upb[i1][s[0] + ioff] - upb[i1][s[0] - ioff - 2];
                upb[i1][s[0] - ioff - 2] = 0.0;
                updb[i1][s[0] + ioff] = updb[i1][s[0] + ioff] - updb[i1][s[0] - ioff - 2];
                updb[i1][s[0] - ioff - 2] = 0.0;
            }
    if (rbc[1])
        for (ioff = k-2; ioff > -1; --ioff)
            for (i0 = e0; i0 > s0-1; --i0) {
                upb[e[1] - ioff][i0] = upb[e[1] - ioff][i0] - upb[e[1] + 2 + ioff][i0];
                upb[e[1] + 2 + ioff][i0] = 0.0;
                updb[e[1] - ioff][i0] = updb[e[1] - ioff][i0] - updb[e[1] + 2 + ioff][i0];
                updb[e[1] + 2 + ioff][i0] = 0.0;
            }
    if (lbc[1])
        for (ioff = k-2; ioff > -1; --ioff)
            for (i0 = e0; i0 > s0-1; --i0) {
                upb[s[1] + ioff][i0] = upb[s[1] + ioff][i0] - upb[s[1] - ioff - 2][i0];
                upb[s[1] - ioff - 2][i0] = 0.0;
                updb[s[1] + ioff][i0] = updb[s[1] + ioff][i0] - updb[s[1] - ioff - 2][i0];
                updb[s[1] - ioff - 2][i0] = 0.0;
            }

    for (i1 = e[1]; i1 > s[1]-1; --i1) {

        for (i0 = s0; i0 < e0+1; ++i0) {
            lapd[i0] = c0*ucd[i1][i0];
            lap[i0] = c0*uc[i1][i0];
        }

        // offset dependent part of stencil
        for (ioff = 1; ioff < k+1; ++ioff)
            for (i0 = s0; i0 < e0+1; ++i0) {
                lapd[i0] = lapd[i0] + c[0][ioff]*(ucd[i1][i0+ioff]+ucd[i1][i0-
                    ioff]) + c[1][ioff]*(ucd[i1+ioff][i0]+ucd[i1-ioff][i0]);
                lap[i0] = lap[i0] + (c[0][ioff]*(uc[i1][i0+ioff]+uc[i1][i0-
                    ioff]) + c[1][ioff]*(uc[i1+ioff][i0]+uc[i1-ioff][i0]));
            }

        for (i0 = e0; i0 > s0-1; --i0) {
            csqb[i1][i0] = csqb[i1][i0] + lap[i0]*upb[i1][i0];
            lapb[i0] = csq[i1][i0]*upb[i1][i0];
	    lapb[i0] += csqd[i1][i0]*updb[i1][i0];
            csqb[i1][i0] = csqb[i1][i0] + lapd[i0]*updb[i1][i0];
            lapdb[i0] = csq[i1][i0]*updb[i1][i0];
        }

        for (ioff = k; ioff > 0; --ioff)
            for (i0 = e0; i0 > s0-1; --i0) {
                tempb = c[0][ioff]*lapb[i0];
                tempb0 = c[1][ioff]*lapb[i0];
                ucb[i1][i0 + ioff] = ucb[i1][i0 + ioff] + tempb;
                ucb[i1][i0 - ioff] = ucb[i1][i0 - ioff] + tempb;
                ucb[i1 + ioff][i0] = ucb[i1 + ioff][i0] + tempb0;
                ucb[i1 - ioff][i0] = ucb[i1 - ioff][i0] + tempb0;
                tempb1 = c[0][ioff]*lapdb[i0];
                tempb2 = c[1][ioff]*lapdb[i0];
                ucdb[i1][i0 + ioff] = ucdb[i1][i0 + ioff] + tempb1;
                ucdb[i1][i0 - ioff] = ucdb[i1][i0 - ioff] + tempb1;
                ucdb[i1 + ioff][i0] = ucdb[i1 + ioff][i0] + tempb2;
                ucdb[i1 - ioff][i0] = ucdb[i1 - ioff][i0] + tempb2;
            }
        for (i0 = e0; i0 > s0-1; --i0) {
            ucb[i1][i0] = ucb[i1][i0] + c0*lapb[i0];
            ucdb[i1][i0] = ucdb[i1][i0] + c0*lapdb[i0];
            ucb[i1][i0] = ucb[i1][i0] + 2.0*upb[i1][i0];
            upb[i1][i0] = -upb[i1][i0];
            ucdb[i1][i0] = ucdb[i1][i0] + 2.0*updb[i1][i0];
            updb[i1][i0] = -updb[i1][i0];
        }
    }
}
