/*        Generated by TAPENADE     (INRIA, Ecuador team)
    Tapenade 3.14 (r7114) - 22 Oct 2018 09:36
*/
/*        Generated by TAPENADE     (INRIA, Ecuador team)
    Tapenade 3.14 (r7114) - 22 Oct 2018 09:36
*/
#include "cstd.h"

/*
  Differentiation of sggrad in forward (tangent) mode:
   variations   of useful results: *y
   with respect to varying inputs: *x
   RW status of diff variables: *x:in *y:out
   Plus diff mem management of: x:in y:in

 staggered grid derivative in dimension idx
 */
void sggrad_d(int *gsx, int *gex, int *gsy, int *gey, int *csy, int *cey, int 
        m, int idx, int dim, float *x, float *xd, float *a, float *y, float *
        yd, int *sggrad) {
    // only one error return, if idx out-of-range
    int ierr = 0;
    // compute offsets, lengths
    int *oxy;
    oxy = (int *)malloc(dim*sizeof(int));
    int *nx;
    nx = (int *)malloc(dim*sizeof(int));
    int *ny;
    ny = (int *)malloc(dim*sizeof(int));
    for (int i = 0; i < dim; ++i) {
        oxy[i] = gsy[i] - gsx[i];
        nx[i] = gex[i] - gsx[i] + 1;
        ny[i] = gey[i] - gsy[i] + 1;
    }
    if (idx == 0) {
        *yd = 0.0;
        for (int i1 = csy[1]-gsy[1]; i1 < cey[1]-gsy[1]+1; ++i1) {
            //#pragma ivdep
            //#pragma GCC ivdep
            for (int i0 = csy[0]-gsy[0]; i0 < cey[0]-gsy[0]+1; ++i0) {
                yd[i0 + i1*ny[0]] = a[0]*(xd[i0+oxy[0]+1+(i1+oxy[1])*nx[0]]-xd
                    [i0+oxy[0]+(i1+oxy[1])*nx[0]]);
                y[i0 + i1*ny[0]] = a[0]*(x[i0+oxy[0]+1+(i1+oxy[1])*nx[0]]-x[i0
                    +oxy[0]+(i1+oxy[1])*nx[0]]);
            }
            for (int j = 1; j < m; ++j)
                //#pragma ivdep
                //#pragma GCC ivdep
                for (int i0 = csy[0]-gsy[0]; i0 < cey[0]-gsy[0]+1; ++i0) {
                    yd[i0 + i1*ny[0]] = yd[i0 + i1*ny[0]] + a[j]*(xd[i0+oxy[0]
                        +1+j+(i1+oxy[1])*nx[0]]-xd[i0+oxy[0]-j+(i1+oxy[1])*nx[
                        0]]);
                    y[i0 + i1*ny[0]] += a[j]*(x[i0+oxy[0]+1+j+(i1+oxy[1])*nx[0
                    ]]-x[i0+oxy[0]-j+(i1+oxy[1])*nx[0]]);
                }
        }
    } else if (idx == 1) {
        *yd = 0.0;
        for (int i1 = csy[1]-gsy[1]; i1 < cey[1]-gsy[1]+1; ++i1) {
            //#pragma ivdep
            //#pragma GCC ivdep
            for (int i0 = csy[0]-gsy[0]; i0 < cey[0]-gsy[0]+1; ++i0) {
                yd[i0 + i1*ny[0]] = a[0]*(xd[i0+oxy[0]+(i1+oxy[1]+1)*nx[0]]-xd
                    [i0+oxy[0]+(i1+oxy[1])*nx[0]]);
                y[i0 + i1*ny[0]] = a[0]*(x[i0+oxy[0]+(i1+oxy[1]+1)*nx[0]]-x[i0
                    +oxy[0]+(i1+oxy[1])*nx[0]]);
            }
            for (int j = 1; j < m; ++j)
                //#pragma ivdep
                //#pragma GCC ivdep
                for (int i0 = csy[0]-gsy[0]; i0 < cey[0]-gsy[0]+1; ++i0) {
                    yd[i0 + i1*ny[0]] = yd[i0 + i1*ny[0]] + a[j]*(xd[i0+oxy[0]
                        +(i1+oxy[1]+1+j)*nx[0]]-xd[i0+oxy[0]+(i1+oxy[1]-j)*nx[
                        0]]);
                    y[i0 + i1*ny[0]] += a[j]*(x[i0+oxy[0]+(i1+oxy[1]+1+j)*nx[0
                    ]]-x[i0+oxy[0]+(i1+oxy[1]-j)*nx[0]]);
                }
        }
    } else {
        ierr = 1;
        *yd = 0.0;
    }
    free(oxy);
    free(nx);
    free(ny);
}
