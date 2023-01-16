#include "cstd.h"
//#include <omp.h>

void foo(float * u, float * v, float ** c, int s, int e, int hmax) {

#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
    int i; int ih; int sh; int eh;
    // i >=s and i+ih >=s - i >= s-ih i => max(s,s-ih)
    // i <=e and i+ih <=e - i <= e-ih i <= min(e,e-ih)
#ifdef _OPENMP
#pragma omp simd private(sh,eh)
    //#pragma omp for schedule(static, OMP_CHUNK_SIZE)
#endif    
    for (ih = -hmax; ih <= hmax; ++ih) {
      sh = (s>=s-ih ? s : s-ih);
      eh = (e<=e-ih ? e : e-ih);
      for (i = sh; i <= eh; ++i) {
#pragma ivdep
	u[i] += c[ih][i]*v[i+ih];
      }
    }

  }
}


