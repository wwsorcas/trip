#include "utils.h"
#include "std_cpp_includes.hh"

extern std::vector<double> fdcoeff(int,
				   double,
				   std::vector<double>);

namespace RVL {
  
// staggered grid coefficients
//#define OLD
#undef OLD

#ifdef OLD

float * sgcoeffs(int k) {
  int n,i;
  float * c = (float *)usermalloc_(k*sizeof(float));
  for (n=1;n<=k;n++) {
    c[n-1]=1.0f;
    for (i=1;i<=k;i++) {
      if (i != n) {
	c[n-1] *= ((2.0f*i-1.0f)*(2.0f*i-1.0f))
	  /(((2.0f*n-1.0f)*(2.0f*n-1.0f)) - ((2.0f*i-1.0f)*(2.0f*i-1.0f)));
      }
    }
    c[n-1]=fabs(c[n-1])/(2.0f*n-1.0f);
    if (((n-(n/2)*2) == 0)) c[n-1]=-c[n-1];
  }
  return c;
}

#else

// version using Leveque's general computation of FD coeffs
// use target point xbar=0, nodes
// -k+1/2,...,k-1/2 for normalized coefficients
// note that sgcoeffs only returns coeffs for positive offsets

float * sgcoeffs(int k) {
  float * c = (float *)usermalloc_(k*sizeof(float));
  std::vector<double> d(2*k);
  std::vector<double> x(2*k);
  double xbar=0.0;
  x[0]=-double(k)+0.5;
  for (int i=1;i<2*k;i++) x[i]=x[i-1]+1.0;
  d = fdcoeff(1,xbar,x);
  for (int i=0;i<k;i++) c[i]=float(d[k+i]);
  return c;
}

#endif

}
