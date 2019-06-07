/ Compatibility file for C99 and C++ complex.  This header
// can be included by either C99 or ANSI C++ programs to
// allow complex arithmetic to be written in a common subset.
// Note that overloads for both the real and complex math
// functions are available after this header has been
// included.

#ifndef MATH_COMPLEX_H_INCLUDED
#define MATH_COMPLEX_H_INCLUDED

#ifdef __cplusplus

#include <cmath>
#include <complex>

using namespace std;

typedef complex<float> float_complex;
typedef complex<double> double_complex;
typedef complex<long double> long_double_complex;

#else

// Note that <tgmath.h> includes <math.h> and <complex.h>
#include <tgmath.h>

typedef float complex float_complex;
typedef double complex double_complex;
typedef long double complex long_double_complex;

#define float_complex(r,i) ((float)(r) + ((float)(i))*I)
#define double_complex(r,i) ((double)(r) + ((double)(i))*I)
#define long_double_complex(r,i) ((long double)(r) + ((long double)(i))*I)

#define real(x) creal(x)
#define imag(x) cimag(x)
#define abs(x) fabs(x)
#define arg(x) carg(x)

#endif  // #ifdef __cplusplus

#endif  // #ifndef MATH_COMPLEX_H_INCLUDED
