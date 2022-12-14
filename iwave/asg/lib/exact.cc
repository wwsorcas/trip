#include "utils.h"
#include "std_cpp_includes.hh"

float lininterp(int nt, float dt, float t, float * f) {
  float s = t/dt;
  int it = int(s);
  if (it < 0 || it > nt-2) return 0.0f;
  else return (s-float(it))*f[it+1]+(float(it)+1.0f-s)*f[it];
}

/* Poisson's formula for point trace of 2D point source pressure
   field, after T. Vdovina, TRIP Annual Report presentation 07

dp/dt = bulk[-div v + f(t) delta(x-x_s)]
dv/dt = -buoy grad p

output is p(x,t) where |x-x_s| = r

limitation: f, p defined on same interval

*/

void exact(int nt, float dt, float r,
	   float bulk, float buoy,
	   float * f,
	   float * p) {
  // compute df
  float * df = (float *)(usermalloc_(nt*sizeof(float)));
  df[0]=(f[1]-f[0])/dt;
  df[nt-1]=(f[nt-1]-f[nt-2])/dt;
  for (int i=1;i<nt-1;i++)
    df[i]=(f[i+1]-f[i-1])/(2.0f*dt);

  // sound speed
  float c = sqrt(bulk*buoy);

  // traveltime
  float xoc = r/c;

  // what again
  float pi = 4.0f*atan(1.0);

  for (int i=0;i<nt;i++) {
    // initialize
    p[i]=0.0f;
    // time
    float t = i*dt;
    // only nonzero after arrival time
    if (t >= xoc) {
      // limit of tau integration
      float taulim = sqrt(t-xoc);
      // greatest multiple of dt below
      int lim = int(taulim/dt);
      // buffer for integral
      float a = 0.0f;
      for (int j=0; j<=lim; j++) {
	// integration variable
	float tau = j*dt;
	// poisson denom
	float denom=sqrt(tau*tau+2.0f*xoc);
	// argument of df/dt
	float targ = t-xoc-tau*tau;
	// accumulate
	a += lininterp(nt,dt,targ,df)/denom;
      }
      // scale by dt and other factors
      p[i]=dt*a/(pi*c*c);
    }
  }
  // cleanup
  userfree_(df);
}
