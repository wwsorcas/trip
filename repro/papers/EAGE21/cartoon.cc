#include <sys/types.h>
#include <sys/uio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <time.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <exception>
#include <stdexcept>
#include <typeinfo>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <valarray>
#include <complex>
#include <memory>

using namespace std;

#define NV 101
#define VE 2.0
#define V0 0.5
#define DV 0.01
#define C0 1.0
#define C1 0.1
#define C2 4.0
#define C3 0.05
#define C4 0.2
#define C5 7.0
#define C6 50.0
#define C7 500.0

int main(int argc, char ** argv) {

  float * v   = new float[NV];
  float * j   = new float[NV];
  float * h   = new float[NV];
  float * t   = new float[NV];
  float * f   = new float[NV];

  for (int i=0; i< NV; i++) {
    v[i] = i*DV - V0;
    h[i] = C0*v[i]*v[i];
    j[i] = h[i] - C3*sin(C2*v[i])*sin(C2*v[i])*sin(C2*v[i])*sin(C2*v[i])
      + C1*sin(C5*v[i])*sin(C5*v[i])*sin(C5*v[i])*C3*sin(C5*v[i])*sin(C5*v[i])*sin(C5*v[i]);
    t[i]= h[i] + C4*sin(C2*v[i])*sin(C2*v[i])*sin(C2*v[i])*sin(C2*v[i])
      - C1*sin(C2*v[i])*sin(C2*v[i])*sin(C2*v[i])*C1*sin(C2*v[i])*sin(C2*v[i])*sin(C2*v[i]);
    f[i]= 0.5*(C7*h[i]+C0*sin(C6*v[i])*sin(C6*v[i]))/(1.0+C7*h[i]);
  }

  FILE * dat = fopen("cartoon.rsf@","w");
  fwrite(h, sizeof(float), NV, dat);
  fwrite(j, sizeof(float), NV, dat);
  fwrite(t, sizeof(float), NV, dat);
  fwrite(f, sizeof(float), NV, dat);    
  fclose(dat);

  ofstream hdr("cartoon.rsf");
  hdr<<"n1="<<NV<<" o1="<<VE-V0<<" d1="<<DV<<endl;
  hdr<<"n2=4 o2=0 d2=1"<<endl;
  hdr<<"data_format=\"native_float\""<<endl;
  hdr<<"esize=4"<<endl;
  hdr<<"in=\"cartoon.rsf@\""<<endl;
  hdr.close();

  delete [] v;
  delete [] h;
  delete [] j;
  delete [] t;
  delete [] f;
}  
