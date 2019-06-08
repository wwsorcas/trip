/* C translation of fdcoeff.m, from R. J. Leveque, downloaded from 
   faculty.washington.edu/rjl/fdmbook/matlab/fdcoeff.m on 2015.08.22
   by WWS.

   to compile main program below to test, define LOCALTEST and something like
   
   g++ -I../../base/include -I../../../rvl/rvl/include fdcoeff.cc
*/
#include "utils.h"
#include "std_cpp_includes.hh"

//#define LOCALTEST

//function c = fdcoeffF(k,xbar,x)

//double * fdcoeff(int k, double xbar, std::vector<double> x) {
std::vector<double> fdcoeff(int k, double xbar, std::vector<double> x) {
  /*
    % Compute coefficients for finite difference approximation for the
    % derivative of order k at xbar based on grid values at points in x.
    %
    % This function returns a row vector c of dimension 1 by n, where n=length(x),
    % containing coefficients to approximate u^{(k)}(xbar), 
    % the k'th derivative of u evaluated at xbar,  based on n values
    % of u at x(1), x(2), ... x(n).  
    %
    % If U is a column vector containing u(x) at these n points, then 
    % c*U will give the approximation to u^{(k)}(xbar).
    %
    % Note for k=0 this can be used to evaluate the interpolating polynomial 
    % itself.
    %
    % Requires length(x) > k.  
    % Usually the elements x(i) are monotonically increasing
    % and x(1) <= xbar <= x(n), but neither condition is required.
    % The x values need not be equally spaced but must be distinct.  
    %
    % This program should give the same results as fdcoeffV.m, but for large
    % values of n is much more stable numerically.
    %
    % Based on the program "weights" in 
    %   B. Fornberg, "Calculation of weights in finite difference formulas",
    %   SIAM Review 40 (1998), pp. 685-691.
    %
    % Note: Forberg's algorithm can be used to simultaneously compute the
    % coefficients for derivatives of order 0, 1, ..., m where m <= n-1.
    % This gives a coefficient matrix C(1:n,1:m) whose k'th column gives
    % the coefficients for the k'th derivative.
    %
    % In this version we set m=k and only compute the coefficients for
    % derivatives of order up to order k, and then return only the k'th column
    % of the resulting C matrix (converted to a row vector).  
    % This routine is then compatible with fdcoeffV.   
    % It can be easily modified to return the whole array if desired.
    %
    % From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)
  */
  
  // n = length(x);
  // assumes that size of x can be represented in 4 bytes
  int n = x.size();
  if (k >= n) {
    //   error('*** length(x) must be larger than k')
    /*
    RVLException e;
    e<<"ERROR: fdcoeff\n";
    e<<"  order k="<<k<<" must be less than number of evaluation points = "<<n<<"\n";
    throw e;
    */
    exit(1);
  }

  int m = k;   
  //% change to m=n-1 if you want to compute coefficients for all
  //% possible derivatives.  Then modify to output all of C.
  double c1 = 1.0;
  double c4 = x[0] - xbar;
  // this appears to be a mistake - row index n is allowed in the matlab code
  // C = zeros(n-1,m+1)
  //  int nn = (n-1)*(m+1);
  int nn = n*(m+1);
  std::vector<double> CapC(nn);
  for (int i=0; i< nn; i++) CapC[i]=0.0;

  //C = zeros(n-1,m+1);
  //C(1,1) = 1;
  CapC[0]=1.0;
  //  for i=1:n-1
  for (int i=1; i<n; i++) {
    int i1 = i+1;
    int mn = min(i,m);
    double c2 = 1.0;
    double c5 = c4;
    c4 = x[i1-1] - xbar;
    for (int j=0;j<i;j++) {
      int j1 = j+1;
      double c3 = x[i1-1] - x[j1-1];
      c2 = c2*c3;
      //      cout <<"i="<<i<<" j="<<j<<" c2="<<c2<<" c3="<<c3<<" c4="<<c4<<endl;

      if (j==i-1) {
	//	for s=mn:-1:1
	for (int s=mn; s>0; s--) {
	  int s1 = s+1;
	  //        C(i1,s1) = c1*(s*C(i1-1,s1-1) - c5*C(i1-1,s1))/c2;
	  CapC[i1-1 + (s1-1)*n] = c1*(s*CapC[i1-2 + (s1-2)*n] - c5*CapC[i1-2 + (s1-1)*n])/c2;
        }
	CapC[i1-1] = -c1*c5*CapC[i1-2]/c2;
      }
      for (int s=mn; s>0; s--) {
	int s1 = s+1;
	//	C(j1,s1) = (c4*C(j1,s1) - s*C(j1,s1-1))/c3;
	CapC[j1-1 + (s1-1)*n] = (c4*CapC[j1-1+(s1-1)*n] - s*CapC[j1-1+(s1-2)*n])/c3;
      }
      CapC[j1-1] = c4*CapC[j1-1]/c3;
    }
    c1 = c2;
  }

  //c = C(:,end)';            % last column of c gives desired row vector
  std::vector<double> RetC(n);
  //double * RetC = (double *)usermalloc_(n*sizeof(double));
  for (int i=0; i<n; i++) {
    RetC[i]=CapC[i+n*m];
  }
  return RetC;
}

#ifdef LOCALTEST

// test main

int main(int argc, char ** argv) {
  cout<<"Levque's Code for Fornberg's Difference Coefficient Formulas\n";
  cout<<"  derivative order=1\n";
  int k=1;
  for (int ord=1; ord<5; ord++) {
    cout<<"    fd rder "<<2*ord<<":\n  ";
    std::vector<double> x(2*ord);
    for (int i=0;i<2*ord;i++) x[i]=double(i-ord)+0.5f;
    std::vector<double> c(2*ord);
    c = fdcoeff(k,0.0f,x);
    //double * c = fdcoeff(k,0.0f,x);
    cout<<"      ";
    for (int i=0;i<2*ord;i++) cout<<c[i]<<"  ";
    cout<<"\n";
    //    userfree_(c);
  }
  cout<<"  derivative order=2\n";
  k=2;
  for (int ord=1; ord<5; ord++) {
    cout<<"    fd rder "<<2*ord<<":\n  ";
    std::vector<double> x(2*ord+1);
    for (int i=0;i<2*ord+1;i++) x[i]=double(i-ord);
    std::vector<double> c(2*ord+1); 
    c = fdcoeff(k,0.0f,x);
    //    double * c = fdcoeff(k,0.0f,x);
    cout<<"      ";
    for (int i=0;i<2*ord+1;i++) cout<<c[i]<<"  ";
    cout<<"\n";
    //    userfree_(c);
  }
}
#endif
