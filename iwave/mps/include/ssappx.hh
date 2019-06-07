// ssappx.hh
// Author: Mario J. Bencomo
// last modified 01/26/16

/**
 * \file pntsrc_appx.hh
 * For computing weights associated with singular source appx.
 */

#ifndef __SSAPPX_HH_
#define __SSAPPX_HH_

#include "MPS_includes.hh"

namespace TSOpt{
  
  /**
   * Truncated factorial: m!/(m-s)! = m*(m-1)*...(m-s+1)
   */
  int factorial(int m);
  
  void midx_update( vector< vector<int> > &MIDX,
		    vector<int> const &midx, 
		    int N, int i, int ell );

  vector< vector<int> > midx_set( int q, int s, int ell );

  float poly_ell( float x, int ell, float h, int s, int q);

  float ssappx( float x, float h, int s, int q);

}
#endif
