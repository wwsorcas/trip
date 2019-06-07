// MPS_spread.hh
// Author: Mario J. Bencomo
// last modified 02/20/17

/**
 * \file MPS_spread.hh
 * Spread operator for MPS's.
 */

#ifndef __MPS_SPREAD_HH_
#define __MPS_SPREAD_HH_

#include "MPS_includes.hh"
#include "MPS_Space.hh"

namespace TSOpt {

  /**
   * Low level kernel for stacking.
   *
   * @param [in] N_src - number of sources
   * @param [in] MPS_dim - dimension of MPS space
   * @param [in] x_file - input SEGY file
   * @param [in] y_file - output SEGY file
   */
  void MPS_spread_kern( int N_src,
			int MPS_dim,
			string x_file,
			string y_file );

  void MPS_stack_kern( int N_src,
		       int MPS_dim,
		       string x_file,
		       string y_file );

  //==============================================================================//
  // MPS_spread class
  //==============================================================================//
  /**
   * Stacking operator for MPS spaces.
   */
  class MPS_spread : public LinearOp<ireal> {
  
  private:

    MPS_Space const &dom; /**< domain MPS space */
    MPS_Space const &rng; /**< range MPS space */
    PARARRAY *pars; 

    void apply( const Vector<ireal> &x,
		Vector<ireal> &y ) const;

    void applyAdj( const Vector<ireal> &y,
		   Vector<ireal> &x ) const;

    void param_set( Vector<ireal> const &x, 
		    Vector<ireal> const &y ) const;

    Operator<ireal> * clone() const{ return new MPS_spread(*this);}
    
    void sanity_checks() const;
    
  public:

    /**
     * Copy constructor.
     */
    MPS_spread( MPS_spread const & op )
      : dom(op.dom), 
	rng(op.rng),
	pars(NULL){
      pars = ps_new();
      ps_copy(&pars,*op.pars);
    }

    /**
     * Main constructor.
     */
    MPS_spread( MPS_Space const &_dom, MPS_Space const &_rng )
      : dom(_dom), rng(_rng) {
      pars = ps_new();
    }

    /**
     * Destructor.
     */
    ~MPS_spread(){ ps_delete(&pars);}

    /**
     * Handy function retrieving MPS domain space.
     */
    const MPS_Space & getMPSDomain() const{ return dom;}
    const MPS_Space & getMPSRange() const{ return rng;}
    const Space<ireal> & getDomain() const {return dom;}
    const Space<ireal> & getRange() const {return rng;}
    ostream & write(ostream & str) const;
  };

}//end TSOpt

#endif
