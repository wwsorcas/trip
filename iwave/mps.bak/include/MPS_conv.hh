// MPS_conv.hh
// Author: Mario J. Bencomo
// last modified 11/10/16

/**
 * \file MPS_conv.hh
 * Convolution operator for MPS's.
 */

#ifndef __MPS_CONV_HH_
#define __MPS_CONV_HH_

#include "MPS_includes.hh"
#include "MPS_Space.hh"

namespace TSOpt {

  /**
   * Low level kernel for convolution.
   *
   * @param [in] s_pos - source positions
   * @param [in] MPS_dim - dimension of MPS space
   * @param [in] MPS_idx - MPS index
   * @param [in] x_file - input SEGY file
   * @param [in] y_files - vector of output SEGY files
   * @param [in] k_files - vector of kernel SEGY files
   */
  void MPS_conv_kern( vector<Rtuple> s_pos,
		      int MPS_dim,
		      int MPS_idx,
		      string x_file,
		      vector<string> y_files,
		      vector<string> k_files );

  /**
   * Low level kernel for cross-correlation.
   *
   * @param [in] s_pos - source positions
   * @param [in] MPS_dim - dimension of MPS space
   * @param [in] MPS_idx - MPS index
   * @param [in] x_file - output SEGY file
   * @param [in] y_files - vector of input SEGY files
   * @param [in] k_files - vector of kernel SEGY files
   */
  void MPS_corr_kern( vector<Rtuple> s_pos,
		      int MPS_dim,
		      int MPS_idx,
		      string x_file,
		      vector<string> y_files,
		      vector<string> k_files );

  // MPI versions
#ifdef IWAVE_USE_MPI
  void MPS_conv_kern_MPI( vector<Rtuple> s_pos,
			  int MPS_dim,
			  int MPS_idx,
			  string x_file,
			  vector<string> y_files,
			  vector<string> k_files );
  
  void MPS_corr_kern_MPI( vector<Rtuple> s_pos,
			  int MPS_dim,
			  int MPS_idx,
			  string x_file,
			  vector<string> y_files,
			  vector<string> k_files );
#endif
  

  //==============================================================================//
  // MPS_conv class
  //==============================================================================//
  /**
   * Convolution operator for MPS spaces.
   */
  class MPS_conv : public LinearOp<ireal> {
  
  private:

    MPS_Space const &dom;                /**< domain MPS space */
    IWaveSpace const &rng;               /**< range IWaveSpace */
    shared_ptr< Vector<ireal> > ker_ptr; /**< Green's function kernels */

    PARARRAY *pars;          /**< parameter array */
    vector<string> ker_keys; /**< keys related to kernels */
    int MPS_idx;             /**< MPS index */

    void apply( const Vector<ireal> &x,
		Vector<ireal> &y ) const;

    void applyAdj( const Vector<ireal> &y,
		   Vector<ireal> &x ) const;

    /**
     * Handy function call for sanity checks.
     */
    void sanity_checks();

    /**
     * Handy function call for setting parameters.
     */
    void param_set( Vector<ireal> const &x, 
		    Vector<ireal> const &y ) const;

    Operator<ireal> * clone() const{ return new MPS_conv(*this);}

  public:

    /**
     * Copy constructor.
     */
    MPS_conv( MPS_conv const &conv );

    /**
     * Constructor with specified Green's function kernels.
     */
    MPS_conv( PARARRAY const &_pars,
	      Vector<ireal> &k, 
	      MPS_Space const &d,
	      IWaveSpace const &r);
    
    /**
     * Constructor without specified Green's function kernels.
     */
    MPS_conv( PARARRAY const &_pars,
	      MPS_Space const &d,
	      IWaveSpace const &r );

    /**
     * Destructor.
     */
    ~MPS_conv();

    static shared_ptr<MPS_conv> newPtr(MPS_conv const &conv_op){
      shared_ptr<MPS_conv> p(new MPS_conv(conv_op));
      return p;
    }

    /**
     * Function for setting MPS index in convolution operator.
     */
    void set_MPS_idx(int idx){MPS_idx=idx;}

    /**
     * Function for setting Green's function kernels.
     */
    void set_kernel(Vector<ireal> const &k);
    
    /**
     * Function for changing Green's function kernels.
     */
    void change_kernel(Vector<ireal> const &k);
    
    /**
     * Handy function retrieving MPS domain space.
     */
    const MPS_Space & get_MPS_Domain() const{ return dom;}

    /**
     * Handy function retrieving IWaveSpace range space.
     */
    const IWaveSpace & get_IWave_Range() const{ return rng;}
    //    const SubIWaveSpace & get_SubIWave_Green() const;

    const Space<ireal> & getDomain() const {return dom;}
    const Space<ireal> & getRange() const {return rng;}
    ostream & write(ostream & str) const;
    void print_info() const{}
  };

}//end TSOpt

#endif
