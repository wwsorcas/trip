// MPS_frac_cal.hh
// Author: Mario J. Bencomo
// last modified 04/10/18

/**
 * \file MPS_frac_cal.hh
 * Fractional calculus operators for MPS spaces.
 */

#ifndef __MPS_FRAC_CAL_HH_
#define __MPS_FRAC_CAL_HH_

#define MY_MAX(x,y) ((x)>(y))?(x):(y)
#define MY_MIN(x,y) -(MY_MAX(-(x),-(y)))
#define MINUS_TO_POW(p) ((p)%2)?(-1):(1)

#include "MPS_includes.hh"
#include "MPS_Space.hh"
#include "ssappx.hh"

namespace TSOpt {

  /**
   * Computing binomial coeff via recursion (for small n,k<10).
   */
  int nchoosek(int n, int k);


  /**
   * Low level fractional derivative kernel.
   *
   * @param [in] q - derivative order
   * @param [in] scal - float scaling factor
   * @param [in] N - size of data
   * @param [in] x - input data
   * @param [out] y - output data
   * @param [in] dt - time step size
   * @param [in] adj - boolean adjoint flag
   */
  void frac_deriv( float q, 
		   float scal, 
		   int N, 
		   const float *x, 
		   float *y, 
		   float dt, 
		   bool adj );

  /**
   * Low level fractional integral kernel.
   *
   * @param [in] q - integration order
   * @param [in] scal - float scaling factor
   * @param [in] N - size of data
   * @param [in] x - input data
   * @param [out] y - output data
   * @param [in] dt - time step size
   * @param [in] adj - boolean adjoint flag
   */
  void frac_integ( float q, 
		   float scal, 
		   int N, 
		   const float *x, 
		   float *y, 
		   float dt, 
		   bool adj );

  /**
   * Low level kernels for MPS frac derivative/integrals and adjoints.
   *
   * @param [in] order - fractional order of derivative/integral
   * @param [in] scalars - scaling factors
   * @param [in] N_src - number of sources
   * @param [in] in_file - filename of input SEGY
   * @param [in] out_file - filename of output SEGY
   * @param [in] adj - boolean adjoint flag
   */
  void MPS_frac_cal_kern( vector<float> order,
			  vector<float> scalars,
			  int N_src,
			  string in_file,
			  string out_file,
			  bool adj,
			  bool inv );

#ifdef IWAVE_USE_MPI
  /** 
   * MPI version of MPS_frac_cal_kern.
   */
  void MPS_frac_cal_kern_MPI( vector<float> order,
			      vector<float> scalars,
			      int N_src,
			      string in_file,
			      string out_file,
			      bool adj,
			      bool inv );
#endif


  //==============================================================================//
  // MPS_frac_cal class
  //==============================================================================//
  /**
   * Fractional derivative/integral operators for MPS's.
   */
  class MPS_frac_cal : public RVL::LinearOpWithInverse<ireal> {
  
  private:

    MPS_Space const &mps_sp; /**< Domain and range spaces.*/
    vector<float> orders;    /**< fractional derivative/integral orders */
    PARARRAY *pars;          /**< array of key-val pairs */
    vector<float> scalars;   /**< scaling factors */
    string in_key;           /**< input key for op*/
    string out_key;          /**< output key for op*/

    //private methods
    void apply( const Vector<ireal> &x,
		Vector<ireal> &y ) const;

    void applyAdj( const Vector<ireal> &y,
		   Vector<ireal> &x ) const;
    
    void applyInv( const Vector<ireal> &y,
		   Vector<ireal> &x ) const;
    
    void applyInvAdj( const Vector<ireal> &x,
		      Vector<ireal> &y ) const;
    
    void sanity_checks();
    void param_set( Vector<ireal> const &x, 
		    Vector<ireal> const &y ) const;

    Operator<ireal> * clone() const{ return new MPS_frac_cal(*this);}

  public:

    /**
     * Copy constructor.
     */
    MPS_frac_cal( MPS_frac_cal const &op );

    /**
     * Main constructor.
     *
     * @param [in] _pars - array of key-val pairs
     * @param [in] _sp - domain and range space of op
     * @param [in] c - wave speed for scaling op
     * @param [in] order_0 - inital order of fractional op
     * @param [in] order_d - increment of order of fractional op
     * @param [in] inv - boolean flag for inverse op
     * @param [in] _in_key - key for input SEGY file
     * @param [in] _out_key - key for output SEGY file
     */
    MPS_frac_cal( PARARRAY const &_pars,
		  MPS_Space const &_sp,
		  float c,
		  float order_0,
		  float order_d,
		  string _in_key="MPS_frac_cal_in",
		  string _out_key="MPS_frac_cal_out");
    
    /**
     * Destructor.
     */
    ~MPS_frac_cal();

    /*
    static std::shared_ptr<MPS_conv> newPtr(MPS_conv const &conv_op){
      std::shared_ptr<MPS_conv> p(new MPS_conv(conv_op));
      return p;
    }
    */

    /**
     * Handy function for retrieving MPS domain space.
     */
    const MPS_Space & get_MPS_Domain() const{ return mps_sp;}

    /**
     * Handy function for retrieving MPS domain space.
     */
    const MPS_Space & get_MPS_Range() const{ return get_MPS_Domain();}    

    const Space<ireal> & getDomain() const {return get_MPS_Domain();}
    const Space<ireal> & getRange() const {return get_MPS_Range();}
    ostream & write(ostream & str) const;
    void print_info() const{}
  };








}//end TSOpt
#endif
