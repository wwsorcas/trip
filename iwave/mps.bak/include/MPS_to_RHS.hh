// MPS_to_RHS.hh
// Author: Mario J. Bencomo
// last modified: 12/12/16

/**
 * \file MPS_to_RHS.hh
 * Implementation of MPS-to-RHS map.
 *
 * Defines MPS_to_RHS, derived class of LinearOp.
 * Domain and range spaces are initialized when operator is constructed.
 * The following MPS_KEYS will be need to be initialized with their values 
 * supplied by user (see main constructor for more info):
 *   - MPS_KEYS.appx_ord
 *   - MPS_KEYS.grid_file
 *   - MPS_KEYS.RHS_files **option, val='empty' for make_RHS=true**
 * Other MPS_KEYS do not need to be initialized.
 */

#ifndef __MPS_TO_RHS_HH_
#define __MPS_TO_RHS_HH_

#include "MPS_includes.hh"
#include "MPS_Space.hh"
#include "ssappx.hh"

namespace TSOpt {

  /**
   * Struct for source stencil points.
   */
  typedef struct s_MPS_point{
    RPNT glo;   /**< global spatial location */
    RPNT loc;   /**< local  spatial location */
    IPNT ext;   /**< -1= stencil point is out of bounds, must apply some extension 
		      0= boundary point,
		      1= interior point */

  } MPS_point;

  /**
   * Struct for source stencil related to singular source approx.
   */
  typedef struct s_MPS_stencil{
    vector<MPS_point> points; /**< stencil points*/
    IPNT n;                   /**< number of stencil points in each dim*/
  } MPS_stencil;

  /**
   * This function writes out a SEGY file associated with an RHS sources from scratch, 
   * by duplicating a reference trace.
   *
   * @param [in] stencils - source stencil data struct
   * @param [in] s_pos - source positions
   * @param [in] ref_file - reference SEGY filename
   * @param [in] RHS_files - filenames of RHS SEGY files to be written
   */
  void write_RHS_SEGY( vector<MPS_stencil> stencils,
		       vector<Rtuple> s_pos,
		       string ref_file,
		       vector<string> RHS_files );


  /**
   * Low level kernel for MPS_to_RHS forward action.
   *
   * @param [in] a_ord - approximation order
   * @param [in] d_ord - derivative order
   * @param [in] g - spatial grid
   * @param [in] s_pos - source position
   * @param [in] stencils - source stencil
   * @param [in] MPS_basis - MPS basis
   * @param [in] MPS_file - input MPS SEGY 
   * @param [in] RHS_files - output RHS SEGY files
   */
  void MPS_to_RHS_kern( int a_ord,
			Ituple d_ord,
			grid g,
			vector<Rtuple> s_pos,
			vector<MPS_stencil> stencils,
			vector<MPS_base> MPS_basis,
			string MPS_file,
			vector<string> RHS_files );

  /**
   * Low level kernel for MPS_to_RHS adjoint action.
   *
   * @param [in] a_ord - approximation order
   * @param [in] d_ord - derivative order
   * @param [in] g - spatial grid
   * @param [in] s_pos - source positions
   * @param [in] stencils - source stencil
   * @param [in] MPS_basis - MPS basis
   * @param [in] MPS_file - input MPS SEGY 
   * @param [in] RHS_files - output RHS SEGY files
   */
  void RHS_to_MPS_kern( int a_ord,
			Ituple d_ord,
			grid g,
			vector<Rtuple> s_pos,
			vector<MPS_stencil> stencils,
			vector<MPS_base> MPS_basis,
			string MPS_file,
			vector<string> RHS_files );

  // The MPS_to_RHS class, a purely virtual class since it depends on MPS_Space,
  // is an LinearOp in ways similar to IWaveOp. It provides a link between the MPS 
  // representation and an IWAVE RHS source input.
  // An MPS_to_RHS subclass must be implemented for a given subclass of MPS_Space.
  // Examples are provided in MPS_Can.cc and MPS_Ex.cc.

  // DATA MEMBERS:
  //                        int a_ord = order of MPS approximation
  //                           grid g = spatial grid info
  // vector<MPS_stencil> MPS_stencils = MPS stencils for src approx
  //                  PARARRAY * pars = parameter table
  //                     MPS_KEYS mks = MPS keys
  //                    bool make_RHS = bool flag for making RHS space from scratch 
  //    shared_ptr<MPS_Space> dom_ptr = MPS space domain
  //   shared_ptr<IWaveSpace> rng_ptr = RHS space range

  //=============================================================================//
  // MPS_to_RHS
  //=============================================================================//
  /**
   * Linear map from MPS to RHS sources.
   */
  class MPS_to_RHS : public LinearOp<ireal> {

  private:

    MPS_Space const &dom; /**< domain MPS space */
    shared_ptr<IWaveSpace> rng_ptr; /**< range space */

    IWaveInfo ic;    /**< IWaveInfo */
    int a_ord;       /**< approximation order */
    grid g;          /**< spatial grid */
    PARARRAY *pars;  /**< key-val pairs */
    MPS_KEYS mks;    /**< MPS keys */
    bool make_RHS;   /**< boolean flag for making RHS SEGY files from scratch */
    vector<MPS_stencil> stencils; /**< source stencils */

    //methods for LinearOp
    void apply( const Vector<ireal> &x,
		Vector<ireal> &y ) const;
    void applyAdj( const Vector<ireal> &y,
		   Vector<ireal> &x ) const;

    //other methods
    void param_set( PARARRAY &locpars, 
		    Vector<ireal> const &x,
		    Vector<ireal> const &y ) const;

    /**
     * Handy function for setting range space.
     */
    void set_rng_ptr();

    /**
     * Handy function for setting source stencil
     */
    void set_stencils();

    /**
     * Handy function for sanity checks
     */
    void sanity_checks() const;

    Operator<ireal> * clone() const{ return new MPS_to_RHS(*this);}

  public:

    /**
     * Empty constructor.
     */
    /*
    MPS_to_RHS()
      : dom(), 
	rng_ptr(), 
	pars(NULL), 
	s_idx(0), 
	s_res(0),
	stencils(0) {}
    */

    /**
     * Copy constructor.
     */
    MPS_to_RHS( MPS_to_RHS const &MtoR );   

    /**
     * Main constructor.
     *
     * MPS_to_RHS main constructor requires the following MPS_KEYS be 
     * defined:
     *   - MPS_KEYS.appx_ord - singular source approximation order 
     *   - MPS_KEYS.grid_file - filename to RSF file containing spatial grid info
     *   - MPS_KEYS.RHS_files - values can be set to 'empty' or a given filename
     *
     * @param [in] _mks - keys related to initialization of MPS spaces, see \ref MPS_KEYS
     * @param [in] pars - input key-val pairs
     * @param [in] mps - MPS domain space 
     * @param [in] _make_RHS - boolean flag for making RSH files related to range from scratch
     */
    MPS_to_RHS( MPS_KEYS mks, 
		PARARRAY const &_pars, 
		MPS_Space const & mps,
		bool _make_RHS=true );
    
    /**
     * Desctructor.
     */
    ~MPS_to_RHS();
    
    /**
     * Handy function for retrieving MPS space domain.
     */
    const MPS_Space & get_MPS_Domain() const{ return dom;}

    /**
     * Handy function for retrieving IWaveSpace range. 
     */
    const IWaveSpace & get_IWave_Range() const {return *(rng_ptr);}

    const Space<ireal> & getDomain() const {return get_MPS_Domain();}
    const Space<ireal> & getRange() const {return get_IWave_Range();}

    ostream & write(ostream & str) const{ 
      str<<"MPS_to_RHS operator\n";
      str<<"domain:\n";
      dom.write(str);
      str<<"range:\n";
      if( !(rng_ptr) )
	str<<"null pointer\n";
      else
	rng_ptr->write(str);
      return str;
    }

    MPS_KEYS get_MPS_keys() const{return mks;}
    PARARRAY & get_pars(){ return *pars;}
    PARARRAY const & get_pars() const{ return *pars;}    
    grid get_grid() const{ return g;}
    int get_a_ord() const{ return a_ord;}
    vector<MPS_stencil> get_stencils() const{ return stencils;}
    string get_dom_filename() const;
    vector<string> get_rng_filenames() const;
    void print_info(ostream &str) const;
    void print_info() const{ print_info(cerr);}
  };



}//end TSOpt

#endif
