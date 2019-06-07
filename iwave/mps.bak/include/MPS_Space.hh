// MPS_Space.hh
// Author: Mario J. Bencomo
// last modified 01/27/16

/**
 * \file MPS_Space.hh
 * MPS space base class.
 *
 * Defines MPS space base class and related structs and functions.
 *
 * The MPS_Space constructor will need the following MPS_KEYS initialized
 * with their values supplied by user (see main constructor for more info):
 *   - MPS_KEYS.MPS_file
 *   - MPS_KEYS.MPS_grid
 *   - MPS_KEYS.MPS_ref   **needed for make=false**
 *   - MPS_KEYS.hdr_files **needed for make=false**
 * Other MPS_KEYS do not need to be initialized.
 */

#ifndef __MPS_SPACE_HH_
#define __MPS_SPACE_HH_

#include "MPS_includes.hh"

#ifndef MAX_MPS
 #define MAX_MPS 1+2+9+27
#endif

#ifndef MAX_D_ORD
 #define MAX_D_ORD 4
#endif

namespace TSOpt {

  template< typename T >
  void print_vec(vector<T> v){
    for(int i=0; i<v.size(); i++)
      cerr << v[i] << "\n";
  }

  /**
   * Useful wrapper struct for IPNT, used in templating C++ vectors.
   */
  typedef struct s_Ituple{
    IPNT coor;
  } Ituple;

  /**
   * Useful wrapper struct for RPNT, used in templating C++ vectors.
   */
  typedef struct s_Rtuple{
    RPNT coor;
  } Rtuple;


  /**
   * Struct containing all info pertinent to an MPS basis.
   */
  typedef struct s_MPS_base{
    int  NT;                  /**< number of terms */              
    vector<Ituple> derv;      /**< derivative info in multi-index form*/
    vector<int> RHS_flags;    /**< flags related to components fo RHS sources */
  } MPS_base;


  /**
   * Struct containing all keys related to MPS_Space.
   */     
  typedef struct s_MPS_KEYS{
    string appx_ord;          /**< singular source approximation order */
    string grid_file;         /**< rsf filename containing spatial grid info */
    string MPS_file;          /**< segy filename for initializing MPS_Space */
    string MPS_ref;           /**< reference segy filename for generating MPS file from scratch */
    vector<string> RHS_files; /**< vector of segy filenames for initializing RHS IWaveSpace's */
    vector<string> hdr_files; /**< vector of segy filenames associated with data, containing src positions */
    vector<string> G_files;   /**< vector of segy filenames for initializing Green data */
    string delta_file;        /**< segy filename associated with input MPS for generating Green data */
  } MPS_KEYS;


  /** \defgroup add_to_pars
   * Clean function call to add key-val pairs to a PARARRAY.
   * @{
   */
  void add_to_pars(PARARRAY const &pars, vector<string> keys, vector<string> vals);
  void add_to_pars(PARARRAY const &pars, string key, string val);
  void add_to_pars(PARARRAY const &pars, string key, int val);
  /** @}*/

  /** \defgroup get_ntr 
   * Convenient functional call to get total number of traces from SEGY file.
   * @{
   */
  int get_ntr(string sufile);
  int get_ntr_serial(string sufile);
  /** @}*/

  /**
   * This function extracts source/receiver positions of given sufile.
   *
   * @param [in] sufile - SEGY filename to extract source/receiver positions from.
   * @param [in] rcv - boolean flag, =true if receiver positions are to be extracted,
   * else extract source positions.
   * @return - vector of source positions
   */
  vector<Rtuple> ex_pos( string sufile, bool rcv );
  vector<Rtuple> ex_pos_serial( string sufile, bool rcv );
  /**
   * This function writes out a SEGY file associated with an MPS space from scratch, 
   * by duplicating a reference trace.
   *
   * @param [in] s_pos - source positions
   * @param [in] N_basis - number of MPS bases (i.e., dimension of MPS space)
   * @param [in] ref_file - reference SEGY filename
   * @param [in] MPS_file - filename of created SEGY
   */
  void write_MPS_SEGY( vector<Rtuple> s_pos,
		       int N_basis,
		       string ref_file,
		       string MPS_file );

  //==================================================================================//
  // The MPS_Space class is an abstract class representing a space of multipole 
  // point sources (MPS's) with respect to a set "basis" of multipole terms.
  //
  // DATA MEMBERS:
  //               MPS_KEYS mks = PARARRAY keys for everything needed in MPS code
  //  Space<ireal> const * core = core SEGYSpace
  //       Vector<Rtuple> s_pos = source positions
  //                  int s_dim = spatial dimension
  //
  // VIRTUAL FUNCTIONS:
  //    clone()
  //    get_type()
  //    get_MPS_basis()
  
  //==================================================================================//
  /**
   * MPS space base class.
   */
  class MPS_Space: public StdSpace<ireal,ireal>{ //ProductSpace<ireal> {

  protected:
    
    MPS_KEYS mks; /**< MPS keys struct, see \ref MPS_KEYS */
    
    //to be set in class constructor or by gen_SEGY_core()
    Space<ireal> const *core; /**< core SEGY space */
    vector<Rtuple> s_pos;     /**< source positions */
    int s_dim;                /**< spatial dimension */
    
    /**
     * Virtual function to be supplied by derived classes. 
     * Supplies a way for copying MPS_Spaces as generic Space pointers.
     */
    virtual Space<ireal> * clone() const =0;

    /**
     * Function for generating MPS SEGY file and setting core.
     */
    void gen_SEGY_core( PARARRAY const& pars );

    /**
     * Convenient function call for running sanity checks post-construction.
     */
    void sanity_checks() const;

  public:
    
    /**
     * Empty constructor.
     */
    MPS_Space()
      : core(NULL),
	s_dim(0),
	s_pos(0) {}

    /**
     * Copy constructor.
     */
    MPS_Space( MPS_Space const &sp );

    /**
     * Main Constructor.
     * MPS_Space main constructor requires the following MPS_KEYS be 
     * initialized, depending on the following cases:
     *
     *  1. MPS file is to be made generated from scratch (make=true);
     *     - MPS_KEYS.MPS_file  - key for filename of MPS to be generated for 
     *                            setting core SEGYSpace
     *     - MPS_KEYS.grid_file - key for RSF file containing spatial grid info
     *     - MPS_KEYS.MPS_ref   - key for SU file used are refence for generating MPS file
     *     - MPS_KEYS.hdr_files - keys for SU files used to extract source positions for
     *                            generating MPS file
     *  2. MPS file exists (make=false);
     *     - MPS_KEYS.MPS_file  - key for filename of MPS used in setting core SEGYSpace
     *     - MPS_KEYS.grid_file - key for RSF file containing spatial grid info
     *    
     * @param [in] _mks - keys related to initialization of MPS spaces, see \ref MPS_KEYS
     * @param [in] pars - input key-val pairs
     * @param [in] make - boolean flag for making core SEGY space from scratch
     */
    MPS_Space( MPS_KEYS _mks, 
	       PARARRAY const& pars,
	       bool make);

    /**
     * Destructor.
     */
    ~MPS_Space();

    //StdSpace functions
    DataContainerFactory const & getDCF() const;
    LinearAlgebraPackage<ireal> const & getLAP() const;

    /**
     * Virtual accessor function to be supplied by subclass.      
     * Returns string containing description of MPS space.
     */
    virtual string get_type() const =0;

    /**
     * Virtual accessor function to be supplied by subclass.
     * Return MPS basis, see \ref MPS_basis.
     */
    virtual vector<MPS_base> get_MPS_basis() const =0;
    
    /**
     * Handy function for retrieving filename of core SEGY file.
     */
    string get_coreFilename() const;

    /**
     * Handy function for accessing derivative order, namely max
     * derivative order for each axis.
     */
    Ituple get_d_ord() const;
    int get_MPS_ord() const;
    
    /**
     * Handy function for accessing dimension (# of MPS bases).
     */
    int get_dim() const{ return this->get_MPS_basis().size(); }

    /**
     * Handy function for accessing spatial dimension.
     */
    int get_s_dim() const{ return s_dim;}

    /**
     * Handy function for accessing source positions.
     */
    vector<Rtuple> get_s_pos() const{ return s_pos;}

    /**
     * Handy function for accessing MPS keys.
     */
    MPS_KEYS get_MPS_keys() const{return mks;}
    
    /**
     * Function for printing extra info of MPS space.
     */
    void print_info(ostream &str) const;
    void print_info() const{ print_info(cerr);};

    /**
     * Function for writing out info of MPS space objects.
     */
    ostream & write(ostream &str) const{
      str << "MPS_Space of type ("<<get_type()<<") with core:\n";
      if(core==NULL)
	str << "NULL core\n";
      else
	core->write(str);
      return str;
    }

  };    

  
}//end TSOpt

#endif
