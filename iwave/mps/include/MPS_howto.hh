// MPS_howto.hh
// Author: Mario J. Bencomo
// last modified 12/07/16

/**
 * \file MPS_howto.hh
 * Detailed example of how to define an MPS space.
 */

#ifndef __MPS_HOWTO_HH_
#define __MPS_HOWTO_HH_

#include "MPS_Space.hh"

namespace TSOpt {

  /**
   * Example scalar (acoustic pressure source) MPS space in 2D/3D:
   *
   * f(t,x) = w_0(t) delta(x-x^*) 
   *        + w_1(t) { (d/dx_0) + (d/dx_1) } delta(x-x^*)
   */
  class ExScal_MPS_Space : public MPS_Space {

  private:

    string type; /**< description of MPS space */
    vector<MPS_base> MPS_basis; /**< MPS basis */ 

    /**
     * Clone function required by base class.
     */
    Space<ireal> * clone() const{return new ExScal_MPS_Space(*this);}

    /**
     * Handy function for building MPS basis.
     */
    void build_MPS_basis();

  public:

    /**
     * Empty constructor.
     */
    ExScal_MPS_Space()
      : MPS_Space(),
	MPS_basis(0){}

    /**
     * Copy constructor.
     */
    ExScal_MPS_Space( ExScal_MPS_Space const &sp)
      : MPS_Space(sp),
	type(sp.type),
	MPS_basis(sp.MPS_basis){}

    /**
     * Main constructor.
     */
    ExScal_MPS_Space( MPS_KEYS _mks,
		      PARARRAY const &pars,
		      bool make=false );
    
    /**
     * Destructor.
     */
    ~ExScal_MPS_Space(){}
    
    /**
     * Accessor function required from base class.
     * Returns tyoe, description of MPS space.
     */
    string get_type() const{ return type;}

    /**
     * Accessor function required from base class.
     * Returns MPS basis.
     */
    vector<MPS_base> get_MPS_basis() const{ return MPS_basis;}

  };


  //----------------------------------------------------------------------//
  ExScal_MPS_Space::ExScal_MPS_Space( MPS_KEYS _mks,
				      PARARRAY const& pars,
				      bool make )
    : MPS_Space(_mks,pars,make) {
  //----------------------------------------------------------------------//
    try{
      
      //setting type
      std::stringstream ss;
      ss << "Example scalar MPS_Space in " << s_dim << "D.";
      type = ss.str();

      //setting MPS_basis and d_ord
      build_MPS_basis();
      	
      //Case where SEGY core is made from scratch.
      //Note, if make==false then base class constructor takes care of 
      //setting SEGYSpace core.
      if( make ){
	gen_SEGY_core(pars);
      }

      sanity_checks();
    }
    catch(RVLException &e){
      e << "ERROR from ExScal_MPS_Space constructor\n";
      throw e;
    }
  }

  //----------------------------------------------------------------------//
  void ExScal_MPS_Space::build_MPS_basis(){
  //----------------------------------------------------------------------//
    try{
      
      //monopole term: b1 = delta(x)
      MPS_base b1; 
      b1.NT = 1; // 1 term
      
      b1.RHS_flags.resize(1); // 1 RHS flag (scalar)
      b1.RHS_flags[0]=1;      // RHS_flags[0]=1 (active in pressure field)
      
      b1.derv.resize(b1.NT); 
      for( int k=0; k<s_dim; k++)
	b1.derv[0].coor[k]=0; // term 0, 0 derivatives in x_k coordinate 
      MPS_basis.push_back(b1);

      //dipole term: b2 = { d/dx_0 + d/dx_1 } delta(x)
      MPS_base b2;
      b2.NT = 2; // 2 terms, (d/dx_0)delta and (d/dx1)delta
      
      b2.RHS_flags.resize(1); // 1 RHS flag (scalar)
      b2.RHS_flags[0]=1;      // RHS_flags[0]=1 (active in pressure field)

      b2.derv.resize(b2.NT);
      b2.derv[0].coor[0]=1;  // term 0, 1 derivative  in x_0 or z coordinate
      b2.derv[0].coor[1]=0;  // term 0, 0 derivatives in x_1 or x coordinate
      b2.derv[0].coor[2]=0;  // term 0, 0 derivatives in x_2 or y coordinate      
      b2.derv[1].coor[0]=0;  // term 1, 0 derivatives in x_0 or z coordinate
      b2.derv[1].coor[1]=1;  // term 1, 1 derivative  in x_1 or x coordinate
      b2.derv[1].coor[2]=0;  // term 1, 0 derivatives in x_2 or y coordinate
      
      MPS_basis.push_back(b2);

    }
    catch(RVLException &e){
      e << "ERROR from ExScal_MPS_Space::build_MPS_basis()!\n";
      throw e;
    }
  }

}//end TSOpt

#endif
