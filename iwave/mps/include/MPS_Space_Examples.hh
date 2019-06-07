// MPS_Space_Examples.hh
// Author: Mario J. Bencomo
// last modified 12/12/16

/**
 * \file MPS_Space_Examples.hh
 * Example MPS spaces.
 */

#ifndef __MPS_CAN_HH_
#define __MPS_CAN_HH_

#include "MPS_Space.hh"

namespace TSOpt {

  //=============================================================================//
  // CanScal_MPS_Space class
  //=============================================================================//
  /**
   * Canonical scalar MPS space.
   *
   * MPS order 0 (2-D), f(t,x) = w_0(t) delta(x-x^*)
   * MPS order 1 (2-D), f(t,x) = w_0(t) delta(x-x^*) 
   *                           + w_1(t) (d/dx_0) delta(x-x^*) 
   *                           + w_2(t) (d/dx_1) delta(x-x^*)
   */
  class CanScal_MPS_Space : public MPS_Space {

  private:

    int order; /**< MPS order */
    string type; /**< description of MPS space */
    vector<MPS_base> MPS_basis; /**< MPS basis */ 
    
    /**
     * Clone function required by base class.
     */
    Space<ireal> * clone() const{return new CanScal_MPS_Space(*this);}
    
    /**
     * Handy function for building MPS basis.
     */
    void build_MPS_basis();

    
  public:

    /**
     * Empty constructor.
     */
    CanScal_MPS_Space()
      : MPS_Space(), 
	order(-1),
	MPS_basis(0){}

    /**
     * Copy constructor.
     */
    CanScal_MPS_Space( CanScal_MPS_Space const &sp)
      : MPS_Space(sp), 
	order(sp.order),
	type(sp.type),
	MPS_basis(sp.MPS_basis){}

    /**
     * Main constructor.
     */
    CanScal_MPS_Space( MPS_KEYS _mks,
		       PARARRAY const &pars,
		       bool make=false );
    
    /**
     * Destructor.
     */
    ~CanScal_MPS_Space(){}

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


  //=============================================================================//
  // CanVec_MPS_Space class
  //=============================================================================//
  /**
   * Canonical vector MPS space.
   *
   * MPS order 0 (2-D), f(x,t) = delta(x-x^*) [w_0(t),w_1(t)]
   * MPS order 1 (2-D), f(x,t) = delta(x-x^*) [w_0(t),w_1(t)]
   *                           + (d/dx_0) delta(x-x^*) [w_2(t),w_3(t)]
   *                           + (d/dx_1) delta(x-x^*) [w_4(t),w_5(t)]
   */
  class CanVec_MPS_Space : public MPS_Space {

  private:

    int order; /**< MPS order */
    string type; /**< description of MPS space */
    vector<MPS_base> MPS_basis; /**< MPS basis */ 

    /**
     * Clone function required by base class.
     */
    Space<ireal> * clone() const{return new CanVec_MPS_Space(*this);}

    /**
     * Handy function for building MPS basis.
     */
    void build_MPS_basis();

  public:

    /**
     * Empty constructor.
     */
    CanVec_MPS_Space()
      : MPS_Space(),
	order(-1),
	MPS_basis(0){}

    /**
     * Copy constructor.
     */
    CanVec_MPS_Space( CanVec_MPS_Space const &sp)
      : MPS_Space(sp), 
	order(sp.order),
	type(sp.type),
	MPS_basis(sp.MPS_basis){}

    /**
     * Main constructor.
     */
    CanVec_MPS_Space( MPS_KEYS mks,
		      PARARRAY const &pars,
		      bool make=false );

    /**
     * Destructor.
     */
    ~CanVec_MPS_Space(){}

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


  //=============================================================================//
  // ExVec_MPS_Space class
  //=============================================================================//
  /**
   * Example 2-D vector MPS space.
   *
   * f(x,t) = delta(x-x^*) [w_0(t),w_1(t)]
   *        + (d/dx_0) delta(x-x^*) [w_2(t),w_2(t)] 
   *        + (d/dx_1) delta(x-x^*) [w_2(t),w_2(t)]
   */
  class ExVec_MPS_Space : public MPS_Space {

  private:

    string type; /**< description of MPS space */
    vector<MPS_base> MPS_basis; /**< MPS basis */ 

    /**
     * Clone function required by base class.
     */
    Space<ireal> * clone() const{return new ExVec_MPS_Space(*this);}

    /**
     * Handy function for building MPS basis.
     */
    void build_MPS_basis();

  public:

    /**
     * Empty constructor.
     */
    ExVec_MPS_Space()
      : MPS_Space(),
	MPS_basis(0){}

    /**
     * Copy constructor.
     */
    ExVec_MPS_Space( ExVec_MPS_Space const &sp)
      : MPS_Space(sp),
	type(sp.type),
	MPS_basis(sp.MPS_basis){}

    /**
     * Main constructor.
     */
    ExVec_MPS_Space( MPS_KEYS mks, 
		     PARARRAY const &pars,
		     bool make=false );

    /**
     * Destructor.
     */
    ~ExVec_MPS_Space(){}

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


  //=============================================================================//
  // Scal_MPS_Space class
  //=============================================================================//
  /**
   * Scalar MPS space.
   *
   * MPS of the form (in 2-D)
   *     f(t,x) = w_0(t) (d/dx_0)^{s_0} (d/dx_1)^{s_1} delta(x-x^*)
   * for a given multi-index s = {s_0,s_1}. Similar for 3-D.
   */
  class Scal_MPS_Space : public MPS_Space {

  private:

    Ituple order; /**< multi-index derivative order*/
    string type; /**< description of MPS space */
    vector<MPS_base> MPS_basis; /**< MPS basis */ 

    /**
     * Clone function required by base class.
     */
    Space<ireal> * clone() const{return new Scal_MPS_Space(*this);}

    /**
     * Handy function for building MPS basis.
     */
    void build_MPS_basis();

  public:

    /**
     * Empty constructor.
     */
    Scal_MPS_Space()
      : MPS_Space(),
	MPS_basis(0){ 
      order.coor[0]=-1;
      order.coor[1]=-1;
      order.coor[2]=-1;
    }

    /**
     * Copy constructor.
     */
    Scal_MPS_Space( Scal_MPS_Space const &sp)
      : MPS_Space(sp),
	order(sp.order),
	type(sp.type),
	MPS_basis(sp.MPS_basis){}

    /**
     * Main constructor.
     */
    Scal_MPS_Space( MPS_KEYS mks, 
		    PARARRAY const &pars,
		    bool make=false );

    /**
     * Destructor.
     */
    ~Scal_MPS_Space(){}

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


}//end TSOpt

#endif
