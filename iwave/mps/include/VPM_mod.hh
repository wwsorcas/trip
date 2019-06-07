// VPM_mod.hh
// Author: Mario J. Bencomo
// last modified 11/10/16

/**
 * \file VPM_mod.hh
 * Modified variable projection method to include preconditioning.
 */

#ifndef __VPM_MOD_HH_
#define __VPM_MOD_HH_

#include "op.hh"
#include "alg.hh"
#include "cgalg.hh"
#include "cgnealg.hh"
#include "terminator.hh"
#include "table.hh"

using namespace RVLAlg;

namespace RVLUmin {
    
  using namespace RVL;
  using namespace RVLAlg;
    
  /** Given a LinOpValOp F and a Vector d in the range of op,
      implements the function
      \f$$
      f(x) = \inf_{dx} \|F(x)dx - d\|^2
      \f$$
      as an RVL::Functional. The linear least squares solver is specified by
      policy. This modified version allows for preconditioning.
  */
  template<typename Scalar, typename LSPolicy, typename LSPolicyData>
  class VPM_mod: public Functional<Scalar>, public LSPolicy {
        
    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;
        
  private:
      
    LinOpValOp<Scalar> const & op; /**< operator */
    LinearOp<Scalar> const & M;    /**< pre-conditioner op */
    Vector<Scalar> const & d;      /**< data */
    mutable Vector<Scalar> dx;     /**< linear solution */
    mutable bool applied;          /**< applied boolean flag */
    ostream & str;
        
  protected:
        
    void apply(const Vector<Scalar> & x,
	       Scalar & val) const {
      try {
	atype rnorm;
	atype nrnorm;
	// access linear part
	LinearRestrictOp<Scalar> lop(op,x);
	dx.zero();
          
	Algorithm * solver
	  = LSPolicy::build(dx,lop,M,d,rnorm,nrnorm,str);
	  
	solver->run();
          
	// get the value of objective function
	val = 0.5*rnorm*rnorm;

	applied = true;
	delete solver;
      }
      catch (RVLException & e) {
	e<<"\ncalled from VPM_mod::apply\n";
	throw e;
      }
    }
        
    void applyGradient(const Vector<Scalar> & x,
		       Vector<Scalar> & g) const {
      try {
	// compute dx, value
	if(!applied){
	  Scalar val;
	  this->apply(x,val);
	}

	// access linear part
	LinearRestrictOp<Scalar> lop(op,x);

	// compute residual dd = F(x)dx-d
	Vector<Scalar> dd(lop.getRange());
	lop.applyOp(dx,dd);	  
	dd.linComb(-1.0,d);

	// access nonlinear part - input product vector with linear
	// part set
	Vector<Scalar> xx(op.getDomain());
	Components<Scalar> cxx(xx);
	cxx[1].copy(dx);
	RestrictOp<Scalar> rop(op,xx,0);

	// access evaluation, linearization
	OperatorEvaluation<Scalar> ropeval(rop,x);

	// naive computation of gradient
	ropeval.getDeriv().applyAdjOp(dd,g);
      }
      catch (RVLException & e) {
	e<<"\ncalled from VPM:_mod:applyGradient\n";
	throw e;
      }   
    }
        
    void applyHessian(const Vector<Scalar> & x,
		      const Vector<Scalar> & dx,
		      Vector<Scalar> & dy) const {
      RVLException e;
      e<<"ERROR: VPM_mod::applyHessian\n";
      e<<"  Hessian not defined\n";
      throw e;
    }
        
    Functional<Scalar> * clone() const {
      return new VPM_mod<Scalar,LSPolicy,LSPolicyData>(*this);
    }
        
  public:
      
    VPM_mod(LinOpValOp<Scalar> const & _op,
	    LinearOp<Scalar> const & _M,
	    Vector<Scalar> const & _d,
	    LSPolicyData const & s,
	    ostream & _str=cerr)
      : LSPolicy(), op(_op), M(_M), d(_d),
        dx(op.getProductDomain()[1]),
        applied(false), str(_str) {
      try{
	dx.zero();
	LSPolicy::assign(s);
	if (s.verbose) {
	  str<<"\n";
	  str<<"==============================================\n";
	  str<<"VPM_mod constructor - ls policy data = \n";
	  s.write(str);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from VPM_mod::Constructor\n";
	throw e;
      }
    }
      
    VPM_mod(VPM_mod<Scalar,LSPolicy,LSPolicyData> const & f) 
      : LSPolicy(f), op(f.op), M(f.M), d(f.d), dx(f.dx),
	applied(f.applied), str(f.str) {}
        
    const Space<Scalar> & getDomain() const { 
      return (op.getProductDomain())[0];
    }
      
    Scalar getMaxStep(const Vector<Scalar> & x,
		      const Vector<Scalar> & dx) const {
      try {
	return op.getMaxStep(x,dx);
      }
      catch (RVLException & e) {
	e<<"\ncalled from VPM_mod::getMaxStep\n";
	throw e;
      }
    }
      
    Vector<Scalar> const & getLSSoln() const { return dx; }
      
    ostream & write(ostream & str) const {
      str<<"VPM_mod: \n";
      str<<"*** LinearOp valued operator:\n";
      op.write(str);
      str <<"*** Preconditioner operator:\n";
      M.write(str);
      str<<"*** data vector:\n";
      d.write(str);
      str<<"*** linear LS solution dx:\n";
      dx.write(str);
      str<<"*** applied = " << applied << "\n";
      return str;
    }
  };

}
#endif
