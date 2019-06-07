#ifndef __RVLALG_LINFIT_L2_H
#define __RVLALG_LINFIT_L2_H

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
      implements the function \f$
      f(x) = F(x, \mbox{argmin}_{z} \|F(x)y - d\|^2) \f$
      as an RVL::Operator. The linear least squares solver is specified by
      policy.
  */
  template<typename Scalar, typename LSPolicy, typename LSPolicyData>
  class ReduceOp_mod: public Operator<Scalar>, public LSPolicy {
        
    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;
        
  private:

    LinearOp<Scalar> const & M;         // preconditioner operator
    LinOpValOp<Scalar> const & op;      // operator
    Vector<Scalar> const & d;           // data
    mutable Vector<Scalar> xx;          // domain vec
    mutable Components<Scalar> cxx;     // components thereof
    mutable bool applied;               // set on evaluation
    ostream & str;

    void reduce(Vector<Scalar> const & x) const {
      try {
	if (!applied) {
	  cxx[0].copy(x);
	  atype rnorm;
	  atype nrnorm;
	  // access linear part
	  LinearRestrictOp<Scalar> lop(op,cxx[0]);
	  cxx[1].zero();
	  Algorithm * solver
	    = LSPolicy::build(cxx[1],lop,M,d,rnorm,nrnorm,str);
	  solver->run();          
	  // get the value of objective function
	  applied = true;
	  delete solver;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from ReduceOp_mod::reduce\n";
	throw e;
      }
    }
    
  protected:
    
    void apply(Vector<Scalar> const & x,
	       Vector<Scalar> & y) const {
      try {
	reduce(x);
	RVL::Operator<Scalar>::export_apply(op, xx, y);
	y.linComb(-ScalarFieldTraits<Scalar>::One(),d);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ReduceOp_mod::apply\n";
	throw e;
      }
    }
        
    void applyDeriv(Vector<Scalar> const & x,
		    Vector<Scalar> const & dx,
		    Vector<Scalar> & dy) const {
      try {
	reduce(x);
	//RestrictOp<Scalar> rop(op,cxx[1],1); 
	RestrictOp<Scalar> rop(op,xx,0); //MJB edit
	RVL::Operator<Scalar>::export_applyDeriv(rop,x,dx,dy);
	//	RVL::OperatorProductDomain<Scalar>
	//	  ::export_applyPartialDeriv(op, 0, xx, dx, dy);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ReduceOp_mod::applyDeriv\n";
	throw e;
      }
    }

    void applyAdjDeriv(Vector<Scalar> const & x,
		       Vector<Scalar> const & dy,
		       Vector<Scalar> & dx) const {
      try {
	reduce(x);
	//RestrictOp<Scalar> rop(op,cxx[1],1); 
	RestrictOp<Scalar> rop(op,xx,0); //MJB edit
	RVL::Operator<Scalar>::export_applyAdjDeriv(rop,x,dy,dx);
	//RVL::Operator<Scalar>::export_applyAdjDeriv(rop,x,dy,dx);
	//	RVL::OperatorProductDomain<Scalar>
	//	  ::export_applyAdjPartialDeriv(op, 0, xx, dy, dx);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ReduceOp_mod::applyAdjDeriv\n";
	throw e;
      }
    }

    Operator<Scalar> * clone() const {
      return new ReduceOp_mod<Scalar,LSPolicy,LSPolicyData>(*this);
    }
        
  public:
      
    ReduceOp_mod(LinOpValOp<Scalar> const & _op,
		 LinearOp<Scalar> const & _M,
		 Vector<Scalar> const & _d,
		 LSPolicyData const & s,
		 ostream & _str=cerr)
      : LSPolicy(), op(_op), M(_M), d(_d),
	xx(op.getProductDomain()),
	cxx(xx),
	applied(false),
	str(_str) {
      try{
	LSPolicy::assign(s);
	if (s.verbose) {
	  str<<"\n";
	  str<<"==============================================\n";
	  str<<"ReduceOp_mod constructor - ls policy data = \n";
	  s.write(str);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from ReduceOp_mod::Constructor\n";
	throw e;
      }
    }
      
    ReduceOp_mod(ReduceOp_mod<Scalar,LSPolicy,LSPolicyData> const & f) 
      : LSPolicy(f), op(f.op), M(f.M), d(f.d), xx(op.getProductDomain()), cxx(xx),
	applied(false), str(f.str) {}
        
    Space<Scalar> const & getDomain() const { return (op.getProductDomain())[0]; }
    Space<Scalar> const & getRange() const { return op.getRange(); }
      
    Scalar getMaxStep(const Vector<Scalar> & x,
		      const Vector<Scalar> & y) const {
      try {
        reduce(x);
	Vector<Scalar> yy(op.getProductDomain());
	Components<Scalar> cyy(yy);
	cyy[0].copy(y);
	cyy[1].zero();
	return op.getMaxStep(xx,yy);
      }
      catch (RVLException & e) {
	e<<"\ncalled from VPM::getMaxStep\n";
	throw e;
      }
    }
      
    Vector<Scalar> const & getLSSoln() const {
      if (!applied) {
	RVLException e;
	e<<"Error: ReduceOp_mod::getLSSoln()\n";
	e<<"  operator not applied, no reduction performed\n";
	throw e;
      }
      return cxx[1];
    }
      
    ostream & write(ostream & str) const {
      str<<"ReduceOp_mod: \n";
      str<<"*** LinearOp valued operator:\n";
      op.write(str);
      str<<"*** Preconditioner:\n";
      M.write(str);
      str<<"*** data vector:\n";
      d.write(str);
      str<<"*** nonlinear domain vector:\n";
      cxx[0].write(str);
      str<<"*** linear LS solution:\n";
      cxx[1].write(str);
      str<<"*** applied = " << applied << "\n";
      return str;
    }
  };

}
#endif
