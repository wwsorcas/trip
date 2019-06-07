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
      f(x) = \inf_{dx} \|F(x)dx - d\|^2 \f$
      as an RVL::Functional. The linear least squares solver is specified by
      policy.
  */
  template<typename Scalar, typename LSPolicy, typename LSPolicyData>
  class VPM: public Functional<Scalar>, public LSPolicy {
        
    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;
        
  private:
      
    LinOpValOp<Scalar> const & op;      // operator
    Vector<Scalar> const & d;           // data
    mutable Vector<Scalar> dx;         // linear solution
    mutable bool applied;
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
	  = LSPolicy::build(dx,lop,d,rnorm,nrnorm,str);
	  
	solver->run();
          
	// get the value of objective function
	val = 0.5*rnorm*rnorm;

	applied = true;
	delete solver;
      }
      catch (RVLException & e) {
	e<<"\ncalled from VPM::apply\n";
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
	e<<"\ncalled from VPM::applyGradient\n";
	throw e;
      }   
    }
        
    void applyHessian(const Vector<Scalar> & x,
		      const Vector<Scalar> & dx,
		      Vector<Scalar> & dy) const {
      RVLException e;
      e<<"ERROR: VPM::applyHessian\n";
      e<<"  Hessian not defined\n";
      throw e;
    }
        
    Functional<Scalar> * clone() const {
      return new VPM<Scalar,LSPolicy,LSPolicyData>(*this);
    }
        
  public:
      
    VPM(LinOpValOp<Scalar> const & _op,
	Vector<Scalar> const & _d,
	LSPolicyData const & s,
	ostream & _str=cerr)
      : LSPolicy(), op(_op), d(_d),
        dx(op.getProductDomain()[1]),
        applied(false), str(_str) {
      try{
	dx.zero();
	LSPolicy::assign(s);
	if (s.verbose) {
	  str<<"\n";
	  str<<"==============================================\n";
	  str<<"VPM constructor - ls policy data = \n";
	  s.write(str);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from VPM::Constructor\n";
	throw e;
      }
    }
      
    VPM(VPM<Scalar,LSPolicy,LSPolicyData> const & f) 
      : LSPolicy(f), op(f.op), d(f.d), dx(f.dx),
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
	e<<"\ncalled from VPM::getMaxStep\n";
	throw e;
      }
    }
      
    Vector<Scalar> const & getLSSoln() const { return dx; }
      
    ostream & write(ostream & str) const {
      str<<"VPM: \n";
      str<<"*** LinearOp valued operator:\n";
      op.write(str);
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
