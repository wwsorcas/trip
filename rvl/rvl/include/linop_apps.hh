/*************************************************************************

Copyright Rice University, 2004, 2005, 2006
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/

#ifndef __RVL_LINOP_APPS
#define __RVL_LINOP_APPS

#include "op.hh"

namespace RVL {

  /** Composition of linear operators \f$ Op_1, Op_2 \mapsto Op_2
      \circ Op_1 \f$ (so subscripts indicate order of evaluation -
      that's how the constructor is organized). Alignment of domains
      and ranges checked as part of construction. */

  template<class Scalar>
  class CompLinearOp: public LinearOp<Scalar> {
  
  private:
  
    std::vector<LinearOp<Scalar> const *> ops;
    mutable bool applied;

  protected:

    virtual LinearOp<Scalar> * clone() const {
      return new CompLinearOp<Scalar>(*this);
    }
      
    void apply(const Vector<Scalar> & x, 
	       Vector<Scalar> & y) const {
      try {
	applied=true;
	if (ops.size()==0) {
	  RVLException e;
	  e<<"ERROR: CompLinearOp::apply\n";
	  e<<"  attempted with unitialized factors\n";
	  throw e;
	}
	else if (ops.size()==1) ops[0]->applyOp(x,y);
	else if (ops.size()==2) {
	  Vector<Scalar> z(ops[0]->getRange());
	  ops[0]->applyOp(x,z);
	  ops[1]->applyOp(z,y);
	}
	else {
	  shared_ptr<Vector<Scalar> > z;
	  shared_ptr<Vector<Scalar> > w = RVL::Vector<Scalar>::newPtr(ops[0]->getRange());
	  ops[0]->applyOp(x,*w);
	  for (size_t i=1; i<ops.size()-1; i++) {
	    z=RVL::Vector<Scalar>::newPtr(ops[i]->getRange());
	    ops[i]->applyOp(*w,*z);
	    w=z;
	  }
	  ops[ops.size()-1]->applyOp(*w,y);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from CompLinearOp::apply\n";
	throw e;
      }
    }

    void applyAdj(const Vector<Scalar> & x, 
		  Vector<Scalar> & y) const {
      try {
	applied=true;
	if (ops.size()==0) {
	  RVLException e;
	  e<<"ERROR: CompLinearOp::applyAdj\n";
	  e<<"  attempted with unitialized factors\n";
	  throw e;
	}
	else if (ops.size()==1) ops[0]->applyAdjOp(x,y);
	else if (ops.size()==2) {
	  Vector<Scalar> z(ops[0]->getRange());
	  ops[1]->applyAdjOp(x,z);
	  ops[0]->applyAdjOp(z,y);
	}
	else {
	  shared_ptr<Vector<Scalar> > z = RVL::Vector<Scalar>::newPtr(ops[ops.size()-1]->getDomain());
	  shared_ptr<Vector<Scalar> > w; 
	  ops[ops.size()-1]->applyAdjOp(x,*z);
	  for (size_t i=ops.size()-2; i>=1; i--) {
	    w=z;
	    z=RVL::Vector<Scalar>::newPtr(ops[i]->getDomain());
	    ops[i]->applyAdjOp(*w,*z);
	  }
	  ops[0]->applyAdjOp(*z,y);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from CompLinearOp::applyAdj\n";
	throw e;
      }
    }

  public:

    CompLinearOp(): ops(0), applied(false) {}
    
    CompLinearOp(LinearOp<Scalar> const & _op1,
		 LinearOp<Scalar> const & _op2): applied(false) {
      try {
	
	ops.push_back(dynamic_cast<RVL::LinearOp<float> const *>
		      (RVL::Operator<float>::export_clone(_op1)));
	ops.push_back(dynamic_cast<RVL::LinearOp<float> const *>
		      (RVL::Operator<float>::export_clone(_op2)));
	if (!ops[0] || !ops[1]) {
	  RVLException e;
	  e<<"ERROR: CompLinearOp constructor\n";
	  e<<"  something seriously wrong - failed to record inputs as LinearOps\n";
	  throw e;
	}

	if (ops[1]->getDomain() != ops[0]->getRange()) {
	  RVLException e;
	  e<<"Error: CompLinearOp constructor\n";
	  e<<"  incompatible domains or ranges\n";
	  e<<"\n";
	  e<<"  **** operator 1:\n";
	  ops[0]->write(e);
	  e<<"\n";
	  e<<"  **** operator 2:\n";
	  ops[1]->write(e);
	  e<<"\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from CompLinearOp::CompLinearOp\n";
	throw e;
      }
    }

    CompLinearOp(CompLinearOp const & op): applied(op.applied) {
      for (size_t i=0; i<op.ops.size();i++) { 
	ops.push_back(dynamic_cast<LinearOp<Scalar> const *>
		      (RVL::Operator<float>::export_clone(*(op.ops[i]))));
      }
    }

    ~CompLinearOp() {
      for (size_t i=0; i<ops.size(); i++) {
	delete ops[i];
      }
    }
    
    void setNext(LinearOp<Scalar> const & lop) {
      if (applied) {
	RVLException e;
	e<<"ERROR: CompLinearOp::setNext\n";
	e<<"  already initialized, used - effectively const\n";
	throw e;
      }
      if (ops.size()>0) {
	if (ops[ops.size()-1]->getRange() != lop.getDomain()) {
	  RVLException e;
	  e<<"ERROR: CompLinearOp::setNext\n";
	  e<<"  input operator domain not same as last added operator range\n";
	  e<<"  domain:\n";
	  lop.getDomain().write(e);
	  e<<"  range:\n";
	  ops[ops.size()-1]->getRange().write(e);
	  throw e;
	}
      }
      ops.push_back(dynamic_cast<RVL::LinearOp<float> const *>
		    (RVL::Operator<float>::export_clone(lop)));
    }
    
    const Space<Scalar> & getDomain() const {
      return ops[0]->getDomain();
    }

    const Space<Scalar> & getRange() const {
      return ops[ops.size()-1]->getRange();
    }
  
    ostream & write(ostream & str) const {
      str<<"CompLinearOp - composition of linear ops\n";
      str<<"\n";
      str<<"in order of application:\n";
      for (size_t i=0; i<ops.size()-1; i++) {
	str<<"**** operator "<<i<<"\n";
	ops[i]->write(str);
	str<<"\n";
      }
      return str;
    }
  };
}

#endif
