/*************************************************************************

Copyright Rice University, 2004, 2005, 2006, 2007, 2008, 2009, 2010
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

#ifndef __RVL_BLOCKOP
#define __RVL_BLOCKOP

#include "op.hh"

namespace RVL {

  /** Operator defined with product domain and range. Partial
      derivatives are defined componentwise, as is access to the
      domain and range as ProductSpaces.  As for the parent class, all
      functions which may change the internal state are protected,
      accessed only by the corresponding OperatorEvaluation objects,
      which act on independent captive instances.
  */
  template<class Scalar> 
  class BlockOperator: public Operator<Scalar> {

    friend class OperatorEvaluation<Scalar>;

  protected:

    virtual void applyComponent(int i, 
				const Vector<Scalar> & x,
				Vector<Scalar> & yi) const = 0;

    virtual void apply(Vector<Scalar> const & x,
		       Vector<Scalar> & y) const {
      try {
	Components<Scalar> yc(y);
	for (int i=0;i<(int)yc.getSize();i++) 
	  applyComponent(i,x,yc[i]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from BlockOperator::apply\n";
	throw e;
      }
    }

    /** \f$dy = \partial_jF_i(x)dx_j\f$, where \f$dx_j \in X_j\f$, \f$dy_i \in Y_i\f$ */
    virtual void applyComponentDeriv(int i, 
				     const Vector<Scalar> & x, 
				     const Vector<Scalar> & dx,
				     Vector<Scalar> & dyi) const = 0;
  
    /** applyDeriv() is implemented in terms of
	applyComponentDeriv(). Default implementation supplied, which
	may be overridden. 
    */
    virtual void applyDeriv(const Vector<Scalar> & x, 
			    const Vector<Scalar> & dx,
			    Vector<Scalar> & dy) const {
      try {
	Components<Scalar> dyc(dy);
	for (int i=0;i<(int)dyc.getSize();i++) 
	  applyComponentDeriv(i,x,dx,dyc[i]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from BlockOperator::applyDeriv\n";
	throw e;
      }
    }

    /** \f$dx_j = \partial_jF_i(x)^*dy_i\f$, where \f$dx_j \in X_j\f$ */
    virtual void applyComponentAdjDeriv(int i, 
					const Vector<Scalar> & x, 
					const Vector<Scalar> & dyi,
					Vector<Scalar> & dx) const = 0;

    /** applyAdjDeriv() is implemented in terms of
	applyComponentAdjDeriv(). Default implementation supplied, which
	may be overridden. */
    virtual void applyAdjDeriv(const Vector<Scalar> & x, 
			       const Vector<Scalar> & dy,
			       Vector<Scalar> & dx) const {
      try {
	Components<Scalar> dyc(dy);
	applyComponentAdjDeriv(0,x,dyc[0],dx);
	if (dyc.getSize()>0) {
	  Vector<Scalar> tmp(this->getDomain(),true);
	  for (int i=1; i<(int)dyc.getSize(); i++) {
	    applyComponentAdjDeriv(i,x,dyc[i],tmp);
	    dx.linComb(1.0,tmp);
	  }
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from BlockOperator::applyAdjDeriv\n";
	throw e;
      }
    }

    /** \f$dy_i = D^2F_i(x)[dx0, dx1]\f$, where \f$dx0, dx1 \in X\f$, \f$dy_i \in Y_i\f$ */
    virtual void applyComponentDeriv2(int i,
				      const Vector<Scalar> & x,
				      const Vector<Scalar> & dx0,
				      const Vector<Scalar> & dx1,
				      Vector<Scalar> & dyi) const = 0;
      
    /** applyDeriv2() is implemented in terms of
	applyComponentDeriv(). Default implementation supplied, which
	may be overridden.
    */
    virtual void applyDeriv2(const Vector<Scalar> & x,
			     const Vector<Scalar> & dx0,
			     const Vector<Scalar> & dx1,
			     Vector<Scalar> & dy) const {
      try {
	Components<Scalar> dyc(dy);
	for (int i=0;i<(int)dyc.getSize();i++)
	  applyComponentDeriv2(i,x,dx0,dx1,dyc[i]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from BlockOperator::applyDeriv2\n";
	throw e;
      }
    }
    /** \f$dx = \sum_i D^2F_i(x)^*[dx,dy_i]\f$, where \f$dy_i \in Y_i\f$ */
    virtual void applyComponentAdjDeriv2(int i,
					 const Vector<Scalar> & x,
					 const Vector<Scalar> & dx0,
					 const Vector<Scalar> & dyi,
					 Vector<Scalar> & dx1) const = 0;
      
    /** applyAdjDeriv() is implemented in terms of
	applyComponentAdjDeriv(). Default implementation supplied, which
	may be overridden. */
    virtual void applyAdjDeriv2(const Vector<Scalar> & x,
				const Vector<Scalar> & dx0,
				const Vector<Scalar> & dy,
				Vector<Scalar> & dx1) const {
      try {
	Components<Scalar> dyc(dy);
	applyComponentAdjDeriv2(0,x,dx0,dyc[0],dx1);
	if (dyc.getSize()>0) {
	  Vector<Scalar> tmp(this->getDomain(),true);
	  for (int i=1; i<(int)dyc.getSize(); i++) {
	    applyComponentAdjDeriv2(i,x,dx0,dyc[i],tmp);
	    dx1.linComb(1.0,tmp);
	  }
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from BlockOperator::applyAdjDeriv2\n";
	throw e;
      }
    }
    /** Primary clone method returns object of this type;
	parent clone method delegates. */
    virtual BlockOperator<Scalar> * cloneBlockOp() const = 0;
    Operator<Scalar> * clone() const { return cloneBlockOp(); }

  public:

    BlockOperator() {}
    BlockOperator(const BlockOperator<Scalar> &) {}
    virtual ~BlockOperator() {}

    /** access to range as ProductSpace */
    virtual const ProductSpace<Scalar> & getProductRange() const = 0;
    /** access to range as Space - delegates to getProductRange */
    const Space<Scalar> & getRange() const { 
      return getProductRange(); 
    }
  
  };

  /** Explicit BlockOp construction for two range components */

  template<typename Scalar>
  class TensorOp: public BlockOperator<Scalar> {

  private:

    Operator<Scalar> const & op1;
    Operator<Scalar> const & op2;
    StdProductSpace<Scalar> rng;

    TensorOp();

  protected:

    void applyComponent(int i, 
			const Vector<Scalar> & x,
			Vector<Scalar> & yi) const {
      try {
	if (i==0) this->export_apply(op1,x,yi);
	else if (i==1) this->export_apply(op2,x,yi);
	else {
	  RVLException e;
	  e<<"Error: TensorOp::applyComponent\n";
	  e<<"index "<<i<<" out of range [0,1]\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TensorOp::applyComponent\n";
	throw e;
      }
    }

    void applyComponentDeriv(int i,
			     const Vector<Scalar> & x, 
			     const Vector<Scalar> & dx,
			     Vector<Scalar> & dyi) const {
      try {
	if (i==0) this->export_applyDeriv(op1,x,dx,dyi);
	else if (i==1) this->export_applyDeriv(op2,x,dx,dyi);
	else {
	  RVLException e;
	  e<<"Error: TensorOp::applyComponentDeriv\n";
	  e<<"index "<<i<<" out of range [0,1]\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TensorOp::applyComponentDeriv\n";
	throw e;
      }
    }

    void applyComponentAdjDeriv(int i, 
				const Vector<Scalar> & x, 
				const Vector<Scalar> & dyi,
				Vector<Scalar> & dx) const {
      try {
	if (i==0) this->export_applyAdjDeriv(op1,x,dyi,dx);
	else if (i==1) this->export_applyAdjDeriv(op2,x,dyi,dx);
	else {
	  RVLException e;
	  e<<"Error: TensorOp::applyComponentAdjDeriv\n";
	  e<<"index "<<i<<" out of range [0,1]\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TensorOp::applyComponentAdjDeriv\n";
	throw e;
      }
    }
  
    void applyComponentDeriv2(int i,
			      const Vector<Scalar> & x,
			      const Vector<Scalar> & dx0,
			      const Vector<Scalar> & dx1,
			      Vector<Scalar> & dyi) const {
      try {
	if (i==0) this->export_applyDeriv2(op1,x,dx0,dx1,dyi);
	else if (i==1) this->export_applyDeriv2(op2,x,dx0,dx1,dyi);
	else {
	  RVLException e;
	  e<<"Error: TensorOp::applyComponentDeriv2\n";
	  e<<"index "<<i<<" out of range [0,1]\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TensorOp::applyComponentDeriv2\n";
	throw e;
      }
    }
      
    void applyComponentAdjDeriv2(int i,
				 const Vector<Scalar> & x,
				 const Vector<Scalar> & dx0,
				 const Vector<Scalar> & dyi,
				 Vector<Scalar> & dx1) const{
      try {
	if (i==0) this->export_applyAdjDeriv2(op1,x,dx0,dyi,dx1);
	else if (i==1) this->export_applyAdjDeriv2(op2,x,dx0,dyi,dx1);
	else {
	  RVLException e;
	  e<<"Error: TensorOp::applyComponentAdjDeriv2\n";
	  e<<"index "<<i<<" out of range [0,1]\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TensorOp::applyComponentAdjDeriv2\n";
	throw e;
      }
    }
      
    TensorOp<Scalar> * cloneTensorOp() const { return new TensorOp(*this); }
    BlockOperator<Scalar> * cloneBlockOp() const { return cloneTensorOp(); }

  public:

    TensorOp(Operator<Scalar> const & _op1,
	     Operator<Scalar> const & _op2)
      : op1(_op1), op2(_op2), rng(op1.getRange(),op2.getRange()) {
      try {
	if (op1.getDomain() != op2.getDomain()) {
	  RVLException e;
	  e<<"Error: TensorOp constructor\n";
	  e<<"input operators do not have same domain\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TensorOp constructor\n";
	throw e;
      }
    }

    TensorOp(TensorOp<Scalar> const & op)
      : op1(op.op1), op2(op.op2), rng(op1.getRange(),op2.getRange()) {}

    ~TensorOp() {}

    Space<Scalar> const & getDomain() const { return op1.getDomain(); }
    ProductSpace<Scalar> const & getProductRange() const { return rng; }

    ostream & write(ostream & str) const {
      str<<"TensorOp: 2x1 block operator\n";
      str<<"*** operator[0,0]:\n";
      op1.write(str);
      str<<"*** operator[1,0]:\n";
      op2.write(str);
      return str;
    }
  };

  /** Explicit BlockOp construction for 2x2 diagonal. Delivered as Op
      rather than BlockOp to avoid any inefficient access to off-diag
      componennts.
  */

  template<typename Scalar>
  class DiagOp: public Operator<Scalar> {

  private:

    std::vector<Operator<Scalar> const * > op;
    StdProductSpace<Scalar> dom;
    StdProductSpace<Scalar> rng;

    bool init() const {
      for (int i=0; i<op.size(); i++)
	if (!(op[i])) return false;
      return true;
    }
    
    DiagOp();

  protected:

    void apply(const Vector<Scalar> & x,
	       Vector<Scalar> & y) const {
      try {
	if (!init()) {
	  RVLException e;
	  e<<"ERROR: DiagOp::apply\n";
	  e<<"  not initialized\n";
	  throw e;
	}
	//	cerr<<"DiagOp::apply\n";
	
	Components<Scalar> cx(x);
	Components<Scalar> cy(y);
	if ((cx.getSize() != op.size()) ||
	    (cy.getSize() != op.size())) {
	  RVLException e;
	  e<<"ERROR: DiagOp::apply\n";
	  e<<"  input or output vector has wrong number of args\n";
	  throw e;
	}
	for (int i=0;i<op.size();i++)
	  RVL::Operator<Scalar>::export_apply(*(op[i]),cx[i],cy[i]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from DiagOp::apply\n";
	throw e;
      }
    }

    void applyDeriv(const Vector<Scalar> & x,
		    const Vector<Scalar> & dx,
		    Vector<Scalar> & dy) const {
      try {
	if (!init()) {
	  RVLException e;
	  e<<"ERROR: DiagOp::applyDeriv\n";
	  e<<"  not initialized\n";
	  throw e;
	}
	Components<Scalar> cx(x);
	Components<Scalar> cdx(dx);
	Components<Scalar> cdy(dy);
	if ((cx.getSize() != op.size()) ||
	    (cdx.getSize() != op.size()) ||
	    (cdy.getSize() != op.size())) {
	  RVLException e;
	  e<<"ERROR: DiagOp::applyDeriv\n";
	  e<<"  input or output vector has wrong number of args\n";
	  throw e;
	}
	for (int i=0; i<op.size(); i++) 
	  RVL::Operator<Scalar>::export_applyDeriv(*(op[i]),cx[i],cdx[i],cdy[i]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from DiagOp::applyDeriv\n";
	throw e;
      }
    }

    void applyAdjDeriv(const Vector<Scalar> & x,
		       const Vector<Scalar> & dy,
		       Vector<Scalar> & dx) const {
      try {
	if (!init()) {
	  RVLException e;
	  e<<"ERROR: DiagOp::applyAdjDeriv\n";
	  e<<"  not initialized\n";
	  throw e;
	}
		
	Components<Scalar> cx(x);
	Components<Scalar> cdx(dx);
	Components<Scalar> cdy(dy);
	if ((cx.getSize() != op.size()) ||
	    (cdx.getSize() != op.size()) ||
	    (cdy.getSize() != op.size())) {
	  RVLException e;
	  e<<"ERROR: DiagOp::applyAdjDeriv\n";
	  e<<"  input or output vector has wrong number of args\n";
	  throw e;
	}	
	for (int i=0; i<op.size(); i++) 
	  RVL::Operator<Scalar>::export_applyAdjDeriv(*(op[i]),cx[i],cdy[i],cdx[i]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from DiagOp::applyAdjDeriv\n";
	throw e;
      }
    }

    void applyDeriv2(const Vector<Scalar> & x,
		     const Vector<Scalar> & dx0,
		     const Vector<Scalar> & dx1,
		     Vector<Scalar> & dy) const {
      try {
	if (!init()) {
	  RVLException e;
	  e<<"ERROR: DiagOp::applyDeriv2\n";
	  e<<"  not initialized\n";
	  throw e;
	}

	Components<Scalar> cx(x);
	Components<Scalar> cdx0(dx0);
	Components<Scalar> cdx1(dx1);
	Components<Scalar> cdy(dy);
	if ((cx.getSize() != op.size()) ||
	    (cdx0.getSize() != op.size()) ||
	    (cdx1.getSize() != op.size()) ||
	    (cdy.getSize() != op.size())) {
	  RVLException e;
	  e<<"ERROR: DiagOp::applyDeriv2\n";
	  e<<"  input or output vector has wrong number of args\n";
	  throw e;
	}
	for (int i=0; i<op.size(); i++) 
	  RVL::Operator<Scalar>::export_applyDeriv2(*(op[i]),cx[i],cdx0[i],cdx1[i],cdy[i]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from DiagOp::applyDeriv2\n";
	throw e;
      }
    }

    void applyAdjDeriv2(const Vector<Scalar> & x,
			const Vector<Scalar> & dx0,
			const Vector<Scalar> & dy,
			Vector<Scalar> & dx1) const {
      try {
	if (!init()) {
	  RVLException e;
	  e<<"ERROR: DiagOp::applyAdjDeriv2\n";
	  e<<"  not initialized\n";
	  throw e;
	}

	Components<Scalar> cx(x);
	Components<Scalar> cdx0(dx0);
	Components<Scalar> cdx1(dx1);
	Components<Scalar> cdy(dy);
	if ((cx.getSize() != op.size()) ||
	    (cdx0.getSize() != op.size()) ||
	    (cdx1.getSize() != op.size()) ||
	    (cdy.getSize() != op.size())) {
	  RVLException e;
	  e<<"ERROR: DiagOp::applyAdjDeriv2\n";
	  e<<"  input or output vector has wrong number of args\n";
	  throw e;
	}
	for (int i=0; i<op.size(); i++) 
	  RVL::Operator<Scalar>::export_applyAdjDeriv2(*(op[i]),cx[i],cdx0[i],cdy[i],cdx1[i]);

      }
      catch (RVLException & e) {
	e<<"\ncalled from DiagOp::applyAdjDeriv2\n";
	throw e;
      }
    }

    Operator<Scalar> * clone() const { return new DiagOp(*this); }

  public:

    DiagOp(size_t nfac)
      : op(nfac, NULL), dom(nfac), rng(nfac) {}

    DiagOp(DiagOp<Scalar> const & x)
      : op(x.op), dom(x.dom), rng(x.rng) {}

    ~DiagOp() {}

    void set(Operator<Scalar> const & a, size_t i) {
      if (op[i]) {
	RVLException e;
	e<<"ERROR: DiagOp::set\n";
	e<<"  entry "<<i<<" already set\n";
	throw e;
      }
      op[i] = RVL::Operator<Scalar>::export_clone(a);
      dom.set(a.getDomain(),i);
      rng.set(a.getRange(),i);
    }
    
    Space<Scalar> const & getDomain() const { return dom; }
    Space<Scalar> const & getRange() const { return rng; }

    ostream & write(ostream & str) const {
      str<<"DiagOp: diagonal block operator\n";
      if (this->init()) {
	for (int i=0; i<op.size(); i++) {
	  str<<"*** operator["<<i<<","<<i<<"]:\n";
	  op[i]->write(str);
	}
      }
      else {
	str<<"  incompletely initialized\n";
      }
      return str;
    }
  };

  // forward declaration
  template<typename Scalar>
  class BlockLinearOpBlock;

  /** Linear Operator defined with product range.
      Y.H. at Oct 22, 2014 
      Since there is no LinOpProdDom, need class with 
      product domain as well - WWS 03.03.15
  */
  template<class Scalar>
  class BlockLinearOp: public LinearOp<Scalar> {
        
    friend class BlockLinearOpBlock<Scalar>;

  protected:
        
    /** apply A_{i,j} */
    virtual void apply(int i, int j,
		       const Vector<Scalar> & xj,
		       Vector<Scalar> & yi) const = 0;
      
    virtual void apply(Vector<Scalar> const & x,
		       Vector<Scalar> & y) const {
      try {
	Components<Scalar> yc(y);
	Components<Scalar> xc(x);
	for (int i=0;i<(int)yc.getSize();i++)  {
	  yc[i].zero();
	  Vector<Scalar> tmp(yc[i].getSpace());
	  for (int j=0; j< (int)xc.getSize(); j++) {
	    apply(i,j,xc[j],tmp);
	    yc[i].linComb(ScalarFieldTraits<Scalar>::One(),tmp);
	  }
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from BlockLinearOp::apply\n";
	throw e;
      }
    }
   
    /** apply A_{i,j}^* */
    virtual void applyAdj(int i, int j,
			  const Vector<Scalar> & yi,
			  Vector<Scalar> & xj) const = 0;
        
    /** applyAdj() is implemented in terms of
	applyComponentAdj(). Default implementation supplied, which
	may be overridden. */
    virtual void applyAdj(const Vector<Scalar> & y,
			  Vector<Scalar> & x) const {
      try {
	Components<Scalar> xc(x);
	Components<Scalar> yc(y);
	for (int j=0; j<(int)xc.getSize();j++) {
	  xc[j].zero();
	  Vector<Scalar> tmp(xc[j].getSpace());
	  for (int i=0; i<(int)yc.getSize(); i++) {
	    applyAdj(i,j,yc[i],tmp);
	    xc[j].linComb(ScalarFieldTraits<Scalar>::One(),tmp);
	  }
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from BlockLinearOp::applyAdjDeriv\n";
	throw e;
      }
    }

    /** Primary clone method returns object of this type;
	parent clone method delegates. */
    virtual BlockLinearOp<Scalar> * cloneBlockLinearOp() const = 0;
    LinearOp<Scalar> * clone() const { return cloneBlockLinearOp(); }
        
  public:
        
    BlockLinearOp() {}
    BlockLinearOp(const BlockLinearOp<Scalar> &) {}
    virtual ~BlockLinearOp() {}
        
    /** access to domain as ProductSpace */
    virtual const ProductSpace<Scalar> & getProductDomain() const = 0;
    /** access to domain as Space - delegates to getProductDomain */
    const Space<Scalar> & getDomain() const { 
      return getProductDomain(); 
    }
    /** access to range as ProductSpace */
    virtual const ProductSpace<Scalar> & getProductRange() const = 0;
    /** access to range as Space - delegates to getProductRange */
    const Space<Scalar> & getRange() const { 
      return getProductRange(); 
    }
        
  };

  /** i,j block of BlockLinearOp, as LinearOp
  */
  template<class Scalar>
  class BlockLinearOpBlock: public LinearOp<Scalar> {

  private:

    BlockLinearOp<Scalar> const & blk;
    int row;
    int col;

  protected:
        
    virtual void apply(const Vector<Scalar> & xj,
		       Vector<Scalar> & yi) const {
      try {
	blk.apply(row, col, xj, yi);
      }      
      catch (RVLException & e) {
	e<<"\ncalled from BlockLinearOpBlock::apply\n";
	throw e;
      }
    }
      
    virtual void applyAdj(const Vector<Scalar> & yi,
			  Vector<Scalar> & xj) const {
      try {
	blk.applyAdj(row, col, yi, xj);
      }
      catch (RVLException & e) {
	e<<"\ncalled from BlockLinearOpBlock::applyAdj\n";
	throw e;
      }
    }

    LinearOp<Scalar> * clone() const { return BlockLinearOpBlock(*this); }
        
  public:
        
    BlockLinearOpBlock(BlockLinearOp<Scalar> const & _blk, int _row, int _col): blk(_blk), row(_row), col(_col) {}
    BlockLinearOpBlock(const BlockLinearOpBlock<Scalar> & b): blk(b.blk), row(b.row), col(b.col) {}
    ~BlockLinearOpBlock() {}
        
    /** access to domain as Space - delegates to getProductDomain */
    const Space<Scalar> & getDomain() const { 
      return blk.getProductDomain()[col]; 
    }
    /** access to range as Space - delegates to getProductRange */
    const Space<Scalar> & getRange() const { 
      return blk.getProductRange()[row]; 
    }
        
  };

  /** Linear Operator defined with product range.
      Y.H. at Oct 22, 2014 
      Renamed - WWS 03.03.15
  */
  template<class Scalar>
  class ColumnLinearOp: public LinearOp<Scalar> {
        
    friend class OperatorEvaluation<Scalar>;
        
  protected:
        
    virtual void applyComponent(int i,
		       const Vector<Scalar> & x,
		       Vector<Scalar> & yi) const = 0;
      
    virtual void apply(Vector<Scalar> const & x,
		       Vector<Scalar> & y) const {
      try {
	Components<Scalar> yc(y);
	for (int i=0;i<(int)yc.getSize();i++) 
	  applyComponent(i,x,yc[i]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ColumnLinearOp::apply\n";
	throw e;
      }
    }

        
    /** \f$x_j = F_i^*y_i\f$, where \f$x_j \in X_j\f$ */
    virtual void applyComponentAdj(int i,
				   const Vector<Scalar> & yi,
				   Vector<Scalar> & x) const = 0;
        
    /** applyAdj() is implemented in terms of
	applyComponentAdj(). Default implementation supplied, which
	may be overridden. */
    virtual void applyAdj(const Vector<Scalar> & x,
			  Vector<Scalar> & y) const {
      try {
	Components<Scalar> xc(x);
	applyComponentAdj(0,xc[0],y);
	if (xc.getSize()>0) {
	  Vector<Scalar> tmp(this->getDomain(),true);
	  for (int i=1; i<(int)xc.getSize(); i++) {
	    applyComponentAdj(i,xc[i],tmp);
	    y.linComb(1.0,tmp);
	  }
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from ColumnLinearOp::applyAdjDeriv\n";
	throw e;
      }
    }

    /** Primary clone method returns object of this type;
	parent clone method delegates. */
    virtual ColumnLinearOp<Scalar> * cloneColumnLinearOp() const = 0;
    LinearOp<Scalar> * clone() const { return cloneColumnLinearOp(); }
        
  public:
        
    ColumnLinearOp() {}
    ColumnLinearOp(const ColumnLinearOp<Scalar> &) {}
    virtual ~ColumnLinearOp() {}
        
    /** access to range as ProductSpace */
    virtual const ProductSpace<Scalar> & getProductRange() const = 0;
    /** access to range as Space - delegates to getProductRange */
    const Space<Scalar> & getRange() const { 
      return getProductRange(); 
    }
        
  };
    
  /** Explicit ColumnLinearOp construction for two range components
      Y.H. at Oct 22, 2014
  */
    
  template<typename Scalar>
  class TensorLinearOp: public ColumnLinearOp<Scalar> {
        
  private:
        
    LinearOp<Scalar> const & op1;
    LinearOp<Scalar> const & op2;
    StdProductSpace<Scalar> rng;
        
    TensorLinearOp();
        
  protected:
        
    void applyComponent(int i,
			const Vector<Scalar> & x,
			Vector<Scalar> & yi) const {
      try {
	if (i==0) this->export_apply(op1,x,yi);
	else if (i==1) this->export_apply(op2,x,yi);
	else {
	  RVLException e;
	  e<<"Error: TensorLinearOp::applyComponent\n";
	  e<<"index "<<i<<" out of range [0,1]\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TensorLinearOp::applyComponent\n";
	throw e;
      }
    }
        
    void applyComponentAdj(int i,
			   const Vector<Scalar> & yi,
			   Vector<Scalar> & x) const {
      try {
	if (i==0) op1.applyAdjOp(yi,x); //this->export_applyAdj(op1,yi,x);
	else if (i==1) op2.applyAdjOp(yi,x); //this->export_applyAdj(op2,yi,x);
	else {
	  RVLException e;
	  e<<"Error: TensorLinearOp::applyComponentAdj\n";
	  e<<"index "<<i<<" out of range [0,1]\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TensorLinearOp::applyComponentAdj\n";
	throw e;
      }
    }
        
    TensorLinearOp<Scalar> * cloneTensorLinearOp() const { return new TensorLinearOp(*this); }
    ColumnLinearOp<Scalar> * cloneColumnLinearOp() const { return cloneTensorLinearOp(); }
        
  public:
        
    TensorLinearOp(LinearOp<Scalar> const & _op1,
		   LinearOp<Scalar> const & _op2)
      : op1(_op1), op2(_op2), rng(op1.getRange(),op2.getRange()) {
      try {
	if (op1.getDomain() != op2.getDomain()) {
	  RVLException e;
	  e<<"Error: TensorLinearOp constructor\n";
	  e<<"input operators do not have same domain\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TensorLinearOp constructor\n";
	throw e;
      }
    }
        
    TensorLinearOp(TensorLinearOp<Scalar> const & op)
      : op1(op.op1), op2(op.op2), rng(op1.getRange(),op2.getRange()) {}
        
    ~TensorLinearOp() {}
        
    Space<Scalar> const & getDomain() const { return op1.getDomain(); }
    ProductSpace<Scalar> const & getProductRange() const { return rng; }
        
    ostream & write(ostream & str) const {
      str<<"TensorLinearOp: 2x1 block operator\n";
      str<<"*** LinearOp[0,0]:\n";
      op1.write(str);
      str<<"*** LinearOp[1,0]:\n";
      op2.write(str);
      return str;
    }
  };

  /** LinOpValOp defined with product range.
      By MJB, 03.08.18 
  */
  template<class Scalar>
  class ColumnLinOpValOp: public LinOpValOp<Scalar> {
        
    friend class OperatorEvaluation<Scalar>;
        
  protected:

    virtual void apply0Component(int i,
				 const Vector<Scalar> & x0,
				 const Vector<Scalar> & x1,
				 Vector<Scalar> & yi) const = 0;

    virtual void apply0(const Vector<Scalar> & x0,
			const Vector<Scalar> & x1,
			Vector<Scalar> & y) const {
      try {
	Components<Scalar> yc(y);
	for (int i=0;i<(int)yc.getSize();i++) 
	  apply0Component(i,x0,x1,yc[i]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ColumnLinOpValOp::apply0\n";
	throw e;
      }
    }

    virtual void applyAdj0Component(int i,
				    const Vector<Scalar> & x0,
				    const Vector<Scalar> & yi,
				    Vector<Scalar> & x1) const = 0;
        
    virtual void applyAdj0(const Vector<Scalar> & x0,
			   const Vector<Scalar> & y,
			   Vector<Scalar> & x1) const {
      try {
	Components<Scalar> yc(y);
	applyAdj0Component(0,x0,yc[0],x1);
	if (yc.getSize()>0) {
	  Vector<Scalar> tmp(this->getDomain(),true);
	  for (int i=1; i<(int)yc.getSize(); i++) {
	    applyAdj0Component(i,x0,yc[i],tmp);
	    x1.linComb(1.0,tmp);
	  }
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from ColumnLinOpValOp::applyAdj0\n";
	throw e;
      }
    }

    virtual void applyPartialDeriv0Component(int i,
					     const Vector<Scalar> & x0,
					     const Vector<Scalar> & x1,
					     const Vector<Scalar> & dx0,
					     Vector<Scalar> & dyi) const = 0;

    virtual void applyPartialDeriv0(const Vector<Scalar> & x0,
				    const Vector<Scalar> & x1,
				    const Vector<Scalar> & dx0,
				    Vector<Scalar> & dy) const  {
      try {
	Components<Scalar> dyc(dy);
	for (int i=0;i<(int)dyc.getSize();i++) 
	  applyPartialDeriv0Component(i,x0,x1,dx0,dyc[i]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ColumnLinOpValOp::applyPartialDeriv0\n";
	throw e;
      }
    }

    virtual void applyAdjPartialDeriv0Component(int i,
						const Vector<Scalar> & x0,
						const Vector<Scalar> & x1,
						const Vector<Scalar> & yi,
						Vector<Scalar> & dx0) const = 0;
        
    virtual void applyAdjPartialDeriv0(const Vector<Scalar> & x0,
				       const Vector<Scalar> & x1,
				       const Vector<Scalar> & y,
				       Vector<Scalar> & dx0) const {
      try {
	Components<Scalar> yc(y);
	applyAdjPartialDeriv0Component(0,x0,x1,yc[0],dx0);
	if (yc.getSize()>0) {
	  Vector<Scalar> tmp(this->getDomain(),true);
	  for (int i=1; i<(int)yc.getSize(); i++) {
	    applyAdjPartialDeriv0Component(i,x0,x1,yc[i],tmp);
	    dx0.linComb(1.0,tmp);
	  }
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from ColumnLinOpValOp::applyAdjPartialDeriv0\n";
	throw e;
      }
    }
        
  public:
        
    ColumnLinOpValOp() {}
    ColumnLinOpValOp(const ColumnLinOpValOp<Scalar> &) {}
    virtual ~ColumnLinOpValOp() {}
        
    /** access to range as ProductSpace */
    virtual const ProductSpace<Scalar> & getProductRange() const = 0;
    /** access to range as Space - delegates to getProductRange */
    const Space<Scalar> & getRange() const { 
      return getProductRange(); 
    }
        
  };
    
  /** Explicit ColumnLinOpVal construction for two range components
      By MJB, 03.08.18
  */
    
  template<typename Scalar>
  class TensorLinOpValOp: public ColumnLinOpValOp<Scalar> {
    
  private:
        
    LinOpValOp<Scalar> const & op1;
    LinOpValOp<Scalar> const & op2;
    StdProductSpace<Scalar> rng;
        
    TensorLinOpValOp();
        
  protected:
        
    void apply0Component(int i,
			 const Vector<Scalar> & x0,
			 const Vector<Scalar> & x1,
			 Vector<Scalar> & yi) const {
      try {
	if (i==0) this->export_apply0(op1,x0,x1,yi);
	else if (i==1) this->export_apply0(op2,x0,x1,yi);
	else {
	  RVLException e;
	  e<<"Error: TensorLinOpValOp::apply0Component\n";
	  e<<"index "<<i<<" out of range [0,1]\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TensorLinOpValOp::apply0Component\n";
	throw e;
      }
    }
    
    void applyAdj0Component(int i,
			    const Vector<Scalar> & x0,
			    const Vector<Scalar> & yi,
			    Vector<Scalar> & x1) const {
      try {
	if (i==0) this->export_applyAdj0(op1,x0,yi,x1);
	else if (i==1) this->export_applyAdj0(op2,x0,yi,x1);
	else {
	  RVLException e;
	  e<<"Error: TensorLinOpValOp::applyAdj0Component\n";
	  e<<"index "<<i<<" out of range [0,1]\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TensorLinOpValOp::applyAdj0Component\n";
	throw e;
      }
    }

   void applyPartialDeriv0Component(int i,
				    const Vector<Scalar> & x0,
				    const Vector<Scalar> & x1,
				    const Vector<Scalar> & dx0,
				    Vector<Scalar> & dyi) const {
     try {
       if (i==0) this->export_applyPartialDeriv0(op1,x0,x1,dx0,dyi);
       else if (i==1) this->export_applyPartialDeriv0(op2,x0,x1,dx0,dyi);
       else {
	 RVLException e;
	 e<<"Error: TensorLinOpValOp::applyPartialDeriv0Component\n";
	 e<<"index "<<i<<" out of range [0,1]\n";
	 throw e;
       }       
     }
     catch (RVLException & e) {
	e<<"\ncalled from TensorLinOpValOp::applyPartialDeriv0Component\n";
	throw e;       
     }
   }

    void applyAdjPartialDeriv0Component(int i,
					const Vector<Scalar> & x0,
					const Vector<Scalar> & x1,
					const Vector<Scalar> & dyi,
					Vector<Scalar> & dx0) const {
      try {
	if (i==0) this->export_applyAdjPartialDeriv0(op1,x0,x1,dyi,dx0);
	else if (i==1) this->export_applyAdjPartialDeriv0(op2,x0,x1,dyi,dx0);
	else {
	  RVLException e;
	  e<<"Error: TensorLinOpValOp::applyAdjPartialDeriv0Component\n";
	  e<<"index "<<i<<" out of range [0,1]\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TensorLinOpValOp::applyAdjPartialDeriv0Component\n";
	throw e;
      }
    }

    OperatorProductDomain<Scalar> * clonePD() const{
      return new TensorLinOpValOp(*this);
    }
        
  public:
        
    TensorLinOpValOp(LinOpValOp<Scalar> const & _op1,
		     LinOpValOp<Scalar> const & _op2)
      : op1(_op1), op2(_op2), rng(op1.getRange(),op2.getRange()) {
      try {
	if (op1.getDomain() != op2.getDomain()) {
	  RVLException e;
	  e<<"Error: TensorLinOpValOp constructor\n";
	  e<<"input operators do not have same domain\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TensorLinOpValOp constructor\n";
	throw e;
      }
    }
        
    TensorLinOpValOp(TensorLinOpValOp<Scalar> const & op)
      : op1(op.op1), op2(op.op2), rng(op1.getRange(),op2.getRange()) {}
        
    ~TensorLinOpValOp() {}
        
    Space<Scalar> const & getDomain() const { return op1.getDomain();}
    ProductSpace<Scalar> const & getProductDomain() const{ return op1.getProductDomain();}
    ProductSpace<Scalar> const & getProductRange() const { return rng;}
        
    ostream & write(ostream & str) const {
      str<<"TensorLinOpValOp: 2x1 block operator\n";
      str<<"*** LinOpValOp[0,0]:\n";
      op1.write(str);
      str<<"*** LinOpValOp[1,0]:\n";
      op2.write(str);
      return str;
    }
  };

  /* FO implementation. 

     FOs arranged in BlockFO form. For an n x m BlockOpFO, the apply
     FOs are stored in a BlockFO of length n. The nm partial derivs
     are stored in a length-n std::vector of BlockFOs of length m, and
     adjoint partial derivatives in a length-m std::vector of BlockFOs
     of length n (i.e. the transpose structure). Because BlockFOs may
     themselves store persistent state, the vectors for deriv and adj
     deriv store pointers.

     Each FO in the apply BlockFO should define evaluation for argument
     vectors of length m.

     Deriv and adj deriv FOs should define evaluation for vectors of
     length m+1, and arg vector should take form
     
     (x[0],...,x[m-1],dx)

     in which dx is a perturbation of one of the jth input component
     for applyPartialDeriv(i,j,...), or of the ith output component
     for applyAdjPartialDeriv(i,j,...).

     Additionally provide for extra parameters in the form of
     RVL::Vectors, supplied in the form of a std::vector of
     pointers. These extra Vectors do not participate in domain or
     range definition, but simply parametrize the action of the
     FOs. Non-Vector parameters should be passed as member data of the
     FOs.
  */

  /*
    template<typename Scalar>
    class BlockOpFO: public BlockOperator<Scalar> {

    private:

    ProductSpace<Scalar> const & dom;
    ProductSpace<Scalar> const & rng;
    std::vector<FunctionObject *> & f; // F_i
    std::vector<std::vector<FunctionObject *> *> & dff; // DF_i/Dx_j
    std::vector<std::vector<FunctionObject *> *> & dfa; // DF_i/Dx_j^T
    std::vector<RVL::Vector<Scalar> const *> & par; // parameters
    BlockOpFO();

    protected:

    virtual void applyComponent(int i, 
    const Vector<Scalar> & x,
    Vector<Scalar> & yi) const {
    try {
    // standard sanity check
    if (x.getSpace() != dom) {
    RVLException e;
    e<<"Error: BlockFO::applyComponent\n";
    e<<"input reference arg not in domain\n";
    throw e;
    }
    if (yi.getSpace() != rng[i]) {
    RVLException e;
    e<<"Error: BlockFO::applyComponent\n";
    e<<"input arg not in component "<<i<<" of range\n";
    throw e;
    }
    Components<Scalar> xc(x);
    std::vector< RVL::Vector<Scalar> const * > xv(xc.getSize()+par.size());
    for (int j=0;j<xc.getSize();j++) xv[j]=&xc[j];
    for (int j=0;j<par.size();j++) xv[j+xc.getSize()]=par[j];
    yi.eval(*(f.at(i)),xv);
    }
    catch (RVLException & e) {
    e<<"\ncalled from BlockOpFO::apply\n";
    throw e;
    }
    catch (out_of_range) {
    RVLException e;
    e<<"Error: out-of-range exception in BlockOpFO::apply\n";
    throw e;
    }
    }

    virtual void applyPartialDeriv(int i, int j,
    const Vector<Scalar> & x, 
    const Vector<Scalar> & dxj,
    Vector<Scalar> & dyi) const {
    try {
    // standard sanity check
    if (x.getSpace() != dom) {
    RVLException e;
    e<<"Error: BlockFO::applyPartialDeriv\n";
    e<<"input reference arg not in domain\n";
    throw e;
    }
    if (dyi.getSpace() != rng[i]) {
    RVLException e;
    e<<"Error: BlockFO::applyPartialDeriv\n";
    e<<"output pert arg not in component "<<i<<" of range\n";
    throw e;
    }
    if (dxj.getSpace() != dom[j]) {
    RVLException e;
    e<<"Error: BlockFO::applyPartialDeriv\n";
    e<<"input pert arg not in component "<<j<<" of domain\n";
    throw e;
    }

    Components<Scalar> xc(x);
    std::vector<RVL::Vector<Scalar> const * > xv(xc.getSize()+1+par.size());
    for (int k=0;k<xc.getSize();k++) xv[k]=&xc[k];
    xv[xc.getSize()]=&dxj;
    for (int k=0;k<par.size();k++) xv[k+xc.getSize()+1]=par[k];
    dyi.eval(*((dff.at(i))->at(j)),xv);
    }
    catch (RVLException & e) {
    e<<"\ncalled from BlockOpFO::applyPartialDeriv\n";
    throw e;
    }
    catch (out_of_range) {
    RVLException e;
    e<<"Error: out-of-range exception in BlockOpFO::applyPartialDeriv\n";
    throw e;
    }
    }

    virtual void applyAdjPartialDeriv(int i, int j,
    const Vector<Scalar> & x, 
    const Vector<Scalar> & dyi,
    Vector<Scalar> & dxj) const {
    try {
    // standard sanity check
    if (x.getSpace() != dom) {
    RVLException e;
    e<<"Error: BlockFO::applyAdjPartialDeriv\n";
    e<<"input reference arg not in domain\n";
    throw e;
    }
    if (dyi.getSpace() != rng[i]) {
    RVLException e;
    e<<"Error: BlockFO::applyAdjPartialDeriv\n";
    e<<"input pert arg not in component "<<i<<" of range\n";
    throw e;
    }
    if (dxj.getSpace() != dom[j]) {
    RVLException e;
    e<<"Error: BlockFO::applyAdjPartialDeriv\n";
    e<<"output pert arg not in component "<<j<<" of domain\n";
    throw e;
    }

    Components<Scalar> xc(x);
    std::vector<RVL::Vector<Scalar> const * > xv(xc.getSize()+1+par.size());
    for (int k=0;k<xc.getSize();k++) xv[k]=&xc[k];
    xv[xc.getSize()]=&dyi;
    for (int k=0;k<par.size();k++) xv[k+xc.getSize()+1]=par[k];
    dxj.eval(*((dfa.at(j))->at(i)),xv);
    }
    catch (RVLException & e) {
    e<<"\ncalled from BlockOpFO::applyAdjPartialDeriv\n";
    throw e;
    }
    catch (out_of_range) {
    RVLException e;
    e<<"Error: out-of-range exception in BlockOpFO::applyAdjPartialDeriv\n";
    throw e;
    }
    }

    virtual BlockOperator<Scalar> * cloneBlockOp() const {
    return new BlockOpFO<Scalar>(*this);
    }

    public:

    BlockOpFO(ProductSpace<Scalar> const & _dom,
    ProductSpace<Scalar> const & _rng,
    std::vector< FunctionObject *> const & _f,
    std::vector< std::vector< FunctionObject *> *> const & _dff,
    std::vector< std::vector< FunctionObject *> *> const & _dfa,
    std::vector< Vector<Scalar> const *> const & _par)
    : dom(_dom), rng(_rng), f(_f), dff(_dff), dfa(_dfa), par(_par) {
    // sanity tests
    int m=dom.getSize();
    int n=rng.getSize();

    bool testdim = true;
    testdim = testdim && (n == f.size());
    testdim = testdim && (n == dff.size());
    for (int i=0;i<n;i++) 
    testdim = testdim && (m == (dff.at(i))->size());
    testdim = testdim && (m == dfa.size());
    for (int j=0;j<m;j++) 
    testdim = testdim && (n == (dfa.at(j))->size());

    if (!testdim) {
    RVLException e;
    e<<"Error: BlockOpFO construction failed, rows="<<n<<" cols="<<m<<"\n";
    e<<"incompatible input FO vectors\n";
    throw e;
    }
    }

    BlockOpFO(const BlockOpFO<Scalar> & a)
    : dom(a.dom), rng(a.rng), f(a.f), dff(a.dff), dfa(a.dfa), par(a.par) {
    // copy parameter vectors
    //      for (int i=0;i<a.par.size();i++) par[i]=a.par[i];
    }

    ~BlockOpFO() {}
   
    virtual const ProductSpace<Scalar> & getProductDomain() const { return dom; }
    virtual const ProductSpace<Scalar> & getProductRange() const { return rng; }
    };

  */

  /** Affine Injection operator. Takes vector in a product of component 
      (sub)spaces as input, replaces corresponding components of 
      output with input. Replaces all other components with constant
      data members.

      Usual linear injection and projection operators are the derivative
      and adjoint deriv of this op.

      For convenience of definition, takes a const reference to a vector 
      in the product range space as arg to constructor, initializes a 
      copy as data member. The apply method copies this vector to the output,
      then replaces the indicated components with corresponding input 
      components. 

      Relies on error trapping in component subspace operations - implements
      no additional traps.

      NOTE added 04.09.15 WWS: because of the current limitations of
      ProductSpace construction, number of components limited to 4,
      with clunky constructor interfaces.
  */

  template<typename Scalar>
  class InjectOp: public Operator<Scalar> {
    
  private:

    std::shared_ptr<Vector<Scalar> > ref;
    Components<Scalar> cref;
    StdProductSpace<Scalar> dom;
    std::vector<int> icvec;
    
  protected:
    
    void apply(Vector<Scalar> const & x,
	       Vector<Scalar> & y) const {
      // rely on built-in error trapping
      try {
	//	cerr<<"InjectOp::apply\n";
	y.copy(*ref);
	Components<Scalar> cy(y);
	Components<Scalar> cx(x);
	for (int i=0; i<icvec.size(); i++) 
	  cy[icvec[i]].copy(cx[i]);
      }
      catch (RVLException &e) {
	e<<"\ncalled from InjectOp::apply\n";
	throw e;
      }
    }

    void applyDeriv(Vector<Scalar> const & x,
		    Vector<Scalar> const & dx,
		    Vector<Scalar> & dy) const {
      try {
	dy.zero();
	Components<Scalar> cx(dx);
	Components<Scalar> cy(dy);
	for (int i=0; i<icvec.size(); i++) 
	  cy[icvec[i]].copy(cx[i]);
      }
      catch (RVLException &e) {
	e<<"\ncalled from InjectOp::applyDeriv\n";
	throw e;
      }
    }
    
    void applyAdjDeriv(Vector<Scalar> const & x,
		       Vector<Scalar> const & dy,
		       Vector<Scalar> & dx) const {
      try {
	Components<Scalar> cy(dy);
	Components<Scalar> cx(dx);
	for (int i=0; i<icvec.size(); i++) 
	  cx[i].copy(cy[icvec[i]]);
      }
      catch (RVLException &e) {
	e<<"\ncalled from InjectOp::applyAdjDeriv\n";
	throw e;
      }
    }

    Operator<Scalar> * clone() const { return new InjectOp(*this); }

  public:

    InjectOp(std::shared_ptr<Vector<Scalar> > _ref,
	     int _icomp)
      : ref(_ref), cref(*ref),
	dom(cref[_icomp].getSpace()),
	icvec(1,_icomp) {}

    InjectOp(std::shared_ptr<Vector<Scalar> > _ref,
	     std::vector<int> _icvec)
      : ref(_ref), cref(*ref),
	dom(_icvec.size()), 
	icvec(_icvec)
    {
      for (int i=0; i<icvec.size(); i++)
	dom.set(cref[i].getSpace(),i);
    }
    
    InjectOp(InjectOp<Scalar> const & f)
      : ref(f.ref), cref(*ref), dom(f.dom), icvec(f.icvec) {}

    ~InjectOp() {}

    Space<Scalar> const & getDomain() const {
      return dom;
    }
    Space<Scalar> const & getRange() const {
      return ref->getSpace();
    }

    ostream & write(ostream & str) const {
      str<<"InjectOp object\n";
      str<<"domain:\n***********\n";
      dom.write(str);
      str<<"range:\n***********\n";
      this->getRange().write(str);
      str<<"components of ref vector:\n***********\n";
      for (int i=0; i<icvec.size(); i++) str<<"  index = "<<i<<" value = "<<icvec[i]<<"\n";
      str<<"reference vector in output product space:\n";
      ref->write(str);
      return str;
    }
  };

  /** convenience specialization for cases where ref vector plays no
      role. Application of InjectOp is equivalent to setting range
      vector to ref data member and applying LinInjectOp.
  */
  template<typename Scalar>
  class LinInjectOp: public LinearOp<Scalar> {
    
  private:

    Space<Scalar> const & dom;
    Space<Scalar> const & rng;
    std::vector<int> icvec;
    
  protected:
    
    void apply(Vector<Scalar> const & x,
	       Vector<Scalar> & y) const {
      // rely on built-in error trapping
      try {
	y.zero();
	Components<Scalar> cy(y);
	Components<Scalar> cx(x);
	for (int i=0; i<icvec.size(); i++) 
	  cy[icvec[i]].copy(cx[i]);
      }
      catch (RVLException &e) {
	e<<"\ncalled from LinInjectOp::apply\n";
	throw e;
      }
    }

    void applyAdj(Vector<Scalar> const & y,
		  Vector<Scalar> & x) const {
      // rely on built-in error trapping
      try {
	Components<Scalar> cy(y);
	Components<Scalar> cx(x);
	for (int i=0; i<icvec.size(); i++) 
	  cx[i].copy(cy[icvec[i]]);
      }
      catch (RVLException &e) {
	e<<"\ncalled from LinInjectOp::apply\n";
	throw e;
      }
    }

    Operator<Scalar> * clone() const { return new LinInjectOp(*this); }

  public:

    LinInjectOp(Space<Scalar> const & _dom,
		Space<Scalar> const & _rng,
		int _icomp)
      : dom(_dom), rng(_rng), icvec(1,_icomp) {}

    LinInjectOp(Space<Scalar> const & _dom,
		  Space<Scalar> const & _rng,
		  std::vector<int> _icvec)
      : dom(_dom), rng(_rng), icvec(_icvec) {}
  
    
    LinInjectOp(LinInjectOp<Scalar> const & f)
      : dom(f.dom), rng(f.rng), icvec(f.icvec) {}

    ~LinInjectOp() {}

    Space<Scalar> const & getDomain() const { return dom; }
    Space<Scalar> const & getRange() const { return rng; }

    ostream & write(ostream & str) const {
      str<<"LinInjectOp object\n";
      str<<"components = \n";
      for (int i=0; i<icvec.size(); i++) str<<"  index = "<<i<<" value = "<<icvec[i]<<"\n";
      str<<"domain = \n";
      dom.write(str);
      str<<"range = \n";
      rng.write(str);
      return str;
    }
  };
  
}
    
#endif
