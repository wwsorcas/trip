/*************************************************************************

Copyright Rice University, 2004, 2005, 2006.
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

#ifndef __RVL_PVEC
#define __RVL_PVEC

#include "space.hh"
#include "productdata.hh"

namespace RVL {

  /** Abstract base class expressing behaviour of Cartesian products of
      vector spaces. Uses the ProductDataContainer class to define data
      component access. Provides basic data container methods (getSize,
      operator[]) for a Cartesian product. All other methods of Space
      are implemented in terms of these. 

      A default implementation is provided for inner product, expressing
      the canonical inner product in a product Hilbert space. This method
      is virtual and may be overridden in subclasses.  */

  template<class Scalar>
  class ProductSpace: public Space<Scalar>,
		      public ROProduct< Space<Scalar> > {

  public:

    virtual ~ProductSpace() {}

    virtual bool isCompatible(DataContainer const & dc) const {
      try {
	ProductDataContainer const & pdc =
	  dynamic_cast<ProductDataContainer const &>(dc);
	if (pdc.getSize() != this->getSize()) return 0;
	for (size_t i=0;i<this->getSize();i++) {
	  if (!((*this)[i]->isCompatible(*(pdc[i])))) return 0;
	}
	return true;
      }
      catch (bad_cast) {
	return false;
      }
    }

    virtual bool operator ==(const Space<Scalar> & sp) const {
      if (this == &sp) return true;
      try {
	const ProductSpace<Scalar> & psp 
	  = dynamic_cast<const ProductSpace<Scalar> &>(sp);
	if (psp.getSize() != this->getSize()) return false;
	for (size_t i=0;i<this->getSize();i++) {
	  if (*(psp[i]) != *((*this)[i])) return false;
	}
	return true;
      }
      catch (bad_cast) {
	return false;
      }
    }

    /** inner product - may be overridden by supplying non-block-diag
	Gram matrix. */
    virtual Scalar inner(std::shared_ptr<DataContainer const> x, 
			 std::shared_ptr<DataContainer const> y) const {
      Scalar res = ScalarFieldTraits<Scalar>::Zero();
      try {
	std::shared_ptr<ProductDataContainer const> xp = 
	  dynamic_pointer_cast<ProductDataContainer const>(x);
	std::shared_ptr<ProductDataContainer const> yp = 
	  dynamic_pointer_cast<ProductDataContainer const>(y);
	if (xp->getSize() != this->getSize() ||
	    yp->getSize() != this->getSize()) {
	  RVLException e; e<<"Error: ProductSpace::inner\n";
	  e<<"input data containers not of same size as space\n";
	  throw e;
	}
	for (size_t i=0;i<this->getSize();i++) {
	  try {
	    res += (*this)[i]->inner((*xp)[i],(*yp)[i]);
	  }
	  catch (RVLException & e) {
	    e<<"\ncalled from ProductSpace::inner\n";
	    throw e;
	  }
	}
	return res;
      }
      catch (bad_cast) {
	RVLException e; e<<"Error: ProductSpace::inner\n";
	e<<"inputs not (both) product vectors\n";
	throw e;
      }
    }

    // zero vector
    virtual void zero(std::shared_ptr<DataContainer> x) const {
      try {
	std::shared_ptr<ProductDataContainer> xp = 
	  dynamic_pointer_cast<ProductDataContainer>(x);
	if (xp->getSize() != this->getSize()) {
	  RVLException e; e<<"Error: ProductSpace::zero\n";
	  e<<"input data container not of same size as space\n";
	  throw e;
	}
	for (size_t i=0;i<this->getSize();i++) {
	  (*this)[i]->zero((*xp)[i]);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from ProductSpace::zero\n";
	throw e;
      }
      catch (bad_cast) {
	RVLException e; e<<"Error: ProductSpace::zero\n";
	e<<"input not product vector\n";
	throw e;
      }
    }

    // linear combination
    virtual void linComb(Scalar a, std::shared_ptr<DataContainer const> x,
			 Scalar b, std::shared_ptr<DataContainer> y) const {
      try {
	std::shared_ptr<ProductDataContainer const> xp = 
	  dynamic_pointer_cast<ProductDataContainer const>(x);
	std::shared_ptr<ProductDataContainer> yp = 
	  dynamic_pointer_cast<ProductDataContainer>(y);
	if (xp->getSize() != this->getSize() || 
	    yp->getSize() != this->getSize()) {
	  RVLException e; e<<"Error: ProductSpace::linComb\n";
	  e<<"input data containers not of same size as space\n";
	  e<<"**** FIRST INPUT PDC:\n";
	  x->write(e);
	  e<<"**** SECOND INPUT PDC:\n";
	  y->write(e);
	  e<<"**** SPACE:\n";
	  Space<Scalar>::write(e);
	  throw e;
	}
	for (size_t i=0;i<this->getSize();i++) {
	  (*this)[i]->linComb(a,(*xp)[i],b,(*yp)[i]);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from ProductSpace::linComb\n";
	throw e;
      }
      catch (bad_cast) {
	RVLException e; e<<"Error: ProductSpace::linComb\n";
	e<<"inputs not all product vectors\n";
	throw e;
      }
    }

    /** Write method - may be overwritten if useful. */
    virtual ostream & write(ostream & str) const {
      str<<"Product Space"<<endl;
      str<<"size = "<<this->getSize()<<"\n";
      for (size_t i=0; i<this->getSize(); i++) {
	str<<"***factor "<<i<<":"<<endl;
	(*this)[i]->write(str);
      }
      return str;
    }
  };

  /** Standard static construction of product spaces, via an (STL) vector of
      dynamically allocated Spaces. 
  */

  template<class Scalar>
  class StdProductSpace: public ProductSpace<Scalar> {

  private: 

    /** primary data vector */
    mutable std::vector<std::shared_ptr<Space<Scalar> > > s;

    /** Default construction - disallowed */
    StdProductSpace();

  protected:

    /** checks that factor Spaces are initialized */
    bool init() const {
      bool res = true;
      for (size_t i=0; i< s.size(); i++) res = res && (s[i]);
      return res;
    }

  public:

    /** virtual copy construction */
    std::shared_ptr<Space<Scalar> > clone() const {
      return make_shared<StdProductSpace<Scalar> >(*this);
    }

    /** assignment of ith factor */
    void set(size_t i, std::shared_ptr<Space<Scalar> > p) const {
      try {
	if (!(s.at(i))) s.at(i)=p;
	else {
	  RVLException e;
	  e<<"Error: StdProductSpace::set\n";
	  e<<"  attempt to set factor "<<i<<" which is already set\n";
	  throw e;
	}
      }
      catch (out_of_range) {
	RVLException e;
	e<<"Error: StdProductSpace::set\n";
	e<<"  attempt to set factor "<<i<<" index out of range [0,"
	 <<this->getSize()<<"]\n";
	throw e;	
      }
    }
    
    /** Copy construction. */
    StdProductSpace(const StdProductSpace<Scalar> & sp):
      s(sp.s.size(),NULL) {
      try {
	if (sp.init()) {
	  for (size_t i=0;i<s.size();i++) {
	    if ((sp.s)[i]) 
	      s[i]=((sp.s)[i])->clone();
	  }
	}
	else {
	  RVLException e;
	  e<<"Error: StdProductSpace copy construction\n";
	  e<<"  prototype not initialized\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from StdProductSpace copy constructor\n";
	throw e;
      }
    }
    
    /** construction for nfac factors, without initialization */
    StdProductSpace(size_t nfac): s(nfac,NULL) {
      //cerr<<"StdProductSpace constructed on "<<nfac<<" factors, not initialized\n";
    }

    ~StdProductSpace() {}

    /** implements virtual DataContainer build method via
	StdProductDataContainer class. Made virtual so that
	subclasses can use specialized subtypes of StdProdDC */
    virtual std::shared_ptr<DataContainer> buildDataContainer() const {
      try {
	if (init()) {
	  std::shared_ptr<StdProductDataContainer> d = 
	    make_shared<StdProductDataContainer>(this->getSize());
	  for (size_t i=0; i<this->getSize(); i++) {
	    cerr<<"SPSP: assign dc "<<i<<"\n";
	    d->set(i,s[i]->buildDataContainer());
	  }
	  return d;
	}
	else {
	  RVLException e;
	  e<<"ERROR: StdProductSpace::buildDataContainer\n";
	  e<<"  called on uninitialized StdProductSpace object\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from StdProductSpace::buildDataContainer\n";
	throw e;
      }
    }

    size_t getSize() const { return s.size(); }
    std::shared_ptr<Space<Scalar> const> operator[](size_t i) const {
      if (init()) return s.at(i);
      else {
	RVLException e;
	e<<"ERROR: StdProductSpace::operator[]\n";
	e<<"  called on uninitialized StdProductSpace object\n";
	throw e;
      }	
    }
  };

  /** Provides indexed access to components (factors) of a Vector in 
      ProductSpace.  
      
      Set up to work properly for non-product vectors also - these have
      one component!

      Note however that the singleton component is not the same as the
      parent vector - these are different objects, referring to the
      same data.
  */
  template<class Scalar>
  class Components: public Product< Vector<Scalar> > {

  private:

    size_t size;
    size_t & ver;

    std::shared_ptr<DataContainer> dc;
    std::shared_ptr<ProductDataContainer> pdc;
    const Space<Scalar> & sp;
    const ProductSpace<Scalar> * psp;

    Components();
    Components(const Components<Scalar> &);

  public:

    /* Constructor. Uses protected Vector constructor to build
       (STL) vector of component Vector(s). */
    Components(const Vector<Scalar> & v)
      : size(1), ver(v.getVersionRef()),
	dc(v.getDataContainer()), sp(v.getSpace()),
	psp(NULL) {
      try {
	if (pdc = 
	    dynamic_pointer_cast<ProductDataContainer>(dc)) {
	  if (psp =
	      dynamic_cast<ProductSpace<Scalar> const *>(&sp)) {
	    size = psp->getSize();
	  }
	  else {
	    RVLException e;
	    e<<"Error: Components constructor\n";
	    e<<"  weird one: dc is pdc, but sp is not psp\n";
	    e<<"  input vector:\n";
	    v.write(e);
	  }
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from Components constructor\n";
	throw e;
      }
    }

    /** return number of components */
    size_t getSize() const { return size; }
    /** return ith component; throw exception if index is out of bounds. */
    std::shared_ptr<Vector<Scalar> > operator[](size_t i) {
      if (i>size-1) {
	RVLException e; e<<"Error: Components::operator[]\n"; 
	e<<"index "<<i<<" out of range [0, "<<size-1<<"]\n";
	throw e;
      }
      if (pdc) {
	//	return make_shared<Vector<Scalar> >(*((*psp)[i]),(*pdc)[i],ver);
	return std::shared_ptr<Vector<Scalar> >(new Vector<Scalar>(*((*psp)[i]),(*pdc)[i],ver));
      }
      else {
	return std::shared_ptr<Vector<Scalar> >(new Vector<Scalar>(sp,dc,ver));
	//	return make_shared<Vector<Scalar> >(sp,dc,ver);
      }
    }

    std::shared_ptr<Vector<Scalar> const> operator[](size_t i) const {
      if (i>size-1) {
	RVLException e; e<<"Error: Components::operator[]\n"; 
	e<<"index "<<i<<" out of range [0, "<<size-1<<"]\n";
	throw e;
      }
      if (pdc) {
	return std::shared_ptr<Vector<Scalar> >(new Vector<Scalar>(*((*psp)[i]),(*pdc)[i],ver));	
	//	return make_shared<Vector<Scalar> const>(*((*psp)[i]),(*pdc)[i],ver);
      }
      else {
	return std::shared_ptr<Vector<Scalar> >(new Vector<Scalar>(sp,dc,ver));
	//	return make_shared<Vector<Scalar> const>(sp,dc,ver);
      }
    }

    ostream & write(ostream & str) const {
      str<<"Components of a Vector"<<"\n";
      str<<"  number = "<<size<<"\n";
      str<<"  parent space = "<<"\n";
      sp.write(str);
      return str;
    }
  };

}

#endif
