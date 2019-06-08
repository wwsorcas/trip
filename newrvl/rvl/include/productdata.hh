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

#ifndef __RVL_PDC
#define __RVL_PDC

#include "data.hh"
#include "product.hh"

namespace RVL {

  /** ProductDataContainers are DataContainers equipped with an
      indexing operator[], which returns a reference to a
      DataContainer when supplied with an in-range int index.
      Since ProductDataContainers act in this way like arrays of
      DataContainers, they also have a getSize() method. Both the
      indexing operator and getSize() are supplied by a mixin Product
      interface.

      This class is abstract, to permit a variety of schemes for
      storing and retrieving references to the component
      DataContainers. Evaluation of Function Objects is effectively
      equivalent to the simple loop algorithms of the Standard
      Library, but defined at a more abstract level, not requiring
      copy semantics of the items over which the iteration takes
      place.
  */

  class ProductDataContainer: public DataContainer, 
			      public Product<DataContainer> {

  public:

    ProductDataContainer() {}
    ProductDataContainer(const ProductDataContainer &) {}
    virtual ~ProductDataContainer() {}

    // can be overridden 
    virtual void eval( FunctionObject & f,
		       std::vector<std::shared_ptr<DataContainer const> > x) {
      try {
	size_t nx = x.size();
	std::vector<std::shared_ptr<ProductDataContainer const> > xp(nx);
	for (size_t i=0;i<nx;i++) {
	  if (!(xp[i]=dynamic_pointer_cast<ProductDataContainer const>(x[i]))) {
	    RVLException e;
	    e<<"Error: ProductDataContainer::eval\n";
	    e<<"argument "<<i<<" is not PDC\n";
	    throw e;
	  }
	  if (xp[i]->getSize() != getSize()) {
	    RVLException e;
	    e<<"Error: ProductDataContainer::eval \n";
	    e<<"input ProductDataContainer arg "<<i
	     <<" does not have same number\n";
	    e<<"of factors\n";
	    throw e;
	  }      
	}
	std::vector<std::shared_ptr<DataContainer const> > tmp(nx);
	size_t nc = this->getSize();
	for (size_t i=0;i<nc;i++) {
	  for (size_t j=0;j<nx;j++) {
	    tmp[j] = (*(xp[j]))[i];
	  }
	  //	  (*this)[i].eval(f,tmp);
	  ((*this)[i])->eval(f,tmp);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from ProductDataContainer::eval\n";
	throw e;
      }
      catch (...) {
	throw;
      }
    }

    virtual void eval(FunctionObjectConstEval & f,
		      std::vector<std::shared_ptr<DataContainer const> > x) const {
      try {
	size_t nx = x.size();
	std::vector<std::shared_ptr<ProductDataContainer const > > xp(nx);
	for (size_t i=0;i<nx;i++) {
	  if (!(xp[i]=dynamic_pointer_cast<ProductDataContainer const> (x[i]))) {
	    RVLException e;
	    e<<"Error: ProductDataContainer::eval\n";
	    e<<"argument "<<i<<" is not PDC\n";
	    throw e;
	  }
	  if (xp[i]->getSize() != getSize()) {
	    RVLException e;
	    e<<"Error: ProductDataContainer::eval \n";
	    e<<"input ProductDataContainer arg "<<i
	     <<" does not have same number\n";
	    e<<"of factors\n";
	    throw e;
	  }      
	}
	std::vector<std::shared_ptr<DataContainer const> > tmp(nx);
	size_t nc = this->getSize();
	for (size_t i=0;i<nc;i++) {
	  for (size_t j=0;j<nx;j++) {
	    //	    cerr<<"i="<<i<<"j="<<j<<"\n";
	    tmp[j] = (*(xp[j]))[i];
	  }
	  //	  cerr<<"eval(f,tmp)\n";
	  ((*this)[i])->eval(f,tmp);
	  //	  cerr<<"exit i="<<i<<"\n";
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from ProductDataContainer::eval (const)\n";
	throw e;
      }
      catch (...) {
	throw;
      }
    }

    /** report to ostream 
	make possible to override - WWS 2016.12.09 */
    virtual ostream & write(ostream & str) const {
      str<<"Product Data Container"<<endl;
      for (size_t i=0; i<getSize(); i++) {
	str<<"***factor "<<i<<":"<<endl;
	(*this)[i]->write(str);
      }
      return str;
    }

  };

  /** Standard static implementation of ProductDataContainer. Stores
      components by std::shared_ptr in a std::vector. The main
      constructor creates a StdProductDataContainer with size()
      uninitialized pointers, which must be initialized before use,
      from externally initialized shared pointers to appropriate
      objects, via the set() method. Both const and non-const
      operator[] versions throw an exception if the indexed shared_ptr
      is not initialized.
  */

  class StdProductDataContainer: public ProductDataContainer {

  private:
    
    mutable std::vector<std::shared_ptr<DataContainer> > data;
    
    /** Copy constructor. Deep copy - creates independent component
	objects, but does not copy data. Really a clone method. 
	Implementation removed due to dependence on DataContainer::clone(),
	which has been deprecated
    */
    StdProductDataContainer(const StdProductDataContainer & p);
    
  public:
    
    /** Default constructor yields size=n object */
    StdProductDataContainer(size_t n): data(n) {}

    ~StdProductDataContainer() {}

    size_t getSize() const { return data.size(); }

    void set(size_t i, std::shared_ptr<DataContainer> p) const {
      try {
	data.at(i) = p;
      }
      catch (out_of_range) {
	RVLException e;
	e<<"Error: StdProductDataContainer::set\n";
	e<<"  attempt to access component "<<i<<" of StdProductDC of size "<<data.size()<<"\n";
	throw e;
      }
    }
  
    std::shared_ptr<DataContainer> operator[](size_t i) {
      try {
	cerr<<"SPDC: op[] i="<<i<<"\n";
	return data.at(i);
      }
      catch (out_of_range) {
	RVLException e;
	e<<"Error: StdProductDataContainer::operator[]\n";
	e<<"attempt to access component "<<i<<" of StdProductDC of size "<<data.size()<<"\n";
	throw e;
      }
    }
  
    std::shared_ptr<DataContainer const> operator[](size_t i) const {
      try {
	if (data.at(i)) {
	  cerr<<"SPDC: op[] const i="<<i<<"\n";
	  return data.at(i);
	}
	else {
	  RVLException e;
	  e<<"Error: StdProductDataContainer::operator[], const\n";
	  e<<"  read access impossible - component "<<i<<" not initialized\n";
	  throw e;
	}
      }
      catch (out_of_range) {
	RVLException e;
	e<<"attempt to access component "<<i<<" of StdProductDC of size "<<this->getSize()<<"\n";
	throw e;
      }
    }  

  };

}

#endif










