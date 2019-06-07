/*************************************************************************

Copyright Rice University, 2004 - 2015.
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

#ifndef __RVL_VLOp
#define __RVL_VLOp

#include "space.hh"

namespace RVL {

  /**  This BFO implements the scalar version of the logistic 
       function 
       \f$ f(x) = x (1 + s^2 x^2)^{-1/2} + m \f$, where
       \f$ s=2/(f_{\rm max} - f_{\rm min})\f$ and
       \f$ m=(f_{\rm max} + f_{\rm min})/2\f$. This function has 
       the properties
       <ul>
       <li>\f$ f(x) \rightarrow f_{\rm max},\,\,x \rightarrow \infty\f$</li>
       <li>\f$ f(x) \rightarrow f_{\rm min},\,\,x \rightarrow -\infty\f$</li>
       <li>\f$ f(0) = m,\,\, f'(0)=1\f$</li>
       </ul>
       The function is applied to each component of the input LDC, value 
       written to the corresponding component of the output LDC. Input 
       and output must have the same length.
       
  */
  template<class Scalar>
  class RVLScalarLogistic: public BinaryLocalFunctionObject<Scalar> {
  private:
    Scalar s;
    Scalar m;
    RVLScalarLogistic(const RVLScalarLogistic<Scalar> &) {}
  public:
    RVLScalarLogistic(Scalar fmin=ScalarFieldTraits<Scalar>::Zero(),
		      Scalar fmax=ScalarFieldTraits<Scalar>::One()) {
      try {
	testRealOnly<Scalar>();
	if (!(fmin<fmax)) {
	  RVLException e;
	  e<<"ERROR: RVLScalarLogistic constructor\n";
	  e<<"fmin="<<fmin<<" not less than fmax="<<fmax<<"\n";
	  throw e;
	}
	s = 2.0 /(fmax-fmin);
	m = 0.5 *(fmax+fmin);
      }
      catch (RVLException & e) {
	e<<"\ncalled from RVLScalarLogistic constructor\n";
	throw e;
      }
    }
    ~RVLScalarLogistic() {}
  
    /** PRECONDITIONS:  Scalar field is real, fmin < fmax, 
        x.getSize() = y.getSize()
        <p>
	POSTCONDITION:  x[i] == f(y[i],fmin,fmax) for all i in range
    */
    using RVL::BinaryLocalEvaluation<Scalar>::operator();
    void operator()(LocalDataContainer<Scalar> & x,
		    LocalDataContainer<Scalar> const & y) {
      if (x.getSize() != y.getSize()) {
	RVLException e; e<<"ERROR: RVLScalarLogistic\n";
	e<<"input size = "<<x.getSize()<<" output size = "<<y.getSize()<<"\n";
	e<<"required to be same\n";
	throw e;
      }
      else {
	size_t n = x.getSize();
	Scalar * px = x.getData();
	Scalar const * py = y.getData();
	for (size_t i=0;i<n;i++) {
	  px[i]=m + py[i]/sqrt(1.0 + s*s*py[i]*py[i]);
	}
      }
    }
    string getName() const  { return "RVLScalarLogistic"; }
  };

  /**  This BFO implements the inverse of the scalar version
       of the logistic 
       function 
       \f$ f(x) = x (1 + s^2 x^2)^{-1/2} + m \f$, where
       \f$ s=2/(f_{\rm max} - f_{\rm min})\f$ and
       \f$ m=(f_{\rm max} + f_{\rm min})/2\f$. 
       The inverse is
       \f$ f^{-1}(y) = (y-m) (1 - s^2(y-m)^2))^{-1/2}\f$,
       well-defined when \f$f_{\rm min} < y < f_{\rm max}\f$.
       
  */
  template<class Scalar>
  class RVLScalarLogisticInverse: public BinaryLocalFunctionObject<Scalar> {
  private:
    Scalar s;
    Scalar m;
    RVLScalarLogisticInverse(const RVLScalarLogisticInverse<Scalar> &) {}
  public:
    RVLScalarLogisticInverse(Scalar fmin=ScalarFieldTraits<Scalar>::Zero(),
			     Scalar fmax=ScalarFieldTraits<Scalar>::One()) {
      try {
	testRealOnly<Scalar>();
	if (!(fmin<fmax)) {
	  RVLException e;
	  e<<"ERROR: RVLScalarLogisticInverse constructor\n";
	  e<<"fmin="<<fmin<<" not less than fmax="<<fmax<<"\n";
	  throw e;
	}
	s = 2.0 /(fmax-fmin);
	m = 0.5 *(fmax+fmin);
      }
      catch (RVLException & e) {
	e<<"\ncalled from RVLScalarLogisticInverse constructor\n";
	throw e;
      }
    }
    ~RVLScalarLogisticInverse() {}
  
    /** PRECONDITIONS:  Scalar field is real, fmin < y[i] < fmax, 
        x.getSize() = y.getSize()
        <p>
	POSTCONDITION:  x[i] == finv(y[i],fmin,fmax) for all i in range
    */
    using RVL::BinaryLocalEvaluation<Scalar>::operator();
    void operator()(LocalDataContainer<Scalar> & x,
		    LocalDataContainer<Scalar> const & y) {
      if (x.getSize() != y.getSize()) {
	RVLException e; e<<"ERROR: RVLScalarLogisticInverse\n";
	e<<"input size = "<<x.getSize()<<" output size = "<<y.getSize()<<"\n";
	e<<"required to be same\n";
	throw e;
      }
      else {
	size_t n = x.getSize();
	Scalar * px = x.getData();
	Scalar const * py = y.getData();
	for (size_t i=0;i<n;i++) {
	  if ((py[i] < m-1/s + numeric_limits<Scalar>::epsilon()) || 
	      (py[i] > m+1/s - numeric_limits<Scalar>::epsilon())) {
	    RVLException e;
	    e<<"ERROR: RVLScalarLogisticInverse\n";
	    e<<"input value = "<<py[i]<<" too close to \n";
	    e<<"fmin = "<<m-1/s<<" or fmax = "<<m+1/s<<"\n";
	    throw e;
	  }
	  px[i]=(py[i]-m)/sqrt(1.0 - s*s*(py[i]-m)*(py[i]-m));
	}
      }
    }
    string getName() const  { return "RVLScalarLogisticInverse"; }
  };

  /**  This TFO implements the scalar version of the logistic 
       function derivative
       \f$ df(x)dx = (1 + s^2 x^2)^{-3/2} dx\f$, where
       \f$ s=2/(f_{\rm max} - f_{\rm min})\f$
       <ul>
       <li>\f$ df(x) \rightarrow 0,\,\,x \rightarrow \pm \infty\f$</li>
       <li>\f$ df(x) > 0\f$</li>
       <li>\f$ df(0)=1\f$</li>
       </ul>
       The function is applied to each component of the input LDC, value 
       written to the corresponding component of the output LDC. Input 
       and output must have the same length.
       
  */
  template<class Scalar>
  class RVLScalarLogisticDeriv: public TernaryLocalFunctionObject<Scalar> {
  private:
    Scalar s;
    Scalar t;
    RVLScalarLogisticDeriv(const RVLScalarLogisticDeriv<Scalar> &) {}
  public:
    RVLScalarLogisticDeriv(Scalar fmin=ScalarFieldTraits<Scalar>::Zero(),
			   Scalar fmax=ScalarFieldTraits<Scalar>::One()) {
      try {
	testRealOnly<Scalar>();
	if (!(fmin<fmax)) {
	  RVLException e;
	  e<<"ERROR: RVLScalarLogisticDeriv constructor\n";
	  e<<"fmin="<<fmin<<" not less than fmax="<<fmax<<"\n";
	  throw e;
	}
	s = 2.0 /(fmax-fmin);
      }
      catch (RVLException & e) {
	e<<"\ncalled from RVLScalarLogisticDeriv constructor\n";
	throw e;
      }
    }
    ~RVLScalarLogisticDeriv() {}
  
    /** PRECONDITIONS:  Scalar field is real, fmin < fmax, 
        x.getSize() = y.getSize()
        <p>
	POSTCONDITION:  x[i] == df(y[i],fmin,fmax)dy[i] for all i in range
    */
    using RVL::TernaryLocalEvaluation<Scalar>::operator();
    void operator()(LocalDataContainer<Scalar> & x,
		    LocalDataContainer<Scalar> const & y, 
		    LocalDataContainer<Scalar> const & dy) {
      if ((x.getSize() != y.getSize()) ||
	  (x.getSize() != dy.getSize())) {
	RVLException e; e<<"ERROR: RVLScalarLogisticDeriv\n";
	e<<"input size = "<<y.getSize()<<"\n";
	e<<"input pert size = "<<dy.getSize()<<"\n";
	e<<"output size = "<<x.getSize()<<"\n";
	e<<"required to be same\n";
	throw e;
      }
      else {
	size_t n = x.getSize();
	Scalar * px = x.getData();
	Scalar const * py = y.getData();
	Scalar const * pdy = dy.getData();
	for (size_t i=0;i<n;i++) {
	  t=1.0/sqrt(1.0 + s*s*py[i]*py[i]);
	  px[i]=pdy[i]*t*t*t;
	}
      }
    }
    string getName() const  { return "RVLScalarLogisticDeriv"; }
  };

  /**  This QFO implements the vector version of the logistic 
       function 
       \f$ f(x) = x (1 + s^2 x^2)^{-1/2} + m \f$, where
       \f$ s=2/(f_{\rm max} - f_{\rm min})\f$ and
       \f$ m=(f_{\rm max} + f_{\rm min})/2\f$. This function has 
       the properties
       <ul>
       <li>\f$ f(x) \rightarrow f_{\rm max},\,\,x \rightarrow \infty\f$</li>
       <li>\f$ f(x) \rightarrow f_{\rm min},\,\,x \rightarrow -\infty\f$</li>
       <li>\f$ f(0) = m,\,\, f'(0)=1\f$</li>
       </ul>
       The function is applied to each component of the input LDC, value 
       written to the corresponding component of the output LDC. Input 
       and output must have the same length. Values of \f$f_{\rm min}\f$
       and \f$f_{\rm max}\f$ taken from two other LDCs of same length, 
       third and fourth arguments respectively.
       
  */
  template<class Scalar>
  class RVLVectorLogistic: public QuaternaryLocalFunctionObject<Scalar> {
  private:
    RVLVectorLogistic(const RVLVectorLogistic<Scalar> &) {}
  public:
    RVLVectorLogistic() {
      try {
	testRealOnly<Scalar>();
      }
      catch (RVLException & e) {
	e<<"\ncalled from RVLVectorLogistic constructor\n";
	throw e;
      }
    }
    ~RVLVectorLogistic() {}
  
    /** PRECONDITIONS:  Scalar field is real, fmin < fmax, 
        x.getSize() = y.getSize()
        <p>
	POSTCONDITION:  x[i] == f(y[i],fmin,fmax) for all i in range
    */
    using RVL::QuaternaryLocalEvaluation<Scalar>::operator();
    void operator()(LocalDataContainer<Scalar> & x,
		    LocalDataContainer<Scalar> const & y,
		    LocalDataContainer<Scalar> const & lb,
		    LocalDataContainer<Scalar> const & ub) {
      if ((x.getSize() != y.getSize()) ||
	  (x.getSize() != lb.getSize()) ||
	  (x.getSize() != ub.getSize())) {
	RVLException e; e<<"ERROR: RVLVectorLogistic\n";
	e<<"input size = "<<x.getSize()<<"\n";
	e<<"output size = "<<y.getSize()<<"\n";
	e<<"lb size = "<<lb.getSize()<<"\n";
	e<<"ub size = "<<ub.getSize()<<"\n";
	e<<"required to be same\n";
	throw e;
      }
      else {
	size_t n = x.getSize();
	Scalar * px = x.getData();
	Scalar const * py = y.getData();
	Scalar const * fmin = lb.getData();
	Scalar const * fmax = ub.getData();
	for (size_t i=0;i<n;i++) {
	  if (!(fmin[i]<fmax[i])) {
	    RVLException e;
	    e<<"ERROR: RVLVectorLogistic::operator()\n";
	    e<<"array index = "<<i<<"\n";
	    e<<"fmin="<<fmin[i]<<" not less than fmax="<<fmax[i]<<"\n";
	    throw e;
	  }
	  Scalar s = 2.0/(fmax[i]-fmin[i]);
	  Scalar m = 0.5 *(fmax[i]+fmin[i]);
	  px[i]=m + py[i]/sqrt(1.0 + s*s*py[i]*py[i]);
	}
      }
    }
    string getName() const  { return "RVLVectorLogistic"; }
  };

  /**  This O implements the vector version of the logistic 
       function derivative
       \f$ df(x)dx = (1 + s^2 x^2)^{-3/2} dx\f$, where
       \f$ s=2/(f_{\rm max} - f_{\rm min})\f$
       <ul>
       <li>\f$ df(x) \rightarrow 0,\,\,x \rightarrow \pm \infty\f$</li>
       <li>\f$ df(x) > 0\f$</li>
       <li>\f$ df(0)=1\f$</li>
       </ul>
       The function is applied to each component of the input LDC, value 
       written to the corresponding component of the output LDC. Input 
       and output must have the same length.
       
  */
  template<class Scalar>
  class RVLVectorLogisticDeriv: public LocalFunctionObject<Scalar> {
  public:
    RVLVectorLogisticDeriv() {}
    RVLVectorLogisticDeriv(const RVLVectorLogisticDeriv<Scalar> &) {}
    ~RVLVectorLogisticDeriv() {}
  
    using RVL::LocalEvaluation<Scalar>::operator();
    void operator()(LocalDataContainer<Scalar> & x,
		    std::vector<LocalDataContainer<Scalar> const *> & v) {
      if (v.size() != 4) {
	RVLException e;
	e<<"ERROR: RVLVectorLogisticDeriv::operator()\n";
	e<<"  input std::vector must have size = 4\n";
	e<<"  layout: v[0]=x, v[1]=dx, v[2]=xmin, v[3]=xmax\n";
	throw e;
      }
      if ((x.getSize() != v[0]->getSize()) ||
	  (x.getSize() != v[1]->getSize()) ||
	  (x.getSize() != v[2]->getSize()) ||
	  (x.getSize() != v[3]->getSize())) {	
	RVLException e; e<<"ERROR: RVLVectorLogisticDeriv\n";
	e<<"  incompatible array lengths\n";
	throw e;
      }
      else {
	size_t n = x.getSize();
	Scalar * px = x.getData();
	Scalar const * py = v[0]->getData();
	Scalar const * pdy = v[1]->getData();
	Scalar const * fm = v[2]->getData();
	Scalar const * fp = v[3]->getData();
	for (size_t i=0;i<n;i++) {
	  if (!(fm[i]<fp[i])) {
	    RVLException e;
	    e<<"ERROR: RVLVectorLogistic::operator()\n";
	    e<<"array index = "<<i<<"\n";
	    e<<"fmin="<<fm[i]<<" not less than fmax="<<fp[i]<<"\n";
	    throw e;
	  }
	  Scalar s = 2.0/(fp[i]-fm[i]);
	  Scalar t=1.0/sqrt(1.0 + s*s*py[i]*py[i]);
	  px[i]=pdy[i]*t*t*t;
	}
      }
    }
    string getName() const  { return "RVLVectorLogisticDeriv"; }
  };

  /**  This FO implements the inverse of the vector version
       of the logistic function 
       \f$ f(x) = x (1 + s^2 x^2)^{-1/2} + m \f$, where
       \f$ s=2/(f_{\rm max} - f_{\rm min})\f$ and
       \f$ m=(f_{\rm max} + f_{\rm min})/2\f$. 
       The inverse is
       \f$ f^{-1}(y) = (y-m) (1 - s^2(y-m)^2))^{-1/2}\f$,
       well-defined when \f$f_{\rm min} < y < f_{\rm max}\f$.
       
  */
  template<class Scalar>
  class RVLVectorLogisticInverse: public LocalFunctionObject<Scalar> {
  public:
    RVLVectorLogisticInverse() {}
    RVLVectorLogisticInverse(const RVLVectorLogisticInverse<Scalar> &) {}
    ~RVLVectorLogisticInverse() {}
  
    using RVL::LocalEvaluation<Scalar>::operator();
    void operator()(LocalDataContainer<Scalar> & x,
		    std::vector<LocalDataContainer<Scalar> const *> & v) {
      if (v.size() != 3) {
	RVLException e;
	e<<"ERROR: RVLVectorLogisticInverse::operator()\n";
	e<<"  input std::vector must have size = 4\n";
	e<<"  layout: v[0]=x, v[1]=xmin, v[2]=xmax\n";
	throw e;
      }
      if ((x.getSize() != v[0]->getSize()) ||
	  (x.getSize() != v[1]->getSize()) ||
	  (x.getSize() != v[2]->getSize())) {	
	RVLException e; e<<"ERROR: RVLVectorLogisticInverse\n";
	e<<"  incompatible array lengths\n";
	throw e;
      }
      else {
	size_t n = x.getSize();
	Scalar * px = x.getData();
	Scalar const * py = v[0]->getData();	
	Scalar const * fm = v[1]->getData();
	Scalar const * fp = v[2]->getData();
	//	cerr<<"****** beg RVLVectorLogisticInverse::operator()"<<endl;
	for (size_t i=0;i<n;i++) {
	  Scalar s = 2.0/(fp[i]-fm[i]);
	  Scalar m = 0.5*(fp[i]+fm[i]);
	  if (s*s*(py[i]-m)*(py[i]-m) > 1.0 - numeric_limits<Scalar>::epsilon()) {
	    //	  if ((py[i] < (m-1/s) + 
	    //	      (py[i] > m+1/s - numeric_limits<Scalar>::epsilon())) {
	    RVLException e;
	    e<<"ERROR: RVLVectorLogisticInverse\n";
	    e<<"  input value = "<<py[i]<<" too close to \n";
	    e<<"  fmin = "<<m-1/s<<" or fmax = "<<m+1/s<<"\n";
	    e<<"  or exceeds these bounds - cannot invert logistic map\n";
	    e<<"  supplied function bounds are "<<fm[i]<<" and "<<fp[i]<<"\n";
	    throw e;
	  }
	  px[i]=(py[i]-m)/sqrt(1.0 - s*s*(py[i]-m)*(py[i]-m));
	  //	  cerr<<"i="<<i<<" input="<<py[i]<<" lb="<<fm[i]<<" ub="<<fp[i]<<" s="<<s<<" m="<<m<<" output="<<px[i]<<endl;
	}
	//	cerr<<"****** end RVLVectorLogisticInverse::operator()"<<endl;	
      }
    }
    string getName() const  { return "RVLVectorLogisticInverse"; }
  };

  template<typename Scalar>
  class VectorLogisticOp: public Operator<Scalar> {

  private:

    // these should be shared_ptr objects - to do
    Vector<Scalar> const & lb;
    Vector<Scalar> const & ub;

    VectorLogisticOp();
    
  protected:

    void apply(Vector<float> const & x,
	       Vector<float> & y) const {
      try {
	RVLVectorLogistic<Scalar> f;
	y.eval(f,x,lb,ub);
      }
      catch (RVLException & e) {
	e<<"\ncalled from VectorLogisticOp::apply\n";
	throw e;
      }
    }
      
    void applyDeriv(Vector<float> const & x,
		    Vector<float> const & dx,
		    Vector<float> & dy) const {
      try {
	RVLVectorLogisticDeriv<Scalar> f;
	std::vector<RVL::Vector<Scalar> const *> v;
	v.push_back(&x);
	v.push_back(&dx);
	v.push_back(&lb);
	v.push_back(&ub);
	dy.eval(f,v);
      }
      catch (RVLException & e) {
	e<<"\ncalled from VectorLogisticOp::applyDeriv\n";
	throw e;
      }
    }
    
    void applyAdjDeriv(Vector<float> const & x,
		       Vector<float> const & dy,
		       Vector<float> & dx) const {
      try {
	RVLVectorLogisticDeriv<Scalar> f;
	std::vector<RVL::Vector<Scalar> const *> v;
	v.push_back(&x);
	v.push_back(&dy);
	v.push_back(&lb);
	v.push_back(&ub);
	dx.eval(f,v);
      }
      catch (RVLException & e) {
	e<<"\ncalled from VectorLogisticOp::applyDeriv\n";
	throw e;
      }
    }

    /* note: not implemented yet

    void applyDeriv2(const Vector<float> & x,
		     const Vector<float> & dx1,
		     const Vector<float> & dx2,
		     Vector<float> & dy) const {
      RVLException e;
      e<<"ERROR: VectorLogisticOp::applyDeriv2 not defined\n";
      e<<"  complain to management\n";
      throw e;
    }
    void applyAdjDeriv2(const Vector<float> & x,
			const Vector<float> & dy,
			const Vector<float> & dx2,
			Vector<float> & dx1) const {
      RVLException e;
      e<<"ERROR: VectorLogisticOp::applyAdjDeriv2 not defined\n";
      e<<"  complain to management\n";
      throw e;
    }
    */
    
    Operator<float> * clone() const { return new VectorLogisticOp<Scalar>(*this); }

  public:

    VectorLogisticOp(Vector<Scalar> const & _lb,
		     Vector<Scalar> const & _ub)
      : lb(_lb), ub(_ub) {}

    VectorLogisticOp(VectorLogisticOp<Scalar> const & op) 
      : lb(op.lb), ub(op.ub) {}

    ~VectorLogisticOp() {}

    Space<Scalar> const & getDomain() const { return lb.getSpace(); }
    Space<Scalar> const & getRange() const { return lb.getSpace(); }

    ostream & write(ostream & str) const {
      str<<"Vector Logistic Operator implementing field box constraints\n";
      str<<"*** lower bound vector: \n";
      lb.write(str);
      str<<"*** upper bound vector: \n";
      ub.write(str);
      return str;
    }
  };

}

#endif
