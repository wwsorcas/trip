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

#ifndef __RVL_FCNS
#define __RVL_FCNS

#include "utility.hh"
#include "local.hh"

namespace RVL {

  /**  This BFO copies the second LDC onto the first
  */
  template<class Scalar>
  class RVLCopy: public BinaryLocalFunctionObject<Scalar> {
  private:
    RVLCopy(const RVLCopy<Scalar> &) {}
  public:
    RVLCopy() {}
    ~RVLCopy() {}
  
    /** PRECONDITION:  x.getSize() <= y.getSize()
	POSTCONDITION:  x[i] == y[i] for all i in range
    */
    using RVL::BinaryLocalEvaluation<Scalar>::operator();
    void operator()(LocalDataContainer<Scalar> & x,
		    LocalDataContainer<Scalar> const & y) {
      if (x.getSize() > y.getSize()) {
	RVLException e; e<<"Error: RVLCopy\n";
	e<<"input shorter than output - copy does not cover\n";
	throw e;
      }
      else {
	size_t n = x.getSize();
	Scalar * px = x.getData();
	Scalar const * py = y.getData();
	for (size_t i=0;i<n;i++) {
	  px[i]=py[i];
	}
      }
    }
    string getName() const  { return "RVLCopy"; }
  };

  template<>
  class RVLCopy<float>: public BinaryLocalFunctionObject<float> {
  private:
    RVLCopy(const RVLCopy<float> &) {}
  public:
    RVLCopy() {}
    ~RVLCopy() {}
  
    /** PRECONDITION:  x.getSize() <= y.getSize()
	POSTCONDITION:  x[i] == y[i] for all i in range
    */
    using RVL::BinaryLocalEvaluation<float>::operator();
    void operator()(LocalDataContainer<float> & x,
		    LocalDataContainer<float> const & y) {
      if (x.getSize() > y.getSize()) {
	RVLException e; e<<"Error: RVLCopy\n";
	e<<"input shorter than output - copy does not cover\n";
	throw e;
      }
      else {
	float * px = x.getData();
	float const * py = y.getData();
	memcpy(px, py, x.getSize()*sizeof(float));
      }
    }
    string getName() const  { return "RVLCopy<float>"; }
  };

  template<>
  class RVLCopy<double>: public BinaryLocalFunctionObject<double> {
  private:
    RVLCopy(const RVLCopy<double> &) {}
  public:
    RVLCopy() {}
    ~RVLCopy() {}
  
    /** PRECONDITION:  x.getSize() <= y.getSize()
	POSTCONDITION:  x[i] == y[i] for all i in range
    */
    using RVL::BinaryLocalEvaluation<double>::operator();
    void operator()(LocalDataContainer<double> & x,
		    LocalDataContainer<double> const & y) {
      if (x.getSize() > y.getSize()) {
	RVLException e; e<<"Error: RVLCopy\n";
	e<<"input shorter than output - copy does not cover\n";
	throw e;
      }
      else {
	double * px = x.getData();
	double const * py = y.getData();
	memcpy(px, py, x.getSize()*sizeof(double));
      }
    }
    string getName() const  { return "RVLCopy<double>"; }
  };

  /**  This UFO scales an LDC
  */
  template<class Scalar>
  class RVLScale: public UnaryLocalFunctionObject<Scalar> {
  private:
    Scalar c;
    RVLScale();
    RVLScale(const RVLScale<Scalar> &);
  public:
    RVLScale(Scalar _c): c(_c) {}
    ~RVLScale() {}

    /** PRECONDITION:  x is an LDC
	POSTCONDITION:  x[i] == c * x[i] for all i in range
    */
    using RVL::UnaryLocalEvaluation<Scalar>::operator();
    void operator()(LocalDataContainer<Scalar> & x) {
      size_t n = x.getSize();
      Scalar * px = x.getData();
      for (size_t i=0;i<n;i++) {
	px[i]=c*(px[i]);
      }
    }
    string getName() const  { return "RVLScale"; }
  };

  /**  This UFOSR and Accumulation finds the maximum element of an LDC.
       Works of course only for Scalar types that have a max function.
  */
  template<class Scalar>
  class RVLMax
    : public UnaryLocalFunctionObjectScalarRedn<Scalar,Scalar> {
  public:
    RVLMax()
      : UnaryLocalFunctionObjectScalarRedn<Scalar,Scalar>(-numeric_limits<Scalar>::max()) {}
    RVLMax(Scalar _res) 
      : UnaryLocalFunctionObjectScalarRedn<Scalar,Scalar>(_res) {}
    RVLMax(const RVLMax<Scalar> & m)
      : UnaryLocalFunctionObjectScalarRedn<Scalar,Scalar>(m) {}
    ~RVLMax() {}

    // override base class
    void setValue() { ScalarRedn<Scalar>::setValue(-numeric_limits<Scalar>::max()); }
    using RVL::UnaryLocalConstEval<Scalar>::operator();
    void operator() (LocalDataContainer<Scalar> const & x) {
      Scalar maxr = ScalarRedn<Scalar>::getValue();
      size_t n=x.getSize();
      Scalar const * px = x.getData();
      for (size_t i=0;i<n;i++) maxr = max<Scalar>(maxr,px[i]);
      ScalarRedn<Scalar>::setValue(maxr);
    }

    string getName() const  { return "RVLMax"; }
    
  };
  
  /**  This UFOSR and Accumulation finds the minimum element of an LDC
       Works of course only for Scalar types that have a min function..
   */
  template<class Scalar>
  class RVLMin: public UnaryLocalFunctionObjectScalarRedn<Scalar,Scalar> {
  public:
    RVLMin()
      : UnaryLocalFunctionObjectScalarRedn<Scalar,Scalar>(numeric_limits<Scalar>::max()) {}
    RVLMin(Scalar _res)
      : UnaryLocalFunctionObjectScalarRedn<Scalar,Scalar>(_res) { }
    RVLMin(const RVLMin<Scalar> & m)
      : UnaryLocalFunctionObjectScalarRedn<Scalar,Scalar>(m) {}
    ~RVLMin() {}
  
    void setValue() { ScalarRedn<Scalar>::setValue(numeric_limits<Scalar>::max()); }
    using RVL::UnaryLocalConstEval<Scalar>::operator();
    void operator() (LocalDataContainer<Scalar> const & x) {
      Scalar minr = ScalarRedn<Scalar>::getValue();
      size_t n=x.getSize();
      Scalar const * px = x.getData();
      for (size_t i=0;i<n;i++) { minr=min<Scalar>(minr,px[i]); }
      ScalarRedn<Scalar>::setValue(minr);
    }

    string getName() const  { return "RVLMin"; }
  };

  /** This BFOSR and Accumulation computes the inner product of two LDCs.
  */
  template<class Scalar>
  class RVLL2innerProd: public BinaryLocalFunctionObjectScalarRedn<Scalar,Scalar> {

  private:
    mutable Scalar scale;
  public:
    /** \param _scale a scale factor applied to the inner product.
	Defaults to 1 \param _init uses 0 as the default initial
	value, but can be set otherwise.
     */
    RVLL2innerProd(Scalar _scale = ScalarFieldTraits<Scalar>::One(), 
		   Scalar _init = ScalarFieldTraits<Scalar>::Zero())
      :  BinaryLocalFunctionObjectScalarRedn<Scalar,Scalar>(_init),
	 scale(_scale) {}
    RVLL2innerProd(const RVLL2innerProd<Scalar> & ipfo) 
      :  BinaryLocalFunctionObjectScalarRedn<Scalar,Scalar>(ipfo),
	 scale(ipfo.scale) {}
    ~RVLL2innerProd() {}

    void setValue() { ScalarRedn<Scalar>::setValue(ScalarFieldTraits<Scalar>::Zero()); }

    using RVL::BinaryLocalConstEval<Scalar>::operator();
    void operator() (LocalDataContainer<Scalar> const & v, 
		     LocalDataContainer<Scalar> const & w) {
      try {
	size_t n=v.getSize();
	if (n != w.getSize()) {
	  RVLException e; e<<"Error: RVLL2innerProd::operator()\n";
	  e<<"operands do not have same dimension\n";
	  e<<"\noperand 1:\n";
	  v.write(e);
	  e<<"\noperand 2:\n";
	  w.write(e);
	  throw e;
	}
	Scalar const * pv = v.getData();
	Scalar const * pw = w.getData();
	Scalar raw = ScalarFieldTraits<Scalar>::Zero();
	Scalar ip = this->getValue();
	for (size_t i=0;i<n;i++) {
	  raw += pv[i]*pw[i];
	}
	ip += scale*raw;
	ScalarRedn<Scalar>::setValue(ip);
      }
      catch (RVLException & e) {
	e<<"\ncalled from RVLL2innerProduct::operator()\n";
	throw e;
      }
    }

    /** added to enable faithful copy in LinAlg package copy constructors */
    Scalar getScale() const { return scale; }

    /** added to separate instantiation from initialization. */
    void setScale(Scalar newscale) { scale = newscale; }

    string getName() const  { return "basic rvl L2innerProd"; }

  };

  /** Complex specialization of L2 Inner Product. Note: scale is REAL.
  */
  template<class Scalar>
  class RVLL2innerProd<complex<Scalar> >: 
    public BinaryLocalFunctionObjectScalarRedn<complex<Scalar>, complex<Scalar> > {
  private:
    Scalar scale;
  public:
    /** \param _scale a scale factor applied to the inner product.
	Defaults to 1 \param _init uses 0 as the default initial
	value, but can be set otherwise.
     */
    RVLL2innerProd(Scalar _scale = ScalarFieldTraits<Scalar>::One(), 
		   complex<Scalar> _init = complex<Scalar>(ScalarFieldTraits<Scalar>::Zero()))
      : BinaryLocalFunctionObjectScalarRedn<complex<Scalar>, complex<Scalar> > (_init),
	 scale(_scale) {}
    RVLL2innerProd(const RVLL2innerProd<complex<Scalar> > & ipfo)
      : BinaryLocalFunctionObjectScalarRedn<complex<Scalar>, complex<Scalar> > (ipfo), 
	 scale(ipfo.scale) {}
    ~RVLL2innerProd() {}

    void setValue() { ScalarRedn<complex<Scalar> >::setValue(complex<Scalar>(ScalarFieldTraits<Scalar>::Zero())); }

    using RVL::BinaryLocalConstEval< complex<Scalar> >::operator();
    void operator() (LocalDataContainer<complex<Scalar> > const & v, 
		     LocalDataContainer<complex<Scalar> > const & w) {
      try {
	size_t n=v.getSize();
	if (n != w.getSize()) {
	  RVLException e; e<<"Error: RVLL2innerProd<complex>::operator()\n";
	  e<<"operands do not have same dimension\n";
	  e<<"\noperand 1:\n";
	  v.write(e);
	  e<<"\noperand 2:\n";
	  w.write(e);
	  throw e;
	}
	complex<Scalar> const * pv = v.getData();
	complex<Scalar> const * pw = w.getData();
	complex<Scalar> raw = complex<Scalar>(ScalarFieldTraits<Scalar>::Zero());
	complex<Scalar> ip = BinaryLocalFunctionObjectScalarRedn<complex<Scalar>, complex<Scalar> >::getValue();
	for (size_t i=0;i<n;i++) {
	  raw += pv[i]*conj(pw[i]);
	}
	ip += scale*raw;
	ScalarRedn<complex<Scalar> >::setValue(ip);
      }
      catch (RVLException & e) {
	e<<"\ncalled from RVLL2innerProduct::operator()\n";
	throw e;
      }
    }

    /** added to enable faithful copy in LinAlg package copy constructors */
    Scalar getScale() const { return scale; }

    /** added to separate instantiation from initialization. */
    void setScale(Scalar newscale) { scale = newscale; }

    string getName() const  { return "basic rvl L2innerProd <complex>"; }

  };

  /** This BFO does v += w
   */
  template<class Scalar>
  class RVLAddAccumulate: public BinaryLocalFunctionObject<Scalar> {
    
  public:
    RVLAddAccumulate() {}
    RVLAddAccumulate(const RVLAddAccumulate<Scalar> &) {}
    virtual ~RVLAddAccumulate() {}

    using RVL::BinaryLocalEvaluation<Scalar>::operator();
    void operator() (LocalDataContainer<Scalar> & v,
		     LocalDataContainer<Scalar> const & w) {
      size_t n=v.getSize();
      if (n != w.getSize()) {
	RVLException e;
	e<<"Error: RVLAddAccumulate::operator()\n";
	e<<"inputs of different sizes\n";
	e<<"first input:\n";
	v.write(e);
	e<<"second input:\n";
	w.write(e);
	throw e;
      }
      Scalar * pv = v.getData();
      Scalar const * pw = w.getData();
      for (size_t i=0;i<n;i++) {
	pv[i] += pw[i];
      }
    }
    string getName() const  { string sname =  "RVLAddAccumulate"; return sname; }
  };

  /**  This UFO sets all elements of a LDC to a constant value.

  */
  template<class Scalar>
  class RVLAssignConst: public UnaryLocalFunctionObject<Scalar> {
  private:
    Scalar c;
    RVLAssignConst();
    RVLAssignConst(const RVLAssignConst<Scalar> &);
  public:
    RVLAssignConst(Scalar _c): c(_c) {}
    ~RVLAssignConst() {}

    using RVL::UnaryLocalEvaluation<Scalar>::operator();
    void operator() (LocalDataContainer<Scalar> & v) {
      size_t n = v.getSize();
      Scalar * pv = v.getData();
      for (size_t i=0;i<n;i++) {
	pv[i]=c;
      }
    }

    string getName() const { return "RVLAssignConst"; }

  };

  /**  Assigns random values in [a, b] to each entry in the LDC
       where the defaults are a = 0, b = 1.
   */
  template<class Scalar>
  class RVLRandomize: public UnaryLocalFunctionObject<Scalar> {
  private:
    Scalar a;
    Scalar w;

    RVLRandomize(const RVLRandomize<Scalar> &) {}
    
  public:
    RVLRandomize(): a(0), w(1) {}
    RVLRandomize( long seed, Scalar _a = 0, Scalar _b = 1 )
      : a(_a), w(_b-_a) { PlantSeeds(seed); }
    ~RVLRandomize() {}
  
    virtual bool readsData(size_t i=0) { return false; }

    using RVL::UnaryLocalEvaluation<Scalar>::operator();
    void operator() (LocalDataContainer<Scalar> & v) {
      size_t n = v.getSize();
      Scalar * pv = v.getData();
      if( a == ScalarFieldTraits<Scalar>::Zero()) {
	if( w == ScalarFieldTraits<Scalar>::One()) {
	  for (size_t i=0;i<n;i++) {
	    //	    pv[i]=rand()/(RAND_MAX+ScalarFieldTraits<Scalar>::One());
	    pv[i]=Random();
	  }
	}
	else {
	  for (size_t i=0;i<n;i++) {
	    //	    pv[i]=w*rand()/(RAND_MAX+ScalarFieldTraits<Scalar>::One());
	    pv[i]=w*Random();
	  }
	}
      } else if (w == 1) {
	for (size_t i=0;i<n;i++) {
	  //	  pv[i]=rand()/(RAND_MAX+ScalarFieldTraits<Scalar>::One()) + a;
	  pv[i]=Random() + a;
	}
      } else {
	for (size_t i=0;i<n;i++) {
	  //	  pv[i]=w*rand()/(RAND_MAX+ScalarFieldTraits<Scalar>::One()) + a;
	  pv[i]=w*Random() + a;
	}
      }
    }
    
    string getName() const { return "RVLRandomize"; }
    
  };

  /**  Assigns random values in [a, b] to each entry in the LDC
       where the defaults are a = 0, b = 1.
   */
  template<class Scalar>
  class RVLRandomize<complex<Scalar> >: public UnaryLocalFunctionObject<complex<Scalar> >{
  private:
    Scalar a;
    Scalar w;

    RVLRandomize(const RVLRandomize<complex<Scalar> >&) {}
    
  public:
    RVLRandomize()
      : a(ScalarFieldTraits<Scalar>::Zero()), w(ScalarFieldTraits<Scalar>::One()) {}
    RVLRandomize( //unsigned int seed, 
		 long seed, 
		 Scalar _a = ScalarFieldTraits<Scalar>::Zero(), 
		 Scalar _b = ScalarFieldTraits<Scalar>::One() )
      : a(_a), w(_b-_a) { //srand(seed); }
      PlantSeeds(seed); }
    ~RVLRandomize() {}
  
    using RVL::UnaryLocalEvaluation< complex<Scalar> >::operator();
    void operator() (LocalDataContainer<complex<Scalar> > & v) {
      size_t n = v.getSize();
      complex<Scalar> * pv = v.getData();
      Scalar rex;
      Scalar imx;
      if( a == ScalarFieldTraits<Scalar>::Zero()) {
	if( w == ScalarFieldTraits<Scalar>::One()) {
	  for (size_t i=0;i<n;i++) {
	    //rex = rand()/(RAND_MAX+ScalarFieldTraits<Scalar>::One());
	    rex = Random();
	    //imx = rand()/(RAND_MAX+ScalarFieldTraits<Scalar>::One());
	    imx = Random();
	    pv[i]=complex<Scalar>(rex,imx);
	  }
	}
	else {
	  for (size_t i=0;i<n;i++) {
	    //rex = w*rand()/(RAND_MAX+ScalarFieldTraits<Scalar>::One());
	    rex = w*Random();
	    //imx = w*rand()/(RAND_MAX+ScalarFieldTraits<Scalar>::One());
	    imx = w*Random();
	    pv[i]=complex<Scalar>(rex,imx);
	  }
	}
      } else if (w == 1) {
	for (size_t i=0;i<n;i++) {
	  //rex = a + rand()/(RAND_MAX+ScalarFieldTraits<Scalar>::One());
	  rex = a + Random();
	  //imx = a + rand()/(RAND_MAX+ScalarFieldTraits<Scalar>::One());
	  imx = a + Random();
	  pv[i]=complex<Scalar>(rex,imx);
	}
      } else {
	for (size_t i=0;i<n;i++) {
	  //rex = a + w*rand()/(RAND_MAX+ScalarFieldTraits<Scalar>::One());
	  rex = a + w*Random();
	  //imx = a + w*rand()/(RAND_MAX+ScalarFieldTraits<Scalar>::One());
	  imx = a + w*Random();
	  pv[i]=complex<Scalar>(rex,imx);
	}
      }
    }
    
    string getName() const { return "RVLRandomize<complex>"; }
    
  };

  /** Specialization for assigning random integers in [0, RAND_MAX]
   */
  template<>
  class RVLRandomize<int>: public UnaryLocalFunctionObject<int> {
  private:
    RVLRandomize(const RVLRandomize<int> &) {}
  public:
    RVLRandomize() {}
    ~RVLRandomize() {}

    using RVL::UnaryLocalEvaluation<int>::operator();
    virtual bool readsData(size_t i=0) { return false; }

    void operator() (LocalDataContainer<int> & v) {
      size_t n = v.getSize();
      int * pv = v.getData();
      for (size_t i=0;i<n;i++) {
	pv[i]=rand();
      }
    }

    string getName() const  { return "RVLRandomize<int>"; }
  };

  /** Read a vector from a text file, where the first line is expected to
      contain data on the dimension of the vector.
  */
  template<class Scalar>
  class ASCIIReader: public UnaryLocalFunctionObject<Scalar> {
  private:
    ifstream from;
    string file;
  
    ASCIIReader() {}
    ASCIIReader(const ASCIIReader &) {}
  public:
    ASCIIReader(string fname): from(fname.c_str(),ios::in),file(fname) {
      if (!from) {
	RVLException e; e<<"Error: ASCIIReader constructor\n";
	e<<"failed to open file "<<fname<<"\n";
	throw e;
      }
      // skip the first line as it is simply used to indicate dimension
      string l;
      getline(from,l);
    }
    ~ASCIIReader() { from.close(); }

    using RVL::UnaryLocalEvaluation<Scalar>::operator();
    void operator() (LocalDataContainer<Scalar> & v) {
      size_t n = v.getSize();
      Scalar * pv = v.getData();
      size_t i=0;
      while (i<n && from>>pv[i]) {
	// cout<<"i = "<<i<<" v = "<<pv[i]<<endl;
	i++;
      }
      if (i!=n) {
	RVLException e; e<<"Error: ASCIIReader::operator()\n";
	e<<"failed to read "<<n<<" scalars from file \""<<file<<"\"\n";
	throw e;
      }
    }

    string getName() const  { return "ASCIIReader"; }
  };

  /** Write an LDC to a file, where the first line is expected to
      contain data on the dimension. Note that the argument to operator()
      should really be const. However that would require either adding 
      another FO type or inventing a void RetType and subclassing this 
      from Reduction, neither of which is appealing at all.
  */
  template<class Scalar>
  class ASCIIWriter: public UnaryLocalFunctionObjectConstEval<Scalar> {
  private:
    ofstream to;
    string file;
    ASCIIWriter() {}
    ASCIIWriter(const ASCIIWriter &) {}
  public:
    ASCIIWriter(string fname): to(fname.c_str(),ios::out), file(fname) {
      if (!to) {
	RVLException e; e<<"Error: ASCIIWriter constructor\n";
	e<<"failed to open file "<<fname<<"\n";
	throw e;
      }
    }
    ~ASCIIWriter() { to.close(); }

    using RVL::UnaryLocalConstEval<Scalar>::operator();
    void operator() (LocalDataContainer<Scalar> const & v) {
      size_t n = v.getSize();
      to<<n<<endl;
      Scalar const * const pv = v.getData();
      to.setf(ios::right,ios::adjustfield);
      size_t i=0;
      while (i<n && to<<pv[i]<<endl) i++;
      if (i!=n) {
	RVLException e; e<<"Error: ASCIIWriter::operator()\n";
	e<<"failed to write "<<n<<" scalars to file \""<<file<<"\"\n";
	throw e;
      }
      to.flush();
    }

    string getName() const  { return "ASCIIWriter"; }
  };

  /** Read a vector from a binary file, where the first 
      entry in the file is element first (in bytes) from SEEK_SET.
      16.12.06: first is setable, via function seek(long).
  */
  template<class Scalar>
  class BinaryReader: public UnaryLocalFunctionObjectConstEval<Scalar> {
  private:
    FILE * fp;
    string file;
    mutable long first;
    BinaryReader(){}
    BinaryReader(const BinaryReader<Scalar> &){}
  public:

    BinaryReader(char const * fname, long _first=0L)
      : file(fname), first(_first) {
      if (!(fp=fopen(fname,"r+"))) {
	RVLException e; e<<"Error: BinaryReader constructor\n";
	e<<"failed to open file "<<fname<<"\n";
	throw e;
      }
      if (fseek(fp,first,SEEK_SET)) {
	RVLException e; e<<"Error: BinaryReader constructor\n";
	e<<"failed to seek to specified first position = "<<first<<" in file "<<fname<<"\n";
	throw e;
      }
    }
    BinaryReader(const string fname, long _first=0L)
      : file(fname), first(_first) {
      if (!(fp=fopen(fname.c_str(),"r+"))) {
	RVLException e; e<<"Error: BinaryReader constructor\n";
	e<<"failed to open file "<<fname<<"\n";
	throw e;
      }
      if (fseek(fp,first,SEEK_SET)) {
	RVLException e; e<<"Error: BinaryReader constructor\n";
	e<<"failed to seek to specified first position = "<<first<<" in file "<<fname<<"\n";
	throw e;
      }
    }
    ~BinaryReader() {fclose(fp);}

    /** seek to specified word in file - note: not offset in bytes! */
    bool seek(size_t firstword) {
      first = firstword * sizeof(Scalar);
      if (fseek(fp,first,SEEK_SET)) return false;
      return true;
    }

    using RVL::UnaryLocalEvaluation<Scalar>::operator();
    void operator() (LocalDataContainer<Scalar> & v) {
      size_t n = v.getSize();
      Scalar * pv = v.getData();
      if (n != fread(pv,sizeof(Scalar),n,fp)) {
	RVLException e; e<<"Error: BinaryReader::operator()\n";
	e<<"failed to read "<<n<<" "<<sizeof(Scalar)<<"-byte words starting at position "<<first<<" from file "<<file<<"\n";
	  throw e;
      }
    }

    string getName() const  { return "BinaryReader"; }

  };

  /** Dumps binary data in memory to a binary file
      16.12.06: starting at file position first, setable via
      seek(int).
   */
  template<class Scalar>
  class BinaryWriter: public UnaryLocalConstEval<Scalar> {
  private:
    FILE * fp;
    string file;
    long first;
    BinaryWriter(){}
    BinaryWriter(const BinaryWriter<Scalar> &){}
  public:

    BinaryWriter(char const * fname, long _first=0L)
      : file(fname), first(_first) {
      if (!(fp=fopen(fname,"w+"))) {
	RVLException e; e<<"Error: BinaryWriter constructor\n";
	e<<"failed to open file "<<fname<<"\n";
	throw e;
      }
      if (fseek(fp,first,SEEK_SET)) {
	RVLException e; e<<"Error: BinaryWriter constructor\n";
	e<<"failed to seek to specified first position = "<<first<<" in file "<<fname<<"\n";
	throw e;
      }
    }

    BinaryWriter(const string fname, long _first=0L): 
      file(fname), first(_first) {
      if (!(fp=fopen(fname.c_str(),"w+"))) {
	RVLException e; e<<"Error: BinaryWriter constructor\n";
	e<<"failed to open file "<<fname<<"\n";
	throw e;
      }
      if (fseek(fp,first,SEEK_SET)) {
	RVLException e; e<<"Error: BinaryWriter constructor\n";
	e<<"failed to seek to specified first position = "<<first<<" in file "<<fname<<"\n";
	throw e;
      }
    }
    ~BinaryWriter() {fclose(fp);}

    /** seek to specified word in file - note: not offset in bytes! */
    bool seek(size_t firstword) {
      first = firstword * sizeof(Scalar);
      if (fseek(fp,first,SEEK_SET)) return false;
      return true;
    }

    using RVL::UnaryLocalConstEval<Scalar>::operator();
    void operator() (LocalDataContainer<Scalar> const & v) {
      size_t n = v.getSize();
      Scalar const * pv = v.getData();
      if (n != fwrite(pv,sizeof(Scalar),n,fp)) {
	RVLException e; e<<"Error: BinaryWriter::operator()\n";
	e<<"failed to write "<<n<<" "<<sizeof(Scalar)<<"-byte words starting at position "<<first<<" to file "<<file<<"\n";
	throw e;
      }
      fflush(fp);
    }

    string getName() const  { return "BinaryWriter"; }

  };

  /** Given a maximum and minimum vectors which define a box constraint, find
      the largest scalar so that \f$x+\alpha dx\f$ is not outside the box.
  */
  template<class Scalar>
  class RVLBoxMaxStep: 
    public QuaternaryLocalFunctionObjectScalarRedn<Scalar, Scalar> {

  public:
    RVLBoxMaxStep() {}
    RVLBoxMaxStep(const RVLBoxMaxStep<Scalar> & b) {}
    ~RVLBoxMaxStep() {}

    using RVL::QuaternaryLocalConstEval<Scalar>::operator();    
    void operator()(LocalDataContainer<Scalar> const & x,
		    LocalDataContainer<Scalar> const & dx,
		    LocalDataContainer<Scalar> const & xmin,
		    LocalDataContainer<Scalar> const & xmax) {
      try {
	size_t n = x.getSize();
	if ((n != dx.getSize()) ||
	    (n != xmin.getSize()) ||
	    (n != xmax.getSize())) {
	  RVLException e;
	  e<<"Error: GridMaxStep::operator()\n";
	  e<<"incompatible inputs\n";
	  throw e;
	}
	
	Scalar step = numeric_limits<Scalar>::max();
	Scalar const * px    = x.getData();
	Scalar const * pdx   = dx.getData();
	Scalar const * pxmin = xmin.getData();
	Scalar const * pxmax = xmax.getData();
	
	for (size_t i=0;i<n;i++) {
	  Scalar tmp;
	  if (pdx[i] > ScalarFieldTraits<Scalar>::Zero() && 
	      !ProtectedDivision<Scalar>(pxmax[i]-px[i],pdx[i],tmp)) {
	    if (tmp > ScalarFieldTraits<Scalar>::Zero()) { step = min(step,tmp); }
	    else { step = ScalarFieldTraits<Scalar>::Zero(); }
	  }
	  if (pdx[i] < ScalarFieldTraits<Scalar>::Zero() &&
	      !ProtectedDivision<Scalar>(pxmin[i]-px[i],pdx[i],tmp)) {
	    if (tmp > ScalarFieldTraits<Scalar>::Zero()) { step = min(step,tmp); }
	    else { step = ScalarFieldTraits<Scalar>::Zero(); }
	  }
	}
	setValue(step);
      }
      catch (RVLException & e) {
	e<<"\ncalled from RVLBoxMaxStep::operator()\n";
	throw e;
      }
    }
   
    string getName() const { string tmp = "RVLBoxMaxStep"; return tmp; }
  };

  /**  Performs an elementwise multiplication.  \f$u = v.*w\f$ */
  template<class Scalar>
  class ElementwiseMultiply: public TernaryLocalFunctionObject<Scalar> {

  public:
    ElementwiseMultiply() {}
    ~ElementwiseMultiply() {}
    
    using RVL::TernaryLocalEvaluation<Scalar>::operator();
    void operator() (LocalDataContainer<Scalar> & u, 
		     LocalDataContainer<Scalar> const & v, 
		     LocalDataContainer<Scalar> const & w) {
      Scalar * up = u.getData(); 
      Scalar const * vp = v.getData();
      Scalar const * wp = w.getData();
      size_t size = u.getSize();
      if( (size > v.getSize())||(size > w.getSize())) {
	RVLException e;
	e << "Error in ElementwiseMultiply::operator() - target LDC has more elements than source LDC.\n";
	throw e;
      }
      for(size_t i = 0; i < size; i++)
	up[i] = vp[i]*wp[i];
    }

    string getName() const  { string s = "ElementwiseMultiply"; return s; }
  };

/**  Performs an elementwise safeguarded division. \f$u = v./w\f$ */
  template<class Scalar>
  class ElementwiseDivision: public TernaryLocalFunctionObject<Scalar> {

  public:
    ElementwiseDivision() {}
    ~ElementwiseDivision() {}
    
    using RVL::TernaryLocalEvaluation<Scalar>::operator();
    void operator() (LocalDataContainer<Scalar> & u, 
		     LocalDataContainer<Scalar> const & v, 
		     LocalDataContainer<Scalar> const & w) {
      Scalar * up = u.getData();
      Scalar const * vp = v.getData(); 
      Scalar const * wp = w.getData();
      size_t size = u.getSize();
      if( (size > v.getSize())||(size > w.getSize())) {
	RVLException e;
	e << "Error: ElementwiseDivision::operator()\n";
	e << "target LDC has more elements than source LDC.\n";
	throw e;
      }
      for(size_t i = 0; i < size; i++) {
	if (ProtectedDivision<Scalar>(vp[i],wp[i],up[i])) {
	  RVLException e;
	  e<<"Error: ElementwiseDivision::operator()\n";
	  e<<"zerodivide in element "<<i<<"\n";
	  throw e;
	}
      }
    }

    string getName() const  { string s = "ElementwiseDivision"; return s; }
  };

/**  Performs an elementwise square root of absolute value \f$u =
     sqrt(abs(v))\f$ */
  template<class Scalar>
  class ElementwiseSqrtAbs: public BinaryLocalFunctionObject<Scalar> {

  public:
    ElementwiseSqrtAbs() {}
    ~ElementwiseSqrtAbs() {}
    
    using RVL::BinaryLocalEvaluation<Scalar>::operator();
    void operator() (LocalDataContainer<Scalar> & u, 
		     LocalDataContainer<Scalar> const & v) {
      Scalar * up = u.getData();
      Scalar const * vp = v.getData(); 
      size_t size = u.getSize();
      if (size > v.getSize()) {
	RVLException e;
	e << "Error: ElementwiseSqrtAbs::operator()\n";
	e << "target LDC has more elements than source LDC.\n";
	throw e;
      }
      for(size_t i = 0; i < size; i++) {
	up[i]=sqrt(abs(vp[i]));
      }
    }

    string getName() const  { string s = "ElementwiseSqrtAbs"; return s; }
  };

}
      
#endif



  
