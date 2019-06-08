/*************************************************************************

Copyright Rice University, 2004-2015.
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
#ifndef __RVL_TS
#define __RVL_TS

#include "space.hh"
#include "op.hh"
#include "productspace.hh"
#include "functions.hh"
#include "write.hh"

#define TOL 0.1

namespace RVL {

  /** A sampling operator has the attributes of a LinearOp, but in
      addition an action time range and current time Relies on
      post-construction initialization via "set" functions. The time attributes permit the action of the operator to be time-dependent. For example, setting the initial data at t=t_0 to would require a   */

  template<typename T>
  class TSSample:
    public LinearOp<T> {

  private:
    mutable T t;
    mutable T tmin;
    mutable T tmax;

  public:
    // copy constructor is default (bytewise)
    TSSample() {
      t = ScalarFieldTraits<T>::Zero();
      tmin = ScalarFieldTraits<T>::Zero();
      tmax = ScalarFieldTraits<T>::Zero();
    }

    T const & getMinTime() const { return tmin;}
    T const & getMaxTime() const { return tmax;}
    void setMinTime(T _tmin) const { tmin=_tmin;}
    void setMaxTime(T _tmax) const { tmax=_tmax;}  
    virtual void setTestTime() = 0;
  };

  // forward declaration - see docs for type below
  class TSDC;

  /** Basic class for FunctionObjects that act on a TSDC - interface
      for all time step functions. Unary action because all time step
      functions in their basic form operate in overwrite mode, i.e.
      x=f(x). Subclass must define action on a TSDC - note that
      parent class has no nontrivial methods, so this is significant!
  */
  class TSFO: public FunctionObject {
  public:
 
    // unary operator interface - overwrite mode 
    virtual void operator()(TSDC & y) const = 0;

  };

  /** Basic data container for time stepping control and state
      vectors. Almost always multicomponent, so derived from a product
      container class. The only nontrivial function of this class is
      to provide an interface to intercept unary evaluation of TSFO's,
      which implement time step functions, and usually interface in
      turn to low-level code. Evaluations of other FunctionObjects
      treated via inherited product container implicit loop.
  */
  class TSDC: public StdProductDataContainer {
  public:
    void eval( FunctionObject & f,
	       std::vector<DataContainer const *> & x) {
      try {
	// cerr<<"TSDC::eval\n";
	// cerr<<"  FO = "<<f.getName()<<endl;
	TSFO * g = NULL;
	if ((x.size()==0) && (g = dynamic_cast<TSFO *>(&f))) {
	  (*g)(*this);
	}
	else {
	  // cerr<<"  -->StdProductDataContainer::eval\n";
	  StdProductDataContainer::eval(f,x);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TSDC::eval(FO,...)\n";
	throw e;
      }
    }

    ostream & write(ostream & str) const {
      str<<"TSDC: specialization of \n";
      StdProductDataContainer::write(str);
      return str;
    }
  };

  /** Space whose DataContainer class is a StdProductDC with TSDC
      components. Only attribute distinguishing this from StdProductDC
      is that the buildDataContainer method is guaranteed to return a
      TSDC.

  */
  template<typename T>
  class TSSpace: public StdProductSpace<T> {
  protected:
    Space<T> * clone() const { return new TSSpace<T>(*this); }
  public:
    TSSpace(size_t nfac): StdProductSpace<T>(nfac) {}
    TSSpace(TSSpace<T> const & sp): StdProductSpace<T>(sp) {}
    
    /** overrides method of StdProductDataContainer class.  */
    DataContainer * buildDataContainer() const {
      if (StdProductSpace<T>::init()) {
	TSDC * d = new TSDC();
	for (size_t i=0; i<StdProductSpace<T>::getSize(); i++) {
	  // this construction uses the buildDC method of each
	  // component space to create a factory, then builds a new DC
	  // and pushes it onto the std::vector of DCs underlying d.
	  SpaceDCF<T> f((*this)[i]);
	  d->push(f);
	}
	// cerr<<"*** from TSSpace::buildDataContainer\n";
	//	d->write(cerr);
	return d;
      }
      else {
	RVLException e;
	e<<"ERROR: TSSpace::buildDataContainer\n";
	e<<"  called on uninitialized StdProductSpace object\n";
	throw e;
      }
    }
    
    ostream & write(ostream & str) const {
      str<<"TSSpace: product space with TSDC components\n";
      ProductSpace<T>::write(str);
      return str;
    }
  };

  // forward declarations of friendship
  template<typename T>
  class TimeStepOp;

  template<typename T>
  class LinRestrictTSStep;

  template<typename T>
  class TanRestrictTSStep;

  /** Builds time steps out of special function objects (TSFOs), which
      access DataContainers directly hence provide convenient wrappers
      for low-level numerical kernels in C. TSFOs have unary
      evaluation interface consistent with overwrite style of typical
      time step kernel. The binary interface for the main protected
      virtual methods (apply,...) is required by the parent class
      (Operator) and used in standard tests (DerivTest, AdjointTest)
      so must be implemented, but only for test purposes. The binary
      interface works by copying the input over the output then
      evaluating the TSFO on the output - inefficient, but OK for
      tests. The unnecessary copy does not take place if the input and
      output are aliased (that is, the same object) which is the main
      use case.

      Mathematically, a TSStep implements an update rule on the pair
      (c,s) consisting of a control c and state s, and relative
      derivative and adjoint calculations. The rule is linear in s and
      parametrized by c, so c is not updated:

      (c,s) -> (c, H(c)s)

      That is, the input and output of apply are both pairs: in 

      apply(x,y);

      any implementation should expect x=(c,s), y=(c',s'). The domain
      space of a TSStep consists of vectors of this form, so the space
      membership checked by OperatorEvaluation, TimeStepOp, and so on
      assures this structure.

      Unlike the generic operator, TSStep has a natural adjoint mode,
      namely 

      (c,s) -> (c, H(c)^T s)

      This method (applyAdj) is thus not inherited from Operator.

      Protected methods have no type-checking because they are only
      accessed by parent Operator, by OperatorEvaluation (through
      export calls to Operator), and by friends TanRestrictTSStep,
      LinRestrictTSStep, and TimeStepOp, which sanity-check arguments
      before call.

      To conform to typical structure of time-stepping
      implementations, the virtual access to derivative is through a
      tangent map (that is, tangent vector to tangent vector)

      ((c,s), (dc,ds)) -> ((c,H(c)s),(dc,(DH(c)dc)s + H(c)ds))

      rather than the derivative as required by Operator, which would
      amount to

      ((c,s), (dc,ds)) -> (dc,(dc,(DH(c)dc)s + H(c)ds))

      that is, the second component of the tangent map. Similar for
      adjoints. TSStep provides an obvious implementation of
      applyDeriv and applyAdjDeriv in terms of applyTangent and
      applyTangentAdj for test purposes. The main use case, in
      implementation of TimeStepOp and related, uses only the
      applyTangent and applyTangentAdj interfaces.

      This is an abstract base class: instantiable subclass requires
      definitions of

      * FOs for apply methods
      * getDomain and getRange: domain and range space access
      * stepTimeFwd and stepTimeBwd methods
      * constructor/destructor
  */
  template<typename T>
  class TSStep:
    public Operator<T> {
    
    friend class TimeStepOp<T>;
    friend class LinRestrictTSStep<T>;
    friend class TanRestrictTSStep<T>;
    
  private:
    virtual TSFO & applyFO() const = 0;
    virtual TSFO & applyAdjFO() const = 0;
    virtual TSFO & applyTangentFO() const = 0;
    virtual TSFO & applyAdjTangentFO() const = 0;
    
  protected:

    /** inherited from Operator. In this an other apply methods, the
     interface is formally binary (or with more arguments, for the
     other methods), with two intended use cases: <ul> <li>x and y are
     same object (checked by pointer), then the applyFO method is
     called to update the output y - this is the typical use case, as
     applyFO operates in overwrite mode by (usually) calling low-level
     code </li> <li>x and y are different, in which case x is copied
     onto y then applyFO is called on y. This mode is provided for
     test purposes, eg. to check applyDeriv by finite difference
     comparison (DerivCheck function).</li> </ul>*/
    void apply(Vector<T> const & x,
	       Vector<T> & y) const {
      try {
	//	cerr<<"  AStep::apply -> copy\n";
	if (&y != &x) y.copy(x);
	// cerr<<"  AStep::apply -> eval(applyFO()) FO = "<<applyFO().getName()<<"\n";
	//	cerr<<"  input vector:\n";
	//	y.write(cerr);
	y.eval(applyFO());
	// cerr<<"  AStep::apply -> stepTimeFwd()\n";	
      }
      catch (RVLException & e) {
	e<<"\ncalled from TSStep::apply\n";
	throw e;
      }
    }

    void applyAdj(Vector<T> const & x,
		  Vector<T> & y) const {
      try {
	if (&y != &x) y.copy(x);
	y.eval(applyAdjFO());
      }
      catch (RVLException & e) {
	e<<"\ncalled from TSStep::apply\n";
      }
    }

    /** tangent map: (c,s,dc,ds) -> (c,F[c,s],dc,DF[c,s][dc,ds]) 
	= (c,H[c]s, dc, (DH[c]dc)s + H[c]ds)
	combine args into product space: (x,dx) = ((c,s),(dc,ds))
	ydy = (y,dy) = (s,ds)

	only TimeStepOp accesses the Tangent methods, so sanity checking
	deferred to there. 
    */
    void applyTangent(Vector<T> const & xdx,
		      Vector<T> & ydy) const {
      try {
	if (&ydy != &xdx) ydy.copy(xdx);
	ydy.eval(applyTangentFO());
      }
      catch (RVLException & e) {
	e<<"\ncalled from TSStep::apply\n";
      }
    }

    /** adjoint tangent map: (c,s,dc,ds) -> (c,H[c]^{-1}s,dc + ((DH[c]
	*)H[c]^{-1}s)^T dc, H[c]^{T} ds) cannot implement H[c]^{-1} in
	general, since H[c] need not be invertible. Instead, assert
	that (c,s) is value of base trajectory at t+dt, and interpret
	(c,H[c]^{-1}s) as state at time t. Therefore this computation
	requires two stages: (1) time reverse, t -> t-dt, and
	extraction of state - method for random access to trajectory,
	TSFO/TSClock; and stationary adjoint tangent map, implementing
	(c,s,dc,ds) -> (c,s,dc + ((DH[c]*)s)^T dc, H[c]^{T} ds). This
	method implements the second step; the first is built into
	TimeStepOp, since it naturally requires access to the time
	loop to build up a cache.
    */
    void applyAdjTangent(Vector<T> const & xdx,
		      Vector<T> & ydy) const {
      try {
	if (&ydy != &xdx) ydy.copy(xdx);
	ydy.eval(applyAdjTangentFO());
      }
      catch (RVLException & e) {
	e<<"\ncalled from TSStep::apply\n";
      }
    }
    
    /** inherited from Operator

	implement via tangent map 
	x   = (c,s)^T
	dx  = (dc,ds)^T
	dy  = (dc,DF[c,s][dc,ds]) = (dc, (DH[c]dc)s + H[c]ds)
	= M[c,s][dc ds]^T

	M[c,s] = [ I            0   ]
                 [ (DH[c] *)s  H[c] ]

        this method and applyAdjDeriv are useful mostly for test purposes

    */
    void applyDeriv(Vector<T> const & x,
		    Vector<T> const & dx,
		    Vector<T> & dy) const {
      try {
	StdProductSpace<T> TX(this->getDomain(), this->getDomain());
	Vector<T> xdx(TX);
	Vector<T> ydy(TX);
	Components<T> cxdx(xdx);
	Components<T> cydy(ydy);
	cxdx[0].copy(x);
	cxdx[1].copy(dx);
	applyTangent(xdx,ydy);
	dy.copy(cydy[1]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from TSStep::applyDeriv\n";
      }
    }
    
    /** inherited from Operator

	x = (c,s)^T
	dx  = (dc,ds)^T
	dy  = M[c,s]^Tdx
	M[c,s]^T = [ I           ((DH[c] *)s)^T ]
                   [ 0                   H[c]^T ]
    */
    void applyAdjDeriv(Vector<T> const & x,
		       Vector<T> const & dx,
		       Vector<T> & dy) const {
      try {
	StdProductSpace<T> TX(this->getDomain(), this->getDomain());
	Vector<T> xdx(TX);
	Vector<T> ydy(TX);
	Components<T> cxdx(xdx);
	Components<T> cydy(ydy);
	cxdx[0].copy(x);
	cxdx[1].copy(dx);
	applyAdjTangent(xdx,ydy);
	dy.copy(cydy[1]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from TSStep::applyAdjDeriv\n";
      }
    }

    // does not need to be actual time step could be total step for
    // method with substeps - in any case is the step between cache
    // time levels for adjoint state loops
    virtual T getTimeStep() const = 0;

    /** move t forward by appropriate step (possibly partial) */
    virtual void stepTimeFwd(T & t) const = 0;
    /** move t backward by appropriate step (possibly partial) */
    virtual void stepTimeBwd(T & t) const = 0;
    
  public:

    TSStep() {}
    TSStep(TSStep<T> & ts) {}
    virtual ~TSStep() {}

    TSSpace<T> const & getTSDomain() const {
      try {
	//	cerr<<"in TSStep::getTSDomain\n";
       	TSSpace<T> const & pdom =
	  dynamic_cast<TSSpace<T> const &>(this->getDomain());
	//////////////////////////////////
	//pdom.write(cerr);
	//////////////////////////////////	
	return pdom;
      }
      catch (bad_cast) {
	RVLException e;
	e<<"Error: TSStep::getTSDomain\n";
	e<<"  it aint a PD\n";
	this->getDomain().write(e);
	throw e;
      }
    }
    
  };

  /** restriction to state component of (ctrl,state) TSSpace pair,
      presumed to be linear
  */
  template<typename T>
  class LinRestrictTSStep: public LinearOp<T> {
  private:
    TSStep<T> const & step;
    mutable Vector<T> xbuf; mutable Components<T> cxbuf;
    mutable Vector<T> ybuf; mutable Components<T> cybuf;    
    
  protected:

    void apply(const Vector<T> & x,
	       Vector<T> & y) const {
      try {
	// input x0 is already stored in cxbuf[0]
	cxbuf[1].copy(x);
	step.apply(xbuf,ybuf);
	y.copy(cybuf[1]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinRestrictTSStep::apply\n";
	throw e;
      }
    }
    
    virtual void applyAdj(const Vector<T> & x,
			  Vector<T> & y) const {
      try {
	// input x0 is already stored in cxbuf[0]
	cxbuf[1].copy(x);
	step.applyAdj(xbuf,ybuf);
	y.copy(cybuf[1]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinRestrictTSStep::applyAdj\n";
	throw e;
      }
    }

    Operator<T> * clone() const {
      return new LinRestrictTSStep<T>(*this);
    }

  public:

    LinRestrictTSStep(TSStep<T> const & _step,
		      Vector<T> const & x0)
      : step(_step),
	xbuf(step.getTSDomain()),
	cxbuf(xbuf),
	ybuf(step.getTSDomain()),
      	cybuf(ybuf) {
      try {
	// this catches any mismatch between x0 and
	// dom[0] - other mismatches caught by LinearOp interface
	cxbuf[0].copy(x0);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinRestrictTSStep constructor\n";
	e<<"input x0:\n";
	x0.write(e);
	e<<"output comp[0]:\n";
	cxbuf[0].write(e);
	throw e;
      }
    }

    LinRestrictTSStep(LinRestrictTSStep<T> const & a)
      : step(a.step),
	xbuf(step.getTSDomain()),
	cxbuf(xbuf),
	ybuf(step.getTSDomain()),
	cybuf(ybuf) {	
      try {
	cxbuf[0].copy(a.cxbuf[0]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinRestrictTSStep constructor\n";
	throw e;
      }
    }

    Space<T> const & getDomain() const { return step.getTSDomain()[1]; }
    Space<T> const & getRange() const { return step.getTSDomain()[1]; }

    ostream & write(ostream & str) const {
      str<<"LinRestrictTSStep linear op: restriction to comp 1\n";
      str<<"  of TSStep operator\n";
      step.write(str);
      str<<"  with comp 0 restricted to\n";
      cxbuf[0].write(str);
      return str;
    }
  };

  /** linear map defined by tangent action - for tests */
  template<typename T>
  class TanRestrictTSStep: public LinearOp<T> {
  private:
    TSStep<T> const & step;
    TSSpace<T> tdom;
    std::shared_ptr<Vector<T> > xbuf;
    std::shared_ptr<Vector<T> > ybuf;
    
  protected:

    void apply(const Vector<T> & x,
	       Vector<T> & y) const {
      try {
	// input x0 is already stored in comp[0]
	Components<T> cxbuf(*xbuf);
	Components<T> cybuf(*ybuf);
	cxbuf[1].copy(x);
	step.applyTangent(*xbuf,*ybuf);
	y.copy(cybuf[1]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from TanRestrictTSStep::apply\n";
	throw e;
      }
    }
    
    virtual void applyAdj(const Vector<T> & x,
			  Vector<T> & y) const {
      try {
	// input x0 is already stored in comp[0]
	Components<T> cxbuf(*xbuf);
	Components<T> cybuf(*ybuf);
	cxbuf[1].copy(x);
	step.applyAdjTangent(*xbuf,*ybuf);
	y.copy(cybuf[1]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from TanRestrictTSStep::applyAdj\n";
	throw e;
      }
    }

    Operator<T> * clone() const {
      return new TanRestrictTSStep<T>(*this);
    }

  public:

    TanRestrictTSStep(TSStep<T> const & _step,
		      Vector<T> const & x0)
      : step(_step), tdom(2) {
      try {
	tdom.set(step.getTSDomain(),0);
	tdom.set(step.getTSDomain(),1);
	xbuf = make_shared<Vector<T> >(tdom);
	ybuf = make_shared<Vector<T> >(tdom);
	Components<T> cxbuf(*xbuf);
	cxbuf[0].copy(x0);
	cxbuf[1].zero();
      }
      catch (RVLException & e) {
	e<<"\ncalled from TanRestrictTSStep constructor\n";
	e<<"input x0:\n";
	x0.write(e);
	e<<"output:\n";
	xbuf.write(e);
	throw e;
      }
    }

    TanRestrictTSStep(TanRestrictTSStep<T> const & a)
      : step(a.step),
	tdom(a.tdom) {
      try {
	xbuf = make_shared<Vector<T> >(tdom);
	ybuf = make_shared<Vector<T> >(tdom);
	xbuf->copy(a.xbuf);
      }
      catch (RVLException & e) {
	e<<"\ncalled from TanRestrictTSStep constructor\n";
	throw e;
      }
    }

    Space<T> const & getDomain() const { return step.getTSDomain(); }
    Space<T> const & getRange() const { return step.getTSDomain(); }

    ostream & write(ostream & str) const {
      str<<"TanRestrictTSStep linear op: restriction to comp 1\n";
      str<<"  of TSStep operator\n";
      step.write(str);
      str<<"  with comp 0 restricted to\n";
      Components<T> cxbuf(*xbuf);
      cxbuf[0].write(str);
      return str;
    }
  };

  /** 
      LinOpValOp interface to generic time loop based on a TSStep
      object and input, output sampling operators. 

      LinOpValOp defines a linear operator valued (nonlinear)
      function, or operator, with explicit access both the the
      argument of the nonlinear function and to the linear operator
      which is the value of the nonlinear function. 

      In this application, I call the argument of the nonlinear
      function the control, or control vector. It typically represents
      the coefficients of a linear partial differential operator. The
      linear operator, obtained by evaluating the TimeStepOp at a
      control, implements the source-to-data relation: the input is
      the source, or right-hand side, represented by the source
      sampling operator, and the output is the data, or sampling of
      the solution, represented by the data sampling operator. Both
      sampling operators act on a state field - in this application,
      the solution field of the partial differential equation. 

      The source and data sampling operators are LinearOps. Each also
      maintains a current time, controlled through a setTime function,
      and also max and min times. The action of these operators may
      vary with the time: for example, initial data is provided by a
      source operator that acts only at the minimum time, and is a
      no-op at other times.

      The source sampling operator (and in adjoint mode the
      data sampling operator) typically increments the state, that is,

      state += source(source_vector)

      RVL::LinearOp supplies a applyPlusOp interface with the obvious
      default implementation, which may be inefficient: specific
      source operator classes can override the default to minimize the
      use of intermediate storage and unnecessary passes through data
      (for example, by only visiting those state field values actually
      updated by a source with limited extent).

      The state is not visible in the main class methods of
      TimeStepOp: it is an attribute of the TSStep data member, which
      implements the time step evolution of the state. The source and
      data sampling operators are linear operators, whose range,
      respectively domain, is the state space.

      The basic time-stepping loop is the body of the apply0 method -
      see the doc for this method for a detailed explanation, and the
      docs for LinOpValOp for an overview of the class, its structure
      and rationale.
*/
  template<typename T>
  class TimeStepOp: public LinOpValOp<T> {
  private:

    TSSample<T> & data;
    StdProductSpace<T> pdom;   // product domain of LOVOp
    TSSpace<T> tdom;   // ctrl-state domain of step tangent map
    // for tests
    mutable RVLRandomize<T> rnd;
    mutable ofstream report;
    bool testflag;

  protected:

    // must be available to children
    TSStep<T> & step;
    TSSample<T> & src;
    
    /** time advance method, provided for use in 
	caches */
    void nowandthen(Vector<T> const & ctrlstate,
		    Vector<T> const & x1,
		    T tstart, T tstop) const {
      try {
	// sanity test should be implicit
	Components<T> comp_ctrlstate(ctrlstate);
	T t = tstart;
	T tlast = t;
	step.setTime(t);
	while (t + TOL*(t-tlast) < tstop) {
	  src.setTime(t);
	  src.applyPlusOp(x1,comp_ctrlstate[1]);
	  step.apply(ctrlstate,ctrlstate);
	  tlast = t;
	  step.stepTimeFwd(t);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TimeStepOp::apply0\n";
	throw e;
      }
    }

    /** Save state at time t. default implementation of save does
	nothing */
    virtual void save(Vector<T> const & x, T t) const {}

    /** Load state at time t. default implementation of load throws
	exception */
    virtual void load(Vector<T> & x, T t) const {
      RVLException e;
      e<<"Error: TimeStepOp::load\n";
      e<<"  no default implementation provided\n";
      e<<"  to employ must implement in child class\n";
      throw e;
    }
    
    /** This function implements the basic timestepping loop in the
	TimeStepOp framework. Glossary:
	<ul>
	<li>x0 = control vector (nonlinear input)</li>
	<li>x1 = source vector (linear input)</li>
	<li>y = data vector (linear output)
	
	<li>ctrlstate = product vector (control,state), domain/range
	vector of TSStep data member. </li>
	
	<li>comp_ctrlstate = components of ctrlstate - x0 is
	compatible with comp_ctrlstate[0] (that is, a member of the
	same RVL::Space) as is automatically checked by copying x0
	onto comp_ctrlstate[0] - an exception will be thrown if these
	vectors are not compatible.</li>

	</ul>

	The state vector is always initialized to zero: further
	initialization may be supplied by the source operator, for
	instance.

	The time interval of the simulation begins with the lesser of
	the min times from source and data sampling operators, and
	ends at the max time of the data sampling operator. This works
	because the sampling operators implement no-ops outside of
	their time intervals.
    */
    void apply0(Vector<T> const & x0,
		Vector<T> const & x1,
		Vector<T> & y) const {
      try {
	// vector in domain of TSStep - (ctrl,state)
	Vector<T> ctrlstate(step.getTSDomain());
	// access to components
	Components<T> comp_ctrlstate(ctrlstate);
	// input x0 is control - copy it in - implicit sanity test
	// is it possible to SAFELY reference outside data and avoid
	// a copy? YES! if components are recorded in smart pointers,
	// as in newrvl
	comp_ctrlstate[0].copy(x0);
	// always zero state vector
	comp_ctrlstate[1].zero();
	// safety - zero output
	y.zero();
	// simulation begins at min of source, data min times
	T t = min(src.getMinTime(), data.getMinTime());
	T tend = data.getMaxTime();
	T tlast = t;
	while (t + TOL*(t-tlast) < tend) {
	  // save before update
	  save(comp_ctrlstate[1]);
	  src.setTime(t);
	  src.applyPlusOp(x1,comp_ctrlstate[1]);
	  step.apply(ctrlstate,ctrlstate);
	  tlast = t;
	  step.stepTimeFwd(t);
	  data.setTime(t);
	  data.applyOp(comp_ctrlstate[1],y);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TimeStepOp::apply0\n";
	throw e;
      }
    }

    void applyAdj0(const Vector<T> & x0,
		   const Vector<T> & x1, 
		   Vector<T> & y) const {
      try {
	Vector<T> ctrlstate(step.getTSDomain());
	Components<T> comp_ctrlstate(ctrlstate);
	comp_ctrlstate[0].copy(x0);
	comp_ctrlstate[1].zero();

	T t = data.getMaxTime();
	step.setTime(t);
	T tbeg = min(src.getMinTime(), data.getMinTime());
	T tnext = t;
	while (t - TOL*(tnext-t) > tbeg) {
	  data.setTime(t);
	  data.applyPlusAdjOp(x1,comp_ctrlstate[1]);
	  step.stepTimeBwd();
	  step.applyAdj(ctrlstate,ctrlstate);
	  tnext = t;
	  t = step.getTime();
	  src.setTime(t);
	  src.applyPlusAdjOp(comp_ctrlstate[1],y);
	}
      }	
      catch (RVLException & e) {
	e<<"\ncalled from TimeStepOp::applyAdj0\n";
	throw e;
      }
    }

    // here x0 is control, x1 is external source, y is external data
    // state does not appear in arg list!
    /** partial deriv wrt control. Dynamics:
	s_{n+1} = H[c]s_n + src_n
	differentiated wrt c:
	ds_{n+1} = (DH[c]dc)s_n + H[c]ds_n
	Two options:

	(1) use LinOpValOp methods:
	ds = H[c]ds
	ds += (DH[c]dc)s ropeval.getDeriv().applyPlusOp(s,ds), 
	rop = RestrictOp(step,s) = (c \rightarrow H[c]s)
	however you slice it, this requires two passes through ds data

	(2) use Op methods
	view op as F[c,s] = H[c]s
	DF[c,s](dc,ds) = (DH[c]dc)s + H[c]ds
	this preserves opportunities for loop fusion, but because (c,s) is
	updated, must be careful that aux data does not get recomputed 
	every step

	Also note that F[c,s_n] = s_{n+1}, that is step map is from 
	(c,s) -> s.

    */
    void applyPartialDeriv0(const Vector<T> & x0,
			    const Vector<T> & x1,
			    const Vector<T> & dx0,
			    Vector<T> & dy) const {
      try {

	Vector<T> tangstate(tdom);
	Components<T> comp_tangstate(tangstate);
	Vector<T> & ctrlstate  = comp_tangstate[0];
	Vector<T> & dctrlstate = comp_tangstate[1];
	Components<T> comp_ctrlstate(ctrlstate);
	comp_ctrlstate[0].copy(x0);
	comp_ctrlstate[1].zero();
	Components<T> comp_dctrlstate(dctrlstate);
	comp_dctrlstate[0].copy(dx0);
	comp_dctrlstate[1].zero();

	dy.zero();

	T t = min(src.getMinTime(), data.getMinTime());
	step.setTime(t);
	T tend = data.getMaxTime();
	T tlast = t;
	while (t + TOL*(t-tlast) < tend) {
	  save(comp_ctrlstate[1]);
	  src.setTime(t);
	  src.applyPlusOp(x1,comp_ctrlstate[1]);
	  step.applyTangent(tangstate,tangstate);
	  step.stepTimeFwd();
	  tlast = t;
	  t = step.getTime();
	  data.setTime(t);
	  data.applyOp(comp_dctrlstate[1],dy);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TimeStepOp::applyPartialDeriv0\n";
	throw e;
      }
    }
    
    void applyAdjPartialDeriv0(const Vector<T> & x0,
			       const Vector<T> & x1,
			       const Vector<T> & dy,
			       Vector<T> & dx0) const {
      try {

	Vector<T> tangstate(tdom);
	Components<T> comp_tangstate(tangstate);
	Vector<T> & ctrlstate  = comp_tangstate[0];
	Vector<T> & dctrlstate = comp_tangstate[1];
	Components<T> comp_ctrlstate(ctrlstate);
	comp_ctrlstate[0].copy(x0);
	Components<T> comp_dctrlstate(dctrlstate);
	comp_dctrlstate[0].zero();
	comp_dctrlstate[1].zero();
	dx0.zero();

	T t = data.getMaxTime();
	step.setTime(t);
	T tbeg = data.getMinTime();
	T tnext = t;
	while (t - TOL*(tnext-t) > tbeg) {
	  data.setTime(t);
	  data.applyPlusAdjOp(dy,comp_dctrlstate[1]);
	  step.stepTimeBwd();
	  load(comp_ctrlstate[1]);
	  step.applyAdjTangent(tangstate,tangstate);	  
	  tnext = t;
	  t = step.getTime();
	}
	dx0.copy(comp_dctrlstate[0]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from TimeStepOp::applyPartialDeriv0\n";
	throw e;
      }
    }
    
    OperatorProductDomain<T> * clonePD() const {
      return new TimeStepOp(*this);
    }
    
  public:
    TimeStepOp(TSStep<T> & _step,
	       TSSample<T> & _src, 
	       TSSample<T> & _data,
	       bool _testflag = false)
      : step(_step), src(_src), data(_data),
	pdom(2),
	tdom(2),
	rnd(getpid(),-ScalarFieldTraits<T>::One(),
	    ScalarFieldTraits<T>::One()),
	testflag(_testflag) {
      pdom.set((step.getTSDomain())[0],0);
      pdom.set(src.getDomain(),1);
      tdom.set(step.getTSDomain(),0);
      tdom.set(step.getTSDomain(),1);
      
      if (testflag) {
	std::stringstream tmp;
	tmp << getpid();
	std::string fn = "TimeStepOp_TestReportPID" + tmp.str() + ".txt";
	report.open(fn,std::ofstream::app);
      }
    }
    TimeStepOp(TimeStepOp<T> const & op)
      : step(op.step), src(op.src), data(op.data),
	pdom(2), tdom(2),
	rnd(getpid(),-ScalarFieldTraits<T>::One(),
	    ScalarFieldTraits<T>::One()),
	testflag(op.testflag) {
      pdom.set((step.getTSDomain())[0],0);
      pdom.set(src.getDomain(),1);
      tdom.set(step.getTSDomain(),0);
      tdom.set(step.getTSDomain(),1);
      if (testflag) {
	std::stringstream tmp;
	tmp << getpid();
	std::string fn = "TimeStepOp_TestReportPID" + tmp.str() + ".txt";
	report.open(fn,std::ofstream::app);
      }
    }
    
    ~TimeStepOp() { report.close(); }

    ProductSpace<float> const & getProductDomain() const { return pdom; }
    Space<float> const & getRange() const { return data.getRange(); }

    // adjoint tests - x0 is arg so that any necessary constraints
    // may be obeyed - design assumes that 
    bool testSrcAdj() {
      if (testflag) {
	report<<"\n\n*** TimeStepOp: Src Adjoint Test\n\n";
	src.setTestTime();
	return AdjointTest(src,rnd,report);
      }
      else {
	RVLException e;
	e<<"Error: TimeStepOp::testSrcAdj\n";
	e<<"  attempt to perform test with testflag not set\n";
	e<<"  testflag must be set on construction\n";
	throw e;
      }
    }

    bool testDataAdj() {
      if (testflag) {
	report<<"\n\n*** TimeStepOp: Data Adjoint Test\n\n";
	data.setTestTime();
	return AdjointTest(data,rnd,report);
      }
      else {
	RVLException e;
	e<<"Error: TimeStepOp::testSrcAdj\n";
	e<<"  attempt to perform test with testflag not set\n";
	e<<"  testflag must be set on construction\n";
	throw e;
      }
    }

    bool testStepAdj0(Vector<T> const & ctrl) {
      if (testflag) {
	report<<"\n\n*** TimeStepOp: Step Adjoint Test Order 0\n\n";
	//	cerr<<"******* before call to LinRestrictOp constr\n";
	//	cerr<<"******* ctrl:\n";
	//	ctrl.write(cerr);

	LinRestrictTSStep<T> linstep(step,ctrl);
	return AdjointTest(linstep,rnd,report);
      }
      else {
	RVLException e;
	e<<"Error: TimeStepOp::testStepAdj0\n";
	e<<"  attempt to perform test with testflag not set\n";
	e<<"  testflag must be set on construction\n";
	throw e;
      }
    }

    bool testStepAdj1(Vector<T> const & ctrl) {
      if (testflag) {
	report<<"\n\n*** TimeStepOp: Step Adjoint Test Order 0\n\n";
	//	cerr<<"******* before call to LinRestrictOp constr\n";
	//	cerr<<"******* ctrl:\n";
	//	ctrl.write(cerr);
	Vector<T> cs(step.getTSDomain());
	Components<T> ccs(cs);
	ccs[0].copy(ctrl);
	ccs[1].eval(rnd);
	TanRestrictTSStep<T> tanstep(step,cs);
	return AdjointTest(tanstep,rnd,report);
      }
      else {
	RVLException e;
	e<<"Error: TimeStepOp::testStepAdj0\n";
	e<<"  attempt to perform test with testflag not set\n";
	e<<"  testflag must be set on construction\n";
	throw e;
      }
    }

    bool testAll(Vector<T> const & x0) {
      if (testflag) {
	return (this->testSrcAdj() &&
		this->testDataAdj() &&
		this->testStepAdj0(x0));
      }
      else {
	RVLException e;
	e<<"Error: TimeStepOp::testAll\n";
	e<<"  attempt to perform test with testflag not set\n";
	e<<"  testflag must be set on construction\n";
	throw e;
      }
    }
    
    ostream & write(ostream & str) const {
      str<<"TimeStepOp: \n";
      str<<"*** time step \n";
      step.write(str);
      str<<"*** source sampler (input) \n";
      src.write(str);
      str<<"*** data sampler (output)\n";
      data.write(str);
      return str;
    }
  };

  template<typename T>
  class TimeStepOpAllCache: public TimeStepOp<T> {

  private:

    size_t count(T t,T dt) {
      return size_t((t+.1*dt)/dt);
    }
    size_t nt;
    T dt;
    T tmin; T tmax;
    StdProductSpace<T> cachesp;
    mutable std::shared_ptr<Vector<T> > cache;
    mutable std::shared_ptr<Components<T> > comp_cache;
    mutable bool init;
    
  protected:

    void save(Vector<T> const & x, T t) const {
      try {
	int it = count(t-tmin,dt);	
	(*comp_cache)[it].copy(x);
      }
      catch (RVLException & e) {
	e<<"\ncalled from TimeStepOpAllCache::apply\n";
	throw e;
      }
    }

    void load(Vector<T> & x, T t) const {
      try {
	int it = count(t-tmin,dt);	
	x.copy((*comp_cache)[it]);
      }
      catch (RVLException & e) {
	e<<"\ncalled from TimeStepOpAllCache::apply\n";
	throw e;
      }
    }

    LinOpValOp<T> * clone() const {
      return new TimeStepOpAllCache(*this);
    }
    
  public:
    TimeStepOpAllCache(TSStep<T> & _step,
		       TSSample<T> & _src, 
		       TSSample<T> & _data,
		       bool _testflag = false)	      
      : TimeStepOp<T>(_step, _src, _data, _testflag),
      dt(_step.getTimeStep()),
      tmin(min(_src.getMinTime(),_data.getMinTime())),
      tmax(_data.getMaxTime()),
      nt(this->count(tmax-tmin,_step.getTimeStep())), cachesp(nt) {
      for (size_t i; i<nt; i++) cachesp.set(_step.getTSDomain()[1],i);
      cache = make_shared<Vector<T> >(cachesp);
      comp_cache = make_shared<Components<T> >(*cache);
    }
    TimeStepOpAllCache(TimeStepOpAllCache<T> const & x)
      : TimeStepOp<T>(x), nt(x.nt), tmin(x.tmin), tmax(x.tmax), cachesp(x.cachesp) {
      cache = make_shared<Vector<T> >(cachesp);
      comp_cache = make_shared<Components<T> >(*cache);
    }

    ostream & write(ostream & str) const {
      TimeStepOp<T>::write(str);
      str<<"TimeStepOpAllCache: TimeStepOp subtype, random initialization and\n";
      str<<"access to vectors in power of domain\n";
      str<<"  number of copies = "<<nt<<"\n";
      return str;
    }
  };

}

#endif
