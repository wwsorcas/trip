/** \mainpage The RVL Time Stepping Operator Class 

Author: William W. Symes<br>
<p>
<h2>Mathematical Setting</h2>
<p>

RVL::TimeStepOp and related classes define simulation operators
derived from time stepping loops. 
<p>
The mathematical setting is based on a several Hilbert spaces:
<ul>
<li> a control space \f$C\f$</li>
<li> a state space \f$S\f$</li>
<li> a source space \f$W\f$</li>
<li> a data space \f$D\f$</li>
<li> time, represented by the real line \f${\bf R}\f$
</ul>
and several operators connecting them:
<ul>
<li> time step \f$ H: C \times S \rightarrow S \f$</li>
<li> time update rule: \f$ h: {\bf R} \rightarrow {\bf R}\f$</li>
<li> source sampling operator \f$R: W \rightarrow S\f$</li>
<li> data sampling operator \f$P: S \rightarrow D\f$</li>
</ul>
The data and source sampling operators are also time-dependent. For instance, a source operator could be constructed that initializes the state vector at time \f$t_0\f$, and at other times is a no-op. Similarly, a seismic trace data sampling operator updates different parts of the traces at different times. 
<p>
In view of the interactions with time, we revise the operator definitions as follows:
<ul>
<li> source sampling operator \f$R: W \times {\bf R} \rightarrow S\f$</li>
<li> data sampling operator \f$P: S \times {\bf R} \rightarrow D\f$</li>
</ul>
<p>
The current version of TSOp relies on several simplifying assumptions that are natural for seismic simulation. The action of the time step \f$H\f$ on the state vector is presumed independent of time, that is, the dynamics modeled by TSOp are autonomous. The only interaction of \f$H\f$ is to update it. \f$H\f$ is presumed to be linear in \f$S\f$, (thus representing linear dynamics), so appropriate notation is
\f$ (c, s) \mapsto H(c)s \f$. The sampling operators are also presumed linear.
<p>
The time loop at the heart of this design creates a sequence \f$ (s_i \in S, t_i \in {\bf R}): i=0,...N\f$, along with a sequence of modifications of a data vector \f$d \in D\f$ via the recursion:
<p>
\f$ s_{i+1} = H(c)s_i + R(t_i) w\f$
<p>
\f$ t_{i+1} = h(t_i)\f$
<p>
\f$ d += P(t_{i+1})s_{i+1}\f$
<p>

The simulation operator implicitly defined by this recursion is \f$F: C \times W \rightarrow D\f$. Since this operator is manifestly linear in \f$ W \f$ given the assumptions above, write \f$ F(c)w \f$.
<p>
If \f$ H \f$ is differentiable, then so is \f$ F \f$. The TSOp class includes methods implementing first derivatives and their adjoints, as well as the adjoint of \f$ w \mapsto F(c)w \f$. Second derivatives should probably be added at some point, and derivatives of arbitrary order could certainly be coded in a similar style.
<p>
<h2>Class Design</h2>
<p>
The main design objectives are:

<ol>

<li> RVL::TimeStepOp defines the simulation operator via a time
step loop, from inputs (coefficients = control vector, right hand side
= source vector) to outputs (data vector), and its derivatives of
orders 1 and 2 and their adjoints</li>

<li>the composition of RVL::TimeStepOp in terms of component computations
(time step, sampling operators) should be explicit and visible in the
code</li>

<li>every major component  has an RVL interface (Operator, LinearOp)
and therefore can be tested separately.</li>

<li>the main class provides a set of tests for its principal
components. Passing these <i>guarantees</i> that the simulation
operators created from the components stand in correct relations, that
is, pass similar tests: the linearized operators are actually
derivatives of the basic operators, the adjoints satisfy the dot
product test, and so on. Including the component tests as methods of
the main class assures in any test, compatible
instantiations of component objects are tested. </li>

<li>admit more or less arbitrary source and data sampling operators</li>

<li>stencil calculations for the time step operator may be implemented
by optimized low level code</li> </ol>

Top level design decisions flowing from these goals:

<ol>
<li>major component classes are RVL::TSSample (sampling operators,
subclass of RVL::LinearOp) and RVL::TSStep (time step operator, subclass of
RVL::Operator)</li>

<li>RVL::TSSample adds to the RVL::LinearOp interface "set" and "get"
methods for current time and time range (min and max times). Since the
action of the RVL::LinearOp methods (applyOp etc.) may depend on time
and time range, RVL::TSSample must be a subclass of RVL::LinearOp with the
"set" and "get" attributes, rather than being defined through mix-in
multiple inheritance.</li>

<li>To achieve arbitrary combinations of sampling operators, create a
composition subclass for RVL::TSSample, essentially an array of
RVL::TSSamples. The time range of the composition is the convex hull of the
union of component time ranges [that is, if P1 has time range [tmin1,
tmax1] and P2 has time range [tmin2, tmax2], then the composition of
P1 and P2 has time range [min(tmin1, tmin2), max(tmax1, tmax2)]]. The
current time is synchronized between components.</li>

<li>the class RVL::TSStep realizes the time step operator \f$ H \f$ by
means of RVL::FunctionObjects encapsulating (possibly) optimized
low-level implementation of stencil operators. Because
RVL::FunctionObjects have direct access to RVL::DataContainer objects
hence to raw data, this design keeps the time step as close as
possible to the underlying low-level code yet allows a generic
implementation via virtual access functions returning references to
the necessary FOs. [in RVL3, the access functions will be replaced by
shared pointer member data] </li>

<li>RVL::TSStep is an RVL::Operator subclass, and the corresponding test
suite (DerivTest, AdjointTest) provides sufficiently direct testing of
the low-level stencil code through the RVL::FunctionObject
interfaces.</li>

<li>the time update function (\f$ h \f$ in the preceding section) is
an attribute of the RVL::TSStep - so \f$ H \f$ and \f$ h \f$ are packaged
together in RVL::TSStep. Implementation: \f$ H \f$ via RVL::LinearOp apply
methods, \f$ h \f$ via stepTimeFwd() and stepTimeBwd() methods </li>

<li>the domain (and range) space of RVL::TSStep, and range of the source
RVL::TSSample and domain of the data RVL::TSSample, are instances of RVL::TSSpace, an
RVL::ProductSpace distinguished by the type of the underlying
RVL::DataContainer it constructs (RVL::TSDC, see below). A RVL::TSSpace has two
components, representing control and state spaces respectively, and
the components are themselves RVL::TSSpaces hence ProductSpaces and may
have multiple components. An simple example is 3D acoustics, in which
the control space has bulk modulus and buoyancy components, and the
state space has pressure and three velocity components.</li>

<li>the main class is RVL::TimeStepOp, a subclass of RVL::LinOpValOp. It
representes the simulation operator \f$ F \f$ described in the
preceding section. \f$ F \f$ is a function on the control space \f$
C\f$, whose values are linear maps from the source space \f$ W \f$ to
the state space \f$ S\f$. Thus \f$ F\f$ conforms conceptually to the
RVL::LinOpValOp interface.</li>

<li>RVL::TimeStepOp has data members src (RVL::TSSample), data (RVL::TSSample), stp
(RVL::TSStep), and t (scalar).</li>

<li>The domain space of RVL::TimeStepOp is an RVL::ProductSpace with two
components: the control space, and thesource space. The range space of
RVL::TimeStepOp is the data space (which may also be a product space).</li>

<li>The apply0 method of RVL::TimeStepOp (the basic apply method, to which
all others are related) acts by creating a vector in the domain of the
RVL::TSStep data member, copying the first argument onto the control
(first) component, then running the time loop.</li>

<li>The time loop is implemented in terms of the data member
operators: slightly simplified, for x0 = control vector, x1 = source
vector, and y = data vector, the function body of apply(x0, x1, y) is

<ul>
<li>	Components<T> comp_ctrlstate(ctrlstate);</li>
<li>	comp_ctrlstate[0].copy(x0);</li>
<li>	comp_ctrlstate[1].zero();</li>
<li>	y.zero();</li>
<li>	T t = min(src.getMinTime(), data.getMinTime());</li>
<li>	T tend = data.getMaxTime();</li>
<li>	T tlast = t;</li>
<li>	while (t + TOL*(t-tlast) < tend) {
<p>
	  src.setTime(t);
<p>
	  src.applyPlusOp(x1,comp_ctrlstate[1]);
<p>
	  step.apply(ctrlstate,ctrlstate);
<p>
	  tlast = t;
<p>
	  step.stepTimeFwd(t);
<p>
	  data.setTime(t);
<p>
	  data.applyOp(comp_ctrlstate[1],y);
<p>	}</li>
</ul>
</li>

</ol>
<p>
<h2>Structure of RVL::TSStep</h2>
<p>
<ol>
<li>RVL::TSStep is an RVL::Operator with domain=range instance of RVL::TSSpace</li>

<li>The domain/range RVL::TSSpace has two component RVL::TSSpaces, representing control
and state vectors (\f$c\f$ and \f$s\f$ in the mathematical
description). Each component subspace is a RVL::TSSpace.</li>

<li>The main distinguishing characteristic of a RVL::TSSpace is that its
buildDataContainer method produces a RVL::TSDC - otherwise it is simply a
RVL::StdProductSpace, i.e. a product space with in-core
components.</li>

<li>RVL::TSDC is a subclass of RVL::StdProductDataContainer, distinguished by an implementation of eval(FO,...) that intercepts a RVL::FunctionObject of the subtype RVL::TSFO (by type checking) and provides a special unary evaluation interface, otherwise delegating to the implemented evaluation loop of the parent class. This interception permits RVL::TSFOs encoding time step functions to use the special structure of RVL::TSDCs.</li>

<li>The intercepted unary evaluation of a RVL::TSFO checks that the RVL::TSDC on which the RVL::TSFO being evaluated has two components, each RVL::TSDCs, then calls the unary evaluation interface. This check is implemented in RVL::TSDC, relieving subclasses of RVL::TSFO the necessity to make this check - they may be written on the assumption that their RVL::TSDC arguments have two components, and that the two components are RVL::TSD Cs (in particular, StdPDCs with an indexing operator)</li>

<li>The control and state components of the RVL::TSDC s belonging to vectors in the domain/range RVL::TSSpace of a RVL::TSStep are LocalDataContainers, providing access to data arrays through the getData() function. Use of the ContentPackage subclass of LDC also provides access to metadata objects for storage of grid geometry and other useful stuff. [In RVL3, LocalContentPackage]. This interface provides all data necessary to set up calls to stencil functions implemented in C, and is used in implementation of the RVL::TSFO s to which RVL::TSStep delegates the stencil computations in its apply() methods.</li>
 
</ol>
*/
