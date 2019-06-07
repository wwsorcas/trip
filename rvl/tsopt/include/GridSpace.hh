/*************************************************************************

Copyright Rice University, 2004-2017.
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
#ifndef __RVL_GS
#define __RVL_GS

#include "space.hh"
#include "productspace.hh"
#include "functions.hh"
#include "contentpackage.hh"
#include "locallinalg.hh"
#include "write.hh"
#include "except.hh"
#include "utils.h"
#include "parser.h"
#include "TSOp.hh"

/* tolerance parameter - use: two grid locations
   x and y are identified if abs(x-y) < TOL*numeric_limits<float>::epsilon()*d.
*/
#define TOLFAC 100.0

/* max dimension of Grid - used to hard-truncate loops */
#define GRID_MAXDIM 10

// define boundary index between external, internal extended axis indices
#define EXTINT 100

namespace RVL {

  /** Axis: Basic Metadata element for regular grids. Defines a uniformly
      sampled axis, after the fashion of SEPlib77 or RSF. Struct part of
      RVLGrid Axis class. 
    
      Innovative feature: each axis carries an integer token (axis.id)
      meant to signify its position in a global ordering of axes (hence
      its physical meaning)..

      Commensurable constraint: all axes may either subgrids of the
      integer (primal) lattice scaled by d, or the dual integer+1/2
      lattice scaled by d, or neither. This status is hecked on
      construction, and a flag set to signify lattice type: atype=0
      for primal, atype=1 for dual, atype=-1 for neither.

      For indexed operations, it is convenient to refer to a global
      latice, with origin in the base grid cell [0, d). The first and
      last grid indices with respect to this lattice are gs0 and
      ge0. It is also possible to defined an axis subinterval or
      window, with limits gs >= gs0 and ge <= ge0. These are defaulted
      to gs0 and ge0 on construction, and may be reset
      post-construction.
  */
  
  class Axis: public Writeable {

  private:
    /** number of gridpoints on axis */
    size_t n; 
    /** step between gridpoints */
    float d;
    /** coordinate of first gridpoint */
    float o;
    /** axis index */
    int id;
    /** axis type (primal or dual), computed */
    int atype;
    /** unit of d */
    std::string dunit;
    /** indices of axis start and end in global grid with 
	origin at zero vector */
    int gs0; int ge0;
    /** active window limits - default to ends of axis */
    mutable int gs;
    mutable int ge;
    
  public:
    
    Axis(size_t _n = 0,
	 float _d = 1.0f,
	 float _o = 0.0f,
	 int _id = 0,
	 std::string _dunit = "" );
    Axis(Axis const & a);
    bool operator==(Axis const & a) const;
    void operator=(Axis const & a);
    size_t getLength() const { return n; }
    float getStep() const { return d; }
    float getOrigin() const { return o; }
    int getID() const { return id; }
    std::string getUnit() const { return dunit; }
    int getAxisType() const { return atype; }
    int getAxisStartIndex() const { return gs0; }
    int getAxisEndIndex() const { return ge0; }
    int getSubAxisStartIndex() const { return gs; }
    int getSubAxisEndIndex() const { return ge; }
    void setSubAxisStartIndex(int _gs) const;
    void setSubAxisEndIndex(int _ge) const;
    // used in implementation of getMetaSize helper function for ContentPackage class
    size_t getAxisSize() const;
    ostream & write(ostream & str) const;
  };

  /** returns true if 
      1. a.d = b.d
      2. a.o - b.o is divisible by a.d
      These conditions imply that a and b are sub-axes of a global axis
  */
  bool areCompatible(Axis const & a, Axis const & b);
  
  class Grid: public Writeable {
    
  private:
    
    size_t gdim;
    mutable std::vector<Axis> axes;
    
    mutable bool init;

    mutable IPNT gs0;
    mutable IPNT gs;
    mutable IPNT ge0;
    mutable IPNT ge;
    
  public:
    
    Grid(size_t _gdim=0);
    Grid(Grid const & g);
    void addAxis(Axis const & a);
    size_t getDimension() const { return gdim; }
    bool isInit() const { return init; }
    Axis const & getAxis(size_t i) const;
    float getCellVolume() const;
    void getAxisLimits(IPNT gs0, IPNT ge0) const;
    void getSubAxisLimits(IPNT gs, IPNT ge) const;
    void setSubAxisLimits(IPNT gs, IPNT ge) const;
    bool operator==(Grid const & g) const;
    bool operator!=(Grid const & g) const { return !((*this)==g); }
    ostream & write(ostream & str) const;
  };

  /** returns true if all axes are compatible, pairwise */
  bool areCompatible(Grid const & g, Grid const & h);
  
  template<> ostream & writeMeta<Grid>(Grid const & g, ostream & e);

  template<> size_t getDataSize<Grid>(Grid const & g);

  template<> size_t getMetaSize<Grid>(Grid const & g);
  
  template<> char * serialize<Axis>(Axis const & a, size_t & len);
  
  template<> Axis * deserialize<Axis>(char * cbuf, size_t len);

  template<> char * serialize<Grid>(Grid const & g, size_t & len);
    
  template<> Grid * deserialize<Grid>(char * cbuf, size_t len);

  std::shared_ptr<Grid> make_padded_grid(Grid const & phys,
					 std::vector<int> nlsloc,
					 std::vector<int> nrsloc);

  /** divvy up an array of n T's into an array of arrays of Ts, each 
      chunklen long - like matlab reshape of vector into matrix */
  template<typename T> T** dimView(T* x, size_t n, size_t chunklen) {
    size_t nchunks = n/chunklen;
    if (n != nchunks*chunklen) {
      RVLException e;
      e<<"Error: dimView\n";
      e<<"  length of input array = "<<n<<" not divisible by chunk length = "
       <<chunklen<<"\n";
      throw e;
    }
    T ** xchunks = (T**)usermalloc_(nchunks*sizeof(T*));
    for (size_t i=0; i< nchunks; i++) xchunks[i] = &(x[i*chunklen]);
    return xchunks;
  }

  class GridDC: public ContentPackage<float,Grid> {
  private:
    mutable float * data1d;
    mutable float * data1ds;
    mutable float const * data1dsc;
    mutable float ** data2d;
    mutable float ** data2ds;
    mutable float const ** data2dsc;
    mutable float *** data3d;
    mutable float *** data3ds;
    mutable float const *** data3dsc;
    mutable float **** data4d;
    mutable float ***** data5d;
  public:
    GridDC(Grid const & g);
    GridDC(GridDC const & d);
    ~GridDC();
    Grid const & getGrid() const {
      Grid const & g = dynamic_cast<Grid const &>(this->getMetadata());
      return g;
    }
    
    float * getData1D();
    float const * getData1D() const;    
    float * getData1DGlobal();
    float const * getData1DGlobal() const;
    float ** getData2D();
    float const ** getData2D() const;
    float ** getData2DGlobal();
    float const ** getData2DGlobal() const;
    float *** getData3D();
    float const *** getData3D() const;
    float *** getData3DGlobal();
    float const *** getData3DGlobal() const ;
    float **** getData4D();
    float const **** getData4D() const;
    float ***** getData5D();
    float const ***** getData5D() const;
  };

  class GridDCFactory: public DataContainerFactory {
  private:
    Grid g;
  protected:
    GridDC * buildDC() const {
      try {
	GridDC * dc = new GridDC(g);
	return dc;
      }
      catch (RVLException & e) {
	e<<"\ncalled from GridDCFactory::buildDC\n";
	throw e;
      }
    }
    DataContainer * build() const { return buildDC(); }
  public:
    GridDCFactory(Grid const & _g): g(_g) {}
    GridDCFactory(GridDCFactory const & gf): g(gf.g) {}
    bool compare(DataContainerFactory const & dcf) const {
      try {
	GridDCFactory const & gdcf = dynamic_cast<GridDCFactory const &>(dcf);
	return (gdcf.g == g);
      }
      catch (bad_cast) {
	return false;
      }
    }
    bool isCompatible(DataContainer const & dc) const {
      try {
	GridDC const & gdc = dynamic_cast<GridDC const &>(dc);
	Grid const & gdcg = gdc.getMetadata();
	return (gdcg == g);
      }
      catch (bad_cast) {
	return false;
      }
    }

    Grid const & getGrid() const { return g; }
    
    ostream & write(ostream & str) const {
      str<<"GridDCFactory based on Grid:\n";
      g.write(str);
      return str;
    }
  };

  class GridSpace: public StdSpace<float,float> {
    
  private:
    
    GridDCFactory f; 
    RVLLinearAlgebraPackage<float> l;

  protected:

    Space<float> * clone() const { return new GridSpace(*this); }
    
  public:
    GridSpace(Grid const & g)
      : f(g), l(g.getCellVolume()) {}
    GridSpace(GridSpace const & sp)
      : f(sp.f), l(sp.l) {}
    LinearAlgebraPackage<float> const & getLAP() const { return l; }
    DataContainerFactory const & getDCF() const { return f; }
    // convenience
    Grid const & getGrid() const { return f.getGrid(); }
    ostream & write(ostream & str) const {
      str<<"GridSpace: vector space class for gridded data\n";
      f.getGrid().write(str);
      return str;
    }
  };

  class GridDomain: public TSSpace<float> {
  private:

    std::map<std::string,int> ind;

    size_t gdim;
    std::vector<float> dxs;
    mutable IPNT gs; // workspace for axis starts
    mutable IPNT ge; // workspace for axis ends
    
  public:

    GridDomain(std::vector<Grid *> ctrlgrids,
	       std::vector<Grid *> stategrids,
	       std::map<std::string,int> _ind);
    GridDomain(GridDomain const & gd);

    void get_ga(IPNT gsa, IPNT gea, int i) const;
    void get_gc(IPNT gsc, IPNT gec, int i) const;
    std::map<std::string,int> get_ind() const { return ind; }
    size_t getDimension() const { return gdim; }
    std::vector<float> getSteps() const { return dxs; }
    
    ostream & write(ostream & str) const;
  };

  /** extracts grid from RSF file format */
  std::shared_ptr<Grid> GridfromFile(std::string fname);

  /** check that calling object is GridDC with same grid as
      defined by file, then read or write data according to 
      value of load flag */
  class GridDCIOFO: public UnaryLocalFunctionObject<float> {
  private:
    std::string fname;
    bool load;
    shared_ptr<Grid> g;
  public:
    GridDCIOFO(std::string _fname, bool _load=true, bool _extend=false)
      : fname(_fname), load(_load), g(GridfromFile(fname)) {}
    GridDCIOFO(GridDCIOFO const & f)
      : fname(f.fname), load(f.load), g(f.g) {}
    void operator()(LocalDataContainer<float> & x);
    std::string getName() const { std::string tmp="GridDCIOFO"; return tmp; }
  };
}

#endif

