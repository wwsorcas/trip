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
#ifndef __RVL_TSSAMPLE__
#define __RVL_TSSAMPLE__

#include "GridSpace.hh"

namespace RVL {
  
  class GridCopyOverlapFO: public BinaryLocalFunctionObject<float> {
  private:
    bool plus;
  public:
    GridCopyOverlapFO(bool _plus=false): plus(_plus) {}
    GridCopyOverlapFO(GridCopyOverlapFO const &) {}
    ~GridCopyOverlapFO() {}
    void operator()(LocalDataContainer<float> & inside,
		    LocalDataContainer<float> const & outside);
    std::string getName() const { std::string tmp = "GridCopyOverlapFO"; return tmp; }
  };

  class GridExtendFO: public UnaryLocalFunctionObject<float> {
  private:
    Grid const & g;
  public:
    GridExtendFO(Grid const & _g)
      : g(_g) {}
    GridExtendFO(GridExtendFO const & f)
      : g(f.g) {}
    void operator()(LocalDataContainer<float> & x);
    std::string getName() const { std::string tmp = "GridExtendFO"; return tmp; }
  };
  
  class GridtoTSOp: public TSSample<float> {
  private:

    float tsamp;
    std::vector<size_t> ind; 
    GridSpace const & dom;
    Space<float> const & rng;
    // only first two indices count: 0 = ctrl or state, 1 = subvector index
    bool extend;
    float tol;
   
  protected:

    void apply(Vector<float> const & x,
	       Vector<float> & y) const;
    void applyPlus(Vector<float> const & x,
		   Vector<float> & y) const;
    void applyAdj(Vector<float> const & x,
		  Vector<float> & y) const;
    void applyAdjPlus(Vector<float> const & x,
		      Vector<float> & y) const;
    LinearOp<float> * clone() const { return new GridtoTSOp(*this); }

  public:

    GridtoTSOp(float _tsamp, std::vector<size_t> _ind,
	       GridSpace const & _gsp,
	       Space<float> const & _tsp,
	       bool _extend = false,
	       float _tol = 100*numeric_limits<float>::epsilon());
    GridtoTSOp(GridtoTSOp const & gop)
      : tsamp(gop.tsamp), ind(gop.ind), dom(gop.dom), rng(gop.rng), extend(gop.extend), tol(gop.tol){
      this->setMinTime(tsamp);
      this->setMaxTime(tsamp);
    }
    void setTestTime() { this->setTime(tsamp); }
    Space<float> const & getDomain() const { return dom; }
    Space<float> const & getRange() const { return rng; }
    void applyPlusOp(Vector<float> const & x,
		     Vector<float> & y) const;
    void applyPlusAdjOp(Vector<float> const & x,
			Vector<float> & y) const;
    ostream & write(ostream & str) const {
      str<<"GridtoTSOp\n";
      str<<"sample time = "<<tsamp<<"\n";
      str<<"sample tolerance = "<<tol<<"\n";
      str<<"control (0) or state (1): "<<ind[0]<<"\n";
      return str;
    }
    
  };
    
}  
    
#endif
