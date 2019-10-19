#include "gridpp_top.hh"
#include "op.hh"

using TSOpt::GridSpace;
using namespace RVL;

class mvaop: public Operator<float> {
private:

  // wavelet file input name
  string wname;

  // data file input name
  string dname;

  // model space : bvel
  GridSpace<float> const & dom;

  // data space : image
  GridSpace<float> const & rng;

  // elevation of sources and receivers
  float r_z;
  float s_z;
    
  // default construction disabled
  mvaop();
    
protected:
    
  void apply(const Vector<float> & x, 
	     Vector<float> & y) const {
    try {

      // extract const Grid reference from range (data) space -
      Grid<float> const & g = rng.getGrid();

      // the rsf/grid data is in the various axes of
      // the Grid g; assign to variables named in your
      // usual way
      int h_x_n = g.getAxis(2).n;
      float h_x_d = g.getAxis(2).d;
      float h_x_o = g.getAxis(2).o;

      // note that r_z, s_z are not properties of the
      // rsf grid of the data - so they must be stored in 
      // data of the operator class. This is wrong, as these
      // numbers are really properties of the data. To include
      // them, we will need to start using the tfile construction.

      // next, get the filenames associated with the input and 
      // output vectors. The AssignParams function does this:
      // create it, then evaluate it on the input and output
      // vectors. The effect of AssignParam evaluation is to 
      // copy the data_type and filename properties of a GridSpace
      // vector into a "data_type = filename" pair in the PARARRAY
      // referenced by the AssignParams object.
 
      // create empty param lists
      PARARRAY * xpar = ps_new();
      PARARRAY * ypar = ps_new();

      // create two AssignParams functions, make them refer
      // to the two PARARRAYs.
      AssignParams xap(*xpar,stdout);
      AssignParams yap(*ypar,stdout);

      // transfer filenames from input and output vectors to 
      // PARARRAYS. If x is a vector in a ProductSpace, i.e. 
      // has more than one component, the AssignParam function
      // gets evaluated once an each component - so you get 
      // all of the data_type = filename pairs recorded in the
      // PARARRAY.
      x.eval(xap); // xap now has both velo and refl names
      y.eval(yap);

      // extract filenames from PARARRAYs
      string vname;
      string rname;
      if (!parse<string>(*xpar,"velocity",vname) ||
	  !parse<string>(*ypar,"reflectivity",rname)) {	
	RVLException e;
	e<<"Error: mvaop::apply - failed to extract filenames\n";
	throw e;
      }

      // build command line - C++, not python!
      // use std::stringstream to translate types
      stringstream cmdss;
      cmdss << "< " << vname << " sfrtm2d adj=1 add=0 wavfile=" << wname << " data=" << dname << " r_z=" << r_z << " s_z=" << s_z << " h_x_n=" << h_x_n << " h_x_d=" << h_x_d << " h_x_o=" << h_x_o << " rho=1.0 bndz=150.0 bndz=150.0 jsnap=50 wfldout=0 verb=1 > " << rname;

      // execute - should monitor for failure
      cout<<"mvaop: executing command "<<endl;
      cout<<cmdss.str()<<endl;
      //      system(cmdss.str().c_str());
	
      // clean up
      ps_delete(&xpar);
      ps_delete(&ypar);
    }
    catch (RVLException & e) {
      e<<"\ncalled in fop::apply\n";
      throw e;
    }
  }

  void applyDeriv(const Vector<float> & x, 
		  const Vector<float> & dx,
		  Vector<float> & dy) const {
    try {

      // extract const Grid reference from range (data) space -
      Grid<float> const & g = rng.getGrid();

      // the rsf/grid data is in the various axes of
      // the Grid g; assign to variables named in your
      // usual way
      int h_x_n = g.getAxis(2).n;
      float h_x_d = g.getAxis(2).d;
      float h_x_o = g.getAxis(2).o;

      // note that r_z, s_z are not properties of the
      // rsf grid of the data - so they must be stored in 
      // data of the operator class. This is wrong, as these
      // numbers are really properties of the data. To include
      // them, we will need to start using the tfile construction.

      // next, get the filenames associated with the input and 
      // output vectors. The AssignParams function does this:
      // create it, then evaluate it on the input and output
      // vectors. The effect of AssignParam evaluation is to 
      // copy the data_type and filename properties of a GridSpace
      // vector into a "data_type = filename" pair in the PARARRAY
      // referenced by the AssignParams object.
 
      // create empty param lists
      PARARRAY * xpar = ps_new();
      PARARRAY * dxpar = ps_new();
      PARARRAY * dypar = ps_new();

      // create two AssignParams functions, make them refer
      // to the two PARARRAYs.
      AssignParams xap(*xpar,stdout);
      AssignParams dxap(*dxpar,stdout);
      AssignParams dyap(*dypar,stdout);

      // transfer filenames from input and output vectors to 
      // PARARRAYS. If x is a vector in a ProductSpace, i.e. 
      // has more than one component, the AssignParam function
      // gets evaluated once an each component - so you get 
      // all of the data_type = filename pairs recorded in the
      // PARARRAY.
      x.eval(xap); // xap now has both velo and refl names
      dx.eval(dxap);
      dy.eval(dyap);

      // extract filenames from PARARRAYs
      string vname;
      string dvname;
      string drname;
      if (!parse<string>(*xpar,"velocity",vname) ||
	  !parse<string>(*dxpar,"velocity",dvname) ||
          !parse<string>(*dypar,"reflectivity",drname)) {	
	RVLException e;
	e<<"Error: mvaop::applyderiv - failed to extract filenames\n";
	throw e;
      }

      // build command line - C++, not python!
      // use std::stringstream to translate types
      stringstream cmdss;
      cmdss << "< " << vname << " sfwemva2d adj=0 add=0 wavfile=" << wname << " data_obs=" << dname << " deltaS=" << dvname << " r_z=" << r_z << " s_z=" << s_z << " h_x_n=" << h_x_n << " h_x_d=" << h_x_d << " h_x_o=" << h_x_o << " rho=1.0 bndz=150.0 bndz=150.0 jsnap=50 wfldout=0 verb=1 > " << drname;

      // execute - should monitor for failure
      cout<<"mvaop: executing command "<<endl;
      cout<<cmdss.str()<<endl;
      system(cmdss.str().c_str());
	
      // clean up
      ps_delete(&xpar);
      ps_delete(&dxpar);
      ps_delete(&dypar);
    }
    catch (RVLException & e) {
      e<<"\ncalled in fmva::applyderiv\n";
      throw e;
    }
  }
  
  void applyAdjDeriv(const Vector<float> & x, 
		     const Vector<float> & dy,
		     Vector<float> & dx) const {
    try {

      // extract const Grid reference from range (data) space -
      Grid<float> const & g = rng.getGrid();

      // the rsf/grid data is in the various axes of
      // the Grid g; assign to variables named in your
      // usual way
      int h_x_n = g.getAxis(2).n;
      float h_x_d = g.getAxis(2).d;
      float h_x_o = g.getAxis(2).o;

      // note that r_z, s_z are not properties of the
      // rsf grid of the data - so they must be stored in 
      // data of the operator class. This is wrong, as these
      // numbers are really properties of the data. To include
      // them, we will need to start using the tfile construction.

      // next, get the filenames associated with the input and 
      // output vectors. The AssignParams function does this:
      // create it, then evaluate it on the input and output
      // vectors. The effect of AssignParam evaluation is to 
      // copy the data_type and filename properties of a GridSpace
      // vector into a "data_type = filename" pair in the PARARRAY
      // referenced by the AssignParams object.
 
      // create empty param lists
      PARARRAY * xpar = ps_new();
      PARARRAY * dxpar = ps_new();
      PARARRAY * dypar = ps_new();

      // create two AssignParams functions, make them refer
      // to the two PARARRAYs.
      AssignParams xap(*xpar,stdout);
      AssignParams dxap(*dxpar,stdout);
      AssignParams dyap(*dypar,stdout);

      // transfer filenames from input and output vectors to 
      // PARARRAYS. If x is a vector in a ProductSpace, i.e. 
      // has more than one component, the AssignParam function
      // gets evaluated once an each component - so you get 
      // all of the data_type = filename pairs recorded in the
      // PARARRAY.
      x.eval(xap); // xap now has vel name
      dx.eval(dxap);
      dy.eval(dyap);

      // extract filenames from PARARRAYs
      string vname;
      string dvname;
      string drname;
      if (!parse<string>(*xpar,"velocity",vname) ||
	  !parse<string>(*dxpar,"velocity",dvname) ||
          !parse<string>(*dypar,"reflectivity",drname)) {	
	RVLException e;
	e<<"Error: mvaop::applyadjder- failed to extract filenames\n";
	throw e;
      }

      // build command line - C++, not python!
      // use std::stringstream to translate types
      stringstream cmdss;
      cmdss << "< " << vname << " sfwemva2d adj=1 add=0 wavfile=" << wname << " data_obs=" << dname << " deltaI=" << drname << " r_z=" << r_z << " s_z=" << s_z << " h_x_n=" << h_x_n << " h_x_d=" << h_x_d << " h_x_o=" << h_x_o << " rho=1.0 bndz=150.0 bndz=150.0 jsnap=50 wfldout=0 verb=1 > " << dvname;

      // execute - should monitor for failure
      cout<<"mvaop: executing command "<<endl;
      cout<<cmdss.str()<<endl;
      system(cmdss.str().c_str());
	
      // clean up
      ps_delete(&xpar);
      ps_delete(&dxpar);
      ps_delete(&dypar);
    }
    catch (RVLException & e) {
      e<<"\ncalled in fmva::applyadjder\n";
      throw e;
    }
  }

public:

  mvaop(mvaop const & A)
  : dom(A.dom), rng(A.rng), wname(A.wname), dname(A.dname),
    s_z(A.s_z), r_z(A.r_z) {}
      
  mvaop(GridSpace<float> const & _dom,
	GridSpace<float> const & _rng,
	string _wname, string _dname,
	float _s_z, float _r_z):
        dom(_dom), rng(_rng), wname(_wname), dname(_dname),
        s_z(_s_z), r_z(_r_z) {}

  ~mvaop() {}

  // this class is considered terminal, with no overrides foreseen,
  // so clone method is not virtual
  Operator<float> * clone() const { return new mvaop(*this); }

  // access to domain, range
  const Space<float> & getDomain() const { return dom; }
  const Space<float> & getRange() const { return rng; }

  ostream & write(ostream & str) const {
    str<<"test operator for system call op design\n";
    return str;
  }

};

int main(int argc, char ** argv) {

  try {

  // strings for names, data type
  string vname = "v.rsf";
  string rname = "r.rsf";
  string dname = "data.rsf";
  string wname = "wavelet.rsf";

  // make sure files are built
//  system("./buildrtmtest.x");

  cerr<<"1\n";
  
  // create domain and range space using rsf file name and data_type keyword
  GridSpace<float> vsp(vname,"velocity");
  GridSpace<float> rsp(rname,"reflectivity");
//  GridSpace<float> dsp(dname,"data");



  // merge two grid space to form a new space
//  StdProductSpace<float> msp(vsp,rsp);

  // create class using constructor: msp and dsp are model and data space
  mvaop op(vsp,rsp,wname,dname,0.0,0.0);

  // create input, output vectors defined in the corresponding space
  Vector<float> x(vsp);
  Vector<float> y(rsp);

  // access to components of model separately:
  // cx[0] = first component = velocity
  // cx[1] = second component = reflectivity
//  Components<float> cx(x);

  // couple vectors to files
  // looks like the concrete definition of model and data vector
  AssignFilename vaf(vname);
  x.eval(vaf);
  AssignFilename raf(rname);
  y.eval(raf);

  // evaluate operator at input
  OperatorEvaluation<float> opeval(op,x);
  
  // copy output to y
  y.copy(opeval.getValue());

  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
