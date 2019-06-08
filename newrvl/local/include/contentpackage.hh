/************************************************************************

Copyright Rice University, 2006-2017
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

#ifndef __RVL_TOOLS_NEW_CONTENT_PACKAGE_
#define __RVL_TOOLS_NEW_CONTENT_PACKAGE_

#include "localdata.hh"
#include "productdata.hh"
#include "server.hh"

namespace RVL {

  // forward declaration
  template<typename T, typename M>
  class LocalContentPackage;

  /** Helper function template to extract the size of a data array
      implicit in a MetaType object. No sensible default, so throws
      exception. Specializations provided for intrinsic integer
      MetaType (size_t, int, long, uint,...). Note - returns size
      of data in WORDS - DataType not being specified cannot convert
      to bytes here. */

  template<typename MetaType>
  size_t getDataSize(MetaType const & md) {
    size_t t = md;
    return t;
  }

  /** Metadata types are not constant sized, so serialization requires
      a size function. Obvious specializations supplied for int types.
      note that return is size in BYTES.
  */

  template<class MetaType>
  size_t getMetaSize(MetaType const & md) {
    return sizeof(MetaType);
  }

  // rewrite of ADP's LocalContentPackage - now stores only a single 
  // metadata object, which is supposed to be collective
  
  // getDataSize and getMetaSize moved to rvl/productdata.hh 03.17 WWS
  
  /** works for any fixed size data, eg. a struct with default-constructed
      members of fixed size */
  template<typename MetaType>
  char * serialize(MetaType const & mt, size_t & len) {
    len = getMetaSize<MetaType>(mt);
    char * cbuf = new char[len];
    memcpy(cbuf, (char *)(&mt), len);
    return cbuf;
  }

  /** Obvious implementation for int types - will need to be specialized
      for others. Will not compile for variable size metatypes for which
      default construction is not defined. */
  template<typename MetaType>
  MetaType * deserialize(char * cbuf, size_t len) {
    // functional only for data of fixed size and 
    // default construction - test this first
    MetaType * mt = new MetaType;
    if (len < getMetaSize<MetaType>(*mt)) {
      RVLException e;
      e<<"ERROR: deserialize[LocalContentPackage MetaType]\n";
      e<<"  asserted size of char buf too small to hold\n";
      e<<"  metadata object of this type\n";
      throw e;
    }
    memcpy((char *)mt,cbuf,getMetaSize<MetaType>(*mt));
    return mt;
  }

  /** serializes CP into two contiguous blocks: (1) serialized
      metadata, followed by (2) data array.
  */
  template<typename DataType, typename MetaType>
  char * serialize(LocalContentPackage<DataType,MetaType> const & cp, size_t & len) {
    size_t mlen = 0;
    char * mbuf = serialize<MetaType>(cp.getMetadata(), mlen);
    size_t dlen = getDataSize<MetaType>(cp.getMetadata())*sizeof(DataType);
    len = mlen + dlen;
    char * cbuf = new char[len];
    memcpy(cbuf,mbuf,mlen);
    memcpy(&(cbuf[mlen]),(char *)(cp.getData()),dlen);
    delete [] mbuf;
    return cbuf;
  }

  /** Assumes that serialized CP consists of two contiguous blocks:
      (1) serialized metadata, followed by (2) data array. 
  */
  template<typename DataType, typename MetaType>
  LocalContentPackage<DataType,MetaType> * deserialize(char * cbuf, size_t len) {
    try {
      MetaType * mt = deserialize<MetaType>(cbuf,len);
      int clen = getMetaSize<MetaType>(*mt) + getDataSize<MetaType>(*mt)*sizeof(DataType);
      if (len != clen) {
	RVLException e;
	e<<"ERROR: deserialize[LocalContentPackage]\n";
	e<<"  asserted size of char buffer appears to differ from that needed\n";
	e<<"  to store CP object of given type\n";
	e<<"  asserted size of char buffer = "<<len<<"\n";
	e<<"  require buffer size for CP   = "<<clen<<"\n";
	throw e;
      }
      LocalContentPackage<DataType,MetaType> * cp = new LocalContentPackage<DataType,MetaType>;
      cp->initialize(*mt);
      delete mt;
      memcpy((char *)(cp->getData()), &(cbuf[getMetaSize<MetaType>(*mt)]),getDataSize<MetaType>(*mt)*sizeof(DataType) );
      return cp;
    }
    catch (RVLException & e) {
      e<<"\ncalled from deserialize[LocalContentPackage]\n";
      throw e;
    }
  }

  /** Helper template for implementation of standard write method. Will
      give sensible if minimal output. */
  template<class MetaType>
  ostream & writeMeta(MetaType const & md, ostream & e) { 
    e<<md<<"\n"; return e; }

  /** A normal pairing of an array of data along with Metadata
      object. 

      The only requirement on DataType: 

      * must have static size - and object of type DataType must have
      * a definite size known at compile time - and possess a default
      * constructor which allocates a memory block of constant,
      * definite size - sizeof(DataType) must make sense, and actually
      * give the size in bytes of a word of type DataType. Essentially
      * the only possibilities are an intrinsic type or a struct made
      * out of intrinsic types or a struct made out of such structs
      * or...

      The only requirements on MetaType: 

      * a deep copy constructor; 

      * correct behaviour of the getDataSize<MetaType>() helper
      * function, if necessary by template specialization, returning
      * the size of the DataType array specified by a MetaType object.

      These assumptions are implicit in the presumption that all
      information necessary to initialize and define the behaviour of
      an object of this type is present in the MetaType data member
      and the DataType type specification..

      This class follows the "initialization separate from
      instantiation" paradigm. The added flexibility is useful eg. in
      distributed computing - a recipient object can be constructed
      before the data it is to receive has been transmitted. Since the
      MetaType data member determines the size of the DataType array,
      the alternative paradigm would require some kind of auxiliary
      factory object.
  */
  
  template<class DataType, class MetaType >
  class LocalContentPackage: public LocalDataContainer<DataType> {
    
  private:

    std::shared_ptr<Server<DataType,MetaType> > fptr;
    MetaType md;
    DataType * d;
    mutable bool readonly;

    LocalContentPackage();
    LocalContentPackage(const LocalContentPackage<DataType,MetaType> & bin);
    
  public:    

    /** copy constructor - deep copy, assumes const word length for
	DataType returned by sizeof. */
    // should there be a copy constructor? fptr copied?
    /*
    LocalContentPackage(const LocalContentPackage<DataType,MetaType> & bin)
      : fptr(bin.fptr), md(bin.md), d(NULL), readonly(true) {
      d=(DataType *)usermalloc_(getDataSize<MetaType>(md)
				*sizeof(DataType));
      memcpy(d,bin.d,getDataSize<MetaType>(md)*sizeof(DataType));
    } 
    */
    /** main constructor */
    LocalContentPackage(std::shared_ptr<Server<DataType,MetaType> > _fptr(),
			bool initzero=false) 
      : fptr(_fptr), md(fptr->getMetaData()), d(NULL), readonly(true) {
      d=(DataType *)usermalloc_(getDataSize<MetaType>(md)
				*sizeof(DataType));
      if (initzero) memset(d,0,getDataSize<MetaType>(md)
			   *sizeof(DataType));
      else fptr->load(d);
    }

    virtual ~LocalContentPackage() {
      if (!readonly) fptr->save(d);
      userfree_(d);
    }
    
    /** access metadata, mutable */
    MetaType & getMetadata() { return md; }

    /** access metadata, const */
    MetaType const & getMetadata() const { return md; }

    /** access data array size */
    size_t getSize() const { 
      return getDataSize<MetaType>(md); 
    }

    /** access data, mutable */
    DataType * getData() { readonly=false; return d; } 
    /** access data, const */
    DataType const * getData() const { return d; } 

    ostream & write(ostream & str) const {
      str<<"LocalContentPackage LDC, MetaData = \n";
      writeMeta<MetaType>(md,str);
      str<<"  and data array\n";
      for (size_t i=0;i<this->getSize(); i++) {
	str<<"  d["<<i<<"]="<<d[i]<<"\n";
      }
      return str;
    }
  };

  /** FO template for transfer of connect info - intercepted by
      ContentPackage::eval(FO) */
  template<typename C>
  class AssignConnectInfo: public FunctionObject {
  private:
    C ci;
  public:
    AssignConnectInfo(C _ci): ci(_ci) {}
    AssignConnectInfo(AssignConnectInfo<C> const & aci): ci(aci.ci) {}
    ~AssignConnectInfo() {}
    C get() const { return ci; }
    std::string getName() const {
      std::string tmp="AssignConnectInfo"; return tmp;
    }
  };

  /* product DC - factors are LCPs, initialized via data servers */
  template<typename T, typename C, typename M>
  class ContentPackage: public ProductDataContainer {
  private:

    std::shared_ptr<ServerFactory<T,C,M> const> fptr;
    mutable std::vector<std::shared_ptr<Server<T,M> > > dptrs;
    mutable std::vector<std::weak_ptr<DataContainer> > wptrs;
    bool verbose;
    // can't think of any sensible copy construction
    
    ContentPackage(ContentPackage<T,C,M> const &);

  public:
    
    ContentPackage(std::shared_ptr<ServerFactory<T,C,M> const> _fptr,
		   bool _verbose = false) 
      : fptr(_fptr), dptrs(0), wptrs(0), verbose(_verbose) {}
      
    ~ContentPackage() {}

    // inherited from product
    size_t getSize() const { return fptr->getSize(); }

    // expose connect info
    C getConnectInfo() const { return fptr->getConnectInfo(); }
	
    // override PDC::eval to intercept connect info assignment
    void eval(FunctionObject & f,
	      std::vector<std::shared_ptr<DataContainer const> > x) {
      try {

	AssignConnectInfo<C> * aci = NULL;
	if (aci=dynamic_cast<AssignConnectInfo<C> *>(&f)) {
	  fptr->setConnectInfo(aci->get());
	}
	else {
	  this->ProductDataContainer::eval(f,x);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from ContentPackage::eval(FO)\n";
	throw e;
      }
    }

    // FOCE eval inherited from PDC
    /*
    void eval(FunctionObjectConstEval & f,
	      std::vector<std::shared_ptr<DataContainer const> > x) const {
      try {
	this->ProductDataContainer::eval(f,x);
      }
      catch (RVLException & e) {
	e<<"\ncalled from ContentPackage::eval(FOCE)\n";
	throw e;
      }
    }
    */
    
    std::shared_ptr<DataContainer> operator[](size_t i) {
      try {
	// initialize data server vector if not already initialized
	if (dptrs.size()==0) {
	  dptrs = fptr->build();
	  wptrs.resize(dptrs.size());
	}
	shared_ptr<DataContainer> spt;
	if (spt = (wptrs.at(i)).lock()) {
	  if (verbose) {
	    cerr<<"ContentPackage: component "<<i<<" exists ref count = "
		<<(wptrs.at(i)).use_count()<<"\n";
	  }
	}
	else {
	  if (verbose) {
	    cerr<<"ContentPackage: component "<<i<<" new\n";
	  }
	  spt = make_shared<LocalContentPackage<T,M> >(dptrs[i]);
	  wptrs.at(i) = spt;
	}
	return spt;
      }
      catch (out_of_range) {
	RVL::RVLException e;
	e<<"Error: ContentPackage::operator[]\n";
	e<<"  index "<<i<<" out of range\n";
	throw e;
      }
      catch (RVLException & e) {
	e<<"\ncalled from ContentPackage::operator[]\n";
	throw e;
      }
    }
    
    std::shared_ptr<DataContainer const> operator[](size_t i) const {
      try {
	// initialize data server vector if not already initialized
	if (dptrs.size()==0) {
	  dptrs = fptr->build();
	  wptrs.resize(dptrs.size());
	}
	
	shared_ptr<DataContainer> spt;
	if (spt = (wptrs.at(i)).lock()) {
	  if (verbose) {
	    cerr<<"ContentPackage: component "<<i<<" exists ref count = "
		<<(wptrs.at(i)).use_count()<<"\n";
	  }
	}
	else {
	  if (verbose) {
	    cerr<<"ContentPackage: const component "<<i<<" new\n";
	  }
	  spt = make_shared<LocalContentPackage<T,M> >(dptrs[i]);
	  wptrs.at(i) = spt;
	}
	return spt;
      }
      catch (out_of_range) {
	RVL::RVLException e;
	e<<"Error: ContentPackage::operator[] const\n";
	e<<"  index "<<i<<" out of range\n";
	throw e;
      }
      catch (RVLException & e) {
	e<<"\ncalled from ContentPackage::operator[] const\n";
	throw e;
      }      
    }
    
    ostream & write(ostream & str) const {
      str<<"ContentPackage PDC with LocalContentPackage<float,size_t> chunks\n";
      str<<"  data server:\n";
      fptr->write(str);
      return str;
    }
    
  };

  void readMetadata(std::vector<size_t> lens,
		    std::string fname,
		    bool verbose=false) {

    FILE * fp=NULL;
    if (!(fp=fopen(fname.c_str(),"r"))) {
      RVLException e;
      e<<"Error: readMetadata\n";
      e<<"  failed to open archival stream - error from fopen\n";
      throw e;
    }
    if (verbose) cerr<<"readMetadata: file "<<fname<<" opened\n";

    // next, determine file length
    fseeko(fp,0L,SEEK_END);
    size_t flen=ftello(fp);
    if (verbose) cerr<<"file length = "<<flen<<"\n";
    fseeko(fp,0L,SEEK_SET);
    size_t fpos=0;
    size_t ichunk=0;

    // read through file, reading metadata objects as they are
    // encountered
    while (fpos < flen) {
      if (verbose) cerr<<"read nchunk\n";
      size_t nchunk;
      fread(&nchunk,sizeof(size_t),1,fp);
      if (verbose) cerr<<"ichunk="<<ichunk<<" nchunk="<<nchunk<<"\n";
      if (nchunk <= 0) {
	RVLException e;
	e<<"Error: readMetadata\n";
	e<<"  chunk size <= 0\n";
	throw e;
      }
      // if lens nontrivial, then chunk lengths MUST match
      if (lens.size() > ichunk) {
	if (nchunk != lens[ichunk]) {
	  RVLException e;
	  e<<"Error: readMetadata\n";
	  e<<"  archival file, chunk lengths specified\n";
	  e<<"  chunk="<<ichunk<<" nchunk="
	   <<nchunk<<" lens="<<lens[ichunk]<<"\n";
	  throw e;
	}
      }
      // else store it
      else {
	lens.push_back(nchunk);
      }
      if (verbose) cerr<<"chunk length = "<<nchunk<<" floats\n";
      fpos += sizeof(size_t) + nchunk*sizeof(float);
      
      fseeko(fp,fpos,SEEK_SET);
      ichunk++;
    }
    if (fpos != flen) {
      RVLException e;
      e<<"Error: readMetadata\n";
      e<<"  failed to reach EOF: flen="<<flen<<" fpos="<<fpos<<"\n";
      throw e;
    }
    fclose(fp);
  }
  
  template<typename T, typename C, typename M>
  class ContentPackageFactory: public DataContainerFactory {

  private:
    bool verbose;
    ContentPackageFactory();
    ContentPackageFactory(const ContentPackageFactory<T,C,M> &);

  public:
    ContentPackageFactory(bool _verbose=false)
      : verbose(_verbose) {}
    
    ~ContentPackageFactory() {}

    size_t getSize() const 
 
    virtual bool compare( DataContainerFactory const & dcf) const {
      ContentPackageFactory<T,C,M> const * f = NULL;
      f = dynamic_cast<ContentPackageFactory<T,C,M> const *>(&dcf);
      if (f) {
	if (this->lens.size() == f->lens.size()) {
	  for (size_t i=0; i<this->lens.size();i++) {
	    if (this->lens[i] != f->lens[i]) return false;
	  }
	  return true;
	}
      }
      return false;
    }
    
    virtual bool isCompatible(DataContainer const & dc) const {
      ContentPackage const * g = NULL;
      g = dynamic_cast<ContentPackage const *>(&dc);
      if (g) {
	if (this->lens.size() == g->getMetadata().size()) {
	  for (size_t i=0; i<this->lens.size();i++) {
	    if (this->lens[i] != (g->getMetadata())[i]) return false;
	  }
	  return true;
	}
      }
      return false;
    }

    std::shared_ptr<DataContainer> build() const {
      if (verbose) cerr<<"ContentPackageFactory::build\n";
      return make_shared<ContentPackage>(lens,verbose);
    }

    ostream & write(ostream & str) const {
      str<<"ContentPackageFactory with metadata\n";
      for (size_t i=0;i<lens.size(); i++) 
	str<<"  lens["<<i<<"] = "<<lens[i]<<"\n";
      return str;
    }

  };

  class ContentPackageSpace: public ProductSpace<float> {

  private:

    ContentPackageFactory dcf;

    std::vector<std::shared_ptr<ocldcSpace> > sp;
    
    ContentPackageSpace();
    ContentPackageSpace(ContentPackageSpace const &);
    
  public:

    ContentPackageSpace(std::vector<size_t> lens, bool verbose=false)
      : dcf(lens,verbose) {
      sp.resize(0);
      for (size_t i=0; i<lens.size(); i++) 
	sp.push_back(std::make_shared<ocldcSpace>(lens[i]));
    }

    std::shared_ptr<DataContainer> 
    buildDataContainer() const { return dcf.build(); }

    std::shared_ptr<Space<float> const> operator[](size_t i) const {
      try {
	return sp.at(i);
      }
      catch (out_of_range) {
	RVLException e;
	e<<"Error: ContentPackageSpace::operator[]\n";
	e<<"  index = "<<i<<" out of range [0,"<<dcf.getSize()<<"]\n";
	throw e;
      }
    }
    
    virtual size_t getSize() const { return dcf.getSize(); }
    
    ostream & write(ostream & str) const {
      str<<"ContentPackageSpace defined by DataContainerFactory\n";
      dcf.write(str);
      return str;
    }
  };

}
#endif
