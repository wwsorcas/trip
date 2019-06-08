/************************************************************************

Copyright Rice University, 2017
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

#ifndef __RVL_TOOLS_SIMPLE_FILE_SERVER_
#define __RVL_TOOLS_SIMPLE_FILE_SERVER_

#include "server.hh"

namespace RVL {
  
  template<typename T>
  class SimpleFileServer: public Server<T,size_t> {
  private:
    size_t n;
    FILE * fp;
    size_t offset;
    bool verbose;
  public:
    SimpleFileServer(size_t _n, FILE * _fp, size_t _offset, bool _verbose)
      : n(_n), fp(_fp), offset(_offset), verbose(_verbose) {
      // sanity check
      if (!fp) {
	RVLException e;
	e<<"Error: SimpleServer constructor\n";
	e<<"  fp=NULL\n";
	throw e;
      }
      fseeko(fp,0L,SEEK_END);
      size_t flen=ftello(fp);
      if (verbose) cerr<<"file length = "<<flen<<"\n";
      if (offset + n*sizeof(T) - 1 > flen) {
	RVLException e;
	e<<"Error: SimpleServer constructor\n";
	e<<"  offset = "<<offset<<" n = "<<n
	 <<" sum greater than file length = "<<flen<<"\n";
	throw e;
      }
      size_t nn;
      fseeko(fp,offset,SEEK_SET);
      fread(&nn,sizeof(size_t),1,fp);
      if (nn!=n) {
        RVLException e;
	e<<"Error: SimpleFileServer constructor\n";
	e<<"  file data incompatible with spec array len\n";
	e<<"  for chunk at offset "<<offset<<"\n";
	e<<"  chunk length from constructor arg = "<<_n<<"\n";
	e<<"  chunk length from file = "<<n<<"\n";
	throw e;
      }
    }
    SimpleFileServer(SimpleFileServer<T> const & sfm)
      : n(sfm.n), fp(sfm.fp), offset(sfm.offset), verbose(sfm.verbose) {}
    ~SimpleFileServer() {}
    void load(T * d) const {
      fseeko(fp,offset+sizeof(size_t),SEEK_SET);
      if (n != fread(d,sizeof(T),n,fp)) {
        RVLException e;
	e<<"Error: SimpleFileServer::load\n";
	e<<"  failed to read "<<n<<" floats at offset "<<offset<<"\n";
	throw e;
      }
    }
    void save(T const * d) const {
      fseeko(fp,offset+sizeof(size_t),SEEK_SET);
      if (n != fwrite(d,sizeof(T),n,fp)) {
        RVLException e;
	e<<"Error: SimpleFileServer::load\n";
	e<<"  failed to read "<<n<<" floats at offset "<<offset<<"\n";
	throw e;
      }
    }
    // expose metadata for CPs
    size_t getMetaData() { return n; }
  };

  template<typename T>
  SimpleFileServerFactory<std::string, size_t>
    : public ServerFactory<T,std::string,size_t> {
  private:
    std::vector<size_t> lens;
    mutable std::string fname;
    mutable FILE * fp;
    mutable std::vector<size_t> offsets;

    void writeFilefromMetadata() {
      if (!fp) {
	RVLException e;
	e<<"ERROR: SimpleFileServerFactory::writeFilefromMetadata\n";
	e<<"  file not open yet\n";
	throw e;
      }
      //    cerr<<"writeFilefromMetadata\n";
      // fill up file writing chunks
      fseeko(fp,0L,SEEK_SET);
      size_t len;
      size_t lenmax=0;
      for (size_t i=0; i<lens.size(); i++) {
	lenmax=max(lenmax,lens[i]);
      }
      float * x = (float *)malloc(lenmax*sizeof(float));
      memset((char *)x,0,lenmax*sizeof(float));
      size_t fpos =0;
      
      for (size_t i=0; i<lens.size(); i++) {
	len=lens[i];
	//if (verbose) cerr<<"write chunk length = "<<len<<"\n";
	fwrite(&len,sizeof(size_t),1,fp);
	//if (verbose) cerr<<"write "<<lens[i]<<" floats\n";
	fwrite(x,sizeof(float),lens[i],fp);
	fpos+=sizeof(size_t)+lens[i]*sizeof(float);
	//if (verbose) cerr<<"chunk len="<<lens[i]<<"\n";
	fflush(fp);
      }
      free(x);
    }

    // private function to initiate data connection
    void initialize() const {
      try {
	if (verbose) cerr<<"SimpleFileServerFactory::initialize: lens.size="<<lens.size()<<"\n";
	
	// branch 1: for empty filename, open temp file
	if (fname.size()==0) {
	  if (verbose) cerr<<"SimpleFileServerFactory::initialize - temp\n";
	  // sanity: if no archival file, then metadata (vector of
	  // chunk lengths, in this case) must be nontrivial
	  if (lens.size()==0) {
	    RVLException e;
	    e<<"Error: SimpleFileServerFactory::initialize\n";
	    e<<"  temp case: metadata (chunk lengths) not supplied\n";
	    e<<"  lens.size=0\n";
	    throw e;
	  }
	  //first identify DATAPATH, if it exists

	  std::string dpath;
	  std::string dpathraw = getenv("DATAPATH");

	  // create full pathname as std::string
	  if (dpathraw.size() > 0) {
	    if (dpathraw[dpathraw.size()-1] != '/') 
	      dpath=dpathraw + "/tmp.XXXXXX";
	    else
	      dpath=dpathraw + "tmp.XXXXXX";
	  }
	  else
	    dpath="./tmp.XXXXXX"; 

	  // overwrite using mkstemp, get file descriptor
	  if (verbose) cerr<<"dpath="<<dpath;
	  char * dpathc = (char *)malloc((dpath.size()+1)*sizeof(char));
	  strcpy(dpathc,dpath.c_str());
	  int fd = mkstemp(dpathc);

	  // initialize fname, open file pointer
	  fname = dpathc;
	  if (verbose) cerr<<"after mkstemp: fname="<<fname<<"\n";
	  if (fd<0) {
	    RVLException e;
	    e<<"Error: SimpleFileServerFactory::initialize\n";
	    e<<"  failed to open temp file - error from mkstemp\n";
	    throw e;
	  }
	  if (!(fp=fdopen(fd,"w+"))) {
	    RVLException e;
	    e<<"Error: SimpleFileServerFactory::initialize\n";
	    e<<"  failed to open temp stream - error from fdopen\n";
	    throw e;
	  }
	  if (verbose) cerr<<"temp file opened\n";
	  if (verbose) cerr<<"initialize file from metadata\n";
	  writeFilefromMetadata(lens,fp);
	  if (verbose) cerr<<"reposition file pointer\n";
	  fseeko(fp,0L,SEEK_SET);
	}

	// branch 2: archival file name - ARCHIVAL FILES MUST EXIST
	// note: metadata presumed to be initialized in
	// any case - so if file is readable, must be checked
	// against contents
	else {

	  if (!(fp=fopen(fname.c_str(),"r+"))) {
	    RVLException e;
	    e<<"Error: SimpleFileServerFactory::initialize\n";
	    e<<"  failed to open archival stream - error from fopen\n";
	    throw e;
	  }
	  if (verbose) cerr<<"file "<<fname<<" opened\n";
	  
	  // next, check consistency of existing file
	  fseeko(fp,0L,SEEK_END);
	  size_t flen=ftello(fp);
	  if (verbose) cerr<<"file length = "<<flen<<"\n";
	  fseeko(fp,0L,SEEK_SET);
	  size_t fpos=0;
	  size_t ichunk=0;
	  while (fpos < flen) {
	    if (verbose) cerr<<"read loop fpos="<<fpos<<" ichunk="<<ichunk<<" lens.size="<<lens.size()<<"\n";
	    if (ichunk > lens.size()) {
	      RVLException e;
	      e<<"Error: SimpleFileServerFactory::initialize\n";
	      e<<"  attempt to access index out of bounds [0,"
	       <<lens.size()<<"]\n";
	      throw e;
	    }
	    size_t nchunk;
	    fread(&nchunk,sizeof(size_t),1,fp);
	    if (nchunk <= 0) {
	      RVLException e;
	      e<<"Error: SimpleFileServerFactory::initialize\n";
	      e<<"  chunk size <= 0\n";
	      throw e;
	    }
	    // if lens nontrivial, then chunk lengths MUST match

	    if (nchunk != lens[ichunk]) {
	      RVLException e;
	      e<<"Error: SimpleFileServerFactory::initialize\n";
	      e<<"  archival file, chunk lengths specified\n";
	      e<<"  chunk="<<ichunk<<" nchunk="<<nchunk<<" lens="<<lens[ichunk]<<"\n";
	      throw e;
	    }
	    if (verbose) cerr<<"chunk length = "<<nchunk<<" floats\n";
	    fpos += sizeof(size_t) + nchunk*sizeof(float);
	    fseeko(fp,fpos,SEEK_SET);
	    ichunk++;
	  }
	  if (fpos != flen) {
	    RVLException e;
	    e<<"Error: SimpleFileServerFactory::initialize\n";
	    e<<"  failed to reach EOF: flen="<<flen
	     <<" fpos="<<fpos<<"\n";
	    throw e;
	  }  
	}
	
      }
      catch (RVLException & e) {
	e<<"\ncalled from SimpleFileServerFactory::initialize\n";
	throw e;
      }     
    }

  public:
    // constructor should get CP Metadata; connection to data
    // source deferred
    SimpleFileServerFactory(std::vector<size_t> _lens,
			    bool _verbose=false)
      : lens(_lens), fp(NULL), fname(""), verbose(_verbose) {
      if (verbose) cerr<<"SimpleFileServerFactory constructor\n";
            // initialize offset vector
      offsets.push_back(0);
      for (size_t i=1; i<lens.size();i++) {
	offsets.push_back(offsets[i-1]+sizeof(size_t) +
			  lens[i-1]*sizeof(float));
      }
    }

    // supply address for data source
    void setConnectInfo(std::string _fname) {
      if ((fp) || (fname.size() > 0)) {
	RVLException e;
	e<<"ERROR: SimpleFileServerFactory::setConnectInfo\n";
	e<<"  cannot set filename to "<<_fname<<"\n";
	e<<"  because it is already set to "<<fname<<"\n";
	throw e;
      }
      fname = _fname;
    }

    // connect and create Server vector
    virtual std::vector<std::shared_ptr<Server<T,size_t> > > build() const {
      try {
	// if not connected but nontrivial name supplied, attempt connect
	if (!fp && fname.size() > 0) initialize();
	// if connected, create and return Server vector
	if (fp) {
	  std::vector<std::shared_ptr<Server<T,size_t> > > retval;
	  for (int i=0; i<lens.size(); i++) {
	    std::shared_ptr<Server<T> > tmp
	      = make_shared<SimpleFileServer<T>(lens[i],fp,offsets[i], verbose);
	    retval.push_back(tmp);
	  }
	  return retval;
	}
	else {
	  RVLException e;
	  e<<"ERROR: SimpleFileServerFactory::build\n";
	  e<<"  source file not opened\n";
	  e<<"  fp = "<<fp<<" fname = "<<fname<<"\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from SimpleFileServerFactory::build\n";
	throw e;
      }
    }
  };

}

#endif

  
