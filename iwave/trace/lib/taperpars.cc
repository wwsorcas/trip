#include "taperpars.hh"

TaperPars::TaperPars(std::string arg) { 

  char sep=';'; // axis data separator
  size_t pos=0;
  //  size_t posc=0;
  std::string arg1 = ""; // axis data
  int ierr;

  while ((arg1=findit<std::string>(arg,sep,pos,arg1,ierr)) != "") {
    //    cerr<<"pos = "<<posc<<" arg1 = "<<arg1<<" ierr="<<ierr<<"\n";
    //    posc=pos+1;
    char sep1=':'; // key-corner separator
    size_t pos1=0;
    size_t pidx=0;
    std::string arg2 = ""; // key, corner data
    TAPERDATA * td = new TAPERDATA; // to be pushed onto tdata
    while ((arg2=findit<std::string>(arg1,sep1,pos1,arg2,ierr)) != "") {
      //      cerr<<"pos1 = "<<posc1<<" arg2 = "<<arg2<<" ierr="<<ierr<<"\n";
      if (pidx==0) {
	if (arg2 != "") td->key = arg2;
	else {
	  RVL::RVLException e;
	  e<<"ERROR: read_taperdata\n";
	  e<<"  failed to read axis key\n";
	  throw e;
	}
      }
      else if (pidx==1) {
	if (arg2=="") {
	  RVL::RVLException e;
	  e<<"ERROR: read_taperdata\n";
	  e<<"  failed to read corner coordinates string\n";
	  throw e;
	}
	else {
	  char sep2=',';
	  size_t pos2=0;
	  std::string arg3="";
	  int tidx =  0;
	  while ((tidx < NCORNER) &&
		 ((arg3=findit<float>(arg2,sep2,pos2,td->t[tidx],ierr)) != "")) tidx++;
	  if (tidx != NCORNER) {
	    RVL::RVLException e;
	    e<<"ERROR: read_taperdata\n";
	    e<<"  failed to read "<<NCORNER<<" corner coordinates\n";
	    throw e;
	  }
	}
      }
      else {
	RVL::RVLException e;
	e<<"ERROR: read_taperdata\n";
	e<<"  too many strings separated by "<<sep1<<"\n";
	e<<"  should be precisely 2 - key and coordinates\n";
	throw e;
      }
      pidx++;
    }
    tdata.push_back(td);
  }
}

TaperPars::TaperPars(TaperPars const & t) {
  for (int i=0; i< t.tdata.size(); i++) {
    TAPERDATA * td = new TAPERDATA;
    td->key = t.tdata[i]->key;
    for (int j=0; j<NCORNER; j++) td->t[j] = t.tdata[i]->t[j];
    this->tdata.push_back(td);
  }
}

TaperPars::~TaperPars() {
  for (int i=0; i<this->tdata.size(); i++) delete tdata[i];
}

TAPERDATA const & TaperPars::operator[](size_t i) const {
  if (i < tdata.size()) return *(tdata[i]);
  else {
    RVL::RVLException e;
    e<<"ERROR: TaperPars::operator[]\n";
    e<<"  index "<<i<<" out of bounds [0,"<<tdata.size()<<"\n";
    throw e;
  }
}
ostream & TaperPars::write(ostream & str) const {
  str<<"TaperPars: object specifying taper data\n";
  for (int i=0; i< tdata.size(); i++) {
    str<<"TaperPars: axis key = "<<tdata[i]->key<<" corners = ";
    for (int k=0;k<NCORNER-1;k++) str<<(tdata[i]->t)[k]<<",";
    str<<(tdata[i]->t)[NCORNER-1];
    str<<"\n";
  }
  return str;
}


