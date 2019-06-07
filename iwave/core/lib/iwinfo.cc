#include "iwinfo.hh"

//#define IWAVE_VERBOSE

size_t pow2(int n) {
  size_t m = 1;
  while (n>0) {
    m *=2;
    n--;
  }
  return m;
}

int IWaveInfo::get_num_fields() const {
  int num=0;
  while ((get_iwave_fields()[num].field != "") && num<IWAVEMAXDATA) num++;
  if (num >= IWAVEMAXDATA) {
    RVLException e;
    e<<"Error: get_num_fields\n";
    e<<"  over limit for number - probably left off last entry with field=\"\"\n";
    throw e;
  }
  return num;
}
  
int IWaveInfo::get_num_iokeys() const {
  int num=0;
  while ((get_iwave_iokeys()[num].keyword != "") && num<IWAVEMAXDATA) num++;
  if (num >= IWAVEMAXDATA) {
    RVLException e;
    e<<"Error: get_num_iokeys\n";
    e<<"  over limit for number - probably left off last entry with keyword=\"\"\n";
    throw e;
  }
  return num;
}

// NOTE: legit return values are non-negative
int IWaveInfo::get_property_iokeys(std::string keyword,
				   std::string property) const {
  try {
    for (int i=0; i<get_num_iokeys(); i++) {
      if (keyword == get_iwave_iokeys()[i].keyword) {
	if (property == "index") return i;
	else if (property == "rarrindex") return get_iwave_iokeys()[i].rarrindex;
	else if (property == "input") return get_iwave_iokeys()[i].input;
	else if (property == "active") return get_iwave_iokeys()[i].active;
	else {
	  RVLException e;
	  e<<"Error: get_property_iokeys\n";
	  e<<"  only properties by keyword are index, rarrindex, input, and active,\n";
	  e<<"  but you asked for "<<property<<"\n";
	  throw e;
	}
      }
    }
    /*
    RVLException e;
    e<<"Error: get_property_iokeys\n";
    e<<"  input keyword = "<<keyword<<" not present in this iokeys struct\n";
    e<<"  model name = "<<iwave_model<<"\n";
    e<<"  available iokey keywords = \n";
    for (int i=0; i<get_num_iokeys(); i++) {
      e<<"    "<<get_iwave_iokeys()[i].keyword<<"\n";
    }
    throw e;
    */
    return -1;
  }
  catch (RVLException & e) {
    e<<"\ncalled from IWaveInvo::get_property_iokeys\n";
    throw e;
  }
}
    
ostream & IWaveInfo::write_iwave_fields(ostream & str) const {
  str <<"Field Definition: name = "<<get_iwave_model()<<"\n";
  for (int i=0;i<get_num_fields();i++) {
    str<<"  field["<<i<<"]="<<get_iwave_fields()[i].field
       <<" dynamic="<<get_iwave_fields()[i].dynamic
       <<" substep="<<get_iwave_fields()[i].substep
       <<" gtype=[";
    for (int j=0; j<RARR_MAX_NDIM-1; j++) 
      str<<get_iwave_fields()[i].gtype[j]<<",";
    str<<get_iwave_fields()[i].gtype[RARR_MAX_NDIM-1]<<"]\n";
  }
  return str;
}
  
ostream & IWaveInfo::write_iwave_iokeys(ostream & str) const {
  str <<"IO Definition: name = "<<get_iwave_model()<<"\n";
  for (int i=0;i<get_num_iokeys();i++) {
    str<<"  keyword["<<i<<"]="<<get_iwave_iokeys()[i].keyword
       <<" index="<<get_iwave_iokeys()[i].rarrindex
       <<" input="<<get_iwave_iokeys()[i].input
       <<" active="<<get_iwave_iokeys()[i].active
       <<"\n";
  }
  return str;
}

namespace TSOpt {

  void IOTask(std::vector<TASK_RELN *> & tr,int order, bool fwd, IWaveInfo const & ic) {
    // loop over 0,...order, adj flag
#ifdef IWAVE_VERBOSE
    cerr<<"IOTask constructor input\n";
    cerr<<"order="<<order<<" fwd="<<fwd<<endl;
    for (int i=0;i<ic.get_num_iokeys();i++) {
      cerr<<"i="<<i<<" keyword="<<ic.get_iwave_iokeys()[i].keyword;
      cerr<<" input="<<ic.get_iwave_iokeys()[i].input;
      cerr<<" active="<<ic.get_iwave_iokeys()[i].active<<endl;
    }
#endif
    if (fwd) {
      for (int n = 0; n<=order; n++) {
	if (n==0) {
	  // all inputs kept for n=0, also all passive outputs, and active if order==0
	  for (int i=0;i<ic.get_num_iokeys();i++) {
	    if (ic.get_iwave_iokeys()[i].input || 
		(!(ic.get_iwave_iokeys()[i].input) && !(ic.get_iwave_iokeys()[i].active))) {
	      TASK_RELN * p = new TASK_RELN;
	      p->iwaveindex = 0;
	      p->keyword = ic.get_iwave_iokeys()[i].keyword;
	      p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	      p->input = ic.get_iwave_iokeys()[i].input;
	      tr.push_back(p);
	    }
	    if (order==0 && !(ic.get_iwave_iokeys()[i].input) && ic.get_iwave_iokeys()[i].active) {
	      TASK_RELN * p = new TASK_RELN;
	      p->iwaveindex = 0;
	      p->keyword = ic.get_iwave_iokeys()[i].keyword;
	      p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	      p->input = ic.get_iwave_iokeys()[i].input;
	      tr.push_back(p);
	    }
	  }
	}
	else if (0<n && n<order) {
	  // only active inputs; keywords altered
	  for (int i=0;i<ic.get_num_iokeys();i++) {
	    if (ic.get_iwave_iokeys()[i].input && ic.get_iwave_iokeys()[i].active) {
	      TASK_RELN * p = new TASK_RELN;
	      p->iwaveindex = pow2(n-1);
	      std::ostringstream t;
	      t<<n;
	      p->keyword = ic.get_iwave_iokeys()[i].keyword + "_d"+t.str();
	      p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	      p->input = ic.get_iwave_iokeys()[i].input;
	      tr.push_back(p);
	    }
	  }
	}
	// active output only in comp 2^order-1
	else {
	  for (int i=0;i<ic.get_num_iokeys();i++) {
	    if (ic.get_iwave_iokeys()[i].input && ic.get_iwave_iokeys()[i].active) {
	      TASK_RELN * p = new TASK_RELN;
	      p->iwaveindex = pow2(n-1);
	      std::stringstream t;
	      t<<n;
	      p->keyword = ic.get_iwave_iokeys()[i].keyword + "_d"+t.str();
	      p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	      p->input = ic.get_iwave_iokeys()[i].input;
	      tr.push_back(p);
	    }
	    if (!(ic.get_iwave_iokeys()[i].input) && ic.get_iwave_iokeys()[i].active) {
	      TASK_RELN * p = new TASK_RELN;
	      p->iwaveindex = pow2(n)-1;
	      p->keyword = ic.get_iwave_iokeys()[i].keyword;
	      p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	      p->input = ic.get_iwave_iokeys()[i].input;
	      tr.push_back(p);
	    }
	  }
	}
      }
    }
    // adjoint branch 
    else {
      // order 0 doesn't make sense for adjoint
      // 16.02.16: nah-ah - this is the source branch
      if (order==0) {
	for (int i=0;i<ic.get_num_iokeys();i++) {
	  // nonlinear input or inactive output
	  if ((ic.get_iwave_iokeys()[i].input && 
	       (ic.get_iwave_iokeys()[i].active ==1)) || 
	      (!(ic.get_iwave_iokeys()[i].input) && 
	       !(ic.get_iwave_iokeys()[i].active))) {
	    TASK_RELN * p = new TASK_RELN;
	    p->iwaveindex = 0;
	    p->keyword = ic.get_iwave_iokeys()[i].keyword;
	    p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	    p->input = ic.get_iwave_iokeys()[i].input;
	    tr.push_back(p);
	  }
	  // linear input turns into output
	  if (ic.get_iwave_iokeys()[i].input && 
	      (ic.get_iwave_iokeys()[i].active>1)) {
	    TASK_RELN * p = new TASK_RELN;
	    p->iwaveindex = 0;
	    p->keyword = ic.get_iwave_iokeys()[i].keyword;
	    p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	    p->input = 0; //ic.get_iwave_iokeys()[i].input;
	    tr.push_back(p);
	  }
	  // active output turns into input
	  if (!(ic.get_iwave_iokeys()[i].input) && 
	      ic.get_iwave_iokeys()[i].active) {
	    TASK_RELN * p = new TASK_RELN;
	    p->iwaveindex = 0;
	    p->keyword = ic.get_iwave_iokeys()[i].keyword;
	    p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
	    p->input = 1;// ic.get_iwave_iokeys()[i].input;
	    tr.push_back(p);
	  }
	}
      }
      else {
	
	for (int n = 0; n<=order; n++) {
	  if (n==0) {

	    // all inputs kept for n=0, also all passive outputs
	    for (int i=0;i<ic.get_num_iokeys();i++) {
	      if (ic.get_iwave_iokeys()[i].input || 
		  (!(ic.get_iwave_iokeys()[i].input) && 
		   !(ic.get_iwave_iokeys()[i].active))) {
		TASK_RELN * p = new TASK_RELN;
		p->iwaveindex = 0;
		p->keyword = ic.get_iwave_iokeys()[i].keyword;
		p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
		p->input = ic.get_iwave_iokeys()[i].input;
		tr.push_back(p);
	      }
	    }

	  }
	  else if (0<n && n<order) {
	    // only active inputs; keywords altered
	    for (int i=0;i<ic.get_num_iokeys();i++) {
	      if (ic.get_iwave_iokeys()[i].input && 
		  (ic.get_iwave_iokeys()[i].active==1)) {
		TASK_RELN * p = new TASK_RELN;
		p->iwaveindex = pow2(n-1);
		std::ostringstream t;
		t<<n;
		p->keyword = ic.get_iwave_iokeys()[i].keyword + "_d"+t.str();
		p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
		p->input = ic.get_iwave_iokeys()[i].input;
		tr.push_back(p);
	      }
	    }
	  }

	  // active output only in comp = order
	  else {
	    // only active inputs; keywords altered
	    for (int i=0;i<ic.get_num_iokeys();i++) {
	      if (ic.get_iwave_iokeys()[i].input && 
		  (ic.get_iwave_iokeys()[i].active==1)) {
		TASK_RELN * p = new TASK_RELN;
		p->iwaveindex = pow2(n-1);
		std::ostringstream t;
		t<<n;
		p->keyword = ic.get_iwave_iokeys()[i].keyword + "_b"+t.str();
		p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
		p->input = !ic.get_iwave_iokeys()[i].input;
		tr.push_back(p);
	      }
	      if (!(ic.get_iwave_iokeys()[i].input) && 
		  ic.get_iwave_iokeys()[i].active) {
		TASK_RELN * p = new TASK_RELN;
		p->iwaveindex = pow2(n)-1;
		p->keyword = ic.get_iwave_iokeys()[i].keyword;
		p->rarrindex = ic.get_iwave_iokeys()[i].rarrindex;
		p->input = !ic.get_iwave_iokeys()[i].input;
		tr.push_back(p);	  
	      }
	    }
	  }
	}
      }
    }
  }

  void IOTaskWriter(std::vector<TASK_RELN *> const & tr, ostream & str) {
    for (size_t i=0;i<tr.size();i++) {
      str<<"index="<<tr[i]->iwaveindex<<" keyword="<<tr[i]->keyword<<" rarrindex="<<tr[i]->rarrindex<<" input=" << tr[i]->input<<"\n";
    }
  }
}

