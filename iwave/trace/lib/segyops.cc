#include "segyops.hh"

// precondition for all of these FOs: input, output have same geom

namespace TSOpt {
    
  using RVL::LocalDataContainer;
  using RVL::RVLException;

  void uconv(int outbeg, int outend, float * restrict out,
	     int inbeg, int inend, const float * restrict in,
	     int kerbeg, int kerend, const float * restrict ker,
	     float scal) {
    for (int i=outbeg; i<= outend; i++) {
      int start = iwave_max(inbeg, i-kerend);
      int stop = iwave_min(inend, i-kerbeg);
      out[i]=0.0f;
      for (int j=start; j<=stop; j++) {
	out[i] += in[j]*ker[i-j];
      }
    }
  }
  
  void conv(int shift, int nout, int nin, int nker,
	    float * restrict out,
	    const float * restrict in,
	    const float * restrict ker,
	    float scal ) {
    int start;
    int stop;
    for (int i=0; i<nout; i++) {
      start = iwave_max(i+shift-nker+1,0);
      stop  = iwave_min(i+shift+1,nin);
      out[i]=0.0f;

      //#pragma ivdep
      for (int j=start; j<stop; j++) {
	out[i] += ker[i-j+shift]*in[j];
      }
      out[i] *= scal;
    }
  }

  void corr(int shift, int nout, int nin, int nker,
	    float * restrict out,
	    const float * restrict in,
	    const float * restrict ker,
	    float scal ) {
    int start;
    int stop;
    for (int i=0; i<nout; i++) {
      start = iwave_max(i-shift,0);
      stop  = iwave_min(i-shift+nker,nin);
      out[i]=0.0f;
      //#pragma ivdep
      for (int j=start; j<stop; j++) {
	out[i] += ker[j-i+shift]*in[j];
      }
      out[i] *= scal;
    }
  }


  /* in physical discrete coordinates with sample rate dt, convolution
     computes out(ip), outbeg <= ip <= outend, as

     sum from jp = max(ip-kerend, inbeg) to min(ip-kerbeg, inend) 
     of ker(ip-jp)in(jp) = ker(kp)in(jp), 

     array indices: ip=i+outbeg, jp=j+inbeg, kp=k+kerbeg=i-j+outbeg-inbeg
     so kernel index = k = i-j+outbeg-inbeg-kerbeg
   
     so 
     max(i+outbeg-kerend,inbeg) <= j+inbeg <= min(i+outbeg-kerbeg,inend)
     max(i+oubbeg-kerend,inbeg)-inbeg <= j <= min(i+outbeg-kerbeg,inend)-inbeg 

     0 <= i <= outend-outbeg

     check: k=i-j+outbeg-inbeg-kerbeg >= i - min(i+outbeg-kerbeg,inend)+inbeg+outbeg-inbeg-kerbeg
     >= i+outbeg-kerbeg - min(i+outbeg-kerbeg,inend) >=0 

     sim k <= kerend-kerbeg

     in terms of start and length
     nout = outend-outbeg+1  : outend = nout+outbeg-1
     nin  = inend-inbeg+1    : inend  = nin +inbeg -1
     nker = kerend-kerbeg+1  : kerend = nker+kerbeg-1

     so: 0 <= i < nout
     max(i+outbeg-nker-kerbeg+1,inbeg)-inbeg <= j < min(i+outbeg-kerbeg,nin+inbeg-1)-inbeg+1
     max(i+outbeg-kerbeg-inbeg+1-nker,0)     <= j < min(i+outbeg-kerbeg-inbeg+1,nin)
  */

  void SEGYConvolve::operator()(LocalDataContainer<float> & x,
				LocalDataContainer<float> const & y,
				LocalDataContainer<float> const & z) {
    try {
      segytrace & cpout = dynamic_cast<segytrace &>(x);
      segy & trout = cpout.getMetadata();
      segytrace const & cpin = dynamic_cast<segytrace const &>(y);
      segy const & trin = cpin.getMetadata();
      segytrace const & cpker = dynamic_cast<segytrace const &>(z);
      segy const & trker = cpker.getMetadata();

      Value val;
      float dt;
      int nin;
      int nout;
      int nker;
      
      int inbeg;
      int outbeg;
      int kerbeg;
      
      string dtstr="dt";
      string nsstr="ns";
      string rtstr="delrt";
      float tmp;

      gethdval(&trout,  (char*)(dtstr.c_str()), &val);
      dt=vtof(hdtype((char*)(dtstr.c_str())), val);
      dt*=0.001;

      gethdval(&trin,  (char*)(nsstr.c_str()), &val);
      nin=vtoi(hdtype(nsstr.c_str()), val);
      gethdval(&trin,  (char*)(rtstr.c_str()), &val);
      tmp=vtof(hdtype(rtstr.c_str()), val);
      if (tmp<0.0f) inbeg=int(-0.01 + tmp/dt);
      else inbeg=int(0.01 + tmp/dt);
      
      gethdval(&trout,  (char*)(nsstr.c_str()), &val);
      nout=vtoi(hdtype(nsstr.c_str()), val);
      gethdval(&trout,  (char*)(rtstr.c_str()), &val);
      tmp=vtof(hdtype(rtstr.c_str()), val);
      if (tmp<0.0f) outbeg=int(-0.01 + tmp/dt);
      else outbeg=int(0.01 + tmp/dt);

      gethdval(&trker,  (char*)(nsstr.c_str()), &val);
      nker=vtoi(hdtype(nsstr.c_str()), val);
      gethdval(&trker,  (char*)(rtstr.c_str()), &val);
      tmp=vtof(hdtype(rtstr.c_str()), val);
      if (tmp<0.0f) kerbeg=int(-0.01 + tmp/dt);
      else kerbeg=int(0.01 + tmp/dt);

      if (adj==0) {
	int ishift=outbeg-inbeg-kerbeg;
	conv(ishift,nout,nin,nker,
	     trout.data,trin.data,trker.data,dt);
	//	const float * in     = &(trin.data[-inbeg]);
	//	float * out    = &(trout.data[-outbeg]);
	//	const float * ker    = &(trker.data[-kerbeg]);

	//	cerr<<"inbeg="<<inbeg<<" outbeg="<<outbeg<<" kerbeg="<<kerbeg<<"\n";
	//	cerr<<"nin="<<nin<<" nout="<<nout<<" nker="<<nker<<"\n";
	
	//	uconv(outbeg,outbeg+nout-1,out,
	//	      inbeg,inbeg+nin-1,in,
	//	      kerbeg,kerbeg+nker-1,ker,
	//	      dt);
      }
      else {
	int ishift=inbeg-kerbeg-outbeg;
	corr(ishift,nout,nin,nker,
	     trout.data,trin.data,trker.data,dt);	
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"ERROR: SEGYConvolve::operator()\n";
      e<<"  at least one input not segytrace type\n";
      throw e;
    }
  }

  void SEGYTaper::operator()(LocalDataContainer<float> & x,
			     LocalDataContainer<float> const & y) {
    try {
      segytrace & cpout = dynamic_cast<segytrace &>(x);
      segy & trout = cpout.getMetadata();
      segytrace const & cpin = dynamic_cast<segytrace const &>(y);
      segy const & trin = cpin.getMetadata();
      Value val;
      float fac=1.0f;
      int dontdoit = 1;
      // fetch coordinates for each taper axis
      float xcoord;
      
      std::string tmp = "scalco";
      gethdval(&trout,(char *)(tmp.c_str()),&val);
      float scalco = vtof(hdtype("scalco"),val);
      for (int i=0; i<tap.getSize(); i++) {
	// special case: offset, computed from source/receiver posns
	if (tap[i].key=="offset") {
	  tmp="sx";
	  gethdval(&trout,(char *)(tmp.c_str()),&val);
	  float sx = vtof(hdtype("sx"),val);	  
	  tmp="sy";
	  gethdval(&trout,(char *)(tmp.c_str()),&val);
	  float sy = vtof(hdtype("sy"),val);	  
	  tmp="gx";
	  gethdval(&trout,(char *)(tmp.c_str()),&val);
	  float gx = vtof(hdtype("gx"),val);	  
	  tmp="gy";
	  gethdval(&trout,(char *)(tmp.c_str()),&val);
	  float gy = vtof(hdtype("gy"),val);	  
	  xcoord = sqrt(((sx-gx)*(sx-gx))+((sy-gy)*(sy-gy)));
	}
	else {
	  gethdval(&trout,(char *)(tap[i].key.c_str()),&val);
	  xcoord = vtof(hdtype(tap[i].key.c_str()),val);
	}
	if (scalco > 0) { xcoord *=  scalco; }
	if (scalco < 0) { xcoord /= -scalco; }
	dontdoit*=taperfun(0,tap[i].t[0],tap[i].t[1],xcoord,fac); // returns 1 if x to r of t[1] 
	dontdoit*=taperfun(1,tap[i].t[2],tap[i].t[3],xcoord,fac); // returns 1 if x to l of t[2]
	//	cerr<<"xcoord="<<xcoord<<" t0="<<tap[i].t[0]<<" t1="<<tap[i].t[1]<<" t2="<<tap[i].t[2]<<" t3="<<tap[i].t[3]<<" dontdoit="<<dontdoit<<" fac="<<fac<<endl;
      }
      tmp="ns";
      gethdval(&trout,(char *)(tmp.c_str()),&val);
      int ns = vtoi(hdtype("ns"),val);
      if (dontdoit) {
	for (int k=0;k<ns;k++) (trout.data)[k] = (trin.data)[k];
      }
      else {
	for (int k=0;k<ns;k++) (trout.data)[k] = fac*((trin.data)[k]);
      }
    }
    catch (bad_cast) {
      RVL::RVLException e;
      e<<"ERROR: SEGYTaper::operator[]\n";
      e<<"  input not of segytrace type\n";
      throw e;
    }
  }
      
  void SEGYLinMute::operator()(LocalDataContainer<float> & x,
			       LocalDataContainer<float> const & y) {
        
    try {
            
      segytrace & cpout = dynamic_cast<segytrace &>(x);
      segytrace const & cpin = dynamic_cast<segytrace const &>(y);
            
      segy & trout = cpout.getMetadata();
      segy const & trin = cpin.getMetadata();
            
      Value val;
      float gx;
      float x;
      float dt;
      float t0;
      float wm;
      int nt;
            
      string name="ns";
      gethdval(&trin,(char *)(name.c_str()),&val);
      nt=vtoi(hdtype((char *)(name.c_str())),val);
            
      name="dt";
      gethdval(&trin,(char *)(name.c_str()),&val);
      dt=0.001*vtof(hdtype((char *)(name.c_str())),val);
            
      name="delrt";
      gethdval(&trin,(char *)(name.c_str()),&val);
      t0=vtof(hdtype((char *)(name.c_str())),val);
            
      name="gx";
      gethdval(&trin,(char *)(name.c_str()),&val);
      gx = vtof(hdtype((char *)(name.c_str())),val);
            
      name="sx";
      gethdval(&trin,(char *)(name.c_str()),&val);
      x = gx - vtof(hdtype((char *)(name.c_str())),val);
            
      // adjust mute so that w=0 works
      if (mute_type) {
	wm = MAX(dt/s,w);
	float wtmp = mutefun((gx - s)/wm) * mutefun((tm - gx)/wm);
	for (int i=0;i<nt;i++)
	  trout.data[i] = trin.data[i]*wtmp;
      }
      else {
	wm=MAX(dt,w);
	//	cerr<<"mute at offset="<<x<<" slope="<<s<<" zotime="<<tm<<" width="<<wm<<endl;
	for (int i=0;i<nt;i++)
	  //	trout.data[i] = trin.data[i]*mutefun((t0+i*dt-s*fabs(x)-tm)/wm);
	  // this version: tm = time AFTER first sample
	  trout.data[i] = trin.data[i]*mutefun((i*dt+t0-s*fabs(x)-tm)/wm);
      }
    }
    catch (bad_cast) {
      cerr<<"Warning: SEGYLinMute::operator()\n";
      cerr<<"input LDC not segytrace, therefore copy op\n";
      RVL::RVLCopy<float> cp;
      cp(x,y);
    }
    catch (RVLException & e) {
      e<<"\ncalled from SEGYLinMute::operator()\n";
      throw e;
    }
  }
    
    
  void SEGYTaperMute::operator()(LocalDataContainer<float> & x,
				 LocalDataContainer<float> const & y) {
        
    try {
            
      segytrace & cpout = dynamic_cast<segytrace &>(x);
      segytrace const & cpin = dynamic_cast<segytrace const &>(y);
            
      segy & trout = cpout.getMetadata();
      segy const & trin = cpin.getMetadata();
            
      Value val;
      float gx;
      float sx;
      float x;
      float dt;
      float t0;
      float wm;
      float wt;
      int nt;
            
      string name="ns";
      gethdval(&trin,(char *)(name.c_str()),&val);
      nt=vtoi(hdtype((char *)(name.c_str())),val);
            
      name="dt";
      gethdval(&trin,(char *)(name.c_str()),&val);
      dt=0.001*vtof(hdtype((char *)(name.c_str())),val);
            
      name="delrt";
      gethdval(&trin,(char *)(name.c_str()),&val);
      t0=vtof(hdtype((char *)(name.c_str())),val);
            
      name="gx";
      gethdval(&trin,(char *)(name.c_str()),&val);
      gx = vtof(hdtype((char *)(name.c_str())),val);
            
      name="sx";
      gethdval(&trin,(char *)(name.c_str()),&val);
      sx = vtof(hdtype((char *)(name.c_str())),val);
      x = gx - sx; 
            
      int itw = (int) (tw/dt + 0.1);
      itw = MAX(1,itw);

      float sxt_wt = 1.0f;
      wt = MAX(1.0,sx_width);
      if (sx < sx_min) sxt_wt = cosfun((sx_min-sx)/wt);
      if (sx > sx_max) sxt_wt = cosfun((sx-sx_max)/wt);
      //cerr << " sxt_wt = " << sxt_wt << endl;
      // adjust mute so that w=0 works
      if (mute_type) {
	wm = MAX(dt/s,w);
	float wtmp = mutefun((gx - s)/wm) * mutefun((tm - gx)/wm);
	wtmp = wtmp * sxt_wt;
	wt=MAX(1.0,width);
	float wttmp=1.0;
	if (taper_type){ // offset
	  if (x < taper_min) wttmp = cosfun((taper_min-x)/wt);
	  else if (x > taper_max) wttmp = cosfun((x - taper_max)/wt);
#pragma ivdep                    
	  for (int i=0;i<nt-itw;i++)
	    trout.data[i] = trin.data[i]*wtmp*wttmp;
#pragma ivdep
	  for (int i=nt-itw;i<nt;i++)
	    trout.data[i] = trin.data[i]*wtmp*wttmp*cosfun(float(nt-i)/(itw+0.0f));

	}
	else { // taper geophone position
	  if (gx < taper_min) wttmp = cosfun((taper_min-gx)/wt);
	  else if (gx > taper_max) wttmp = cosfun((gx - taper_max)/wt);
#pragma ivdep                    
	  for (int i=0;i<nt-itw;i++)
	    trout.data[i] = trin.data[i]*wtmp*wttmp;
#pragma ivdep
	  for (int i=nt-itw;i<nt;i++)
	    trout.data[i] = trin.data[i]*wtmp*wttmp*cosfun(float(nt-i)/(itw+0.0f));
	}
      }
      else {
	wm=MAX(dt,w);
	float wttmp=1.0;
	wt=MAX(1.0,width);
	if (taper_type){ // offset
	  if (x < taper_min) wttmp = cosfun((taper_min-x)/wt);
	  else if (x > taper_max) wttmp = cosfun((x - taper_max)/wt);
	  wttmp = wttmp * sxt_wt;
#pragma ivdep                    
	  for (int i=0;i<nt-itw;i++)
	    trout.data[i] = trin.data[i]*wttmp*mutefun((i*dt+t0-s*fabs(x)-tm)/wm);
#pragma ivdep
	  for (int i=nt-itw;i<nt;i++)
	    trout.data[i] = trin.data[i]*wttmp*mutefun((i*dt+t0-s*fabs(x)-tm)/wm)*cosfun(float(nt-i)/(itw+0.0f));
	}
	else { // taper geophone position
	  if (gx < taper_min) wttmp = cosfun((taper_min-gx)/wt);
	  else if (gx > taper_max) wttmp = cosfun((gx - taper_max)/wt);
	  wttmp = wttmp * sxt_wt;
#pragma ivdep
	  for (int i=0;i<nt-itw;i++)
	    trout.data[i] = trin.data[i]*wttmp*mutefun((i*dt+t0-s*fabs(x)-tm)/wm);
#pragma ivdep
	  for (int i=nt-itw;i<nt;i++)
	    trout.data[i] = trin.data[i]*wttmp*mutefun((i*dt+t0-s*fabs(x)-tm)/wm)*cosfun(float(nt-i)/(itw+0.0f));
	}
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: SEGYTaperMute::operator()\n";
      e<<"input LDC not segytrace\n";
      throw e;
    }
  }

  void SEGYFwdDiffIP::operator()(LocalDataContainer<float> & x) {
	 
    try {
            
      segytrace & cpout = dynamic_cast<segytrace &>(x);
            
      segy & trout = cpout.getMetadata();
            
      Value val;
      float dt;
      int nt;
            
      string name="ns";
      gethdval(&trout,(char *)(name.c_str()),&val);
      nt=vtoi(hdtype((char *)(name.c_str())),val);
      name="dt";
      gethdval(&trout,(char *)(name.c_str()),&val);
      dt=0.001*vtof(hdtype((char *)(name.c_str())),val);
            
      for (int j=0;j<nint;j++) {
#pragma ivdep
	for (int i=nt;i>0;i--) trout.data[i]
				 = (trout.data[i] - trout.data[i-1])/dt;
	trout.data[0] /= dt;
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: SEGYFwdDiffIP::operator()\n";
      e<<"input LDC not segytrace\n";
      throw e;
    }
  }  
    
  void SEGYFwdIntIP::operator()(LocalDataContainer<float> & x) {
	 
    try {
            
      segytrace & cpout = dynamic_cast<segytrace &>(x);
            
      segy & trout = cpout.getMetadata();
            
      Value val;
      float dt;
      int nt;
            
      string name="ns";
      gethdval(&trout,(char *)(name.c_str()),&val);
      nt=vtoi(hdtype((char *)(name.c_str())),val);
      name="dt";
      gethdval(&trout,(char *)(name.c_str()),&val);
      dt=0.001*vtof(hdtype((char *)(name.c_str())),val);
            
      for (int j=0;j<nint;j++) {
	float prev = trout.data[0];
	float curr = 0.0f;
	trout.data[0] = 0.0f;
#pragma ivdep
	for (int i=1;i<nt;i++) {
	  curr = trout.data[i];
	  trout.data[i] = 0.5*dt*(trout.data[i]+prev) + trout.data[i-1];
	  prev = curr;
	}
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: SEGYFwdIntIP::operator()\n";
      e<<"input LDC not segytrace\n";
      throw e;
    }
  }
    
  void SEGYAdjIntIP::operator()(LocalDataContainer<float> & x) {        
    try {
            
      segytrace & cpout = dynamic_cast<segytrace &>(x);
      
      segy & trout = cpout.getMetadata();
            
      Value val;
      float dt;
      int nt;
            
      string name="ns";
      gethdval(&trout,(char *)(name.c_str()),&val);
      nt=vtoi(hdtype((char *)(name.c_str())),val);
      name="dt";
      gethdval(&trout,(char *)(name.c_str()),&val);
      dt=0.001*vtof(hdtype((char *)(name.c_str())),val);
            
      for (int j=0;j<nint;j++) {
	float prev = trout.data[nt-1];
	float curr = 0.0f;
	trout.data[nt-1] = 0.5*dt*prev;
#pragma ivdep
	// prev = i+1
	for (int i=nt-2;i>0;i--) {
	  curr = trout.data[i];
	  trout.data[i] = 0.5*dt*(trout.data[i]+prev) + trout.data[i+1];
	  prev = curr;
	}
	trout.data[0]=0.5*(0.5*dt*prev + trout.data[1]);
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: SEGYAdjIntIP::operator()\n";
      e<<"input LDC not segytrace\n";
      throw e;
    }
  }

  void readTraceSampleFO::operator()(LocalDataContainer<float> const & x) {
    try {
      segytrace const & cpout = dynamic_cast<segytrace const &>(x);
      segy const & trout = cpout.getMetadata();
      if (!init) {
	
	Value val;
	
	std::string str="dt";
	gethdval(&trout,(char*)(str.c_str()), &val);
	if (RVL::ProtectedDivision<float>(1000.0f,vtof(hdtype(str.c_str()), val),dtr)) {
	  RVL::RVLException e;
	  e<<"Error: readTraceSample::operator()\n";
	  e<<"  zerodivide by dt\n";
	  throw e;
	}
	str="ns";
	gethdval(&trout,  (char*)(str.c_str()), &val);
	n=vtoi(hdtype(str.c_str()), val);
	str="delrt";
	gethdval(&trout,  (char*)(str.c_str()), &val);
	
	beg=int(0.01 + vtof(hdtype(str.c_str()), val)*dtr);
	
	if (beg>0 || beg<-(n-1)) {
	  RVLException e;
	  e<<"Error: readTraceSampleFO::operator()\n";
	  e<<"  trace header with dt="<<1.0f/dtr<<" ns="<<n<<" beg="<<beg<<"\n";
	  throw e;
	}
	init=true;
      }

      // fwd map scales by 1/dt re discrete delta - adjoint is just
      // copy because of dt factor in norm on trace side
      if (idx<v.getSize()) {
	v.getData()[idx]=trout.data[-beg];
	idx++;
      }
      else {
	RVLException e;
	e<<"Error: readTraceSample::operator()\n";
	e<<"  too many input traces - more than "<<v.getSize()<<"\n";
	throw e;
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"ERROR: readTraceSample::operator()\n";
      e<<"  at least one input not segytrace type\n";
      throw e;
    }
  }

  void writeTraceSampleFO::operator()(LocalDataContainer<float> & x) {
    try {
      segytrace & cpout = dynamic_cast<segytrace &>(x);
      segy & trout = cpout.getMetadata();
      if (!init) {
	
	Value val;
	
	std::string str="dt";
	gethdval(&trout,(char*)(str.c_str()), &val);
	if (RVL::ProtectedDivision<float>(1000.0f,vtof(hdtype(str.c_str()), val),dtr)) {
	  RVL::RVLException e;
	  e<<"Error: writeTraceSample::operator()\n";
	  e<<"  zerodivide by dt\n";
	  throw e;
	}
	str="ns";
	gethdval(&trout,  (char*)(str.c_str()), &val);
	n=vtoi(hdtype(str.c_str()), val);
	str="delrt";
	gethdval(&trout,  (char*)(str.c_str()), &val);
	
	beg=int(0.01 + vtof(hdtype(str.c_str()), val)*dtr);
	
	if (beg>0 || beg<-(n-1)) {
	  RVLException e;
	  e<<"Error: writeTraceSampleFO::operator()\n";
	  e<<"  trace header with dt="<<1.0f/dtr<<" ns="<<n<<" beg="<<beg<<"\n";
	  throw e;
	}
	init=true;
      }

      // fwd map scales by 1/dt re discrete delta - adjoint is just
      // copy because of dt factor in norm on trace side
      if (idx<v.getSize()) {
	trout.data[-beg]=dtr*v.getData()[idx];
	idx++;
      }
      else {
	RVLException e;
	e<<"Error: writeTraceSample::operator()\n";
	e<<"  too many input traces - more than "<<v.getSize()<<"\n";
	throw e;
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"ERROR: writeTraceSample::operator()\n";
      e<<"  at least one input not segytrace type\n";
      throw e;
    }
  }

  void OpZFwdFO::operator()(LocalDataContainer<float> & x,
			    LocalDataContainer<float> const & y) {
    try {
      segytrace & cpout = dynamic_cast<segytrace &>(x);
      segy & trout = cpout.getMetadata();
      segytrace const & cpin = dynamic_cast<segytrace const &>(y);
      segy const & trin = cpin.getMetadata();
	
      Value val;
      int nin;
      int nout;
      string nsstr="ns";      
      float fscalel;
      float fselev;
      float fgelev;
      float fac=1.0f;
      char * key = (char *)usermalloc_(128*sizeof(char));

      // verify that output trace has zero src/rec depth
      strcpy(key,"scalel");
      gethdval(&trout,key,&val);
      fscalel = vtof(hdtype(key),val);
      strcpy(key,"selev");
      gethdval(&trout,key,&val);
      fselev = vtof(hdtype(key),val);
      strcpy(key,"gelev");
      gethdval(&trout,key,&val);
      fgelev = vtof(hdtype(key),val);	
      if (fscalel > 0) { fselev *=  fscalel; fgelev *=  fscalel; }
      if (fscalel < 0) { fselev /= -fscalel; fgelev /= -fscalel;}
      if (((fabs(fselev)+1.0f) > 1.0f) || ((fabs(fgelev)+1.0f) > 1.0f)) {
	RVL::RVLException e;
	e<<"ERROR: OpZFwdFO\n";
	e<<"  output trace has nonzero src/rec depths\n";
	throw e;
      }

      // extract selev, gelev from input trace
      strcpy(key,"scalel");
      gethdval(&trin,key,&val);
      fscalel = vtof(hdtype(key),val);
      strcpy(key,"selev");
      gethdval(&trin,key,&val);
      fselev = vtof(hdtype(key),val);
      strcpy(key,"gelev");
      gethdval(&trin,key,&val);
      fgelev = vtof(hdtype(key),val);	
      if (fscalel > 0) { fselev *=  fscalel; fgelev *=  fscalel; }
      if (fscalel < 0) { fselev /= -fscalel; fgelev /= -fscalel;}

      if (RVL::ProtectedDivision<float>(1.0f,4.0f*fselev*fgelev,fac)) {
	RVLException e;
	e<<"Error: OpZFwdFO\n";
	e<<"  one or both free surface data src/rec depths vanishes\n";
	throw e;
      }
	
      gethdval(&trin,  (char*)(nsstr.c_str()), &val);
      nin=vtoi(hdtype(nsstr.c_str()), val);
	
      gethdval(&trout,  (char*)(nsstr.c_str()), &val);
      nout=vtoi(hdtype(nsstr.c_str()), val);
	
      if (nin != nout) {
	RVL::RVLException e;
	e<<"Error: OpZFwdFO\n";
	e<<"  traces of unequal length\n";
	throw e;
      }

#pragma ivdep
      for (int i=0; i<nin; i++) { (trout.data)[i] = fac*(trin.data)[i]; }

      userfree_(key);
    }
    catch (bad_cast) {
      RVL::RVLException e;
      e<<"ERROR: OpZFwdFO::operator()\n";
      e<<"  at least one input not segytrace type\n";
      throw e;
    }
  }

  void OpZAdjFO::operator()(LocalDataContainer<float> & x,
			    LocalDataContainer<float> const & y) {
    try {
      segytrace & cpout = dynamic_cast<segytrace &>(x);
      segy & trout = cpout.getMetadata();
      segytrace const & cpin = dynamic_cast<segytrace const &>(y);
      segy const & trin = cpin.getMetadata();
	
      Value val;
      int nin;
      int nout;
      string nsstr="ns";      
      float fscalel;
      float fselev;
      float fgelev;
      float fac=1.0f;
      char * key = (char *)usermalloc_(128*sizeof(char));
	
      // verify that input trace has zero src/rec depths
      strcpy(key,"scalel");
      gethdval(&trin,key,&val);
      fscalel = vtof(hdtype(key),val);
      strcpy(key,"selev");
      gethdval(&trin,key,&val);
      fselev = vtof(hdtype(key),val);
      strcpy(key,"gelev");
      gethdval(&trin,key,&val);
      fgelev = vtof(hdtype(key),val);	
      if (fscalel > 0) { fselev *=  fscalel; fgelev *=  fscalel; }
      if (fscalel < 0) { fselev /= -fscalel; fgelev /= -fscalel;}
      if (((fabs(fselev)+1.0f) > 1.0f) || ((fabs(fgelev)+1.0f) > 1.0f)) {
	RVL::RVLException e;
	e<<"ERROR: OpZAdjFO\n";
	e<<"  input trace has nonzero src/rec depths\n";
	throw e;
      }

      // extract src/rec depths from output trace
      strcpy(key,"scalel");
      gethdval(&trout,key,&val);
      fscalel = vtof(hdtype(key),val);
      strcpy(key,"selev");
      gethdval(&trout,key,&val);
      fselev = vtof(hdtype(key),val);
      strcpy(key,"gelev");
      gethdval(&trout,key,&val);
      fgelev = vtof(hdtype(key),val);	
      if (fscalel > 0) { fselev *=  fscalel; fgelev *=  fscalel; }
      if (fscalel < 0) { fselev /= -fscalel; fgelev /= -fscalel;}

      if (RVL::ProtectedDivision<float>(1.0f,4.0f*fselev*fgelev,fac)) {
	RVLException e;
	e<<"Error: OpZAdjFO\n";
	e<<"  one or both free surface data src/rec depths vanishes\n";
	throw e;
      }
	
      gethdval(&trin,  (char*)(nsstr.c_str()), &val);
      nin=vtoi(hdtype(nsstr.c_str()), val);
	
      gethdval(&trout,  (char*)(nsstr.c_str()), &val);
      nout=vtoi(hdtype(nsstr.c_str()), val);
	
      if (nin != nout) {
	RVL::RVLException e;
	e<<"Error: OpZAdjFO\n";
	e<<"  traces of unequal length\n";
	throw e;
      }

#pragma ivdep	
      for (int i=0; i<nin; i++) { (trout.data)[i] = fac*(trin.data)[i]; }
    }
    catch (bad_cast) {
      RVL::RVLException e;
      e<<"ERROR: OpZAdjFO::operator()\n";
      e<<"  at least one input not segytrace type\n";
      throw e;
    }
  }

  void SRCellVolFO::operator()(LocalDataContainer<float> const & y) {
    try {
      segytrace const & cpin = dynamic_cast<segytrace const &>(y);
      segy const & trin = cpin.getMetadata();
	
      Value val;
      float fscalco;
      float fsx;
      float fsy;      
      float fgx;
      float fgy;
      
      char * key = (char *)usermalloc_(128*sizeof(char));
	
      strcpy(key,"scalco");
      gethdval(&trin,key,&val);
      fscalco = vtof(hdtype(key),val);
      strcpy(key,"sx");
      gethdval(&trin,key,&val);
      fsx = vtof(hdtype(key),val);
      strcpy(key,"sy");
      gethdval(&trin,key,&val);
      fsy = vtof(hdtype(key),val);	
      strcpy(key,"gx");
      gethdval(&trin,key,&val);
      fgx = vtof(hdtype(key),val);
      strcpy(key,"gy");
      gethdval(&trin,key,&val);
      fgy = vtof(hdtype(key),val);	
      if (fscalco > 0) {
	fsx *=  fscalco; fsy *=  fscalco; 
	fgx *=  fscalco; fgy *=  fscalco;
      }
      if (fscalco < 0) {
	fsx /=  -fscalco; fsy /=  -fscalco; 
	fgx /=  -fscalco; fgy /=  -fscalco;
      }

      // subsequent traces
      // if any of the coordinates does not increase
      // at all, then its cell factor remains 1.
      if (init) {
	if (fgy > gylast) 
	  dgy = iwave_max(dgy,fgy-gylast);
	if (fgx > gxlast) 
	  dgx = iwave_max(dgx,fgx-gxlast);
	if (fsy > sylast) 
	  dsy = iwave_max(dsy,fsy-sylast);
	if (fsx > sxlast) 
	  dsx = iwave_max(dsx,fsx-sxlast);
      }
      // first trace
      else {
	init=true;
      }

      sxlast=fsx;
      sylast=fsy;
      gxlast=fgx;
      gylast=fgy;

      //      if (retrieveGlobalRank()==0) 
      //	cerr<<"dgx="<<dgx<<" dgy="<<dgy<<" dsx="<<dsx<<" dsy="<<dsy<<"\n";
      // now set cell volume estimate from this trace
      RVL::ScalarRedn<float>::setValue(dgx*dgy*dsx*dsy);
    }
    catch (bad_cast) {
      RVL::RVLException e;
      e<<"ERROR: SRCellVolFO::operator()\n";
      e<<"  at least one input not segytrace type\n";
      throw e;
    }
  }
  
  void TraceScaleFO::operator()(LocalDataContainer<float> & xout,
				LocalDataContainer<float> const & xin) {
    try {
      segytrace const & cpin = dynamic_cast<segytrace const &>(xin);
      segy const & trin = cpin.getMetadata();
      segytrace & cpout = dynamic_cast<segytrace &>(xout);
      segy & trout = cpout.getMetadata();
	
      Value val;
      int nin;
      int nout;
      float hfac;
      string nsstr="ns";      
	
      gethdval(&trin,  (char*)(nsstr.c_str()), &val);
      nin=vtoi(hdtype(nsstr.c_str()), val);
      gethdval(&trout,  (char*)(nsstr.c_str()), &val);
      nout=vtoi(hdtype(nsstr.c_str()), val);
      
      if (nin != nout) {
	RVL::RVLException e;
	e<<"ERROR: TraceScaleFO::operator()\n";
	e<<"  input length = "<<nin<<" != output length = "<<nout<<"\n";
	throw e;
      }
	
      if (key.size() > 0) {
	gethdval(&trout, (char*)(key.c_str()), &val);
	hfac=vtof(hdtype(key.c_str()),val);
	if (f) {
	  hfac = (*f)(hfac);
	}
	else {
	  hfac *= fac;
	}
      }
      else {
	hfac=fac;
      }
      
#pragma ivdep
      for (int i=0; i<nout; i++) { (trout.data)[i] = hfac*(trin.data)[i]; }
    }
    catch (bad_cast) {
      RVL::RVLException e;
      e<<"ERROR: TraceScaleFO::operator()\n";
      e<<"  at least one input not segytrace type\n";
      throw e;
    }
    catch (RVL::RVLException & e) {
      e<<"\ncalled from TraceScaleFO::operator()\n";
      throw e;
    }
  }

  void SEGYTFTScaleFO::operator()(LocalDataContainer<float> & x) {
    try {
      segytrace & cpout = dynamic_cast<segytrace &>(x);
      segy & trout = cpout.getMetadata();
	
      Value val;
      string dtstr="dt";      
	
      gethdval(&trout,  (char*)(dtstr.c_str()), &val);
      // convert to ms
      float dt = 0.001*vtof(hdtype(dtstr.c_str()), val);

      FTScaleFO ft(dt,power,ab,locut,lopas,hipas,hicut);

      ft(x);
      
    }
    catch (bad_cast) {
      RVL::RVLException e;
      e<<"ERROR: SEGYTFTScaleFO::operator()\n";
      e<<"  at least one input not segytrace type\n";
      throw e;
    }
    catch (RVL::RVLException & e) {
      e<<"\ncalled from SEGYTFTScaleFO::operator()\n";
      throw e;
    }    
  }
}


