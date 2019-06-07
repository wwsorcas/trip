#include "gridtrace.hh"

void createExplReflTraces(std::string gridname,
			  std::string tracename,
			  std::string dataname,
			  std::string cwproot,
			  FILE * stream) {
  try {
    grid g;
    init_default_grid(&g);
    if (read_grid(&g, gridname.c_str(),stream)) {
      RVL::RVLException e;
      e<<"Error: createExplReflTraces\n";
      e<<"  failed to open grid file "<<gridname<<"\n";
      throw e;
    }

    // for the moment limited to 2D
    if (get_dimension_grid(g) !=2) {
      RVL::RVLException e;
      e<<"Error: createExplReflTraces\n";
      e<<"  only 2D case currently supported\n";
      throw e;
    }
    IPNT gs0;
    IPNT ge0;
    IPNT n;
    get_gs(gs0,g);
    get_ge(ge0,g);
    get_n(n,g);
    RPNT d;
    RPNT o;
    get_d(d,g);
    get_o(o,g);

    FILE * fpin = NULL;
    Value val;
    float dt;
    segy tr;
    
#ifdef IWAVE_USE_FMGR
    /* iwave_fopen version */
    if (!(fpin=iwave_const_fopen(dataname.c_str(),"r+",NULL,stream))) {
      RVL::RVLException e;
      e<<"Error: createExplReflTraces\n";
      e<<"  failed to open data file\n";
      throw e;
    }
#else
    /* stdio version */
    if (!(fpin=fopen(dataname.c_str(),"r+"))) {
      RVL::RVLException e;
      e<<"Error: createExplReflTraces\n";
      e<<"  failed to open data file\n";
      throw e;
    }
#endif
    fgettr(fpin,&tr);
    std::string sdt = "dt";
    gethdval(&tr,(char*)(sdt.c_str()),&val);
    dt = vtof(hdtype(sdt.c_str()),val);
    cerr<<"dt="<<dt<<endl;
    
#ifdef IWAVE_USE_FMGR
    iwave_fclose(fpin);
#else
    fclose(fpin);
#endif

    std::stringstream cmd;
    cmd << cwproot<<"/bin/suplane npl=1 nt=11 ntr=" << n[0]*n[1] << " ix1=0 nspk=1 ct1=6 dt="
	<< 1.e-6*dt << " offset=0 len1=" << n[0]*n[1]
	<< " | "<<cwproot<<"/bin/sushw key=delrt,gelev,gx a=" << -0.005*dt<<","
	<<-gs0[0]*d[0]<<","<<gs0[1]*d[1]
	<< " b=0,"<<-d[0]<<",0 c=0,0,"<<d[1]<<" j="<<n[0]<<","<<n[0]<<","<<n[0]
	<<"  > "<< tracename <<" \n";

    cerr<<"cmd = "<<cmd.str().c_str();
    
    if (system(cmd.str().c_str())) {
      RVL::RVLException e;
      e<<"Error: createExplReflTraces\n";
      e<<"  system call failed\n";
      throw e;
    }
  }
  catch (RVL::RVLException & e) {
    e<<"\ncalled from createExplReflTraces\n";
    throw e;
  }
}

