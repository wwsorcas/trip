#include <su.h>
#include <segy.h>
#include <header.h>
#include <parser.h>

using RVL::valparse;
using RVL::RVLException;

#define TOL 0.001

const char * sdoc[] = { 
  "Usage: wiw.x data= jvpm= alpha= mmin= mmax= nm=",
  "",
  "Synopsis: Implements formulas (38)-(41) in \"Solution of an Acoustic",
  "Transmission Inverse Problem by Source Extension\". ",
  "The source-receiver distance r is fixed at 1 km.",
  "Output is the reduced DSO objective J_VPM, the mean square error e,",
  "and the penalty term p, sampled over the",
  "interval [mmin, mmax] at sample interval dm = (mmax-mmin)/(nm-1). The",
  "data is read from an SU file, the output J_VPM, e, and p samples are ",
  "organized in an RSF file pair.",
  "",
  "Required inputs:",
  "  data [string]     = input data trace",
  "  jvpm [string]     = output jvpm samples",
  "  alpha [float]     = penalty parameter",
  "  mmin [float]      = min slowness",
  "  mmax [float]      = max slowness",
  "  nm [int]          = number of slowness samples",
  "",
  "NB: working units are s, km. File units must be scaled: SEGY stores dt in",
  "microsec, delrt (= t0) in ms.",
  NULL};

int main(int argc, char ** argv) {

  try {
  
    PARARRAY * par;    /* param array */
 
    std::string data;     /* data file name */
    std::string jvpm;     /* output file name */
    float alpha;     /* */
    FILE * fpw;      /* data file pointer */
    FILE * fph;      /* jvpm header file pointer */
    FILE * fpd;      /* jvpm data file pointer */
    segy tr;         /* data trace workspace */
    Value val;       /* header word workspace */
    float dt;        /* time step */
    float t0;        /* start time of wavelet series */
    int nt;          /* number of samples per trace */
    float dm;        /* slowness interval */
    float mmin;      /* min slowness */
    float mmax;      /* max slowness */
    int nm;          /* number of slowness samples */
    float * e;       /* array of scaled data error values */
    float * p;       /* array of scaled annihilator values */
    float * ap;      /* array of scaled annihilator values */
    float * jval;    /* array of jvpm values */
    float eplus;     /* upper residual bound for alpha update */
    float eminus;    /* lower residual bound for alpha update */
    float efactor;   /* factor by which lower residual bound is less than upper */

    xargc=argc; xargv=argv;
    requestdoc(1);

    /* extract input parameters */
    par=ps_new();
    if ( ps_createargs(par, argc - 1, argv + 1) ) {
      RVLException e;
      e<<"ERROR: wiw from ps_creatargs \n";
      e<<"  called with args:\n";
      e<<"  argc = "<<argc-1<<"\n";
      for (int i=0;i<argc-1;i++) 
	e<<"  argv["<<i<<"] = "<<argv[i+1]<<"\n";
      throw e;
    }
  
    data   = RVL::valparse<std::string>(*par,"data");
    jvpm   = RVL::valparse<std::string>(*par,"jvpm");
    alpha  = RVL::valparse<float>(*par,"alpha");
    mmin   = RVL::valparse<float>(*par,"mmin");
    mmax   = RVL::valparse<float>(*par,"mmax");
    nm     = RVL::valparse<int>(*par,"nm");
    eplus  = RVL::valparse<float>(*par,"eplus",0.0f);
    efactor= RVL::valparse<float>(*par,"efactor",0.9f);
    if ((efactor<0.01) || (efactor>0.99f)) {
      RVLException e;
      e<<"Error: WIW\n";
      e<<"  efactor="<<efactor<<" should be in [0.01, 0.99]\n";
      throw e;
    }
    eminus = efactor*eplus;
    if (nm<2) {
      RVL::RVLException e;
      e<<"ERROR: wiw\n";
      e<<"  nm =  "<<nm<<" must be at least 2\n";
      throw e;
    }
    dm = (mmax-mmin)/(nm-1);

    /* open data files */
    if (!(fpw=fopen(data.c_str(),"r"))) {
      RVL::RVLException e;
      e<<"ERROR: wiw\n";
      e<<"  failed to open input file "<<data<<"\n";
      throw e;
    }

    if (!(fph=fopen(jvpm.c_str(),"w"))) {
      RVL::RVLException e;
      e<<"ERROR: wiw\n";
      e<<"  failed to open input file "<<jvpm<<"\n";
      throw e;
    }

    std::string jvpd = "/var/tmp/"+jvpm + "@";
    if (!(fpd=fopen(jvpd.c_str(),"w"))) {
      RVL::RVLException e;
      e<<"ERROR: wiw\n";
      e<<"  failed to open input file "<<jvpm<<"\n";
      throw e;
    }

    /* read data */
    if (!(fgettr(fpw,&tr))) {
      RVL::RVLException e;
      e<<"ERROR: wiw\n";
      e<<"  failed to read trace from "<<data<<"\n";
      throw e;
    }

    /* workspace to make up for gethdval not taking const args */
    char * str = new char[128];
    
    /* read time sample info */
    strcpy(str,"ns");
    gethdval(&tr,str,&val);
    nt = vtoi(hdtype(str),val);

    /* file dt in \mu s; scale to s */
    strcpy(str,"dt");
    gethdval(&tr,str,&val);
    dt = 0.000001*vtof(hdtype(str),val);

    /* file t0 in ms; scale to s */
    strcpy(str,"delrt");
    gethdval(&tr,str,&val);
    t0 = 0.001*vtof(hdtype(str),val);
    delete [] str;
    
    /* workspace for jvpm values */
    e = new float[nm];
    p = new float[nm];
    ap = new float[nm];
    jval = new float[nm];

    /* equations 38-40 in WIW paper */
    /* loop over slowness */

    float jopt=std::numeric_limits<float>::max();
    float eopt=std::numeric_limits<float>::max();
    float popt=std::numeric_limits<float>::max();
    float mopt=mmin;
    
    float pi = 3.1415927;
    float fac = (4.0*pi)*(4.0*pi);
    float alphasq = alpha*alpha;
    float dns = 0.0f;
    for (int i=0;i<nt;i++) {
      dns += 0.5*dt*(tr.data[i])*(tr.data[i]);
    }
    // cout<<"0.5*data norm squared = "<<dns<<endl;
    for (int k=0; k<nm; k++) {
      float m = mmin + k*dm;
      e[k]=0.0f;
      p[k]=0.0f;
      for (int i=0;i<nt;i++) {
	/* scale times to s */
	float t = t0 + i*dt;
	// (4pi(t-m))^2
	float amp0 = alphasq*fac*(t-m)*(t-m);
	float amp1 = amp0/(1.0f + amp0);
	float amp2 = fac*(t-m)*(t-m)/((1.0f + amp0)*(1.0f + amp0));
	float wsq  = (tr.data[i])*(tr.data[i]);
	// accumulate (4pialpha(t-m))^4/(1+(4pialpha(t-m))^2)^2
	// * d(t)^2
	e[k]+=amp1*amp1*wsq;
	// accumulate (t-m)^2/(1+(4pialpha(t-m))^2)^2 * d(t)^2
	p[k]+=amp2*wsq;
      }
      e[k]*=0.5*dt;
      p[k]*=0.5*dt;
      ap[k]=alphasq*p[k];
      jval[k]=e[k]+ap[k];
      // normalize by 1/2|d|^2 = lim alpha -> \infty
      e[k]/=dns;
      p[k]/=dns;
      ap[k]/=dns;
      jval[k]/=dns;
      // update minimizer
      if (jval[k] <= jopt) {
	jopt = jval[k];
	eopt = e[k];
	popt = p[k];
	mopt = m;
      }
      
    }
    if (eplus > 0) {
      cout<<"jopt="<<jopt<<" eopt="<<eopt<<" popt="<<popt<<" mopt="<<mopt<<endl;
      if (eopt > eplus) {
	cout<<"eopt="<<eopt<<" > eplus="<<eplus<<": reduce alpha\n";
      }
      else if (eopt < eminus) {
	float alphaplus = sqrt(alphasq + (eplus-eopt)/(2.0f*popt));
	cout<<"recommended next alpha = "<<alphaplus<<"\n";
      }
      else {
	cout<<"alpha looks good, carry on\n";
      }
    }
      
    /* write jvpm data file */
    fwrite(jval,sizeof(float),nm,fpd);
    fwrite(e,sizeof(float),nm,fpd);
    fwrite(ap,sizeof(float),nm,fpd);

    delete [] jval;
    delete [] e;
    delete [] p;
    delete [] ap;
    
    /* write jvpm header file */
    fprintf(fph,"n1=%d ",nm);
    fprintf(fph,"d1=%g ",dm);
    fprintf(fph,"o1=%g\n",mmin);
    fprintf(fph,"n2=%d ",3);
    fprintf(fph,"d2=%g ",1.0f);
    fprintf(fph,"o2=%g\n",0.0f);    
    fprintf(fph,"esize=4 data_format=\"native_float\"\n");
    fprintf(fph,"in=\"%s\"\n",jvpd.c_str());

    /* close files */
    fclose(fpw);
    fclose(fph);
    fclose(fpd);
    ps_delete(&par);

  }
  catch (RVL::RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
  
