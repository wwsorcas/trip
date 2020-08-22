#include <su.h>
#include <segy.h>
#include <header.h>
#include <parser.h>

using RVL::valparse;
using RVL::RVLException;

#define TOL 0.001

const char * sdoc[] = { 
  "Usage: wiw.x wavelet= jvpm= lambda= alpha= mex= mmin= mmax= nm=",
  "",
  "Synopsis: Implements formula (18) in \"Full waveform inversion by source",
  "extension: why it works\". The source-receiver distance r is fixed at 1 km.",
  "The interval on which the wavelet is defined should include the subinterval",
  "[-lambda,lambda] - this is checked. Output is J_VPM sampled over the",
  "interval [mmin, mmax] at sample interval dm = (mmax-mmin)/(nm-1). The",
  "wavelet is read from an SU file, the output J_VPM samples are organized",
  "in an RSF file pair.",
  "",
  "Required inputs:",
  "  wavelet [string]  = input source wavelet, root name (missing .su)",
  "  jvpm [string]     = output jvpm samples, root name (missing .rsf)",
  "  lambda [float]    = radius of assumed wavelet support",
  "  alpha [float]     = penalty parameter",
  "  mex [float]       = target (\"exact\") slowness",
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
 
    std::string wav;      /* wavelet file name */
    std::string jvpm;     /* output file name */
    float lambda;    /* */
    float alpha;     /* */
    FILE * fpw;      /* wavelet file pointer */
    FILE * fph;      /* jvpm header file pointer */
    FILE * fpd;      /* jvpm data file pointer */
    segy tr;         /* wavelet trace workspace */
    Value val;       /* header word workspace */
    float dt;        /* time step */
    float t0;        /* start time of wavelet series */
    int nt;          /* number of samples per trace */
    float dm;        /* slowness interval */
    float mmin;      /* min slowness */
    float mmax;      /* max slowness */
    float mex;       /* target ("exact") slowness */
    int nm;          /* number of slowness samples */
    float * e;       /* array of scaled data error values */
    float * p;       /* array of scaled annihilator values */
    float * jval;    /* array of jvpm values */

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
  
    wav    = RVL::valparse<std::string>(*par,"wavelet");
    jvpm   = RVL::valparse<std::string>(*par,"jvpm");
    lambda = RVL::valparse<float>(*par,"lambda");
    alpha  = RVL::valparse<float>(*par,"alpha");
    mex    = RVL::valparse<float>(*par,"mex");
    mmin   = RVL::valparse<float>(*par,"mmin");
    mmax   = RVL::valparse<float>(*par,"mmax");
    nm     = RVL::valparse<int>(*par,"nm");
    if (nm<2) {
      RVL::RVLException e;
      e<<"ERROR: wiw\n";
      e<<"  nm =  "<<nm<<" must be at least 2\n";
      throw e;
    }
    dm = (mmax-mmin)/(nm-1);

    /* open data files */
    if (!(fpw=fopen(wav.c_str(),"r"))) {
      RVL::RVLException e;
      e<<"ERROR: wiw\n";
      e<<"  failed to open input file "<<wav<<"\n";
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

    /* read wavelet */
    if (!(fgettr(fpw,&tr))) {
      RVL::RVLException e;
      e<<"ERROR: wiw\n";
      e<<"  failed to read trace from "<<wav<<"\n";
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
    
    /* check time interval */
    if ((-lambda < t0) || (lambda > t0+(nt-1)*dt)) {
      RVL::RVLException e;
      e<<"ERROR: wiw\n";
      e<<"  target support interval [-"<<lambda<<", "<<lambda<<"] not \n";
      e<<"  subset of wavelet definition interval ";
      e<<"["<<t0<<", "<<t0+(nt-1)*dt<<"\n";
      throw e;
    }      

    /* workspace for jvpm values */
    e = new float[nm];
    p = new float[nm];
    jval = new float[nm];

    /* equation 18 in WIW paper */
    /* loop over slowness */
    
    float pi = 3.1415927;
    float fac = (4.0*pi)*(4.0*pi)*alpha*alpha;
    float scl = dt/(32*pi*pi);
    float dns = 0.0f;
    for (int i=0;i<nt;i++) {
      dns += (tr.data[i])*(tr.data[i]);
    }
    dns *= scl;
    cout<<"0.5*data norm squared = "<<dns<<endl;
    for (int k=0; k<nm; k++) {
      float m = mmin + k*dm;
      e[k]=0.0f;
      p[k]=0.0f;
      for (int i=0;i<nt;i++) {
	/* scale times to s */
	float t = t0 + i*dt;
	// (4pialpha(t+(m*-m)))^2
	float amp0 = (fac*(t+(mex-m))*(t+(mex-m)));
	float amp1 = amp0/(1.0f + amp0);
	float amp2 = (t+(mex-m))/(1.0f + amp0);
	float wsq  = (tr.data[i])*(tr.data[i]);
	// accumulate (4pialpha(t+(m*-m)))^4/(1+(4pialpha(t+(m*-m)))^2)^2
	// * w(t)^2
	e[k]+=amp1*amp1*wsq;
	// accumulate (t+(m*-m))^2/(1+(4pialpha(t+(m*-m)))^2)^2 * w(t)^2
	p[k]+=amp2*amp2*wsq;
      }
      e[k]*=scl/dns;
      p[k]*=0.5*dt*alpha*alpha/dns;
      jval[k]=e[k]+p[k];
      //      cout<<"k="<<k<<" jval="<<jval[k]<<endl;
    } 

    /* write jvpm data file */
    fwrite(jval,sizeof(float),nm,fpd);
    fwrite(e,sizeof(float),nm,fpd);
    fwrite(p,sizeof(float),nm,fpd);

    delete [] jval;
    delete [] e;
    delete [] p;
    
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
  
