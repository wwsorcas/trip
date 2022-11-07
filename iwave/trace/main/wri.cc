#include <su.h>
#include <segy.h>
#include <header.h>
#include <parser.h>

using RVL::valparse;
using RVL::RVLException;

#define TOL 0.001

const char * sdoc[] = { 
  "Usage: wri.x wavelet= residual= jvpm= alpha= mex= mmin= mmax= nm= offset=",
  "",
  "Synopsis: Implements formulas (13), (35), and (36) in \"Wavefield ",
  "reconstruction inversion: an example\". Output is fwi and wri objectives",
  "sampled over the interval [mmin, mmax] at sample interval dm = ",
  "(mmax-mmin)/(nm-1), and recorded in two columns (n2=2). The",
  "wavelet and residual traces are read from SU files, the output",
  "samples are organized in an RSF file pair.",
  "",
  "Required inputs:",
  "  wavelet [string]  = input source wavelet filename (.su)",
  "  data [string]     = output data trace filename (.su)",
  "  residual [string] = output residual trace filename (.su)",
  "  jvpm [string]     = output jvpm samples, filename (.rsf)",
  "  alpha [float]     = penalty parameter for WRI",
  "  mex [float]       = target (\"exact\") slowness",
  "  mmin [float]      = min slowness",
  "  mmax [float]      = max slowness",
  "  nm [int]          = number of slowness samples",
  "  offset [float]    = source-receiver distance",
  "",
  "NB: working units are s, km. File units must be scaled: SEGY stores dt in",
  "microsec, delrt (= t0) in ms.",
  NULL};

int main(int argc, char ** argv) {

  try {
  
    PARARRAY * par;    /* param array */
 
    std::string wav;      /* wavelet file name */
    std::string res;      /* residual trace file name */
    std::string dat;      /* data trace file name */
    std::string jvpm;     /* output file name */
    float alpha;     /* */
    FILE * fpw;      /* wavelet file pointer */
    FILE * fpr;      /* residual file pointer */
    FILE * fps;      /* data file pointer */
    FILE * fph;      /* jvpm header file pointer */
    FILE * fpd;      /* jvpm data file pointer */
    segy trw;        /* wavelet trace workspace */
    segy trr;        /* residual/data trace workspace */
    Value val;       /* header word workspace */
    float dtw;       /* time step */
    float t0w;       /* start time of wavelet series */
    int ntr;         /* number of samples per residual trace */
    float dtr;       /* time step */
    float t0r;       /* start time of residual series */
    int ntw;         /* number of samples per trace */
    float dm;        /* slowness interval */
    float mmin;      /* min slowness */
    float mmax;      /* max slowness */
    float mex;       /* target ("exact") slowness */
    int nm;          /* number of slowness samples */
    float * fwi;     /* array of fwi objective values */
    float * wri;     /* array of wri objective values */
    float * r;       /* residual array */
    float * w;       /* wavelet array */
    float m;         /* trial slowness */
    float t;         /* time */
    float tint;      /* shift time */
    int jint;        /* shift index */
    float frac;      /* interpolation time, normalized */
    float offs;      /* source-receiver offset */

    
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
    res    = RVL::valparse<std::string>(*par,"residual");
    dat    = RVL::valparse<std::string>(*par,"data");
    jvpm   = RVL::valparse<std::string>(*par,"jvpm");
    alpha  = RVL::valparse<float>(*par,"alpha");
    mex    = RVL::valparse<float>(*par,"mex");
    mmin   = RVL::valparse<float>(*par,"mmin");
    mmax   = RVL::valparse<float>(*par,"mmax");
    nm     = RVL::valparse<int>(*par,"nm");
    offs   = RVL::valparse<float>(*par,"offset");
      
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

    if (!(fpr=fopen(res.c_str(),"r+"))) {
      RVL::RVLException e;
      e<<"ERROR: wiw\n";
      e<<"  failed to open input file "<<res<<"\n";
      throw e;
    }

    if (!(fps=fopen(dat.c_str(),"w"))) {
      RVL::RVLException e;
      e<<"ERROR: wiw\n";
      e<<"  failed to open input file "<<dat<<"\n";
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
    if (!(fgettr(fpw,&trw))) {
      RVL::RVLException e;
      e<<"ERROR: wiw\n";
      e<<"  failed to read trace from "<<wav<<"\n";
      throw e;
    }

    /* read residual - at this point, only header is significant */
    if (!(fgettr(fpr,&trr))) {
      RVL::RVLException e;
      e<<"ERROR: wiw\n";
      e<<"  failed to read trace from "<<wav<<"\n";
      throw e;
    }

    /* workspace to make up for gethdval not taking const args */
    char * str = new char[128];
    
    /* read time sample info */
    strcpy(str,"ns");
    gethdval(&trw,str,&val);
    ntw = vtoi(hdtype(str),val);

    /* file dt in \mu s; scale to s */
    strcpy(str,"dt");
    gethdval(&trw,str,&val);
    dtw = 0.000001*vtof(hdtype(str),val);

    /* wavelet file t0 in ms; scale to s */
    strcpy(str,"delrt");
    gethdval(&trw,str,&val);
    t0w = 0.001*vtof(hdtype(str),val);

    /* residual file t0 in ms; scale to s */
    strcpy(str,"delrt");
    gethdval(&trr,str,&val);
    t0r = 0.001*vtof(hdtype(str),val);
    delete [] str;
    
    /* read time sample info */
    strcpy(str,"ns");
    gethdval(&trr,str,&val);
    ntr = vtoi(hdtype(str),val);

    /* file dt in \mu s; scale to s */
    strcpy(str,"dt");
    gethdval(&trr,str,&val);
    dtr = 0.000001*vtof(hdtype(str),val);

    /* workspace for jvpm values */
    fwi = new float[nm];
    wri = new float[nm];
    r   = trr.data;
    w   = trw.data;
    
    /* equation 18 in WIW paper */
    /* loop over slowness */
    for (int k=0; k<nm; k++) {
      m = mmin + k*dm;
      fwi[k]=0.0f;
      wri[k]=0.0f;
      /* zero residual trace */
      for (int i=0;i<ntr;i++) {
	r[i]=0;
      }
      for (int i=0;i<ntr;i++) {      
	t = t0r + i*dtr;
	/* shift wavelet by mex*r (target data) */
        tint = (t-mex*offs-t0w)/dtw;
	jint = tint;
	frac = tint-jint;
	if ((jint >= 0) && (jint<ntw-1)) {
	  r[i] = 0.5*mex*(w[jint+1]*frac + w[jint]*(1.0f-frac));
	}
      }
      /* first time through, write out data trace */
      if (k==0) {
	fputtr(fps,&trr);
      }
      for (int i=0;i<ntr;i++) {      
	t = t0r + i*dtr;      
	/* shift wavelet by m*r (predicted data) */
        tint = (t-m*offs-t0w)/dtw;
	jint = tint;
	frac = tint-jint;
	if ((jint >= 0) && (jint<ntw-1)) {
	  r[i] -= 0.5*m*(w[jint+1]*frac + w[jint]*(1.0f-frac));
	}
	fwi[k]+=0.5*dtr*r[i]*r[i];
      }
      // equations 35, 36
      wri[k]=0.5*alpha*alpha*fwi[k]/(0.25*m*m*offs + alpha*alpha);
      /* write residual trace */
      fputtr(fpr,&trr);
    } 

    /* write fwi, wri data files */
    fwrite(fwi,sizeof(float),nm,fpd);
    fwrite(wri,sizeof(float),nm,fpd);

    delete [] fwi;
    delete [] wri;
    
    /* write jvpm header file */
    fprintf(fph,"n1=%d ",nm);
    fprintf(fph,"d1=%g ",dm);
    fprintf(fph,"o1=%g\n",mmin);
    fprintf(fph,"n2=%d ",2);
    fprintf(fph,"d2=%g ",1.0f);
    fprintf(fph,"o2=%g\n",0.0f);    
    fprintf(fph,"esize=4 data_format=\"native_float\"\n");
    fprintf(fph,"in=\"%s\"\n",jvpd.c_str());

    /* close files */
    fclose(fpw);
    fclose(fpr);
    fclose(fps);
    fclose(fph);
    fclose(fpd);
    ps_delete(&par);

  }
  catch (RVL::RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
  
