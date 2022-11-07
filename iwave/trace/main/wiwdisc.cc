#include <su.h>
#include <segy.h>
#include <header.h>
#include <parser.h>

using RVL::valparse;
using RVL::RVLException;

#define TOL 0.001

const char * sdoc[] = { 
  "Usage: wiwdisc.x data= out= maxit= mmin= mmax= nm=",
  "",
  "Synopsis: Implements discrepancy principle updates of alpha",
  "and m using formulas (38)-(41) in \"Solution of an Acoustic",
  "Transmission Inverse Problem by Source Extension\". ",
  "The source-receiver distance r is fixed at 1 km.",
  "Output is a record of the update loop. Slowness update by",
  "minimizing VPM objective by examining all values sampled over the",
  "interval [mmin, mmax] at sample interval dm = (mmax-mmin)/(nm-1). The",
  "data is read from an SU file.",
  ""
  "Penalty weight alpha computed to place residual between eminus and ",
  "eplus, on either side of etarget. etarget = 1/(2(*n*sn), eminus = ",
  "etarget*efactor, eplus = etarget/efactor, and efactor < 1 is fudge ",
  "factor providing nonzero tolerance for data misfit.",
  "",
  "Required inputs:",
  "  data [string]     = input data trace",
  "  res [string]      = output residual trace",
  "  out [string]      = output text file",
  "  maxit             = max iterations (alpha and m)",
  "  maxalpha          = max iterations (alpha, given m)",
  //  "  maxm              = max iterations (m, given alpha)",
  "  mmin [float]      = min slowness",
  "  mmax [float]      = max slowness",
  "  nm [int]          = number of slowness samples",
  "  sn [float]        = target signal-to-noise ratio",
  "  efactor [float]   = fudge factor < 1",
  "  rule [string]     = update rule",
  "",
  "NB: working units are s, km. File units must be scaled: SEGY stores dt in",
  "microsec, delrt (= t0) in ms.",
  NULL};

int main(int argc, char ** argv) {

  try {
  
    PARARRAY * par;  /* param array */
 
    std::string data;/* data file name */
    std::string res; /* residual file name */
    std::string out; /* output file name */
    float alpha;     /* penalty weight */
    FILE * fpd;      /* data file pointer */
    FILE * fpr;      /* res file pointer */
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
    int maxit;       /* max number of m/alpha updates */
    int maxalpha;    /* max number of alpha updates */
    //    int maxm;        /* max number of m updates */
    float sn;        /* target s/n */
    float etarget;   /* 1/2 square of n/s */
    float eplus;     /* upper 1/2 ms residual bound for alpha update */
    float eminus;    /* lower 1/2 ms residual bound for alpha update */
    float efactor;   /* factor for eplus, eminus vs etarget */
    string rule;     /* update rule */

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
  
    data     = RVL::valparse<std::string>(*par,"data");
    res      = RVL::valparse<std::string>(*par,"res");
    out      = RVL::valparse<std::string>(*par,"out");
    maxit    = RVL::valparse<int>(*par,"maxit");
    maxalpha = RVL::valparse<int>(*par,"maxalpha");
    //    maxm     = RVL::valparse<int>(*par,"maxm");
    mmin     = RVL::valparse<float>(*par,"mmin");
    mmax     = RVL::valparse<float>(*par,"mmax");
    nm       = RVL::valparse<int>(*par,"nm");
    sn       = RVL::valparse<float>(*par,"sn",0.0f);
    efactor  = RVL::valparse<float>(*par,"efactor",0.9f);
    rule     = RVL::valparse<std::string>(*par,"rule","pbound");

    std::ofstream outf;
    outf.open(out);

    if ((efactor<0.01) || (efactor>0.99f)) {
      RVLException e;
      e<<"ERROR: wiwdisc\n";
      e<<"  efactor="<<efactor<<" should be in [0.01, 0.99]\n";
      throw e;
    }
    etarget = 1/(2.0*sn*sn);
    eplus = etarget/efactor;
    eminus = etarget*efactor;
    outf<<"rule="<<rule<<endl;      
    outf<<"s/n="<<sn<<endl;
    outf<<"etarget="<<etarget<<endl;
    outf<<"eplus="<<eplus<<endl;
    outf<<"eminus="<<eminus<<endl;
    
    if (nm<2) {
      RVL::RVLException e;
      e<<"ERROR: wiwdisc\n";
      e<<"  nm =  "<<nm<<" must be at least 2\n";
      throw e;
    }
    dm = (mmax-mmin)/(nm-1);

    /* open data file */
    if (!(fpd=fopen(data.c_str(),"r"))) {
      RVL::RVLException e;
      e<<"ERROR: wiw\n";
      e<<"  failed to open input file "<<data<<"\n";
      throw e;
    }

    /* read data */
    if (!(fgettr(fpd,&tr))) {
      RVL::RVLException e;
      e<<"ERROR: wiw\n";
      e<<"  failed to read trace from "<<data<<"\n";
      throw e;
    }

    fclose(fpd);

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

    /* workspace for residual */
    float pi = 3.1415927;
    float fac = (4.0*pi)*(4.0*pi);
    float dns = 0.0f;
    for (int i=0;i<nt;i++) {
      dns += 0.5*dt*(tr.data[i])*(tr.data[i]);
    }
    /* equations 38-40 in WIW paper */
    /* loop over slowness */

    // initialize alpha
    alpha=0.0f;
    float eoptc=0.0f;
    float alphasq=0.0f;
    float alphasqc=0.0f;
    float ehi=0.0f;
    float elo=0.0f;
    float alphahi=0.0f;
    float alphalo=0.0f;
    float mopt=mmin;
    float eopt=0.0;
    float popt=0.0;
    float jopt=0.0;

    int termcode=-1;
    int it = 0;
    
    while ((it<maxit) && (termcode<0)) {
      
      outf<<"-------------------------------"<<endl;
      outf<<"  iteration="<<it<<endl;
      if (termcode > 0) {
	outf<<"  termination: termcode="<<termcode<<endl;
      }
      outf<<"-------------------------------"<<endl;
      outf<<"  ---- update alpha"<<endl;
      outf<<"-------------------------------"<<endl;

      alpha=0.0f;

      for (int a=0; a<maxalpha; a++) {

	alphasq = alpha*alpha;
	
	// (re)compute j, e, p
	eopt=0.0f;
	popt=0.0f;
	for (int i=0;i<nt;i++) {
	  /* scale times to s */
	  float t = t0 + i*dt;
	  // (4pi(t-mopt))^2
	  float amp0 = alphasq*fac*(t-mopt)*(t-mopt);
	  float amp1 = amp0/(1.0f + amp0);
	  float amp2 = fac*(t-mopt)*(t-mopt)/((1.0f + amp0)*(1.0f + amp0));
	  float wsq  = (tr.data[i])*(tr.data[i]);
	  // accumulate (4pialpha(t-m))^4/(1+(4pialpha(t-m))^2)^2
	  // * d(t)^2
	  eopt+=amp1*amp1*wsq;
	  // accumulate (t-m)^2/(1+(4pialpha(t-m))^2)^2 * d(t)^2
	  popt+=amp2*wsq;
	}
	eopt*=0.5*dt/dns;
	popt*=0.5*dt/dns;
	jopt=eopt+alphasq*popt;

	// alpha update tree
	outf<<"-------------------------------"<<endl;
	outf<<"alpha update # "<<a<<endl;
	outf<<"alpha="<<alpha<<endl;
	outf<<"j="<<jopt<<endl;
	outf<<"e="<<eopt<<endl;
	outf<<"p="<<popt<<endl;
	outf<<"m="<<mopt<<endl;
	// successful termination
	if ((eopt >= eminus) && (eopt <= eplus)) {
	  outf<<"-------------------------------"<<endl;
	  outf<<"bounds test passed"<<endl;
	  outf<<"mest="<<mopt<<endl;
	  outf<<"eest="<<eopt<<endl;
	  outf<<"rest="<<sqrt(2.0*eopt)<<endl;
	  outf<<"aest="<<alpha<<endl;
	  outf<<"-------------------------------"<<endl;
	  termcode=0;
	  break;
	}
	else if (rule=="pbound") {
	  if (eopt > eplus) {
	    //	    outf<<"eopt="<<eopt<<" > eplus="<<eplus<<": should not happen!\n";
	    termcode=1;
	    outf<<"-------------------------------"<<endl;
	    outf<<"termcode=1"<<endl;
	    outf<<"-------------------------------"<<endl;	    
	    break;
	  }
	  else if (eopt < eminus) {
	    alpha = sqrt(alphasq + (eplus-eopt)/(2.0f*popt));
	  }
	}
	else if (rule=="fixedsecant") {
	  alpha = sqrt(alphasq*etarget/eopt);
	}
	else if (rule=="secant") {
	  if (a==0) {
	    outf<<"first step by p bound"<<endl;
	    outf<<"eminus="<<eminus<<" eplus="<<eplus<<endl;
	    //	    outf<<"eopt="<<eopt<<" old elo="<<elo<<" old ehi="<<ehi<<endl;
	    alpha = sqrt(alphasq + (eplus-eopt)/(2.0f*popt));
	    alphalo=sqrt(alphasq);
	    elo=eopt;
	    alphahi=alphalo;
	    ehi=elo;
	    outf<<"alphalo="<<alphalo<<" alphahi="<<alphahi<<endl;
	    outf<<"eest="<<eopt<<" elo="<<elo<<" ehi="<<ehi<<endl;
	  }
	  else {
	    // increase alpha by power law rule until it exceeds plus
	    if (eopt<eminus && ehi<eminus) {
	      outf<<"power law step"<<endl;
	      outf<<"eminus="<<eminus<<" eplus="<<eplus<<endl;
	      outf<<"eest="<<eopt<<" old elo="<<elo<<" old ehi="<<ehi<<endl;
	      alpha = 2.0*sqrt(alphasq*sqrt(etarget/eopt));
      	      alphasqc = alphasq;
	      eoptc = eopt;
	      elo=ehi;
	      alphalo=alphahi;
	      ehi=eopt;
	      alphahi=sqrt(alphasq);
	      outf<<"elo="<<elo<<" ehi="<<ehi<<endl;
	      outf<<"alphalo="<<alphalo<<" alphahi="<<alphahi<<endl;
	    }
	    // if overshoot eplus, then eoptc<eminus so take
	    // first secant step aiming for etarget
	    // once ehi > eplus, will never be less
	    else if (eopt>eplus && ehi<eminus) {
	      outf<<"first secant step"<<endl;
	      outf<<"eminus="<<eminus<<" eplus="<<eplus<<endl;
	      outf<<"eest="<<eopt<<" old elo="<<elo<<" old ehi="<<ehi<<endl;  
	      alpha = sqrt(alphasqc
			   + (alphasq-alphasqc)*(etarget-eoptc)/(eopt-eoptc));
	      alphasqc = alphasq;
	      eoptc = eopt;
	      elo=ehi;
	      alphalo=alphahi;
	      ehi=eopt;
	      alphahi=sqrt(alphasq);
	      outf<<"elo="<<elo<<" ehi="<<ehi<<endl;
	      outf<<"alphalo="<<alphalo<<" alphahi="<<alphahi<<endl;
	    }
	    // subsequent secant steps keep elo<eminus, ehi>eplus
	    else {
	      outf<<"secant step"<<endl;
	      outf<<"eminus="<<eminus<<" eplus="<<eplus<<endl;
	      outf<<"eest="<<eopt<<" old elo="<<elo<<" old ehi="<<ehi<<endl;  
	      if ((eopt<eminus) || (eopt>eplus)) {
		if (eopt<eminus) {
		  elo=eopt;
		  alphalo=sqrt(alphasq);
		}
		else {
		  ehi=eopt;
		  alphahi=sqrt(alphasq);
		}
   		alpha = sqrt(sqrt(alphalo*alphalo*alphalo*alphalo 
		   + (alphahi*alphahi*alphahi*alphahi-
		      alphalo*alphalo*alphalo*alphalo)*
				  (etarget-elo)/(ehi-elo)));
		//   		alpha = sqrt(alphalo*alphalo 
		//		   + (alphahi*alphahi-alphalo*alphalo)*
		//		     (etarget-elo)/(ehi-elo));
      		//alpha = alphalo 
		//   + (alphahi-alphalo)*
		//     (etarget-elo)/(ehi-elo);
		alphasqc = alphasq;
		eoptc = eopt;
		outf<<"secant step"<<endl;
		outf<<"elo="<<elo<<" ehi="<<ehi<<endl;
		outf<<"alphalo="<<alphalo<<" alphahi="<<alphahi<<endl;
	      }
	      else {
		outf<<"-------------------------------"<<endl;
		outf<<"mest="<<mopt<<endl;
		outf<<"eest="<<eopt<<endl;
		outf<<"rest="<<sqrt(2.0*eopt)<<endl;
		outf<<"aest="<<alpha<<endl;
		outf<<"-------------------------------"<<endl;
		break;
	      }
	    }
	  }
	}
	else if (rule=="bisect") {
	  if (a==0) {
	    alpha = sqrt(alphasq + (eplus-eopt)/(2.0f*popt));
	  }
	  else {
	    if (eopt<eminus) alpha *= 2.0;
	    else if (eopt>eplus) alpha /= 1.8;
	    else {
	      outf<<"-------------------------------"<<endl;
	      outf<<"mest="<<mopt<<endl;
	      outf<<"eest="<<eopt<<endl;
	      outf<<"rest="<<sqrt(2.0*eopt)<<endl;
	      outf<<"aest="<<alpha<<endl;
	      outf<<"-------------------------------"<<endl;
	    }
	  }
	}
	else {
	  RVLException e;
	  e<<"ERROR: wiwdisc\n";
	  e<<"  rule = "<<rule<<" not known update rule\n";
	  throw e;
	}
      
        if (a==maxalpha-1) {
	  termcode=2;
	  outf<<"-------------------------------"<<endl;
	  outf<<"termcode=2"<<endl;
	  outf<<"-------------------------------"<<endl;
	}
      } // end e tests

      // failure to update alpha
      if (termcode>0) {
	outf<<"-------------------------------"<<endl;
	outf<<"  alpha update failed, terminate algorithm"<<endl;
	outf<<"-------------------------------"<<endl;
      }

      // alpha update successful, proceed to m update
      if (termcode==0) {
	outf<<"-------------------------------"<<endl;
	outf<<"  ---- update m"<<endl;
	outf<<"-------------------------------"<<endl;
	
	eopt=std::numeric_limits<float>::max();
	popt=std::numeric_limits<float>::max();
	jopt=std::numeric_limits<float>::max();
	alphasq = alpha*alpha;
	
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
	} // end m search loop
	
	outf<<"-------------------------------"<<endl;
	outf<<"mest="<<mopt<<endl;
	outf<<"-------------------------------"<<endl;
	
      }// end m, alpha update loop

      if ((eopt >= eminus) && (eopt <= eplus)) {
         outf<<"-------------------------------"<<endl;
	 outf<<"bounds test passed for optimal model"<<endl;
	 outf<<"mopt="<<mopt<<endl;
	 outf<<"eopt="<<eopt<<endl;
	 outf<<"ropt="<<sqrt(2.0*eopt)<<endl;
	 outf<<"aopt="<<alpha<<endl;
	 outf<<"terminate"<<endl;
	 outf<<"-------------------------------"<<endl;
	 break;
      }
      else {
	termcode=-1;
	it++;
      }
    }

    
    if (termcode==0) {
      // write residual
      for (int i=0;i<nt;i++) {
	/* scale times to s */
	float t = t0 + i*dt;
	// (4pi(t-m))^2
	float amp0 = alpha*alpha*fac*(t-mopt)*(t-mopt);
	tr.data[i] *= -amp0/(1.0f + amp0);
      }
      
      /* open residual file */
      if (!(fpr=fopen(res.c_str(),"w"))) {
	RVL::RVLException e;
	e<<"ERROR: wiw\n";
	e<<"  failed to open input file "<<data<<"\n";
	throw e;
      }
      
      /* write residual */
      fputtr(fpr,&tr);
      
      fclose(fpr);
    }
    
    delete [] jval;
    delete [] e;
    delete [] p;
    delete [] ap;
    
    /* close files */
    outf.flush();
    outf.close();
    ps_delete(&par);
    
  }
  catch (RVL::RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
  
