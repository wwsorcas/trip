#include "fftops.hh"

//*******************************************************************************************************

/* Weight vector: v=exp(-i*2PI/2N) */
int ffttrig(int nfft,float_complex *trigs)
{
  int i;
  float_complex cc;
  float dinc=2.f*acos(-1.0)/nfft;
  // WWS
  float_complex sm1(0.0f,1.0f);
  for(i=0;i<nfft;i++)
    {
      cc=dinc*i;
      //      trigs[i]=cexp(I*cc);
      trigs[i]=std::exp(sm1*cc);      
    }
  return 0;
}

//*******************************************************************************************************

/* DFT function (decimation in time, base 2) */
void fft(float_complex *a,float_complex *b,float_complex *trigs,int nfft,int isign)
{
  int la=1,ii=1,i,j,m,k,l;
  float_complex t;

  m=nfft/2;
  while(la<nfft)
    {
      i=0;
      j=0;
      if(ii>0)
	{
	  for(k=0;k<=m-la;k+=la)
	    {
	      t=trigs[k];
	      if(isign<0)
		t=conj(t);
	      for(l=1;l<=la;l++)
		{
		  b[j]=a[i]+a[i+m];
		  b[j+la]=t*(a[i]-a[i+m]);
		  i++;
		  j++;
		}
	      j+=la;
	    }
	}
      else
	{
	  for(k=0;k<=m-la;k+=la)
	    {
	      t=trigs[k];
	      if(isign<0)
		t=conj(t);
	      for(l=1;l<=la;l++)
		{
		  a[j]=b[i]+b[i+m];
		  a[j+la]=t*(b[i]-b[i+m]);
		  i++;
		  j++;
		}
	      j+=la;
	    }
	}
      ii=-ii;
      la*=2;
    }
    
  if(ii<0)
    {
      for(j=0;j<nfft;j++)
	a[j]=b[j];
    }
}

//*******************************************************************************************************

/* Multiplication by |freq|^power */
void MultFreq(float *Trace, int Nt, float dt, float_complex *Aux, int nfft, int power, float *band) {

	int   i, isign;
	int   if1, if2, if3, if4;
	float df, arg, taper;
	float fact;
	float_complex *a, *b, *trigs;
	
	a     = &Aux[0];
	b     = &Aux[nfft];
	trigs = &Aux[2*nfft];
    
    //For band-pass filter
   	float PI  = acos (-1.0);
	df  = 2.f*PI / dt / nfft;
	if1 = floor(2.f*PI * band[0] / df);
	if2 = floor(2.f*PI * band[1] / df);
	if3 = ceil (2.f*PI * band[2] / df);
	if4 = ceil (2.f*PI * band[3] / df);

	//Copy into a complex array
	for(i = 0; i < Nt; i++)
		a[i] = Trace[i];

   	//Forward complex DFT
	ffttrig(nfft, trigs);
	isign = -1;
	fft(a, b, trigs, nfft, isign);

	//Loop over all frequencies in the band, scale by (I*f)^power
	for(i = 1; i < nfft/2; i++) {
			
		fact = powf(i*df, power);
		if ((i < if1) || (i >= if4)) {
			a[i]      = (float_complex)0.f;
			a[nfft-i] = (float_complex)0.f;
		} else if ((i >= if1) && (i < if2)) {
			arg   = (float) (if2-i)/(if2-if1) * PI/2.f;
			taper = pow(cosf(arg), 2);
			a[i]      *= fact * taper;
			a[nfft-i] *= fact * taper;
		} else if ((i >= if2) && (i < if3)) {
			a[i]      *= fact;
			a[nfft-i] *= fact;
		} else if ((i >= if3) && (i < if4)) {
			arg = (float) (i-if3)/(if4-if3) * PI/2.f;
			taper = pow(cosf(arg), 2);
			a[i]      *= fact * taper;
			a[nfft-i] *= fact * taper;
		}

	}
	a[0]      = (float_complex)0.f;
	a[nfft/2] = (float_complex)0.f;

	//Inverse complex DFT
	isign = 1;
	fft(a, b, trigs, nfft, isign);
	
	//Copy back into output trace
	for(i = 0; i < Nt; i++)
	  //		Trace[i] = creal(a[i]) / nfft;
	  Trace[i] = std::real(a[i]) / nfft;	  

}

//*******************************************************************************************************

/* Differentiation / Integration - Multiplication by (I*freq)^power */
void DiffInt (float *Trace, int Nt, float dt, float_complex *Aux, int nfft, int power, float *band) {

	int   i, isign;
	int   if1, if2, if3, if4;
	float df, arg, taper;
	float_complex fact1, fact2;
	float_complex *a, *b, *trigs;

	// WWS
	float_complex sm1(0.0f,1.0f);
	
	a     = &Aux[0];
	b     = &Aux[nfft];
	trigs = &Aux[2*nfft];
    
    //For band-pass filter
   	float PI  = acos (-1.0);
	df  = 2.f*PI / dt / nfft;
	if1 = floor(2.f*PI * band[0] / df);
	if2 = floor(2.f*PI * band[1] / df);
	if3 = ceil (2.f*PI * band[2] / df);
	if4 = ceil (2.f*PI * band[3] / df);

	//Copy into a complex array
	for(i = 0; i < Nt; i++)
		a[i] = Trace[i];

   	//Forward complex DFT
	ffttrig(nfft, trigs);
	isign = -1;
	fft(a, b, trigs, nfft, isign);

	//Loop over all frequencies in the band, scale by (I*f)^power
	for(i = 1; i < nfft/2; i++) {
	    //		fact1 = cpowf(i*df*I, power);
	  fact1 = std::pow(i*df*sm1, power);
	  //		fact2 = cpowf(-i*df*I, power);
	  fact2 = std::pow(-i*df*sm1, power);		
		if ((i < if1) || (i >= if4)) {
			a[i]      = (float_complex)0.f;
			a[nfft-i] = (float_complex)0.f;
		} else if ((i >= if1) && (i < if2)) {
			arg   = (float) (if2-i)/(if2-if1) * PI/2.f;
			taper = powf(cosf(arg), 2);
			a[i]      *= fact1 * taper;
			a[nfft-i] *= fact2 * taper;
		} else if ((i >= if2) && (i < if3)) {
			a[i]      *= fact1;
			a[nfft-i] *= fact2;
		} else if ((i >= if3) && (i < if4)) {
			arg = (float) (i-if3)/(if4-if3) * PI/2.f;
			taper = powf(cosf(arg), 2);
			a[i]      *= fact1 * taper;
			a[nfft-i] *= fact2 * taper;
		}

	}
	a[0]      = (float_complex)0.f;
	a[nfft/2] = (float_complex)0.f;

	//Inverse complex DFT
	isign = 1;
	fft(a, b, trigs, nfft, isign);
	
	//Copy back into output trace
	for(i = 0; i < Nt; i++)
	  //		Trace[i] = creal(a[i]) / nfft;
	  Trace[i] = std::real(a[i]) / nfft;	

}




//*******************************************************************************************************


namespace TSOpt {

  void FTScale(float *Trace, int Nt, float dt,
	       int ab, int power, float *band) {
    int nfft=1;
    while (nfft < Nt) { nfft *= 2; }
    nfft *= 2;
    float_complex * Aux =
      (float_complex*)calloc(3*nfft,sizeof(float_complex));
    if (ab)    
      MultFreq(Trace, Nt, dt,
	       Aux, nfft, power, band);
    else
      DiffInt(Trace, Nt, dt,
	      Aux, nfft, power, band);
    free(Aux);
  }
  
  void FTScaleFO::operator()(RVL::LocalDataContainer<float> & x) {

    int Nt = x.getSize();
    float * Trace = x.getData();
    FTScale(Trace, Nt, dt, ab, power, band);
  }

}
