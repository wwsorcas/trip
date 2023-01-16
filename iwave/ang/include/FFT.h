//#######################################################################################################
// HEADER:	FFT.h
// DESCRIPTION:	Fourier Transform functions (real and complex)
// CONTENT:	DFT (Discrete complex Fourier Transform)
//		ffttrig Function of the weight vector: v=exp(-i*2PI/2N)
//		Real FFT (Real to Complex)
//		Real INVERSE FFT (Complex to Real) 
//#######################################################################################################

#define  PI 3.141592654

#include <complex.h>

/* FFT auxiliary arrays */
typedef struct {

	int   n1_nfft;
	int   n2_nfft;
	float complex ***M;
	float complex *a;
	float complex *b;
	float complex *c;
	float complex *d;
	float complex *trigs1;
	float complex *trigs2;
	float complex *trigs3;
	float *k1;
	float *k2;

} stFFT;

//*******************************************************************************************************
//The FFT function (decimation in time)
void fft(float complex *a,float complex *b,float complex *trigs,int nfft,int isign)
{
  int la=1,ii=1,i,j,m,k,l;
  float complex t;

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

//Function of the weight vector: v=exp(-i*2PI/2N)
int ffttrig(int nfft,float complex *trigs)
{
  int i;
  float complex cc;
  float dinc=2.f*acos(-1.0)/nfft;
  for(i=0;i<nfft;i++)
    {
      cc=dinc*i;
      trigs[i]=cexp(I*cc);
    }
  return 0;
}


//*******************************************************************************************************

//The real FFT function (real to complex)
void rfftrc(float complex *a, float complex *b, float complex *trigs, float complex *trigs1, int nfft)
{
  int i,n,isign;
  n=nfft/2;
  isign=-1;
  float complex A;

  fft(a,b,trigs,n,isign);

  a[n]=a[0];
  for(i=0;i<=n;i++)
    {
      A=conj(a[n-i]);
      b[i]=a[i]+A-I*conj(trigs1[i])*(a[i]-A);
    }
  for(i=0;i<=n;i++)
    a[i]=(float complex)0.5*b[i];
}


//*******************************************************************************************************

//The real FFT function (complex to real)
void rfftcr(float complex *a, float complex *b, float complex *trigs, float complex *trigs1, int nfft)
{
   int i,n=nfft/2,isign=1;
   float con=1./nfft;
   float complex A;

   for(i=0;i<n;i++)
     {
       A=conj(a[n-i]);
       b[i]=a[i]+A+I*trigs1[i]*(a[i]-A);
     }
   for(i=0;i<n;i++)
     a[i]=(float complex)con*b[i];

   fft(a,b,trigs,n,isign);
}


//*******************************************************************************************************

//Allocate auxiliary arrays for FFTs functions
void AllocAuxFFT(stFFT *pFFT , int n3 , int n2 , int n1)
{
	int i;

	//Calculate n1_nfft and n2_nfft size (power of 2)
    for(i = 0 ; pow(2, i) < n1 ; i++)
      ;
    pFFT->n1_nfft = (int)pow(2, i);
    for(i = 0 ; pow(2 , i) < n2 ; i++)
      ;
    pFFT->n2_nfft=(int)pow(2, i);
    //Double the model size by zero-padding - avoid rap-around effect 
    pFFT->n1_nfft *= 2;
    pFFT->n2_nfft *= 2;
    
    //Memory allocation
    pFFT->M = alloc_complex_cube (n3, pFFT->n2_nfft, pFFT->n1_nfft);
    pFFT->a = alloc_complex_array(pFFT->n1_nfft);
    pFFT->b = alloc_complex_array(pFFT->n1_nfft);
    pFFT->c = alloc_complex_array(pFFT->n2_nfft);
    pFFT->d = alloc_complex_array(pFFT->n2_nfft);
    pFFT->trigs1 = alloc_complex_array(pFFT->n1_nfft);
    pFFT->trigs2 = alloc_complex_array(pFFT->n1_nfft/2);
    pFFT->trigs3 = alloc_complex_array(pFFT->n2_nfft);
	
    //frequency allocation
    pFFT->k1 = alloc_array(pFFT->n1_nfft/2);
    pFFT->k2 = alloc_array(pFFT->n2_nfft);
    
}


//*******************************************************************************************************

//Zero-out auxiliary arrays of FFTs functions
void ZeroOutAuxFFT(stFFT *pFFT  , int n3)
{
	int i , m;

	//Zero-out arrays
	for(i = 0 ; i < n3 ; i++)
		for(m = 0 ; m < pFFT->n2_nfft ; m++)
			memset(pFFT->M[i][m] , 0. , pFFT->n1_nfft*sizeof(float complex));
	memset(pFFT->a , 0. , pFFT->n1_nfft*sizeof(float complex));
	memset(pFFT->b , 0. , pFFT->n1_nfft*sizeof(float complex));
	memset(pFFT->c , 0. , pFFT->n2_nfft*sizeof(float complex));
	memset(pFFT->d , 0. , pFFT->n2_nfft*sizeof(float complex));
	memset(pFFT->trigs1 , 0. , pFFT->n1_nfft*sizeof(float complex));
	memset(pFFT->trigs2 , 0. , pFFT->n1_nfft/2*sizeof(float complex));
	memset(pFFT->trigs3 , 0. , pFFT->n2_nfft*sizeof(float complex));
	
}


//*******************************************************************************************************

//Copy into a complex array for FFT
void CopyForFFT(stFFT *pFFT , float ***Vol , int n3 , int n2 , int n1)
{
	int i , j , m , n;

	for(i = 0 ; i < n3 ; i++)
		for(m = 0 ; m < n2 ; m++)
	   		for(n = 0 , j = 0 ; n < n1/2 ; n++, j+=2)
				pFFT->M[i][m][n] = Vol[i][m][j] + Vol[i][m][j+1]*I;

}


//*******************************************************************************************************

//Copy complex array of FFT back to original volume 
void CopyFromFFT(stFFT *pFFT , float ***Vol , int n3 , int n2 , int n1)
{
	int i , m , n;

	for(i = 0 ; i < n3 ; i++)
		for(m = 0 ; m < n2 ; m++)
	   		for(n = 0 ; n < n1 ; n++)
				Vol[i][m][n] = creal(pFFT->M[i][m][n]);

}

//*******************************************************************************************************

//Copy into a complex array for FFT (transpose axis 2 ans 3)
void CopyForFFT_transp(stFFT *pFFT , float ***Vol , int n3 , int n2 , int n1)
{
	int i , j , m , n;

	for(i = 0 ; i < n3 ; i++)
		for(m = 0 ; m < n2 ; m++)
	   		for(n = 0 , j = 0 ; n < n1/2 ; n++, j+=2)
				pFFT->M[m][i][n] = Vol[i][m][j] + Vol[i][m][j+1]*I;

}


//*******************************************************************************************************

//Copy complex array of FFT back to original volume (transpose axis 2 ans 3)
void CopyFromFFT_transp(stFFT *pFFT , float ***Vol , int n3 , int n2 , int n1)
{
	int i , m , n;

	for(i = 0 ; i < n3 ; i++)
		for(m = 0 ; m < n2 ; m++)
	   		for(n = 0 ; n < n1 ; n++)
				Vol[m][i][n] = creal(pFFT->M[i][m][n]);

}


//*******************************************************************************************************

//1D-RealFFT of volume projections
void FFT1DRealVol(stFFT *pFFT , int n3)
{
	int i , m , n;
	
	//Real FFT into Fourier space
    ffttrig(pFFT->n1_nfft   , pFFT->trigs1);
    ffttrig(pFFT->n1_nfft/2 , pFFT->trigs2);
	for(i = 0 ; i < n3 ; i++) {
    	for(m = 0 ; m < pFFT->n2_nfft ; m++) {
			for(n = 0 ; n < pFFT->n1_nfft/2 ; n++)
				pFFT->a[n] = pFFT->M[i][m][n];
			rfftrc(pFFT->a , pFFT->b , pFFT->trigs2 , pFFT->trigs1 , pFFT->n1_nfft);
			for(n = 0 ; n < pFFT->n1_nfft/2 ; n++)
				pFFT->M[i][m][n] = pFFT->a[n];
		}
	}
	
}


//*******************************************************************************************************

//1D-Inverse RealFFT of volume projections
void InvFFT1DRealVol(stFFT *pFFT , int n3)
{
	int i , j , m , n;

	//Real Inverse FFT from Fourier space into axis 1
    ffttrig(pFFT->n1_nfft   , pFFT->trigs1);
    ffttrig(pFFT->n1_nfft/2 , pFFT->trigs2);
	for(i = 0 ; i < n3 ; i++)  
    	for(m = 0 ; m < pFFT->n2_nfft ; m++) {
			for(n = 0 ; n < pFFT->n1_nfft/2 ; n++)
				pFFT->a[n] = pFFT->M[i][m][n];
			rfftcr(pFFT->a , pFFT->b , pFFT->trigs2 , pFFT->trigs1 , pFFT->n1_nfft);
			for(n = 0 , j = 0 ; n < pFFT->n1_nfft/2 ; n++, j+=2) {
				pFFT->M[i][m][j]   = creal(pFFT->a[n]);
				pFFT->M[i][m][j+1] = cimag(pFFT->a[n]);
			}
		}

}


//*******************************************************************************************************

//2D-FFT of volume projections
void FFT2DVol(stFFT *pFFT , int n3)
{
	int i , m , n , isign;

	//Real FFT of axis 1 into Fourier space
	FFT1DRealVol (pFFT , n3);

    //Complex DFT of axis 2 into Fourier space
    ffttrig(pFFT->n2_nfft , pFFT->trigs3);
    isign = -1;
	for(i = 0 ; i < n3 ; i++)
    	for(n = 0 ; n < pFFT->n1_nfft/2 ; n++) {
			for(m = 0 ; m < pFFT->n2_nfft ; m++)
				pFFT->c[m] = pFFT->M[i][m][n];
			fft(pFFT->c , pFFT->d , pFFT->trigs3 , pFFT->n2_nfft , isign);
			for(m = 0 ; m < pFFT->n2_nfft ; m++)
				pFFT->M[i][m][n] = pFFT->c[m];
		}

}


//*******************************************************************************************************

//2D-Inverse FFT of volume projections
void InvFFT2DVol(stFFT *pFFT , int n3)
{
	int i , m , n , isign;

	//Complex Inverse DFT from Fourier space into axis 2
   	ffttrig(pFFT->n2_nfft , pFFT->trigs3);
	isign = 1;
	for(i = 0 ; i < n3 ; i++)
		for(n = 0 ; n < pFFT->n1_nfft/2 ; n++) {
	    	for(m = 0 ; m < pFFT->n2_nfft ; m++)
				pFFT->c[m] = pFFT->M[i][m][n];
			fft(pFFT->c , pFFT->d , pFFT->trigs3 , pFFT->n2_nfft , isign);
			//Scale factor
			for(m = 0 ; m < pFFT->n2_nfft ; m++)
				pFFT->M[i][m][n] = pFFT->c[m] / pFFT->n2_nfft;
		}

	//Real Inverse FFT from Fourier space into axis 1
	InvFFT1DRealVol (pFFT , n3);

}


//*******************************************************************************************************

//3D-FFT of volume (Not3: n3 must be equal to pFFT->n2_nfft)
void FFT3DVol(stFFT *pFFT , int n3)
{
	int i , m , n , isign;

	//2D-FFT of volume projections
	FFT2DVol(pFFT , n3);
	
	//Zero-out arrays
	memset(pFFT->c , 0. , pFFT->n2_nfft*sizeof(float complex));
	memset(pFFT->d , 0. , pFFT->n2_nfft*sizeof(float complex));
	memset(pFFT->trigs3 , 0. , pFFT->n2_nfft*sizeof(float complex));
	
	//Complex DFT of axis 3 into Fourier space
    ffttrig(pFFT->n2_nfft , pFFT->trigs3);
    isign = -1;
   	for(n = 0 ; n < pFFT->n1_nfft/2 ; n++)
		for(m = 0 ; m < pFFT->n2_nfft ; m++) {
			for(i = 0 ; i < n3 ; i++)
				pFFT->c[i] = pFFT->M[i][m][n];
			fft(pFFT->c , pFFT->d , pFFT->trigs3 , pFFT->n2_nfft , isign);
			for(i = 0 ; i < n3 ; i++)
				pFFT->M[i][m][n] = pFFT->c[i];
		}
}


//*******************************************************************************************************

//3D-Inverse FFT of volume (Not3: n3 must be equal to pFFT->n2_nfft)
void InvFFT3DVol(stFFT *pFFT , int n3)
{
	int i , m , n , isign;

	//Zero-out arrays
	memset(pFFT->c , 0. , pFFT->n2_nfft*sizeof(float complex));
	memset(pFFT->d , 0. , pFFT->n2_nfft*sizeof(float complex));
	memset(pFFT->trigs3 , 0. , pFFT->n2_nfft*sizeof(float complex));

	//Complex Inverse DFT from Fourier space into axis 3
   	ffttrig(pFFT->n2_nfft , pFFT->trigs3);
	isign = 1;
	for(n = 0 ; n < pFFT->n1_nfft/2 ; n++)
	   	for(m = 0 ; m < pFFT->n2_nfft ; m++) {
	   		for(i = 0 ; i < n3 ; i++)
				pFFT->c[i] = pFFT->M[i][m][n];
			fft(pFFT->c , pFFT->d , pFFT->trigs3 , pFFT->n2_nfft , isign);
			//Scale factor
	   		for(i = 0 ; i < n3 ; i++)
				pFFT->M[i][m][n] = pFFT->c[i] / n3;
		}

	//Zero-out arrays
	memset(pFFT->c , 0. , pFFT->n2_nfft*sizeof(float complex));
	memset(pFFT->d , 0. , pFFT->n2_nfft*sizeof(float complex));
	memset(pFFT->trigs3 , 0. , pFFT->n2_nfft*sizeof(float complex));

	//2D-Inverse FFT of volume projections
	InvFFT2DVol (pFFT , n3);

}


//*******************************************************************************************************

//Frequency vectors
void Calc2DFreq(stFFT *pFFT , float d2 , float d1)
{
	int i;
	int n1 = pFFT->n1_nfft;
	int n2 = pFFT->n2_nfft;
    
	//dimension1
	for(i = 0 ; i < n1/2 ; i++)
		pFFT->k1[i] = 2.f*PI*i/d1/n1;

	//dimension2
	for(i = 0 ; i < n2/2 ; i++)
		pFFT->k2[i] = 2.f*PI*i/d2/n2;
	for(i = n2/2 ; i < n2 ; i++)
		pFFT->k2[i] = 2.f*PI*(i-n2)/d2/n2;

}


//*******************************************************************************************************

//Laplacian operator 2D
float Laplac2D(stFFT *pFFT , int i2 , int i1)
{
    
	return pFFT->k2[i2]*pFFT->k2[i2] + pFFT->k1[i1]*pFFT->k1[i1];

}


//*******************************************************************************************************

//Band Pass Filter
void BandPassFilter2D(stFFT *pFFT , int n3 , float band2 , float band1)
{
    int   i, n, m;
    float pass;
    
		for(m = 0 ; m < pFFT->n2_nfft ; m++)
	    	for(n = 0 ; n < pFFT->n1_nfft/2 ; n++) {
				pass = ((pFFT->k1[n] > band1) && (fabsf(pFFT->k2[m]) > band2)) ? 0.f : 1.f;
			   	for(i = 0 ; i < n3 ; i++)
			   		pFFT->M[i][m][n] *= pass;
			}

}


//*******************************************************************************************************

