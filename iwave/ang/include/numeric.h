//#######################################################################################################
// HEADER:		numeric.h
// DESCRIPTION:	numeric utilities
// CONTENT:		first derivetive (FD central)
//				median filter
//				differentiation / integration operator (Fourier domain)
//
//#######################################################################################################

#include "FFT.h"

#define PI 3.141592654


/* Finite Difference scheme - 1st derivetive */
void FD1st (float* f, float* g, int nx, float dx_inv) {

	int i;
	
	//Forward FD
	g[0] = (f[1] - f[0]) * dx_inv;
	
	//Central FD
	for (i = 0; i < nx-1; i++)
		g[i] = (f[i+1]-f[i-1]) * dx_inv * 0.5;
	
	//Backward FD
	g[nx-1] = (f[nx-1] - f[nx-2]) * dx_inv;

}


//*******************************************************************************************************

/* Median Filter - 1D */

int cmpfunc (const void* a, const void* b) {
	return (*(int*)a - *(int*)b);
}

void Median1D (float *f, float *g, int wind, int nx) {

	int i, j;
	float* buff = &g[nx];
	
	//Median Filter
	for (i = 0; i < wind; i++)
		g[i] = f[i];
	for (i = wind; i < nx-wind; i++) {
		for (j = -wind; j <= wind; j++)
			buff[j+wind] = f[i+j];
		qsort(buff, wind*2+1, sizeof(float), cmpfunc);
		g[i] = buff[wind];
	}
	for (i = nx-wind; i < nx; i++)
		g[i] = f[i];

}


//*******************************************************************************************************

/* Time Differentiation / Integration of siesmic traces */
void DerivIntegral(float *IN , float *OUT , int n3, int n2, int n1, float dt, int power, float *band) {

	int   i, j, k, nfft;
	int   isign, if1, if2, if3, if4;
	float df, arg, taper;
	float complex fact1, fact2;
	float complex *a, *b, *trigs;
	
	//Find nfft
    for(i = 0 ;pow(2, i) < n1 ;i++);
    	nfft = (int)pow(2, i);
	nfft *= 2;
	
	//Memory allocation
	a 	  = alloc_complex_array  (nfft);
	b     = alloc_complex_array  (nfft);
	trigs = alloc_complex_array  (nfft);

	printf("Internal Memory allocated\n");

	ffttrig(nfft,   trigs);
    
    //For band pass filter
	df  = 2.f*PI / dt / nfft;
	if1 = floor(2.f*PI * band[0] / df);
	if2 = floor(2.f*PI * band[1] / df);
	if3 = ceil (2.f*PI * band[2] / df);
	if4 = ceil (2.f*PI * band[3] / df);
	
			    
	// Loop over input traces
	for(k = 0; k < n3; k++)
		for(j = 0; j < n2; j++) {

			int i0 = k*n1*n2 + j*n1;
			memset(a, 0., nfft*sizeof(float complex));
			
//########################################### Forward FFT #########################################

			//Copy into a complex array
			for(i = 0; i < n1; i++)
				a[i] = IN[i0+i];
		   	//Complex DFT
			isign = -1;
			fft(a, b, trigs, nfft, isign);
	
//			printf("Fwd FFT\n");
		
//################################# Differentiate or Integrate ####################################	
			
			//Loop over all frequencies in the band, scale by (I*f)^power
			for(i = 1; i < nfft/2; i++) {
				
				fact1 = cpowf(i*df*I, power);
				fact2 = cpowf(-i*df*I, power);
				if ((i < if1) || (i >= if4)) {
					a[i]      = (float complex)0.f;
					a[nfft-i] = (float complex)0.f;
				} else if ((i >= if1) && (i < if2)) {
					arg   = (float) (if2-i)/(if2-if1) * PI/2.f;
					taper = pow(cosf(arg), 2);
					a[i]      *= fact1 * taper;
					a[nfft-i] *= fact2 * taper;
				} else if ((i >= if2) && (i < if3)) {
					a[i]      *= fact1;
					a[nfft-i] *= fact2;
				} else if ((i >= if3) && (i < if4)) {
					arg = (float) (i-if3)/(if4-if3) * PI/2.f;
					taper = pow(cosf(arg), 2);
					a[i]      *= fact1 * taper;
					a[nfft-i] *= fact2 * taper;
				}
		
			}
			a[0]      = (float complex)0.f;
			a[nfft/2] = (float complex)0.f;


//########################################### Backward FFT ########################################

			//Inverse complex DFT
			isign = 1;
			fft(a, b, trigs, nfft, isign);
	
			// copy to output traces
			for(i = 0; i < n1; i++)
				OUT[i0+i] = creal(a[i]) / nfft;

//			printf("Bwd FFT\n");

		}
	
	
	free(a);
	free(b);
	free(trigs);

}


//*******************************************************************************************************

