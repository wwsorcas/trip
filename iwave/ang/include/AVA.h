//#######################################################################################################
// HEADER:		AVA.h
// DESCRIPTION:	AVA Inversion
// CONTENT:		Build AVA Matrix
//				Inversion
//#######################################################################################################

#define PI 3.141592654


/* AVA struct */
typedef struct {

	int   nz;		//depth samples
	int   nx;		//horizontal samples
	int	  ng;		//angle
	float dg;		//angle increment
	float og;		//angle origin
	float *sing;	//sin^2(angle)
	float *x;
	float *d;
	float *Atd;
	float *r;
	float *v;
	float *Av;
	float **A;		//AVA matrix
	float **AtA;

} stAVA;


//*******************************************************************************************************

/*  Build wavelet - FIRST derivative of a Gaussian */
void GaussDeriv1Wavelet (float *trace , float dt , int Nt , float fpeak) {

	int i;
	
	float t0 =  Nt/2*dt; 					  //zero crossing of the wavelet
	float s  =  1.0 / (sqrt(2) * PI * fpeak); //sigma (sqrt of the variance)
	float a  =  1.0 / (2*s*s);
	float b  = -1.0 / (sqrt(2*PI)*pow(s,3)) * 1.e6;
	
	for(i = 0; i < Nt; i++)
		trace[i] = b * (i*dt-t0) * exp(-a*pow(i*dt-t0,2));

}


//*******************************************************************************************************

/*  Build wavelet - SECOND derivative of a Gaussian (Ricker Wavelet) */
void GaussDeriv2Wavelet (float *trace , float dt , int Nt , float fpeak) {

	int i;
	
	float t0 =  Nt/2*dt; 					  //zero crossing of the wavelet
	float s  =  1.0 / (sqrt(2) * PI * fpeak); //sigma (sqrt of the variance)
	float a  =  1.0 / (2*s*s);
	float b  =  1.0 / (sqrt(2*PI)*pow(s,3));
	
	for(i = 0; i < Nt; i++)
		trace[i] = b * (1 - pow((i*dt-t0)/s,2)) * exp(-a*pow(i*dt-t0,2));

}


//*******************************************************************************************************

/* Normalization */
void Normalization (float *trace , int Nt) {

	int   i;
	float Max = 0;

	//Find Max value in the trace
	for(i = 0; i < Nt; i++)
		Max = fabsf(trace[i]) > Max ? fabsf(trace[i]) : Max;
	//Normalization
	for(i = 0; i < Nt; i++)
		trace[i] /= Max; 

}


//*******************************************************************************************************

/* Convolution function */
void conv (float*f,float*g,float*h,float dt,int N,int M) {
	int i, j, k;
	int k0 = N/2;

	memset (h, 0, M*sizeof(float));
  
	for (i = 0; i < N; i++) {
		for (j = 0; j < M; j++) {
			k = i + j;
		  	//Only size of M is needed for the convolution
		  	if ((k >= k0) && (k < k0+M))
		  		h[k-k0] += f[i] * g[j] * dt;
	  	}
    }
}


//*******************************************************************************************************

/* Build Convolution Matrix */
void BuildConvMat (float **A, float *x, float *f, float *g, float dt, int Nf, int Ng) {
	int i, j, k, k0;
	
	int N = Ng;
	//	int M = Ng+Nf-1;
	k0 = Nf/2;
	
	//Build A
	for (i = 0; i < N; i++)
		for (j = Nf-1, k = i; j >=0; j--, k++)
			A[i][k] = f[j] * dt;
	//Build x
	if (g != NULL)
		for (j = 0; j < Ng; j++)
			x[Nf-1-k0+j] = g[j];
}


//*******************************************************************************************************

/* Build AVA Matrix */
void BuildAVAMat (stAVA *pAVA) {

	int i, j, n0, m0;

	//Array of sin^2(gamma)
	for (j = 0; j < pAVA->ng; j++) {
		float ang = (pAVA->og + j*pAVA->dg) * PI/180.f;
		pAVA->sing[j] = cosf(2.f*ang); //sinf(ang) * sinf(ang);
	}

	//AVA matrix
	for (j = 0; j < pAVA->ng; j++) 
		for (i = 0; i < pAVA->nz; i++) {
			n0 = (j*pAVA->nz) + i;
			m0 =  i;
			pAVA->A[n0][2*m0]   = 1.f;
			pAVA->A[n0][2*m0+1] = 1.f * pAVA->sing[j];
		}

}


//*******************************************************************************************************

/* Build AVA Matrix  - Shuey (1985)*/
void BuildAVAMat_Shuey (stAVA *pAVA) {

	int i, j, n0, m0;

	//Array of sin^2(gamma)
	for (j = 0; j < pAVA->ng; j++) {
		float ang = (pAVA->og + j*pAVA->dg) * PI/180.f;
		pAVA->sing[j] = sinf(ang) * sinf(ang);
	}

	//AVA matrix
	for (j = 0; j < pAVA->ng; j++)
		for (i = 0; i < pAVA->nz; i++) {
			n0 = (j*pAVA->nz) + i;
			m0 =  i;
			pAVA->A[n0][2*m0]   =  0.5 / (1.f - pAVA->sing[j]);
			pAVA->A[n0][2*m0+1] = -0.5 * pAVA->sing[j] / (1.f - pAVA->sing[j]);
		}

}

//*******************************************************************************************************

