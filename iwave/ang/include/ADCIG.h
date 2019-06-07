//#######################################################################################################
// HEADER:		ADCIG.h
// DESCRIPTION:	Decompose angle domain CIGs
// CONTENT:		HOCIG -> ADCIG scattering-angle
//				HOCIG -> ADCIG dip-angle
//				HOCIG -> ADCIG multi-angle (MPCIG)
//				HOCIG -> ADCIG dip-angle of h
//				ADCIG scattering-angle -> HOCIG
//				ADCIG dip-angle of h -> HOCIG
//#######################################################################################################

#define PI 3.141592654


/* ADCIG struct */
typedef struct {

	int   nz;
	int   nx;
	int   nh;
	int	  np;	/* scattering slope */
	int	  nq;	/* dip slope */
	int   nxWin;
	float dz;
	float dx;
	float dh;
	float dp;
	float dq;
	float ***HOCIG;
	float ***ADCIG;
	float *taper;
	float complex ***phaseShift;
	float complex **M1;
	float complex ***M2;
	float complex ***Maux;
	float complex *a;
	float complex *b;
	float complex *trigs1;
	float complex *trigs2;

} stADCIG;


//*******************************************************************************************************

/* scattering-angle CIG decomposition */
void HO2ADscat (stADCIG *pADCIG, int nthread, int taperFact) {

	int   i, j, k, n, m;
	int   nfft;
	int   thread_id = 0;
	float *taper;
	float complex ***phaseShift;
	float complex **Maux;
	float complex ***M;
	float complex *a, *b;
	float complex *trigs1, *trigs2;

	//Find nfft
    for(i = 0 ;pow(2, i) < pADCIG->nz ;i++);
    	nfft = (int)pow(2, i);
	nfft *= 2;

	//Memory allocation
	taper      = alloc_array          (pADCIG->nh);
	phaseShift = alloc_complex_cube   (nfft/2, pADCIG->np, pADCIG->nh);
	M		   = alloc_complex_cube   (pADCIG->nx, (int)fmax(pADCIG->np, pADCIG->nh), nfft);
	Maux       = alloc_complex_matrix (nthread, pADCIG->np);
	a 		   = alloc_complex_array  (nfft);
	b          = alloc_complex_array  (nfft);
	trigs1     = alloc_complex_array  (nfft);
	trigs2     = alloc_complex_array  (nfft/2);

	printf("Internal Memory allocated\n");

	//Build the phase shift table
	for(n = 0; n < nfft/2; n++) {
		float Kz = n * PI/pADCIG->dz/(0.5*nfft);
		for(m = 0; m < pADCIG->np; m++) {
			float Slope = (m - pADCIG->np/2) * pADCIG->dp;
			for(j = 0; j < pADCIG->nh; j++)
				phaseShift[n][m][j] = cexp(I*Kz*Slope*(j-pADCIG->nh/2)*pADCIG->dh) * pADCIG->dh;
		}
	}

	//Build the taper table
	for(j = 0; j < pADCIG->nh; j++) {
		float arg = (float) (j-pADCIG->nh/2)/(pADCIG->nh/2) * PI/2;
		taper[j] = pow(cosf(arg),taperFact);
	}

	//Copy into a complex array
	for(k = 0; k < pADCIG->nx; k++)
		for(j = 0; j < pADCIG->nh; j++)
	   		for(n = 0, i = 0; n < pADCIG->nz/2; n++, i+=2)
				M[k][j][n] = pADCIG->HOCIG[k][j][i] + pADCIG->HOCIG[k][j][i+1]*I;

   	//1D Real FFT (z->Kz)	
    ffttrig(nfft,   trigs1);
    ffttrig(nfft/2, trigs2);
	for(k = 0; k < pADCIG->nx; k++) {
    	for(j = 0; j < pADCIG->nh; j++) {
			for(n = 0; n < nfft/2; n++)
				a[n] = M[k][j][n];
			rfftrc(a, b, trigs2, trigs1, nfft);
			for(n = 0; n < nfft/2; n++)
				M[k][j][n] = a[n];
		}
	}

	int count = 0;
#ifdef _OPENMP
#pragma omp parallel default(shared) private(thread_id, n, k, m, j)
#endif
{
#ifdef _OPENMP
	thread_id = omp_get_thread_num();
#pragma omp for
#endif
	//Phase shift to angles
	for(k = 0; k < pADCIG->nx; k++) {
		for(n = 0; n < nfft/2; n++) {
		
			memset (Maux[thread_id], 0.f, pADCIG->np * sizeof(float complex));
			
			for(m = 0; m < pADCIG->np; m++)	
				for(j = 0; j < pADCIG->nh; j++)
					Maux[thread_id][m] += M[k][j][n] * phaseShift[n][m][j] * taper[j];

			for(m = 0; m < pADCIG->np; m++)
				M[k][m][n] = Maux[thread_id][m];
		}
		count++;
//		printf ("thread %d did step %d\n", thread_id, count);
	}
} // End of parallel region

	//Inverse 1D Real FFT (Kz->z)
	for(k = 0; k < pADCIG->nx; k++) {
    	for(m = 0; m < pADCIG->np; m++) {
			for(n = 0; n < nfft/2; n++)
				a[n] = M[k][m][n];
			rfftcr(a, b, trigs2, trigs1, nfft);
			for(n = 0, i = 0; n < nfft/2; n++, i+=2) {
				M[k][m][i]   = creal(a[n]);
				M[k][m][i+1] = cimag(a[n]);
			}
		}
	}

	//Copy to ADCIG
	for(k = 0; k < pADCIG->nx ; k++)
		for(m = 0; m < pADCIG->np; m++)
			for(i = 0; i < pADCIG->nz; i++)
				pADCIG->ADCIG[k][m][i] = creal(M[k][m][i]);
				
	printf("ADCIG calculated\n");
		

	free(taper);
	free(phaseShift);
	free(M);
	free(Maux);
	free(a);
	free(b);
	free(trigs1);
	free(trigs2);

}


//*******************************************************************************************************

/* dip-angle CIG decomposition */
void HO2ADdip (stADCIG *pADCIG, int nthread, int taperFact, int hfst, int hlst) {

	int   i, j, k, n, m;
	int   nfft;
	//	int   thread_id = 0;
	float *taper;
	float complex ***phaseShift;
	float complex ***Maux;
	float complex ***M;
	float complex *a, *b;
	float complex *trigs1, *trigs2;

	//Find nfft
    for(i = 0 ;pow(2, i) < pADCIG->nz ;i++);
    	nfft = (int)pow(2, i);
	nfft *= 2;

	//Memory allocation
	taper      = alloc_array          (pADCIG->nxWin);
	phaseShift = alloc_complex_cube   (nfft/2, pADCIG->nq, pADCIG->nxWin);
	M		   = alloc_complex_cube   (pADCIG->nh, pADCIG->nx, nfft);
	Maux       = alloc_complex_cube   (pADCIG->nx, nfft, pADCIG->nq);
	a 		   = alloc_complex_array  (nfft);
	b          = alloc_complex_array  (nfft);
	trigs1     = alloc_complex_array  (nfft);
	trigs2     = alloc_complex_array  (nfft/2);

	printf("Internal Memory allocated\n");

	//Build the phase shift table
	for(n = 0; n < nfft/2; n++) {
		float Kz = n * PI/pADCIG->dz/(0.5*nfft);
		for(m = 0; m < pADCIG->nq; m++) {
			float Slope = (m - pADCIG->nq/2) * pADCIG->dq;
			for(k = 0; k < pADCIG->nxWin; k++)
				phaseShift[n][m][k] = cexp(I*Kz*Slope*(k-pADCIG->nxWin/2)*pADCIG->dx) * pADCIG->dx;
		}
	}

	//Build the taper table
	for(k = 0; k < pADCIG->nxWin; k++) {
		float arg = (float) (k-pADCIG->nxWin/2)/(pADCIG->nxWin/2) * PI/2;
		taper[k] = pow(cosf(arg),taperFact);
	}

	//Copy into a complex array
	for(j = 0; j < pADCIG->nh; j++)
		for(k = 0; k < pADCIG->nx; k++)
	   		for(n = 0, i = 0; n < pADCIG->nz/2; n++, i+=2)
				M[j][k][n] = pADCIG->HOCIG[k][j][i] + pADCIG->HOCIG[k][j][i+1]*I;

   	//1D Real FFT (z->Kz)	
    ffttrig(nfft,   trigs1);
    ffttrig(nfft/2, trigs2);
   	for(j = 0; j < pADCIG->nh; j++) {
   		for(k = 0; k < pADCIG->nx; k++) {
			for(n = 0; n < nfft/2; n++)
				a[n] = M[j][k][n];
			rfftrc(a, b, trigs2, trigs1, nfft);
			for(n = 0; n < nfft/2; n++)
				M[j][k][n] = a[n];
		}
	}

	//Loop over subsurface offsets
	for(j = hfst; j <= hlst; j++) {

#ifdef _OPENMP
#pragma omp parallel default(shared) private(thread_id, i, k, m, n)
#endif
{
#ifdef _OPENMP
		thread_id = omp_get_thread_num();
#pragma omp for
#endif
		//Phase shift to angles
		for(k = 0; k < pADCIG->nx; k++) {
			for(i = 0; i < nfft/2; i++) {
		
				memset (Maux[k][i], 0.f, pADCIG->nq * sizeof(float complex));
			
				for(m = 0; m < pADCIG->nq; m++) {
					for(n = 0; n < pADCIG->nxWin; n++) {
						int xi = k + (n-pADCIG->nxWin/2);
						if (xi >=0 && xi < pADCIG->nx)
							Maux[k][i][m] += M[j][xi][i] * phaseShift[i][m][n] * taper[n];
					}
				//Scaling Factor!!
				Maux[k][i][m] /= pADCIG->nxWin;
				}
			}
		}
} // End of parallel region

		//Inverse 1D Real FFT (Kz->z)
		for(k = 0 ; k < pADCIG->nx ; k++) {
    		for(m = 0 ; m < pADCIG->nq ; m++) {
				for(i = 0 ; i < nfft/2 ; i++)
					a[i] = Maux[k][i][m];
				rfftcr(a , b , trigs2 , trigs1 , nfft);
				for(i = 0 , n = 0 ; i < nfft/2 ; i++, n+=2) {
					Maux[k][n][m]   = creal(a[i]);
					Maux[k][n+1][m] = cimag(a[i]);
				}
			}
		}

		//Copy to ADCIG (sum over all h)
		for(k = 0; k < pADCIG->nx ; k++)
			for(m = 0; m < pADCIG->nq; m++)
				for(i = 0; i < pADCIG->nz; i++)
					pADCIG->ADCIG[k][m][i] += creal(Maux[k][i][m]);
		
		printf ("End of step %d\n", j);

	} //End of subsurface offset loop

	//Scaling factor!!
	for(k = 0; k < pADCIG->nx ; k++)
		for(m = 0; m < pADCIG->nq; m++)
			for(i = 0; i < pADCIG->nz; i++)
				pADCIG->ADCIG[k][m][i] /= (hlst-hfst+1);//pADCIG->nh;

	printf("ADCIG calculated\n");


	free(taper);
	free(phaseShift);
	free(M);
	free(Maux);
	free(a);
	free(b);
	free(trigs1);
	free(trigs2);

}


//*******************************************************************************************************

/* Multi-Parameter CIG decomposition */
void HO2ADMP (stADCIG *pADCIG, int nthread, int fsth, int lsth, int xCIG) {

	int   i, j, k, n, m;
	int   nfft;
	int   thread_id = 0;
	float *taper;
	float complex ***phaseShift;
	float complex **Mauxq, **Mauxp;
	float complex ***M;
	float complex *a, *b;
	float complex *trigs1, *trigs2;

	//Find nfft
    for(i = 0 ;pow(2, i) < pADCIG->nz ;i++);
    	nfft = (int)pow(2, i);
	nfft *= 2;

	//Memory allocation
	taper      = alloc_array          ((int)fmax(pADCIG->nxWin, pADCIG->nh));
	phaseShift = alloc_complex_cube   (nfft/2, (int)fmax(pADCIG->nq, pADCIG->np), (int)fmax(pADCIG->nxWin, pADCIG->nh));
	M		   = alloc_complex_cube   ((int)fmax(pADCIG->nh, pADCIG->np), (int) fmax(pADCIG->nx, pADCIG->nq), nfft);
	Mauxq      = alloc_complex_matrix (nfft, pADCIG->nq);
	Mauxp	   = alloc_complex_matrix (nfft, pADCIG->np);
	a 		   = alloc_complex_array  (nfft);
	b          = alloc_complex_array  (nfft);
	trigs1     = alloc_complex_array  (nfft);
	trigs2     = alloc_complex_array  (nfft/2);

	printf("Internal Memory allocated\n");

	//Build the phase shift table
	for(n = 0; n < nfft/2; n++) {
		float Kz = n * PI/pADCIG->dz/(0.5*nfft);
		for(m = 0; m < pADCIG->nq; m++) {
			float Slope = (m - pADCIG->nq/2) * pADCIG->dq;
			for(k = 0; k < pADCIG->nxWin; k++)
				phaseShift[n][m][k] = cexp(I*Kz*Slope*(k-pADCIG->nxWin/2)*pADCIG->dx) * pADCIG->dx;
		}
	}

	//Build the taper table
	for(k = 0; k < pADCIG->nxWin; k++) {
		float arg = (float) (k-pADCIG->nxWin/2)/(pADCIG->nxWin/2) * PI/2;
		taper[k] = pow(cosf(arg),1);
	}

	//Copy into a complex array
	for(j = 0; j < pADCIG->nh; j++)
		for(k = 0; k < pADCIG->nx; k++)
	   		for(n = 0, i = 0; n < pADCIG->nz/2; n++, i+=2)
				M[j][k][n] = pADCIG->HOCIG[k][j][i] + pADCIG->HOCIG[k][j][i+1]*I;

   	//1D Real FFT (z->Kz)	
    ffttrig(nfft,   trigs1);
    ffttrig(nfft/2, trigs2);
   	for(j = 0; j < pADCIG->nh; j++) {
   		for(k = 0; k < pADCIG->nx; k++) {
			for(n = 0; n < nfft/2; n++)
				a[n] = M[j][k][n];
			rfftrc(a, b, trigs2, trigs1, nfft);
			for(n = 0; n < nfft/2; n++)
				M[j][k][n] = a[n];
		}
	}

	//Loop over all subsurface offsets
	for(j = fsth; j <= lsth; j++) {

#ifdef _OPENM
#pragma omp parallel default(shared) private(thread_id, i, m, n)
#endif
{
#ifdef _OPENM
		thread_id = omp_get_thread_num();
#pragma omp for
#endif
		//Phase shift to dip angle
		for(i = 0; i < nfft/2; i++) {
		
			memset (Mauxq[i], 0.f, pADCIG->nq * sizeof(float complex));
			
			for(m = 0; m < pADCIG->nq; m++) {
				for(n = 0; n < pADCIG->nxWin; n++) {
					int xi = xCIG + (n-pADCIG->nxWin/2);
					if (xi >=0 && xi < pADCIG->nx)
						Mauxq[i][m] += M[j][xi][i] * phaseShift[i][m][n] * taper[n];
				}
			}
			//Scaling Factor!!
			Mauxq[i][m] /= pADCIG->nxWin;
		}
} // End of parallel region

		//Copy back M
		for(m = 0; m < pADCIG->nq; m++)
			for(i = 0; i < nfft/2; i++)
				M[j][m][i] = Mauxq[i][m];
			
		printf ("End of step %d\n", j);
		
	} //End of subsurface offset loop

	printf("dips of h calculated\n");


	//Build the phase shift table
	for(n = 0; n < nfft/2; n++) {
		float Kz = n * PI/pADCIG->dz/(0.5*nfft);
		for(m = 0; m < pADCIG->np; m++) {
			float Slope = (m - pADCIG->np/2) * pADCIG->dp;
			for(j = 0; j < pADCIG->nh; j++)
				phaseShift[n][m][j] = cexp(I*Kz*Slope*(j-pADCIG->nh/2)*pADCIG->dh) * pADCIG->dh;
		}
	}
	
	//Build the taper table
	for(j = 0; j < pADCIG->nh; j++) {
		float arg = (float) (j-pADCIG->nh/2)/(pADCIG->nh/2) * PI/2;
		taper[j] = pow(cosf(arg),2);
	}

	int count = 0;
#ifdef _OPENMP
#pragma omp parallel default(shared) private(thread_id, n, k, m, j)
#endif
{
#ifdef _OPENMP
	thread_id = omp_get_thread_num();
#pragma omp for
#endif
	//Phase shift to scattering angle
	for(k = 0; k < pADCIG->nq; k++) {
		for(n = 0; n < nfft/2; n++) {

			memset (Mauxp[thread_id], 0.f, pADCIG->np * sizeof(float complex));

			for(m = 0; m < pADCIG->np; m++)	
				for(j = 0; j < pADCIG->nh; j++)
					Mauxp[thread_id][m] += M[j][k][n] * phaseShift[n][m][j] * taper[j];

			for(m = 0; m < pADCIG->np; m++)
				M[m][k][n] = Mauxp[thread_id][m];

		}

		count++;
		printf ("thread %d did step %d\n", thread_id, count);
	}
} // End of parallel region


	//Inverse 1D Real FFT (Kz->z)
	for(k = 0 ; k < pADCIG->nq ; k++) {
    	for(m = 0 ; m < pADCIG->np ; m++) {
			for(i = 0 ; i < nfft/2 ; i++)
				a[i] = M[m][k][i];
			rfftcr(a , b , trigs2 , trigs1 , nfft);
			for(i = 0 , n = 0 ; i < nfft/2 ; i++, n+=2) {
				M[m][k][n]   = creal(a[i]);
				M[m][k][n+1] = cimag(a[i]);
			}
		}
	}

	//Copy to ADCIG
	for(k = 0; k < pADCIG->nq ; k++)
		for(m = 0; m < pADCIG->np; m++)
			for(i = 0; i < pADCIG->nz; i++)
				pADCIG->ADCIG[k][m][i] = creal(M[m][k][i]);
   		
	printf("ADCIG calculated\n");
	
	
	free(taper);
	free(phaseShift);
	free(M);
//	free(Mauxq);
	free(Mauxp);
	free(a);
	free(b);
	free(trigs1);
	free(trigs2);
		
}


//*******************************************************************************************************

/* dip-angle CIG of h decomposition */
void HO2ADdipOfh (stADCIG *pADCIG, int nfft, int nthread, int taperFact, int ih) {

  int   i, k, n, m;
	//	int   thread_id = 0;
	float *taper;
	float complex ***phaseShift;
	float complex ***Maux;
	float complex **M;
	float complex *a, *b;
	float complex *trigs1, *trigs2;

	//Assign buffers
	taper      = pADCIG->taper;
	phaseShift = pADCIG->phaseShift;
	M		   = pADCIG->M1;
	Maux       = pADCIG->Maux;
	a 		   = pADCIG->a;
	b          = pADCIG->b;
	trigs1     = pADCIG->trigs1;
	trigs2     = pADCIG->trigs2;

	//Build the phase shift table
	for(n = 0; n < nfft/2; n++) {
		float Kz = n * PI/pADCIG->dz/(0.5*nfft);
		for(m = 0; m < pADCIG->nq; m++) {
			float Slope = (m - pADCIG->nq/2) * pADCIG->dq;
			for(k = 0; k < pADCIG->nxWin; k++)
				phaseShift[n][m][k] = cexp(I*Kz*Slope*(k-pADCIG->nxWin/2)*pADCIG->dx) * pADCIG->dx;
		}
	}

	//Build the taper table
	for(k = 0; k < pADCIG->nxWin; k++) {
		float arg = (float) (k-pADCIG->nxWin/2)/(pADCIG->nxWin/2) * PI/2;
		taper[k] = pow(cosf(arg),taperFact);
	}

	//Copy into a complex array
	for(k = 0; k < pADCIG->nx; k++) {
		memset(M[k], 0.f, nfft * sizeof(complex float));
		for(n = 0, i = 0; n < pADCIG->nz/2; n++, i+=2)
			M[k][n] = pADCIG->HOCIG[k][ih][i] + pADCIG->HOCIG[k][ih][i+1]*I;
	}

   	//1D Real FFT (z->Kz)	
    ffttrig(nfft,   trigs1);
    ffttrig(nfft/2, trigs2);
   	for(k = 0; k < pADCIG->nx; k++) {
		for(n = 0; n < nfft/2; n++)
			a[n] = M[k][n];
		rfftrc(a, b, trigs2, trigs1, nfft);
		for(n = 0; n < nfft/2; n++)
			M[k][n] = a[n];
	}

	//for current subsurface offset ih
#ifdef _OPENMP
#pragma omp parallel default(shared) private(thread_id, i, k, m, n)
#endif
{
#ifdef _OPENMP
	thread_id = omp_get_thread_num();
#pragma omp for
#endif
	//Phase shift to angles
	for(k = 0; k < pADCIG->nx; k++) {
		for(i = 0; i < nfft/2; i++) {
	
			memset (Maux[k][i], 0.f, pADCIG->nq * sizeof(float complex));
		
			for(m = 0; m < pADCIG->nq; m++) {
				for(n = 0; n < pADCIG->nxWin; n++) {
					int xi = k + (n-pADCIG->nxWin/2);
					if (xi >=0 && xi < pADCIG->nx)
						Maux[k][i][m] += M[xi][i] * phaseShift[i][m][n] * taper[n];
				}
				//Scaling Factor!!
				Maux[k][i][m] /= pADCIG->nxWin;
			}
		}
	}
} // End of parallel region

	//Inverse 1D Real FFT (Kz->z)
	for(k = 0 ; k < pADCIG->nx ; k++) {
   		for(m = 0 ; m < pADCIG->nq ; m++) {
			for(i = 0 ; i < nfft/2 ; i++)
				a[i] = Maux[k][i][m];
			rfftcr(a , b , trigs2 , trigs1 , nfft);
			for(i = 0 , n = 0 ; i < nfft/2 ; i++, n+=2) {
				Maux[k][n][m]   = creal(a[i]);
				Maux[k][n+1][m] = cimag(a[i]);
			}
		}
	}

	//Copy to ADCIG (current h)
	for(k = 0; k < pADCIG->nx ; k++)
		for(m = 0; m < pADCIG->nq; m++)
			for(i = 0; i < pADCIG->nz; i++)
				pADCIG->ADCIG[k][m][i] = creal(Maux[k][i][m]);
	
	printf("ADCIG calculated for ih %d\n", ih);
	
}


//*******************************************************************************************************

/* Inverse scattering-angle CIG decomposition */
void ADscat2HO (stADCIG *pADCIG, int nthread, int taperFact) {

	int   i, j, k, n, m;
	int   nfft;
	int   thread_id = 0;
	float *taper;
	float complex ***phaseShift;
	float complex **Maux;
	float complex ***M;
	float complex *a, *b;
	float complex *trigs1, *trigs2;

	//Find nfft
    for(i = 0 ;pow(2, i) < pADCIG->nz ;i++);
    	nfft = (int)pow(2, i);
	nfft *= 2;

	//Memory allocation
	taper      = alloc_array          (pADCIG->np);
	phaseShift = alloc_complex_cube   (nfft/2, pADCIG->np, pADCIG->nh);
	M		   = alloc_complex_cube   (pADCIG->nx, (int)fmax(pADCIG->np,pADCIG->nh), nfft);
	Maux       = alloc_complex_matrix (nthread, pADCIG->nh);
	a 		   = alloc_complex_array  (nfft);
	b          = alloc_complex_array  (nfft);
	trigs1     = alloc_complex_array  (nfft);
	trigs2     = alloc_complex_array  (nfft/2);

	printf("Internal Memory allocated\n");

	//Build the phase shift table
	for(n = 0; n < nfft/2; n++) {
		float Kz = n * PI/pADCIG->dz/(0.5*nfft);
		for(m = 0; m < pADCIG->np; m++) {
			float Slope = (m - pADCIG->np/2) * pADCIG->dp;
			for(j = 0; j < pADCIG->nh; j++)
				phaseShift[n][m][j] = cexp(-I*Kz*Slope*(j-pADCIG->nh/2)*pADCIG->dh) * Kz * pADCIG->dp / (2*PI);
		}
	}

	//Build the taper table
	for(m = 0; m < pADCIG->np; m++) {
		float arg = (float) (m-pADCIG->np/2)/(pADCIG->np/2) * PI/2;
		taper[m] = pow(cosf(arg),taperFact);
	}
	
	//Copy into a complex array
	for(k = 0; k < pADCIG->nx; k++)
		for(m = 0; m < pADCIG->np; m++)
	   		for(n = 0, i = 0; n < pADCIG->nz/2; n++, i+=2)
				M[k][m][n] = pADCIG->ADCIG[k][m][i] + pADCIG->ADCIG[k][m][i+1]*I;

   	//1D Real FFT (z->Kz)
    ffttrig(nfft,   trigs1);
    ffttrig(nfft/2, trigs2);
	for(k = 0; k < pADCIG->nx; k++) {
    	for(m = 0; m < pADCIG->np; m++) {
			for(n = 0; n < nfft/2; n++)
				a[n] = M[k][m][n];
			rfftrc(a, b, trigs2, trigs1, nfft);
			for(n = 0; n < nfft/2; n++)
				M[k][m][n] = a[n];
		}
	}

	int count = 0;
#ifdef _OPENMP
#pragma omp parallel default(shared) private(thread_id, i, k, m, j)
#endif
{
#ifdef _OPENMP
	thread_id = omp_get_thread_num();
#pragma omp for
#endif
	//Phase shift to subsurface offset
	for(k = 0; k < pADCIG->nx; k++) {
		for(i = 0; i < nfft/2; i++) {
		
			memset (Maux[thread_id], 0.f, pADCIG->nh * sizeof(float complex));
			
			for(j = 0; j < pADCIG->nh; j++)
				for(m = 0; m < pADCIG->np; m++)	
					Maux[thread_id][j] += M[k][m][i] * phaseShift[i][m][j] * taper[m];

			for(j = 0; j < pADCIG->nh; j++)
				M[k][j][i] = Maux[thread_id][j];
		}
		count++;
		printf ("thread %d did step %d\n", thread_id, count);
	}
} // End of parallel region

	//Inverse 1D Real FFT (Kz->z)
	for(k = 0; k < pADCIG->nx; k++) {
    	for(j = 0; j < pADCIG->nh; j++) {
			for(n = 0; n < nfft/2; n++)
				a[n] = M[k][j][n];
			rfftcr(a, b, trigs2, trigs1, nfft);
			for(n = 0, i = 0; n < nfft/2; n++, i+=2) {
				M[k][j][i]   = creal(a[n]);
				M[k][j][i+1] = cimag(a[n]);
			}
		}
	}

	//Copy to HOCIG
	for(k = 0; k < pADCIG->nx ; k++)
		for(j = 0; j < pADCIG->nh; j++)
			for(i = 0; i < pADCIG->nz; i++)
				pADCIG->HOCIG[k][j][i] = creal(M[k][j][i]);
				
	printf("HOCIG calculated\n");
		

	free(taper);
	free(phaseShift);
	free(M);
	free(Maux);
	free(a);
	free(b);
	free(trigs1);
	free(trigs2);

}


//*******************************************************************************************************

/* Inverse dip-angle CIG of h decomposition */
void ADdipOfh2HO (stADCIG *pADCIG, int nfft, int nthread, int taperFact, int ih) {

  int   i, k, n, m;
  //	int   thread_id = 0;
	float *taper;
	float complex ***phaseShift;
	float complex ***Maux;
	float complex **M;
	float complex *a, *b;
	float complex *trigs1, *trigs2;

	//Assign buffers
	taper      = pADCIG->taper;
	phaseShift = pADCIG->phaseShift;
	M		   = pADCIG->M1;
	Maux       = pADCIG->Maux;
	a 		   = pADCIG->a;
	b          = pADCIG->b;
	trigs1     = pADCIG->trigs1;
	trigs2     = pADCIG->trigs2;

	//Build the phase shift table
	for(n = 0; n < nfft/2; n++) {
		float Kz = n * PI/pADCIG->dz/(0.5*nfft);
		for(m = 0; m < pADCIG->nq; m++) {
			float Slope = (m - pADCIG->nq/2) * pADCIG->dq;
			for(k = 0; k < pADCIG->nxWin; k++)
				phaseShift[n][m][k] = cexp(-I*Kz*Slope*(k-pADCIG->nxWin/2)*pADCIG->dx) * Kz * pADCIG->dq / (2*PI);
		}
	}

	//Build the taper table
	for(m = 0; m < pADCIG->nq; m++) {
		float arg = (float) (m-pADCIG->nq/2)/(pADCIG->nq/2) * PI/2;
		taper[m] = pow(cosf(arg),taperFact);
	}

	//Copy into a complex array
	for(k = 0 ; k < pADCIG->nx ; k++)
		for(m = 0 ; m < pADCIG->nq ; m++) {
			memset(Maux[k][m] , 0. , nfft * sizeof(float complex));
	   		for(n = 0 , i = 0 ; n < pADCIG->nz/2 ; n++, i+=2)
				Maux[k][m][n] = pADCIG->ADCIG[k][m][i] + pADCIG->ADCIG[k][m][i+1]*I;
		}

   	//1D Real FFT (z->Kz)	
    ffttrig(nfft,   trigs1);
    ffttrig(nfft/2, trigs2);
	for(k = 0 ; k < pADCIG->nx ; k++) {
    	for(m = 0 ; m < pADCIG->nq ; m++) {
			for(i = 0 ; i < nfft/2 ; i++)
				a[i] = Maux[k][m][i];
			rfftrc(a, b ,trigs2, trigs1, nfft);
			for(i = 0 ; i < nfft/2 ; i++)
				Maux[k][m][i] = a[i];
		}
	}

#ifdef _OPENMP
#pragma omp parallel default(shared) private(thread_id, i, k, m, n)
#endif
{
#ifdef _OPENMP
		thread_id = omp_get_thread_num();
#pragma omp for
#endif
		//Phase shift to subsurface offset
		for(i = 0; i < nfft/2; i++) {
		
			memset(M[i], 0., pADCIG->nx * sizeof(float complex));
			for(k = 0; k < pADCIG->nx; k++) {

				for(n = 0; n < pADCIG->nxWin; n++) {
					int xi = k + (n-pADCIG->nxWin/2);
					if (xi >=0 && xi < pADCIG->nx) {
						for(m = 0; m < pADCIG->nq; m++)
							M[i][xi] += Maux[k][m][i] * phaseShift[i][m][n] * taper[m];
					}
				}

			}
		}
} // End of parallel region

	//Inverse 1D Real FFT (Kz->z)
	for(k = 0 ; k < pADCIG->nx ; k++) {
		for(i = 0 ; i < nfft/2 ; i++)
			a[i] = M[i][k];
		rfftcr(a , b , trigs2 , trigs1 , nfft);
		for(i = 0 , n = 0 ; i < nfft/2 ; i++, n+=2) {
			M[n][k]   = creal(a[i]);
			M[n+1][k] = cimag(a[i]);
		}
	}

	//Copy to HOCIG (current h)
	for(k = 0; k < pADCIG->nx ; k++)
		for(i = 0; i < pADCIG->nz; i++)
			pADCIG->HOCIG[k][ih][i] = creal(M[i][k]);
	
	printf("HOCIG calculated for ih %d\n", ih);
	
}


//*******************************************************************************************************
			
