//#######################################################################################################
// FILE       :      DeconLS.c
// WRITER     :      Raanan Dafny
// EXERCISE   :      Least-Square Deconvolution.
// DESCRIPTION:      Input data(x,gamma,z).
//					 Compute Convolution matrix.
//					 Conjugate Gradient solver for least-square inversion.
// HEADER FILES:     memory_alloc.h
//					 matrix.h
//					 AVA.h
//					 Interp.h
// INPUT FILES :     Data
// OUTPUT FILE :     DeconData
//#######################################################################################################

#define PI 3.141592654

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#include "memory_alloc.h"
#include "matrix.h"
#include "AVA.h"
#include "Interp.h"


int main (int argc, char *argv[]) {

	int     i, j, n0, m0;
	int	    Nt, Ntz;
	//    int     nthread = 1;
    float   dz, dt, tmax;
    float   *t, *wvlt, *wvlt_tz;
    stAVA   *pAVA;
    FILE    *fin1, *fin2, *fout;
    
   	//Allocate AVA struct
	pAVA = (stAVA*) calloc (1, sizeof(stAVA));
	
    //Image volume parameters
    pAVA->nz            = atoi (argv[1]);
    pAVA->ng            = atoi (argv[2]);
    pAVA->nx            = atoi (argv[3]);
	dz					= atof (argv[4]);
	Nt                  = atoi (argv[5]);
	dt  				= atof (argv[6]);
    char  *fileNameIN1  = argv[7];
    char  *fileNameIN2  = argv[8];
    char  *fileNameOUT  = argv[9];
    
	tmax = dt*(Nt-1);
	Ntz  = (int)(tmax/dz) + 1;
	
	//Convolution matrix dimensions
	int Ntraces = pAVA->nx*pAVA->ng;
	int NofA    = pAVA->nz;
	int MofA    = pAVA->nz+Ntz-1;
	int N       = NofA*Ntraces;
	int M       = MofA*Ntraces;
	
	//Conjugate Gradient parameters
	int   itrmaxCg = 120;
	float eps      = 0.1e-3;//0.1e-4;
	    
    //I/O files
    if ((fin1 = fopen(fileNameIN1,"rb")) == NULL) {
      	printf("Error opening input file 1!!\n");
      	exit(1);
    }
    if ((fin2 = fopen(fileNameIN2,"rb")) == NULL) {
      	printf("Error opening input file 2!!\n");
      	exit(1);
    }
   	fout = fopen(fileNameOUT,"wb");

    //omp setting
#ifdef _OPENMP
	omp_set_num_threads(8);
	nthread = omp_get_max_threads();
#endif

//####################################### Memory Allocation #######################################

	//Memory allocation
	t          = alloc_array  (Nt);
	wvlt       = alloc_array  (Nt);
	wvlt_tz    = alloc_array  (Ntz);
	pAVA->x    = alloc_array  (M);
	pAVA->d    = alloc_array  (N);
	pAVA->Atd  = alloc_array  (M);
	pAVA->r    = alloc_array  (M);
	pAVA->v    = alloc_array  (M);
	pAVA->Av   = alloc_array  (M);
	pAVA->A    = alloc_matrix (NofA, MofA);
	pAVA->AtA  = alloc_matrix (MofA, MofA);
	
	printf("Memory allocated\n");

//########################################## Input data ###########################################

	//Input the data
	if (fread(&wvlt[0], sizeof(float), Nt, fin1) != Nt) {
	   	 	printf("Error reading wavelet trace !!\n");
	   	 	exit(1);
	}
	for(j = 0; j < Ntraces; j++) {
		n0 = j*pAVA->nz;
	  	if (fread(&pAVA->d[n0], sizeof(float), pAVA->nz, fin2) != pAVA->nz) {
	   	 	printf("Error reading data trace %d !!\n", j);
	   	 	exit(1);
	    }
	}

	//Interpolate the source wavelet trace
	for(i = 0; i < Nt; i++)
		t[i] = i*dt;
	for(i = 0; i < Ntz; i++)
		wvlt_tz[i] = LG4PInterp (wvlt, t, i*dz, Nt);
	
	printf("Data input\n");

//################################## Compute Convolution matrix ###################################

	//Build convolution matrix
	BuildConvMat (pAVA->A, pAVA->x, wvlt_tz, NULL, dz, Ntz, pAVA->nz);
	Mat_ATA (pAVA->A , pAVA->AtA , NofA , MofA);
	printf ("Built AtA\n");

#ifdef _OPENMP
#pragma omp parallel default(shared) private(j, n0, m0) num_threads(8)
#endif
{
#ifdef _OPENMP
#pragma omp for
#endif
	//Compute Atd
	for (j = 0; j < Ntraces; j++) {
		n0 = j*NofA;
		m0 = j*MofA;
		Mat_ATd (pAVA->A , &pAVA->d[n0], &pAVA->Atd[m0], NofA , MofA);
	}
} // End of parallel region
	printf ("Built Atd\n");
		
//####################################### LS Deconvolution ########################################

	//Conjugate Gradient solver
	if (!CG_AVA (pAVA->x, pAVA->r, pAVA->v, pAVA->Av, pAVA->Atd, pAVA->AtA, M, MofA, itrmaxCg, eps))
		printf ("CG Error: residue too large after %d iterations!\n", itrmaxCg);

	printf ("CG completed\n");

//########################################## Output data ##########################################

	//Output deconvolved data
	for (j = 0; j < Ntraces; j++) {
		m0 = j*MofA;
		if (fwrite(&pAVA->x[m0+Ntz/2], sizeof(float), NofA, fout) != NofA) {
	   		printf ("Error writing trace %d !!\n", j);
    		exit(1);
    	}
    }
	
	printf ("Data Output\n");


	fclose(fin1);
	fclose(fin2);
	fclose(fout);
	
	free(t);
	free(wvlt);
	free(wvlt_tz);
	free(pAVA);

}
