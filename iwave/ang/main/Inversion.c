//#######################################################################################################
// FILE       :      Inversion.c
// WRITER     :      Raanan Dafny
// EXERCISE   :      AVA inversion of angle-domain CIGs.
// DESCRIPTION:      Input CIGs d(x, gamma, z).
//					 Compute AtA and Atd.
//					 Conjugate Gradient solver for Intercept and Gradient.
// HEADER FILES:     memory_alloc.h
//					 matrix.h
//					 AVA.h
// INPUT FILES :     ADCIG
// OUTPUT FILE :     Intercept, Gradient
//#######################################################################################################

#define PI 3.141592654

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#include "memory_alloc.h"
#include "matrix.h"
#include "AVA.h"


int main (int argc, char *argv[]) {

	int     i, j, k, n0, m0;
	//    int     nthread = 1;
    float   *trace;
    stAVA   *pAVA;
    FILE    *fin, *fout1, *fout2;
    
   	//Allocate AVA struct
	pAVA = (stAVA*) calloc (1, sizeof(stAVA));
	
    //Image volume parameters
    pAVA->ng            = atoi (argv[1]);
    pAVA->nx            = atoi (argv[2]);
    pAVA->nz            = atoi (argv[3]);
    pAVA->dg            = atof (argv[4]);
	pAVA->og            = atof (argv[5]);
    char  *fileNameIN   = argv[6];
    char  *fileNameOUT1 = argv[7];
    char  *fileNameOUT2 = argv[8];
    
	//AVA matrix dimensions
	int NofA = pAVA->nz*pAVA->ng;
	int MofA = 2*pAVA->nz;
	int N    = NofA*pAVA->nx;
	int M    = MofA*pAVA->nx;
	
	//Conjugate Gradient parameters
	int   itrmaxCg = 60;
	float eps      = 0.1e-5;
	    
    //I/O files
    if ((fin = fopen(fileNameIN,"rb")) == NULL) {
      	printf("Error opening input file !!\n");
      	exit(1);
    }    
   	fout1 = fopen(fileNameOUT1,"wb");
   	fout2 = fopen(fileNameOUT2,"wb");

    //omp setting
#ifdef _OPENMP
	omp_set_num_threads(8);
	nthread = omp_get_max_threads();
#endif

//####################################### Memory Allocation #######################################

	//Memory allocation
	trace      = alloc_array  (pAVA->nz);
	pAVA->sing = alloc_array  (pAVA->ng);
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

	//Input the angle-domain CIG
   	for (k = 0; k < pAVA->nx; k++)
		for (j = 0; j < pAVA->ng; j++) {
			n0 = (k*pAVA->nz*pAVA->ng) + (j*pAVA->nz);
			if (fread(&pAVA->d[n0], sizeof(float), pAVA->nz, fin) != pAVA->nz) {
		   		printf("Error reading angle %d in x %d !!\n", j, k);
		   		exit(1);
			}
		}

	printf("Data input\n");

//###################################### Compute AVA matrix #######################################

	//Build AVA matrix
	BuildAVAMat (pAVA);	
	printf("Built AVA matrix\n");

	//Compute AtA
	Mat_ATA   (pAVA->A , pAVA->AtA , NofA , MofA);
	printf("Built AtA\n");

#ifdef _OPENMP
#pragma omp parallel default(shared) private(k, n0, m0) num_threads(8)
#endif
{
#ifdef _OPENMP
#pragma omp for
#endif
	//Compute Atd
	for (k = 0; k < pAVA->nx; k++) {
		n0 = (k*pAVA->nz*pAVA->ng);
		m0 = (k*pAVA->nz);
		Mat_ATd (pAVA->A , &pAVA->d[n0], &pAVA->Atd[2*m0], NofA , MofA);
	}
} // End of parallel region
	printf("Built Atd\n");
		
//######################################### AVA inversion #########################################

	//Conjugate Gradient solver
	if (!CG_AVA (pAVA->x, pAVA->r, pAVA->v, pAVA->Av, pAVA->Atd, pAVA->AtA, M, MofA, itrmaxCg, eps))
		printf("CG Error: residue too large after %d iterations!\n", itrmaxCg);

	printf("CG completed\n");

//######################################### Output model ##########################################

	//Output inverted Intercept
   	for (k = 0; k < pAVA->nx; k++) {
		for (i = 0; i < pAVA->nz; i++) {
			m0 = (k*pAVA->nz) + i;
			trace[i] = pAVA->x[2*m0];
		}
		if (fwrite(&trace[0], sizeof(float), pAVA->nz, fout1) != pAVA->nz) {
			printf("Error writing Intercept at x %d !!\n", k);
			exit(1);
		}
	}
	//Output inverted Gradient
   	for (k = 0; k < pAVA->nx; k++) {
		for (i = 0; i < pAVA->nz; i++) {
			m0 = (k*pAVA->nz) + i;
			trace[i] = pAVA->x[2*m0+1];
		}
		if (fwrite(&trace[0], sizeof(float), pAVA->nz, fout2) != pAVA->nz) {
    		printf("Error writing Gradient at x %d !!\n", k);
    		exit(1);
    	}
    }

	printf("Model output\n");


	fclose(fin);
	fclose(fout1);
	fclose(fout2);
	
	free(trace);
	free(pAVA);

}
