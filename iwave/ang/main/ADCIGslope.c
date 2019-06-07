//#######################################################################################################
// FILE       :      ADCIGslope.c
// WRITER     :      Raanan Dafny
// EXERCISE   :      Generate scattering-angle CIGs in Fourier domain.
// DESCRIPTION:      Input Extended Image I(z,x,h).
//					 1D Fourier Transform to I(kz,x,h).
//					 Phase shift.
//					 1D inverse Fourier Transform to I(z,x,ang).
// HEADER FILES:     memory_alloc.h
//					 FFT.h
//					 ADCIG.h
// INPUT FILES :     multiRef
// OUTPUT FILE :     ADCIG
//#######################################################################################################


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <complex.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#include "memory_alloc.h"
#include "FFT.h"
#include "ADCIG.h"


int main (int argc, char *argv[]) {

	int     i, j, k, m;
    int     nthread = 1;
    float   *trace;
    stADCIG *pADCIG;
    FILE    *fin, *fout;

	int taperFact = 0;

	//Allocate ADCIG struct
	pADCIG = (stADCIG*) calloc (1, sizeof(stADCIG));
	
	//Extended Image Volume parameters
	pADCIG->nh          = atoi (argv[1]);
	pADCIG->nx          = atoi (argv[2]);
    pADCIG->nz          = atoi (argv[3]);
	pADCIG->np          = atoi (argv[4]);
    pADCIG->dh          = atof (argv[5]);
    pADCIG->dx          = atof (argv[6]);
    pADCIG->dz          = atof (argv[7]);
	pADCIG->dp          = atof (argv[8]);
    char  *fileNameIN  = argv[9];
    char  *fileNameOUT = argv[10];

	//I/O files
    if((fin = fopen(fileNameIN,"rb")) == NULL) {
      	printf("Error opening input file 1!!\n");
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
	trace 		  = alloc_array (pADCIG->nz);
	pADCIG->HOCIG = alloc_cube  (pADCIG->nx, pADCIG->nh, pADCIG->nz);
	pADCIG->ADCIG = alloc_cube  (pADCIG->nx, pADCIG->np, pADCIG->nz);
	
	printf("Memory allocated\n");


//########################################## HOCIG Input ##########################################
	
	//Input the Extended Image volume
   	for (j = 0; j < pADCIG->nh; j++) {
	   	for(k = 0; k < pADCIG->nx; k++) {
	  		if (fread(&trace[0], sizeof(float), pADCIG->nz, fin) != pADCIG->nz) {
	      		printf("Error reading x %d in h %d !!\n",j,k);
	      		exit(1);
	    	}
		    for(i = 0; i < pADCIG->nz; i++)
		    	pADCIG->HOCIG[k][j][i] = trace[i];
    	}
    }

	printf("HOCIG input\n");

	
//####################################### HOCIG to ADCIG ##########################################

	//Transform HOCIG to scattering-angle ADCIG
	HO2ADscat (pADCIG, nthread, taperFact);
	
	
//########################################### Output ADCIG ########################################

	//Output the ADCIG
	for(k = 0; k < pADCIG->nx ; k++)
		for(m = 0; m < pADCIG->np; m++) {
   			if(fwrite(&pADCIG->ADCIG[k][m][0], sizeof(float), pADCIG->nz, fout) != pADCIG->nz) {
     			printf("Error writing x %d in h %d !!\n",m,k);
    			exit(1);
    	}
   	}	
	
	printf("ADCIG Output\n");


	fclose(fin);
	fclose(fout);
	
	free(trace);
	free(pADCIG);
	
}

