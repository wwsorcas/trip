//#######################################################################################################
// FILE       :      Wd_2017.c
// WRITER     :      Raanan Dafny
// EXERCISE   :      Approximate inverse - VERSION_2 (Hou and Symes, 2017)
// DESCRIPTION:      Data side weight: Wd = It^3.
//					 Intergration via FFT.
// HEADER FILES:     memory_alloc.h
//					 numeric.h
// INPUT FILES :     myborn
// OUTPUT FILE :     mybornWd
//#######################################################################################################

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <complex.h>
#include "memory_alloc.h"
#include "numeric.h"


int main (int argc, char *argv[]) {

	int   j, k;
 	float *IN, *OUT;
	FILE  *fin, *fout;
	
	int    ns          = atoi (argv[1]);
	int    nr          = atoi (argv[2]);
    int    nt          = atoi (argv[3]);
    float  dt          = atof (argv[4]);
    float  f1          = atof (argv[5]);
    float  f2          = atof (argv[6]);
    float  f3          = atof (argv[7]);
    float  f4          = atof (argv[8]);
    char  *fileNameIN  = argv[9];
    char  *fileNameOUT = argv[10];
 	
 	//Band Pass filter
	float band[4] = {f1, f2, f3, f4};
 	
	//Differentiate (>0) or integrate (<0)
 	int   power = -3;

   	//I/O files
	if((fin = fopen(fileNameIN,"rb")) == NULL) {
      	printf("Error opening input file 1!!\n");
      	exit(1);
    }    
	fout = fopen(fileNameOUT,"wb");

//####################################### Memory Allocation #######################################
	
	//Memory allocation
	IN  = alloc_array (ns*nr*nt);
	OUT = alloc_array (ns*nr*nt);

	printf("Memory allocated\n");

//########################################## Input Data ##########################################
	
	//Input data traces
	for (k = 0; k < ns; k++)
		for (j = 0; j < nr; j++) {
			int i0 = k*nt*nr + j*nt;
			if (fread(&IN[i0], sizeof(float), nt, fin) != nt) {
				printf("Error reading trace !!\n");
	    		exit(1);
			}
		}

	printf("Data input\n");
		
//################################# Differentiate or Integrate ####################################	

	//Differentiation or Integration via Fourier transform
	//Band-Pass filter is applied
	DerivIntegral(IN, OUT, ns, nr, nt, dt, power, band);
	
	printf("Derivation/Integration completed\n");
	
//########################################## Output data #########################################
	
	for (k = 0; k < ns; k++)
		for (j = 0; j < nr; j++) {
			int i0 = k*nt*nr + j*nt;
   			if(fwrite(&OUT[i0], sizeof(float), nt, fout) != nt) {
				printf("Error writing trace!!\n");
    			exit(1);
    		}
		}
		
	printf("Data output\n");

	fclose (fin);
	fclose (fout);

	free(IN);
	free(OUT);
	

}
