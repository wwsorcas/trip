//#######################################################################################################
// FILE       :      Wm_2017.c
// WRITER     :      Raanan Dafny
// EXERCISE   :      Approximate inverse - VERSION_2 (Hou and Symes, 2017)
// DESCRIPTION:      Model side weight: Wm = 32*DEN^2*VEL^6*Dz.
//					 Differentiation via FFT.
// HEADER FILES:     memory_alloc.h
//					 numeric.h
// INPUT FILES :     mymig
// OUTPUT FILE :     mymigWm
//#######################################################################################################

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <complex.h>
#include "memory_alloc.h"
#include "numeric.h"


int main (int argc, char *argv[]) {

	int   i, j, k;
 	float *IN, *OUT, *VEL, *DEN;
	FILE  *fin1, *fin2, *fin3, *fout;
	
	int    nh          = atoi (argv[1]);
	int    nx          = atoi (argv[2]);
    int    nz          = atoi (argv[3]);
    float  dz          = atof (argv[6]);
    float  f1          = atof (argv[7]);
    float  f2          = atof (argv[8]);
    float  f3          = atof (argv[9]);
    float  f4          = atof (argv[10]);
    char  *fileNameIN  = argv[11];
    char  *fileNameOUT = argv[12];
    char  *fileNameVEL = argv[13];
    char  *fileNameDEN = argv[14];

 	//Band Pass filter
	float band[4] = {f1, f2, f3, f4};
 	
	//Differentiate (>0) or integrate (<0)
 	int   power = 1;

   	//I/O files
	if((fin1 = fopen(fileNameIN,"rb")) == NULL) {
      	printf("Error opening input file 1!!\n");
      	exit(1);
    }
    if((fin2 = fopen(fileNameVEL,"rb")) == NULL) {
        printf("Error opening input file 2!!\n");
        exit(1);
    }
    if((fin3 = fopen(fileNameDEN,"rb")) == NULL) {
        printf("Error opening input file 3!!\n");
        exit(1);
    }


	fout = fopen(fileNameOUT,"wb");

//####################################### Memory Allocation #######################################
	
	//Memory allocation
	IN  = alloc_array (nh*nx*nz);
	OUT = alloc_array (nh*nx*nz);
        VEL = alloc_array (nx*nz);
        DEN = alloc_array (nx*nz);

	printf("Memory allocated\n");

//########################################## Input Volumes ########################################
	
	//Input image
	for (k = 0; k < nh; k++)
		for (j = 0; j < nx; j++) {
			int i0 = k*nz*nx + j*nz;
			if (fread(&IN[i0], sizeof(float), nz, fin1) != nz) {
				printf("Error reading trace from file 1!!\n");
	    		exit(1);
			}
		}

	//Mute shallow section (with taper)
	int z_mute = 20;
        for (k = 0; k < nh; k++)
                for (j = 0; j < nx; j++) {
                        int i0 = k*nz*nx + j*nz;
                        for (i = 0; i < z_mute; i++) {
                                float arg = (float) (z_mute - i)/z_mute * PI/2;
                                float taper = pow(cosf(arg),4);
                                IN[i0+i] *= taper;
                        }
		}

	//Input velocity
        for (j = 0; j < nx; j++) {
                int j0 = j*nz;		
                if (fread(&VEL[j0], sizeof(float), nz, fin2) != nz) {
                        printf("Error reading trace from file 2!!\n");
                        exit(1);
                }
	}
	//Input density
        for (j = 0; j < nx; j++) {
                int j0 = j*nz;
                if (fread(&DEN[j0], sizeof(float), nz, fin3) != nz) {
                        printf("Error reading trace from file 3!!\n");
                        exit(1);
                }
        }

	printf("Volumes In\n");
		
//################################# Differentiate or Integrate ####################################	

	//Differentiation or Integration via Fourier transform
	//Band-Pass filter is applied
	DerivIntegral(IN, OUT, nh, nx, nz, dz, power, band);
	
	printf("Derivation/Integration completed\n");

//################################### multiply by 32*DEN^2*VEL^6 ##################################

        // Multiplication by a physical scale factor
        for (k = 0; k < nh; k++)
                for (j = 0; j < nx; j++) {
                        int i0 = k*nz*nx + j*nz;
                	int j0 = j*nz;
                        for (i = 0; i < nz; i++)
//                                OUT[i0+i] *= 32.0 * powf(DEN[j0+i], 2) * powf(VEL[j0+i], 6);
                                  OUT[i0+i] *= 32.0 * powf(VEL[j0+i], 6);
//                                OUT[i0+i] *= 32.0 * powf(DEN[j0+i], 2) * powf(VEL[j0+i], 10);

		}

        printf("Scalar Multiply\n");
	
//########################################## Output Image #########################################
	
	for (k = 0; k < nh; k++)
		for (j = 0; j < nx; j++) {
			int i0 = k*nz*nx + j*nz;
   			if(fwrite(&OUT[i0], sizeof(float), nz, fout) != nz) {
				printf("Error writing trace!!\n");
    			exit(1);
    			}
		}
		
	printf("Volume Out\n");

	fclose (fin1);
	fclose (fin2);
	fclose (fin3);
	fclose (fout);

	free(IN);
        free(OUT);
	free(VEL);
	free(DEN);

}
