//#######################################################################################################
// FILE       :      Wm.c
// WRITER     :      Raanan Dafny
// EXERCISE   :      Approximate inverse - VERSION_1 (Hou and Symes, 2015)
// DESCRIPTION:      Model side weight: Wm = -16*DEN^2*VEL^7*sqrt(-Laplac(xz))*sqrt(-Laplac(hz)).
//					 Laplacian operator via FFT.
// HEADER FILES:     memory_alloc.h
//					 FFT.h
// INPUT FILES :     mymig
// OUTPUT FILE :     mymigWm
//#######################################################################################################

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "memory_alloc.h"
#include "FFT.h"


int main (int argc, char *argv[]) {

	int   i, j, k, m, n;
 	float ***IN;
 	float **VEL, **DEN;
    stFFT *pFFT;
	FILE  *fin1, *fin2, *fin3, *fout;
	
	int    nh          = atoi (argv[1]);
	int    nx          = atoi (argv[2]);
    int    nz          = atoi (argv[3]);
    float  dh          = atof (argv[4]);
    float  dx          = atof (argv[5]);
    float  dz          = atof (argv[6]);
    //    float  f1          = atof (argv[7]);
    //    float  f2          = atof (argv[8]);
    //    float  f3          = atof (argv[9]);
    //    float  f4          = atof (argv[10]);
    char  *fileNameIN  = argv[11];
    char  *fileNameOUT = argv[12];
    char  *fileNameVEL = argv[13];
    char  *fileNameDEN = argv[14];

        //Band Pass filter
    //        float band[4] = {f1, f2, f3, f4};
        
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
	IN  = alloc_cube   (nh, nx, nz);
	VEL = alloc_matrix (nx, nz);
        DEN = alloc_matrix (nx, nz);
	
	//Allocate FFT struct
	pFFT = (stFFT*) calloc (1, sizeof(stFFT));

	//Allocate auxiliary arrays for FFTs functions
	int nbigger = (nx > nh) ? nx : nh;
	AllocAuxFFT(pFFT ,nbigger ,nbigger ,nz);

	printf("Memory Allocated\n");

//######################################### Input Volumes #########################################
	
	//Input image
	for (k = 0; k < nh; k++)
		for (j = 0; j < nx; j++)
			if (fread(&IN[k][j][0], sizeof(float), nz, fin1) != nz) {
				printf("Error reading trace from file 1!!\n");
	    		exit(1);
			}
	//Mute shallow section (with taper)
	int z_mute = 20;
	for (k = 0; k < nh; k++)
		for (j = 0; j < nx; j++)
			for (i = 0; i < z_mute; i++) {
				float arg = (float) (z_mute - i)/z_mute * PI/2;
				float taper = pow(cosf(arg),4);
				IN[k][j][i] *= taper;
			}
					
	//Input velocity
	for (j = 0; j < nx; j++)
		if (fread(&VEL[j][0], sizeof(float), nz, fin2) != nz) {
			printf("Error reading trace from file 2!!\n");
	   		exit(1);
		}
	//Input density
        for (j = 0; j < nx; j++)
                if (fread(&DEN[j][0], sizeof(float), nz, fin3) != nz) {
                        printf("Error reading trace from file 3!!\n");
                        exit(1);
                }

	printf("Volumes In\n");
		
//####################################### sqrt(Dx^2+Dz^2) #########################################

	//Copy into a complex array for FFT
	CopyForFFT(pFFT, IN, nh, nx, nz);
	
	//2D-FFT of volume projections
	FFT2DVol(pFFT, nh);
	
	//Frequency vectors
	Calc2DFreq(pFFT, dx, dz);

	//Multiply by sqrt(kx^2+kz^2)
   	for(m = 0; m < pFFT->n2_nfft; m++)
		for(n = 0; n < pFFT->n1_nfft/2; n++) {
			float sqrtLaplac = sqrt(Laplac2D(pFFT, m, n));
			for (k = 0; k < nh; k++)
				pFFT->M[k][m][n] *= sqrtLaplac;
		}

	//2D-Inverse FFT of volume projections
	InvFFT2DVol(pFFT, nh);
	
	//Copy complex array of FFT back to original volume 
	CopyFromFFT(pFFT, IN, nh, nx, nz);

	//zero-out arrays
	ZeroOutAuxFFT(pFFT, nh);

	printf("Laplac xz\n");

//####################################### sqrt(Dh^2+Dz^2) #########################################

	//Copy into a complex array for FFT (transpose axis 2 ans 3)
	CopyForFFT_transp(pFFT, IN, nh, nx, nz);
	
	//2D-FFT of volume projections
	FFT2DVol(pFFT, nx);
	
	//Frequency vectors
	Calc2DFreq(pFFT, dh, dz);
	
	//Multiply by sqrt(kh^2+kz^2)
   	for(m = 0; m < pFFT->n2_nfft; m++)
		for(n = 0; n < pFFT->n1_nfft/2; n++) {
			float sqrtLaplac = sqrt(Laplac2D(pFFT, m, n));
			for (k = 0; k < nx; k++)
				pFFT->M[k][m][n] *= sqrtLaplac;
		}
	
	//2D-Inverse FFT of volume projections
	InvFFT2DVol(pFFT, nx);
	
	//Copy complex array of FFT back to original volume (transpose axis 2 ans 3)
	CopyFromFFT_transp(pFFT, IN, nx, nh, nz);

	printf("Laplac hz\n");

//################################# multiply by -16*DEN^2*VEL^7 ###################################

	// Multiplication by a physical scale factor
	for (k = 0; k < nh; k++)
		for (j = 0; j < nx; j++)
			for (i = 0; i < nz; i++)
				IN[k][j][i] *= -16.0 * powf(DEN[j][i], 2) * powf(VEL[j][i], 7);

	printf("Scalar Multiply\n");

//########################################## Output Image #########################################

	for (k = 0; k < nh; k++)
		for (j = 0; j < nx; j++)
   			if(fwrite(&IN[k][j][0], sizeof(float), nz, fout) != nz) {
				printf("Error writing trace!!\n");
    			exit(1);
    		}

	printf("Volume Out\n");


	fclose (fin1);
	fclose (fin2);
        fclose (fin3);
	fclose (fout);

	free(IN);
	free(VEL);
        free(DEN);
	free(pFFT);
	

}
