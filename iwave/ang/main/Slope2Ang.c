//#######################################################################################################
// FILE       :      Slope2Ang.c
// WRITER     :      Raanan Dafny
// EXERCISE   :      Interpolate slope to angle axis.
// DESCRIPTION:      Input ADCIG slope gathers.
//					 Interpolate slope -> angle.
//				     Output ADCIG angle gathers.
// HEADER FILES:     memory_alloc.h
//					 Interp.h
// INPUT FILES :     slope gathers
// OUTPUT FILE :     angle gathers
//#######################################################################################################

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "memory_alloc.h"
#include "Interp.h"

#define PI 3.141592654


int main (int argc, char *argv[]) {

	int    i, j, k;
    float  ***IN, ***OUT;
    float  *trace, *trace_p;
    FILE   *fin, *fout;

	float  deg2rad = PI / 180.0;
    
    //Angle Volume parameters
   	int    np          = atoi (argv[1]);
	int    nx          = atoi (argv[2]);
    int    nz          = atoi (argv[3]);
	int    nang        = atoi (argv[4]);
	float  dp          = atof (argv[5]);
    float  dang        = atof (argv[6]);
    char  *fileNameIN  = argv[7];
    char  *fileNameOUT = argv[8];
	
	//I/O files
    if((fin = fopen(fileNameIN,"rb")) == NULL) {
      	printf("Error opening input file !!\n");
      	exit(1);
    }
    
	fout = fopen(fileNameOUT,"wb");
	
	
//####################################### Memory Allocation #######################################

	//Memory allocation
	trace   = alloc_array(nz);
	trace_p = alloc_array(np);
	IN      = alloc_cube  (nx, nz, np);
	OUT     = alloc_cube  (nx, nang, nz);
	
	printf("Memory allocated\n");
	
	
//######################################### Input gathers #########################################

	//Input slope gathers
	for(k = 0; k < nx; k++)
		for(j = 0; j < np; j++) {
	  		if (fread(&trace[0], sizeof(float), nz, fin) != nz) {
	   	 		printf("Error reading x %d p %d !!\n",k,j);
	   	 		exit(1);
	    	}
	    	for(i = 0; i < nz; i++)
	    		IN[k][i][j] = trace[i];
	    }

	printf("Slope gathers input\n");


//######################################### Interpolation #########################################

	//Build slope trace
	for(j = 0; j < np; j++)
		trace_p[j] = atanf((-1*np/2 + j) * dp);
	
	//Interpulate
	for(k = 0; k < nx; k++)
		for(i = 0; i < nz; i++)
			for(j = 0; j < nang; j++) {
				float ang = (-1*nang/2 + j) * dang * deg2rad;
				OUT[k][j][i] = LG4PInterp(IN[k][i] , trace_p , ang , np);
			}


//########################################### Output HOCIG ########################################

	for(k = 0; k < nx; k++)
		for(j = 0; j < nang; j++)
   			if(fwrite(&OUT[k][j][0], sizeof(float), nz, fout) != nz) {
	    		printf("Error writing x %d ang %d !!\n",k,j);
    			exit(1);
    		}
	
	printf("Angle gathers Output\n");


	fclose(fin);
	fclose(fout);
	
	free(trace);
	free(trace_p);
	free(IN);
	free(OUT);
	
}

