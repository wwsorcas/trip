//#######################################################################################################
// HEADER:	memory_alloc.h
// DESCRIPTION:	Memory allocation functions (real and complex)
// CONTENT:	1-D real array allocation
//		2-D real array allocation
//		3-D real array allocation
//		1-D complex array allocation
//		2-D complex array allocation
//		3-D complex array allocation
//#######################################################################################################

#include<complex.h>

//*******************************************************************************************************

	float* alloc_array(int size)
	{
		float *array;
	
		array=(float*)calloc(size,sizeof(float));
    		if(array==0) return 0;
	
		return array;
	}


//*******************************************************************************************************

	float** alloc_matrix(int ROW,int COL)
	{
		int i;
		float **matrix;
	
		matrix=(float**)calloc(ROW,sizeof(float*));
    		if(matrix==0) return 0;
    		for(i=0;i<ROW;i++)
		{
		    matrix[i]=(float*)calloc(COL,sizeof(float));
		    if (matrix[i]==0) return 0;
		}

		return matrix;	
	}


//*******************************************************************************************************

	float*** alloc_cube(int DEPTH,int ROW,int COL)
	{
		int i,j;
		float ***cube;
	
		cube=(float***)calloc(DEPTH,sizeof(float**));
  		if(cube==0) return 0;
  		for(i=0;i<DEPTH;i++)
    		{
      			cube[i]=(float**)calloc(ROW,sizeof(float*));
      			if(cube[i]==0) return 0;
      			for(j=0;j<ROW;j++)
			{
	  			cube[i][j]=(float*)calloc(COL,sizeof(float));
	  			if(cube[i][j]==0) return 0;
			}
    		}	

		return cube;	
	}


//*******************************************************************************************************	

	float**** alloc_4D(int d4,int DEPTH,int ROW,int COL)
	{
		int i,j,k;
		float ****cube;
	
		cube=(float****)calloc(d4,sizeof(float***));
  		//if(cube==0) return 0;
		for(k=0;k<d4;k++)
		{
			    cube[k]=(float***)calloc(DEPTH,sizeof(float**));
	       		    //if(cube[k]==0) return 0;
  			    for(i=0;i<DEPTH;i++)
    			    {
      			    	cube[k][i]=(float**)calloc(ROW,sizeof(float*));
      				//if(cube[k][i]==0) return 0;
      				for(j=0;j<ROW;j++)
				{
	  				cube[k][i][j]=(float*)calloc(COL,sizeof(float));
	  				//if(cube[k][i][j]==0) return 0;
				}
	    	  	    }
		}	

		return cube;	
	}
	

//*******************************************************************************************************	

	float complex* alloc_complex_array(int size)
	{
		float complex *array;
	
		array=(float complex*)calloc(size,sizeof(float complex));
    		if(array==0) return 0;
	
		return array;
	}


//*******************************************************************************************************

	float complex** alloc_complex_matrix(int ROW,int COL)
	{
		int i;
		float complex **matrix;
	
		matrix=(float complex**)calloc(ROW,sizeof(float complex*));
    		if(matrix==0) return 0;
    		for(i=0;i<ROW;i++)
		{
		    matrix[i]=(float complex*)calloc(COL,sizeof(float complex));
		    if (matrix[i]==0) return 0;
		}

		return matrix;	
	}


//*******************************************************************************************************

	float complex*** alloc_complex_cube(int DEPTH,int ROW,int COL)
	{
		int i,j;
		float complex ***cube;
	
		cube=(float complex***)calloc(DEPTH,sizeof(float complex**));
  		if(cube==0) return 0;
  		for(i=0;i<DEPTH;i++)
    		{
      			cube[i]=(float complex**)calloc(ROW,sizeof(float complex*));
      			if(cube[i]==0) return 0;
      			for(j=0;j<ROW;j++)
			{
	  			cube[i][j]=(float complex*)calloc(COL,sizeof(float complex));
	  			if(cube[i][j]==0) return 0;
			}
    		}	

		return cube;	
	}


//*******************************************************************************************************	

	float complex**** alloc_complex_4D(int d4,int DEPTH,int ROW,int COL)
	{
		int i,j,k;
		float complex ****cube;
	
		cube=(float complex****)calloc(d4,sizeof(float complex***));
  		//if(cube==0) return 0;
		for(k=0;k<d4;k++)
		{
			    cube[k]=(float complex***)calloc(DEPTH,sizeof(float complex**));
	       		    //if(cube[k]==0) return 0;
  			    for(i=0;i<DEPTH;i++)
    			    {
      			    	cube[k][i]=(float complex**)calloc(ROW,sizeof(float complex*));
      				//if(cube[k][i]==0) return 0;
      				for(j=0;j<ROW;j++)
				{
	  				cube[k][i][j]=(float complex*)calloc(COL,sizeof(float complex));
	  				//if(cube[k][i][j]==0) return 0;
				}
	    	  	    }
		}	

		return cube;	
	}
	
	
//*******************************************************************************************************	

