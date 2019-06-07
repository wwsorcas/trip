//#######################################################################################################
// HEADER:		Intrep.h
// DESCRIPTION:	Interpolation functions
// CONTENT:		Linear interpolation
//				Inverse linear interpolation
//				La-Grange 4 points interpolation function
//				Cubic Spline interpulation
//#######################################################################################################


//*******************************************************************************************************

//Linear interpolation
void LinearInterp (float* v, float* d, float* delta, int N, int M , int* NM) {

	int i;
	for(i =0; i < N; i++) {

		float del   = delta[i];
		float delm1 = 1. - del;
		int   nm    = NM[i];

		d[i] = delm1 * v[nm] + del * v[nm + 1];
		
	}

}

//*******************************************************************************************************

//Inverse linear interpolation
void InvLinearInterp (float* d, float* v, float* delta, int N, int M, int* NM) {

	int i;
	memset((char *)v, 0, M * sizeof(float));

	for(i = 0; i < N; i++){

		float del   = delta[i];
		float delm1 = 1. - del;
		int   nm = NM[i];

		v[nm]     += delm1 * d[i];
		v[nm + 1] += del   * d[i];

	}

}

//*******************************************************************************************************

//La-Grange 4 points interpolation function
float LG4PInterp (float*P, float*x, float t, int samples) {

	int   i, j, n;
	float sum=0, multy, help[4], x_help[4];

	// If out of boundry
	if(t > x[samples-1])
		return P[samples-1];
	else if (t < x[0])
		return P[0];
	//If not, choose the 4 closest points to the interpulation point (t)
	else {
		
		n=0;
		while ((x[n] < t) && (n < samples-1))
			n++;

		if (n == 0) {
			for (i = 0; i < 4; i++) {
				help[i]   = P[i];
				x_help[i] = x[i];
			}
		}
		else if ((n > 0) && (n < samples-1)) {
			n -= 1;
			for (i = 0; i < 4; i++) {
				help[i]   = P[n+i];
				x_help[i] = x[n+i];
			}
		}
		else if (n == samples-1) {
			n -= 3;
			for (i = 0; i < 4; i++) {
				help[i]   = P[n+i];
				x_help[i] = x[n+i];
			}
		}

		//Interpolate with the 4 closest points
		for (i = 0; i < 4; i++) {
			multy = 1.f;
			for (j = 0; j < 4; j++) {
				if(j != i)
					multy *= ((t - x_help[j]) / (x_help[i] - x_help[j]));
			}
			sum += multy * help[i];
		}
		return sum;
		
	}
	
}

//*******************************************************************************************************

//La-Grange 4 points interpolation complex function
float complex LG4PInterpComplex (float complex*P, float*x, float t, int samples) {

	int   i, j, n;
	float multy, x_help[4];
	float complex sum=0, help[4];

	// choose the 4 closest points to the interpulation point (t)
	if (t > x[0] && t < x[samples-1]) {
		
		n=0;
		while ((x[n] < t) && (n < samples-1))
			n++;
		if (n == 0) {
			for (i = 0; i < 4; i++) {
				help[i]   = P[i];
				x_help[i] = x[i];
			}
		}
		else if ((n > 0) && (n < samples-1)) {
			n -= 1;
			for (i = 0; i < 4; i++) {
				help[i]   = P[n+i];
				x_help[i] = x[n+i];
			}
		}
		else if (n == samples-1) {
			n -= 3;
			for (i = 0; i < 4; i++) {
				help[i]   = P[n+i];
				x_help[i] = x[n+i];
			}
		}

		//Interpolate with the 4 closest points
		for (i = 0; i < 4; i++) {
			multy = 1.f;
			for (j = 0; j < 4; j++) {
				if(j != i)
					multy *= ((t - x_help[j]) / (x_help[i] - x_help[j]));
			}
			sum += multy * help[i];
		}
		
	}
	
		return sum;
	
}

//*******************************************************************************************************

//Cubic spline interpulation

void cubicSpline (int p, float *x, float *fx, float **spline) {

	int i, j, k;
	float h[p-1], b[p-1], u[p-2], v[p-2], matrix[p-2][p-1], help1[p-1], help2, m[p];

	for (i = 0; i < p-1 ; i++) {
		h[i] = x[i+1] - x[i];
		b[i] = (fx[i+1]-fx[i]) / h[i];
	}

	for (i = 1; i < p-1; i++) {
		u[i-1] = 2*(h[i-1] + h[i]);
		v[i-1] = 6*(b[i] - b[i-1]);
	}

	//build the tri-diagonal matrix
	for (i = 0; i < p-2; i++) {
		for (j = 0; j < p-1; j++) {
			if (i == j)
				matrix[i][j] = u[i];
			else if (j == i-1)
				matrix[i][j] = h[i];
			else if ((j == i+1) && (j != p-2))
				matrix[i][j] = h[i+1];
			else if(j == p-2)
				matrix[i][j] = v[i];
			else
				matrix[i][j] = 0.f;
		}
	}

	//solve the matrix by the Gauss elimination
	for (i = 0; i < p-2; i++) {
		for (j = 0; j < p-1; j++)
			help1[j] = matrix[i][j];
		for (j = 0; j < p-1; j++)
			matrix[i][j] /= help1[i];
		k = i+1;
		while (k < p-2) {
			help2 = matrix[k][i];
			for (j = 0; j < p-1; j++)
				matrix[k][j] -= help2 * matrix[i][j];
			k++;
		}
	}

	//the matrix solution(m)
	m[p-2] = matrix[p-3][p-2];
	m[p-1] = 0;
	for (i = p-3; i > 0; i--)
		m[i] = matrix[i-1][p-2] - (matrix[i-1][i] * m[i+1]);
	m[0] = 0;

	//calculate the spline coefficients (each row in the spline matrix is a different spline)
	for(i = 0; i < p-1; i++) {
		spline[i][0] = (m[i+1]-m[i])/(6*h[i]);
		spline[i][1] = (m[i]*x[i+1]-m[i+1]*x[i])/(2*h[i]);
		spline[i][2] = ((pow(x[i],2)*m[i+1]-pow(x[i+1],2)*m[i])/(2*h[i]))+((fx[i+1]-fx[i])/h[i])+((m[i]*h[i]-m[i+1]*h[i])/6);
		spline[i][3] = ((pow(x[i+1],3)*m[i]-pow(x[i],3)*m[i+1])/(6*h[i]))+(x[i]*(m[i+1]*h[i]/6-fx[i+1]/h[i]))+(x[i+1]*(fx[i]/h[i]-m[i]*h[i]/6));
	}
	
}

