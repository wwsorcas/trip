//#######################################################################################################
// HEADER:	matrix.h
// DESCRIPTION:	Matrix manipulations functions
// CONTENT:     Matrix transpose
//		Build ATA (including ATdt vector)
//              Matrix multiplication solver 
// 		Matrix-Vector multiplication solver
//	        Pivotint	
//		Gauss elimination function (square matrix solver)
//              LU Decomposition 
//              QR Decomposition (square matrix)
//              Eigen-Values (square matrix)
//              SVD Singular-Value-Decomposition
//		Conjugate gradient least squares solver
//#######################################################################################################

#include <cstd.h>

//*******************************************************************************************************

/*
 *Matrix transpose
 */
void Mat_Transpose(float **A , float **AT , int d1 , int d2)
{
    int i , j;
    
    for(i = 0 ; i < d2 ; i++)
	for(j = 0 ; j < d1 ; j++)
	    AT[i][j] = A[j][i];
}


//*******************************************************************************************************

/*
 *Matrix ATA Building
 */
void Mat_ATA(float **A, float **ATA , int d1 , int d2)
{
    int i , j , k;
    
	for(j = 0 ; j < d2 ; j++)
	  memset((char *)ATA[j], 0, d2*sizeof(float));

#ifdef _OPENMP
#pragma omp parallel default(shared) private(j, k, i) num_threads(8)
#endif
{
#ifdef _OPENMP
#pragma omp for
#endif
    //ATA
    for(j = 0 ; j < d2 ; j++) {
		for(k = 0 ; k < d2 ; k++) {
			for(i = 0 ; i < d1 ; i++) {
		   		ATA[j][k]+= A[i][j] * A[i][k];
		   	}
		}
	}
		   		
} // End of parallel region

}


//*******************************************************************************************************

/*
 *Spray ATA Matrix into d3 (construct block diagonal matrix)
 */
void Spray_ATA(float **ATA, int MofA, int d3)
{
    int i , j , k, l, i0;
	
	//Spray
	for(l = 1 ; l < d3 ; l++)
		for(k = 1 ; k < d3 ; k++)
			if (l == k) {
				i0 = k*MofA;
			    for(j = 0 ; j < MofA ; j++)
					for(i = 0 ; i < MofA ; i++)
						ATA[i0+i][i0+j] = ATA[i][j];
			}

}

//*******************************************************************************************************

/*
 *Vector ATd Building
 */
void Mat_ATd(float **A , float *d, float *ATd, int d1 , int d2)
{
    int i , j;
    
    memset((char *)ATd, 0, d2*sizeof(float));

    //ATd
    for(j = 0 ; j < d2 ; j++)
		for(i = 0 ; i < d1 ; i++)
		    ATd[j]+= A[i][j] * d[i];

}


//*******************************************************************************************************

/*
 *Matrix multiplication solver
 */
void Mat_Multi(float **A , float **B , float **C , int d1 , int d2 , int d3)
{
    int i , j , k;
  
    for(j = 0 ; j < d1 ; j++)
	for(k = 0 ; k < d3 ; k++)
	    {
		C[j][k] = 0.f;
		for(i = 0 ; i < d2 ; i++)
		    C[j][k]+= A[j][i] * B[i][k];
	    }

}


//*******************************************************************************************************

/*
 *Matrix-Vector multiplication solver
 */
void Mat_Vec_Multi(float **A , float *B , float *C , int d1 , int d2)
{
    int i , j;
  
    for(j = 0 ; j < d1 ; j++)
	    {
		C[j]=0.f;
		for(i = 0 ; i < d2 ; i++)
		    C[j]+= A[j][i] * B[i];
	    }

}


//*******************************************************************************************************

/*
 *2D Matrix inversion
 */
void Mat_Inv_2D(float **A , float **InvA)
{
	float det;
	
	det = A[0][0]*A[1][1] - A[0][1]*A[1][0];

	InvA[0][0] = A[1][1] / det;
	InvA[0][1] = -A[0][1] / det;
	InvA[1][0] = -A[1][0] / det;
	InvA[1][1] = A[0][0] / det;
}

//*******************************************************************************************************

/*
 *Pivotint -  biggest coefficients are on the main diagonal
 */
void Pivoting(float **matrix , float **Pivot , int d1 , int d2)
{
    int i , j , k , l;
    float n;

    for(j = 0 ; j < d2 ; j++)
	{
	    n = 0;
	    for(i = 0 ; i < d1 ; i++)
		{
		    if(fabs(matrix[i][j]) > n)
			{
			    n = fabs(matrix[i][j]);
			    l = i;
			}
		}
	    for(k = 0 ; k < d2 ; k++)
		{
		    Pivot[j][k] = matrix[l][k];
		    matrix[l][k] = 0;
		}
	}

} 


//*******************************************************************************************************

/*
 *Gauss elimination function (square matrix solver)
 */
void GE(float **matrix , float *sol_x , int r , int c)
{
    int i , j , k , l , m;
    float n , sum = 0;
    
    for(i = 0 ; i < r ; i++)
	{
	    for(j = 0 ; j < c ; j++)
		sol_x[j] = matrix[i][j];
	    for(j = 0 ; j < c ; j++)
		matrix[i][j]/= sol_x[i];
	    k = i + 1;
	    while(k < r)
		{
		    n = matrix[k][i];
		    for(j = 0 ; j < c ; j++)
			matrix[k][j]-= n * matrix[i][j];
		    k++;
		}
	}
    /*	
     *backward substitution
     */
    k = r - 1;
    for(i = 0 ; i < r ; i++)
	{
	    sum = 0;
	    l = 0;
	    m = r - 1;
	    for(j = 0 ; j < i ; j++)
		{
		    sum+= sol_x[l] * matrix[k][m];
		    l++;
		    m--;
		}
	    sol_x[i] = matrix[k][c-1] - sum;
	    k--;
	}
    
}


//*******************************************************************************************************

/*
 *LU Decoposition (lower triangular-upper triangular)
 */
void LU(float **matrix , int d1 , int d2)
{
    int i , j , k;
    float sum;

    for(j = 0; j < d2 ; j++)
	{
	    for(i = 0 ; i < d1 ; i++)
		{
		    sum = 0;
		    if(i<= j)
			{
			    for(k = 0 ; k < i ; k++)
				sum+= matrix[i][k]*matrix[k][j];
			    matrix[i][j]-= sum;
			}
		    else
			{
			    for(k = 0 ; k < j ; k++)
				sum+= matrix[i][k] * matrix[k][j];
			    matrix[i][j]=(1 / matrix[j][j]) * (matrix[i][j] - sum);
			}
		}
	}
    

}


//*******************************************************************************************************

/*
 *QR Decomposition (orthogonal-upper triangular)
 */
int QR_Decomp(float **matrix , float **Q , float **R , int d1)
{
    int i , j , k;
    float sum;
    
    for(j = 0 ; j < d1 ; j++)
	for(i = 0 ; i < d1 ; i++) 
	    Q[i][j] = matrix[i][j];
    
    for(j = 0 ; j < d1 ; j++)
	{
	    sum = 0;
	    for(i = 0 ; i < d1 ; i++) 
		sum+= (Q[i][j]*Q[i][j]);
	    R[j][j] = sqrt(sum);

	    if (R[j][j] == 0)
		{
		    printf("The matrix has lineary dependent columns!\n");
		    return(1);
		}
	    else
		{
		    for(i = 0 ; i < d1 ; i++) 
			Q[i][j]/= R[j][j];
		    for(k = j+1 ; k < d1 ; k++)
			{
			    sum=0;
			    for(i = 0 ; i < d1 ; i++) 
				sum+= (Q[i][j]*Q[i][k]);
			    R[j][k] = sum;
			    for(i = 0 ; i < d1 ; i++)  
				Q[i][k]-= (Q[i][j]*R[j][k]);
			}
		}
	}

    return 0;
}		    


//*******************************************************************************************************

/*
 *Eigen-Values
 */
int Eig_Val(float **matrix , float **V , float *eig , int d1)
{
    int i , j , k , iter , col_depend;
    float **Q , **QT , **R , **help;

    Q = alloc_matrix(d1,d1);
    QT = alloc_matrix(d1,d1);
    R = alloc_matrix(d1,d1);
    help = alloc_matrix(d1,d1);

    for(i = 0 ; i < d1 ; i++)
	help[i][i] = 1.;

    iter = 10;
    for(k = 0 ; k < iter ; k++) 
	{
	    col_depend = QR_Decomp(matrix , Q , R , d1);
	    if(col_depend == 0)
		{
		    Mat_Multi(help , Q , V , d1 , d1 , d1);
		    for(i = 0 ; i < d1 ; i++)
			for(j = 0 ; j < d1 ; j++)
				help[i][j] = V[i][j];
		    Mat_Transpose(Q , QT , d1 , d1);
		    Mat_Multi(QT , matrix , R , d1 , d1 , d1);
		    Mat_Multi(R , Q , matrix , d1 , d1 , d1);
		}
	    else break;
	}
    
    for(i = 0 ; i < d1 ; i++) 
	eig[i] = matrix[i][i];


    free(Q);
    free(QT);
    free(R);
    free(help);

    return col_depend; 

}


//*******************************************************************************************************

/*
 *SVD (Singular Value Decomposition)
 */
void SVD(int m , int n , int withu , int withv , float eps , float tol , float **a , float *q , float **u , float **v)
{
  int i,j,k,l,l1,iter;//,retval;
	float c,f,g,h,s,x,y,z;
	float *e;

	e = (float *)calloc(n,sizeof(float));
	//	retval = 0;

/* Copy 'a' to 'u' */    
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++)
			u[i][j] = a[i][j];
	}
/* Householder's reduction to bidiagonal form. */
	g = x = 0.0;    
	for (i=0;i<n;i++) {
		e[i] = g;
		s = 0.0;
		l = i+1;
		for (j=i;j<m;j++)
			s += (u[j][i]*u[j][i]);
		if (s < tol)
			g = 0.0;
		else {
			f = u[i][i];
			g = (f < 0) ? sqrt(s) : -sqrt(s);
			h = f * g - s;
			u[i][i] = f - g;
			for (j=l;j<n;j++) {
				s = 0.0;
				for (k=i;k<m;k++)
					s += (u[k][i] * u[k][j]);
				f = s / h;
				for (k=i;k<m;k++)
					u[k][j] += (f * u[k][i]);
			} /* end j */
		} /* end s */
		q[i] = g;
		s = 0.0;
		for (j=l;j<n;j++)
			s += (u[i][j] * u[i][j]);
		if (s < tol)
			g = 0.0;
		else {
			f = u[i][i+1];
			g = (f < 0) ? sqrt(s) : -sqrt(s);
			h = f * g - s;
			u[i][i+1] = f - g;
			for (j=l;j<n;j++) 
				e[j] = u[i][j]/h;
			for (j=l;j<m;j++) {
				s = 0.0;
				for (k=l;k<n;k++) 
					s += (u[j][k] * u[i][k]);
				for (k=l;k<n;k++)
					u[j][k] += (s * e[k]);
			} /* end j */
		} /* end s */
		y = fabs(q[i]) + fabs(e[i]);                         
		if (y > x)
			x = y;
	} /* end i */

/* accumulation of right-hand transformations */
	if (withv) {
		for (i=n-1;i>=0;i--) {
			if (g != 0.0) {
				h = u[i][i+1] * g;
				for (j=l;j<n;j++)
					v[j][i] = u[i][j]/h;
				for (j=l;j<n;j++) {
					s = 0.0;
					for (k=l;k<n;k++) 
						s += (u[i][k] * v[k][j]);
					for (k=l;k<n;k++)
						v[k][j] += (s * v[k][i]);

				} /* end j */
			} /* end g */
			for (j=l;j<n;j++)
				v[i][j] = v[j][i] = 0.0;
			v[i][i] = 1.0;
			g = e[i];
			l = i;
		} /* end i */
 
	} /* end withv, parens added for clarity */

/* accumulation of left-hand transformations */
	if (withu) {
		for (i=n;i<m;i++) {
			for (j=n;j<m;j++)
				u[i][j] = 0.0;
			u[i][i] = 1.0;
		}
	}
	if (withu) {
		for (i=n-1;i>=0;i--) {
			l = i + 1;
			g = q[i];
			for (j=l;j<m;j++)  /* upper limit was 'n' */
				u[i][j] = 0.0;
			if (g != 0.0) {
				h = u[i][i] * g;
				for (j=l;j<m;j++) { /* upper limit was 'n' */
					s = 0.0;
					for (k=l;k<m;k++)
						s += (u[k][i] * u[k][j]);
					f = s / h;
					for (k=i;k<m;k++) 
						u[k][j] += (f * u[k][i]);
				} /* end j */
				for (j=i;j<m;j++) 
					u[j][i] /= g;
			} /* end g */
			else {
				for (j=i;j<m;j++)
					u[j][i] = 0.0;
			}
			u[i][i] += 1.0;
		} /* end i*/
	} /* end withu, parens added for clarity */

/* diagonalization of the bidiagonal form */
	eps *= x;
	for (k=n-1;k>=0;k--) {
		iter = 0;
test_f_splitting:
		for (l=k;l>=0;l--) {
			if (fabs(e[l]) <= eps) goto test_f_convergence;
			if (fabs(q[l-1]) <= eps) goto cancellation;
		} /* end l */

/* cancellation of e[l] if l > 0 */
cancellation:
		c = 0.0;
		s = 1.0;
		l1 = l - 1;
		for (i=l;i<=k;i++) {
			f = s * e[i];
			e[i] *= c;
			if (fabs(f) <= eps) goto test_f_convergence;
			g = q[i];
			h = q[i] = sqrt(f*f + g*g);
			c = g / h;
			s = -f / h;
			if (withu) {
				for (j=0;j<m;j++) {
					y = u[j][l1];
					z = u[j][i];
					u[j][l1] = y * c + z * s;
					u[j][i] = -y * s + z * c;
				} /* end j */
			} /* end withu, parens added for clarity */
		} /* end i */
test_f_convergence:
		z = q[k];
		if (l == k) goto convergence;

/* shift from bottom 2x2 minor */
		iter++;
		if (iter > 30) {
		  //			retval = k;
			break;
		}
		x = q[l];
		y = q[k-1];
		g = e[k-1];
		h = e[k];
		f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2*h*y);
		g = sqrt(f*f + 1.0);
		f = ((x-z)*(x+z) + h*(y/((f<0)?(f-g):(f+g))-h))/x;
/* next QR transformation */
		c = s = 1.0;
		for (i=l+1;i<=k;i++) {
			g = e[i];
			y = q[i];
			h = s * g;
			g *= c;
			e[i-1] = z = sqrt(f*f+h*h);
			c = f / z;
			s = h / z;
			f = x * c + g * s;
			g = -x * s + g * c;
			h = y * s;
			y *= c;
			if (withv) {
				for (j=0;j<n;j++) {
					x = v[j][i-1];
					z = v[j][i];
					v[j][i-1] = x * c + z * s;
					v[j][i] = -x * s + z * c;
				} /* end j */
			} /* end withv, parens added for clarity */
			q[i-1] = z = sqrt(f*f + h*h);
			c = f/z;
			s = h/z;
			f = c * g + s * y;
			x = -s * g + c * y;
			if (withu) {
				for (j=0;j<m;j++) {
					y = u[j][i-1];
					z = u[j][i];
					u[j][i-1] = y * c + z * s;
					u[j][i] = -y * s + z * c;
				} /* end j */
			} /* end withu, parens added for clarity */
		} /* end i */
		e[l] = 0.0;
		e[k] = f;
		q[k] = x;
		goto test_f_splitting;
convergence:
		if (z < 0.0) {
/* q[k] is made non-negative */
			q[k] = - z;
			if (withv) {
				for (j=0;j<n;j++)
					v[j][k] = -v[j][k];
			} /* end withv, parens added for clarity */
		} /* end z */
	} /* end k */
	
	free(e);
}



//*******************************************************************************************************

/*
 *Conjugate gradient least squares solver: A*x=d
 */
int CG(float* x, float* r, float* v, float* Av, float* d, float** A, int N, int itrmaxCg, float eps)
{

    int i, j, itr;
    //double rs, rr0, vAv, rr, rAv, alpha, beta;
    double rs, rr0, vAv, rr, alpha, beta;

	//Initiation: r = d-Ax (x=0) ; rr = r*r ; v = r
    rr0 = 0.;
    for(i = 0 ; i < N ; ++i){
    	r[i] = d[i];
		rr0 += (double)(r[i] * r[i]);
		v[i] = r[i];
		x[i] = 0.;
    }
    rr = rr0;

	//If right-hand-side is null
    if(!rr){
		printf("\n from CG: rr = 0 ");
		return 0;
    }

	//CG iterations
    for(itr = 1 ; itr < itrmaxCg ; ++itr){
		vAv = 0.;
#ifdef _OPENMP
#pragma omp parallel default(shared) private(i, j) num_threads(8)
#endif
{
#ifdef _OPENMP
#pragma omp for
#endif
		//Set: Av = A*v ; vAv = v*Av
		for(i = 0 ; i < N ; ++i){
		    Av[i] = 0.;
		    for(j = 0 ; j < N ; ++j){
				Av[i] += A[i][j] * v[j];
	    	}
		}
		
} // End of parallel region

		for(i = 0 ; i < N ; ++i)
	    	vAv += (double)(v[i] * Av[i]);

		//Avoide devision by zero
		if(!vAv){
	    	printf("\n from CG: vAv = 0 ");
	    	return 0;
        }

		// Calculate alpha
		alpha = rr / vAv;
		//Save rr of cuurent iteration
		rs  = rr;
		rr  = 0.;
		//rAv = 0.;
		//Update: x = x+alpha*v ; r = r-alpha*Av
		for(i = 0 ; i < N ; ++i){
		    x[i] += alpha * v[i];
		    r[i] -= alpha * Av[i];
		    rr  += (double)(r[i] * r[i]);
		    //rAv += (double)(r[i] * Av[i]);
		}
		//Reset the residual every ten iterations
		if(itr/10*10 == itr){
		  //rAv = 0.;
		    rr  = 0.;
#ifdef _OPENMP
#pragma omp parallel default(shared) private(i, j) num_threads(8)
#endif
{
#ifdef _OPENMP
#pragma omp for
#endif
		    //Set: r = d-Ax
		    for(i = 0; i < N; ++i){
				r[i] = d[i];
				for(j = 0 ; j < N ; ++j){
				    r[i] -= A[i][j] * x[j];
				}
	    	}
	    	
} // End of parallel region	

		    for(i = 0; i < N; ++i){
				rr  += (double)(r[i] * r[i]);
				//rAv += (double)(r[i] * Av[i]);
	    	}
		}
		printf("CGitr = %d , residue = %lf\n" , itr , rr/rr0);

		//Check if residual is small enough
		if(rr/rr0 < eps) return 1;
		
		//If not, Set beta
		beta = rr / rs;
		//Set for next iteration: v = r+beta*v
		for(i = 0 ; i < N ; ++i)
			v[i] = r[i] + beta * v[i];
    }
    
    return 0;
    
}


//*******************************************************************************************************

/*
 *Conjugate gradient least squares solver: A*x=d (For AVA inversion)
 */
int CG_AVA(float* x, float* r, float* v, float* Av, float* d, float** A, int N, int NofA, int itrmaxCg, float eps)
{

    int i, j, itr;
    int i0, j0, ires, jres;
    double rs, rr0, vAv, rr, alpha, beta;
    //    double rs, rr0, vAv, rr, rAv, alpha, beta;

	//Initiation: r = d-Ax (x=0) ; rr = r*r ; v = r
    rr0 = 0.;
    for(i = 0 ; i < N ; ++i){
    	r[i] = d[i];
		rr0 += (double)(r[i] * r[i]);
		v[i] = r[i];
		x[i] = 0.;
    }
    rr = rr0;

	//If right-hand-side is null
    if(!rr){
		printf("\n from CG: rr = 0 ");
		return 0;
    }

	//CG iterations
    for(itr = 1 ; itr < itrmaxCg ; ++itr){
		vAv = 0.;
#ifdef _OPENMP
#pragma omp parallel default(shared) private(i, j, i0, j0, ires, jres) num_threads(8)
#endif
{
#ifdef _OPENMP
#pragma omp for
#endif
		//Set: Av = A*v ; vAv = v*Av
		for(i = 0 ; i < N ; ++i){
		    Av[i] = 0.;
	    	i0 = (int)(i/NofA);
		    for(j = 0 ; j < N ; ++j){
		    	j0 = (int)(j/NofA);
		    	//A is a block diagonal matrix
		    	if (i0 == j0) {
		    		ires = i%NofA;
		    		jres = j%NofA;
					Av[i] += A[ires][jres] * v[j];
				}
	    	}
	    }

} // End of parallel region

		for(i = 0 ; i < N ; ++i)
	    	vAv += (double)(v[i] * Av[i]);

		//Avoide devision by zero
		if(!vAv){
	    	printf("\n from CG: vAv = 0 ");
	    	return 0;
        }

		// Calculate alpha
		alpha = rr / vAv;
		//Save rr of cuurent iteration
		rs  = rr;
		rr  = 0.;
		//		rAv = 0.;
		//Update: x = x+alpha*v ; r = r-alpha*Av
		for(i = 0 ; i < N ; ++i){
		    x[i] += alpha * v[i];
		    r[i] -= alpha * Av[i];
		    rr  += (double)(r[i] * r[i]);
		    //		    rAv += (double)(r[i] * Av[i]);
		}
		//Reset the residual every ten iterations
		if(itr/10*10 == itr){
		  //		    rAv = 0.;
		    rr  = 0.;
#ifdef _OPENMP
#pragma omp parallel default(shared) private(i, j, i0, j0, ires, jres) num_threads(8)
#endif
{
#ifdef _OPENMP
#pragma omp for
#endif
		    //Set: r = d-Ax
		    for(i = 0; i < N; ++i){
				r[i] = d[i];
				i0 = (int)(i/NofA);
				for(j = 0 ; j < N ; ++j){
		    		j0 = (int)(j/NofA);
		    		//A is a block diagonal matrix
		    		if (i0 == j0) {
						ires = i%NofA;
		    			jres = j%NofA;
				    	r[i] -= A[ires][jres] * x[j];
				    }
				}
			}

} // End of parallel region	

		    for(i = 0; i < N; ++i){
				rr  += (double)(r[i] * r[i]);
				//				rAv += (double)(r[i] * Av[i]);
	    	}
		}
		printf("CGitr = %d , residue = %lf\n" , itr , rr/rr0);

		//Check if residual is small enough
		if(rr/rr0 < eps) return 1;
		
		//If not, Set beta
		beta = rr / rs;
		//Set for next iteration: v = r+beta*v
		for(i = 0 ; i < N ; ++i)
			v[i] = r[i] + beta * v[i];
    }
    
    return 0;
    
}


//*******************************************************************************************************


