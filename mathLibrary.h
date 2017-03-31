#ifndef MATH_LIBRARY_H
#define MATH_LIBRARY_H 1

// http://www.mymathlib.com/

#include <stdlib.h>
#include <float.h>                           // required for DBL_EPSILON
#include <math.h>                           // required for fabs()

//#include "configuration.h"

typedef double scalar;
// typedef float scalar;

#define MIL 1024
#define MILION (MIL*MIL)

				//J:row, K: column, NK: number of columns
#define pos2Dto1D(J, K, NK) ( (J)*(NK) + (K) )
				//I:plano, J:row, K: column, NJ:number of rwos, NK: number of columns
#define pos3Dto1D(I, J, K, NJ, NK) ( (I)*(NJ)*(NK) + (J)*(NK) + (K) )

//#define SIGN(r) (((r)<0.0)?-1.0:1.0)

#ifndef MAX
	#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
	#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef PAYOFF_CALL
	#define PAYOFF_CALL(S,K) ((((S)-(K))>0)?((S)-(K)):0)
#endif
#ifndef PAYOFF_PUT
	#define PAYOFF_PUT(S,K) ((((K)-(S))>0)?((K)-(S)):0)
#endif

#ifndef LEN
	#define LEN(x) (sizeof(x)/sizeof(x[0]))
#endif

#ifndef ABS
	#define ABS(a) (((a) < (0.0)) ? (-a) : (a))
#endif

#ifndef FMAX
    #define FMAX(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef SQR
    #define SQR(a) ((a)*(a))
#endif

#ifndef SGN
    #define SGN(a) ((a)<0?-1:1)
#endif

#ifndef SIGN
    #define SIGN(a,b) (SGN(a)*SGN(b))
#endif

#ifndef RAND_UNIFORM
	#define RAND_UNIFORM ( (scalar)random() / ((scalar)RAND_MAX + 1.0) )
#endif

#ifndef PI
	#define PI 3.141592653589793238462643
#endif

/*******************RANDOM NUMBERS********************/
scalar randBoxMuller(){
	scalar  r = sqrt(-2.0*log(RAND_UNIFORM));
	scalar phi = 2.0*PI*RAND_UNIFORM;
	return r*cos(phi);
}//end randBoxMuller

inline void BoxMuller(scalar *u1, scalar *u2){
	scalar  r = sqrt(-2.0*log(*u1));
	scalar phi = 2.0*PI*(*u2);
	*u1 = r*cos(phi);
	*u2 = r*sin(phi);
}//end BoxMuller

/*******************RANDOM NUMBERS********************/

/*******************QR factorization--least squared solver*********************/
int qrdcmp(scalar *a, int m, int n, scalar *c, scalar *d) {
  int i,j,k;
  scalar scale,sigma,sum,tau;
  int sing = 0;

  for(k = 0; k < n; k++) {
    scale = 0.0;
    for (i = k; i < m; i++)
      scale = FMAX(scale,fabs(a[pos2Dto1D(i, k, n)]));
    if (scale == 0.0) {
      sing = 1;
      c[k] = d[k] = 0.0;
    } else {
      for(i = k; i < m; i++)
        a[pos2Dto1D(i, k, n)] /= scale;
      for(sum = 0.0, i = k; i < m; i++)
    sum += SQR(a[pos2Dto1D(i, k, n)]);
      sigma = sqrt(sum)*SGN(a[pos2Dto1D(k, k, n)]);
      a[pos2Dto1D(k, k, n)] += sigma;
      c[k] = sigma*a[pos2Dto1D(k, k, n)];
      d[k] = -scale*sigma;
      for (j = k+1; j < n; j++) {
        for(sum = 0.0, i = k; i < m; i++)
            sum += a[pos2Dto1D(i, k, n)]*a[pos2Dto1D(i, j, n)];
        tau = sum/c[k];
        for(i = k; i < m; i++)
            a[pos2Dto1D(i, j, n)] -= tau*a[pos2Dto1D(i, k, n)];
      }
    }
  }
  //d[n-1] = a[pos2Dto1D(n-1, n-1, n)];
  if (d[n-1] == 0.0)
    sing = 1;

  return sing;
}

/* solves Rx = b, based upon QR decomp */
void rsolv(scalar *a, int n, scalar *d, scalar *b) {
  int i,j;
  scalar sum;
  b[n-1] /= d[n-1];
  for (i = n-2; i >= 0; i--) {
    for(sum = 0.0, j = i+1; j < n; j++)
      sum += a[pos2Dto1D(i, j, n)]*b[j];
    b[i] = (b[i] - sum)/d[i];
  }
}

/* solves Ax = b, based upon QR decomp */
void qrsolv(scalar *a, int m, int n, scalar *c, scalar *d, scalar *b) {
  int i,j;
  scalar sum, tau;

  for (j = 0; j < n; j++) {
    for (sum = 0.0, i = j; i < m; i++)
      sum += a[pos2Dto1D(i, j, n)]*b[i];
    tau = sum/c[j];
    for (i = j; i < m; i++)
      b[i] -= tau*a[pos2Dto1D(i, j, n)];
  }
  rsolv(a, n, d, b);
}

int leastsq(scalar *a, int m, int n, scalar *b) {

  scalar *d, *c;

  c = (scalar *)malloc(n*sizeof(scalar));
  d = (scalar *)malloc(n*sizeof(scalar));

  if (qrdcmp(a, m, n, c, d))
    return -1;

  qrsolv(a, m, n, c, d, b);
  free(c);
  free(d);

  return 0;
}

/*******************QR factorization--least squared solver*********************/

/*******************USEFUL STATISTIC FUNCTIONS************************/

#define  A1  (-3.969683028665376e+01)
#define  A2   2.209460984245205e+02
#define  A3  (-2.759285104469687e+02)
#define  A4   1.383577518672690e+02
#define  A5  (-3.066479806614716e+01)
#define  A6   2.506628277459239e+00

#define  B1  (-5.447609879822406e+01)
#define  B2   1.615858368580409e+02
#define  B3  (-1.556989798598866e+02)
#define  B4   6.680131188771972e+01
#define  B5  (-1.328068155288572e+01)

#define  C1  (-7.784894002430293e-03)
#define  C2  (-3.223964580411365e-01)
#define  C3  (-2.400758277161838e+00)
#define  C4  (-2.549732539343734e+00)
#define  C5   4.374664141464968e+00
#define  C6   2.938163982698783e+00

#define  D1   7.784695709041462e-03
#define  D2   3.224671290700398e-01
#define  D3   2.445134137142996e+00
#define  D4   3.754408661907416e+00

#define P_LOW   0.02425
#define P_HIGH  0.97575

scalar norminv(scalar p){

	scalar x;
	scalar q, r, u, e;
	if ((0.0 < p )  && (p < P_LOW)){
	q = sqrt(-2.0*log(p));
	x = (((((C1*q + C2)*q + C3)*q + C4)*q + C5)*q + C6) / ((((D1*q + D2)*q + D3)*q + D4)*q + 1.0);
	}
	else{
			if ((P_LOW <= p) && (p <= P_HIGH)){
			q = p - 0.5;
			r = q*q;
			x = (((((A1*r + A2)*r + A3)*r + A4)*r + A5)*r + A6)*q /(((((B1*r + B2)*r + B3)*r + B4)*r + B5)*r + 1.0);
			}
			else{
					if ((P_HIGH < p)&&(p < 1.0)){
					q = sqrt(-2.0*log(1.0 - p));
					x = -(((((C1*q + C2)*q + C3)*q + C4)*q + C5)*q + C6) / ((((D1*q + D2)*q + D3)*q + D4)*q + 1.0);
					}
			}
	}

	if(( 0.0 < p)&&(p < 1.0)){
	e = 0.5*erfc(-x/sqrt(2.0)) - p;
	u = e*sqrt(2.0*PI)*exp(x*x/2.0);
	x = x - u/(1.0 + x*u/2.0);
	}

	return x;
}//end norminv

scalar logninv(scalar u, scalar mu, scalar sigma){

	return exp(sigma*norminv(u) + mu);
}//end logninv

scalar lower_incomplete_gamma(scalar x, scalar s){

	scalar sum = 0.0;

	int j;
	for(j = 0;j < 30; j++){
		sum += pow(x, j)/tgamma(s + j + 1.0);
	}

	return pow(x,s)*tgamma(s)*exp(-x)*sum;
}//end lower_incomplete_gamma

scalar besseli(scalar a, scalar x){

	unsigned int fact = 1;
	scalar sum = 0.0;

	int j;
	for(j=0;j<30;j++){
		sum += pow(x/2.0, 2.0*j + a)/(fact*tgamma(a + j + 1.0));
		fact *= j + 1;
	}
	return sum;
}//end besseli

scalar chi2cdf(scalar x, scalar k){

	return lower_incomplete_gamma(x/2.0, k/2.0)/tgamma(k/2.0);
}//end chi2cdf

scalar ncx2cdf(scalar x, scalar k, scalar l){

	unsigned int fact = 1;
	scalar sum = 0.0;

	int j;
	for(j = 0; j<30; j++){
		sum += pow(l/2.0, j)*chi2cdf(x, k + 2*j)/fact;
		fact *= j + 1;
	}
	return exp(-l/2.0)*sum;
}//end ncx2cdf
/*******************USEFUL STATISTIC FUNCTIONS************************/

////////////////////////////////////////////////////////////////////////////////
// File: jacobi_cyclic_method.c                                               //
// Routines:                                                                  //
//    Jacobi_Cyclic_Method                                                    //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Jacobi_Cyclic_Method                                                 //
//            (scalar eigenvalues[], scalar *eigenvectors, scalar *A, int n)  //
//                                                                            //
//  Description:                                                              //
//     Find the eigenvalues and eigenvectors of a symmetric n x n matrix A    //
//     using the Jacobi method. Upon return, the input matrix A will have     //
//     been modified.                                                         //
//     The Jacobi procedure for finding the eigenvalues and eigenvectors of a //
//     symmetric matrix A is based on finding a similarity transformation     //
//     which diagonalizes A.  The similarity transformation is given by a     //
//     product of a sequence of orthogonal (rotation) matrices each of which  //
//     annihilates an off-diagonal element and its transpose.  The rotation   //
//     effects only the rows and columns containing the off-diagonal element  //
//     and its transpose, i.e. if a[i][j] is an off-diagonal element, then    //
//     the orthogonal transformation rotates rows a[i][] and a[j][], and      //
//     equivalently it rotates columns a[][i] and a[][j], so that a[i][j] = 0 //
//     and a[j][i] = 0.                                                       //
//     The cyclic Jacobi method considers the off-diagonal elements in the    //
//     following order: (0,1),(0,2),...,(0,n-1),(1,2),...,(n-2,n-1).  If the  //
//     the magnitude of the off-diagonal element is greater than a treshold,  //
//     then a rotation is performed to annihilate that off-diagnonal element. //
//     The process described above is called a sweep.  After a sweep has been //
//     completed, the threshold is lowered and another sweep is performed     //
//     with the new threshold. This process is completed until the final      //
//     sweep is performed with the final threshold.                           //
//     The orthogonal transformation which annihilates the matrix element     //
//     a[k][m], k != m, is Q = q[i][j], where q[i][j] = 0 if i != j, i,j != k //
//     i,j != m and q[i][j] = 1 if i = j, i,j != k, i,j != m, q[k][k] =       //
//     q[m][m] = cos(phi), q[k][m] = -sin(phi), and q[m][k] = sin(phi), where //
//     the angle phi is determined by requiring a[k][m] -> 0.  This condition //
//     on the angle phi is equivalent to                                      //
//               cot(2 phi) = 0.5 * (a[k][k] - a[m][m]) / a[k][m]             //
//     Since tan(2 phi) = 2 tan(phi) / (1.0 - tan(phi)^2),                    //
//               tan(phi)^2 + 2cot(2 phi) * tan(phi) - 1 = 0.                 //
//     Solving for tan(phi), choosing the solution with smallest magnitude,   //
//       tan(phi) = - cot(2 phi) + sgn(cot(2 phi)) sqrt(cot(2phi)^2 + 1).     //
//     Then cos(phi)^2 = 1 / (1 + tan(phi)^2) and sin(phi)^2 = 1 - cos(phi)^2 //
//     Finally by taking the sqrts and assigning the sign to the sin the same //
//     as that of the tan, the orthogonal transformation Q is determined.     //
//     Let A" be the matrix obtained from the matrix A by applying the        //
//     similarity transformation Q, since Q is orthogonal, A" = Q'AQ, where Q'//
//     is the transpose of Q (which is the same as the inverse of Q).  Then   //
//         a"[i][j] = Q'[i][p] a[p][q] Q[q][j] = Q[p][i] a[p][q] Q[q][j],     //
//     where repeated indices are summed over.                                //
//     If i is not equal to either k or m, then Q[i][j] is the Kronecker      //
//     delta.   So if both i and j are not equal to either k or m,            //
//                                a"[i][j] = a[i][j].                         //
//     If i = k, j = k,                                                       //
//        a"[k][k] =                                                          //
//           a[k][k]*cos(phi)^2 + a[k][m]*sin(2 phi) + a[m][m]*sin(phi)^2     //
//     If i = k, j = m,                                                       //
//        a"[k][m] = a"[m][k] = 0 =                                           //
//           a[k][m]*cos(2 phi) + 0.5 * (a[m][m] - a[k][k])*sin(2 phi)        //
//     If i = k, j != k or m,                                                 //
//        a"[k][j] = a"[j][k] = a[k][j] * cos(phi) + a[m][j] * sin(phi)       //
//     If i = m, j = k, a"[m][k] = 0                                          //
//     If i = m, j = m,                                                       //
//        a"[m][m] =                                                          //
//           a[m][m]*cos(phi)^2 - a[k][m]*sin(2 phi) + a[k][k]*sin(phi)^2     //
//     If i= m, j != k or m,                                                  //
//        a"[m][j] = a"[j][m] = a[m][j] * cos(phi) - a[k][j] * sin(phi)       //
//                                                                            //
//     If X is the matrix of normalized eigenvectors stored so that the ith   //
//     column corresponds to the ith eigenvalue, then AX = X Lamda, where     //
//     Lambda is the diagonal matrix with the ith eigenvalue stored at        //
//     Lambda[i][i], i.e. X'AX = Lambda and X is orthogonal, the eigenvectors //
//     are normalized and orthogonal.  So, X = Q1 Q2 ... Qs, where Qi is      //
//     the ith orthogonal matrix,  i.e. X can be recursively approximated by  //
//     the recursion relation X" = X Q, where Q is the orthogonal matrix and  //
//     the initial estimate for X is the identity matrix.                     //
//     If j = k, then x"[i][k] = x[i][k] * cos(phi) + x[i][m] * sin(phi),     //
//     if j = m, then x"[i][m] = x[i][m] * cos(phi) - x[i][k] * sin(phi), and //
//     if j != k and j != m, then x"[i][j] = x[i][j].                         //
//                                                                            //
//  Arguments:                                                                //
//     scalar  eigenvalues                                                    //
//        Array of dimension n, which upon return contains the eigenvalues of //
//        the matrix A.                                                       //
//     scalar* eigenvectors                                                   //
//        Matrix of eigenvectors, the ith column of which contains an         //
//        eigenvector corresponding to the ith eigenvalue in the array        //
//        eigenvalues.                                                        //
//     scalar* A                                                              //
//        Pointer to the first element of the symmetric n x n matrix A. The   //
//        input matrix A is modified during the process.                      //
//     int     n                                                              //
//        The dimension of the array eigenvalues, number of columns and rows  //
//        of the matrices eigenvectors and A.                                 //
//                                                                            //
//  Return Values:                                                            //
//     Function is of type void.                                              //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     scalar A[N][N], scalar eigenvalues[N], scalar eigenvectors[N][N]       //
//                                                                            //
//     (your code to initialize the matrix A )                                //
//                                                                            //
//     Jacobi_Cyclic_Method(eigenvalues, (scalar*)eigenvectors,               //
//                                                          (scalar *) A, N); //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Jacobi_Cyclic_Method(scalar eigenvalues[], scalar *eigenvectors, scalar *A, int n)
{
   int /*row,*/ i, j, k, m;
   scalar *pAk, *pAm, *p_r, *p_e;
   scalar threshold_norm;
   scalar threshold;
   scalar tan_phi, sin_phi, cos_phi, tan2_phi, sin2_phi, cos2_phi;
   scalar sin_2phi, /*cos_2phi,*/ cot_2phi;
   scalar dum1;
   scalar dum2;
   scalar dum3;
   //scalar r;
   scalar max;

                  // Take care of trivial cases

   if ( n < 1) return;
   if ( n == 1) {
      eigenvalues[0] = *A;
      *eigenvectors = 1.0;
      return;
   }

          // Initialize the eigenvalues to the identity matrix.

   for (p_e = eigenvectors, i = 0; i < n; i++)
      for (j = 0; j < n; p_e++, j++)
         if (i == j) *p_e = 1.0; else *p_e = 0.0;
  
            // Calculate the threshold and threshold_norm.
 
   for (threshold = 0.0, pAk = A, i = 0; i < ( n - 1 ); pAk += n, i++) 
      for (j = i + 1; j < n; j++) threshold += *(pAk + j) * *(pAk + j);
   threshold = sqrt(threshold + threshold);
   threshold_norm = threshold * DBL_EPSILON;
   max = threshold + 1.0;
   while (threshold > threshold_norm) {
      threshold /= 10.0;
      if (max < threshold) continue;
      max = 0.0;
      for (pAk = A, k = 0; k < (n-1); pAk += n, k++) {
         for (pAm = pAk + n, m = k + 1; m < n; pAm += n, m++) {
            if ( fabs(*(pAk + m)) < threshold ) continue;

                 // Calculate the sin and cos of the rotation angle which
                 // annihilates A[k][m].

            cot_2phi = 0.5 * ( *(pAk + k) - *(pAm + m) ) / *(pAk + m);
            dum1 = sqrt( cot_2phi * cot_2phi + 1.0);
            if (cot_2phi < 0.0) dum1 = -dum1;
            tan_phi = -cot_2phi + dum1;
            tan2_phi = tan_phi * tan_phi;
            sin2_phi = tan2_phi / (1.0 + tan2_phi);
            cos2_phi = 1.0 - sin2_phi;
            sin_phi = sqrt(sin2_phi);
            if (tan_phi < 0.0) sin_phi = - sin_phi;
            cos_phi = sqrt(cos2_phi); 
            sin_2phi = 2.0 * sin_phi * cos_phi;
            //cos_2phi = cos2_phi - sin2_phi;

                     // Rotate columns k and m for both the matrix A 
                     //     and the matrix of eigenvectors.

            p_r = A;
            dum1 = *(pAk + k);
            dum2 = *(pAm + m);
            dum3 = *(pAk + m);
            *(pAk + k) = dum1 * cos2_phi + dum2 * sin2_phi + dum3 * sin_2phi;
            *(pAm + m) = dum1 * sin2_phi + dum2 * cos2_phi - dum3 * sin_2phi;
            *(pAk + m) = 0.0;
            *(pAm + k) = 0.0;
            for (i = 0; i < n; p_r += n, i++) {
               if ( (i == k) || (i == m) ) continue;
               if ( i < k ) dum1 = *(p_r + k); else dum1 = *(pAk + i);
               if ( i < m ) dum2 = *(p_r + m); else dum2 = *(pAm + i);
               dum3 = dum1 * cos_phi + dum2 * sin_phi;
               if ( i < k ) *(p_r + k) = dum3; else *(pAk + i) = dum3;
               dum3 = - dum1 * sin_phi + dum2 * cos_phi;
               if ( i < m ) *(p_r + m) = dum3; else *(pAm + i) = dum3;
            }
            for (p_e = eigenvectors, i = 0; i < n; p_e += n, i++) {
               dum1 = *(p_e + k);
               dum2 = *(p_e + m);
               *(p_e + k) = dum1 * cos_phi + dum2 * sin_phi;
               *(p_e + m) = - dum1 * sin_phi + dum2 * cos_phi;
            }
         }
         for (i = 0; i < n; i++)
            if ( i == k ) continue;
            else if ( max < fabs(*(pAk + i))) max = fabs(*(pAk + i));
      }
   }
   for (pAk = A, k = 0; k < n; pAk += n, k++) eigenvalues[k] = *(pAk + k); 
}

////////////////////////////////////////////////////////////////////////////////
// File: sort_eigenvalues.c                                                   //
// Routines:                                                                  //
//    Sort_Eigenvalues                                                        //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Sort_Eigenvalues(scalar eigenvalues[], scalar *eigenvectors, int n,  //
//                                                            int sort_order) //
//                                                                            //
//  Description:                                                              //
//     Sort the eigenvalues and corresponding eigenvectors. If sort_order is  //
//     positive, then the eigenvalues are sorted in ascending order and if    //
//     sort_order is negative, then the eigenvalues are sorted in descending  //
//     order.  The columns of the matrix eigenvectors are permuted so that    //
//     the ith eigenvalue has corresponding eigenvector in the ith column of  //
//     the matrix eigenvectors.                                               //
//                                                                            //
//  Arguments:                                                                //
//     scalar  eigenvalues  - Array, dimension n, containing the eigenvalues. //
//     scalar* eigenvectors - Matrix of eigenvectors, the ith column of which //
//                            contains an eigenvector corresponding to the    //
//                            ith eigenvalue in the array eigenvalues.        //
//     int     n            - The dimension of the array eigenvalues and the  //
//                            number of columns and rows of the matrix        //
//                            eigenvectors.                                   //
//     int     sort_order   - An indicator used to specify the order in which //
//                            the eigenvalues and eigenvectors are to be      //
//                            returned.  If sort_order > 0, then the          //
//                            eigenvalues are sorted in increasing order.     //
//                            If sort_order < 0, then the eigenvalues are     //
//                            sorted in decreasing order.  If sort_order = 0, //
//                            then the eigenvalues are returned in the order, //
//                            found.                                          //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     scalar scalar eigenvalues[N], scalar eigenvectors[N][N];               //
//     int sort_order = 1;                                                    //
//                                                                            //
//     (your code to calculate the eigenvalues and eigenvectors)              //
//                                                                            //
//     Sort_Eigenvalues(eigenvalues, (scalar*)eigenvectors, N,  sort_order);  //
//     printf(" The Eigenvalues are: \n"); ...                                //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Sort_Eigenvalues(scalar eigenvalues[], scalar *eigenvectors, int n,
                                                                int sort_order)
{
   int i,j;
   int m;
   scalar x, *pm, *pi;

   if (sort_order == 0) return;
   if (sort_order < 0) {
      for (i = 0; i < (n - 1); i++) {
         m = i;
         for (j = i+1; j < n; j++) if (eigenvalues[m] < eigenvalues[j]) m = j;
         if (m == i) continue;
         x = eigenvalues[m];
         eigenvalues[m] = eigenvalues[i];
         eigenvalues[i] = x;
         pm = eigenvectors + m;
         pi = eigenvectors + i;
         for (j = 0; j < n; pm += n, pi += n, j++) {
            x = *pm;
            *pm = *pi;
            *pi = x;
         }
      }
   }
   else {
      for (i = 0; i < (n - 1); i++) {
         m = i;
         for (j = i+1; j < n; j++) if (eigenvalues[m] > eigenvalues[j]) m = j;
         if (m == i) continue;
         x = eigenvalues[m];
         eigenvalues[m] = eigenvalues[i];
         eigenvalues[i] = x;
         pm = eigenvectors + m;
         pi = eigenvectors + i;
         for (j = 0; j < n; pm += n, pi += n, j++) {
            x = *pm;
            *pm = *pi;
            *pi = x;
         }
      }
   }
   return;
}

// A matriz simétrica de entrada de tamaño nxn, L matriz triangular inferior de salida con la descomposición de cholesky
// Valor de retorno; 1 si la matriz es definida positiva, 0 si no lo es
int Cholesky(scalar *A, int n, scalar *L) {
	int i,j,k;

	scalar aux;
 
	//inicializamos la matriz resultado L a ceros, para meter ceros en el triangulo superior
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			L[pos2Dto1D(i,j,n)] = 0.0;

    for (i = 0; i < n; i++)
        for (j = 0; j < (i+1); j++) {
            scalar s = 0;
            for (k = 0; k < j; k++)
                s += L[i * n + k] * L[j * n + k];
            //L[i * n + j] = (i == j) ?
            //               sqrt(A[i * n + i] - s) :
            //               (1.0 / L[j * n + j] * (A[i * n + j] - s));
			if (i==j) {
				aux = A[i * n + i] - s;
				if (aux <= 0.0)
					return 0;
				else
					L[i * n + j] = sqrt(aux);
			} else
				L[i * n + j] = (1.0 / L[j * n + j] * (A[i * n + j] - s));

        }

	return 1;
}

void diagonal(int N, scalar *A, scalar *b){

int i, j, k;
scalar temp = 0.0;
	
	for(i=0; i<N; i++){
		if(A[pos2Dto1D(i,i,N)]==0){
			for(j=0; j<N; j++){
				if(j==i) continue;
				if(A[pos2Dto1D(j,i,N)] != 0 && A[pos2Dto1D(i,j,N)] != 0){
					for(k=0; k<N; k++){
						temp = A[pos2Dto1D(j,k,N)];
						A[pos2Dto1D(j,k,N)] = A[pos2Dto1D(i,k,N)];
						A[pos2Dto1D(i,k,N)] = temp;
					}
					temp = b[j];
					b[j] = b[i];
					b[i] = temp;
					break;
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// File: set_diagonal.c                                                       //
// Routine(s):                                                                //
//    Set_Diagonal                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Set_Diagonal(scalar *A, scalar v[], int nrows, int ncols)            //
//                                                                            //
//  Description:                                                              //
//     Copy the vector v to the diagonal A[i][i], where 0 <= i <              //
//     min( nrows, ncols ).                                                   //
//     Note that v should be declared "scalar v[N]", N >= min( nrows, ncols ) //
//     in the calling routine.                                                //
//                                                                            //
//  Arguments:                                                                //
//     scalar *A    Pointer to the first element of the source matrix A.      //
//     scalar v[]   Source of the new diagonal for the matrix A.              //
//     int    nrows The number of rows matrix A.                              //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     scalar A[M][N],  v[N];                                                 //
//                                                                            //
//     (your code to initialize the matrix A and the vector v)                //
//                                                                            //
//     Set_Diagonal(&A[0][0], v, M, N);                                       //
//     printf("The matrix A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
void Set_Diagonal(scalar *A, scalar v[], int nrows, int ncols)
{
   int n;

   if (nrows < ncols) n = nrows; else n = ncols;

   for (; n > 0 ; A += (ncols + 1), n--)  *A = *v++;
}


void Set_Matrix(scalar *A, int nrows, int ncols, scalar value) {
	unsigned int i, j;
	
	for(i=0;i<nrows;i++)
		for(j=0;j<ncols;j++)
			A[pos2Dto1D(i,j,ncols)] = value;
}

void Set_MatrixInt(int *A, int nrows, int ncols, int value) {
    unsigned int i, j;

    for(i=0;i<nrows;i++)
        for(j=0;j<ncols;j++)
            A[pos2Dto1D(i,j,ncols)] = value;
}

void Set_Matrix3D(scalar *A, int i_size, int j_size, int k_size, scalar value) {
	unsigned int i, j, k;
	
	for(i=0;i<i_size;i++)
		for(j=0;j<j_size;j++)
			for(k=0;k<k_size;k++)
				A[pos3Dto1D(i,j,k, j_size, k_size)] = value;
}

void Set_Vector(scalar *v, int n, scalar value) {
	unsigned int i;
	
	for(i=0;i<n;i++)
		v[i] = value;
}

void Set_VectorInt(int *v, int n, int value) {
	unsigned int i;
	
	for(i=0;i<n;i++)
		v[i] = value;
}


////////////////////////////////////////////////////////////////////////////////
// File: transpose_matrix.c                                                   //
// Routine(s):                                                                //
//    Transpose_Matrix                                                        //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Transpose_Matrix(scalar *At, scalar *A, int nrows, int ncols)        //
//                                                                            //
//  Description:                                                              //
//     Take the transpose of A and store in At, i.e. At = A'.                 //
//     The matrix At should be declared as scalar At[ncols][nrows] in the     //
//     calling routine, and the matrix A declared as scalar A[nrows[ncols].   //
//     In general, At and A should be disjoint i.e. their memory locations    //
//     should be distinct.                                                    //
//                                                                            //
//  Arguments:                                                                //
//     scalar *At   Pointer to the first element of the matrix At.            //
//     scalar *A    Pointer to the first element of the matrix A.             //
//     int    nrows The number of rows of matrix A and number of columns of   //
//                  the matrix At.                                            //
//     int    ncols The number of columns of the matrix A and the number of   //
//                  rows of the matrix At.                                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     scalar A[M][N],  At[N][M];                                             //
//                                                                            //
//     (your code to initialize the matrix A)                                 //
//                                                                            //
//     Transpose_Matrix(&At[0][0], &A[0][0], M, N);                           //
//     printf("The transpose of A is the matrix At \n"); ...                  //
////////////////////////////////////////////////////////////////////////////////
void Transpose_Matrix(scalar *At, scalar *A, int nrows, int ncols) 
{
   scalar *pA;
   scalar *pAt;
   int i,j;

   for (i = 0; i < nrows; At += 1, A += ncols, i++) {
      pAt = At;
      pA = A;
      for (j = 0; j < ncols; pAt += nrows, j++) *pAt = pA[j];
   }
}

void Get2DMatrixFrom3D(scalar *matrix3D, int dim, int dim_index, int i_size, int j_size, int k_size, scalar *R) {
	int i,j,k;
	
	switch(dim){
		case(0):	
			for(j=0;j<j_size;j++)
				for(k=0;k<k_size;k++)
					R[pos2Dto1D(j, k, k_size)] = matrix3D[pos3Dto1D(dim_index, j, k, j_size, k_size)];
		break;
		
		case(1):	
			for(i=0;i<i_size;i++)
				for(k=0;k<k_size;k++)
					R[pos2Dto1D(i, k, k_size)] = matrix3D[pos3Dto1D(i, dim_index, k, j_size, k_size)];
		break;
		
		case(2):	
			for(i=0;i<i_size;i++)
				for(j=0;j<j_size;j++)
					R[pos2Dto1D(i, j, j_size)] = matrix3D[pos3Dto1D(i, j, dim_index, j_size, k_size)];
		break;
		
		default:
			printf("\nError: dim must be 0, 1 or 2 (i, j or k).\n");
		break;
	}
}

void GetVectorFrom2D(scalar *matrix2D, int dim, int dim_index, int i_size, int j_size, scalar *v) {
    int i,j;

    switch(dim){
        case(0):
            for(j=0;j<j_size;j++)
                v[j] = matrix2D[pos2Dto1D(dim_index, j, j_size)];
        break;

        case(1):
            for(i=0;i<i_size;i++)
                v[i] = matrix2D[pos2Dto1D(i, dim_index, j_size)];
        break;

        default:
            printf("\nError: dim must be 0 or 1 (i or j).\n");
        break;
    }
}

void GetVectorFrom3DMatrix(scalar *matrix3D, int dim1, int dim1_index, int dim2, int dim2_index, int i_size, int j_size, int k_size, scalar *v) {
	int i,j,k;
	
	switch(dim1+dim2){
		case(1):
			for(k=0;k<k_size;k++)
				v[k] = matrix3D[pos3Dto1D(dim1_index, dim2_index, k, j_size, k_size)];
		break;
		
		case(2):	
			for(j=0;j<j_size;j++)
				 v[j] = matrix3D[pos3Dto1D(dim1_index, j, dim2_index, j_size, k_size)];
		break;
		
		case(3):	
			for(i=0;i<i_size;i++)
				 v[i] = matrix3D[pos3Dto1D(i, dim1_index, dim2_index, j_size, k_size)];
		break;
		
		default:
			printf("\nError: dim1 and dim2 must be 0, 1 or 2 (i, j or k).\n");
		break;
	}
}

void Set2DMatrixTo3D(scalar *matrix3D, int dim, int dim_index, int i_size, int j_size, int k_size, scalar *R) {
	int i,j,k;
	
	switch(dim){
		case(0):	
			for(j=0;j<j_size;j++)
				for(k=0;k<k_size;k++)
					matrix3D[pos3Dto1D(dim_index, j, k, j_size, k_size)] = R[pos2Dto1D(j, k, k_size)];
		break;
		
		case(1):	
			for(i=0;i<i_size;i++)
				for(k=0;k<k_size;k++)
					matrix3D[pos3Dto1D(i, dim_index, k, j_size, k_size)] = R[pos2Dto1D(i, k, k_size)];
		break;
		
		case(2):	
			for(i=0;i<i_size;i++)
				for(j=0;j<j_size;j++)
					matrix3D[pos3Dto1D(i, j, dim_index, j_size, k_size)] = R[pos2Dto1D(i, j, j_size)];
		break;
		
		default:
			printf("\nError: dim must be 0, 1 or 2 (i, j or k).\n");
		break;
	}
}

void SetVectorTo3DMatrix(scalar *matrix3D, int dim1, int dim1_index, int dim2, int dim2_index, int i_size, int j_size, int k_size, scalar *v) {
	int i,j,k;
	
	switch(dim1+dim2){
		case(1):
			for(k=0;k<k_size;k++)
				matrix3D[pos3Dto1D(dim1_index, dim2_index, k, j_size, k_size)] = v[k];
		break;
		
		case(2):	
			for(j=0;j<j_size;j++)
				matrix3D[pos3Dto1D(dim1_index, j, dim2_index, j_size, k_size)] = v[j];
		break;
		
		case(3):	
			for(i=0;i<i_size;i++)
				matrix3D[pos3Dto1D(i, dim1_index, dim2_index, j_size, k_size)] = v[i];
		break;
		
		default:
			printf("\nError: dim1 and dim2 must be 0, 1 or 2 (i, j or k).\n");
		break;
	}
}

void SetVectorTo2DMatrix(scalar *matrix2D, int dim, int dim_index, int i_size, int j_size, scalar *v) {
    int i,j;

    switch(dim){
        case(0):
            for(j=0;j<j_size;j++)
                matrix2D[pos2Dto1D(dim_index, j, j_size)] = v[j];
        break;

        case(1):
            for(i=0;i<i_size;i++)
                matrix2D[pos2Dto1D(i, dim_index, j_size)] = v[i];
        break;

        default:
            printf("\nError: dim must be 0 or 1 (i or j).\n");
        break;
    }
}

void SetVectorTo2DMatrix_int(int *matrix2D, int dim, int dim_index, int i_size, int j_size, int *v) {
    int i,j;

    switch(dim){
        case(0):
            for(j=0;j<j_size;j++)
                matrix2D[pos2Dto1D(dim_index, j, j_size)] = v[j];
        break;

        case(1):
            for(i=0;i<i_size;i++)
                matrix2D[pos2Dto1D(i, dim_index, j_size)] = v[i];
        break;

        default:
            printf("\nError: dim must be 0 or 1 (i or j).\n");
        break;
    }
}

// http://www.mymathlib.com/matrices/arithmetic/mul_mat.html

// Ares(nrows1 x ncols2) = A1(nrows1 x ncols1) x A2(ncols1 x ncols2)
// Se asume que todas las dimensiones de las matrices son correctas
void Matrix_Multiplication(scalar *M1, int nrows1, int ncols1, scalar *M2, int ncols2, scalar *Mres){
	scalar aux = 0.0;
	int i,j,k;

	for(i=0;i<nrows1;i++) {
		for(j=0;j<ncols2;j++) {
			for(k=0;k<ncols1;k++) {
				aux += M1[pos2Dto1D(i,k,ncols1)] * M2[pos2Dto1D(k,j,ncols2)];
			}
			Mres[pos2Dto1D(i,j,ncols2)] = aux;
			aux = 0.0;
		}
	}

}

scalar norm2(scalar *v, int size) {
	unsigned int i;
	scalar norm = 0.0;

	for(i=0;i<size;i++)
		norm += pow((double)v[i],2.0);
	
	norm = sqrt(norm);

	return norm;
}

scalar scalarProduct(scalar *v1, scalar *v2, int size) {
	unsigned int i;
	scalar sum = 0.0;

	for(i=0;i<size;i++)
		sum += v1[i]*v2[i];

	return sum;
}

void AddScalar(scalar *A, scalar *R, int nrows, int ncols, scalar s){
	int i,j;

	for(i=0;i<nrows;i++)
		for(j=0;j<ncols;j++)
			R[pos2Dto1D(i,j,ncols)] = A[pos2Dto1D(i,j,ncols)] + s;
}

void TimesScalar(scalar *A, scalar *R, int nrows, int ncols, scalar s){
	int i,j;

	for(i=0;i<nrows;i++)
		for(j=0;j<ncols;j++)
			R[pos2Dto1D(i,j,ncols)] = A[pos2Dto1D(i,j,ncols)]*s;
}

void PowScalar(scalar *A, scalar *R, int nrows, int ncols, scalar s){
	int i,j;

	for(i=0;i<nrows;i++)
		for(j=0;j<ncols;j++)
			R[pos2Dto1D(i,j,ncols)] = pow(A[pos2Dto1D(i,j,ncols)], s);
}

void MaxScalar(scalar *A, scalar *R, int nrows, int ncols, scalar s){
	int i,j;

	for(i=0;i<nrows;i++)
		for(j=0;j<ncols;j++)
			R[pos2Dto1D(i,j,ncols)] = MAX(A[pos2Dto1D(i,j,ncols)], s);
}

void AddNbyN(scalar *A, scalar *B, scalar *R, int nrows, int ncols){
	int i,j;

	for(i=0;i<nrows;i++)
		for(j=0;j<ncols;j++)
			R[pos2Dto1D(i,j,ncols)] = A[pos2Dto1D(i,j,ncols)] + B[pos2Dto1D(i,j,ncols)];
}

void SubNbyN(scalar *A, scalar *B, scalar *R, int nrows, int ncols){
	int i,j;

	for(i=0;i<nrows;i++)
		for(j=0;j<ncols;j++)
			R[pos2Dto1D(i,j,ncols)] = A[pos2Dto1D(i,j,ncols)] - B[pos2Dto1D(i,j,ncols)];
}

void TimesNbyN(scalar *A, scalar *B, scalar *R, int nrows, int ncols){
	int i,j;

	for(i=0;i<nrows;i++)
		for(j=0;j<ncols;j++)
			R[pos2Dto1D(i,j,ncols)] = A[pos2Dto1D(i,j,ncols)]*B[pos2Dto1D(i,j,ncols)];
}

void RepMatCols(scalar *v, scalar *R, int nrows, int ncols){
	int i,j;

	for(i=0;i<nrows;i++)
		for(j=0;j<ncols;j++)
			R[pos2Dto1D(i,j,ncols)] = v[i];
}

// salvamento de correlaciones de la matriz simétrica A[size][size], resultado en Ares, en B la matriz que utilizamos en lugar de cholesky(Ares), OJO!! cholesky(Ares) puede dar error pq Ares NO es definida positiva, es SEMIdefinida  positiva!!
void correlationsRescue(scalar *A, int size, scalar *Ares, scalar *B) {
	
	scalar *eigenvalues, *eigenvectors;
	eigenvalues = (scalar*)malloc(size*sizeof(scalar));
	eigenvectors = (scalar*)malloc(size*size*sizeof(scalar));

	Jacobi_Cyclic_Method(eigenvalues, eigenvectors, A, size);
	Sort_Eigenvalues(eigenvalues, eigenvectors, size, -1);

	unsigned int i;
	for(i=0;i<size;i++)
		if (eigenvalues[i] <= 0)
			eigenvalues[i] = 0;


	scalar *eigenvaluesMatrix = (scalar*)malloc(size*size*sizeof(scalar));
	Set_Matrix(eigenvaluesMatrix, size, size, 0.0);
	Set_Diagonal(eigenvaluesMatrix, eigenvalues, size, size);

	for(i=0;i<size;i++)
		eigenvaluesMatrix[pos2Dto1D(i,i,size)] = sqrt(eigenvaluesMatrix[pos2Dto1D(i,i,size)]);

	//scalar *B = (scalar*)malloc(size*size*sizeof(scalar));
	Matrix_Multiplication(eigenvectors, size, size, eigenvaluesMatrix, size, B);

	unsigned int j;	
	scalar norma;
	for(i=0;i<size;i++) {
		norma = norm2(&(B[pos2Dto1D(i,0,size)]), size);
		for(j=0;j<size;j++) {
			B[pos2Dto1D(i,j,size)] = B[pos2Dto1D(i,j,size)]/norma;
		}
	}

	scalar *BTranspose = (scalar*)malloc(size*size*sizeof(scalar));
	Transpose_Matrix(BTranspose, B, size, size);

	Matrix_Multiplication(B, size, size, BTranspose, size, Ares);

	free(eigenvalues); free(eigenvectors); free(eigenvaluesMatrix); free(BTranspose);
}

///////////////////////////////////////////////////////////////////////////////
// Polynomial approximation of cumulative normal distribution function
///////////////////////////////////////////////////////////////////////////////
scalar cnd(scalar d){
    const double       N1 = 0.31938153;
    const double       N2 = -0.356563782;
    const double       N3 = 1.781477937;
    const double       N4 = -1.821255978;
    const double       N5 = 1.330274429;
    const double RSQRT2PI = 0.39894228040143267793994605993438;

    double
        K = 1.0 / (1.0 + 0.2316419 * (double)fabs(d));

    double
        cnd = RSQRT2PI * exp(- 0.5 * (double)d * (double)d) * 
        (K * (N1 + K * (N2 + K * (N3 + K * (N4 + K * N5)))));

    if(d > 0)
        cnd = 1.0 - cnd;

    return cnd;
}

void sumXcols(scalar *matrix, int nrows, int ncols, scalar *vRes){
	int i,j;
	
	for(j=0;j<ncols;j++) {
		vRes[j] = 0.0;
		for(i=0;i<nrows;i++){
			vRes[j] = vRes[j] + matrix[pos2Dto1D(i,j,ncols)];
		}
	}
}

scalar MinVector(scalar *vector, int n){
scalar aux = vector[0];
int i;
	for(i=1;i<n;i++)
		if(vector[i]<aux)
			aux = vector[i];
	return aux;
}

scalar MaxVector(scalar *vector, int n){
scalar aux = vector[0];
int i;
	for(i=1;i<n;i++)
		if(vector[i]>aux)
			aux = vector[i];
	return aux;
}

scalar MeanVector(scalar *vector, int n){

	int i;
	scalar sum;
	sum = 0.0;
	for(i=0;i<n;i++)
		sum += vector[i];

	return sum/n;
}

void printMatrix(scalar *matrix, int nrows, int ncols) {
	int i,j;

	for(i=0;i<nrows;i++) {
		for(j=0;j<ncols;j++){
			printf("%lf ", *(matrix++));
		}
		printf("\n");
	}
	printf("\n");
}

void printVector(scalar *vector, int n){
	int i;

	for(i=0;i<n;i++)
		printf("%lf ", vector[i]);
	printf("\n");
}

void printVectorInt(int *vector, int n){
    int i;

    for(i=0;i<n;i++)
        printf("%d ", vector[i]);
    printf("\n");
}

#endif /* mathLibrary.h  */
