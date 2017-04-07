#include <stdio.h>
#include <stdlib.h>

#include "mathLibrary.h"

#include "mygsllib.h"

#define eSIMS (64*MIL)
#define pSIMS (64*MIL)
#define STEPS 10
#define EXERS 10
#define ASSETS 5
#define BUNDS 4
// #define BUNDS (1*MIL)
#define BASISF 3
#define TRIALS 1

/*********************************************************
Global Variables
**********************************************************/
unsigned int SIMS;

/*********************************************************
Time Measuring functions
**********************************************************/
struct timeval t1, t2, Gt1, Gt2;

void cTic(){
	gettimeofday(&t1, NULL);
}

void cToc(char *c){
	gettimeofday(&t2, NULL);
	printf("Time %s: %f s\n", c, (t2.tv_usec - t1.tv_usec)/1000000.0f  + (t2.tv_sec - t1.tv_sec));
}

void GcTic(){
	gettimeofday(&Gt1, NULL);
}

void GcToc(char *c){
	gettimeofday(&Gt2, NULL);
	printf("Global Time %s: %f s\n", c, (Gt2.tv_usec - Gt1.tv_usec)/1000000.0f  + (Gt2.tv_sec - Gt1.tv_sec));
}

void compute_h(gsl_matrix *SData, gsl_vector *h){

	/***************Arithmetic Basket option**************/
	gsl_vector *ones = gsl_vector_alloc(ASSETS);
	gsl_vector_set_all(ones, 1.0);

	gsl_blas_dgemv(CblasTrans, 1.0/ASSETS, SData, ones, 0.0, h);

	gsl_vector_free(ones);
	/*****************************************************/
}//end compute_h

int Bundling(gsl_matrix **SPaths, unsigned int BundIdx[STEPS+1][SIMS]){

	unsigned int i, k;
	gsl_vector *CData = gsl_vector_alloc(SIMS);
	gsl_permutation *p = gsl_permutation_alloc(SIMS);
	for(k = 0; k<STEPS+1; k++){
		compute_h(SPaths[k], CData);
// 		printVector(CData->data, SIMS);
		gsl_sort_vector_index(p, CData);
// 		for(i = 0; i<SIMS; i++)
// 			printf("%ld ", gsl_permutation_get(p, i));
// 		printf("\n");
		for(i = 0; i<SIMS; i++)
			BundIdx[k][gsl_permutation_get(p, i)] = (unsigned int)i/(SIMS/BUNDS);
// 		for(i = 0; i<SIMS; i++)
// 			printf("%ld ", BundIdx[k][i]);
// 		printf("\n");
	}
	gsl_permutation_free(p);
	gsl_vector_free(CData);
	return 0;
}

void compute_expectation(gsl_matrix *SPaths, scalar r, gsl_vector *sigma, gsl_matrix *rho, scalar dt, gsl_matrix *EPaths){

	int i, j1, j2;

	scalar sum, sig11, sig22, sig12, sum_mu, siga;

	gsl_vector *SData =  gsl_vector_alloc(ASSETS);

	for(i = 0; i<SIMS; i++){
		gsl_matrix_get_col(SData, SPaths, i);

		// First Basis functions
		gsl_matrix_set(EPaths, 0, i, 1.0);

		// Second Basis functions
		sum = 0.0;
		for(j1 = 0; j1<ASSETS; j1++)
			sum += exp(r*dt)*gsl_vector_get(SData,j1)/ASSETS;
		gsl_matrix_set(EPaths, 1, i, sum);

		// Third Basis functions
		sum = 0.0;
		for(j1 = 0; j1<ASSETS; j1++){
			sig11 = gsl_vector_get(sigma, j1)*gsl_vector_get(sigma, j1);
			for(j2 = 0; j2<ASSETS; j2++){
				sig22 = gsl_vector_get(sigma, j2)*gsl_vector_get(sigma, j2);
				sig12 = gsl_vector_get(sigma, j1)*gsl_vector_get(sigma, j2);
				sum_mu = (r - 0.5*sig11)*dt + (r - 0.5*sig22)*dt;
				siga = (sig11 + sig22) + 2.0*gsl_matrix_get(rho, j1, j2)*sig12;
				sum += exp(sum_mu + 0.5*siga*dt)*(gsl_vector_get(SData, j1)/ASSETS)*(gsl_vector_get(SData, j2)/ASSETS);
			}//endfor j2
		}//endfor j1
		gsl_matrix_set(EPaths, 2, i, sum);
	}//endfor SIMS

}//end compute_expectation

int GBM_model(gsl_matrix *S, gsl_vector *sig, scalar mu, scalar dt, gsl_matrix *Z){

	int i, j;
	scalar S_, sig_, Z_;
	for(i = 0; i < S->size1; i++){
		sig_ = gsl_vector_get(sig, i);
		for(j = 0; j < S->size2; j++){
			S_ = gsl_matrix_get(S, i, j);
			Z_ = gsl_matrix_get(Z, i, j);
			gsl_matrix_set(Z, i, j, exp( log(S_) + ((mu - 0.5*sig_*sig_)*dt + sig_*sqrt(dt)*Z_) ));
		}
	}
	return 0;
}

/*********************************************************
Monte Carlo simulation for several assets
- simulations are stored in SPaths
**********************************************************/
void MonteCarlo_MultiAsset(scalar mu,  scalar T, gsl_vector *S0, gsl_vector *sigma, gsl_matrix *rho, gsl_matrix *C, gsl_matrix **SPaths, gsl_vector **ZPaths, gsl_matrix **EPaths){

	int i, j, j2, k;
	scalar dt;
	dt = T/STEPS;

	for(i = 0; i<SIMS; i++)
		gsl_matrix_set_col(SPaths[0], i, S0);
	compute_expectation(SPaths[0], mu, sigma, rho, dt, EPaths[0]);

	
	gsl_matrix *Z = gsl_matrix_alloc(ASSETS, SIMS);
	for(k = 1; k<STEPS+1; k++){

		randomGaussianMatrix(0.0, 1.0, ASSETS, SIMS, Z);
		gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, C, Z);
		
		GBM_model(SPaths[k-1], sigma, mu, dt, Z);
		gsl_matrix_memcpy(SPaths[k], Z);

		compute_h(SPaths[k], ZPaths[k-1]);

		compute_expectation(SPaths[k], mu, sigma, rho, dt, EPaths[k]);

	}//end for STEPS
	gsl_matrix_free(Z);


}//end MonteCarlo_MultiAsset

scalar stdeviation(scalar *v, int n, scalar mean){

	int i;
	scalar sum = 0.0;
	for(i=0;i<n;i++)
		sum += (v[i] - mean)*(v[i] - mean);

	return sqrt(sum/(n-1));

}//end std

/*********************************************************
Compute how many paths are in each bundle
- number of paths for each bundle is stored in pXb
**********************************************************/
void pathsXbundle(int *gIdx, int *pXb){

	int i;
	Set_VectorInt(pXb, BUNDS, 0);

	for(i=0;i<SIMS;i++)
		pXb[gIdx[i]] += 1;

}//end pathsXbundle

scalar directEstimator(gsl_matrix **SPaths, unsigned int BundIdx[][SIMS], gsl_vector **ZPaths, gsl_matrix **EPaths, scalar X, scalar r, scalar T, gsl_matrix **bundle_a){

	int st, i, b, bind, bsize;

	int pXb[BUNDS];

	scalar price;
	scalar dt = T/STEPS;

	gsl_vector *ZData = gsl_vector_alloc(SIMS);
	gsl_matrix *EData = gsl_matrix_alloc(BASISF, SIMS);
	gsl_matrix *RegrMat = gsl_matrix_alloc(SIMS, BASISF);
	gsl_permutation *p = gsl_permutation_alloc(SIMS);
	gsl_vector *tau = gsl_vector_alloc(BASISF);
	gsl_vector *vZ = gsl_vector_alloc(SIMS);
	gsl_vector *alpha = gsl_vector_alloc(BASISF);
	gsl_vector *residual = gsl_vector_alloc(SIMS/BUNDS);

	gsl_vector *ContValue = gsl_vector_calloc(SIMS);

	for(st = STEPS-1; st>-1; st--){
// 		printf("step: %d\n", st);

		/***************** Extract the data for time st ******************/
		gsl_vector_memcpy(ZData, ZPaths[st]);
		gsl_matrix_memcpy(EData, EPaths[st]);
// 		printVectorInt(BundIdx[st], SIMS);

		/***************** Sort according to the bundles ******************/
		/***************** Efficiency ******************/
		gsl_sort_int_index(p->data, BundIdx[st], 1, SIMS);
// 		for(i = 0; i<SIMS; i++)
// 			printf("%ld ", gsl_permutation_get(p, i));
// 		printf("\n");
// 		gsl_permutation_inverse(p_inv, p);
// 		for(i = 0; i<SIMS; i++)
// 			printf("%ld ", gsl_permutation_get(p_inv, i));
// 		printf("\n");

// 		printVector(ZData->data, SIMS);
// 		printMatrix(EData->data, BASISF, SIMS);
		gsl_permute_vector(p, ZData);
		gsl_permute_vector(p, ContValue);
		gsl_permute_matrix(p, EData);
// 		printVector(ZData->data, SIMS);
// 		printMatrix(EData->data, BASISF, SIMS);

		/***************** Regression matrix ******************/
		// Make it general (This is just a particular choice of basis functions)
		gsl_vector_set_all(vZ, 1.0);
		gsl_matrix_set_col(RegrMat, 0, vZ);
		for(b = 1; b<BASISF; b++){
			gsl_vector_mul(vZ, ZData);
			gsl_matrix_set_col(RegrMat, b, vZ);
		}

		/***************** Payoff ******************/
// 		gsl_vector_memcpy(vZ, ZData);
// 		gsl_vector_add_constant(ZData, - X);//Call option
		gsl_vector_scale(ZData, -1.0);
		gsl_vector_add_constant(ZData, X);//Put option
// 		printVector(ZData->data, SIMS);
		gsl_vector_max_value(ZData, 0.0);
		gsl_vector_max_elements(ZData, ContValue);
// 		printVector(ZData->data, SIMS);

		/***************** Regression ******************/
		gsl_matrix_view subRegrMat;
		gsl_vector_view subZData;
		gsl_matrix_view subEData;
		gsl_vector_view subContValue;

		pathsXbundle(BundIdx[st], pXb);
// 		printVectorInt(pXb, BUNDS);
		bind = 0;
		for(b = 0; b<BUNDS; b++){
			bsize = pXb[b];
// 			printf("bind: %ld, bsize: %ld\n", bind, bsize);
			subRegrMat = gsl_matrix_submatrix(RegrMat, bind, 0, bsize, BASISF);
			gsl_linalg_QR_decomp(&subRegrMat.matrix, tau);
			
			subZData = gsl_vector_subvector(ZData, bind, bsize);
			
			gsl_linalg_QR_lssolve(&subRegrMat.matrix, tau, &subZData.vector, alpha, residual);
			
// 			printVector(alpha->data, BASISF);
			gsl_matrix_set_row(bundle_a[st], b, alpha);

			subEData = gsl_matrix_submatrix(EData, 0, bind, BASISF, bsize);
			subContValue = gsl_vector_subvector(ContValue, bind, bsize);

			gsl_blas_dgemv(CblasTrans, exp(-r*dt), &subEData.matrix, alpha, 0.0, &subContValue.vector);
// 			printVector(ContValue->data, SIMS);
			
			bind += bsize;
		}

		gsl_permute_vector_inverse(p, ContValue);

	}//end for steps


	price = gsl_vector_sum(ContValue)/SIMS;
	
// 	printf("\n\nPrice: %lf\n", price);

	gsl_vector_free(ContValue);

	gsl_vector_free(ZData);
	gsl_matrix_free(EData);
	gsl_permutation_free(p);
	gsl_vector_free(tau);
	gsl_vector_free(vZ);
	gsl_vector_free(alpha);
	gsl_vector_free(residual);

	return price;

}//end directEstimator

int main(int argc, char **argv){

	//Input data
	scalar K, r, q, T, rho;

	K = 40;
	r = 0.06;
	q = 0.0;
	T = 1.0;
	rho = 0.25;

	gsl_vector *S0 = gsl_vector_alloc(ASSETS);
	gsl_vector *sigma = gsl_vector_alloc(ASSETS);
	gsl_matrix *Rho = gsl_matrix_alloc(ASSETS, ASSETS);
	gsl_matrix_set_identity(Rho);

	unsigned int i, j;
	for(i=0;i<ASSETS;i++){
		gsl_vector_set(S0, i, 40);
		gsl_vector_set(sigma, i, 0.2);
		for(j=0;j<ASSETS;j++)
			 if(i!=j)
				gsl_matrix_set(Rho, i, j, rho);
	}

	gsl_matrix *C = gsl_matrix_alloc(ASSETS, ASSETS);
	gsl_matrix_memcpy(C, Rho);
	gsl_linalg_cholesky_decomp(C);

	scalar dE[TRIALS];
	scalar pE[TRIALS];
	scalar sum_dE = 0.0;
	scalar sum_pE = 0.0;

	gsl_matrix *SPaths[STEPS+1];
	gsl_vector *ZPaths[STEPS];
	gsl_matrix *EPaths[STEPS+1];
	gsl_matrix *bundle_a[STEPS];

GcTic();

	gen = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(gen, 1234);
	for(i=0;i<TRIALS;i++){

		/**********************DIRECT ESTIMATOR************************/
		SIMS = eSIMS;
		printf("\nSIMS = %d | STEPS = %d | EXERS = %d | ASSETS = %d | BUNDS = %d | RATIO = %d\n", SIMS, STEPS, EXERS, ASSETS, BUNDS, SIMS/BUNDS);

		alloc_array_gsl_matrix(SPaths, STEPS+1, ASSETS, SIMS);
		alloc_array_gsl_vector(ZPaths, STEPS, SIMS);
		alloc_array_gsl_matrix(EPaths, STEPS+1, BASISF, SIMS);
		
		unsigned int BundIdx[STEPS+1][SIMS];
		
		alloc_array_gsl_matrix(bundle_a, STEPS, BUNDS, BASISF);
cTic();
		MonteCarlo_MultiAsset(r-q, T, S0, sigma, Rho, C, SPaths, ZPaths, EPaths);
cToc("MC DE");

cTic();
		Bundling(SPaths, BundIdx);
cToc("Bundling");

cTic();
		dE[i] = directEstimator(SPaths, BundIdx, ZPaths, EPaths, K, r, T, bundle_a);
		sum_dE += dE[i];
cToc("directEst");
		printf("\nDirect Estimator: %lf\n\n", dE[i]);

		free_array_gsl_matrix(SPaths, STEPS+1);
		free_array_gsl_vector(ZPaths, STEPS);
		free_array_gsl_matrix(EPaths, STEPS+1);

		free_array_gsl_matrix(bundle_a	, STEPS);


	}//end for trials

	gsl_rng_free(gen);
	
	gsl_vector_free(S0);
	gsl_vector_free(sigma);
	gsl_matrix_free(Rho);
	gsl_matrix_free(C);

GcToc("Final");

	printf("------------------MEAN-------------------");
	printf("\nDirect Estimator: %lf\n", sum_dE/TRIALS);
	printf("\nPath Estimator: %lf\n", sum_pE/TRIALS);

	printf("------------------STD-------------------");
	printf("\nDirect Estimator: %.10lf\n", stdeviation(dE, TRIALS, sum_dE/TRIALS));
	printf("\nPath Estimator: %.10lf\n", stdeviation(pE, TRIALS, sum_pE/TRIALS));

	printf("\n\n");

}//end main

