#include <stdio.h>
#include <stdlib.h>

#include "mathLibrary.h"

#include "mygsllib.h"

#define bSIMS (16*MIL)
#define eSIMS (64*MIL)
#define STEPS 10
#define EXERS 10
#define ASSETS 5
#define BUNDS 64
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

/*******************RANDOM NUMBERS********************/

/*scalar box_muller(scalar m, scalar s){
	scalar x1, x2, w, y1;
	static scalar y2;
	static int use_last = 0;

	if (use_last){
		y1 = y2;
		use_last = 0;
	}else{
		do{
			x1 = 2.0*gsl_rng_uniform(gen) - 1.0;
			x2 = 2.0*gsl_rng_uniform(gen) - 1.0;
			w = x1*x1 + x2*x2;
		}while ( w >= 1.0 );

		w = sqrt((-2.0*log(w))/w);
		y1 = x1*w;
		y2 = x2*w;
		use_last = 1;
	}

	return( m + y1*s );
}//end box_muller
*/

/*********************************************************
Auxiliary function of k-means
**********************************************************/
void mean(scalar *X, scalar *c, int *gIdx){

	int p,s;
	scalar sum[BUNDS*ASSETS];
	Set_Matrix(sum, BUNDS, ASSETS, 0.0);
	int div[BUNDS];
	for(p=0;p<BUNDS;p++)
		div[p] = 0;

	for(p=0;p<SIMS;p++){
		for(s=0;s<ASSETS;s++){
			sum[pos2Dto1D(gIdx[p],s,ASSETS)] += X[pos2Dto1D(p,s,ASSETS)];
		}
		div[gIdx[p]] += 1;
	}

	//printMatrix(sum, BUNDS, ASSETS);

	for(p=0;p<BUNDS;p++)
		if(div[p]!=0)
			for(s=0;s<ASSETS;s++)
				c[pos2Dto1D(p,s,ASSETS)] = sum[pos2Dto1D(p,s,ASSETS)]/div[p];

}//end means

/*void computeDifferences_l(scalar *X, scalar *c, int *gIdx){

	int p, t, s;
	scalar d, D, vmin;

	for(p=0;p<SIMS;p++){
		vmin = 99999.999;
		for(t=0;t<BUNDS;t++){
			D = 0.0;
			for(s=0;s<ASSETS;s++){
				d = X[pos2Dto1D(p, s, ASSETS)] - c[pos2Dto1D(t, s, ASSETS)];
				D += d*d;
			}//end for ASSETS
			if( D < vmin ){
				gIdx[p] = t;
				vmin = D;
			}//end if
		}//end for BUNDS
	}//end for SIMS
}//end computeDifferences

void k_means_l(scalar *X, int *gIdx, scalar *c){

	computeDifferences_l(SIMS, X, c, gIdx);

	mean(SIMS, X, c, gIdx);

}//end k_means*/

int computeDifferences(scalar *X, scalar *c, int *gIdx){

	int p, t, s, stop;
	scalar d, D, g0, vmin;

	stop = 1;
	for(p=0;p<SIMS;p++){
		g0 = gIdx[p];
		vmin = 99999.999;
		for(t=0;t<BUNDS;t++){
			D = 0.0;
			for(s=0;s<ASSETS;s++){
				d = X[pos2Dto1D(p, s, ASSETS)] - c[pos2Dto1D(t, s, ASSETS)];
				D += d*d;
			}//end for ASSETS
			if( D < vmin ){
				gIdx[p] = t;
				vmin = D;
			}//end if
		}//end for BUNDS
		if(g0 != gIdx[p])
			stop = 0;
	}//end for SIMS

	return stop;

}//end computeDifferences

/*********************************************************
k-means algorithm for bundling
- indexes of groups are stored in gIdx
- centroids of groups are stored in c
**********************************************************/
void k_means(scalar *X, int *gIdx, scalar *c){

	//int t, s;
	int stop, c_ind;

	/*for(t=0;t<BUNDS;t++){
		c_ind = (int)(RAND_UNIFORM*SIMS);
		for(s=0;s<ASSETS;s++)
			c[pos2Dto1D(t,s,ASSETS)] = X[pos2Dto1D(c_ind,s,ASSETS)];
	}*/

	do{

		stop = computeDifferences(X, c, gIdx);

		mean(X, c, gIdx);

	}while(!stop);//end while

}//end k_means

/*********************************************************
Bundling SPaths step by step
- indexes of groups are stored in gIdx
- centroids of groups are stored in ctrs
**********************************************************/
/*void Bundling(scalar *SPaths, int *gIdx, scalar *ctrs){

	int st, b, s, c_ind;
	scalar BunS[SIMS*ASSETS];
	scalar c[BUNDS*ASSETS];

	for(b=0;b<BUNDS;b++){
		c_ind = (int)(RAND_UNIFORM*SIMS);
		for(s=0;s<ASSETS;s++)
			c[pos2Dto1D(b,s,ASSETS)] = SPaths[pos3Dto1D(c_ind, s, 1, ASSETS, STEPS+1)];
	}

	for(st=1;st<STEPS+1;st++){
		Get2DMatrixFrom3D(SPaths, 2, st, SIMS, ASSETS, STEPS+1, BunS);
		//printMatrix(BunS, SIMS, ASSETS);
		k_means(BunS, gIdx, c);
		Set2DMatrixTo3D(ctrs, 2, st-1, BUNDS, ASSETS, STEPS, c);
		//printMatrix(ctrs, BUNDS, ASSETS);
	}//end for STEPS+1

}//end Bundling*/

void reBundling(scalar *SPaths, int *gIdx, scalar *ctrs, int st){

	scalar BunS[SIMS*ASSETS];
	scalar c[BUNDS*ASSETS];

	Get2DMatrixFrom3D(SPaths, 2, st, SIMS, ASSETS, STEPS+1, BunS);
	Get2DMatrixFrom3D(ctrs, 2, st-1, BUNDS, ASSETS, STEPS, c);
	computeDifferences(BunS, c, gIdx);

}//end Bundling

void compute_h_payoff(scalar SData, int path, int step, int asset, scalar *h){

	/***************Arithmetic Basket option**************/
	if(asset == 0)
		h[pos2Dto1D(path, step, STEPS+1)] = 0.0;

	h[pos2Dto1D(path, step, STEPS+1)] += SData/ASSETS;
	/*****************************************************/
}//end compute_h_payoff

void compute_h_gsl(gsl_matrix *SData, gsl_vector *h){

	/***************Arithmetic Basket option**************/
	gsl_vector *ones = gsl_vector_alloc(ASSETS);
	gsl_vector_set_all(ones, 1.0);

	gsl_blas_dgemv(CblasTrans, 1.0/ASSETS, SData, ones, 0.0, h);

	gsl_vector_free(ones);
	/*****************************************************/
}//end compute_h_gsl

int Bundling(gsl_matrix **SPaths, unsigned int BundIdx[STEPS+1][SIMS]){

	unsigned int i, k;
	gsl_vector *CData = gsl_vector_alloc(SIMS);
	gsl_permutation *p = gsl_permutation_alloc(SIMS);
	for(k = 0; k<STEPS+1; k++){
		compute_h_gsl(SPaths[k], CData);
// 		printVector(CData->data, SIMS);
		gsl_sort_vector_index(p, CData);
		for(i = 0; i<SIMS; i++)
			printf("%ld ", gsl_permutation_get(p, i));
		printf("\n");
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

void compute_expectation(scalar *SData, scalar r, scalar *sigma, scalar *rho, scalar dt, scalar *EPaths, int path, int step){

	scalar sum_mu, siga;
	int i1, i2, i3, i4;

	EPaths[pos3Dto1D(path, 0, step, BASISF, STEPS+1)] = 1.0;

	EPaths[pos3Dto1D(path, 1, step, BASISF, STEPS+1)] = 0.0;
	for(i1=0; i1<ASSETS; i1++)
		EPaths[pos3Dto1D(path, 1, step, BASISF, STEPS+1)] += exp(r*dt)*SData[i1]/ASSETS;

	EPaths[pos3Dto1D(path, 2, step, BASISF, STEPS+1)] = 0.0;
	for(i1=0; i1<ASSETS; i1++){
		for(i2=0; i2<ASSETS; i2++){
			sum_mu = (r - 0.5*sigma[i1]*sigma[i1])*dt + (r - 0.5*sigma[i2]*sigma[i2])*dt;
			siga = (sigma[i1]*sigma[i1] + sigma[i2]*sigma[i2]) + 2.0*rho[pos2Dto1D(i1,i2,ASSETS)]*sigma[i1]*sigma[i2];
			EPaths[pos3Dto1D(path, 2, step, BASISF, STEPS+1)] += exp(sum_mu + 0.5*siga*dt)*(SData[i1]/ASSETS)*(SData[i2]/ASSETS);
		}
	}

// 	EPaths[pos3Dto1D(path, 3, step, BASISF, STEPS+1)] = 0.0;
// 	for(i1=0; i1<ASSETS; i1++){
// 		for(i2=0; i2<ASSETS; i2++){
// 			for(i3=0; i3<ASSETS; i3++){
// 				sum_mu = (r - 0.5*sigma[i1]*sigma[i1])*dt + (r - 0.5*sigma[i2]*sigma[i2])*dt + (r - 0.5*sigma[i3]*sigma[i3])*dt;
// 				siga = (sigma[i1]*sigma[i1] + sigma[i2]*sigma[i2] + sigma[i3]*sigma[i3]) + 2.0*(rho[pos2Dto1D(i1,i2,ASSETS)]*sigma[i1]*sigma[i2] + rho[pos2Dto1D(i1,i3,ASSETS)]*sigma[i1]*sigma[i3] + rho[pos2Dto1D(i2,i3,ASSETS)]*sigma[i2]*sigma[i3]);
// 				EPaths[pos3Dto1D(path, 3, step, BASISF, STEPS+1)] += exp(sum_mu + 0.5*siga*dt)*(SData[i1]/ASSETS)*(SData[i2]/ASSETS)*(SData[i3]/ASSETS);
// 			}
// 		}
// 	}
// 
// 	EPaths[pos3Dto1D(path, 4, step, BASISF, STEPS+1)] = 0.0;
// 	for(i1=0; i1<ASSETS; i1++){
// 		for(i2=0; i2<ASSETS; i2++){
// 			for(i3=0; i3<ASSETS; i3++){
// 				for(i4=0; i4<ASSETS; i4++){
// 					sum_mu = (r - 0.5*sigma[i1]*sigma[i1])*dt + (r - 0.5*sigma[i2]*sigma[i2])*dt + (r - 0.5*sigma[i3]*sigma[i3])*dt + (r - 0.5*sigma[i4]*sigma[i4])*dt;
// 					siga = (sigma[i1]*sigma[i1] + sigma[i2]*sigma[i2] + sigma[i3]*sigma[i3] + sigma[i4]*sigma[i4]) + 2.0*(rho[pos2Dto1D(i1,i2,ASSETS)]*sigma[i1]*sigma[i2] + rho[pos2Dto1D(i1,i3,ASSETS)]*sigma[i1]*sigma[i3] + rho[pos2Dto1D(i1,i4,ASSETS)]*sigma[i1]*sigma[i4] + rho[pos2Dto1D(i2,i3,ASSETS)]*sigma[i2]*sigma[i3] + rho[pos2Dto1D(i2,i4,ASSETS)]*sigma[i2]*sigma[i4] + rho[pos2Dto1D(i3,i4,ASSETS)]*sigma[i3]*sigma[i4]);
// 					EPaths[pos3Dto1D(path, 4, step, BASISF, STEPS+1)] += exp(sum_mu + 0.5*siga*dt)*(SData[i1]/ASSETS)*(SData[i2]/ASSETS)*(SData[i3]/ASSETS)*(SData[i4]/ASSETS);
// 				}
// 			}
// 		}
// 	}

}//end compute_expectation

void compute_expectation_gsl(gsl_matrix *SPaths, scalar r, gsl_vector *sigma, gsl_matrix *rho, scalar dt, gsl_matrix *EPaths){

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
// 		EPaths[pos3Dto1D(path, 1, step, BASISF, STEPS+1)] += exp(r*dt)*SData[i1]/ASSETS;

		// Third Basis functions
		sum = 0.0;
		for(j1 = 0; j1<ASSETS; j1++){
			sig11 = gsl_vector_get(sigma, j1)*gsl_vector_get(sigma, j1);
			for(j2 = 0; j2<ASSETS; j2++){
				sig22 = gsl_vector_get(sigma, j2)*gsl_vector_get(sigma, j2);
				sig12 = gsl_vector_get(sigma, j1)*gsl_vector_get(sigma, j2);
// 				sum_mu = (r - 0.5*sigma[i1]*sigma[i1])*dt + (r - 0.5*sigma[i2]*sigma[i2])*dt;
				sum_mu = (r - 0.5*sig11)*dt + (r - 0.5*sig22)*dt;
// 				siga = (sigma[i1]*sigma[i1] + sigma[i2]*sigma[i2]) + 2.0*rho[pos2Dto1D(i1,i2,ASSETS)]*sigma[i1]*sigma[i2];
				siga = (sig11 + sig22) + 2.0*gsl_matrix_get(rho, j1, j2)*sig12;
// 				EPaths[pos3Dto1D(path, 2, step, BASISF, STEPS+1)] += exp(sum_mu + 0.5*siga*dt)*(SData[i1]/ASSETS)*(SData[i2]/ASSETS);
				sum += exp(sum_mu + 0.5*siga*dt)*(gsl_vector_get(SData, j1)/ASSETS)*(gsl_vector_get(SData, j2)/ASSETS);
			}//endfor j2
		}//endfor j1
		gsl_matrix_set(EPaths, 2, i, sum);
	}//endfor SIMS

}//end compute_expectation_gsl

scalar GBM_model(scalar S, scalar sig, scalar mu, scalar dt, scalar dZ){

	return exp( log(S) + ((mu - 0.5*sig*sig)*dt + sig*sqrt(dt)*dZ) );
}

int GBM_model_gsl(gsl_matrix *S, gsl_vector *sig, scalar mu, scalar dt, gsl_matrix *Z){

	int i, j;
	scalar S_, sig_, Z_;
	for(i = 0; i < S->size1; i++){
		sig_ = gsl_vector_get(sig, i);
		for(j = 0; j < S->size2; j++){
			S_ = gsl_matrix_get(S, i, j);
// 			sig_ = gsl_matrix_get(sig, i, j);
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
void MonteCarlo_MultiAsset(scalar mu,  scalar T, /*scalar*/gsl_vector *S0, /*scalar*/gsl_vector *sigma, /*scalar*/gsl_matrix *rho, /*scalar*/gsl_matrix *C, scalar *SPaths, scalar *ZPaths, scalar *EPaths, gsl_matrix **SPaths_g, gsl_vector **ZPaths_g, gsl_matrix **EPaths_g){

	int i, j, j2, k;
	scalar dt;
	dt = T/STEPS;

	gsl_vector *dW = gsl_vector_alloc(ASSETS);
	scalar SData[ASSETS];
	scalar aux_SData;

	gsl_matrix *vZ[STEPS];
	alloc_array_gsl_matrix(vZ, STEPS, ASSETS, SIMS);

// 	gsl_matrix *SPaths_g[STEPS+1];
// 	alloc_array_gsl_matrix(SPaths_g, STEPS+1, ASSETS, SIMS);
	for(i = 0; i<SIMS; i++)
		gsl_matrix_set_col(SPaths_g[0], i, S0);

// 	gsl_vector *ZPaths_g[STEPS];
// 	alloc_array_gsl_vector(ZPaths_g, STEPS, SIMS);
	
// 	gsl_matrix *EPaths_g[STEPS+1];
// 	alloc_array_gsl_matrix(EPaths_g, STEPS+1, BASISF, SIMS);
	compute_expectation_gsl(SPaths_g[0], mu, sigma, rho, dt, EPaths_g[0]);

	for(k = 0; k<STEPS; k++){
		randomGaussianMatrix(0.0, 1.0, ASSETS, SIMS, vZ[k]);
		gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, C, vZ[k]);
	}
	
	gsl_matrix *Z = gsl_matrix_alloc(ASSETS, SIMS);
	for(k = 1; k<STEPS+1; k++){

// 		randomGaussianMatrix(0.0, 1.0, ASSETS, SIMS, Z);
// 		gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, C, Z);
		gsl_matrix_memcpy(Z, vZ[k-1]);
		GBM_model_gsl(SPaths_g[k-1], sigma, mu, dt, Z);
		gsl_matrix_memcpy(SPaths_g[k], Z);

		compute_h_gsl(SPaths_g[k], ZPaths_g[k-1]);

		compute_expectation_gsl(SPaths_g[k], mu, sigma, rho, dt, EPaths_g[k]);

	}//end for STEPS
	gsl_matrix_free(Z);
	
	printf("\nReference value GSL: %lf.\n\n", gsl_matrix_get(SPaths_g[STEPS], 0, 0));
	printf("\nReference value GSL: %lf.\n\n", gsl_vector_get(ZPaths_g[STEPS-1], 0));
	printf("\nReference value GSL: %lf.\n\n", gsl_matrix_get(EPaths_g[STEPS], 1, 0));

	
	for(k = 0; k<STEPS+1; k++){
		for(i = 0; i<SIMS; i++){
// 			for(j = 0; j<ASSETS; j++){
// 				dW[j] = box_muller(0.0, 1.0);
// 				gsl_vector_set(dW, j, gsl_ran_gaussian_ziggurat(gen, 1.0));
// 			}//end for ASSETS

// 			gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, C, dW);

			if(k!=0)
				gsl_matrix_get_col(dW, vZ[k-1], i);

			ZPaths[pos2Dto1D(i,k,STEPS+1)] = 1.0;
			for(j = 0; j<ASSETS; j++){
// 				dZ = 0.0;

// 				for(j2 = 0; j2<ASSETS; j2++)
// 					dZ += dW[j2]*C[pos2Dto1D(j, j2, ASSETS)];
// 				for(j2 = 0; j2<j+1; j2++)
// 					dZ += gsl_vector_get(dW, j2)*gsl_matrix_get(C, j, j2);

				if(k==0){
// 					aux_SData = S0[j];
					aux_SData = gsl_vector_get(S0, j);
					SPaths[pos3Dto1D(i, j, k, ASSETS, STEPS+1)] = aux_SData;
				}else{
// 					aux_SData = GBM_model(SPaths[pos3Dto1D(i, j, k-1, ASSETS, STEPS+1)], sigma[j], mu, dt, dZ);
					aux_SData = GBM_model(SPaths[pos3Dto1D(i, j, k-1, ASSETS, STEPS+1)], gsl_vector_get(sigma, j), mu, dt, gsl_vector_get(dW, j));
					SPaths[pos3Dto1D(i, j, k, ASSETS, STEPS+1)] = aux_SData;
				}

				compute_h_payoff(aux_SData, i, k, j, ZPaths);
			}//end for ASSETS
			GetVectorFrom3DMatrix(SPaths, 0, i, 2, k, SIMS, ASSETS, STEPS+1, SData);
			compute_expectation(SData, mu, sigma->data, rho->data, dt, EPaths, i, k);
		}//end for SIMS

	}//end for STEPS
	
	
// 	scalar error = 0.0;
// 	for(i = 0; i<SIMS; i++)
// 		for(j = 0; j<ASSETS; j++)
// 			for(k = 0; k<STEPS+1; k++)
// 				error += (SPaths[pos3Dto1D(i, j, k, ASSETS, STEPS+1)] - gsl_matrix_get(SPaths_g[k], j, i))*(SPaths[pos3Dto1D(i, j, k, ASSETS, STEPS+1)] - gsl_matrix_get(SPaths_g[k], j, i));
// 	printf("\nError = %lf\n\n", error);
	
	free_array_gsl_matrix(vZ, STEPS);
// 	free_array_gsl_matrix(SPaths_g, STEPS+1);
// 	free_array_gsl_vector(ZPaths_g, STEPS);
// 	free_array_gsl_matrix(EPaths_g, STEPS+1);

}//end MonteCarlo_MultiAsset

/*********************************************************
Compute ZPaths - Geometric mean by assets
- the mean is stored in ZPaths
**********************************************************/
// void computeZPaths(scalar *SPaths, scalar *ZPaths){
// 
// 	int i, j, k;
// 
// 	Set_Matrix(ZPaths, SIMS, STEPS+1, 1.0);
// 
// 	for(i=0;i<SIMS;i++)
// 		for(j=0;j<ASSETS;j++)
// 			for(k=0;k<STEPS+1;k++)
// 				ZPaths[pos2Dto1D(i,k,STEPS+1)] *= pow(SPaths[pos3Dto1D(i, j, k, ASSETS, STEPS+1)], 1.0/ASSETS);
// 
// }//end computeZPaths

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

void sortBundles(int *gIdx, int *old_Idx, scalar *ZData, scalar *EData){

	int aux_old_Idx[2*SIMS];
	scalar aux_ZData[SIMS];
	scalar aux_EData[SIMS*BASISF];

	gsl_sort_int_index(aux_old_Idx, gIdx, 1, SIMS);

	int i,j;
	for(i=0;i<SIMS;i++){
		old_Idx[i] = aux_old_Idx[2*i];
		aux_ZData[i] = ZData[old_Idx[i]];
		for(j=0;j<BASISF;j++)
			aux_EData[pos2Dto1D(i,j,BASISF)] = EData[pos2Dto1D(old_Idx[i],j,BASISF)];
	}

	for(i=0;i<SIMS;i++){
		ZData[i] = aux_ZData[i];
		for(j=0;j<BASISF;j++)
			EData[pos2Dto1D(i,j,BASISF)] = aux_EData[pos2Dto1D(i,j,BASISF)];
// 		ZDataCV[i] = aux_ZDataCV[i];
	}

}//end sortBundles

void computeOV(scalar *ZData, int *old_Idx, int sb, int bsize, int st, scalar *RegrMat, scalar X, scalar *OptionValue, scalar *ContinuationValue){

	scalar IntrinsicValue;
	int p, s;

	for(p=0;p<bsize;p++){

		for(s=0;s<BASISF;s++)
			RegrMat[pos2Dto1D(p, s, BASISF)] = pow(ZData[p], s);

		if( (st+1)%(STEPS/EXERS) == 0 ){
			//printf("\n---- %lf ----\n", ContinuationValue[p]);
			IntrinsicValue = PAYOFF_PUT(ZData[p], X);
			OptionValue[p] = MAX(IntrinsicValue, ContinuationValue[old_Idx[sb+p]]);
		}else{
			//printf("\nit mustn't entry\n");
			OptionValue[p] = ContinuationValue[old_Idx[sb+p]];
		}

	}// end for bsize

}//end computeOV

/*********************************************************
Solve a overdetermined system by means of least squared QR solver
- The solution is stored in res
**********************************************************/
/*void gsl_least_squared_solver(scalar *A_data, scalar *b_data, int rows, int cols, scalar *res){

	gsl_matrix_view A = gsl_matrix_view_array (A_data, rows, cols);
	gsl_vector_view b = gsl_vector_view_array (b_data, rows);

	gsl_vector *tau = gsl_vector_alloc(cols);

	gsl_linalg_QR_decomp(&A.matrix, tau);

	gsl_vector *x = gsl_vector_alloc(cols);
	gsl_vector *residual = gsl_vector_alloc(rows);

	gsl_linalg_QR_lssolve(&A.matrix, tau, &b.vector, x, residual);
	
	int i;
	for(i=0;i<cols;i++)
		res[i] = gsl_vector_get(x, i);

}//end least_squared_solver

void my_least_squared_solver(scalar *A_data, scalar *b_data, int rows, int cols, scalar *res){

	int i;

	leastsq(A_data, rows, cols, b_data);

	for(i=0;i<cols;i++)
		res[i] = b_data[i];

}//end my_least_squared_solver*/

void computeOptionPrice(scalar *EData, int *old_Idx, int sb, int bsize, scalar *a1, scalar r, scalar dt, scalar *ContinuationValue){

	scalar sum_RegrMatCV;
	int p, s;
    
	for(p=0;p<bsize;p++){

		sum_RegrMatCV = 0.0;
		for(s=0;s<BASISF;s++)
			sum_RegrMatCV += a1[s]*EData[pos2Dto1D(p,s,BASISF)];//exp((s*sum_mu + 0.5*s*s*siga)*dt)*pow(ZDataCV[p], s);

		ContinuationValue[old_Idx[sb+p]] = sum_RegrMatCV*exp(-r*dt);

	}// end for bsize
}//end computeOptionPrice

void computeBundles(scalar *ZData, scalar *EData, int *old_Idx, int *pXb, int st, scalar X, scalar r, scalar dt, scalar *mat_a1, scalar *ContinuationValue){

	int b, bb, sum_bb;
	scalar RegrMat[SIMS*BASISF];
	scalar OptionValue[SIMS];

	for(b=0;b<BUNDS;b++){
		sum_bb = 0;
		for(bb=0;bb<b;bb++)
			sum_bb += pXb[bb];

		computeOV(ZData+sum_bb, old_Idx, sum_bb, pXb[b], st, RegrMat+sum_bb*BASISF, X, OptionValue+sum_bb, ContinuationValue);

		leastsq(RegrMat+sum_bb*BASISF, pXb[b], BASISF, OptionValue+sum_bb);

		SetVectorTo2DMatrix(mat_a1, 1, b, BASISF, BUNDS, OptionValue+sum_bb);

		computeOptionPrice(EData+sum_bb*BASISF, old_Idx, sum_bb, pXb[b], OptionValue+sum_bb, r, dt, ContinuationValue);

	}//end for BUNDS

}//end computeBundles

/*********************************************************
Compute de Option value by means of Direct estimator
- 
**********************************************************/
scalar directEstimator(scalar *SPaths, scalar *ctrs, int *gIdx, scalar *ZPaths, scalar *EPaths, scalar X, scalar r, scalar T, scalar *bundle_a1){

	int st, i;

	scalar *ContinuationValue;
	int pXb[BUNDS];
	scalar ZData[SIMS];
	scalar EData[SIMS*BASISF];
	int old_Idx[SIMS];
	scalar mat_a1[BASISF*BUNDS];
	scalar price;
	scalar dt = T/STEPS;

	ContinuationValue = (scalar *)malloc(SIMS*sizeof(scalar));
	Set_Vector(ContinuationValue, SIMS, 0.0);

	for(st = STEPS-1; st>=0; st--){

		if(st>0){
			reBundling(SPaths, gIdx, ctrs, st);
			pathsXbundle(gIdx, pXb);
		}else{
			Set_VectorInt(pXb, BUNDS, 0);
			for(i=0;i<SIMS;i++){
				gIdx[i] = i/(SIMS/BUNDS);
				pXb[i/(SIMS/BUNDS)] += 1;
			}
		}

		GetVectorFrom2D(ZPaths, 1, st+1, SIMS, STEPS+1, ZData);
		Get2DMatrixFrom3D(EPaths, 2, st, SIMS, BASISF, STEPS+1, EData);

		sortBundles(gIdx, old_Idx, ZData, EData);

		computeBundles(ZData, EData, old_Idx, pXb, st, X, r, dt, mat_a1, ContinuationValue);

		Set2DMatrixTo3D(bundle_a1, 2, st, BASISF, BUNDS, STEPS, mat_a1);

	}//end for steps
	//printf("\nCV = %lf\n", ContinuationValue[0]);
	price = MeanVector(ContinuationValue, SIMS);
	free(ContinuationValue);

	return price;

}//end directEstimator

scalar directEstimator_gsl(gsl_matrix **SPaths, unsigned int BundIdx[][SIMS], gsl_vector **ZPaths, gsl_matrix **EPaths, scalar X, scalar r, scalar T, gsl_matrix **bundle_a){

	unsigned int st, i, b, bind, bsize;

// 	scalar *ContinuationValue;
	int pXb[BUNDS];
// 	scalar ZData[SIMS];
// 	scalar EData[SIMS*BASISF];
// 	int old_Idx[SIMS];
// 	scalar mat_a1[BASISF*BUNDS];
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
// 	ContinuationValue = (scalar *)malloc(SIMS*sizeof(scalar));
// 	Set_Vector(ContinuationValue, SIMS, 0.0);

	for(st = STEPS-1; st>0; st--){
		printf("step: %ld\n", st);
// 		GetVectorFrom2D(ZPaths, 1, st+1, SIMS, STEPS+1, ZData);
// 		ZData = ZPaths[st];
		gsl_vector_memcpy(ZData, ZPaths[st]);
		
// 		Get2DMatrixFrom3D(EPaths, 2, st, SIMS, BASISF, STEPS+1, EData);
// 		EData = EPaths[st];
		gsl_matrix_memcpy(EData, EPaths[st]);
// 		printVectorInt(BundIdx[st], SIMS);
		gsl_sort_int_index(p->data, BundIdx[st], 1, SIMS);
// 		for(i = 0; i<SIMS; i++)
// 			printf("%ld ", gsl_permutation_get(p, i));
// 		printf("\n");

// 		printVector(ZData->data, SIMS);
// 		printMatrix(EData->data, BASISF, SIMS);
		gsl_permute_vector(p, ZData);
		gsl_permute_matrix(p, EData);
// 		printVector(ZData->data, SIMS);
// 		printMatrix(EData->data, BASISF, SIMS);

		// Make it general (This is just a particular choice of basis functions)
		gsl_vector_set_all(vZ, 1.0);
		gsl_matrix_set_col(RegrMat, 0, vZ);
		for(b = 1; b<BASISF; b++){
			gsl_vector_mul(vZ, ZData);
			gsl_matrix_set_col(RegrMat, b, vZ);
		}

// 		gsl_vector_memcpy(vZ, ZData);
// 		gsl_vector_add_constant(ZData, - X);//Call option
		gsl_vector_scale(ZData, -1.0);
		gsl_vector_add_constant(ZData, X);//Put option
// 		printVector(ZData->data, SIMS);
		gsl_vector_max_value(ZData, 0.0);
		gsl_vector_max_elements(ZData, ContValue);
// 		printVector(ZData->data, SIMS);

		pathsXbundle(BundIdx[st], pXb);
// 		printVectorInt(pXb, BUNDS);
		
		gsl_matrix_view subRegrMat;
		gsl_vector_view subZData;
		gsl_matrix_view subEData;
		gsl_vector_view subContValue;
		bind = 0;
		for(b = 0; b<BUNDS; b++){
			bsize = pXb[b];
// 			printf("bind: %ld, bsize: %ld\n", bind, bsize);
			subRegrMat = gsl_matrix_submatrix(RegrMat, bind, 0, bsize, BASISF);
			gsl_linalg_QR_decomp(&subRegrMat.matrix, tau);
			
			subZData = gsl_vector_subvector(ZData, bind, bsize);
			
			gsl_linalg_QR_lssolve(&subRegrMat.matrix, tau, &subZData.vector, alpha, residual);
			
			printVector(alpha->data, BASISF);
			gsl_matrix_set_row(bundle_a[st], b, alpha);

			subEData = gsl_matrix_submatrix(EData, 0, bind, BASISF, bsize);
			subContValue = gsl_vector_subvector(ContValue, bind, bsize);

			gsl_blas_dgemv(CblasTrans, exp(-r*dt), &subEData.matrix, alpha, 0.0, &subContValue.vector);
			printVector(ContValue->data, SIMS);
			
			bind += bsize;
		}
		
		

// 		sortBundles(gIdx, old_Idx, ZData, EData);

// 		computeBundles(ZData, EData, old_Idx, pXb, st, X, r, dt, mat_a1, ContinuationValue);

// 		Set2DMatrixTo3D(bundle_a1, 2, st, BASISF, BUNDS, STEPS, mat_a1);

	}//end for steps
	//printf("\nCV = %lf\n", ContinuationValue[0]);
// 	price = MeanVector(ContinuationValue, SIMS);
// 	free(ContinuationValue);

	gsl_vector_max_value(ContValue, 0.0);
	price = gsl_vector_sum(ContValue)/SIMS;
	
	printf("\n\nPrice: %lf\n", price);

	gsl_vector_free(ZData);
	gsl_matrix_free(EData);
	gsl_permutation_free(p);
	gsl_vector_free(tau);
	gsl_vector_free(vZ);
	gsl_vector_free(alpha);
	gsl_vector_free(residual);

	gsl_vector_free(ContValue);

	return price;

}//end directEstimator_gsl

void computeCV(scalar *ZPaths, scalar *EPaths, scalar *CashFlows, int *ExerciseTime, scalar *bundle_a1, int *mat_gIdx, scalar X, scalar r, scalar dt, scalar *ContinuationValue){

	scalar ZData, ZDataCV, IntrinsicValue, sum_RegrMatCV;
	int p, st, s;

	for(p=0;p<SIMS;p++){

		for(st = STEPS-1; st>=0; st--){

			if( (st+1)%(STEPS/EXERS) == 0 ){
				ZData = ZPaths[pos2Dto1D(p, st+1, STEPS+1)];

				IntrinsicValue = PAYOFF_PUT(ZData, X);
				if(IntrinsicValue > ContinuationValue[p]){
					CashFlows[p] = IntrinsicValue;
					ExerciseTime[p] = st+1;
				}
			}

			if(st != 0){

				sum_RegrMatCV = 0.0;
				for(s=0;s<BASISF;s++)
					sum_RegrMatCV += bundle_a1[pos3Dto1D(s, mat_gIdx[pos2Dto1D(p, st-1, STEPS-1)], st, BUNDS, STEPS)]*EPaths[pos3Dto1D(p,s,st,BASISF,STEPS+1)];

				ContinuationValue[p] = sum_RegrMatCV*exp(-r*dt);
			}

		}//end for STEPS

	}// end for SIMS

}//end computeCV

scalar discountCashFlows(scalar *CashFlows, int *ExerciseTime, scalar r, scalar dt){

	scalar sum_CashF = 0.0;
	int p;
	for(p=0;p<SIMS;p++)
		sum_CashF += CashFlows[p]*exp(-r*dt*ExerciseTime[p]);
	return sum_CashF/SIMS;
}//end discountCashFlows

/*********************************************************
Compute de Option value by means of Path estimator
-
**********************************************************/
scalar pathEstimator(scalar *SPaths, scalar *ctrs, int *gIdx, scalar *ZPaths, scalar *EPaths, scalar X, scalar r, scalar T, scalar *bundle_a1){

	int st;
	scalar CashFlows[SIMS];
	int ExerciseTime[SIMS];
	scalar ContinuationValue[SIMS];
	int mat_gIdx[SIMS*(STEPS-1)];
	scalar dt = T/STEPS;

	Set_VectorInt(ExerciseTime, SIMS, STEPS);
	Set_Vector(ContinuationValue, SIMS, 0.0);
	Set_Vector(CashFlows, SIMS, 0.0);

	for(st = STEPS-1; st>0; st--){
		reBundling(SPaths, gIdx, ctrs, st);
		SetVectorTo2DMatrix_int(mat_gIdx, 1, st-1, SIMS, STEPS-1, gIdx);
	}

	computeCV(ZPaths, EPaths, CashFlows, ExerciseTime, bundle_a1, mat_gIdx, X, r, dt, ContinuationValue);

	return discountCashFlows(CashFlows, ExerciseTime, r, dt);

}//end pathEstimator

int main(int argc, char **argv){

	//Number of paths
	/*unsigned int directE_sims, i, j;
	if (argc > 1)
		directE_sims = atoi(argv[1]);
	else
        directE_sims = 2*MILION;*/

	//Input data
	scalar K, r, q, T, rho;

	K = 40;
	r = 0.06;
	q = 0.0;
	T = 1.0;
	rho = 0.25;

	/*
	scalar S0[] = {40, 40, 40, 40, 40};//Can be different
	scalar sigma[] = {0.2, 0.2, 0.2, 0.2, 0.2};

	scalar Rho[] = {1.0, rho, rho, rho, rho,
					rho, 1.0, rho, rho, rho,
					rho, rho, 1.0, rho, rho,
					rho, rho, rho, 1.0, rho,
					rho, rho, rho, rho, 1.0 };
	*/

// 	cTic();
	gsl_vector *S0_gsl = gsl_vector_alloc(ASSETS);
	gsl_vector *sigma_gsl = gsl_vector_alloc(ASSETS);
	gsl_matrix *Rho_gsl = gsl_matrix_alloc(ASSETS, ASSETS);
	gsl_matrix_set_identity(Rho_gsl);

	unsigned int i, j;
	for(i=0;i<ASSETS;i++){
		gsl_vector_set(S0_gsl, i, 40);
		gsl_vector_set(sigma_gsl, i, 0.2);
		for(j=0;j<ASSETS;j++)
			 if(i!=j)
				gsl_matrix_set(Rho_gsl, i, j, rho);
	}
	
// 	gsl_matrix_fprintf (stdout, Rho_gsl, "%g");

	gsl_matrix *C_gsl = gsl_matrix_alloc(ASSETS, ASSETS);
	gsl_matrix_memcpy(C_gsl, Rho_gsl);
	gsl_linalg_cholesky_decomp(C_gsl);
// 	for(i=0;i<ASSETS;i++)
// 		for(j=i+1;j<ASSETS;j++)
// 			gsl_matrix_set(C_gsl, i, j, 0.0);
// 	gsl_matrix_fprintf(stdout, C_gsl, "%g");
// 	printf("\nC = \n");printMatrix(C_gsl->data, ASSETS, ASSETS);
// 	cToc("gsl");

// 	cTic();
// 	scalar S0[ASSETS], sigma[ASSETS], Rho[ASSETS*ASSETS];
// // 	unsigned int i, j;
// 	for(i=0;i<ASSETS;i++){
// 		S0[i] = 40;
// 		sigma[i] = 0.2;
// 		for(j=0;j<ASSETS;j++)
// 			 if(i==j)
// 				 Rho[pos2Dto1D(i,j,ASSETS)] = 1.0;
// 			 else
// 				 Rho[pos2Dto1D(i,j,ASSETS)] = rho;
// 	}
// 	
// 	//printf("S0 = "); printVector(S0, ASSETS);
// 	//printf("sigma = ");printVector(sigma, ASSETS);
// // 	printf("\nsig = \n");printMatrix(Rho, ASSETS, ASSETS);

// 	scalar C[ASSETS*ASSETS];
// 	Cholesky(Rho, ASSETS, C);

// 	printf("\nC = \n");printMatrix(C, ASSETS, ASSETS);
// 	cToc("myMath");

	scalar dE[TRIALS];
	scalar pE[TRIALS];
	scalar sum_dE = 0.0;
	scalar sum_pE = 0.0;
	scalar *SPaths, *ZPaths, *EPaths;
	scalar ctrs[BUNDS*ASSETS*STEPS];
	int *gIdx;
	scalar bundle_a1[BASISF*BUNDS*STEPS];

GcTic();

	gen = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(gen, 1234);
	for(i=0;i<TRIALS;i++){

		/**********************BUNDLING***********************/
		//SIMS = directE_sims;
		SIMS = bSIMS;

		//3D matrix of paths. Store the whole data.
		SPaths = (scalar *)malloc(SIMS*ASSETS*(STEPS+1)*sizeof(scalar));
		ZPaths = (scalar *)malloc(SIMS*(STEPS+1)*sizeof(scalar));
		EPaths = (scalar *)malloc(SIMS*BASISF*(STEPS+1)*sizeof(scalar));
		
		gsl_matrix *SPaths_g[STEPS+1]; alloc_array_gsl_matrix(SPaths_g, STEPS+1, ASSETS, SIMS);
		gsl_vector *ZPaths_g[STEPS]; alloc_array_gsl_vector(ZPaths_g, STEPS, SIMS);
		gsl_matrix *EPaths_g[STEPS+1]; alloc_array_gsl_matrix(EPaths_g, STEPS+1, BASISF, SIMS);
		unsigned int BundIdx[STEPS+1][SIMS];
		gsl_matrix *bundle_a[STEPS]; alloc_array_gsl_matrix(bundle_a, STEPS, BUNDS, BASISF);
cTic();
		//Training paths for bundling
		MonteCarlo_MultiAsset(r-q, T, S0_gsl, sigma_gsl, Rho_gsl, C_gsl, SPaths, ZPaths, EPaths, SPaths_g, ZPaths_g, EPaths_g);
cToc("MCbund");

		printf("\nReference value: %lf.\n\n", SPaths[pos3Dto1D(0, 0, STEPS, ASSETS, STEPS+1)]);
		printf("\nReference value: %lf.\n\n", ZPaths[pos2Dto1D(0, STEPS, STEPS+1)]);
		printf("\nReference value: %lf.\n\n", EPaths[pos3Dto1D(0, 1, STEPS, BASISF, STEPS+1)]);
cTic();
// 		gIdx = (int *)malloc(SIMS*sizeof(int));
// 		Bundling(SPaths, gIdx, ctrs);
// 		Bundling(SPaths_g, gIdx, ctrs_g);
		Bundling(SPaths_g, BundIdx);
// 		printVectorInt(BundIdx[1], SIMS);
cToc("Bund");
		
		directEstimator_gsl(SPaths_g, BundIdx, ZPaths_g, EPaths_g, K, r, T, bundle_a);
		
		free(SPaths);
		free(ZPaths);
		free(EPaths);
		free(gIdx);

		free_array_gsl_matrix(SPaths_g, STEPS+1);
		free_array_gsl_vector(ZPaths_g, STEPS);
		free_array_gsl_matrix(EPaths_g, STEPS+1);
		/**********************BUNDLING***********************/

		/**********************DIRECT ESTIMATOR***********************
		SIMS = eSIMS;

		printf("\nSIMS = %ld | STEPS = %d | EXERS = %d | ASSETS = %d | BUNDS = %d | RATIO = %d\n", SIMS, STEPS, EXERS, ASSETS, BUNDS, SIMS/BUNDS);

		SPaths = (scalar *)malloc(SIMS*ASSETS*(STEPS+1)*sizeof(scalar));
		ZPaths = (scalar *)malloc(SIMS*(STEPS+1)*sizeof(scalar));
		EPaths = (scalar *)malloc(SIMS*BASISF*(STEPS+1)*sizeof(scalar));
cTic();
		//Monte Carlo paths for simulation
		MonteCarlo_MultiAsset(r-q, T, S0_gsl, sigma_gsl, Rho_gsl, C_gsl, SPaths, ZPaths, EPaths);
cToc("MCdirect");
// 		ZPaths = (scalar *)malloc(SIMS*(STEPS+1)*sizeof(scalar));
// cTic();
// 		computeZPaths(SPaths, ZPaths);
// cToc("ZPdirect");
		gIdx = (int *)malloc(SIMS*sizeof(int));
cTic();
		dE[i] = directEstimator(SPaths, ctrs, gIdx, ZPaths, EPaths, K, r, T, bundle_a1);
		sum_dE += dE[i];
cToc("directEst");
		printf("\nDirect Estimator: %lf\n\n", dE[i]);

		//free(SPaths);
		//free(ZPaths);
		//free(gIdx);
		/**********************DIRECT ESTIMATOR***********************/

		/**********************PATH ESTIMATOR***********************
		//SIMS = eSIMS;

		//SPaths = (scalar *)malloc(SIMS*ASSETS*(STEPS+1)*sizeof(scalar));
		//gIdx = (int *)malloc(SIMS*sizeof(int));
cTic();
		//Monte Carlo paths for simulation
		MonteCarlo_MultiAsset(r-q, T, S0_gsl, sigma_gsl, Rho_gsl, C_gsl, SPaths, ZPaths, EPaths);
cToc("MCpath");
// 		ZPaths = (scalar *)malloc(SIMS*(STEPS+1)*sizeof(scalar));
// cTic();
// 		computeZPaths(SPaths, ZPaths);
// cToc("ZPpath");
// 		printMatrix(ZPaths, SIMS, STEPS+1);
cTic();
		pE[i] = pathEstimator(SPaths, ctrs, gIdx, ZPaths, EPaths, K, r, T, bundle_a1);
		sum_pE += pE[i];
cToc("pathEst");
		printf("\nPath Estimator: %lf\n\n", pE[i]);

		free(SPaths);
		free(ZPaths);
		free(EPaths);
		free(gIdx);
		/**********************PATH ESTIMATOR***********************/

	}//end for trials

	gsl_rng_free(gen);

GcToc("Final");

	printf("------------------MEAN-------------------");
	printf("\nDirect Estimator: %lf\n", sum_dE/TRIALS);
	printf("\nPath Estimator: %lf\n", sum_pE/TRIALS);

	printf("------------------STD-------------------");
	printf("\nDirect Estimator: %.10lf\n", stdeviation(dE, TRIALS, sum_dE/TRIALS));
	printf("\nPath Estimator: %.10lf\n", stdeviation(pE, TRIALS, sum_pE/TRIALS));

	printf("\n\n");

}//end main