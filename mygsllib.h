#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_linalg.h>

/*********************************************************
Auxiliary GSL functions
**********************************************************/

gsl_rng *gen;

int randomGaussianMatrix(scalar m, scalar s, int nrows, int ncols, gsl_matrix *Z){

	unsigned int i,j;
	for(i=0;i<nrows;i++)
		for(j=0;j<ncols;j++)
// 			gsl_matrix_set(Z, i, j, m + gsl_ran_gaussian(gen, s));
			gsl_matrix_set(Z, i, j, m + gsl_ran_gaussian_ziggurat(gen, s));

	return 0;
}

int alloc_array_gsl_vector(gsl_vector **arr, int arrsize, int vecsize){

	int i;

	for(i = 0; i < arrsize; i++)
		arr[i] = gsl_vector_alloc(vecsize);

	return 0;
}

int free_array_gsl_vector(gsl_vector **arr, int arrsize){

	int i;

	for(i = 0; i < arrsize; i++)
		gsl_vector_free(arr[i]);

	return 0;
}

int alloc_array_gsl_permutation(gsl_permutation **arr, int arrsize, int vecsize){

	int i;

	for(i = 0; i < arrsize; i++)
		arr[i] = gsl_permutation_alloc(vecsize);

	return 0;
}

int free_array_gsl_permutation(gsl_permutation **arr, int arrsize){

	int i;

	for(i = 0; i < arrsize; i++)
		gsl_permutation_free(arr[i]);

	return 0;
}

int alloc_array_gsl_matrix(gsl_matrix **arr, int arrsize, int matrows, int matcols){

	int i;

	for(i = 0; i < arrsize; i++)
		arr[i] = gsl_matrix_alloc(matrows, matcols);

	return 0;
}

int free_array_gsl_matrix(gsl_matrix **arr, int arrsize){

	int i;

	for(i = 0; i < arrsize; i++)
		gsl_matrix_free(arr[i]);

	return 0;
}

int gsl_matrix_exp_elements(gsl_matrix *a){

	int i,j;

	for(i = 0; i < a->size1; i++)
		for(j = 0; j < a->size2; j++)
			gsl_matrix_set(a, i, j, exp(gsl_matrix_get(a, i, j)));

	return 0;
}

int gsl_matrix_log_elements(gsl_matrix *a){

	int i,j;

	for(i = 0; i < a->size1; i++)
		for(j = 0; j < a->size2; j++)
			gsl_matrix_set(a, i, j, log(gsl_matrix_get(a, i, j)));

	return 0;
}

int gsl_matrix_max_elements(gsl_matrix *a, double x){

	int i,j;

	for(i = 0; i < a->size1; i++)
		for(j = 0; j < a->size2; j++)
			gsl_matrix_set(a, i, j, MAX(gsl_matrix_get(a, i, j), x));

	return 0;
}

scalar gsl_vector_sum(gsl_vector *v){

	int i;

	scalar sum = 0.0;
	for(i = 0; i < v->size; i++)
		sum += gsl_vector_get(v, i);

	return sum;
}

int gsl_vector_max_value(gsl_vector *v, double x){

	int i;

	for(i = 0; i < v->size; i++)
		gsl_vector_set(v, i, MAX(gsl_vector_get(v, i), x));

	return 0;
}

int gsl_vector_max_elements(gsl_vector *v1, gsl_vector *v2){

	int i;

	for(i = 0; i < v1->size; i++)
		gsl_vector_set(v1, i, MAX(gsl_vector_get(v1, i), gsl_vector_get(v2, i)));

	return 0;
}
