#ifndef TRANS_TYPE
#define TRANS_TYPE


#include <gsl/gsl_vector.h>

// corresponding python variable: element of 'list_n_original_diff'
typedef struct {
    gsl_vector *inc1;
    gsl_vector *inc2;
    gsl_vector *skp1;
    gsl_vector *skp2;
    int inclu_len;
    int skip_len;
    int flag;
    char* id;
} odiff;


typedef struct {
    int batch_size;
    void **datum;
} batch_datum;


// used to represent the reture value of
// MLE_marginal_iteration and MLE_marginal_iteration_constrain.
// corresponding python variable: '[current_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]'
typedef struct {
    double sum;
    struct {
        gsl_vector *psi1;
        gsl_vector *psi2;
        double beta0;
        double beta1;
        double var1;
        double var2;
    } params;
} mle_result;


// we use a linked list to represent a python list.
// corresponding python variable: 'list_n_original_diff'
typedef struct dnode{
    odiff* data;
    struct dnode* next;
    struct dnode* end;
} diff_list_node;

// we use fortran l_bfgs_b routine to solve the numerical optimization problem.
// (i.e. f = 0.0 at the optimal solution.)
// only the minimal intersection of C's and Fortran's many data types can be relied on:
// the following are some essential mapping between c type and fortran type
// counterpart of f2c.h
typedef int integer;
typedef float real;
typedef double doublereal;
typedef long int logical;


#endif
