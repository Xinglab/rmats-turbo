#ifndef TRANS_MYFUNC
#define TRANS_MYFUNC


#include <stdarg.h>
#include <gsl/gsl_vector.h>
#include "../include/type.h"


void mle_result_set(mle_result* mle, double sum, gsl_vector* psi1, gsl_vector* psi2,
                    double beta0, double beta1, double var1, double var2);

double myfunc_multivar(const double x[], va_list argv);

void myfunc_multivar_der(const double x[], double res[], va_list argv);

double myfunc_1_2(const double x[], va_list argv);

void myfunc_der_1_2(const double x[], double res[], va_list argv);

double myfunc_marginal_1_2(const double x[], va_list argv);

void myfunc_marginal_1_2_der(const double x[], double res[], va_list argv);

double myfunc_individual(const double x[], va_list argv);

void myfunc_individual_der(const double x[], double res[], va_list argv);

double myfunc_marginal(const double x[], va_list argv);

void myfunc_marginal_der(const double x[], double res[], va_list argv);

int MLE_marginal_iteration(gsl_vector* i1, gsl_vector* i2,
                           gsl_vector* s1, gsl_vector* s2,
                           const int inclu_len, const int skip_len,
                           mle_result* mle);

int MLE_marginal_iteration_constrain(gsl_vector* i1, gsl_vector* i2,
                                     gsl_vector* s1, gsl_vector* s2,
                                     const int inclu_len, const int skip_len,
                                     mle_result* mle);

void* thread_wrapper_for_LT(void* arg);

double likelihood_test(gsl_vector *i1, gsl_vector *i2, gsl_vector *s1, gsl_vector *s2,
                        int inclu_len, int skip_len, int flag, char* id);

int vec2psi(gsl_vector* psi, gsl_vector *inc, gsl_vector *skp,
            int inclu_len, int skip_len);


#endif
