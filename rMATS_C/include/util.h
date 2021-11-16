#ifndef TRANS_UTIL
#define TRANS_UTIL


#define prefetch(x) __builtin_prefetch(x)
#define prefetch_for_r(x) __builtin_prefetch(x, 0, 3)
#define prefetch_for_w(x) __builtin_prefetch(x, 1, 3)


#include <stdarg.h>
#include <gsl/gsl_vector.h>
#include "../include/type.h"

int setulb_(integer *, integer *, doublereal *, 
            doublereal *, doublereal *, integer *, doublereal *, doublereal *,
            doublereal *, doublereal *, doublereal *, integer *, char *, 
            integer *, char *, logical *, integer *, doublereal *, integer *);

double cuscumsum(gsl_vector *vec, double (*fun) (double, va_list), int argc, ...);

int parse_title(char* str, char** output);

int str_to_vector(char* str, gsl_vector** vec);

int parse_line(char* str, char** id, size_t* id_n, gsl_vector** inc1,
               gsl_vector** skp1, gsl_vector** inc2, gsl_vector** skp2,
               int* inclu_len, int* skip_len);

int parse_file(const char* filename, diff_list_node* list, char** title);

double logit(double i);

double cumprod(const gsl_vector* vec);

double sum_for_multivar(const double i, va_list argv);

double sum_for_multivar_der(const double i, va_list argv);

double myfunc_marginal_2der(const double x, const double I, const double S,
                            const double beta, const double var,
                            const int inclu_len, const int skip_len);

double sum_for_marginal(const double i, va_list argv);

double sum_for_marginal_der(const double i, va_list argv);

odiff* diff_alloc(gsl_vector* inc1, gsl_vector* inc2,
                  gsl_vector* skp1, gsl_vector* skp2,
                  int inclu_len, int skip_len, int flag, char* id);

int diff_append(diff_list_node* header, odiff* data);

void mp_threadpool(int nthread, int ntask, void* (*func)(void *), void** datum, void **ret);

double l_bfgs_b_wrapper(integer n, integer m, doublereal x[], doublereal l[],
                        doublereal u[], integer nbd[],
                        double (*fp) (const double x[], va_list argv),
                        void (*gp) (const double x[], double res[], va_list argv),
                        doublereal factr, doublereal pgtol, integer iprint,
                        int maxfun, int maxiter, int argc, ...);


#endif
