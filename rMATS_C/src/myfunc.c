#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics_double.h>
#include "../include/type.h"
#include "../include/util.h"
#include "../include/myfunc.h"
#include "../include/global.h"


extern double rho;
extern double cutoff;
extern clock_t dur;


#define ITER_CUTOFF 0.01


/**
 * @brief invoke this function before return from
 *        MLE_marginal_iteration and MLE_marginal_iteration_constrain.
 *
 * @param mle
 * @param sum
 * @param psi1
 * @param psi2
 * @param beta0
 * @param beta1
 * @param var1
 * @param var2
 */
void mle_result_set(mle_result* mle, double sum, gsl_vector* psi1, gsl_vector* psi2,
                    double beta0, double beta1, double var1, double var2) {
    mle->sum = sum;
    mle->params.psi1 = psi1;
    mle->params.psi2 = psi2;
    mle->params.beta0 = beta0;
    mle->params.beta1 = beta1;
    mle->params.var1 = var1;
    mle->params.var2 = var2;
    return;
}


double myfunc_multivar(const double x[], va_list argv) {
    gsl_vector *psi1 = va_arg(argv, gsl_vector*);
    gsl_vector *psi2 = va_arg(argv, gsl_vector*);
    double var1 = va_arg(argv, double), var2 = va_arg(argv, double);
    double sum1, sum2;

    // TODO  memcpy replaced.
    sum1 = cuscumsum(psi1, sum_for_multivar, 1, x[0]);
    sum2 = cuscumsum(psi2, sum_for_multivar, 1, x[1]);
    sum1 = sum1/(var1*2), sum2 = sum2/(var2*2);

    return sum1 + sum2 + 0.05 * pow(rho, 2) / (1 - pow(rho, 2)) * \
        (pow(gsl_cdf_gaussian_Pinv(x[0], 1), 2) + pow(gsl_cdf_gaussian_Pinv(x[1], 1), 2) - \
         2 * rho * gsl_cdf_gaussian_Pinv(x[0], 1) * gsl_cdf_gaussian_Pinv(x[1], 1));
}

void myfunc_multivar_der(const double x[], double res[], va_list argv) {
    gsl_vector *psi1 = va_arg(argv, gsl_vector*);
    gsl_vector *psi2 = va_arg(argv, gsl_vector*);
    double var1 = va_arg(argv, double), var2 = va_arg(argv, double);
    double sum1, sum2, tmp = pow(rho,2)/(1-pow(rho,2));
    double ppfx0 = gsl_cdf_gaussian_Pinv(x[0], 1);
    double ppfx1 = gsl_cdf_gaussian_Pinv(x[1], 1);

    sum1 = cuscumsum(psi1, sum_for_multivar_der, 1, x[0]);
    sum2 = cuscumsum(psi2, sum_for_multivar_der, 1, x[1]);
    sum1 = sum1/(var1*2), sum2 = sum2/(var2*2);

    // TODO are there any better way to convert cdf to pdf?
    res[0] = sum1 + tmp * 0.1*(ppfx0-rho*ppfx1)/gsl_ran_gaussian_pdf(ppfx0, 1);
    res[1] = sum2 + tmp * 0.1*(ppfx1-rho*ppfx0)/gsl_ran_gaussian_pdf(ppfx1, 1);

    return;
}


double myfunc_1_2(const double x[], va_list argv) {
    int flag = va_arg(argv, int);
    gsl_vector *psi1 = va_arg(argv, gsl_vector*);
    gsl_vector *psi2 = va_arg(argv, gsl_vector*);
    double var1 = va_arg(argv, double), var2 = va_arg(argv, double);
    double sum1, sum2;

    if (flag == 1) {
        sum1 = cuscumsum(psi1, sum_for_multivar, 1, x[0]+cutoff);
        sum2 = cuscumsum(psi2, sum_for_multivar, 1, x[0]);
    } else {
        sum1 = cuscumsum(psi1, sum_for_multivar, 1, x[0]);
        sum2 = cuscumsum(psi2, sum_for_multivar, 1, x[0]+cutoff);
    }
    sum1 = sum1/(var1*2), sum2 = sum2/(var2*2);

    return sum1 + sum2 + 0.05 * pow(rho, 2) / (1 - pow(rho, 2)) * \
        (pow(gsl_cdf_gaussian_Pinv(x[0]+cutoff, 1), 2) + pow(gsl_cdf_gaussian_Pinv(x[0], 1), 2) - \
         2 * rho * gsl_cdf_gaussian_Pinv(x[0]+cutoff, 1) * gsl_cdf_gaussian_Pinv(x[0], 1));
}


void myfunc_der_1_2(const double x[], double res[], va_list argv) {
    int flag = va_arg(argv, int);
    gsl_vector *psi1 = va_arg(argv, gsl_vector*);
    gsl_vector *psi2 = va_arg(argv, gsl_vector*);
    double var1 = va_arg(argv, double), var2 = va_arg(argv, double);
    double sum1, sum2, tmp = pow(rho,2)/(1-pow(rho,2));
    double ppfx0 = gsl_cdf_gaussian_Pinv(x[0], 1);
    double ppfxc = gsl_cdf_gaussian_Pinv(x[0]+cutoff, 1);

    if (flag == 1) {
        sum1 = cuscumsum(psi1, sum_for_multivar_der, 1, x[0]+cutoff);
        sum2 = cuscumsum(psi2, sum_for_multivar_der, 1, x[0]);
    } else {
        sum1 = cuscumsum(psi1, sum_for_multivar_der, 1, x[0]);
        sum2 = cuscumsum(psi2, sum_for_multivar_der, 1, x[0]+cutoff);
    }
    sum1 = sum1/(var1*2), sum2 = sum2/(var2*2);

    // TODO are there any better way to convert cdf to pdf?
    res[0] = sum1 + tmp * 0.1*(ppfxc-rho*ppfx0)/gsl_ran_gaussian_pdf(ppfxc, 1);
    res[1] = sum2 + tmp * 0.1*(ppfx0-rho*ppfxc)/gsl_ran_gaussian_pdf(ppfx0, 1);
    res[0] += res[1];

    return;
}


double myfunc_marginal(const double x[], va_list argv) {
    gsl_vector* I = va_arg(argv, gsl_vector*), *S = va_arg(argv, gsl_vector*);
    gsl_vector* psi = va_arg(argv, gsl_vector*);
    double var = va_arg(argv, double), sum = 0;
    int inclu_len = va_arg(argv, int), skip_len = va_arg(argv, int), idx = 0;

    sum = cuscumsum(psi, sum_for_marginal, 7, &idx, x[0], I, S, var, inclu_len, skip_len);

    return sum;
}


void myfunc_marginal_der(const double x[], double res[], va_list argv) {
    gsl_vector* I = va_arg(argv, gsl_vector*), *S = va_arg(argv, gsl_vector*);
    gsl_vector* psi = va_arg(argv, gsl_vector*);
    double var = va_arg(argv, double), sum;
    int inclu_len = va_arg(argv, int), skip_len = va_arg(argv, int), idx = 0;

    sum = cuscumsum(psi, sum_for_marginal_der, 7, &idx, x[0], I, S, var, inclu_len, skip_len);
    res[0] = sum;

    return;
}


double myfunc_marginal_1_2(const double x[], va_list argv) {
    int flag = va_arg(argv, int);
    double beta1, beta2;
    if (flag == 1) {
        beta1 = x[0] + cutoff, beta2 = x[0];
    } else {
        beta1 = x[0], beta2 = x[0] + cutoff;
    }
    // TODO parameter ordering
    return myfunc_marginal(&beta1, argv) + myfunc_marginal(&beta2, argv);
}


void myfunc_marginal_1_2_der(const double x[], double res[], va_list argv) {
    int flag = va_arg(argv, int);
    double beta1, beta2;
    double tmp1, tmp2;
    if (flag == 1) {
        beta1 = x[0] + cutoff, beta2 = x[0];
    } else {
        beta1 = x[0], beta2 = x[0] + cutoff;
    }
    // TODO parameter ordering
    myfunc_marginal_der(&beta1, &tmp1, argv);
    myfunc_marginal_der(&beta2, &tmp2, argv);
    res[0] = tmp1 + tmp2;

    return;
}


double myfunc_individual(const double x[], va_list argv) {
    double I = va_arg(argv, double), S = va_arg(argv, double);
    double beta = va_arg(argv, double), var = va_arg(argv, double);
    int inclu_len = va_arg(argv, int), skip_len = va_arg(argv, int);
    double new_psi = inclu_len * x[0]/(inclu_len * x[0] + skip_len * (1 - x[0]));

    // TODO This change the result.
    return pow((logit(x[0]) - logit(beta)), 2)/(2 * var) -
           (I * log(new_psi) + S * log(1 - new_psi) - log(x[0]) - log(1-x[0]));
}


void myfunc_individual_der(const double x[], double res[], va_list argv) {
    double I = va_arg(argv, double), S = va_arg(argv, double);
    double beta = va_arg(argv, double), var = va_arg(argv, double);
    int inclu_len = va_arg(argv, int), skip_len = va_arg(argv, int);
    double new_psi = inclu_len * x[0]/(inclu_len * x[0] + skip_len * (1 - x[0]));
    double new_psi_der = inclu_len * skip_len/pow(inclu_len * x[0] + skip_len * (1 - x[0]), 2);

    res[0] = 1/x[0] + S/(1 - new_psi) * new_psi_der + \
             (logit(x[0]) - logit(beta))/(var * x[0] * (1 - x[0])) - (1/(1-x[0]) + I/new_psi * new_psi_der);

    return;
}


#define individual_for_loop() \
    for (idx = 0; idx < psi1->size; ++idx) { \
        x[0] = gsl_vector_get(psi1, idx); \
        cur_sum += l_bfgs_b_wrapper(n, m, x, l, u, nbd, \
                                    myfunc_individual, myfunc_individual_der, \
                                    factr, pgtol, iprint, 15000, 15000, \
                                    6, gsl_vector_get(i1, idx), gsl_vector_get(s1, idx), \
                                    beta0, var1, inclu_len, skip_len); \
        gsl_vector_set(psi1, idx, x[0]); \
    } \
    for (idx = 0; idx < psi2->size; ++idx) { \
        x[0] = gsl_vector_get(psi2, idx); \
        cur_sum += l_bfgs_b_wrapper(n, m, x, l, u, nbd, \
                                    myfunc_individual, myfunc_individual_der, \
                                    factr, pgtol, iprint, 15000, 15000, \
                                    6, gsl_vector_get(i2, idx), gsl_vector_get(s2, idx), \
                                    beta1, var2, inclu_len, skip_len); \
        gsl_vector_set(psi2, idx, x[0]); \
    } \

int MLE_marginal_iteration(gsl_vector* i1, gsl_vector* i2,
                           gsl_vector* s1, gsl_vector* s2,
                           const int inclu_len, const int skip_len,
                           mle_result* mle) {
    gsl_vector* psi1 = gsl_vector_alloc(i1->size); // TODO deallocation.
    gsl_vector* psi2 = gsl_vector_alloc(i2->size);
    vec2psi(psi1, i1, s1, inclu_len, skip_len);
    vec2psi(psi2, i2, s2, inclu_len, skip_len);
    int iter_max = 100, count = 0, i = 0;
    size_t idx = 0;
    double prev_sum = 0, cur_sum = 0;
    double var1 = 10 * gsl_stats_variance(psi1->data, 1, psi1->size) * (psi1->size-1)/(psi1->size);
    double var2 = 10 * gsl_stats_variance(psi2->data, 1, psi2->size) * (psi2->size-1)/(psi2->size);
    double beta0, beta1, iter_cutoff = 1;
    integer n, m = 10, nbd[nmax];
    doublereal x[nmax], l[nmax], u[nmax], factr=1.0e7, pgtol=1.0e-5, iprint=-1;

    if (var1 <= 0.01 || psi1->size == 1) {
        var1 = 0.01;
    }
    if (var2 <= 0.01 || psi2->size == 1) {
        var2 = 0.01;
    }

    // TODO & or && ?
    // According to original python code 'while((iter_cutoff>0.01)&(count<=iter_maxrun)):',
    // it's a bit arithmetic '&'. However,  it should be a logical 'and' in such senario.
    while((iter_cutoff > ITER_CUTOFF) && (count <= iter_max)) {
        ++count, n = 2;
        x[0] = gsl_stats_mean(psi1->data, 1, psi1->size);
        x[1] = gsl_stats_mean(psi2->data, 1, psi2->size);
        for (i = 0; i < n; ++i) {
            nbd[i] = 2;
            l[i] = 0.01;
            u[i] = 0.99;
        }
        l_bfgs_b_wrapper(n, m, x, l, u, nbd,
                         myfunc_multivar, myfunc_multivar_der, factr, pgtol,
                         iprint, 15000, 15000, 4, psi1, psi2, var1, var2);
        beta0 = x[0], beta1 = x[1];
        n = 1, cur_sum = 0;

        individual_for_loop();

        if (count > 1) {
            iter_cutoff = fabs(prev_sum - cur_sum);
        }
        prev_sum = cur_sum;
    }
    if (count > iter_max) {
        mle_result_set(mle, cur_sum, psi1, psi2, 0, 0, var1, var2);
        return 0;
    }

    iter_cutoff = 1, iter_max = 100, count = 0, prev_sum = 0;
    while((iter_cutoff > ITER_CUTOFF) && (count <= iter_max)) {
        ++count, n = 1;
        x[0] = beta0;
        l_bfgs_b_wrapper(n, m, x, l, u, nbd,
                         myfunc_marginal, myfunc_marginal_der,
                         factr, pgtol, iprint, 15000, 15000,
                         6, i1, s1, psi1, var1, inclu_len, skip_len);
        beta0 = x[0];
        x[0] = beta1;
        l_bfgs_b_wrapper(n, m, x, l, u, nbd,
                         myfunc_marginal, myfunc_marginal_der,
                         factr, pgtol, iprint, 15000, 15000,
                         6, i2, s2, psi2, var2, inclu_len, skip_len);
        beta1 = x[0];

        cur_sum = 0;

        individual_for_loop();

        if (count > 1) {
            iter_cutoff = fabs(prev_sum - cur_sum);
        }
        prev_sum = cur_sum;
    }

    gsl_vector_free(psi1);
    gsl_vector_free(psi2);

    if (count > iter_max) {
        mle_result_set(mle, cur_sum, psi1, psi2, 0, 0, var1, var2);
        return 0;
    }

    mle_result_set(mle, cur_sum, psi1, psi2, beta0, beta1, var1, var2);

    return 0;
}


int MLE_marginal_iteration_constrain(gsl_vector* i1, gsl_vector* i2,
                                     gsl_vector* s1, gsl_vector* s2,
                                     const int inclu_len, const int skip_len,
                                     mle_result* mle) {
    gsl_vector* psi1 = gsl_vector_alloc(i1->size); // TODO deallocation.
    gsl_vector* psi2 = gsl_vector_alloc(i2->size);
    vec2psi(psi1, i1, s1, inclu_len, skip_len);
    vec2psi(psi2, i2, s2, inclu_len, skip_len);
    int iter_max = 100, count = 0, i = 0;
    size_t idx = 0;
    double prev_sum = 0, cur_sum = 0, iter_cutoff = 1;
    double var1 = 10 * gsl_stats_variance(psi1->data, 1, psi1->size) * (psi1->size-1)/(psi1->size);
    double var2 = 10 * gsl_stats_variance(psi2->data, 1, psi2->size) * (psi2->size-1)/(psi2->size);
    double beta0, beta1;
    integer n = 1, m = 10, nbd[nmax];
    doublereal x[nmax], l[nmax], u[nmax], factr=1.0e7, pgtol=1.0e-5, iprint=-1;

    if (var1 <= 0.01 || psi1->size == 1) {
        var1 = 0.01;
    }
    if (var2 <= 0.01 || psi2->size == 1) {
        var2 = 0.01;
    }

    // TODO & or && ?
    // According to original python code 'while((iter_cutoff>0.01)&(count<=iter_maxrun)):',
    // it's a bit arithmetic '&'. However,  it should be a logical 'and' in such senario.
    while((iter_cutoff > ITER_CUTOFF) && (count <= iter_max)) {
        ++count;
        x[0] = gsl_stats_mean(psi1->data, 1, psi1->size);
        x[1] = gsl_stats_mean(psi2->data, 1, psi2->size);
        for (i = 0; i < n; ++i) {
            nbd[i] = 2;
            l[i] = 0.001;
            u[i] = 0.999-cutoff;
        }
        if (x[0] > x[1]) {
            l_bfgs_b_wrapper(n, m, x+1, l, u, nbd,
                             myfunc_1_2, myfunc_der_1_2, factr, pgtol,
                             iprint, 15000, 15000, 5, 1, psi1, psi2, var1, var2);
            beta1 = fmax(fmin(x[1], 1-cutoff), 0);
            beta0 = beta1 + cutoff;
        } else {
            l_bfgs_b_wrapper(n, m, x, l, u, nbd,
                             myfunc_1_2, myfunc_der_1_2, factr, pgtol,
                             iprint, 15000, 15000, 5, 2, psi1, psi2, var1, var2);
            beta0 = fmax(fmin(x[0], 1-cutoff), 0);
            beta1 = beta0 + cutoff;
        }
        cur_sum = 0;
        for (i = 0; i < n; ++i) {
            nbd[i] = 2;
            l[i] = 0.01;
            u[i] = 0.99;
        }

        individual_for_loop();

        if (count > 1) {
            iter_cutoff = fabs(prev_sum - cur_sum);
        }
        prev_sum = cur_sum;
    }
    if (count > iter_max) {
        mle_result_set(mle, cur_sum, psi1, psi2, 0, 0, var1, var2);
        return 0;
    }

    iter_cutoff = 1, iter_max = 100, count = 0, prev_sum = 0;
    while((iter_cutoff > ITER_CUTOFF) && (count <= iter_max)) {
        ++count;
        for (i = 0; i < n; ++i) {
            nbd[i] = 2;
            l[i] = 0.001;
            u[i] = 0.999-cutoff;
        }
        x[0] = beta0;
        x[1] = beta1;
        if (gsl_stats_mean(psi1->data, 1, psi1->size) > gsl_stats_mean(psi2->data, 1, psi2->size)) {
            l_bfgs_b_wrapper(n, m, x+1, l, u, nbd,
                             myfunc_marginal_1_2, myfunc_marginal_1_2_der, factr, pgtol,
                             iprint, 15000, 15000, 13, 1, i1, s1, psi1, var1, inclu_len, skip_len,
                             i2, s2, psi2, var2, inclu_len, skip_len);
            beta1 = fmax(fmin(x[1], 1-cutoff), 0);
            beta0 = beta1 + cutoff;
        } else {
            l_bfgs_b_wrapper(n, m, x, l, u, nbd,
                             myfunc_marginal_1_2, myfunc_marginal_1_2_der, factr, pgtol,
                             iprint, 15000, 15000, 13, 2, i1, s1, psi1, var1, inclu_len, skip_len,
                             i2, s2, psi2, var2, inclu_len, skip_len);
            beta0 = fmax(fmin(x[0], 1-cutoff), 0);
            beta1 = beta0 + cutoff;
        }
        cur_sum = 0;
        for (i = 0; i < n; ++i) {
            nbd[i] = 2;
            l[i] = 0.01;
            u[i] = 0.99;
        }

        individual_for_loop();

        if (count > 1) {
            iter_cutoff = fabs(prev_sum - cur_sum);
        }
        prev_sum = cur_sum;
    }

    gsl_vector_free(psi1);
    gsl_vector_free(psi2);

    if (count > iter_max) {
        mle_result_set(mle, cur_sum, psi1, psi2, 0, 0, var1, var2);
        return 0;
    }

    mle_result_set(mle, cur_sum, psi1, psi2, beta0, beta1, var1, var2);

    return 0;
}


void* thread_wrapper_for_LT(void* arg) {
    double *ret = (double*)malloc(sizeof(double));
    odiff *data = (odiff*)arg;
    *ret = likelihood_test(data->inc1, data->inc2,
                           data->skp1, data->skp2,
                           data->inclu_len, data->skip_len,
                           data->flag, data->id);
    return (void*)ret;
}


double likelihood_test(gsl_vector *i1, gsl_vector *i2, gsl_vector *s1, gsl_vector *s2,
                       int inclu_len, int skip_len, int flag, char* id) {
    printf("Testing %s\n", id);
    if (!flag) {
        printf("1 return from: %s\n", id);
        return 1;
    } else {
        mle_result mle;
        MLE_marginal_iteration(i1, i2, s1, s2, inclu_len, skip_len, &mle);
        if (fabs(mle.params.beta0 - mle.params.beta1) <= cutoff) {
            printf("2 return from: %s\n", id);
			return 1;
        } else {
            mle_result mle_constrain;
            MLE_marginal_iteration_constrain(i1, i2, s1, s2, inclu_len, skip_len, &mle_constrain);
            printf("3 return from: %s\n", id);
            return 1 - gsl_cdf_chisq_P(2 * (fabs(mle_constrain.sum - mle.sum)), 1);
        }
    }
}


int vec2psi(gsl_vector* psi, gsl_vector *inc, gsl_vector *skp,
            int inclu_len, int skip_len) {
    size_t idx;
    for(idx = 0; idx < inc->size; ++idx) {
        if (gsl_vector_get(inc, idx) + gsl_vector_get(skp, idx) == 0) {
            gsl_vector_set(psi, idx, 0.5);
        } else {
            gsl_vector_set(psi, idx, gsl_vector_get(inc, idx)/
                    (gsl_vector_get(inc, idx) + inclu_len * gsl_vector_get(skp, idx)/skip_len));
        }
    }
    return 0;
}
