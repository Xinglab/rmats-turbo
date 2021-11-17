#include <math.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <omp.h>
#include "../include/type.h"
#include "../include/util.h"
#include "../include/global.h"


extern clock_t dur;
extern char *optarg;
extern int optind;


/**
 * @brief initiator of odiff.
 *
 * @param inc1
 * @param inc2
 * @param skp1
 * @param skp2
 * @param inclu_len
 * @param skip_len
 * @param flag
 * @param id
 *
 * @return 
 */
odiff* diff_alloc(gsl_vector* inc1, gsl_vector* inc2,
                  gsl_vector* skp1, gsl_vector* skp2,
                  int inclu_len, int skip_len, int flag, char* id) {
    odiff *diff = (odiff*)malloc(sizeof(odiff));
    diff->inc1 = inc1, diff->skp1 = skp1, diff->inc2 = inc2, diff->skp2 = skp2;
    diff->inclu_len = inclu_len, diff->skip_len = skip_len;
    diff->flag = flag, diff->id = (char*)malloc(sizeof(char)*(strlen(id)+1));
    strcpy(diff->id, id);
    return diff;
}


/**
 * @brief cumulative summation of a gsl_vector.
 *
 * @param vec
 * @param fun
 * @param argc
 * @param ...
 *
 * @return 
 */
double cuscumsum(gsl_vector *vec, double (*fun) (double, va_list), int argc, ...) {
    va_list argv, tmp;
    double sum = 0;
    size_t i;
    va_start(argv, argc);
    for (i = 0; i < vec->size; ++i) {
        va_copy(tmp, argv);
        sum += fun(gsl_vector_get(vec, i), tmp);
    }
    va_end(tmp);
    va_end(argv);
    return sum;
}


int parse_title(char* str, char** output) {
    int idx = 0;
    char *token = strtok(str, " \t\r\n\v\f");
    while(token) {
        output[idx] = (char*)malloc(sizeof(char)*(strlen(token) + 1));
        strcpy(output[idx++], token);
        token = strtok(NULL, " \t\r\n\v\f");
    }

    return idx-1;
}


/**
 * @brief 
 *
 * @param input a string in format '1,2,3,4,5,6,...'.
 * @param vec
 *
 * @return 
 */
int str_to_vector(char* input, gsl_vector** vec) {
    int idx = 0, count = 1;
    size_t i = 0;

    for (i = 0; i < strlen(input); ++i) {
        if (input[i] == ',') {
            ++count;
        }
    }

    *vec = gsl_vector_alloc(count);
    char *str = strtok(input, ",");
    while(str) {
        if (strcmp(str, "NA") == 0 || strcmp(str, "NAN") == 0 ||
            strcmp(str, "na") == 0 || strcmp(str, "nan") == 0) {
            return -1;
        }
        gsl_vector_set(*vec, idx++, atof(str));
        str = strtok(NULL, ",");
    }

    return count;
}


int parse_line(char* str, char** id, size_t* id_n, gsl_vector** inc1,
               gsl_vector** skp1, gsl_vector** inc2, gsl_vector** skp2,
               int* inclu_len, int* skip_len) {
    int idx = 0;
    size_t found_id_len = 0;
    char *output[7];

    char *token = strtok(str, " \t\r\n\v\f");
    while(token) {
        output[idx++] = token;
        token = strtok(NULL, " \t\r\n\v\f");
    }

    found_id_len = strlen(output[0]);
    if ((*id_n) <= found_id_len) {
        *id_n = found_id_len + 1;
        if ((*id) != NULL) {
            free(*id);
        }
        *id = (char*)malloc(sizeof(char)*(*id_n));
    }
    strcpy(*id, output[0]);
    if (str_to_vector(output[1], inc1) == -1 ||
        str_to_vector(output[2], skp1) == -1 ||
        str_to_vector(output[3], inc2) == -1 ||
        str_to_vector(output[4], skp2) == -1) {
        printf("An error occured. rMATS cannot handle missing value. Sample Id: %s\n", *id);
        printf("Exiting.\n");
        exit(0);
        
    }
    *inclu_len = atoi(output[5]);
    *skip_len = atoi(output[6]);

    return idx-1;
}


int parse_file(const char* filename, diff_list_node* list, char** title_element_list) {
    FILE *ifp;
    char *str_line = NULL, *id = NULL;
    size_t str_line_n = 0, id_n = 0;
    int row_num=0, inclu_len, skip_len;
    gsl_vector *inc1 = NULL, *skp1 = NULL, *inc2 = NULL, *skp2 = NULL;

    if ((ifp = fopen(filename, "r")) == NULL) {
        printf("Fail to open!");
        return 0;
    }
    if (getline(&str_line, &str_line_n, ifp) == -1) {
        printf("Failed to read first line of input file.\n");
        printf("Exiting\n");
        exit(0);
    }
    parse_title(str_line, title_element_list);

    while (getline(&str_line, &str_line_n, ifp) != -1) {
        ++row_num;
        parse_line(str_line, &id, &id_n, &inc1, &skp1, &inc2, &skp2, &inclu_len, &skip_len);
        if (inc1->size != skp1->size || inc2->size != skp2->size) {
            printf("An error occured. The length of the pair of vector should be equal. Sample Id: %s\n", id);
            printf("Size of vector: %ld, %ld, %ld, %ld\n", inc1->size, skp1->size, inc2->size, skp2->size);
            printf("Exiting\n");
            exit(0);
        }
        gsl_vector_add(inc1, skp1), gsl_vector_add(inc2, skp2);

        // TODO | or || ?
        // According to original python code 'if (vecprod(vecadd(inc1,skp1))==0) | (vecprod(vecadd(inc2,skp2))==0):',
        // it's a bit arithmetic '|'. However,  it should be a logical 'or' in such senario.
        if (cumprod(inc1) == 0 || cumprod(inc2) == 0) {
            gsl_vector_sub(inc1, skp1), gsl_vector_sub(inc2, skp2);
            diff_append(list, diff_alloc(inc1, inc2, skp1, skp2, inclu_len, skip_len, 0, id));
        } else {
            gsl_vector_sub(inc1, skp1), gsl_vector_sub(inc2, skp2);

            // TODO According to original python code, this function will not change anything.
            // original comment: add 1 in both inclusion and skipping counts for robustness in small sample size
            gsl_vector_add_constant(inc1, 0.0), gsl_vector_add_constant(skp1, 0.0);
            gsl_vector_add_constant(inc2, 0.0), gsl_vector_add_constant(skp2, 0.0);
            diff_append(list, diff_alloc(inc1, inc2, skp1, skp2, inclu_len, skip_len, 1, id));
        }
    }
    if (str_line != NULL) {
        free(str_line);
    }
    if (id != NULL) {
        free(id);
    }
    fclose(ifp);

    return row_num;
}


/**
 * @brief handy function.
 *
 * @param i
 * @param argv
 *
 * @return 
 */
double logit(double i) {
    if (i < 0.01) {
        i = 0.01;
    } else if (i > 0.99) {
        i = 0.99;
    }
    return log(i/(1-i));
}


/**
 * @brief cumulative production of a gsl_vector.
 *
 * @param vec
 *
 * @return 
 */
double cumprod(const gsl_vector* vec) {
    double res = 1;
    size_t i;
    for (i = 0; i < vec->size; ++i) {
        res *= gsl_vector_get(vec, i);
    }
    return res;
}


// for multivar, 1, 2
double sum_for_multivar(const double i, va_list argv) {
    double pa = va_arg(argv, double);
    return pow(logit(i)-logit(pa), 2);
}


// for multivar_der, 1_der, 2_der
double sum_for_multivar_der(const double i, va_list argv) {
    double pa = va_arg(argv, double);
    return 2 * (logit(i) - logit(pa))/(pa*pa-pa);
}


double myfunc_marginal_2der(const double x, const double I, const double S,
                            const double beta, const double var,
                            const int inclu_len, const int skip_len) {
    double tmp1, tmp2, tmp3, one_x = 1-x;
    double powx = pow(x,2), pow1_x = pow(one_x,2), pow_len = pow(inclu_len*x + skip_len*one_x,2);

    tmp1 = (((2 * x - 1) * (logit(beta) - logit(x)) - 1)/var - 1)/(powx*pow1_x);
    tmp2 = I * skip_len * (2*inclu_len*x+skip_len)/(powx * pow_len);
    tmp3 = S * inclu_len * (inclu_len+2*skip_len*one_x)/(pow1_x * pow_len);

    return tmp1 - tmp2 - tmp3;
}


// for marginal
double sum_for_marginal(const double i, va_list argv) {
    int *idx = va_arg(argv, int*);
    double beta = va_arg(argv, double);
    double I_ = gsl_vector_get(va_arg(argv, gsl_vector*), *idx);
    double S_ = gsl_vector_get(va_arg(argv, gsl_vector*), *idx);
    double var = va_arg(argv, double), new_psi, f1, f1_2der;
    int inclu_len = va_arg(argv, int), skip_len = va_arg(argv, int);
    *idx += 1;

    new_psi = inclu_len * i/(inclu_len * i + skip_len * (1 - i));
    f1 = I_ * log(new_psi) + S_ * log(1 - new_psi) -
         pow(logit(i) - logit(beta),2)/(2*var) - log(i) - log(1-i);
    f1_2der = fabs(myfunc_marginal_2der(i, I_, S_, beta,
                                        var, inclu_len, skip_len));

    return 0.5 * log(fabs(f1_2der) + 0.00001) - f1;
}


// for marginal_der
double sum_for_marginal_der(const double i, va_list argv) {
    int *idx = va_arg(argv, int*);
    double beta = va_arg(argv, double);
    double I_ = gsl_vector_get(va_arg(argv, gsl_vector*), *idx);
    double S_ = gsl_vector_get(va_arg(argv, gsl_vector*), *idx);
    double var = va_arg(argv, double);
    int inclu_len = va_arg(argv, int), skip_len = va_arg(argv, int);
    double f1_1der, f1_2der, f1_3der, one_i = 1-i;
    *idx += 1;
    double powi = pow(i,2), pow1_i = pow(one_i,2);

    f1_3der = (2 * i - 1)/(powi * pow1_i*beta*(1-beta)*var);
    f1_2der = myfunc_marginal_2der(i, I_, S_, beta, var, inclu_len, skip_len);
    f1_1der = (logit(i) - logit(beta))/(beta * (1-beta) * var);

    return 0.5 * f1_3der/f1_2der - f1_1der;
}


// function used to manipulate our linked list.
int diff_append(diff_list_node* header, odiff* data) {
    diff_list_node *tmp = (diff_list_node*)malloc(sizeof(diff_list_node));
    tmp->data = data;
    tmp->end = tmp;
    tmp->next = NULL;
    header->end->next = tmp;
    header->end->end = tmp;
    header->end = tmp;

    return 0;
}


void mp_threadpool(int nthread, int ntask, void* (*func)(void *), void** datum, void **ret) {
    int i;
    omp_set_num_threads(nthread);

    #pragma omp parallel for private(i)
    for (i = 0; i < ntask; ++i) {
        ret[i] = (*func)(datum[i]);
    }

    return;
}


/**
 * @brief C wrapper of Fortran l_bfgs_b routine.
 *
 * @param n
 * @param m
 * @param x[]
 * @param l[]
 * @param u[]
 * @param nbd[]
 * @param fp
 * @param gp
 * @param factr
 * @param pgtol
 * @param iprint
 * @param maxiter
 * @param argc
 * @param ...
 *
 * @return 
********************************************************************
c     --------------------------------------------------------------
c             DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------
c
c     n is an INTEGER variable that must be set by the user to the
c       number of variables.  It is not altered by the routine.
c
c     m is an INTEGER variable that must be set by the user to the
c       number of corrections used in the limited memory matrix.
c       It is not altered by the routine.  Values of m < 3  are
c       not recommended, and large values of m can result in excessive
c       computing time. The range  3 <= m <= 20 is recommended. 
c
c     x is a DOUBLE PRECISION array of length n.  On initial entry
c       it must be set by the user to the values of the initial
c       estimate of the solution vector.  Upon successful exit, it
c       contains the values of the variables at the best point
c       found (usually an approximate solution).
c
c     l is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the lower bounds on the variables. If
c       the i-th variable has no lower bound, l(i) need not be defined.
c
c     u is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the upper bounds on the variables. If
c       the i-th variable has no upper bound, u(i) need not be defined.
c
c     nbd is an INTEGER array of dimension n that must be set by the
c       user to the type of bounds imposed on the variables:
c       nbd(i)=0 if x(i) is unbounded,
c              1 if x(i) has only a lower bound,
c              2 if x(i) has both lower and upper bounds, 
c              3 if x(i) has only an upper bound.
c
c     f is a DOUBLE PRECISION variable.  If the routine setulb returns
c       with task(1:2)= 'FG', then f must be set by the user to
c       contain the value of the function at the point x.
c
c     g is a DOUBLE PRECISION array of length n.  If the routine setulb
c       returns with taskb(1:2)= 'FG', then g must be set by the user to
c       contain the components of the gradient at the point x.
c
c     factr is a DOUBLE PRECISION variable that must be set by the user.
c       It is a tolerance in the termination test for the algorithm.
c       The iteration will stop when
c
c        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
c
c       where epsmch is the machine precision which is automatically
c       generated by the code. Typical values for factr on a computer
c       with 15 digits of accuracy in double precision are:
c       factr=1.d+12 for low accuracy;
c             1.d+7  for moderate accuracy; 
c             1.d+1  for extremely high accuracy.
c       The user can suppress this termination test by setting factr=0.
c
c     pgtol is a double precision variable.
c       On entry pgtol >= 0 is specified by the user.  The iteration
c         will stop when
c
c                 max{|proj g_i | i = 1, ..., n} <= pgtol
c
c         where pg_i is the ith component of the projected gradient.
c       The user can suppress this termination test by setting pgtol=0.
c
c     wa is a DOUBLE PRECISION  array of length 
c       (2mmax + 4)nmax + 12mmax^2 + 12mmax used as workspace.
c       This array must not be altered by the user.
c
c     iwa is an INTEGER  array of length 3nmax used as
c       workspace. This array must not be altered by the user.
c
c     task is a CHARACTER string of length 60.
c       On first entry, it must be set to 'START'.
c       On a return with task(1:2)='FG', the user must evaluate the
c         function f and gradient g at the returned value of x.
c       On a return with task(1:5)='NEW_X', an iteration of the
c         algorithm has concluded, and f and g contain f(x) and g(x)
c         respectively.  The user can decide whether to continue or stop
c         the iteration. 
c       When
c         task(1:4)='CONV', the termination test in L-BFGS-B has been 
c           satisfied;
c         task(1:4)='ABNO', the routine has terminated abnormally
c           without being able to satisfy the termination conditions,
c           x contains the best approximation found,
c           f and g contain f(x) and g(x) respectively;
c         task(1:5)='ERROR', the routine has detected an error in the
c           input parameters;
c       On exit with task = 'CONV', 'ABNO' or 'ERROR', the variable task
c         contains additional information that the user can print.
c       This array should not be altered unless the user wants to
c          stop the run for some reason.  See driver2 or driver3
c          for a detailed explanation on how to stop the run 
c          by assigning task(1:4)='STOP' in the driver.
c
c     iprint is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     csave  is a CHARACTER working array of length 60.
c
c     lsave is a LOGICAL working array of dimension 4.
c       On exit with task = 'NEW_X', the following information is
c         available:
c       lsave(1) = .true.  the initial x did not satisfy the bounds;
c       lsave(2) = .true.  the problem contains bounds;
c       lsave(3) = .true.  each variable has upper and lower bounds.
c
c     isave is an INTEGER working array of dimension 44.
c       On exit with task = 'NEW_X', it contains information that
c       the user may want to access:
c         isave(30) = the current iteration number;
c         isave(34) = the total number of function and gradient
c                         evaluations;
c         isave(36) = the number of function value or gradient
c                                  evaluations in the current iteration;
c         isave(38) = the number of free variables in the current
c                         iteration;
c         isave(39) = the number of active constraints at the current
c                         iteration;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in isave
c
c     dsave is a DOUBLE PRECISION working array of dimension 29.
c       On exit with task = 'NEW_X', it contains information that
c         the user may want to access:
c         dsave(2) = the value of f at the previous iteration;
c         dsave(5) = the machine precision epsmch generated by the code;
c         dsave(13) = the infinity norm of the projected gradient;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in dsave
c
c     --------------------------------------------------------------
c           END OF THE DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------
c
c     << An example of subroutine 'timer' for AIX Version 3.2 >>
c
c     subroutine timer(ttime)
c     double precision ttime
c     integer itemp, integer mclock
c     itemp = mclock()
c     ttime = dble(itemp)*1.0d-2
c     return
c     end

*********************************************************************
 */
double l_bfgs_b_wrapper(integer n, integer m, doublereal x[], doublereal l[],
                        doublereal u[], integer nbd[],
                        double (*fp) (const double x[], va_list argv),
                        void (*gp) (const double x[], double res[], va_list argv),
                        doublereal factr, doublereal pgtol, integer iprint,
                        int maxfun, int maxiter, int argc, ...) {
    va_list argv, tmp;
    char task[SIXTY], csave[SIXTY];
    logical lsave[4];
    integer iwa[3*n], isave[44];
    doublereal dsave[29], wa[2*m*n+5*n+11*m*m+8*m];
    va_start(argv, argc);

    strcpy(task,"START");
    int i, count = 0, funcount = 0, maxls = 20;
    for(i=5;i<SIXTY;i++) {
        task[i]=' ';
    }

    doublereal f = 0, g[n];
    for (i = 0; i < n; ++i) {
        g[i] = 0;
    }

    // TODO # The minimization routine wants f and g at the current x.
    //      # Note that interruptions due to maxfun are postponed
    //      # until the completion of the current minimization iteration.
    //      # Overwrite f and g:
    do {
        setulb_(&n,&m,x,l,u,nbd,&f,g,&factr,&pgtol,wa,iwa,task,&iprint,
                csave,lsave,isave,dsave, &maxls);

        if (strncmp(task,"FG",2)==0) {

            va_copy(tmp, argv);
            f = fp(x, tmp);
            va_copy(tmp, argv);
            gp(x, g, tmp);
            ++funcount;

        } else if(strncmp(task,"NEW_X",5)==0){
            ++count;
        } else {
            break;
        }
    } while(count < maxiter && funcount < maxfun);

    f = fp(x, argv);
    va_end(argv);
    va_end(tmp);

    return f;
}
