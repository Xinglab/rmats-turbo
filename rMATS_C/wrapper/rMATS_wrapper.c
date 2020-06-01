#include "time.h"
#include "math.h"
#include "stdio.h"
#include "string.h"
#include "../include/rMATS_wrapper.h"
#include "../include/type.h"
#include "../include/util.h"
#include "../include/myfunc.h"
#include "../include/global.h"


clock_t dur = 0;
double cutoff = 0.1;
double rho = 0.9;


int likelihood_test_wrapper(char* inputf, double ctof, int nthread, double* output, int* len) {
    char *title_element_list[100];
    int row_num = 0, i = 0;
    cutoff = ctof;

    printf("number of thread=%d; input file=%s; cutoff=%g;\n",\
           nthread, inputf, cutoff);

    diff_list_node *list = (diff_list_node*)malloc(sizeof(diff_list_node));
    list->data = NULL, list->next = NULL, list->end = list;
    row_num = parse_file(inputf, list, title_element_list);

    clock_t start, diff;
    start = clock();

    odiff *datum[row_num];
    diff_list_node *node = list;
    double **prob = (double**)malloc(sizeof(double*)*row_num);
    for (i = 0; i < row_num; ++i) {
        node = node->next;
        datum[i] = node->data;
    }
    mp_threadpool(nthread, row_num, thread_wrapper_for_LT, (void**)datum, (void**)prob);

    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Total CPU time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("CPU time per thread taken %d seconds %d milliseconds\n", msec/1000/nthread, (msec%1000)/nthread);

    msec = dur * 1000 / CLOCKS_PER_SEC;
    printf("Time for func(single thread): %d seconds %d milliseconds\n", msec/1000, msec%1000);

    for (i = 0; i < row_num; ++i) {
        output[i] = *prob[i];
    }
    *len = row_num;

    return 0;
}
