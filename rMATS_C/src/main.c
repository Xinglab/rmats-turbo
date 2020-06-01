#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include "../include/type.h"
#include "../include/util.h"
#include "../include/myfunc.h"
#include "../include/global.h"


double rho = 0.9;
double cutoff = 0.0001;
clock_t dur = 0;
// char *optarg;
// int optind = 0;


int main(int argc, char *argv[]) {
    char *inputf = "input.txt", *outputf = (char*)malloc(sizeof(char)*MAX_CHAR);
    char str_line[MAX_LINE], *title_element_list[100];
    int nthread = 1, batch_size = 1, opt, row_num = 0, i = 0;

    while ((opt = getopt(argc, argv, "t:bc:o:i:")) != -1) {
        switch (opt) {
        case 't':
            nthread = atoi(optarg);
            break;
        case 'b':
            batch_size = 0;
            break;
        case 'c':
            cutoff = atof(optarg);
            break;
        case 'o':
            // strcpy(outputf, optarg);
            // strcat(outputf, "/rMATS_Result_P.txt");
            outputf = optarg;
            break;
        case 'i':
            inputf = optarg;
            break;
        default: /* '?' */
            fprintf(stderr, "Usage: %s [-t nthread] [-c] cutoff [-o] output folder [-i] input file\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    printf("number of thread=%d; input file=%s; output folder=%s; cutoff=%g;\n",\
           nthread, inputf, outputf, cutoff);

    if (optind > argc) {
        fprintf(stderr, "Expected argument after options\n");
        exit(EXIT_FAILURE);
    }

    if (nthread > MAXTHREAD) {
        printf("nthread (%d) should be less than MAXTHREAD (%d)\n", nthread, MAXTHREAD);
        exit(0);
    }

    diff_list_node *list = (diff_list_node*)malloc(sizeof(diff_list_node));
    list->data = NULL, list->next = NULL, list->end = list;
    row_num = parse_file(inputf, list, title_element_list);

    clock_t start, diff;
    start = clock();

    odiff *datum[row_num];
    diff_list_node *node = list;
    double *prob[row_num];
    for (i = 0; i < row_num; ++i) {
        node = node->next;
        datum[i] = node->data;
    }
    mp_threadpool(nthread, row_num, thread_wrapper_for_LT, (void**)datum, (void**)prob);

    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Total Wallclock time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("Wallclock time per thread taken %d seconds %d milliseconds\n", msec/1000/nthread, (msec%1000)/nthread);

    FILE *ifp = NULL, *ofp = NULL;
    if ((ifp = fopen(inputf, "r")) == NULL) {
        printf("Fail to open %s!", inputf);
        return -1;
    }
    if ((ofp = fopen(outputf, "w")) == NULL) {
        printf("Fail to open %s!", outputf);
        return -1;
    }

    // original comment: analyze the title of the inputed data file to find the
    //                   information of how much simulation are involved
    fgets(str_line, MAX_LINE, ifp);
    str_line[strlen(str_line)-1] = 0;
    fprintf(ofp, "%s\tPValue\n", str_line);
    for (i = 0; i < row_num; ++i) {
        fgets(str_line, MAX_LINE, ifp);
        str_line[strlen(str_line)-1] = 0;
        fprintf(ofp, "%s\t%.12g\n", str_line, *prob[i]);
    }
    fclose(ifp);
    fclose(ofp);

    msec = dur * 1000 / CLOCKS_PER_SEC;
    printf("Time for func(single thread): %d seconds %d milliseconds\n", msec/1000, msec%1000);

    return 0;
}
