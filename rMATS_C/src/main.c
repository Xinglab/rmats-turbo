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


int main(int argc, char *argv[]) {
    char *inputf = "input.txt", *outputf = NULL;
    char *str_line = NULL, *title_element_list[100];
    int nthread = 1, opt, row_num = 0, i = 0;
    size_t str_line_n = 0;

    while ((opt = getopt(argc, argv, "t:bc:o:i:")) != -1) {
        switch (opt) {
        case 't':
            nthread = atoi(optarg);
            break;
        case 'b':
            // ignore legacy batch_size parameter
            break;
        case 'c':
            cutoff = atof(optarg);
            break;
        case 'o':
            outputf = (char*)malloc(sizeof(char)*(strlen(optarg) + 1));
            strcpy(outputf, optarg);
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

    odiff **datum = (odiff**)malloc(row_num*sizeof(odiff*));
    if (datum == NULL) {
      fprintf(stderr, "Failed to allocate datum (%lu)\n", row_num*sizeof(odiff*));
      exit(EXIT_FAILURE);
    }
    diff_list_node *node = list;
    double **prob = (double**)malloc(row_num*sizeof(double*));
    if (prob == NULL) {
      fprintf(stderr, "Failed to allocate prob (%lu)\n", row_num*sizeof(double*));
      exit(EXIT_FAILURE);
    }
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
        printf("Fail to open %s!\n", inputf);
        return -1;
    }
    if ((ofp = fopen(outputf, "w")) == NULL) {
        printf("Fail to open %s!\n", outputf);
        return -1;
    }

    // original comment: analyze the title of the inputed data file to find the
    //                   information of how much simulation are involved
    if (getline(&str_line, &str_line_n, ifp) == -1) {
        printf("Failed to read first line of input file\n");
        return -1;
    }
    str_line[strlen(str_line)-1] = 0; // overwrite the old newline
    fprintf(ofp, "%s\tPValue\n", str_line); // append PValue header to old headers
    for (i = 0; i < row_num; ++i) {
        if (getline(&str_line, &str_line_n, ifp) == -1) {
            printf("Failed to read row %d of input file\n", i);
            return -1;
        }
        str_line[strlen(str_line)-1] = 0; // overwrite old newline
        fprintf(ofp, "%s\t%.12g\n", str_line, *prob[i]); // append PValue
    }
    fclose(ifp);
    fclose(ofp);

    msec = dur * 1000 / CLOCKS_PER_SEC;
    printf("Time for func(single thread): %d seconds %d milliseconds\n", msec/1000, msec%1000);

    free(str_line);
    free(outputf);
    free(datum);
    free(prob);
    return 0;
}
