#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <time.h>

#include "orbel.h"

extern double orbel_zget_(double *q);

static int test_orbel_zget(FILE* logfile, unsigned int num_vals)
{
    fputs("Testing orbel_zget() function:\n", logfile);
    double orbel_in = -100.0;
    int num_failures = 0;
    
    while (--num_vals > 0)
    {
        orbel_in = rand() > (RAND_MAX / 2) ? rand() : -rand();
        const double cresult = corbel_zget(orbel_in);
        const double fresult = orbel_zget_(&orbel_in);
        
#ifdef DEBUG
        fprintf(logfile, "\tcorbel_zget(%f): %f\n", orbel_in, cresult);
        fprintf(logfile, "\torbel_zget(%f): %f\n", orbel_in, fresult);
#endif

        if (fabs(cresult - fresult) > 0.1)
        {
            fprintf(logfile,
                    "FAILED: test_orbel_zget(q = %f) - expected %f, got %f\n",
                    orbel_in, fresult, cresult);
            num_failures++;
        } else {
            fprintf(logfile,
                    "SUCCEESS: test_orbel_zget(q = %f) - both are %f\n",
                    orbel_in, fresult);
        }

    }

    fprintf(logfile, "orbel_zget() failures: %i\n", num_failures);

    return num_failures;
}

int main(int argc, char** argv)
{
    if (argc > 2)
    {
        fputs("main.bin usage: ./main.bin [logfile]", stderr);
        exit(EXIT_FAILURE);
    }
    
#ifndef RANDOM_SEED
    srand(clock());
#else
    srand(RANDOM_SEED);
#endif
    
    FILE* fp = NULL;
    if (argc == 1)
    {
        fp = stderr;
    } else if (argc == 2) {
        fp = fopen(argv[1], "w+");
    }
    
    printf("%s\n", test_orbel_zget(stderr, 100) ? "Failure" : "Success");
    exit(EXIT_SUCCESS);
}
