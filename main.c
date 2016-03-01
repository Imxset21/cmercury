#include <stdio.h>
#include <stdlib.h>

#include "orbel.h"

int main(void)
{
    printf("orbel_zget(%f): %f\n", 1.0, orbel_zget(1.0));
    return EXIT_SUCCESS;
}
