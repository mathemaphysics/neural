// vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab:smarttab
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "config.h"

typedef struct
{
    int nl;
    int nn;
    double *sts;
    double *bss;
    double **wts;
} nnet_t;

