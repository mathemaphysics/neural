// vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab:smarttab
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "config.h"

/**
 * Structure which contains all information
 * pertaining to every neuron in every layer
 * and their connection weights
 */
typedef struct
{
    /**
     * Number of layers
     */
    int nls;

    /**
     * Number of neurons in each layer
     */
    int *nns;

    /**
     * Neuron activation for each neuron
     * in each layer; intended to indicate
     * "states"
     */
    double **sts;

    /**
     * Neron biases for each neuron in each
     * layer; intended to indicate "biases"
     * via "bss"
     */
    double **bss;

    /**
     * Weight connecting each adjacent layer
     * of neurons to the next; intended to
     * indicate "weights"
     */
    double **wts;
} nnet_t;

/**
 * Initializes all data structures required
 * to run an artificial neural network of the
 * given input size
 * @param obj Neural network object to initialize
 * @param 
 */
int nnet_init( nnet_t *obj, int nls, int *nns )
{
    int i;

    /* Make sure good values encountered */
    if( nls < 0 )
        return -1;
    for(i=0;i<nls;i++)
    {
        if( nns[i] < 0 )
            return -2;
    }

    /* Now allocate all layers and neurons */
    obj->nls = nls;
    obj->nns = (int*) malloc( nls * sizeof(int) );
    if( obj->nns == NULL )
        return -3; /* Failure of malloc() */
    for(i=0;i<nls;i++)
    {
        obj->nns[i] = nns[i];
        obj->bss = (double**) malloc( nls * sizeof(double*) );
        if( obj->bss == NULL )
            return -4;
        obj->sts[i] = (double*) malloc( nns[i] * sizeof(double) );
        // Check malloc() worked
        obj->bss[i] = (double*) malloc( nns[i] * sizeof(double) );
        // Check malloc() worked
    }

    return 0;
}

/**
 * Clean up data and functions contained in
 * the neural network object
 * @param obj Neural network object of interest
 */
int nnet_free( nnet_t *obj )
{
    
}

