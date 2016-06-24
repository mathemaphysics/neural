// vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab:smarttab
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "config.h"

#define RBMPROB( E, T ) 1.0 / ( 1.0 + exp( E / T ) )
#define WTS_STDEV 0.001
#define RBM_VERBOSE_LEVEL 0
#define DBN_SAMPLE_INTERVAL 200

typedef struct
{
    long nv;
    long nh;
    char *vis;
    char *hid;
    double *vth;
    double *hth;
    double *wts;
    double temp;
} rbm_t;

typedef struct
{
    long nlayer;
    rbm_t *layers;
    double temp;
} dbn_t;

int rbm_init( rbm_t *obj_in, long nv_in, long nh_in, double temp_in )
{
    long h,i;
    if( obj_in == NULL )
	return -1;
    obj_in->nv = nv_in;
    obj_in->nh = nh_in;
    obj_in->vis = (char*) malloc( nv_in * sizeof(char) );
    obj_in->hid = (char*) malloc( nh_in * sizeof(char) );
    obj_in->vth = (double*) malloc( nv_in * sizeof(double) );
    obj_in->hth = (double*) malloc( nh_in * sizeof(double) );
    obj_in->wts = (double*) malloc( nv_in * nh_in * sizeof(double) );
    obj_in->temp = temp_in; /* Default temperature = 1 */
    srand((unsigned)time(0));
    for(i=0;i<nv_in;i++)
	obj_in->vis[i] = (char) (rand() % 2); /* Should set this to input vector */
    for(i=0;i<nh_in;i++)
	obj_in->hid[i] = (char) (rand() % 2); /* Randomize initial hidden state */
    for(i=0;i<nv_in;i++)
	obj_in->vth[i] = 0.5; /* Initialize all visible to same threshold */
    for(i=0;i<nh_in;i++)
	obj_in->hth[i] = 0.5; /* Initialize all hidden to same threshold */
    for(i=0;i<nv_in*nh_in;i++)
	obj_in->wts[i] = 0.0;
    for(h=0;h<100;h++) /* Generate an initially Gaussian weight distribution with stdev given */
    {
	for(i=0;i<nv_in*nh_in;i++)
	{
	    /* Equal probabilities of going right and left */
	    if( (double) rand() / (double) RAND_MAX < 0.5 )
		obj_in->wts[i] += WTS_STDEV;
	    else
		obj_in->wts[i] -= WTS_STDEV;
	}
    }
    return 0;
}

int rbm_free( rbm_t *obj_in )
{
    if( obj_in == NULL )
	return -1;
    if( obj_in->vis != NULL )
        free( obj_in->vis );
    if( obj_in->hid != NULL )
        free( obj_in->hid );
    if( obj_in->vth != NULL )
        free( obj_in->vth );
    if( obj_in->hth != NULL )
        free( obj_in->hth );
    if( obj_in->wts != NULL )
        free( obj_in->wts );
    return 0;
}

int rbm_set_vis( rbm_t *obj_in, char *vis_in )
{
    long i;
    if( obj_in == NULL )
	return -1;
    for(i=0;i<obj_in->nv;i++)
	obj_in->vis[i] = vis_in[i];
    return 0;
}

int rbm_set_hid( rbm_t *obj_in, char *hid_in )
{
    long i;
    if( obj_in == NULL )
	return -1;
    for(i=0;i<obj_in->nh;i++)
	obj_in->hid[i] = hid_in[i];
    return 0;
}

int rbm_set_temp( rbm_t *obj_in, double temp_in )
{
    if( obj_in == NULL )
	return -1;
    obj_in->temp = temp_in;
    return 0;
}

int rbm_save( rbm_t *obj_in, char *fn_in )
{
    FILE *fp;
    if( obj_in == NULL )
	return -1;
    fp = fopen( fn_in, "wb" );
    if( fp == NULL )
        return -2;
    fwrite( &(obj_in->nv), sizeof(long), (size_t) 1, fp );
    fwrite( &(obj_in->nh), sizeof(long), (size_t) 1, fp );
    fwrite( obj_in->vth, sizeof(double), (size_t) obj_in->nv, fp );
    fwrite( obj_in->hth, sizeof(double), (size_t) obj_in->nh, fp );
    fwrite( obj_in->wts, sizeof(double), (size_t) (obj_in->nv*obj_in->nh), fp );
    fwrite( obj_in->vis, sizeof(char), (size_t) obj_in->nv, fp );
    fwrite( obj_in->hid, sizeof(char), (size_t) obj_in->nh, fp );
    fclose( fp );
    return 0;
}

int rbm_load( rbm_t *obj_in, char *fn_in )
{
    FILE *fp;
    if( obj_in == NULL )
	return -1;
    fp = fopen( fn_in, "rb" );
    if( fp == NULL )
        return -2;
    fread( &(obj_in->nv), sizeof(long), (size_t) 1, fp );
    fread( &(obj_in->nh), sizeof(long), (size_t) 1, fp );
    fread( obj_in->vth, sizeof(double), (size_t) obj_in->nv, fp );
    fread( obj_in->hth, sizeof(double), (size_t) obj_in->nh, fp );
    fread( obj_in->wts, sizeof(double), (size_t) (obj_in->nv*obj_in->nh), fp );
    fread( obj_in->vis, sizeof(char), (size_t) obj_in->nv, fp );
    fread( obj_in->hid, sizeof(char), (size_t) obj_in->nh, fp );
    fclose( fp );
    return 0;
}

double rbm_energy( rbm_t *obj_in )
{
    long i,j; /* Row = vis index, col = hid index */
    double sum = 0.0;
    for(i=0;i<obj_in->nv;i++)
	for(j=0;j<obj_in->nh;j++)
	    if( obj_in->vis[i] == 1 && obj_in->hid[j] == 1 )
		sum -= obj_in->wts[i*obj_in->nh+j];
    for(i=0;i<obj_in->nv;i++)
	if( obj_in->vis[i] == 1 )
	    sum -= obj_in->vth[i];
    for(i=0;i<obj_in->nh;i++)
	if( obj_in->hid[i] == 1 )
	    sum -= obj_in->hth[i];
    return sum;
}

void rbm_state( rbm_t *obj_in, int deg_in )
{
    long i,j;
    if( deg_in & 32 )
    {
	if( deg_in & 16 )
	    printf( "E = " );
	printf( "%10.5f", rbm_energy( obj_in ) );
	printf( "\n" );
    }
    if( deg_in & 1 )
    {
	if( deg_in & 16 )
	    printf( "V:\n" );
	for(i=0;i<obj_in->nv;i++)
	    printf( "%3d", obj_in->vis[i] );
	printf( "\n" );
    }
    if( deg_in & 2 )
    {
	if( deg_in & 16 )
	    printf( "H:\n" );
	for(i=0;i<obj_in->nh;i++)
	    printf( "%3d", obj_in->hid[i] );
	printf( "\n" );
    }
    if( deg_in & 4 )
    {
	if( deg_in & 16 )
	    printf( "W:\n" );
	for(i=0;i<obj_in->nv;i++)
	{
	    for(j=0;j<obj_in->nh;j++)
		printf( "%10.5f", obj_in->wts[i*obj_in->nh+j] );
	    printf( "\n" );
	}
    }
    if( deg_in & 8 )
    {
	if( deg_in & 16 )
	    printf( "VB:\n" );
	for(i=0;i<obj_in->nv;i++)
	    printf( "%10.5f", obj_in->vth[i] );
	printf( "\n" );
	if( deg_in & 16 )
	    printf( "HB:\n" );
	for(i=0;i<obj_in->nh;i++)
	    printf( "%10.5f", obj_in->hth[i] );
	printf( "\n" );
    }
}

/**
 * idx_in is an index in the range [0,nh+nv] which indexes
 * hidden and visible nodes with one index;
 * indexes 0 through nh-1 reference all of the hidden nodes
 * while indexes nh through nh+nv-1 references the visible nodes
*/
double rbm_energy_diff( rbm_t *obj_in, long idx_in )
{
    long i,idx;
    double sum = 0.0;
    if( idx_in < obj_in->nh ) /* If this is a hidden index */
    {
	idx = idx_in;
	for(i=0;i<obj_in->nv;i++)
	    if( obj_in->vis[i] == 1 )
		sum -= obj_in->wts[i*obj_in->nh+idx]; /* wts has visible as first index */
	sum -= obj_in->hth[idx];
    }
    else /* If this is a visible index */
    {
	idx = idx_in - obj_in->nh;
	for(i=0;i<obj_in->nh;i++)
	    if( obj_in->hid[i] == 1 )
		sum -= obj_in->wts[idx*obj_in->nh+i];
	sum -= obj_in->vth[idx];
    }
    return sum;
}

int rbm_run_step( rbm_t *obj_in, long idx_in )
{
    static double p,r,en;
    if( obj_in == NULL )
	return -1;
    if( idx_in > obj_in->nv + obj_in->nh || idx_in < 0 )
	return -2;
    if( idx_in < obj_in->nh ) /* Then update obj_in->hid[idx_in] */
    {
	en = rbm_energy_diff( obj_in, idx_in );
	p = RBMPROB( en, obj_in->temp );
	r = (double) rand() / (double) RAND_MAX;
	if( r < p )
	    obj_in->hid[idx_in] = 1;
	else
	    obj_in->hid[idx_in] = 0;
    }
    else
    {
	en = rbm_energy_diff( obj_in, idx_in );
	p = RBMPROB( en, obj_in->temp );
	r = (double) rand() / (double) RAND_MAX;
	if( r < p )
	    obj_in->vis[idx_in-obj_in->nh] = 1;
	else
	    obj_in->vis[idx_in-obj_in->nh] = 0;
    }
    return 0;
}

int rbm_update_hid( rbm_t *obj_in )
{
    long i;
    if( obj_in == NULL )
	return -1;
    for(i=0;i<obj_in->nh;i++)
	rbm_run_step( obj_in, i );
    return 0;
}

int rbm_update_vis( rbm_t *obj_in )
{
    long i;
    if( obj_in == NULL )
	return -1;
    for(i=0;i<obj_in->nv;i++)
	rbm_run_step( obj_in, obj_in->nh + i );
    return 0;
}

void rbm_zero_stats( rbm_t *obj_in, double *cor_out, double *vav_out, double *hav_out )
{
    long i;
    for(i=0;i<obj_in->nv*obj_in->nh;i++)
	cor_out[i] = 0.0;
    for(i=0;i<obj_in->nv;i++)
	vav_out[i] = 0.0;
    for(i=0;i<obj_in->nh;i++)
	hav_out[i] = 0.0;
}

void rbm_update_weights( rbm_t *obj_in, double *cor_in, double *vav_in, double *hav_in, double lr_in )
{
    long i,j;
    for(i=0;i<obj_in->nv;i++)
	for(j=0;j<obj_in->nh;j++)
	    obj_in->wts[i*obj_in->nh+j] += lr_in * cor_in[i*obj_in->nh+j];
    for(i=0;i<obj_in->nv;i++)
	obj_in->vth[i] += lr_in * vav_in[i];
    for(i=0;i<obj_in->nh;i++)
	obj_in->hth[i] += lr_in * hav_in[i];
}

int rbm_cd_mc( rbm_t *obj_in, char *vis_in, long nstep_in, long nmcs_in, double *cor_out, double *vav_out, double *hav_out )
{
    long i,j,k;

    /* Measure <vh> first by running with v = v0 fixed */
    if( vis_in != NULL )
	rbm_set_vis( obj_in, vis_in );
    for(i=0;i<nstep_in;i++) /* Remember vis is not changed within this loop */
    {
	/* Update for another sample */
	rbm_update_hid( obj_in ); /* Running through this loop does not change vis */
	
	/* Calculate correlations and averages */
	if( cor_out != NULL )
	    for(j=0;j<obj_in->nv;j++)
		for(k=0;k<obj_in->nh;k++)
		    if( obj_in->vis[j] == 1 && obj_in->hid[k] == 1 )
			cor_out[j*obj_in->nh+k] += 1.0 / (double) nstep_in;
	if( vav_out != NULL )
	    for(j=0;j<obj_in->nv;j++)
		if( obj_in->vis[j] == 1 )
		    vav_out[j] += 1.0 / (double) nstep_in;
	if( hav_out != NULL )
	    for(j=0;j<obj_in->nh;j++)
		if( obj_in->hid[j] == 1 )
		    hav_out[j] += 1.0 / (double) nstep_in;
    }
    for(i=0;i<nstep_in;i++) /* nstep_in = the number of independent Markov chains to run */
    {
	/* Run a Markov chain to get a visible vector reproduced via nmcs_in steps */
	if( vis_in != NULL ) /* But running this loop multiple times changes vis; need it to be reset */
	    rbm_set_vis( obj_in, vis_in );
	for(j=0;j<nmcs_in;j++) /* nmcs_in = the number of Markov chain steps to run */
	{
	    rbm_update_hid( obj_in ); /* Sample p(h|v) */
	    rbm_update_vis( obj_in ); /* Sample p(v|h) */
	}

	/* Update hid one more time */
	rbm_update_hid( obj_in );

	/* Calculate (the negative of) correlations and averages to finish calculating the learning signal */
	if( cor_out != NULL )
	    for(j=0;j<obj_in->nv;j++)
		for(k=0;k<obj_in->nh;k++)
		    if( obj_in->vis[j] == 1 && obj_in->hid[k] == 1 )
			cor_out[j*obj_in->nh+k] -= 1.0 / (double) nstep_in;
	if( vav_out != NULL )
	    for(j=0;j<obj_in->nv;j++)
		if( obj_in->vis[j] == 1 )
		    vav_out[j] -= 1.0 / (double) nstep_in;
	if( hav_out != NULL )
	    for(j=0;j<obj_in->nh;j++)
		if( obj_in->hid[j] == 1 )
		    hav_out[j] -= 1.0 / (double) nstep_in;
    }
    return 0;
}

int dbn_init( dbn_t *obj_in, long nlayer_in, long *lsizes_in, double temp_in )
{
    long i;
    obj_in->nlayer = nlayer_in;
    obj_in->layers = (rbm_t*) malloc( ( nlayer_in - 1 ) * sizeof(rbm_t) );
    for(i=0;i<nlayer_in-1;i++)
	rbm_init( &(obj_in->layers[i]), lsizes_in[i], lsizes_in[i+1], temp_in );
    return 0;
}

int dbn_free( dbn_t *obj_in )
{
    long i;
    for(i=obj_in->nlayer-2;i>=0;i--)
        rbm_free( &obj_in->layers[i] );
    return 0;
}

void dbn_state( dbn_t *obj_in )
{
    long i;
    for(i=0;i<obj_in->nlayer-1;i++)
	rbm_state( &(obj_in->layers[i]), 4|8|16 );
}

int dbn_save( dbn_t *obj_in, char *fn_in )
{
    long i;
    FILE *fp;
    if( obj_in == NULL )
	return -1;
    if( obj_in->layers == NULL )
	return -2;
    fp = fopen( fn_in, "wb" );
    if( fp == NULL )
	return -3;
    for(i=0;i<obj_in->nlayer-1;i++)
    {
	fwrite( &(obj_in->layers[i].nv), sizeof(long), (size_t) 1, fp );
	fwrite( &(obj_in->layers[i].nh), sizeof(long), (size_t) 1, fp );
	fwrite( obj_in->layers[i].vth, sizeof(double), (size_t) obj_in->layers[i].nv, fp );
	fwrite( obj_in->layers[i].hth, sizeof(double), (size_t) obj_in->layers[i].nh, fp );
	fwrite( obj_in->layers[i].wts, sizeof(double), (size_t) (obj_in->layers[i].nv*obj_in->layers[i].nh), fp );
    }
    fclose( fp );
    return 0;
}

int dbn_save_init( dbn_t *obj_in, char *fn_in )
{
    long i;
    FILE *fp;
    if( obj_in == NULL )
        return -1;
    if( obj_in->layers == NULL )
        return -2;
    fp = fopen( fn_in, "wb" );
    if( fp == NULL )
        return -3;
    fwrite( &(obj_in->nlayer), sizeof(long), (size_t) 1, fp );
    for(i=0;i<obj_in->nlayer-1;i++)
        fwrite( &(obj_in->layers[i].nv), sizeof(long), (size_t) 1, fp );
    fwrite( &(obj_in->layers[obj_in->nlayer-2].nh), sizeof(long), (size_t) 1, fp );
    for(i=0;i<obj_in->nlayer-1;i++)
    {
        fwrite( &(obj_in->layers[i].nv), sizeof(long), (size_t) 1, fp );
        fwrite( &(obj_in->layers[i].nh), sizeof(long), (size_t) 1, fp );
        fwrite( obj_in->layers[i].vth, sizeof(double), (size_t) obj_in->layers[i].nv, fp );
        fwrite( obj_in->layers[i].hth, sizeof(double), (size_t) obj_in->layers[i].nh, fp );
        fwrite( obj_in->layers[i].wts, sizeof(double), (size_t) (obj_in->layers[i].nv*obj_in->layers[i].nh), fp );
    }
    fclose( fp );
    return 0;
}

int dbn_load( dbn_t *obj_in, char *fn_in )
{
    long i;
    FILE *fp;
    if( obj_in == NULL )
        return -1;
    if( obj_in->layers == NULL )
        return -2;
    fp = fopen( fn_in, "rb" );
    if( fp == NULL )
        return -3;
    for(i=0;i<obj_in->nlayer-1;i++)
    {
        fread( &(obj_in->layers[i].nv), sizeof(long), (size_t) 1, fp );
        fread( &(obj_in->layers[i].nh), sizeof(long), (size_t) 1, fp );
        fread( obj_in->layers[i].vth, sizeof(double), (size_t) obj_in->layers[i].nv, fp );
        fread( obj_in->layers[i].hth, sizeof(double), (size_t) obj_in->layers[i].nh, fp );
        fread( obj_in->layers[i].wts, sizeof(double), (size_t) (obj_in->layers[i].nv*obj_in->layers[i].nh), fp );
    }
    fclose( fp );
    return 0;
}

int dbn_load_init( dbn_t *obj_in, char *fn_in )
{
    long i,nl,*ls;
    double temp;
    FILE *fp;
    if( obj_in == NULL )
        return -1;
    fp = fopen( fn_in, "rb" );
    if( fp == NULL )
        return -2;
    fread( &nl, sizeof(long), (size_t) 1, fp );
    if( nl <= 0 )
        return -3;
    ls = (long*) malloc( nl * sizeof(long) );
    for(i=0;i<nl;i++)
        fread( &ls[i], sizeof(long), (size_t) 1, fp );
    dbn_init( obj_in, nl, ls, temp );
    for(i=0;i<nl-1;i++)
    {
        fread( &(obj_in->layers[i].nv), sizeof(long), (size_t) 1, fp );
        fread( &(obj_in->layers[i].nh), sizeof(long), (size_t) 1, fp );
        fread( obj_in->layers[i].vth, sizeof(double), (size_t) obj_in->layers[i].nv, fp );
        fread( obj_in->layers[i].hth, sizeof(double), (size_t) obj_in->layers[i].nh, fp );
        fread( obj_in->layers[i].wts, sizeof(double), (size_t) (obj_in->layers[i].nv*obj_in->layers[i].nh), fp );
    }
    fclose( fp );
    return 0;
}

/* Sample total posterior at each layer by propagating input distribution of vectors forward through the system */
int dbn_cd_mc( dbn_t *obj_in, long ntrain_in, char *train_in, long nstep_in, long nequil_in, long nmcs_in, long nsamp_in, double *cor_in, double *vav_in, double *hav_in, char *samp_in, double lr_in )
{
    long i,j,k,m;
    double tmp;

    /* First train the top layer */
    for(i=0;i<nstep_in;i++)
    {
        for(j=0;j<ntrain_in;j++)
        {
            rbm_zero_stats( &(obj_in->layers[0]), cor_in, vav_in, hav_in );
            rbm_cd_mc( &(obj_in->layers[0]), train_in + j * obj_in->layers[0].nv, nequil_in, nmcs_in, cor_in, vav_in, hav_in );
            rbm_update_weights( &(obj_in->layers[0]), cor_in, vav_in, hav_in, lr_in );
        }
    }

    /* Now train the other layers */
    for(i=1;i<obj_in->nlayer-1;i++)
    {
        /* Take samples by propagating forward the input data vector distribution from the top */
        for(j=0;j<nsamp_in;j++) /* Take nsamp_in samples of each input data vector and learn them */
        {
            for(k=0;k<ntrain_in;k++)
            {
                rbm_set_vis( &obj_in->layers[0], train_in + k * obj_in->layers[0].nv );
                for(m=1;m<=i;m++)
                {
                    rbm_update_hid( &obj_in->layers[m-1] );
                    rbm_set_vis( &obj_in->layers[m], obj_in->layers[m-1].hid );
                }
                rbm_zero_stats( &obj_in->layers[i], cor_in, vav_in, hav_in );
                rbm_cd_mc( &obj_in->layers[i], obj_in->layers[i-1].hid, nequil_in, nmcs_in, cor_in, vav_in, hav_in );
                rbm_update_weights( &obj_in->layers[i], cor_in, vav_in, hav_in, lr_in );
            }
        }
    }
    return 0;
}

/* Connect the RBMs in obj_in->layers to one another */
int dbn_assemble( dbn_t *obj_in )
{
    long i;
    for(i=0;i<obj_in->nlayer-1;i++)
	obj_in->layers[i+1].vis = obj_in->layers[i].hid;
}

/* Sets the value of vis and propagates the activities up to the top-level associative memory */
int dbn_infer( dbn_t *obj_in, char *vis_in )
{
    long i;

    rbm_set_vis( &(obj_in->layers[0]), vis_in );
    for(i=0;i<obj_in->nlayer-1;i++)
	rbm_update_hid( &(obj_in->layers[i]) );
    return 0;
}

