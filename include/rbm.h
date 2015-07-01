#ifndef RBM_H
#define RBM_H

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

int rbm_init( rbm_t *, long, long, double );
int rbm_free( rbm_t * );
int rbm_set_vis( rbm_t *, char * );
int rbm_set_hid( rbm_t *, char * );
int rbm_set_temp( rbm_t *, double );
int rbm_save( rbm_t *, char * );
int rbm_load( rbm_t *, char * );
double rbm_energy( rbm_t * );
void rbm_state( rbm_t *, int );
double rbm_energy_diff( rbm_t *, long );
int rbm_run_step( rbm_t *, long );
int rbm_update_hid( rbm_t * );
int rbm_update_vis( rbm_t * );
void rbm_zero_stats( rbm_t *, double *, double *, double * );
void rbm_update_weights( rbm_t *, double *, double *, double *, double );
int rbm_cd_mc( rbm_t *, char *, long, long, double *, double *, double * );
int dbn_init( dbn_t *, long, long *, double );
int dbn_free( dbn_t * );
void dbn_state( dbn_t * );
int dbn_save( dbn_t *, char * );
int dbn_save_init( dbn_t *, char * );
int dbn_load( dbn_t *, char * );
int dbn_load_init( dbn_t *, char * );
int dbn_cd_mc( dbn_t *, long, char *, long, long, long, long, double *, double *, double *, char *, double );
int dbn_assemble( dbn_t * );
int dbn_infer( dbn_t *, char * );

#endif

