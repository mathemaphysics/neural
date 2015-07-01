#ifndef CLUSTER_H
#define CLUSTER_H

int cluster_centroid( int, int, int, double *, int, int * );

int cluster_dbscan( int, int, double *, double, int, int **, int *** );

int cluster_dbscan_density( int, int, double *, double (*)(double*), int, int **, int *** );

int cluster_optics(  );

int cluster_fuzzy(  );

#endif

