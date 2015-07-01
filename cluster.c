#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int cluster_centroid( int k_in, int n_in, int dim_in, double *pts_in, int maxit_in, int *cmemb_out )
{
	int i,j,k,m,min,*size;
	double tmp,dist,mdist,*mean;

	/* Allocate space for mean calculations */
	size = (int*) malloc( k_in * sizeof(int) );
	mean = (double*) malloc( dim_in * k_in * sizeof(double) );

	/* Randomly assign points to clusters to begin */
	srand( (unsigned) time( 0 ) );
	for(i=0;i<n_in;i++)
		cmemb_out[i] = rand() % k_in; /* This is stupid but will work for now */

	/* Start iteration */
	for(i=0;i<maxit_in;i++)
	{
		/* Calculate current cluster means */
		for(j=0;j<dim_in*k_in;j++)
			mean[j] = 0.0, size[j] = 0;
		for(j=0;j<n_in;j++)
		{
			for(k=0;k<dim_in;k++)
				mean[dim_in*cmemb_out[j]+k] += pts_in[j*dim_in+k];
			size[cmemb_out[j]]++;
		}
		for(j=0;j<k_in;j++)
			for(k=0;k<dim_in;k++)
				mean[dim_in*j+k] /= (double) size[j];

		/* Assign each point to cluster closest to it */
		for(j=0;j<n_in;j++)
		{
			/* Find the cluster whose mean is closest to point j */
			for(k=0;k<k_in;k++)
			{
				/* Calculate distance between point n and cluster k */
				dist = 0.0;
				for(m=0;m<dim_in;m++)
					tmp = pts_in[dim_in*j+m] - mean[dim_in*k+m], dist += tmp * tmp;
				dist = sqrt( dist );
				if( k == 0 )
				{
					min = 0;
					mdist = dist;
				}
				else if( dist < mdist )
				{
					min = k;
					mdist = dist;
				}
			}

			/* Now assign point j to cluster min */
			cmemb_out[j] = min;
		}
	}

	return 0;
}

#define CLUSTER_INCREMENT 128
#define CLUSTER_MEMBER_INCREMENT 512
#undef CLUSTER_USE_MATRIX

static void distance_matrix( int n_in, int dim_in, double *pts_in, double **mat_out )
{
	int i,j,k;
	double tmp,dist;

	/* Allocate the space as part of the function */
	*mat_out = (double*) malloc( n_in * n_in * sizeof(double) );
	if( *mat_out == NULL )
	{
		fprintf( stderr, "Could not allocate enough memory. Exiting.\n" );
		exit(0);
	}

	/* Calculate pair-wise distances */
	for(i=0;i<n_in;i++)
	{
		/* Only calculate (i,j) = (j,i) once */
		for(j=i+1;j<n_in;j++)
		{
			dist = 0.0;
			for(k=0;k<dim_in;k++)
				tmp = pts_in[i*dim_in+k] - pts_in[j*dim_in+k], dist += tmp * tmp;
			dist = sqrt( dist );
			(*mat_out)[i*n_in+j] = dist;
			(*mat_out)[j*n_in+i] = dist;
		}
	}
}

/**
 * Be careful. The number of clusters changes and the number of members in
 * each cluster can also change. Make sure memory is being managed well.
 * @param n_in Number of input points to cluster
 * @param dim_in Dimension of the points to cluster
 * @param pts_in The list of n_in points of dimension dim_in
 * @param ep_in Epsilon distance to use to search point density
 * @param minpts_in Minimum number of points required to form a cluster
 * @param cs_out Cluster size vector output
 * @param clst_out Cluster list, a pass-by-reference of a int **
 * @return Returns the number of clusters found
 */
int cluster_dbscan( int n_in, int dim_in, double *pts_in, double ep_in, int minpts_in, int **cs_out, int ***clst_out )
{
	/* Variables nc = number of clusters, cs_out = cluster sizes */
	int i,j,k,m,nc,nca,nnb,nba,nnbp,nbap,*nb,*nbp,*csa; /* Variable nb is for storing neighbors at each iteration */
	int *vtd; /* Variables vtd[i] = -3 if pts_in[i] not visited, vtd[i] = -2 if visited, vtd[i] = -1 if noise, vtd[i] >= 0 if i belongs to cluster vtd[i] */
	double tmp,dist,*mat;

	/* Allocate space to store if each vertex is visited or not visited */
	vtd = (int*) malloc( n_in * sizeof(int) );
	for(i=0;i<n_in;i++)
		vtd[i] = -3; /* Initialize not visited */

	/* Allocate space for temporary neighbor storage */
	nba = CLUSTER_MEMBER_INCREMENT;
	nb = (int*) malloc( CLUSTER_MEMBER_INCREMENT * sizeof(int) );
	nbap = CLUSTER_MEMBER_INCREMENT;
	nbp = (int*) malloc( CLUSTER_MEMBER_INCREMENT * sizeof(int) );

	/* Allocate base memory for outputting the clusters and cluster sizes */
	nca = CLUSTER_INCREMENT; /* Initialize an allocated size for cluster data */
	csa = (int*) malloc( CLUSTER_INCREMENT * sizeof(int) ); /* Variable csa stores the amount of space allocated for each cluster individually */
	*cs_out = (int*) malloc( CLUSTER_INCREMENT * sizeof(int) );
	*clst_out = (int**) malloc( CLUSTER_INCREMENT * sizeof(int*) );

#ifdef CLUSTER_USE_MATRIX
	distance_matrix( n_in, dim_in, pts_in, &mat );
#endif

	/* Iterate through all unvisited nodes and their neighbors and mark them visited */
	nc = 0; /* Start with zero clusters */
	for(i=0;i<n_in;i++)
	{
		/* If point i has been visited then continue */
		if( vtd[i] != -3 )
			continue;
		vtd[i] = -2; /* Mark as visited but not necessarily noise */

		/* Get all epsilon-neighbors of pts[i*dim] */
		nnb = 0; /* Initialize number of neighbors found to zero */
		for(j=0;j<n_in;j++)
		{
			if( j == i )
				continue;
#ifdef CLUSTER_USE_MATRIX
			dist = mat[i*n_in+j];
#else
			dist = 0.0;
			for(k=0;k<dim_in;k++)
				tmp = pts_in[i*dim_in+k] - pts_in[j*dim_in+k], dist += tmp * tmp;
			dist = sqrt( dist );
#endif

			/* Tack it on to the end of nb for temporary keeping */
			if( dist < ep_in )
			{
				/* Allocate more neighbor space if not enough already */
				if( nnb + 1 > nba )
				{
					nba += CLUSTER_MEMBER_INCREMENT;
					nb = (int*) realloc( nb, nba * sizeof(int) );
					if( nb == NULL )
						return -1;
				}

				nb[nnb++] = j;
			}
		}

		/* Now have all neighbors counted; figure out what to do */
		if( nnb < minpts_in )
			vtd[i] = -1; /* Mark as noise and move on */
		else /* Expand to create a new cluster */
		{
			/* Create a new cluster */
			if( nc + 1 > nca )
			{
				nca += CLUSTER_INCREMENT; /* Increment nca to show that there is more space allocated */
				csa = (int*) realloc( csa, nca * sizeof(int) ); /* Allocate more space for storing individual allocated cluster sizes */
				*cs_out = (int*) realloc( *cs_out, nca * sizeof(int) ); /* Allocate more space for actual cluster sizes */
				*clst_out = (int**) realloc( *clst_out, nca * sizeof(int*) ); /* Space for the actual cluster members */
			}

			/* Add onto the end at position nc */
			csa[nc] = CLUSTER_MEMBER_INCREMENT;
			(*cs_out)[nc] = 1;
			(*clst_out)[nc] = (int*) malloc( CLUSTER_MEMBER_INCREMENT * sizeof(int) );
			(*clst_out)[nc][0] = i;
			vtd[i] = nc; /* Save the index of the cluster to which it belongs */
			for(j=0;j<nnb;j++)
			{
				/* Do epsilon-density for each neighbor nb[j] */
				if( vtd[nb[j]] == -3 )
				{
					vtd[nb[j]] = -2; /* Mark as visited but not necessarily noise */

					nnbp = 0;
					for(k=0;k<n_in;k++)
					{
						if( k == nb[j] )
							continue;
#ifdef CLUSTER_USE_MATRIX
						dist = mat[nb[j]*n_in+k];
#else
						dist = 0.0;
						for(m=0;m<dim_in;m++)
							tmp = pts_in[nb[j]*dim_in+m] - pts_in[k*dim_in+m], dist += tmp * tmp;
						dist = sqrt( dist );
#endif

						if( dist < ep_in )
						{
							/* Allocate more space if necessary */
							if( nnbp + 1 > nbap )
							{
								nbap += CLUSTER_MEMBER_INCREMENT;
								nbp = (int*) realloc( nbp, nbap * sizeof(int) );
								if( nbp == NULL )
									return -1;
							}

							/* Add it to the neighbors prime list */
							nbp[nnbp++] = k;
						}
					}
					if( nnbp >= minpts_in )
					{
						/* Combine nb with nbp */
						if( nnb + nnbp + 1 > nba )
						{
							nba += CLUSTER_MEMBER_INCREMENT;
							nb = (int*) realloc( nb, nba * sizeof(int) );
							if( nb == NULL )
								return -1;
						}

						/* Combine these two neighbor groups */
						for(k=0;k<nnbp;k++)
							nb[nnb+k] = nbp[k];
						nnb += nnbp; /* Important: Note that because nnb can increase, the loop expands! */
					}
				}

				/* If nb[j] is not yet a part of a cluster, then add it to cluster nc */
				if( vtd[nb[j]] < 0 ) /* Variable vtd[nb[j]] < 0 means taht nb[j] vertex is not part of a cluster yet */
				{
					/* Then add it to the new cluster just created */
					if( (*cs_out)[nc] + 1 > csa[nc] )
					{
						csa[nc] += CLUSTER_MEMBER_INCREMENT;
						(*clst_out)[nc] = (int*) realloc( (*clst_out)[nc], csa[nc] * sizeof(int) );
						if( (*clst_out)[nc] == NULL )
							return -1;
					}
					(*clst_out)[nc][(*cs_out)[nc]++] = nb[j];
					vtd[nb[j]] = nc;
				}
			}

			/* Finally increment after all additions made to the output cluster data structures */
			++nc;
		}
	}

	/* Clean up here; this will leak memory if you don't */
	free( vtd );
	free( csa );
	free( nb );
	free( nbp );

	return nc;
}

/**
 * Be careful. The number of clusters changes and the number of members in
 * each cluster can also change. Make sure memory is being managed well.
 * @param n_in Number of input points to cluster
 * @param dim_in Dimension of the points to cluster
 * @param pts_in The list of n_in points of dimension dim_in
 * @param ep_in Epsilon distance to use to search point density
 * @param minpts_in Minimum number of points required to form a cluster
 * @param cs_out Cluster size vector output
 * @param clst_out Cluster list, a pass-by-reference of a int **
 * @return Returns the number of clusters found
 */
int cluster_dbscan_density( int n_in, int dim_in, double *pts_in, double (*ep_in)(double*), int minpts_in, int **cs_out, int ***clst_out )
{
	/* Variables nc = number of clusters, cs_out = cluster sizes */
	int i,j,k,m,nc,nca,nnb,nba,nnbp,nbap,*nb,*nbp,*csa; /* Variable nb is for storing neighbors at each iteration */
	int *vtd; /* Variables vtd[i] = -3 if pts_in[i] not visited, vtd[i] = -2 if visited, vtd[i] = -1 if noise, vtd[i] >= 0 if i belongs to cluster vtd[i] */
	double tmp,dist,*mat;

	/* Allocate space to store if each vertex is visited or not visited */
	vtd = (int*) malloc( n_in * sizeof(int) );
	for(i=0;i<n_in;i++)
		vtd[i] = -3; /* Initialize not visited */

	/* Allocate space for temporary neighbor storage */
	nba = CLUSTER_MEMBER_INCREMENT;
	nb = (int*) malloc( CLUSTER_MEMBER_INCREMENT * sizeof(int) );
	nbap = CLUSTER_MEMBER_INCREMENT;
	nbp = (int*) malloc( CLUSTER_MEMBER_INCREMENT * sizeof(int) );

	/* Allocate base memory for outputting the clusters and cluster sizes */
	nca = CLUSTER_INCREMENT; /* Initialize an allocated size for cluster data */
	csa = (int*) malloc( CLUSTER_INCREMENT * sizeof(int) ); /* Variable csa stores the amount of space allocated for each cluster individually */
	*cs_out = (int*) malloc( CLUSTER_INCREMENT * sizeof(int) );
	*clst_out = (int**) malloc( CLUSTER_INCREMENT * sizeof(int*) );

#ifdef CLUSTER_USE_MATRIX
	distance_matrix( n_in, dim_in, pts_in, &mat );
#endif

	/* Iterate through all unvisited nodes and their neighbors and mark them visited */
	nc = 0; /* Start with zero clusters */
	for(i=0;i<n_in;i++)
	{
		/* If point i has been visited then continue */
		if( vtd[i] != -3 )
			continue;
		vtd[i] = -2; /* Mark as visited but not necessarily noise */

		/* Get all epsilon-neighbors of pts[i*dim] */
		nnb = 0; /* Initialize number of neighbors found to zero */
		for(j=0;j<n_in;j++)
		{
			if( j == i )
				continue;
#ifdef CLUSTER_USE_MATRIX
			dist = mat[i*n_in+j];
#else
			dist = 0.0;
			for(k=0;k<dim_in;k++)
				tmp = pts_in[i*dim_in+k] - pts_in[j*dim_in+k], dist += tmp * tmp;
			dist = sqrt( dist );
#endif

			/* Tack it on to the end of nb for temporary keeping */
			if( dist < ep_in( pts_in + i * dim_in ) )
			{
				/* Allocate more neighbor space if not enough already */
				if( nnb + 1 > nba )
				{
					nba += CLUSTER_MEMBER_INCREMENT;
					nb = (int*) realloc( nb, nba * sizeof(int) );
					if( nb == NULL )
						return -1;
				}

				nb[nnb++] = j;
			}
		}

		/* Now have all neighbors counted; figure out what to do */
		if( nnb < minpts_in )
			vtd[i] = -1; /* Mark as noise and move on */
		else /* Expand to create a new cluster */
		{
			/* Create a new cluster */
			if( nc + 1 > nca )
			{
				nca += CLUSTER_INCREMENT; /* Increment nca to show that there is more space allocated */
				csa = (int*) realloc( csa, nca * sizeof(int) ); /* Allocate more space for storing individual allocated cluster sizes */
				*cs_out = (int*) realloc( *cs_out, nca * sizeof(int) ); /* Allocate more space for actual cluster sizes */
				*clst_out = (int**) realloc( *clst_out, nca * sizeof(int*) ); /* Space for the actual cluster members */
			}

			/* Add onto the end at position nc */
			csa[nc] = CLUSTER_MEMBER_INCREMENT;
			(*cs_out)[nc] = 1;
			(*clst_out)[nc] = (int*) malloc( CLUSTER_MEMBER_INCREMENT * sizeof(int) );
			(*clst_out)[nc][0] = i;
			vtd[i] = nc; /* Save the index of the cluster to which it belongs */
			for(j=0;j<nnb;j++)
			{
				/* Do epsilon-density for each neighbor nb[j] */
				if( vtd[nb[j]] == -3 )
				{
					vtd[nb[j]] = -2; /* Mark as visited but not necessarily noise */

					nnbp = 0;
					for(k=0;k<n_in;k++)
					{
						if( k == nb[j] )
							continue;
#ifdef CLUSTER_USE_MATRIX
						dist = mat[nb[j]*n_in+k];
#else
						dist = 0.0;
						for(m=0;m<dim_in;m++)
							tmp = pts_in[nb[j]*dim_in+m] - pts_in[k*dim_in+m], dist += tmp * tmp;
						dist = sqrt( dist );
#endif

						if( dist < ep_in( pts_in + j * dim_in ) )
						{
							/* Allocate more space if necessary */
							if( nnbp + 1 > nbap )
							{
								nbap += CLUSTER_MEMBER_INCREMENT;
								nbp = (int*) realloc( nbp, nbap * sizeof(int) );
								if( nbp == NULL )
									return -1;
							}

							/* Add it to the neighbors prime list */
							nbp[nnbp++] = k;
						}
					}
					if( nnbp >= minpts_in )
					{
						/* Combine nb with nbp */
						if( nnb + nnbp + 1 > nba )
						{
							nba += CLUSTER_MEMBER_INCREMENT;
							nb = (int*) realloc( nb, nba * sizeof(int) );
							if( nb == NULL )
								return -1;
						}

						/* Combine these two neighbor groups */
						for(k=0;k<nnbp;k++)
							nb[nnb+k] = nbp[k];
						nnb += nnbp; /* Important: Note that because nnb can increase, the loop expands! */
					}
				}

				/* If nb[j] is not yet a part of a cluster, then add it to cluster nc */
				if( vtd[nb[j]] < 0 ) /* Variable vtd[nb[j]] < 0 means taht nb[j] vertex is not part of a cluster yet */
				{
					/* Then add it to the new cluster just created */
					if( (*cs_out)[nc] + 1 > csa[nc] )
					{
						csa[nc] += CLUSTER_MEMBER_INCREMENT;
						(*clst_out)[nc] = (int*) realloc( (*clst_out)[nc], csa[nc] * sizeof(int) );
						if( (*clst_out)[nc] == NULL )
							return -1;
					}
					(*clst_out)[nc][(*cs_out)[nc]++] = nb[j];
					vtd[nb[j]] = nc;
				}
			}

			/* Finally increment after all additions made to the output cluster data structures */
			++nc;
		}
	}

	/* Clean up here; this will leak memory if you don't */
	free( vtd );
	free( csa );
	free( nb );
	free( nbp );

	return nc;
}

int cluster_optics(  )
{
	
}

int cluster_fuzzy(  )
{
	
}

