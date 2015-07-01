#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ncurses.h>
#include "config.h"

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

void dbn_print_curses( dbn_t *obj_in )
{
	int i,j;
	for(i=0;i<obj_in->nlayer-1;i++)
		for(j=0;j<obj_in->layers[i].nv;j++)
			mvprintw( i, j, "%d", (int) obj_in->layers[i].vis[j] );
	for(j=0;j<obj_in->layers[obj_in->nlayer-2].nh;j++)
		mvprintw( obj_in->nlayer - 1, j, "%d", (int) obj_in->layers[obj_in->nlayer-2].hid[j] );
	refresh();
}

/**
 * This function makes use of the ncurses library to display the
 * neural network in real time as  it is updated. It also allows
 * manual control of time stepping.
 */
void dbn_run( dbn_t *obj_in )
{
	/* Assume that the dbn has been properly set up already */
	long i,j;
	char ch;

	/* Initialize ncurses */
	initscr();
	raw();
	keypad( stdscr, TRUE );
	noecho();

	for(;;)
	{
		/* Print the output first, then update */
		dbn_print_curses( obj_in );

		/* Do the updating here */
		for(i=0;i<8;i++)
		{
			rbm_update_hid( &(obj_in->layers[obj_in->nlayer-2]) );
			rbm_update_vis( &(obj_in->layers[obj_in->nlayer-2]) );
		}
		for(i=obj_in->nlayer-3;i>=0;i--)
			rbm_update_vis( &(obj_in->layers[i]) );

		ch = getch();
		if( ch == 'q' )
			break;
	}
	endwin();
}

inline void usage()
{
    printf( "%s\n", PACKAGE_STRING );
    printf( "Usage: rbmrun [params] [temp]\n" );
}

int main( int argc, char *argv[] )
{
    long i;
    double tmp;
    dbn_t dbn;

    if( argc < 3 )
    {
        usage();
        return 0;
    }

    /* Call the dbn trainer to do all the shit above in one command */
    dbn_load_init( &dbn, argv[1] );
    dbn_assemble( &dbn );

    for(i=0;i<dbn.nlayer-1;i++)
        rbm_set_temp( &dbn.layers[i], 0.1 );

    tmp = atof( argv[2] );
    for(i=0;i<dbn.nlayer-1;i++)
        rbm_set_temp( &dbn.layers[i], tmp );

    /* Initialize to random */
    for(i=0;i<dbn.layers[0].nv;i++)
	if( (double) rand() / (double) RAND_MAX < 0.5 )
	    dbn.layers[0].vis[i] = (char) 1;
	else
	    dbn.layers[0].vis[i] = (char) 0;

    /* Evolve forward */
    for(i=0;i<dbn.nlayer-1;i++)
        rbm_update_hid( &dbn.layers[i] );
    rbm_update_vis( &dbn.layers[dbn.nlayer-2] );

    for(i=0;i<512;i++)
    {
        rbm_update_hid( &dbn.layers[dbn.nlayer-2] );
        rbm_update_vis( &dbn.layers[dbn.nlayer-2] );
    }

    /* Run free via ncurses */
    dbn_run( &dbn );

    return 0;
}

