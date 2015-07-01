#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include <ncurses.h>
#include "parse.h"
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

//void dbn_print_curses( dbn_t *obj_in )
//{
//        int i,j;
//        for(i=0;i<obj_in->nlayer-1;i++)
//                for(j=0;j<obj_in->layers[i].nv;j++)
//                        mvprintw( i, j, "%d", (int) obj_in->layers[i].vis[j] );
//        for(j=0;j<obj_in->layers[obj_in->nlayer-2].nh;j++)
//                mvprintw( obj_in->nlayer - 1, j, "%d", (int) obj_in->layers[obj_in->nlayer-2].hid[j] );
//        refresh();
//}

//void dbn_run( dbn_t *obj_in )
//{
//        /* Assume that the dbn has been properly set up already */
//        long i,j;
//        char ch;
//
//        /* Initialize ncurses */
//        initscr();
//        raw();
//        keypad( stdscr, TRUE );
//        noecho();
//
//        for(;;)
//        {
//                /* Print the output first, then update */
//                dbn_print_curses( obj_in );
//
//                /* Do the updating here */
//                for(i=0;i<4;i++)
//                {
//                        rbm_update_hid( &(obj_in->layers[obj_in->nlayer-2]) );
//                        rbm_update_vis( &(obj_in->layers[obj_in->nlayer-2]) );
//                }
//                for(i=obj_in->nlayer-3;i>=0;i--)
//                        rbm_update_vis( &(obj_in->layers[i]) );
//
//                ch = getch();
//                if( ch == 'q' )
//                       break;
//        }
//        endwin();
//}

inline void usage()
{
    printf( "%s\n", PACKAGE_STRING );
    printf( "Usage: rbmtrain [nlayer] [lsizes] [train-file]\n" );
}

int main( int argc, char *argv[] )
{
    char buf[512],*tok[128],*train,*samp;
    long h,i,j,k,n,nl,nt,max,*lsizes;
    double *cor,*vav,*hav;
    FILE *fp;
    dbn_t dbn;

    if( argc < 4 )
    {
        usage();
        return 0;
    }

    /* Set up here by reading the number of layers and layers sizes in argv[1] and argv[2] */
    nl = atoi( argv[1] );
    strcpy( buf, argv[2] );
    n = parse_stokenize( buf, tok, "," );
    if( n != nl )
    {
        printf( "Number of layers and layer sizes cardinality do not match. Exiting.\n" );
        return 0;
    }
    lsizes = (long*) malloc( nl * sizeof(long) );
    for(i=0;i<n;i++)
    {
        lsizes[i] = atoi( tok[i] );
        if( i == 0 )
            max = lsizes[i];
        else
            if( lsizes[i] > max )
                max = lsizes[i];
    }
    cor = (double*) malloc( max * max * sizeof(double) );
    vav = (double*) malloc( max * sizeof(double) );
    hav = (double*) malloc( max * sizeof(double) );
    samp = (char*) malloc( max * max * sizeof(char) );

    /* Start here by loading the training data and determining the number of layers and their sizes */
    if( argv[3] != NULL )
    {
        fp = fopen( argv[3], "r" );
        if( fp == NULL )
        {
            printf( "Could not open training file. Exiting.\n" );
            return 0;
        }

        /* Read values here */
        n = parse_read_line( fp, buf );
        if( n <= 0 )
        {
            printf( "Could not read header. Exiting.\n" );
            return 0;
        }
        buf[n] = '\0'; /* Make sure null-terminated */
        nt = atoi( buf );
        if( nt <= 0 )
        {
            printf( "Zero or negative number of sample vectors. Exiting.\n" );
            return 0;
        }
        train = (char*) malloc( nt * lsizes[0] * sizeof(char) );
        for(i=0;i<nt;i++)
        {
            /* Read each line as a training vector of characters which must be converted to 0's and 1's */
            n = parse_read_line( fp, buf );
            if( n <= 0 )
            {
                printf( "Number of training vectors specified is not consistent within the file. Exiting.\n" );
                fclose( fp );
                return 0;
            }

            /* Read line character-by-character */
            for(j=0;j<n;j++)
                if( buf[j] == '0' )
                    train[i*lsizes[0]+j] = 0;
                 else
                    train[i*lsizes[0]+j] = 1;
        }
    }

    printf( "Number of layers: %d\n", nl );
    printf( "Layer sizes: " );
    for(i=0;i<nl;i++)
        printf( "%10d", lsizes[i] );
    printf( "\n" );

    /* Call the dbn trainer to do all the shit above in one command */
    dbn_init( &dbn, nl, lsizes, 5.0 );

    printf( "Training set: %d\n", nt );
    for(i=0;i<nt;i++)
    {
        for(j=0;j<lsizes[0];j++)
            printf( "%d", train[i*lsizes[0]+j] );
        printf( "\n" );
    }

    /* Train */
    dbn_cd_mc( &dbn, nt, train, 1000, 100, 1, 5000, cor, vav, hav, samp, 0.15 );

    /* Output */
    dbn_assemble( &dbn );
    dbn_save_init( &dbn, "params.dbn" );

    return 0;
}

