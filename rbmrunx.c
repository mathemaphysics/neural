// vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab:smarttab
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include <X11/Xlib.h>
#include "include/rbm.h"
#include "config.h"

int main( int argc, char **argv )
{
	/* Display variables */
	Display *display;
	Window win;
	int i,j,k;
	int ch,screen_num,win_width,win_height,win_x,win_y;
	int win_border_width = 2;

	/* Check to make sure number of arguments is right */
	if( argc < 3 )
	{
		fprintf( stderr, "Usage: %s [dbn-params] [temp]\n", argv[0] );
		return 0;
	}

	/* Set up the belief net here */
	dbn_t dbn;
	dbn_load_init( &dbn, argv[1] );
	dbn_assemble( &dbn );
	dbn.temp = atof( argv[2] );

	/* Correctly set the temperature */
	for(i=0;i<dbn.nlayer-1;i++)
		dbn.layers[i].temp = dbn.temp;

	/* Initialize to a random vector */
	for(i=0;i<dbn.layers[0].nv;i++)
		if( rand() < ( RAND_MAX / 2 ) )
			dbn.layers[0].vis[i] = (char) 0;
		else
			dbn.layers[0].vis[i] = (char) 1;

	/* Actually open a connection to the X server */
	display = XOpenDisplay( getenv( "DISPLAY" ) );

	/* Find out which screen is the default for this display */
	screen_num = DefaultScreen( display );

	/* Set the width and height of the window to open */
	win_width = DisplayWidth( display, screen_num ) / 5;
	win_height = DisplayHeight( display, screen_num ) / 5;

	/* Set the position of the upper left corner of the window */
	win_x = 0;
	win_y = 0;

	if( display == NULL )
	{
		fprintf( stderr, "Cannot connect to X server\n" );
		exit( 1 );
	}

	/* Create the actual window data structure for mapping to the screen */
	win = XCreateSimpleWindow( display, RootWindow( display, screen_num ), win_x, win_y,
					win_width, win_height, win_border_width,
					BlackPixel( display, screen_num ),
					WhitePixel( display, screen_num ) );

	/* Write the stuff to the X server and sync */
	XMapWindow( display, win );
	XSync( display, False );

	/* Create the GC (graphics context) */
	XGCValues values;
	int line_style,cap_style,join_style;
	unsigned int line_width;
	unsigned long valuemask = 0;

	/* Values defining properties of geometric objects */
	line_width = 2;
	line_style = LineSolid;
	cap_style = CapButt;
	join_style = JoinBevel;

	/* Create the context */
	GC gc = XCreateGC( display, win, valuemask, &values );
	XSetForeground( display, gc, BlackPixel( display, screen_num ) );
	XSetBackground( display, gc, WhitePixel( display, screen_num ) );
	XSetLineAttributes( display, gc, line_width, line_style, cap_style, join_style );
	XSetFillStyle( display, gc, FillSolid );

	/* Draw some stuff into the window */
	XSync( display, False );

	/* Evolve forward */
	for(i=0;i<dbn.nlayer-1;i++)
		rbm_update_hid( &dbn.layers[i] );
	for(i=0;i<1000;i++)
	{
		rbm_update_vis( &dbn.layers[dbn.nlayer-2] );
		rbm_update_hid( &dbn.layers[dbn.nlayer-2] );
	}
	for(i=dbn.nlayer-2;i>=0;i--)
		rbm_update_vis( &dbn.layers[i] );

	/* Iterate through the DBN a few times and output the state to the X server */
	for(;;)
    {
        for(i=0;i<8;i++)
        {
            rbm_update_vis( &dbn.layers[dbn.nlayer-2] );
            rbm_update_hid( &dbn.layers[dbn.nlayer-2] );
        }
        for(i=dbn.nlayer-2;i>=0;i--)
            rbm_update_vis( &dbn.layers[i] );
        XClearWindow( display, win );
        for(j=0;j<dbn.nlayer-1;j++)
            for(i=0;i<dbn.layers[j].nv;i++)
                if( dbn.layers[j].vis[i] == 1 )
                    XDrawLine( display, win, gc, 5+5*i, 20+j*5, 5+5*i+3, 20+j*5 );
        for(i=0;i<dbn.layers[dbn.nlayer-2].nh;i++)
            if( dbn.layers[dbn.nlayer-2].hid[i] == 1 )
                XDrawLine( display, win, gc, 5+5*i, 20+(dbn.nlayer-1)*5, 5+5*i+3, 20+(dbn.nlayer-1)*5 );
        XFlush( display );
        usleep( 40000 );
    }

    /* Clean up the network trash */
    dbn_free( &dbn );

    /* Close the display and exit */
    XCloseDisplay( display );

    return 0;
}

