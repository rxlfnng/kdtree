/*! gcc -std=c89 -pedantic -Wall -g -o test2 test2.c libkdtree.a -lm */
/* Extended test program, contributed by David Underhill */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include "kdtree.h"

#define DEF_NUM_PTS 10


/* get a random double between -10 and 10 */
static double rd( void );

int main(int argc, char **argv) {
  int i, num_pts = DEF_NUM_PTS;
  void *ptree;
  char *data, *pch;
  kdres *presults;
  double pos[3], dist;
  double pt[3] = { 0, 0, 1 };
  double radius = 10;
  int num_near = 10;
  

  srand( time(0) );

  /* create a k-d tree for 3-dimensional points */
  ptree = kd_create( 3 );

  /* add some random nodes to the tree (assert nodes are successfully inserted) */
  for( i=0; i<num_pts; i++ ) {
    data[i] = 'a' + i;
    assert( 0 == kd_insert3( ptree, rd(), rd(), rd(), &data[i] ) );
  }

  /* find points num_near closest to the origin */
  presults = kd_nearest_n( ptree, pt, num_near );

  /* print out all the points found in results */
  printf( "found %d results:\n", kd_res_size(presults) );

  for(i=0; i<num_near; ++i {
    /* get the data and position of the current result item */
    pch = (char*)kd_res_item( presults, i );

    /* compute the distance of the current result from the pt */
    dist = kd_res_dis(presults, i);

    /* print out the retrieved data */
    printf( "node at (%.3f, %.3f, %.3f) is %.3f away and has data=%c\n", 
	    pos[0], pos[1], pos[2], dist, *pch );

  }

  /* free our tree, results set, and other allocated memory */
  free( data );
  kd_res_free( presults );
  kd_free( ptree );

  return 0;
}



static double rd( void ) {
  return (double)rand()/RAND_MAX * 20.0 - 10.0;
}
