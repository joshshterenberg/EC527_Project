/****************************************************************************


   gcc -O1 -std=gnu11 test_SOR.c -lpthread -lrt -lm -o test_SOR

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

#define CPNS 3.0

typedef struct {
  long int* ids;
  double* zs;
  long int* vertex_ids;
} track_soa, *track_soa;

/* Prototypes */
void SOR(arr_ptr v, int *iterations);

double interval(struct timespec start, struct timespec end)
{
  struct timespec temp;
  temp.tv_sec = end.tv_sec - start.tv_sec;
  temp.tv_nsec = end.tv_nsec - start.tv_nsec;
  if (temp.tv_nsec < 0) {
    temp.tv_sec = temp.tv_sec - 1;
    temp.tv_nsec = temp.tv_nsec + 1000000000;
  }
  return (((double)temp.tv_sec) + ((double)temp.tv_nsec)*1.0e-9);
}

/*****************************************************************************/
int main(int argc, char *argv[]) {

  //--------------------------------------------preproc
  int OPTION;
  struct timespec time_start, time_stop;
  double time_stamp[OPTIONS][NUM_TESTS];
  int convergence[OPTIONS][NUM_TESTS];
  int *iterations;
  int i, j;
  int NUM_VERTICES = 10;
  int NUM_TRACKS_PER_VERTEX = 50;
  int SAMPLE_NUM = 12; //related to gaussean generation variance, CLT

  //create list of vertices based on known z values
  double TRUE_Z_VALS[NUM_VERTICES], z_vals[NUM_VERTICES];
  for (i = 0; i < NUM_VERTICES; i++) {
    srand(time(NULL));
    TRUE_Z_VALS[i] = -10 + (rand() * 20 / RAND_MAX);    
  }

  //use vertex list to generate tracks list w associations
  // this can be roughly gaussean distributed for now, will need to match MC later
  track_soa tracks[NUM_VERTICES * NUM_TRACKS_PER_VERTEX];
  for (i = 0; i < NUM_VERTICES * NUM_TRACKS_PER_VERTEX; i++) {
    tracks[i]->ids = i;
    tracks[i]->vertex_idx = i/NUM_VERTICES; //assigned in order
    double track_pos = 0;
    for (j = 0; j < SAMPLE_NUM; j++) {
      srand(time(NULL));
      track_pos += rand() / RAND_MAX; 
    }
    tracks[i]->zs = (track_pos * 2 / (SAMPLE_NUM)) + TRUE_Z_VALS[i/NUM_VERTICES];
  }

  //--------------------------------------------proc
  //get weighted average of all tracks for a given vertex, save
  //loop through tracks again
  //if track too far away (3sigma) influence = 0 else gaussian
  //weighted average of tracks including influence = vertex z
  //pass tracks list to GPU and do it all again

  printf("serial\n");
  for (x=0; x<NUM_TESTS && (n = A*x*x + B*x + C, n<=alloc_size); x++) {
    printf("  iter %d rowlen = %d\n", x, GHOST+n);
    init_array_rand(v0, GHOST+n);
    set_arr_rowlen(v0, GHOST+n);
    clock_gettime(CLOCK_REALTIME, &time_start);
    SOR(v0, iterations);
    clock_gettime(CLOCK_REALTIME, &time_stop);
    time_stamp[OPTION][x] = interval(time_start, time_stop);
    convergence[OPTION][x] = *iterations;
  }

  //--------------------------------------------postproc
  //get differences between calculated z values and real ones (efficiency)
  //get time differences
  //analyze them and stuff
  double errs[NUM_VERTICES];
  for (i = 0; i < NUM_VERTICES; i++) {
    errs[i] = TRUE_Z_VALS[i] - z_vals[i];
  }
  //printf(", %10.4g\n", (double)CPNS * 1.0e9 * time_stamp[OPTION][i]); //time

} /* end main */

/************************************/

/* SOR */
void SOR(arr_ptr v, int *iterations)
{
  long int i, j;
  long int rowlen = get_arr_rowlen(v);
  data_t *data = get_array_start(v);
  double change, total_change = 1.0e10;   /* start w/ something big */
  int iters = 0;

  while ((total_change/(double)(rowlen*rowlen)) > (double)TOL) {
    iters++;
    total_change = 0;
    for (i = 1; i < rowlen-1; i++) {
      for (j = 1; j < rowlen-1; j++) {
        change = data[i*rowlen+j] - .25 * (data[(i-1)*rowlen+j] +
                                          data[(i+1)*rowlen+j] +
                                          data[i*rowlen+j+1] +
                                          data[i*rowlen+j-1]);
        data[i*rowlen+j] -= change * OMEGA;
        if (change < 0){
          change = -change;
        }
        total_change += change;
      }
    }
    if (abs(data[(rowlen-2)*(rowlen-2)]) > 10.0*(MAXVAL - MINVAL)) {
      printf("SOR: SUSPECT DIVERGENCE iter = %ld\n", iters);
      break;
    }
  }
  *iterations = iters;
  printf("    SOR() done after %d iters\n", iters);
}


