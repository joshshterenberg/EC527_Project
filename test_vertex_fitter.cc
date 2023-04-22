/****************************************************************************
 * Mini-application to test performance optimizations for a weighted track
 * filter and vertex fitter.
 *
 * Authors: Joshua Shterenberg, Daniel Wilson
 *
 * This mini-application represents and is derived from code in the
 * https://github.com/cms-sw/cmssw project, which is released with an Apache-2.0 license.
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include "CUDADataFormats/Track/interface/TrackForPVHeterogeneous.h"

struct track_soa_t {
  long int *ids;
  double *zs;
  long int *vertex_ids;
};

int main(int argc, char *argv[]) {
  //--------------------------------------------preproc
  int i, j;
  int NUM_VERTICES = 10;
  int NUM_TRACKS_PER_VERTEX = 50;
  int SAMPLE_NUM = 12;  //related to Gaussian generation variance, CLT

  TrackForPVHeterogeneous CPUtracks(cms::cuda::make_host_unique<TrackForPV::TrackForPVSoA>(cudaStreamDefault));
  TrackForPVHeterogeneous GPUtracks(cms::cuda::make_device_unique<TrackForPV::TrackForPVSoA>(cudaStreamDefault));

  //create list of vertices based on known z values
  double TRUE_Z_VALS[NUM_VERTICES], z_vals[NUM_VERTICES];
  for (i = 0; i < NUM_VERTICES; i++) {
    srand(time(NULL));
    TRUE_Z_VALS[i] = -10 + (rand() * 20.0 / RAND_MAX);
  }

  //use vertex list to generate tracks list w associations
  // this can be roughly Gaussian distributed for now, will need to match MC later
  track_soa_t tracks;
  tracks.ids = new long int[NUM_VERTICES * NUM_TRACKS_PER_VERTEX];
  tracks.zs = new double[NUM_VERTICES * NUM_TRACKS_PER_VERTEX];
  tracks.vertex_ids = new long int[NUM_VERTICES * NUM_TRACKS_PER_VERTEX];
  
  for (i = 0; i < NUM_VERTICES * NUM_TRACKS_PER_VERTEX; i++) {
    tracks.ids[i] = i;
    tracks.vertex_ids[i] = i / NUM_VERTICES;  //assigned in order
    double track_pos = 0;
    for (j = 0; j < SAMPLE_NUM; j++) {
      srand(time(NULL));
      track_pos += (double)rand() / RAND_MAX;
    }
    tracks.zs[i] = (track_pos * 2 / (SAMPLE_NUM)) + TRUE_Z_VALS[i / NUM_VERTICES];
  }

  //--------------------------------------------proc
  //get weighted average of all tracks for a given vertex, save
  //loop through tracks again
  //if track too far away (3sigma) influence = 0 else gaussian
  //weighted average of tracks including influence = vertex z
  //pass tracks list to GPU and do it all again

  //--------------------------------------------postproc
  //get differences between calculated z values and real ones (efficiency)
  //get time differences
  //analyze them and stuff
  double errs[NUM_VERTICES];
  for (i = 0; i < NUM_VERTICES; i++) {
    errs[i] = TRUE_Z_VALS[i] - z_vals[i];
  }

  return 0;
}
