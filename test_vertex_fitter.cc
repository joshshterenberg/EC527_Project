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

struct track_soa_t {
  long int *ids;
  double *zs;
  long int *vertex_ids;
  long int *cluster_ids;
};

int main(int argc, char *argv[]) {
  //--------------------------------------------preproc
  int i, j;
  int NUM_VERTICES = 10; // Typically on the order of a few hundred. Other
                         // real-world applications are on the order of hundreds of thousands. This is all
                         // for one event. Expect many events per second.
  int NUM_TRACKS_PER_VERTEX = 50;
  int SAMPLE_NUM = 12;  //related to Gaussian generation variance, CLT

  //create list of vertices based on known z values
  double TRUE_Z_VALS[NUM_VERTICES], z_vals[NUM_VERTICES];
  // this is the ground truth -- actual vertices
  for (i = 0; i < NUM_VERTICES; i++) {
    srand(time(NULL));
    TRUE_Z_VALS[i] = -10 + (rand() * 20.0 / RAND_MAX); //detector is ~21 meters long
  }

  //use vertex list to generate tracks list w associations
  // this can be roughly Gaussian distributed for now, will need to match MC later
  track_soa_t tracks;
  tracks.ids = new long int[NUM_VERTICES * NUM_TRACKS_PER_VERTEX];
  tracks.zs = new double[NUM_VERTICES * NUM_TRACKS_PER_VERTEX];
  tracks.vertex_ids = new long int[NUM_VERTICES * NUM_TRACKS_PER_VERTEX];
  
  for (i = 0; i < NUM_VERTICES * NUM_TRACKS_PER_VERTEX; i++) {
    // distribute as gaussian around the true points. This generates the data we
    // are allowed to observe in the processing stage.
    tracks.ids[i] = i;
    tracks.vertex_ids[i] = i / NUM_VERTICES;  //assigned in order
    double track_pos = 0;
    for (j = 0; j < SAMPLE_NUM; j++) {
      srand(time(NULL));
      track_pos += (double)rand() / RAND_MAX;
    }
    tracks.zs[i] = (track_pos * 2 / (SAMPLE_NUM)) + TRUE_Z_VALS[i / NUM_VERTICES];
    // TODO: assign clusters. i.e., assign values to cluster_ids.
    //       -- How can we make this separate from vertices?
  }

  // Filtering usually happens here in pre-proc. That gets rid of major
  // outliers (e.g., from instrumentation noise). We don't implement that in
  // this demo since we assume our generated data is "filtered".

  //--------------------------------------------proc
  //1. For each cluster (vertex) get sample mean and sample stddev
  //2. For each track, calculate distance from cluster's sample mean. If GT 3stddev, then
  //   it is an outlier. Have an "influence" scalar in the loop, 0 if outlier.
  //   Inliers get weight == gaussian(cluster mean, cluster stddev)(track's Z position).
  //3. For each track, get influence-weighted mean of z positions. This gives
  //   z-position of cluster. This is our final output, i.e., the vertex position. Assign to z_vals.

  //--------------------------------------------postproc
  //get differences between calculated z values and real ones (efficiency). Use MSE from TRUE_Z_VALS to z_vals as our error calculation.
  //get time differences
  //analyze them and stuff
  double errs[NUM_VERTICES];
  for (i = 0; i < NUM_VERTICES; i++) {
    errs[i] = TRUE_Z_VALS[i] - z_vals[i];
  }

  // Above test data is perfectly clean, so we need to expect epsilon error.
  // TODO: Get testing data from Joshua's summer work: QCD 5, 10, 20, and 40. Use that to validate on non-clean data.

  return 0;
}
