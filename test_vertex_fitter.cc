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
  double *weight;
  long int *vertex_ids;
  long int *cluster_ids;
};

int main(int argc, char *argv[]) {
  //--------------------------------------------preproc
  int i, j;
  static const int NUM_VERTICES = 10;  // Typically on the order of a few hundred. Other
                                       // real-world applications are on the order of hundreds of thousands. This is all
                                       // for one event. Expect many events per second.
  static const int NUM_TRACKS_PER_VERTEX = 50;
  static const int NUM_TRACKS = NUM_VERTICES * NUM_TRACKS_PER_VERTEX;
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
  tracks.ids = std::make_unique<long int[]>(NUM_TRACKS);
  tracks.zs = std::make_unique<double[]>(NUM_TRACKS);
  tracks.weight = std::make_unique<double[]>(NUM_TRACKS);
  tracks.vertex_ids = std::make_unique<long int[]>(NUM_TRACKS);
  tracks.cluster_ids = std::make_unique<long int[]>(NUM_TRACKS);
  
  for (i = 0; i < NUM_TRACKS; i++) {
    // distribute as gaussian around the true points. This generates the data we
    // are allowed to observe in the processing stage.
    tracks.ids[i] = i;
    tracks.vertex_ids[i] = i / NUM_TRACKS_PER_VERTEX;  //assigned in order
    double track_pos = 0;
    for (j = 0; j < SAMPLE_NUM; j++) {
      srand(time(NULL));
      track_pos += (double)rand() / RAND_MAX;
    }
    tracks.zs[i] = (track_pos * 2 / (SAMPLE_NUM)) + TRUE_Z_VALS[i / NUM_TRACKS_PER_VERTEX];
    // TODO: assign clusters. i.e., don't just assign cluster==vertex. Need to
    //       simulate physical behavior that illustrates some of the challenges
    //       we'd see, or we need to load ground truths from actual measurements.
    tracks.cluster_ids[i] = tracks.vertex_ids[i];
  }

  // Filtering usually happens here in pre-proc. That gets rid of major
  // outliers (e.g., from instrumentation noise). We don't implement that in
  // this demo since we assume our generated data is "filtered".

  //--------------------------------------------proc
  //===========================================================================
  // 1. For each cluster (vertex) get sample mean and sample stddev
  // Naive 2-pass sample variance computation: Calculate mean then calculate ssqdiff
  //===========================================================================
  unsigned cluster_track_count[NUM_VERTICES]; // How many tracks are observed in a cluster
  double cluster_track_mean[NUM_VERTICES]; // Mean of track-x-pos observations in a cluster
  double cluster_track_std[NUM_VERTICES]; // Sample standard deviation of track-x-pos observations in a cluster
  for (i = 0; i < NUM_VERTICES; ++i) {
	  cluster_track_count[i] = 0;
	  cluster_track_mean[i] = 0;
	  cluster_track_std[i] = 0;
  }

  // Variance first pass: Calculate the mean.
  for (i = 0; i < NUM_TRACKS; ++i) {
    const long int cluster_id = tracks.cluster_ids[i];
    cluster_track_mean[cluster_id] += tracks.zs[i];
    cluster_track_count[cluster_id] += 1;
  }
  for (i = 0; i < NUM_VERTICES; ++i) {
    cluster_track_mean[i] /= cluster_track_count[i];
  }

  // Variance second pass: get the sum of square differences, divided by n-1
  for (i = 0; i < NUM_TRACKS; ++i) {
    const long int cluster_id = tracks.cluster_ids[i];
    const double diff = (tracks.zs[i] - cluster_track_mean[cluster_id]);
    cluster_track_std[cluster_id] += diff * diff;
  }
  for (i = 0; i < NUM_VERTICES; ++i) {
    // Calculate standard deviation from variance
    cluster_track_std[i] = sqrt(cluster_track_std[i] / (cluster_track_count[i] - 1));
  }

  //===========================================================================
  // 2. For each track, calculate distance from cluster's sample mean. If
  // greater than 3*std, then it is an outlier. Have a weight scalar: 0 if
  // outlier. Inliers get gaussian(cluster mean, cluster stddev)(track's Z pos)
  //===========================================================================
  // TODO: the stddev calculation contains distance from mean as an
  // intermediate step. May reuse here. Profile to see if needed (may
  // complicate other code transformations if we unnecessarily intertwine those
  // steps).
  double sum_of_weights = 0;
  for (i = 0; i < NUM_TRACKS; ++i) {
    const long int cluster_id = tracks.cluster_ids[i];
    const double diff = (tracks.zs[i] - cluster_track_mean[cluster_id]);
    if (diff > cluster_track_std[cluster_id] * 3) {
      tracks.weight[i] = 0;
    } else {
      // TODO: tracks.weight[i] = gaussian(cluster_track_mean[cluster_id], cluster_track_std[cluster_id])(tracks.zs[i])
      tracks.weight[i] = 1;
    }
    sum_of_weights += tracks.weight[i];
  }

  //===========================================================================
  // 3. For each track, get influence-weighted mean of z positions. This gives
  // z-position of cluster. This is our final output, i.e., the vertex
  // position. Assign to z_vals.
  //===========================================================================
  for (i = 0; i < NUM_VERTICES; ++i) {
    z_vals[i] = 0;
  }
  for (i = 0; i < NUM_TRACKS; ++i) {
    const long int cluster_id = tracks.cluster_ids[i];
    z_vals[cluster_id] += tracks.zs[i] * tracks.weight[i];
  }
  for (i = 0; i < NUM_VERTICES; ++i) {
    z_vals[i] /= sum_of_weights;
  }

  //--------------------------------------------postproc
  //get differences between calculated z values and real ones (efficiency). Use MSE from TRUE_Z_VALS to z_vals as our error calculation.
  //get time differences
  //analyze them and stuff
  double errs[NUM_VERTICES];
  double mean_square_error = 0;
  for (i = 0; i < NUM_VERTICES; i++) {
    double err = TRUE_Z_VALS[i] - z_vals[i];
    mean_square_error += err * err;
  }
  mean_square_error /= NUM_VERTICES;
  printf("Mean square error of z positions: %g\n", mean_square_error);

  // Above test data is perfectly clean, so we need to expect epsilon error.
  // TODO: Get testing data from Joshua's summer work: QCD 5, 10, 20, and 40. Use that to validate on non-clean data.

  return 0;
}
