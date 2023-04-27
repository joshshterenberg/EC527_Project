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
#include <memory>

struct track_soa_t {
  std::unique_ptr<long int[]> ids;
  std::unique_ptr<double[]> zs;
  std::unique_ptr<long int[]> vertex_ids;
  std::unique_ptr<long int[]> cluster_ids;
};

// Return the probability of drawing `x` on a Gaussian distribution defined by
// `mean` and `std`.
static inline double gaussian_pdf(double mean, double std, double x)
{
  const double xmstd = ((x - mean) / std);
  const double coeff = std * sqrt(2 * M_PI);
  //return exp(-0.5 * xmstd * xmstd) / coeff;
  return (1 - 0.5*pow(xmstd,2) + 0.25*pow(xmstd,4)) / coeff;
}

// Return the mean square error between two arrays of doubles.
static double mse(size_t length, double* v1, double* v2) {
  double sumsq = 0;
  for (size_t i = 0; i < length; i++) {
    double err = v1[i] - v2[i];
    sumsq += err * err;
  }
  return sumsq / length;
}

// Given a collection of track observations assigned to clusters, estimate the
// z positions of the tracks' vertices.
// Args:
// vertex_count: How many z_vals should be written
// z_vals: Estimated z positions of the vertices
// track_count: How many track observations are given
// tracks: Track observations
static void serial_vertex_fit(size_t vertex_count, double* z_vals, size_t track_count, const track_soa_t* tracks) {
  //===========================================================================
  // 1. For each cluster (vertex) get sample mean and sample stddev
  // Naive 2-pass sample variance computation: Calculate mean then calculate ssqdiff
  //===========================================================================
  // How many tracks are observed in a cluster
  std::unique_ptr<unsigned[]> cluster_track_count(new unsigned[vertex_count]);
  // Mean of track-x-pos observations in a cluster
  std::unique_ptr<double[]> cluster_track_mean(new double[vertex_count]);
  // Sample standard deviation of track-x-pos observations in a cluster
  std::unique_ptr<double[]> cluster_track_std(new double[vertex_count]);
  std::unique_ptr<double[]> track_weights(new double[track_count]);
  std::unique_ptr<double[]> cluster_sum_of_weights(new double[vertex_count]);
  size_t i;

  for (i = 0; i < vertex_count; ++i) {
    cluster_track_count[i] = 0;
    cluster_track_mean[i] = 0;
    cluster_track_std[i] = 0;
    cluster_sum_of_weights[i] = 0;
  }

  // Variance first pass: Calculate the mean.
  for (i = 0; i < track_count; ++i) {
    const long int cluster_id = tracks->cluster_ids[i];
    cluster_track_mean[cluster_id] += tracks->zs[i];
    cluster_track_count[cluster_id] += 1;
  }
  for (i = 0; i < vertex_count; ++i) {
    cluster_track_mean[i] /= cluster_track_count[i];
  }

  // Variance second pass: get the sum of square differences, divided by n-1
  for (i = 0; i < track_count; ++i) {
    const long int cluster_id = tracks->cluster_ids[i];
    const double diff = (tracks->zs[i] - cluster_track_mean[cluster_id]);
    cluster_track_std[cluster_id] += diff * diff;
  }
  for (i = 0; i < vertex_count; ++i) {
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
  for (i = 0; i < track_count; ++i) {
    const long int cluster_id = tracks->cluster_ids[i];
    const double diff = (tracks->zs[i] - cluster_track_mean[cluster_id]);
    if (diff > cluster_track_std[cluster_id] * 3) {
      track_weights[i] = 0;
    } else {
      // TODO: or instead, should I be doing cdf(some_max)-cdf(some_min)?
      track_weights[i] = gaussian_pdf(cluster_track_mean[cluster_id], cluster_track_std[cluster_id], tracks->zs[i]);
    }
    cluster_sum_of_weights[cluster_id] += track_weights[i];
  }

  //===========================================================================
  // 3. For each track, get influence-weighted mean of z positions. This gives
  // z-position of cluster. This is our final output, i.e., the vertex
  // position. Assign to z_vals.
  //===========================================================================
  for (i = 0; i < vertex_count; ++i) {
    z_vals[i] = 0;
  }
  for (i = 0; i < track_count; ++i) {
    const long int cluster_id = tracks->cluster_ids[i];
    z_vals[cluster_id] += tracks->zs[i] * track_weights[i];
  }
  for (i = 0; i < vertex_count; ++i) {
    z_vals[i] /= cluster_sum_of_weights[i];
  }
}

// Generate track positions (for test data) based on the true vertex positions
// this can be roughly Gaussian distributed for now, will need to match MC later
static void generate_tracks_from_vertices(size_t tracks_per_vertex, size_t track_count, track_soa_t* tracks, size_t vertex_count, const double* vertices) {
  static const unsigned SAMPLE_NUM = 12;  //related to Gaussian generation variance, CLT
  for (size_t i = 0; i < track_count; i++) {
    // distribute as gaussian around the true points. This generates the data we
    // are allowed to observe in the processing stage.
    tracks->ids[i] = i;
    tracks->vertex_ids[i] = i / tracks_per_vertex;  //assigned in order
    double track_pos = 0;
    // Transform an Irwin-Hall distribution to approximate a normal distribution
    for (unsigned j = 0; j < SAMPLE_NUM; j++) {
      track_pos += (double)rand() / RAND_MAX;
    }
    track_pos -= 6;  // Necessary to make the distribution central around the mean
    tracks->zs[i] = (track_pos * 2 / (SAMPLE_NUM)) + vertices[i / tracks_per_vertex];
    // TODO: assign clusters. i.e., don't just assign cluster==vertex. Need to
    //       simulate physical behavior that illustrates some of the challenges
    //       we'd see, or we need to load ground truths from actual measurements.
    tracks->cluster_ids[i] = tracks->vertex_ids[i];
  }
}

double interval(struct timespec start, struct timespec end) {
  struct timespec temp;
  temp.tv_sec = end.tv_sec - start.tv_sec;
  temp.tv_nsec = end.tv_nsec - start.tv_nsec;
  if (temp.tv_nsec < 0) {
    temp.tv_sec = temp.tv_sec - 1;
    temp.tv_nsec = temp.tv_nsec + 1000000000;
  }
  return (((double)temp.tv_sec) + ((double)temp.tv_nsec) * 1.0e-9);
}

int main(int argc, char *argv[]) {
  //--------------------------------------------preproc
  static const unsigned NUM_VERTICES = 20000;  // Typically on the order of a few hundred. Other
                                       // real-world applications are on the order of hundreds of thousands. This is all
                                       // for one event. Expect many events per second.
  static const unsigned NUM_TRACKS_PER_VERTEX = 50;
  static const unsigned NUM_TRACKS = NUM_VERTICES * NUM_TRACKS_PER_VERTEX;
  static const unsigned TRIAL_COUNT = 10; // How many times to recalculate z values for profiling.

  double TRUE_Z_VALS[NUM_VERTICES];  // Expected output
  double  z_vals[NUM_VERTICES]; // Actual output

  struct timespec time_start, time_end;
  double serial_elapsed = 0;

  srand(time(NULL));

  //create list of vertices based on known z values
  // this is the ground truth -- actual vertices
  for (unsigned i = 0; i < NUM_VERTICES; i++) {
    TRUE_Z_VALS[i] = -10 + (rand() * 20.0 / RAND_MAX); //detector is ~21 meters long
  }

  track_soa_t tracks;
  tracks.ids = std::make_unique<long int[]>(NUM_TRACKS);
  tracks.zs = std::make_unique<double[]>(NUM_TRACKS);
  tracks.vertex_ids = std::make_unique<long int[]>(NUM_TRACKS);
  tracks.cluster_ids = std::make_unique<long int[]>(NUM_TRACKS);
  generate_tracks_from_vertices(NUM_TRACKS_PER_VERTEX, NUM_TRACKS, &tracks, NUM_VERTICES, TRUE_Z_VALS);

  // Filtering track outliers usually also happens in pre-processing. That gets
  // rid of very strong outliers (e.g., from instrumentation noise). We don't
  // implement that in this demo since we assume our generated data is
  // "filtered".
  for (unsigned trial = 0; trial < TRIAL_COUNT; ++trial) {
    for (size_t i = 0; i < NUM_VERTICES; ++i) {
      z_vals[i] = 0;
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);
    serial_vertex_fit(NUM_VERTICES, z_vals, NUM_TRACKS, &tracks);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_end);
    serial_elapsed += interval(time_start, time_end);
  }
  serial_elapsed /= TRIAL_COUNT;

  //get time differences
  double mean_square_error = mse(NUM_VERTICES, TRUE_Z_VALS, z_vals);
  printf("Serial implementation mean square error of z positions: %g\n", mean_square_error);
  printf("Serial implementation time (s): %g\n", serial_elapsed);

  // Above test data is perfectly clean, so we need to expect epsilon error.
  // TODO: Get testing data from Joshua's summer work: QCD 5, 10, 20, and 40. Use that to validate on non-clean data.
  return 0;
}
