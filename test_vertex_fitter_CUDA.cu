/****************************************************************************
 * Mini-application to test performance optimizations for a weighted track
 * filter and vertex fitter.
 *
 * Authors: Joshua Shterenberg, Daniel Wilson
 *
 * This mini-application represents and is derived from code in the
 * https://github.com/cms-sw/cmssw project, which is released with an Apache-2.0 license.
 ***************************************************************************/
//
//nvcc -arch sm_35 test_vertex_fitter_CUDA.cu -o test_vertex_fitter_CUDA
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <memory>

// Assertion to check for errors
#define CUDA_SAFE_CALL(ans) { gpuAssert((ans), (char *)__FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
  if (code != cudaSuccess)
  {
    fprintf(stderr, "CUDA_SAFE_CALL: %s %s %d\n",
                                       cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                        __longlong_as_double(assumed)));
    } while (assumed != old);

    return __longlong_as_double(old);
}

#define PRINT_TIME		1
#ifndef NUM_VERTICES
#define NUM_VERTICES 1024
#endif
#ifndef NUM_TRACKS_PER_VERTEX
#define NUM_TRACKS_PER_VERTEX 2048
#endif
#define NUM_TRACKS		(NUM_VERTICES * NUM_TRACKS_PER_VERTEX)
#define SAMPLE_NUM		12
#ifndef VERTEX_COUNT
#define VERTEX_COUNT 1024
#endif
#ifndef TRACKS_PER_VERTEX
#define TRACKS_PER_VERTEX 2048
#endif

struct track_soa_t {
  long int* ids;
  double* zs;
  double* weight;
  long int* vertex_ids;
  long int* cluster_ids;
};

__global__ void proc(double* tracks_zs, double* tracks_weight, long int* tracks_cluster_ids, double* z_vals) {
  // current parallel strat: every track gets its own thread
  //   (should be ok because 1024 max threads/block)
  // assuming the grid is 1d
  // THIS IS A BAD IDEA IN PRACTICE BECAUSE NUM_TRACKS_PER_VERTEX ISN'T CONSTANT

  //assuming that first NUM_TRACKS_PER_VERTEX tracks are assoc with vertex 0, etc.

  const int i_track = blockIdx.x * blockDim.x + threadIdx.x; //track
  const int i_vertex = i_track / NUM_TRACKS_PER_VERTEX; //vertex
  __shared__ unsigned cluster_track_count[NUM_VERTICES];
  __shared__ double cluster_track_mean[NUM_VERTICES];
  __shared__ double cluster_track_std[NUM_VERTICES];
  __shared__ double cluster_sum_of_weights[NUM_VERTICES];

  cluster_track_count[i_vertex] = 0;
  cluster_track_mean[i_vertex] = 0;
  cluster_track_std[i_vertex] = 0;
  cluster_sum_of_weights[i_vertex] = 0;

  __syncthreads(); //always sync write/read clusters

  long int cluster_id = tracks_cluster_ids[i_track];
  atomicAdd(&cluster_track_mean[cluster_id], tracks_zs[i_track]);
  atomicAdd(&cluster_track_count[cluster_id], 1);
  
  __syncthreads();

  cluster_track_mean[i_vertex] /= cluster_track_count[i_vertex];

  __syncthreads();

  double diff = -1.0 * abs(tracks_zs[i_track] - cluster_track_mean[cluster_id]);
  atomicAdd(&cluster_track_std[cluster_id], diff*diff);

  __syncthreads();

  cluster_track_std[i_vertex] = sqrt(cluster_track_std[i_vertex] / (cluster_track_count[i_vertex] - 1));

  __syncthreads();

  if (diff <= cluster_track_std[cluster_id] * 3) {
    double xmstd = (tracks_zs[i_track] - cluster_track_mean[cluster_id]) / cluster_track_std[cluster_id];
    tracks_weight[i_track] = exp(-0.5 * xmstd * xmstd) / (cluster_track_std[cluster_id] * sqrt(2 * M_PI));
  } else tracks_weight[i_track] = 0;

  atomicAdd(&cluster_sum_of_weights[cluster_id], tracks_weight[i_track]);

  __syncthreads();
 
  atomicAdd(&z_vals[cluster_id], tracks_zs[i_track] * tracks_weight[i_track]);

  __syncthreads();

  z_vals[i_vertex] /= cluster_sum_of_weights[cluster_id];

}

int main(int argc, char *argv[]) {

  cudaEvent_t start, stop;
  float elapsed;

  //--------------------------------------------preproc
  int i, j;
  srand(time(NULL));

  double TRUE_Z_VALS[NUM_VERTICES], z_vals[NUM_VERTICES];
  for (i = 0; i < NUM_VERTICES; i++) {
    TRUE_Z_VALS[i] = -10 + (rand() * 20.0 / RAND_MAX);
    z_vals[i] = 0;
  }

  track_soa_t tracks;
  
  tracks.ids = (long int*) malloc(NUM_TRACKS * sizeof(long int));
  tracks.zs = (double*) malloc(NUM_TRACKS * sizeof(double));
  tracks.weight = (double*) malloc(NUM_TRACKS * sizeof(double));
  tracks.vertex_ids = (long int*) malloc(NUM_TRACKS * sizeof(long int));
  tracks.cluster_ids = (long int*) malloc(NUM_TRACKS * sizeof(long int));
  
  for (i = 0; i < NUM_TRACKS; i++) {
    tracks.ids[i] = i;
    tracks.vertex_ids[i] = i / NUM_TRACKS_PER_VERTEX;
    double track_pos = 0;
    for (j = 0; j < SAMPLE_NUM; j++) {
      track_pos += (double)rand() / RAND_MAX;
    }
    track_pos -= 6;
    tracks.zs[i] = (track_pos * 2 / (SAMPLE_NUM)) + TRUE_Z_VALS[i / NUM_TRACKS_PER_VERTEX];
    tracks.cluster_ids[i] = tracks.vertex_ids[i];
  }

  /////////////////////////CUDA////////////////////////
  CUDA_SAFE_CALL(cudaSetDevice(0));
  double *tracks_zs_gpu, *tracks_weight_gpu, *z_vals_gpu;
  long int *tracks_cluster_ids_gpu;
  size_t track_size = NUM_TRACKS * sizeof(double), 
         vertex_size = NUM_VERTICES * sizeof(double),
         int_size = NUM_TRACKS * sizeof(long int);

  CUDA_SAFE_CALL(cudaMalloc(&tracks_zs_gpu, track_size));
  CUDA_SAFE_CALL(cudaMalloc(&tracks_weight_gpu, track_size));
  CUDA_SAFE_CALL(cudaMalloc(&tracks_cluster_ids_gpu, int_size));
  CUDA_SAFE_CALL(cudaMalloc(&z_vals_gpu, vertex_size));
  CUDA_SAFE_CALL(cudaMemcpy(tracks_zs_gpu, tracks.zs, track_size, cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(tracks_weight_gpu, tracks.weight, track_size, cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(tracks_cluster_ids_gpu, tracks.cluster_ids, int_size, cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(z_vals_gpu, z_vals, vertex_size, cudaMemcpyHostToDevice));

  dim3 dimGrid(1);
  dim3 dimBlock(NUM_TRACKS);

#if PRINT_TIME
  // Create the cuda events
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  // Record event on the default stream
  cudaEventRecord(start, 0);
#endif
  proc<<<dimGrid, dimBlock>>>
    (tracks_zs_gpu, tracks_weight_gpu, tracks_cluster_ids_gpu, z_vals_gpu);
  CUDA_SAFE_CALL(cudaDeviceSynchronize());
  CUDA_SAFE_CALL(cudaPeekAtLastError());
#if PRINT_TIME
  // Stop and destroy the timer
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsed, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
#endif
  CUDA_SAFE_CALL(cudaMemcpy(z_vals, z_vals_gpu, vertex_size, cudaMemcpyDeviceToHost));
  /////////////////////////CUDA////////////////////////

  //--------------------------------------------postproc

  double mean_square_error = 0;
  for (i = 0; i < NUM_VERTICES; i++) {
    double err = TRUE_Z_VALS[i] - z_vals[i];
    mean_square_error += err * err;
  }
  mean_square_error /= NUM_VERTICES;
  printf("GPU implementation mean square error of z positions: %g\n", mean_square_error);
  printf("GPU implementation time (s): %f (msec)\n", elapsed);

  return 0;
}
