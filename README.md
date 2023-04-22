# EC527 Project

### Joshua Shterenberg, Daniel Wilson

Simplified multicore and GPU versions of Track-Vertex Fitting on the CMSSW. This
is a mini-app version of the weighted fitter for the primary vertex producer
at https://github.com/joshshterenberg/cmssw/blob/from-CMSSW_12_6_0_pre5/RecoVertex/PrimaryVertexProducer/plugins/fitterCUDA.cc

### Dependencies
* [Eigen linear algebra library](https://eigen.tuxfamily.org/index.php?title=Main_Page) used throughout cmssw. E.g., `apt install libeigen3-dev`
* [lib fmt](https://fmt.dev/) used by cmssw logger. E.g., `apt install libfmt-dev`
* libtinyxml2 used by cmssw logger. E.g., `apt install libtinyxml2-dev`
* TBB used by cmssw throughout. E.g., `apt install libtbb-dev`

The cmssw_include directory contains select dependencies copied from the parent
cmssw project.

### Building
To build on SCC, first download and build all of the external dependencies by
running `install_deps.sh`. This only needs to be done once.

Next, set up your current shell's environment to be made aware of the downloaded
dependencies and of the dependencies that are already available in SCC modules
by running `source build_env.scc.sh`. This needs to be done once per shell
session.

Next, run `make -j4` to build this project.

### Running
To run on the SCC, request a node with GPUs. E.g., the following request
asks for a 2-hour session with GPUs: `qrsh -l gpus=1 -P ec527 -l h_rt=2:00:00`

Set up your shell session to be made aware of the project's dependencies and
to load a compatible CUDA runtime by running `source run_env.scc.sh`. This needs
to be done once per shell session in your `qrsh` runs.

Run the application. E.g., `./test_vertex_fitter`
