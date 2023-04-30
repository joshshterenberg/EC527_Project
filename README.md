# EC527 Project

### Joshua Shterenberg, Daniel Wilson

Simplified multicore and GPU versions of Track-Vertex Fitting on the CMSSW. This
is a mini-app version of the weighted fitter for the primary vertex producer
at https://github.com/joshshterenberg/cmssw/blob/from-CMSSW_12_6_0_pre5/RecoVertex/PrimaryVertexProducer/plugins/fitterCUDA.cc

### Building
To build on SCC set up your current shell's environment to be made aware of the
downloaded dependencies and of the dependencies that are already available in
SCC modules by running `source build_env.scc.sh`. This needs to be done once
per shell session.

To build on eng-grid, run `source build_env.eng-grid.sh`. This only needs to be
done once per shell session, before running make.

Next, run `make` to build this project.

To change the problem size, re-build as `make -B CXXFLAGS='-DNUM_VERTICES=1024 -DNUM_TRACKS_PER_VERTEX=1024'`

### Running
To run on the SCC, request a node with GPUs. E.g., the following request
asks for a 2-hour session with GPUs: `qrsh -l gpus=1 -P ec527 -l h_rt=2:00:00`

Set up your shell session to be made aware of the project's dependencies and
to load a compatible CUDA runtime by running `source run_env.scc.sh`. This needs
to be done once per shell session in your `qrsh` runs.

Run the application. E.g., `./test_vertex_fitter`
