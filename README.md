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
