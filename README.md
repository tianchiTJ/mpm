# High-Performance Material Point Method (CB-Geo mpm)
> [CB-Geo Computational Geomechanics Research Group](https://www.cb-geo.com)

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/cb-geo/mpm/develop/license.md)
[![Developer docs](https://img.shields.io/badge/developer-docs-blue.svg)](http://cb-geo.github.io/mpm)
[![User docs](https://img.shields.io/badge/user-docs-blue.svg)](https://mpm.cb-geo.com/)
[![CircleCI](https://circleci.com/gh/cb-geo/mpm.svg?style=svg)](https://circleci.com/gh/cb-geo/mpm)
[![codecov](https://codecov.io/gh/cb-geo/mpm/branch/develop/graph/badge.svg)](https://codecov.io/gh/cb-geo/mpm)
[![](https://img.shields.io/github/issues-raw/cb-geo/mpm.svg)](https://github.com/cb-geo/mpm/issues)
[![Project management](https://img.shields.io/badge/projects-view-ff69b4.svg)](https://github.com/cb-geo/mpm/projects/)

## Documentation

Please refer to [CB-Geo MPM Documentation](https://cb-geo.github.io/mpm-doc) for information on compiling, and running the code. The documentation also include the MPM theory.

## Install dependencies

* Docker image for CB-Geo mpm code [https://hub.docker.com/r/cbgeo/mpm](https://hub.docker.com/r/cbgeo/mpm)

* Instructions for running mpm docker container: [https://github.com/cb-geo/docker-mpm/blob/master/README.md](https://github.com/cb-geo/mpm-container/blob/master/README.md).

### Prerequisite packages
> The following prerequisite packages can be found in the docker image:

* [Boost](http://www.boost.org/)
* [Eigen](http://eigen.tuxfamily.org/)
* [Intel TBB](https://www.threadingbuildingblocks.org/)
* [HDF5](https://support.hdfgroup.org/HDF5/)

#### Optional
* [MPI](https://www.open-mpi.org/)
* [VTK](https://www.vtk.org/)

### Fedora installation

Please run the following command:

```shell
dnf install -y boost boost-devel clang cmake cppcheck eigen3-devel findutils gcc gcc-c++ \
                   git hdf5 hdf5-devel hdf5-openmpi hdf5-openmpi-devel kernel-devel lcov\
                   make openmpi openmpi-devel sqlite sqlite-devel tar tbb tbb-devel valgrind vim \
                   voro++ voro++-devel vtk vtk-devel wget
```

### Ubuntu installation

Please run the following commands to install dependencies:

```
sudo apt-get install -y cmake gcc git libboost-all-dev libeigen3-dev libhdf5-serial-dev libopenmpi-dev \
                        libtbb-dev libvtk7-dev

```

## Compile
> See https://mpm-doc.cb-geo.com/ for more detailed instructions.

0. Run `mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ ..`.

1. Run `make clean && make -jN` (where N is the number of cores).

### Compile mpm or mpmtest

* To compile either `mpm` or `mpmtest` alone, run `make mpm -jN` or `make mpmtest -jN` (where N is the number of cores).

### Compile without tests [Editing CMake options]

To compile without tests run: `mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DMPM_BUILD_TESTING=Off ..`.

### Run tests

0. Run `./mpmtest -s` (for a verbose output) or `ctest -VV`.

### Run MPM
> See https://mpm-doc.cb-geo.com/ for more detailed instructions.

The CB-Geo MPM code uses a `JSON` file for input configuration. To run the mpm code:

```
   ./mpm  [-p <tbb_parallel>] [-i <input_file>] -f <working_dir> [--]
          [--version] [-h]
```

For example:

```
./mpm -f /path/to/input-dir/ -i mpm-usf-3d.json -p 8
```

Where:

```

   -p <tbb_parallel>,  --tbb_parallel <tbb_parallel>
     Number of parallel TBB threads

   -i <input_file>,  --input_file <input_file>
     Input JSON file [mpm.json]

   -f <working_dir>,  --working_dir <working_dir>
     (required)  Current working folder

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```

## Compile with MPI (Running on a cluster)

The CB-Geo mpm code can be compiled with `MPI` to distribute the workload across compute nodes in a cluster.

Additional steps to load `OpenMPI` on Fedora:

```
source /etc/profile.d/modules.sh
export MODULEPATH=$MODULEPATH:/usr/share/modulefiles
module load mpi/openmpi-x86_64
```

Compile with OpenMPI:

```
mkdir build && cd build
export CXX_COMPILER=mpicxx
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_EXPORT_COMPILE_COMMANDS=On ..
make -jN
```

### Running the code with MPI

To run the CB-Geo mpm code on a cluster with MPI:

```
mpirun -N <#-MPI-tasks> ./mpm -f /path/to/input-dir/ -i mpm.json
```

For example to run the code on 4 compute nodes (MPI tasks):

```
mpirun -N 4 ./mpm -f ~/benchmarks/3d/uniaxial-stress -i mpm.json
```

## Input file format

### Input JSON

```
{
    "title" : "model_name",
    "input_files" : {
      "mesh" : "mesh.txt",
      "particles" : "particles.txt",
      "velocity_constraints" : "velocity-constraints.txt",
      "particles_volumes" : "particles-volumes.txt",
      "particles_stresses" : "particles-stresses.txt",
      "entity_sets": "entity-sets.json"
      },
    "mesh" : {
      "isoparametric": false,
      "cell_type": "ED2Q4",
      "mesh_reader": "Ascii2D",
      "node_type": "N2D"
      },
    "particle": {
      "material_id": 0,
      "particle_type": "P2D",
      "particle_sets": [
      {
        "set_id": 0,
        "initialise_material": true,
        "material_id": 1,
        "change_material": false,
        "remove": false
      },
      {
        "set_id": 1,
        "initialise_material": false,
        "change_material": true,
        "new_material_id": 1,
        "cmstep": 1000,
        "remove": false
      },
      {
        "set_id": 2,
        "initialise_material": false,
        "change_material": false,
        "remove": true,
        "rstep":1000
      }
      ]
      },
    "materials" : [
      {
        "id" : 0,
        "type" : "MohrCoulomb2D",
        "density" : 1800.0,
        "youngs_modulus" : 200000.0,
        "poisson_ratio" : 0.3
        "friction" : 0.0,
        "dilation" : 0.0,
        "cohesion" : 50000.0,
        "residual_friction" : 0.0,
        "residual_dilation" : 0.0,
        "residual_cohesion" : 50000.0,
        "peak_epds" : 0.0,
        "critical_epds" : 20000.0,
        "tension_cutoff" : 0.1,
        "softening" : false,
        "tolerance" : 0.1
      },
      {
        "id" : 1,
        "type" : "MohrCoulomb2D",
        "density" : 1800,
        "youngs_modulus" : 300000.0,
        "poisson_ratio" : 0.3,
        "friction" : 0.0,
        "dilation" : 0.0,
        "cohesion" : 50000.0,
        "residual_friction" : 0.0,
        "residual_dilation" : 0.0,
        "residual_cohesion" : 50000.0,
        "peak_epds" : 0.0,
        "critical_epds" : 20000.0,
        "tension_cutoff" : 0.1,
        "softening" : false,
        "tolerance" : 0.1
      }
      ],
    "analysis" : {
      "type" : "MPMExplicitUSF2D",
	    "velocity_update" : false,
      "strain_energy" : false,
	    "uuid": "model_name",
	    "gravity": [ 0.0, -9.81],
	    "dt" : 0.0001,
	    "nsteps" : 50000,
	    "resume" : {
        "resume": false,
        "uuid": "model_name",
        "step" : 40000
      }
      },
    "post_processing" : {
      "path" : "results/",
	    "vtk": [
              "displacements", "velocities",
              "strains", "shear_strains", "plastic_strains", "epds",
              "stresses", "shear_stresses",
              "strain_energy"
             ],
      "output_steps" : 1000
    }
}
```
