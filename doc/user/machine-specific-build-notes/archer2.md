# Building on ARCHER2

HPE Cray EX https://docs.archer2.ac.uk

```
module load cmake/3.21.3
module load PrgEnv-gnu
module swap gcc gcc/11.2.0
module load boost/1.72.0
module load parmetis/4.0.3
module load cray-hdf5-parallel

thisdir=$(readlink -f $(dirname $BASH_SOURCE))
build_dir=$thisdir/build
install_dir=$thisdir/install
source_dir=path/to/hemelb

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$install -DHEMELB_BUILD_RBC=ON -B $build_dir -S $source_dir
cmake --build build
```
