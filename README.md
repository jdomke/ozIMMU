# ozIMMU - DGEMM on Int8 Tensor Core

This library intercepts function calls for cuBLAS DGEMM functions and executes ozIMMU instead

## Prepare
```bash
wget https://github.com/Kitware/CMake/releases/download/v3.29.7/cmake-3.29.7-linux-aarch64.sh
mkdir -p cmake; bash ./cmake-3.29.7-linux-aarch64.sh --prefix=$(pwd)/cmake --skip-license
export PATH="$(pwd)/cmake/bin:${PATH}"
```

## Build
```bash
git clone https://github.com/jdomke/ozIMMU --recursive
cd ozIMMU
git checkout low-prec-gracehopper
sed -i -e 's@CUDA_ARCHITECTURES.*@CUDA_ARCHITECTURES 90)@g' CMakeLists.txt
cmake -B build \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    -DCUDAToolkit_ROOT=/usr/local/cuda \
    -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc \
    -DCMAKE_CUDA_ARCHITECTURES=90
cmake --build build --config Release
OZIMMU_LIB=$(pwd)/build/libozimmu.so
```

## Test example
```bash
git clone https://github.com/NVIDIA/cuda-samples.git
sed -i -e 's@= . \* block_size@= 32 * block_size@g' -e 's@Sgemm@Dgemm@g' -e 's@float@double@g' -e 's@double msecTotal@float msecTotal@g' cuda-samples/Samples/4_CUDA_Libraries/matrixMulCUBLAS/matrixMulCUBLAS.cpp
sed -i -e 's@sdkCompareL2fe(const.*@sdkCompareL2fe(const double *reference, const double *data,@g' cuda-samples/Common/helper_image.h
make -C cuda-samples/Samples/4_CUDA_Libraries/matrixMulCUBLAS/ SMS="90"
# check if it runs
cuda-samples/Samples/4_CUDA_Libraries/matrixMulCUBLAS/matrixMulCUBLAS -device=0 -sizemult=1
# try with ozaki scheme
OZIMMU_COMPUTE_MODE=fp64_int8_4 OZIMMU_INFO=1 OZIMMU_ERROR=1 OZIMMU_ENABLE_CULIP_PROFILING=1 LD_PRELOAD=${OZIMMU_LIB} \
    cuda-samples/Samples/4_CUDA_Libraries/matrixMulCUBLAS/matrixMulCUBLAS -device=0 -sizemult=1
```

## Test RI-MP2 quantum chemistry code
```bash
# compile with -Knolargepage,nohpctag,simd_reg_size=agnostic
# move dependencies to GH
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${PWD}/libs"
# prep mpi env
export OMPI_HOME="/opt/nvidia/hpc_sdk/Linux_aarch64/24.3/comm_libs/openmpi/openmpi-3.1.5"
export PATH="${OMPI_HOME}/bin:${PATH}"
export LD_LIBRARY_PATH="${OMPI_HOME}/lib:/opt/nvidia/hpc_sdk/Linux_aarch64/24.3/compilers/lib:${LD_LIBRARY_PATH}"
ulimit -s hard
MPIEXTRA=("-mca" "orte_base_help_aggregate" "0")
MPIEXTRA+=("-mca" "btl_base_warn_component_unused" "0")
MPIEXTRA+=("--oversubscribe")
MPIEXTRA+=("-x" "OMP_STACKSIZE=512M" "-x" "OMP_NUM_THREADS=1")
[ ! -f INPUT ] && tar xzf data.tar.gz
mpirun ${MPIEXTRA[@]} -np 4 ./rimp2_nohpc_fast.exe
```

## Test RI-MP2 with ozIMMU
```bash
OZIEXTRA=("-x" "LD_PRELOAD=${OZIMMU_LIB}")
OZIEXTRA+=("-x" "OZIMMU_COMPUTE_MODE=fp64_int8_4")
OZIEXTRA+=("-x" "OZIMMU_INFO=1")
OZIEXTRA+=("-x" "OZIMMU_ERROR=1")
OZIEXTRA==("-x" "OZIMMU_ENABLE_CULIP_PROFILING=1")
mpirun ${MPIEXTRA[@]} ${OZIEXTRA[@]} -np 4 ./rimp2_nohpc_fast.exe
```

## Evaluate slice-count w/ RI-MP2 and ozIMMU
```bash
for SC in $(seq 3 18); do
    OZIEXTRA=("-x" "LD_PRELOAD=${OZIMMU_LIB}")
    OZIEXTRA+=("-x" "OZIMMU_COMPUTE_MODE=fp64_int8_${SC}")
    mpirun ${MPIEXTRA[@]} ${OZIEXTRA[@]} -np 4 ./rimp2_nohpc_fast.exe
done
```

The supported compute modes are [here](#supported-compute-mode).

3. Execute the application

### Supported compute mode
| Mode          | Tensor Core type | Num splits |                         |
|:--------------|:-----------------|:-----------|:------------------------|
|dgemm          | --               | --         | Disable hijacking       |
|fp64_int8_3    | Int8 TC          | 3          |                         |
|fp64_int8_4    | Int8 TC          | 4          |                         |
|fp64_int8_5    | Int8 TC          | 5          |                         |
|fp64_int8_6    | Int8 TC          | 6          |                         |
|fp64_int8_7    | Int8 TC          | 7          |                         |
|fp64_int8_8    | Int8 TC          | 8          |                         |
|fp64_int8_9    | Int8 TC          | 9          |                         |
|fp64_int8_10   | Int8 TC          | 10         |                         |
|fp64_int8_11   | Int8 TC          | 11         |                         |
|fp64_int8_12   | Int8 TC          | 12         |                         |
|fp64_int8_13   | Int8 TC          | 13         |                         |
|fp64_int8_14   | Int8 TC          | 14         |                         |
|fp64_int8_15   | Int8 TC          | 15         |                         |
|fp64_int8_16   | Int8 TC          | 16         |                         |
|fp64_int8_17   | Int8 TC          | 17         |                         |
|fp64_int8_18   | Int8 TC          | 18         |                         |
|fp64_int8_auto | Int8 TC          | AUTO       | fp64_int8_3..18 / dgemm |


### Optional environmental variables
```bash
# Show info log
export OZIMMU_INFO=1

# Show error and warning log
export OZIMMU_ERROR=1

# Show CULiP ( https://github.com/enp1s0/CULiP ) log
export OZIMMU_ENABLE_CULIP_PROFILING=1

# Choose malloc mode
export OZIMMU_MALLOC_ASYNC=1

# Set AUTO mode mantissa loss threshold
export OZIMMU_AUTO_AVG_MANTISSA_LOSS_THRESHOLD=1.5
```

## Citation
```bibtex
@article{ootomo2024dgemm,
    author = {Hiroyuki Ootomo and Katsuhisa Ozaki and Rio Yokota},
    title = {DGEMM on integer matrix multiplication unit},
    journal = {The International Journal of High Performance Computing Applications},
    year = {2024},
    doi = {10.1177/10943420241239588},
    URL = {https://doi.org/10.1177/10943420241239588}
}
```

## License
MIT
