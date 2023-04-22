module load cuda intel gcc/12.2.0 eigen
export CUDART_CFLAGS="-I${SCC_CUDA_INCLUDE}"
export CUDART_LDFLAGS="-L${SCC_CUDA_LIB}64" # The SCC-provided variable gives the wrong dir for 64-bit
export EIGEN_CFLAGS="-I${SCC_EIGEN_INCLUDE}/eigen3"
export FMT_CFLAGS="-I${PWD}/thirdparty/usr/local/include"
export FMT_LDFLAGS="-L${PWD}/thirdparty/usr/local/lib64"
export TINYXML_CFLAGS="-I${PWD}/thirdparty/usr/local/include"
export TINYXML_LDFLAGS="-L${PWD}/thirdparty/usr/local/lib64"
export LIBMD_CFLAGS="-I${PWD}/thirdparty/usr/include"
export LIBMD_LDFLAGS="-L${PWD}/thirdparty/usr/lib64"
