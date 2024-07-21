#include <iostream>

// CUDA runtime
#include <cuda_runtime.h>
#include <cublas_v2.h>

static const char *_cudaGetErrorEnum(cudaError_t error)
{
    return cudaGetErrorName(error);
}

template <typename T>
void check(T result, char const *const func, const char *const file,
           int const line)
{
    if (result) {
        fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n", file,
                line, static_cast<unsigned int>(result),
                _cudaGetErrorEnum(result), func);
        exit(EXIT_FAILURE);
    }
}

// This will output the proper CUDA error strings in the event
// that a CUDA host call returns an error
#define checkCudaErrors(val) check((val), #val, __FILE__, __LINE__)


extern "C" void offload_dgemm(const int oLayout,
		              const int oTransA, const int oTransB,
                              const int oM, const int oN, const int oK,
                              const double oAlpha,
                              const double * oA, const int oLda,
                              const double * oB, const int oLdb,
                              const double oBeta,
                              double * oC, const int oLdc);

void offload_dgemm(const int oLayout,
		   const int oTransA, const int oTransB,
                   const int oM, const int oN, const int oK,
                   const double oAlpha,
                   const double * oA, const int oLda,
                   const double * oB, const int oLdb,
                   const double oBeta,
                   double * oC, const int oLdc)
{
    std::cerr << "offload_dgemm: l,ta,tb,m,n,k,lda,ldb,ldc="
        << oLayout << ", " << oTransA << ", " << oTransB << ", "
        << oM << ", " << oN << ", " << oK << ", "
        << oLda << ", " << oLdb << ", " << oLdc << std::endl;

    cudaError_t cudaStat;
    cublasStatus_t stat;
    cublasHandle_t handle;

    double *A = NULL, *B = NULL, *C = NULL;
    int ka = ((!oTransA) ? oK : oM), kb = ((!oTransB) ? oN : oK);
    int memSizeA = oLda * ka * sizeof(*A), memSizeB = oLdb * kb * sizeof(*B);
    int memSizeC = oLdc * oN * sizeof(*C);

    checkCudaErrors(
            cudaMalloc((void**)&A, memSizeA)
            );
    checkCudaErrors(
            cudaMemcpy(A, oA, memSizeA, cudaMemcpyHostToDevice)
            );
    checkCudaErrors(
            cudaMalloc((void**)&A, memSizeB)
            );
    checkCudaErrors(
            cudaMemcpy(B, oB, memSizeB, cudaMemcpyHostToDevice)
            );
    checkCudaErrors(
            cudaMalloc((void**)&A, memSizeC)
            );
    checkCudaErrors(
            cudaMemcpy(C, oC, memSizeC, cudaMemcpyHostToDevice)
            );

    cublasOperation_t transa = (0 == oTransA) ? CUBLAS_OP_N : CUBLAS_OP_T;
    cublasOperation_t transb = (0 == oTransB) ? CUBLAS_OP_N : CUBLAS_OP_T;
    int m = oM, n = oN, k = oK, lda = oLda, ldb = oLdb, ldc = oLdc;
    double alpha[1] = { oAlpha }, beta[1] = { oBeta };
//    checkCudaErrors(
            cublasDgemm_v2(handle, transa, transb, m, n, k,
                           alpha, A, lda, B, ldb, beta, C, oLdc);
//            );

    checkCudaErrors(
            cudaMemcpy(oC, C, memSizeC, cudaMemcpyDeviceToHost)
            );
}
