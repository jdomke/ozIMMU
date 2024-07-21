#include <iostream>

// CUDA runtime
#include <cuda_runtime.h>
#include <cublas_v2.h>

static const char *_cudaGetErrorEnum(cudaError_t error)
{
    return cudaGetErrorName(error);
}

static const char *_cudaGetErrorEnum(cublasStatus_t error)
{
    switch (error) {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";
        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";
        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";
        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";
        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";
        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";
        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";
        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
        case CUBLAS_STATUS_NOT_SUPPORTED:
            return "CUBLAS_STATUS_NOT_SUPPORTED";
        case CUBLAS_STATUS_LICENSE_ERROR:
            return "CUBLAS_STATUS_LICENSE_ERROR";
    }
    return "<unknown>";
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

    double *A = NULL, *B = NULL, *C = NULL;
    int ka = ((!oTransA) ? oK : oM), kb = ((!oTransB) ? oN : oK);
    size_t memSizeA = oLda * ka * sizeof(*A), memSizeB = oLdb * kb * sizeof(*B);
    size_t memSizeC = oLdc * oN * sizeof(*C);

    checkCudaErrors(
            cudaMalloc((void**)&A, memSizeA)
            );
    checkCudaErrors(
            cudaMemcpy(A, oA, memSizeA, cudaMemcpyHostToDevice)
            );
    checkCudaErrors(
            cudaMalloc((void**)&B, memSizeB)
            );
    checkCudaErrors(
            cudaMemcpy(B, oB, memSizeB, cudaMemcpyHostToDevice)
            );
    checkCudaErrors(
            cudaMalloc((void**)&C, memSizeC)
            );
    checkCudaErrors(
            cudaMemcpy(C, oC, memSizeC, cudaMemcpyHostToDevice)
            );

    cublasHandle_t handle;
    checkCudaErrors(
            cublasCreate(&handle)
            );
    checkCudaErrors(
            cublasSetMathMode(handle, CUBLAS_DEFAULT_MATH)
            );

    cublasOperation_t transa = (0 == oTransA) ? CUBLAS_OP_N : CUBLAS_OP_T;
    cublasOperation_t transb = (0 == oTransB) ? CUBLAS_OP_N : CUBLAS_OP_T;
    int m = oM, n = oN, k = oK, lda = oLda, ldb = oLdb, ldc = oLdc;
    double alpha[1] = { oAlpha }, beta[1] = { oBeta };
    checkCudaErrors(
            cublasDgemm_v2(handle, transa, transb, m, n, k,
                           alpha, A, lda, B, ldb, beta, C, oLdc)
            );

    checkCudaErrors(
            cudaMemcpy(oC, C, memSizeC, cudaMemcpyDeviceToHost)
            );
}
