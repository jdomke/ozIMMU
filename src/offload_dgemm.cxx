#include <iostream>

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
    std::cerr << "offload_dgemm: l,ta,tb,m,n,k,lda,ldb,ldc=" << oLayout << oTransA << oTransB << oM << oN << oK << oLda << oLdb << oLdc << std::endl;;

    /*
    cublasHandle_t handle;
    cublasOperation_t transa = (0 == oTransA) ? CUBLAS_OP_N : CUBLAS_OP_T;
    cublasOperation_t transb = (0 == oTransB) ? CUBLAS_OP_N : CUBLAS_OP_T;
    int m = oM, n = oN, k = oK, lda = oLda, ldb = oLdb, ldc = oLdc;
    double alpha[1] = { oAlpha }, beta[1] = { oBeta };
    double *A = oA, *B = oB, *C = oC;
    cublasDgemm_v2(handle, transa, transb, m, n, k,
                   alpha, A, lda, B, ldb, beta, C, oLdc);
    */
}
