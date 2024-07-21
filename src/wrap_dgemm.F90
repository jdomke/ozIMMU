subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta,&
        c, ldc)
#ifndef DONTWRAPGEMM
      use iso_c_binding
      implicit none
      interface
      subroutine offload_dgemm(oLayout, oTransA, oTransB, oM, oN, oK, &
                               oAlpha, oA, oLda, oB, oLdb, oBeta,     &
                               oC, oLdc) bind(c)
      import c_int, c_double
      integer(kind=c_int), intent(in), value :: oLayout
      integer(kind=c_int), intent(in), value :: oTransA, oTransB
      integer(kind=c_int), intent(in), value :: oM, oN, oK
      integer(kind=c_int), intent(in), value :: oLda, oLdb, oLdc
      real(kind=c_double), intent(in)    :: oA(oLda,*), oB(oLdb,*)
      real(kind=c_double), intent(inout) :: oC(oLdc,*)
      end subroutine offload_dgemm
      end interface
#else
      implicit none
#endif
      character :: transa, transb
      integer(kind=4) :: m, n, k, lda, ldb, ldc, lay, ta, tb
      double precision :: alpha, beta
      double precision, dimension(lda,*) :: a
      double precision, dimension(ldb,*) :: b
      double precision, dimension(ldc,*) :: c

#ifndef DONTWRAPGEMM
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
      lay = 0
      ta = 0
      tb = 0
      ! we are in CblasRowMajor, need to revert back to cblas interface
      ! https://netlib.org/lapack/explore-html/dc/d18/cblas__dgemm_8c_source.html#l00102
      if (lsame(transa,'T')) then
          tb = 1
      else if (lsame(transa,'C')) then
          tb = 2
      end if
      if (lsame(transb,'T')) then
          ta = 1
      else if (lsame(transb,'C')) then
          ta = 2
      end if
      call offload_dgemm(lay, ta, tb, n, m, k,                        &
                         alpha, b, ldb, a, lda, beta, c, ldc)
#else
!https://netlib.org/lapack/explore-html/d7/d2b/dgemm_8f_source.html
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
!     ..
!     .. External Subroutines ..
      EXTERNAL xerbla
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC max
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,J,L,NROWA,NROWB
      LOGICAL NOTA,NOTB
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
!     ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA and NROWB  as the number of rows of  A
!     and  B  respectively.
!
      write(*,*) "DGEMM with ", transa, transb, m, n, k, lda, ldb, ldc
      nota = lsame(transa,'N')
      notb = lsame(transb,'N')
      IF (nota) THEN
          nrowa = m
      ELSE
          nrowa = k
      END IF
      IF (notb) THEN
          nrowb = k
      ELSE
          nrowb = n
      END IF
!
!     Test the input parameters.
!
      info = 0
      IF ((.NOT.nota) .AND. (.NOT.lsame(transa,'C')) .AND.            &
          (.NOT.lsame(transa,'T'))) THEN
          info = 1
      ELSE IF ((.NOT.notb) .AND. (.NOT.lsame(transb,'C')) .AND.       &
               (.NOT.lsame(transb,'T'))) THEN
          info = 2
      ELSE IF (m.LT.0) THEN
          info = 3
      ELSE IF (n.LT.0) THEN
          info = 4
      ELSE IF (k.LT.0) THEN
          info = 5
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 8
      ELSE IF (ldb.LT.max(1,nrowb)) THEN
          info = 10
      ELSE IF (ldc.LT.max(1,m)) THEN
          info = 13
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DGEMM ',info)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR.                                 &
          (((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one))) RETURN
!
!     And if  alpha.eq.zero.
!
      IF (alpha.EQ.zero) THEN
          IF (beta.EQ.zero) THEN
              DO 20 j = 1,n
                  DO 10 i = 1,m
                      c(i,j) = zero
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 j = 1,n
                  DO 30 i = 1,m
                      c(i,j) = beta*c(i,j)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
!
!     Start the operations.
!
      IF (notb) THEN
          IF (nota) THEN
!
!           Form  C := alpha*A*B + beta*C.
!
              DO 90 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 50 i = 1,m
                          c(i,j) = zero
   50                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 60 i = 1,m
                          c(i,j) = beta*c(i,j)
   60                 CONTINUE
                  END IF
                  DO 80 l = 1,k
                      temp = alpha*b(l,j)
                      DO 70 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          ELSE
!
!           Form  C := alpha*A**T*B + beta*C
!
              DO 120 j = 1,n
                  DO 110 i = 1,m
                      temp = zero
                      DO 100 l = 1,k
                          temp = temp + a(l,i)*b(l,j)
  100                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  110             CONTINUE
  120         CONTINUE
          END IF
      ELSE
          IF (nota) THEN
!
!           Form  C := alpha*A*B**T + beta*C
!
              DO 170 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 130 i = 1,m
                          c(i,j) = zero
  130                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 140 i = 1,m
                          c(i,j) = beta*c(i,j)
  140                 CONTINUE
                  END IF
                  DO 160 l = 1,k
                      temp = alpha*b(j,l)
                      DO 150 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
  150                 CONTINUE
  160             CONTINUE
  170         CONTINUE
          ELSE
!
!           Form  C := alpha*A**T*B**T + beta*C
!
              DO 200 j = 1,n
                  DO 190 i = 1,m
                      temp = zero
                      DO 180 l = 1,k
                          temp = temp + a(l,i)*b(j,l)
  180                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  190             CONTINUE
  200         CONTINUE
          END IF
      END IF
!
#endif
end subroutine
