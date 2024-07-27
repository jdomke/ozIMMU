subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
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
      real(kind=c_double), intent(in), value :: oAlpha, oBeta
      real(kind=c_double), intent(in)    :: oA(oLda,*), oB(oLdb,*)
      real(kind=c_double), intent(inout) :: oC(oLdc,*)
      end subroutine offload_dgemm
      end interface
#else
      implicit none
#endif
      character :: uplo, trans
      integer(kind=4) :: i, j, n, k, lda, ldc, lay, ta, tb
      double precision :: alpha, beta
      double precision, dimension(lda,*) :: a
      double precision, dimension(ldc,*) :: c
      integer(kind=4) :: ka, istat
      character(len=255) :: errmsg
      double precision, dimension(:,:), allocatable :: pA, pB, pC

#ifndef DONTWRAPGEMM
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
!     .. Local Scalars ..
      LOGICAL upper
      !write(6,*) "DSYRK with ", uplo, trans, n, k, lda, ldc
      !call flush(6)
      lay = 0
      ta = 0
      tb = 1
      ! early return if possible
      if (n .le. 0) then
          return
      else if (k .le. 0) then
          if (upper) then
              do j = 1, n
                  c(1:j,j) = beta*c(1:j,j)
              end do
          else
              do j = 1, n
                  c(j:n,j) = beta*c(j:n,j)
              end do
          end if
      end if
      ! cuBLAS API Reference guide: For maximum compatibility with
      ! existing Fortran [...], the cuBLAS library uses column-major
      ! -> we are in F hence ColMajor, no need to revert back
      if (lsame(trans,'T')) then
          ta = 1
          tb = 0
      else if (lsame(trans,'C')) then
          ta = 2
          tb = 0
      end if
      upper = lsame(uplo,'U')
      ka = k
      if (ta .ne. 0) then
          ka = n
      end if
      ! re-use dgemm, so copy A into B and do some other tricks
      allocate(pA(1:max(1024,lda),1:max(1024,ka)),                    &
               pB(1:max(1024,lda),1:max(1024,ka)),                    &
               pC(1:max(1024,ldc),1:max(1024,n)),                     &
               STAT=istat, ERRMSG=errmsg)
      if (istat .ne. 0) then
          write(*,*) errmsg, " : istat =", istat
          call abort
      end if
      pA = 0
      pB = 0
      pC = 0
      pA(1:lda,1:ka) = a(1:lda,1:ka)
      pB(1:lda,1:ka) = a(1:lda,1:ka)
      !pC(1:ldc,1:n) = c(1:ldc,1:n)
      if (upper) then
          do j = 1, n
              pC(1:j,j) = c(1:j,j)
          end do
      else
          do j = 1, n
              pC(j:n,j) = c(j:n,j)
          end do
      end if
      call offload_dgemm(lay, ta, tb,                                 &
                         max(1024,n), max(1024,n), max(1024,k),       &
                         alpha,                                       &
                         pA, max(1024,lda),                           &
                         pB, max(1024,lda),                           &
                         beta,                                        &
                         pC, max(1024,ldc))
      if (upper) then
          do j = 1, n
              c(1:j,j) = pC(1:j,j)
          end do
      else
          do j = 1, n
              c(j:n,j) = pC(j:n,j)
          end do
      end if
      !write(6,*) "JJ", c(1:min(m,2),1:min(n,2)), "...", c(max(ldc-1,1):ldc,max(n-1,1):n)
      deallocate(pA, pB, pC, STAT=istat, ERRMSG=errmsg)
      if (istat .ne. 0) then
          write(*,*) errmsg, " : istat =", istat
          call abort
      end if
      return
#else
!https://netlib.org/lapack/explore-html/dc/d05/dsyrk_8f_source.html
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
      INTEGER INFO,L,NROWA
      LOGICAL UPPER
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
!     ..
!
!     Test the input parameters.
!
      !write(6,*) "DSYRK with ", uplo, trans, n, k, lda, ldc
      !call flush(6)
      IF (lsame(trans,'N')) THEN
          nrowa = n
      ELSE
          nrowa = k
      END IF
      upper = lsame(uplo,'U')
!
      info = 0
      IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
          info = 1
      ELSE IF ((.NOT.lsame(trans,'N')) .AND.                          &
               (.NOT.lsame(trans,'T')) .AND.                          &
               (.NOT.lsame(trans,'C'))) THEN
          info = 2
      ELSE IF (n.LT.0) THEN
          info = 3
      ELSE IF (k.LT.0) THEN
          info = 4
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 7
      ELSE IF (ldc.LT.max(1,n)) THEN
          info = 10
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('DSYRK ',info)
          RETURN
      END IF
!
!     Quick return if possible.
!
      IF ((n.EQ.0) .OR. (((alpha.EQ.zero).OR.                         &
          (k.EQ.0)).AND. (beta.EQ.one))) RETURN
!
!     And when  alpha.eq.zero.
!
      IF (alpha.EQ.zero) THEN
          IF (upper) THEN
              IF (beta.EQ.zero) THEN
                  DO 20 j = 1,n
                      DO 10 i = 1,j
                          c(i,j) = zero
   10                 CONTINUE
   20             CONTINUE
              ELSE
                  DO 40 j = 1,n
                      DO 30 i = 1,j
                          c(i,j) = beta*c(i,j)
   30                 CONTINUE
   40             CONTINUE
              END IF
          ELSE
              IF (beta.EQ.zero) THEN
                  DO 60 j = 1,n
                      DO 50 i = j,n
                          c(i,j) = zero
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 80 j = 1,n
                      DO 70 i = j,n
                          c(i,j) = beta*c(i,j)
   70                 CONTINUE
   80             CONTINUE
              END IF
          END IF
          RETURN
      END IF
!
!     Start the operations.
!
      IF (lsame(trans,'N')) THEN
!
!        Form  C := alpha*A*A**T + beta*C.
!
          IF (upper) THEN
              DO 130 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 90 i = 1,j
                          c(i,j) = zero
   90                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 100 i = 1,j
                          c(i,j) = beta*c(i,j)
  100                 CONTINUE
                  END IF
                  DO 120 l = 1,k
                      IF (a(j,l).NE.zero) THEN
                          temp = alpha*a(j,l)
                          DO 110 i = 1,j
                              c(i,j) = c(i,j) + temp*a(i,l)
  110                     CONTINUE
                      END IF
  120             CONTINUE
  130         CONTINUE
          ELSE
              DO 180 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 140 i = j,n
                          c(i,j) = zero
  140                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 150 i = j,n
                          c(i,j) = beta*c(i,j)
  150                 CONTINUE
                  END IF
                  DO 170 l = 1,k
                      IF (a(j,l).NE.zero) THEN
                          temp = alpha*a(j,l)
                          DO 160 i = j,n
                              c(i,j) = c(i,j) + temp*a(i,l)
  160                     CONTINUE
                      END IF
  170             CONTINUE
  180         CONTINUE
          END IF
      ELSE
!
!        Form  C := alpha*A**T*A + beta*C.
!
          IF (upper) THEN
              DO 210 j = 1,n
                  DO 200 i = 1,j
                      temp = zero
                      DO 190 l = 1,k
                          temp = temp + a(l,i)*a(l,j)
  190                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  200             CONTINUE
  210         CONTINUE
          ELSE
              DO 240 j = 1,n
                  DO 230 i = j,n
                      temp = zero
                      DO 220 l = 1,k
                          temp = temp + a(l,i)*a(l,j)
  220                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  230             CONTINUE
  240         CONTINUE
          END IF
      END IF
!
#endif
end subroutine
