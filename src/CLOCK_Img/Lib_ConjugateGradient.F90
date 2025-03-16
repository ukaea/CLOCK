
module Lib_ConjugateGradient
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_ConjugateGradient from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
!*    Copyright (C) 2024  Daniel Mason

!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.

!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.

!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!*-----------------------------------------------------------------------------------------------------------------------------------
!*      perform a conjugate gradient relaxation to solve the symmetric equations
!*          A x = b
!*      where A is sparse and symmetric and real

   use OMP_LIB
   use iso_fortran_env
   implicit none
   private

   public      ::      conjgrad

contains
!---^^^^^^^^

   subroutine conjgrad(A, indx, x, b, eps)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      On input A_ik is the kth non-zero column of row i of the sparse symmetric matrix A
      !*      indx_0i is the number of non-zero columns of row i
      !*      indx_ki is the column index of the kth non-zero column
      !*      x is the initial guess, b is the target
      !*      eps is the tolerance for the solution - want |Ax-b|^2 < N eps where N is the number of rows
      !*      on output
      !*      x is the solution, eps is the error |Ax-b|^2/N

      real(kind=real64), dimension(:, :), intent(in)         ::      A
      integer, dimension(0:, :), intent(in)                  ::      indx
      real(kind=real64), dimension(:), intent(in)           ::      b
      real(kind=real64), dimension(:), intent(inout)        ::      x
      real(kind=real64), intent(inout)                     ::      eps

      real(kind=real64), dimension(:), allocatable      ::      rr, pp, qq
      real(kind=real64)                               ::      r2, aa, bb

      integer                     ::      ii, kk, NN, maxSteps

      NN = size(x)
      allocate (rr(NN))
      allocate (pp(NN))
      allocate (qq(NN))

      call sparseMatVec(A, indx, x, rr)
      rr(1:NN) = b(1:NN) - rr(1:NN)
      pp(1:NN) = rr(1:NN)
      maxSteps = ceiling(sqrt(1.0d0*NN))

      do kk = 1, maxSteps

         !---    r2 = r.r
         r2 = 0.0d0; do ii = 1, NN; r2 = r2 + rr(ii)*rr(ii); end do
         if (r2 <= NN*eps) exit

         !---    q = A p
         call sparseMatVec(A, indx, pp, qq)

         !---    a = p.q
         aa = 0.0d0; do ii = 1, NN; aa = aa + pp(ii)*qq(ii); end do

         !---    a = r2/a
         if (abs(aa) < 1.0d-12) exit
         !   weird result: p.Ap = 0 can only be true if p is a zero eigenmode of A, or if p=0
         !   if this is true, the best solution is probably to try again with a different guess vector
         aa = r2/aa

         !---    x = x + a p
         x(1:NN) = x(1:NN) + aa*pp(1:NN)

         !---    r = r - a q
         rr(1:NN) = rr(1:NN) - aa*qq(1:NN)

         !---    aa = r.r/r2
         bb = 0.0d0; do ii = 1, NN; bb = bb + rr(ii)*rr(ii); end do
         bb = bb/r2          !   note r2>0

         !---   p = r + ap
         pp(1:NN) = rr(1:NN) + bb*pp(1:NN)

      end do

      !---    at this point we should have the solution
      call sparseMatVec(A, indx, x, rr)
      rr(1:NN) = b(1:NN) - rr(1:NN)
      r2 = 0.0d0; do ii = 1, NN; r2 = r2 + rr(ii)*rr(ii); end do
      eps = r2/NN

      deallocate (rr)
      deallocate (pp)
      deallocate (qq)

      return
   end subroutine conjgrad

   subroutine sparseMatVec(A, indx, x, y)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      compute sparse matrix multiply Ax = y
      !*      On input A_ik is the kth non-zero column of row i of the sparse symmetric matrix A
      !*      indx_0i is the number of non-zero columns of row i
      !*      indx_ki is the column index of the kth non-zero column

      real(kind=real64), dimension(:, :), intent(in)         ::      A
      integer, dimension(0:, :), intent(in)                  ::      indx
      real(kind=real64), dimension(:), intent(in)           ::      x
      real(kind=real64), dimension(:), intent(out)          ::      y
      real(kind=real64)           ::      yi
      integer                     ::      ii, jj, kk, NN
      NN = size(x)

!$OMP PARALLEL PRIVATE(ii,jj,kk,yi) SHARED(x,y,A,indx,NN)
!$OMP DO
      do ii = 1, NN
         yi = 0                              !   compute y_i = A_ij x_j
         do kk = 1, indx(0, ii)                !   for each non-zero column at row i
            jj = indx(kk, ii)                !   kth non-zero column is indexed j
            yi = yi + A(kk, ii)*x(jj)        !   sum A_ij x_j
         end do
         y(ii) = yi
      end do
!$OMP END DO
!$OMP END PARALLEL

      return
   end subroutine sparseMatVec

end module Lib_ConjugateGradient

