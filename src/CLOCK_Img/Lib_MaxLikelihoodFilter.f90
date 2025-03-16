
module Lib_MaxLikelihoodFilter
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_MaxLikelihoodFilter from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
   use Lib_ConjugateGradient
   use iso_fortran_env
   implicit none
   private

   public      ::      maxLikelihoodFilter

contains
!---^^^^^^^^

   subroutine maxLikelihoodFilter(f, g, s, t)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      compute a maximum likelihood filtering with range given by t
      !*      The probability of finding scene f given that the true solution is g
      !*      is  P = p1 p2
      !*      where   p1 = Exp[ - (f-g)^2/(2 s^2 ) ]
      !*              p2 = Exp[ - ( del^2(g) )^2/(2 t^2) ]
      !*      so the log likelihood is proportional to
      !*          lambda = -(t^2/2) (f-g)^2 - (s^2/2) ( del^2(g) )^2
      !*      which we can maximise wrt g with
      !*          d lambda / dg = 0
      !*      which is a set of linear equations in f,g,s,t

      real(kind=real64), dimension(:, :), intent(in)     ::      f
      real(kind=real64), dimension(:, :), intent(out)    ::      g
      real(kind=real64), intent(in)                    ::      s, t

      real(kind=real64), dimension(:, :), allocatable    ::      AA              !   (1:13,Nx Ny)    -   filter kernel
      integer, dimension(:, :), allocatable              ::      indx            !   (0:13,Nx Ny)    -   index of non-zero rows/cols

      integer         ::      Nx, Ny
      real(kind=real64), dimension(:), allocatable      ::      ff, gg
      real(kind=real64)       ::       ww, xx, eps
      integer         ::      ii, jj, ix, iy, jx, jy, kx, ky, dd, mi
      logical         ::      ok

      real(kind=real64), dimension(0:4), parameter      ::      kernel4 = (/20, -8, 2, 0, 1/)
      real(kind=real64), dimension(0:8), parameter      ::      kernel8 = (/36, 0, -14, 0, 9, -4, 0, 0, 4/)/9.0d0
      integer, parameter                               ::      bandwidth_kernel4 = 13
      integer, parameter                               ::      bandwidth_kernel8 = 25

      logical             ::      useKernel8 = .false.
      integer                                         ::      bandwidth
      real(kind=real64), dimension(:), allocatable      ::      kernel

      if (useKernel8) then
         allocate (kernel(0:size(kernel8) - 1))
         kernel = kernel8
         bandwidth = bandwidth_kernel8
      else
         allocate (kernel(0:size(kernel4) - 1))
         kernel = kernel4
         bandwidth = bandwidth_kernel4
      end if

      Nx = size(f, dim=1)
      Ny = size(f, dim=2)
      allocate (AA(bandwidth, Nx*Ny))
      allocate (indx(0:bandwidth, Nx*Ny))
      allocate (ff(Nx*Ny))
      allocate (gg(Nx*Ny))

      !---    pack f into vector
      do jj = 1, Ny
         ff((jj - 1)*Nx + 1:jj*Nx) = f(1:Nx, jj)
      end do
      gg = ff

      !---    compute laplacian squared kernel

      !kernel(0) = 20*(s*s)/(t*t)      !   note I'm setting central element at end        !   0   0   1   0   0
      !kernel(1) = -8*(s*s)/(t*t)                                                         !   0   2  -8   2   0
      !kernel(2) = 2*(s*s)/(t*t)                                                          !   1  -8  20  -8   1
      !kernel(3) = 0                                                                      !   0   2  -8   2   0
      !kernel(4) = (s*s)/(t*t)                                                            !   0   0   1   0   0
      !

      !kernel(0:4) = (/ -20,-8,2,0,1 /) * (s*s)/(t*t)

      !kernel(0:8) = ((/ 36,0,-14,0,9,-4,0,0,4 /)/9.0d0) * (s*s)/(t*t)

      kernel = kernel*(s*s)/(t*t)

      !---    compute the non-local filter kernel as a sparse matrix
      do iy = 1, Ny
         do ix = 1, Nx

            ii = ix + (iy - 1)*Nx
            mi = 0

            do ky = -2, 2
               jy = ky + iy
               ok = (jy*(Ny + 1 - jy) > 0)             !   inside region
               do kx = -2, 2
                  dd = kx*kx + ky*ky
                  if (dd >= size(kernel)) cycle               !   outside kernel

                  jx = kx + ix
                  if (ok .and. (jx*(Nx + 1 - jx) > 0)) then   !   inside region
                     jj = jx + (jy - 1)*Nx
                     if (ii == jj) cycle           !   note: am going to set central element at end...
                     mi = mi + 1
                     xx = kernel(dd)
                     ww = ww + xx
                     AA(mi, ii) = xx
                     indx(mi, ii) = jj
                  end if
               end do
            end do

            !---    normalise the kernel, add diagonal
            ww = sum(AA(1:mi, ii))
            mi = mi + 1
            indx(mi, ii) = ii
            AA(mi, ii) = 1.0d0 - ww
            indx(0, ii) = mi

         end do
      end do

      !---    now use conjugate gradients to find maximum
      do ii = 1, 3
         eps = 1.0d-8
         call conjgrad(AA, indx, gg, ff, eps)
         if (eps < 1.0d-8) exit
      end do

      g(1:Nx, 1:Ny) = reshape(gg, (/Nx, Ny/))

      return
   end subroutine maxLikelihoodFilter

end module Lib_MaxlikelihoodFilter

