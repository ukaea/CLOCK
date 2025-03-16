
module Lib_DataDoubler
!---^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_DataDoubler from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      module to reduce data set in size by 1/2 and then reexpand it
   use Lib_Splines1
   use iso_fortran_env
   implicit none
   private

   public      ::      halfLine, doubleLine
   public      ::      halfImage, doubleImage
   public      ::      halfField, doubleField

contains
!---^^^^^^^^

   subroutine halfLine(n, y2, y1, pbc)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      take the data from y2(1:n) and interpolate onto the line y1(1:(n+1)/2)
      integer, intent(in)                              ::      n
      real(kind=real64), dimension(:), intent(in)       ::      y2
      real(kind=real64), dimension(:), intent(out)      ::      y1
      logical, intent(in)                              ::      pbc

      integer             ::      ii, mm
      type(Spline1)       ::      ss
      real(kind=real64)   ::      xx

      ss = Spline1_ctor(y2, pbc)

      mm = (n + 1)/2
      y1 = 0.0d0
      if (pbc) then
         do ii = 1, mm
            xx = real((ii - 1)*n, kind=real64)/mm
            y1(ii) = Splint(ss, xx)
         end do
      else
         do ii = 1, mm
            xx = real((ii - 0.5d0)*(n - 1), kind=real64)/mm
            y1(ii) = Splint(ss, xx)
         end do
      end if
      return
   end subroutine halfLine

   subroutine doubleLine(n, y1, y2, pbc)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      take the data from y1(1:(n+1)/2) and interpolate onto the line y2(1:n)
      integer, intent(in)                              ::      n
      real(kind=real64), dimension(:), intent(in)       ::      y1
      real(kind=real64), dimension(:), intent(out)      ::      y2
      logical, intent(in)                              ::      pbc

      integer             ::      ii, mm
      type(Spline1)       ::      ss
      real(kind=real64)   ::      xx

      mm = (n + 1)/2
      ss = Spline1_ctor(y1(1:mm), pbc)
      if (pbc) then
         do ii = 1, n
            xx = real((ii - 1)*mm, kind=real64)/n
            y2(ii) = Splint(ss, xx)
         end do
      else
         do ii = 1, n
            xx = real((ii - 1)*mm, kind=real64)/(n - 1) - 0.5d0
            y2(ii) = Splint(ss, xx)
         end do
      end if
      return
   end subroutine doubleLine

   subroutine halfImage(nx, ny, y2, y1, pbc)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      take the data from y2(1:n) and interpolate onto the line y1(1:(n+1)/2)
      integer, intent(in)                              ::      nx, ny
      real(kind=real64), dimension(:, :), intent(in)     ::      y2
      real(kind=real64), dimension(:, :), intent(out)    ::      y1
      logical, intent(in)                              ::      pbc

      real(kind=real64), dimension(:, :), allocatable    ::      y2_tmp

      integer             ::      ii, mx
      mx = (nx + 1)/2
      allocate (y2_tmp(mx, ny))
      do ii = 1, ny
         call halfLine(nx, y2(:, ii), y2_tmp(:, ii), pbc)
      end do
      do ii = 1, mx
         call halfLine(ny, y2_tmp(ii, :), y1(ii, :), pbc)
      end do
      deallocate (y2_tmp)

      return
   end subroutine halfImage

   subroutine doubleImage(nx, ny, y1, y2, pbc)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      take the data from y1(1:(n+1)/2) and interpolate onto the line y2(1:n)
      integer, intent(in)                              ::      nx, ny
      real(kind=real64), dimension(:, :), intent(in)     ::      y1
      real(kind=real64), dimension(:, :), intent(out)    ::      y2
      logical, intent(in)                              ::      pbc

      real(kind=real64), dimension(:, :), allocatable    ::      y1_tmp

      integer             ::      ii, my
      my = (ny + 1)/2
      allocate (y1_tmp(nx, my))
      do ii = 1, my
         call doubleLine(nx, y1(:, ii), y1_tmp(:, ii), pbc)
      end do
      do ii = 1, nx
         call doubleLine(ny, y1_tmp(ii, :), y2(ii, :), pbc)
      end do
      deallocate (y1_tmp)

      return
   end subroutine doubleImage

   subroutine halfField(nx, ny, nz, y2, y1, pbc)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      take the data from y2(1:n) and interpolate onto the line y1(1:(n+1)/2)
      integer, intent(in)                                  ::      nx, ny, nz
      real(kind=real64), dimension(:, :, :), intent(in)       ::      y2
      real(kind=real64), dimension(:, :, :), intent(out)      ::      y1
      logical, intent(in)                                  ::      pbc

      real(kind=real64), dimension(:, :, :), allocatable      ::      y2_tmp

      integer             ::      ii, jj, mx, my

      mx = (nx + 1)/2
      my = (ny + 1)/2
      allocate (y2_tmp(mx, my, nz))
      do ii = 1, nz
         call halfImage(nx, ny, y2(:, :, ii), y2_tmp(:, :, ii), pbc)
      end do
      do jj = 1, my
         do ii = 1, mx
            call halfLine(nz, y2_tmp(ii, jj, :), y1(ii, jj, :), pbc)
         end do
      end do
      deallocate (y2_tmp)

      return
   end subroutine halfField

   subroutine doubleField(nx, ny, nz, y1, y2, pbc)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      take the data from y1(1:(n+1)/2) and interpolate onto the line y2(1:n)
      integer, intent(in)                              ::      nx, ny, nz
      real(kind=real64), dimension(:, :, :), intent(in)   ::      y1
      real(kind=real64), dimension(:, :, :), intent(out)  ::      y2
      logical, intent(in)                              ::      pbc

      real(kind=real64), dimension(:, :, :), allocatable    ::      y1_tmp

      integer             ::      ii, jj, mz

      mz = (nz + 1)/2
      allocate (y1_tmp(nx, ny, mz))
      do ii = 1, mz
         call doubleImage(nx, ny, y1(:, :, ii), y1_tmp(:, :, ii), pbc)
      end do
      do jj = 1, ny
         do ii = 1, nx
            call doubleLine(nz, y1_tmp(ii, jj, :), y2(ii, jj, :), pbc)
         end do
      end do
      deallocate (y1_tmp)

      return
   end subroutine doubleField

end module Lib_DataDoubler
