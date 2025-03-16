
module Lib_LowPassFilter3d
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_LowPassFilter3d from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
   use iso_fortran_env
   implicit none
   private

   !---

   public      ::      generateFourierCoefficients
   public      ::      interpolateFromFourierCoefficients

   !---

#ifdef DEBUG
   logical, private, parameter         ::      LowPassFilter3d_dbg = .true.
#else
   logical, private, parameter         ::      LowPassFilter3d_dbg = .false.
#endif

   real(kind=real64), private, parameter                     ::      PI = 3.14159265358979d0

   !---

   interface interpolateFromFourierCoefficients
      module procedure interpolateFromFourierCoefficients20
      module procedure interpolateFromFourierCoefficients21
      module procedure interpolateFromFourierCoefficients22
      module procedure interpolateFromFourierCoefficients220
      module procedure interpolateFromFourierCoefficients31
      module procedure interpolateFromFourierCoefficients32
   end interface

   interface generateFourierCoefficients
      module procedure generateFourierCoefficients20
      module procedure generateFourierCoefficients21
      module procedure generateFourierCoefficients3
   end interface

contains
!---^^^^^^^^
   subroutine generateFourierCoefficients20(x, dat, xmax, nFourierCoeffs, qdat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given the data dat(1:n) at 2d points x(1:2,n) ranged from 0:xmax(1:2)
      !*      compute the first few fourier coefficients
      real(kind=real64), dimension(:, :), intent(in)     ::      x
      real(kind=real64), dimension(:), intent(in)       ::      dat
      real(kind=real64), dimension(2), intent(in)       ::      xmax
      integer, dimension(2), intent(in)                 ::      nFourierCoeffs
      real(kind=real64), dimension(:, 0:, 0:), intent(out)      ::      qdat          !   (1:4,0:nFourierCoeffs,0:nFourierCoeffs)

      integer                     ::          ix, iy
      real(kind=real64)           ::          qx, qy

      integer                     ::          ii, nn
      real(kind=real64)           ::          cx, sx, cy, sy, xx, yy
      real(kind=real64)           ::          d2, d1

      nn = size(dat)

      qdat = 0.0d0
      d1 = 2*PI/xmax(1)
      d2 = 2*PI/xmax(2)

      do ii = 1, nn
         xx = x(1, ii)
         yy = x(2, ii)

         do iy = 0, nFourierCoeffs(2)
            qy = yy*d2*iy
            cy = cos(qy)
            sy = sin(qy)

            do ix = 0, nFourierCoeffs(1)
               qx = xx*d1*ix
               cx = cos(qx)
               sx = sin(qx)

               qdat(1, ix, iy) = qdat(1, ix, iy) + cx*cy*dat(ii)
               qdat(2, ix, iy) = qdat(2, ix, iy) + sx*cy*dat(ii)
               qdat(3, ix, iy) = qdat(3, ix, iy) + cx*sy*dat(ii)
               qdat(4, ix, iy) = qdat(4, ix, iy) + sx*sy*dat(ii)

            end do
         end do
      end do

      xx = 1.0d0/nn

      do iy = 0, nFourierCoeffs(2)
         d2 = 2; if (iy == 0) d2 = 1
         do ix = 0, nFourierCoeffs(1)
            d1 = 2; if (ix == 0) d1 = 1
            qdat(:, ix, iy) = qdat(:, ix, iy)*d1*d2*xx
         end do
      end do

      return
   end subroutine generateFourierCoefficients20

   subroutine generateFourierCoefficients21(x, dat, xmax, nFourierCoeffs, qdat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given the data dat(1:d,1:n) at 2d points x(1:2:n) ranged from 0:xmax(1:2)
      !*      compute the first few fourier coefficients
      real(kind=real64), dimension(:, :), intent(in)     ::      x
      real(kind=real64), dimension(:, :), intent(in)     ::      dat
      real(kind=real64), dimension(2), intent(in)       ::      xmax
      integer, dimension(2), intent(in)                 ::      nFourierCoeffs
      real(kind=real64), dimension(:, :, 0:, 0:), intent(out)      ::      qdat          !   (1:d,1:4,0:nFourierCoeffs,0:nFourierCoeffs)

      integer                     ::          ix, iy
      real(kind=real64)           ::          qx, qy

      integer                     ::          ii, nn, dd
      real(kind=real64)           ::          cx, sx, cy, sy, xx, yy
      real(kind=real64)           ::          d2, d1

      dd = size(dat, dim=1)
      nn = size(dat, dim=2)

      qdat = 0.0d0
      d1 = 2*PI/xmax(1)
      d2 = 2*PI/xmax(2)

      do ii = 1, nn
         xx = x(1, ii)
         yy = x(2, ii)

         do iy = 0, nFourierCoeffs(2)
            qy = yy*d2*iy
            cy = cos(qy)
            sy = sin(qy)

            do ix = 0, nFourierCoeffs(1)
               qx = xx*d1*ix
               cx = cos(qx)
               sx = sin(qx)

               qdat(1:dd, 1, ix, iy) = qdat(1:dd, 1, ix, iy) + cx*cy*dat(1:dd, ii)
               qdat(1:dd, 2, ix, iy) = qdat(1:dd, 2, ix, iy) + sx*cy*dat(1:dd, ii)
               qdat(1:dd, 3, ix, iy) = qdat(1:dd, 3, ix, iy) + cx*sy*dat(1:dd, ii)
               qdat(1:dd, 4, ix, iy) = qdat(1:dd, 4, ix, iy) + sx*sy*dat(1:dd, ii)

            end do
         end do
      end do

      xx = 1.0d0/nn

      do iy = 0, nFourierCoeffs(2)
         d2 = 2; if (iy == 0) d2 = 1
         do ix = 0, nFourierCoeffs(1)
            d1 = 2; if (ix == 0) d1 = 1
            qdat(:, :, ix, iy) = qdat(:, :, ix, iy)*d1*d2*xx
         end do
      end do

      return
   end subroutine generateFourierCoefficients21

   subroutine generateFourierCoefficients3(x, dat, xmax, nFourierCoeffs, qdat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given the data dat(1:d,1:n) at 3d points x(1:3,1:n) ranged from 0:xmax(1:3)
      !*      compute the first few fourier coefficients
      real(kind=real64), dimension(:, :), intent(in)     ::      x
      real(kind=real64), dimension(:, :), intent(in)     ::      dat
      real(kind=real64), dimension(3), intent(in)       ::      xmax
      integer, dimension(3), intent(in)                 ::      nFourierCoeffs
      real(kind=real64), dimension(:, :, 0:, 0:, 0:), intent(out)      ::      qdat          !   (1:d,1:8,0:nFourierCoeffs,0:nFourierCoeffs,0:nFourierCoeffs)

      integer                     ::          ix, iy, iz
      real(kind=real64)           ::          qx, qy, qz

      integer                     ::          ii, nn, dd
      real(kind=real64)           ::          cx, sx, cy, sy, cz, sz, xx, yy, zz
      real(kind=real64)           ::          d3, d2, d1
      real(kind=real64), dimension(1:size(dat, dim=1))  ::      dcz, dsz

      dd = size(dat, dim=1)
      nn = size(dat, dim=2)

      qdat = 0.0d0
      d1 = 2*PI/xmax(1)
      d2 = 2*PI/xmax(2)
      d3 = 2*PI/xmax(3)

      do ii = 1, nn
         xx = x(1, ii)
         yy = x(2, ii)
         zz = x(3, ii)

         do iz = 0, nFourierCoeffs(3)
            qz = zz*d3*iz
            cz = cos(qz)
            sz = sin(qz)
            dcz(1:dd) = dat(1:dd, ii)*cz
            dsz(1:dd) = dat(1:dd, ii)*sz

            do iy = 0, nFourierCoeffs(2)
               qy = yy*d2*iy
               cy = cos(qy)
               sy = sin(qy)

               do ix = 0, nFourierCoeffs(1)
                  qx = xx*d1*ix
                  cx = cos(qx)
                  sx = sin(qx)

                  qdat(1:dd, 1, ix, iy, iz) = qdat(1:dd, 1, ix, iy, iz) + cx*cy*dcz(1:dd)
                  qdat(1:dd, 2, ix, iy, iz) = qdat(1:dd, 2, ix, iy, iz) + sx*cy*dcz(1:dd)
                  qdat(1:dd, 3, ix, iy, iz) = qdat(1:dd, 3, ix, iy, iz) + cx*sy*dcz(1:dd)
                  qdat(1:dd, 4, ix, iy, iz) = qdat(1:dd, 4, ix, iy, iz) + sx*sy*dcz(1:dd)
                  qdat(1:dd, 5, ix, iy, iz) = qdat(1:dd, 5, ix, iy, iz) + cx*cy*dsz(1:dd)
                  qdat(1:dd, 6, ix, iy, iz) = qdat(1:dd, 6, ix, iy, iz) + sx*cy*dsz(1:dd)
                  qdat(1:dd, 7, ix, iy, iz) = qdat(1:dd, 7, ix, iy, iz) + cx*sy*dsz(1:dd)
                  qdat(1:dd, 8, ix, iy, iz) = qdat(1:dd, 8, ix, iy, iz) + sx*sy*dsz(1:dd)

               end do
            end do
         end do
      end do

      xx = 1.0d0/nn

      do iz = 0, nFourierCoeffs(3)
         d3 = 2; if (iz == 0) d3 = 1
         do iy = 0, nFourierCoeffs(2)
            d2 = 2; if (iy == 0) d2 = 1
            do ix = 0, nFourierCoeffs(1)
               d1 = 2; if (ix == 0) d1 = 1
               qdat(:, :, ix, iy, iz) = qdat(:, :, ix, iy, iz)*d1*d2*d3*xx
            end do
         end do
      end do

      return
   end subroutine generateFourierCoefficients3

   !---

   function interpolateFromFourierCoefficients20(x, xmax, nFourierCoeffs, qdat) result(dat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(2), intent(in)       ::      x
      real(kind=real64), dimension(2), intent(in)       ::      xmax
      integer, dimension(2), intent(in)                 ::      nFourierCoeffs
      real(kind=real64), dimension(:, 0:, 0:), intent(in)      ::      qdat                 !   (1:d,1:4,0:nFourierCoeffs,0:nFourierCoeffs)
      real(kind=real64)                               ::      dat

      integer                     ::          ix, iy
      real(kind=real64)           ::          qx, qy

      real(kind=real64)           ::          cx, sx, cy, sy
      real(kind=real64)           ::          d2, d1

      d1 = x(1)*2*PI/xmax(1)
      d2 = x(2)*2*PI/xmax(2)

      dat = 0

      do iy = 0, nFourierCoeffs(2)
         qy = d2*iy
         sy = sin(qy)
         cy = cos(qy)

         do ix = 0, nFourierCoeffs(1)
            qx = d1*ix
            sx = sin(qx)
            cx = cos(qx)

            dat = dat + ((qdat(1, ix, iy)*cx &
                          + qdat(2, ix, iy)*sx)*cy &
                         + (qdat(3, ix, iy)*cx &
                            + qdat(4, ix, iy)*sx)*sy)
         end do
      end do

      return
   end function interpolateFromFourierCoefficients20

   function interpolateFromFourierCoefficients21(x, xmax, nFourierCoeffs, qdat) result(dat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(2), intent(in)       ::      x
      real(kind=real64), dimension(2), intent(in)       ::      xmax
      integer, dimension(2), intent(in)                 ::      nFourierCoeffs
      real(kind=real64), dimension(:, :, 0:, 0:), intent(in)      ::      qdat                 !   (1:d,1:4,0:nFourierCoeffs,0:nFourierCoeffs)
      real(kind=real64), dimension(size(qdat, dim=1))   ::      dat

      integer                     ::          ix, iy
      real(kind=real64)           ::          qx, qy

      integer                     ::          dd
      real(kind=real64)           ::          cx, sx, cy, sy
      real(kind=real64)           ::          d2, d1

      dd = size(dat)

      d1 = x(1)*2*PI/xmax(1)
      d2 = x(2)*2*PI/xmax(2)

      dat = 0

      do iy = 0, nFourierCoeffs(2)
         qy = d2*iy
         sy = sin(qy)
         cy = cos(qy)

         do ix = 0, nFourierCoeffs(1)
            qx = d1*ix
            sx = sin(qx)
            cx = cos(qx)

            dat(1:dd) = dat(1:dd) + ((qdat(1:dd, 1, ix, iy)*cx &
                                      + qdat(1:dd, 2, ix, iy)*sx)*cy &
                                     + (qdat(1:dd, 3, ix, iy)*cx &
                                        + qdat(1:dd, 4, ix, iy)*sx)*sy)
         end do
      end do

      return
   end function interpolateFromFourierCoefficients21

   function interpolateFromFourierCoefficients22(nx, ny, nFourierCoeffs, qdat) result(dat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      interpolate at each point of a regularly spaced lattice
      integer, intent(in)                              ::      nx, ny
      integer, dimension(2), intent(in)                 ::      nFourierCoeffs
      real(kind=real64), dimension(:, :, 0:, 0:), intent(in)      ::      qdat                 !   (1:d,1:4,0:nFourierCoeffs,0:nFourierCoeffs)
      real(kind=real64), dimension(size(qdat, dim=1), 0:nx - 1, 0:ny - 1)   ::      dat

      integer                     ::          ix, iy
      integer                     ::          jx, jy
      real(kind=real64)           ::          qx, qy

      integer                     ::          dd
      real(kind=real64)           ::          cx, sx, cy, sy
      real(kind=real64)           ::          d2, d1

      real(kind=real64), dimension(size(qdat, dim=1))       ::      ddat

      dd = size(qdat, dim=1)

      d1 = 2*PI/nx
      d2 = 2*PI/ny

      dat = 0

      do jy = 0, ny - 1
      do iy = 0, nFourierCoeffs(2)
         qy = d2*iy*jy
         sy = sin(qy)
         cy = cos(qy)

         do jx = 0, nx - 1
         do ix = 0, nFourierCoeffs(1)
            qx = d1*ix*jx
            sx = sin(qx)
            cx = cos(qx)

            ddat(1:dd) = ((qdat(1:dd, 1, ix, iy)*cx &
                           + qdat(1:dd, 2, ix, iy)*sx)*cy &
                          + (qdat(1:dd, 3, ix, iy)*cx &
                             + qdat(1:dd, 4, ix, iy)*sx)*sy)

            dat(1:dd, jx, jy) = dat(1:dd, jx, jy) + ddat(1:dd)

         end do
         end do
      end do
      end do

      return
   end function interpolateFromFourierCoefficients22

   function interpolateFromFourierCoefficients220(nx, ny, nFourierCoeffs, qdat) result(dat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      interpolate at each point of a regularly spaced lattice
      integer, intent(in)                              ::      nx, ny
      integer, dimension(2), intent(in)                 ::      nFourierCoeffs
      real(kind=real64), dimension(:, 0:, 0:), intent(in)      ::      qdat                 !   (1:4,0:nFourierCoeffs,0:nFourierCoeffs)
      real(kind=real64), dimension(0:nx - 1, 0:ny - 1)   ::      dat

      integer                     ::          ix, iy
      integer                     ::          jx, jy
      real(kind=real64)           ::          qx, qy

      real(kind=real64)           ::          cx, sx, cy, sy
      real(kind=real64)           ::          d2, d1

      real(kind=real64)           ::      ddat

      d1 = 2*PI/nx
      d2 = 2*PI/ny

      dat = 0

      do jy = 0, ny - 1
      do iy = 0, nFourierCoeffs(2)
         qy = d2*iy*jy
         sy = sin(qy)
         cy = cos(qy)

         do jx = 0, nx - 1
         do ix = 0, nFourierCoeffs(1)
            qx = d1*ix*jx
            sx = sin(qx)
            cx = cos(qx)

            ddat = ((qdat(1, ix, iy)*cx &
                     + qdat(2, ix, iy)*sx)*cy &
                    + (qdat(3, ix, iy)*cx &
                       + qdat(4, ix, iy)*sx)*sy)

            dat(jx, jy) = dat(jx, jy) + ddat

         end do
         end do
      end do
      end do

      return
   end function interpolateFromFourierCoefficients220

   function interpolateFromFourierCoefficients31(x, xmax, nFourierCoeffs, qdat) result(dat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(3), intent(in)       ::      x
      real(kind=real64), dimension(3), intent(in)       ::      xmax
      integer, dimension(3), intent(in)                 ::      nFourierCoeffs
      real(kind=real64), dimension(:, :, 0:, 0:, 0:), intent(in)      ::      qdat                 !   (1:d,1:8,0:nFourierCoeffs,0:nFourierCoeffs,0:nFourierCoeffs)
      real(kind=real64), dimension(size(qdat, dim=1))   ::      dat

      integer                     ::          ix, iy, iz
      real(kind=real64)           ::          qx, qy, qz

      integer                     ::          dd
      real(kind=real64)           ::          cx, sx, cy, sy, cz, sz
      real(kind=real64)           ::          d3, d2, d1

      dd = size(dat)

      d1 = x(1)*2*PI/xmax(1)
      d2 = x(2)*2*PI/xmax(2)
      d3 = x(3)*2*PI/xmax(3)

      dat = 0

      do iz = 0, nFourierCoeffs(3)
         qz = d3*iz
         sz = sin(qz)
         cz = cos(qz)

         do iy = 0, nFourierCoeffs(2)
            qy = d2*iy
            sy = sin(qy)
            cy = cos(qy)

            do ix = 0, nFourierCoeffs(1)
               qx = d1*ix
               sx = sin(qx)
               cx = cos(qx)

               dat(1:dd) = dat(1:dd) + ((qdat(1:dd, 1, ix, iy, iz)*cx &
                                         + qdat(1:dd, 2, ix, iy, iz)*sx)*cy &
                                        + (qdat(1:dd, 3, ix, iy, iz)*cx &
                                           + qdat(1:dd, 4, ix, iy, iz)*sx)*sy)*cz &
                           + ((qdat(1:dd, 5, ix, iy, iz)*cx &
                               + qdat(1:dd, 6, ix, iy, iz)*sx)*cy &
                              + (qdat(1:dd, 7, ix, iy, iz)*cx &
                                 + qdat(1:dd, 8, ix, iy, iz)*sx)*sy)*sz

            end do
         end do
      end do

      return
   end function interpolateFromFourierCoefficients31

   function interpolateFromFourierCoefficients32(nx, ny, nz, nFourierCoeffs, qdat) result(dat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      interpolate at each point of a regularly spaced lattice
      integer, intent(in)                              ::      nx, ny, nz
      integer, dimension(3), intent(in)                 ::      nFourierCoeffs
      real(kind=real64), dimension(:, :, 0:, 0:, 0:), intent(in)      ::      qdat                 !   (1:d,1:8,0:nFourierCoeffs,0:nFourierCoeffs,0:nFourierCoeffs)
      real(kind=real64), dimension(size(qdat, dim=1), 0:nx - 1, 0:ny - 1, 0:nz - 1)   ::      dat

      integer                     ::          ix, iy, iz
      integer                     ::          jx, jy, jz
      real(kind=real64)           ::          qx, qy, qz

      integer                     ::          dd
      real(kind=real64)           ::          cx, sx, cy, sy, cz, sz
      real(kind=real64)           ::          d3, d2, d1

      real(kind=real64), dimension(size(qdat, dim=1))       ::      ddat

      dd = size(qdat, dim=1)

      d1 = 2*PI/nx
      d2 = 2*PI/ny
      d3 = 2*PI/nz

      dat = 0

      do jz = 0, nz - 1
      do iz = 0, nFourierCoeffs(3)
         qz = d3*iz*jz
         sz = sin(qz)
         cz = cos(qz)

         do jy = 0, ny - 1
         do iy = 0, nFourierCoeffs(2)
            qy = d2*iy*jy
            sy = sin(qy)
            cy = cos(qy)

            do jx = 0, nx - 1
            do ix = 0, nFourierCoeffs(1)
               qx = d1*ix*jx
               sx = sin(qx)
               cx = cos(qx)

               ddat(1:dd) = ((qdat(1:dd, 1, ix, iy, iz)*cx &
                              + qdat(1:dd, 2, ix, iy, iz)*sx)*cy &
                             + (qdat(1:dd, 3, ix, iy, iz)*cx &
                                + qdat(1:dd, 4, ix, iy, iz)*sx)*sy)*cz &
                            + ((qdat(1:dd, 5, ix, iy, iz)*cx &
                                + qdat(1:dd, 6, ix, iy, iz)*sx)*cy &
                               + (qdat(1:dd, 7, ix, iy, iz)*cx &
                                  + qdat(1:dd, 8, ix, iy, iz)*sx)*sy)*sz

               dat(1:dd, jx, jy, jz) = dat(1:dd, jx, jy, jz) + ddat(1:dd)

            end do
            end do
         end do
         end do
      end do
      end do

      return
   end function interpolateFromFourierCoefficients32

end module Lib_LowPassFilter3d

