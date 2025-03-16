
module Lib_GaussianBlurs
!---^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_GaussianBlurs from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      perform a pixel-by-pixel blur on a 2D dataset
!*      ignore pixels labelled by GAUSSIANBLUR_IGNORE
   use iso_fortran_env
   implicit none
   private

   !---

   public      ::      GaussianBlur_ctor
   public      ::      delete
   public      ::      report

   public      ::      blur

   !---

#ifdef DEBUG
   logical, private, parameter          ::      GaussianBlur_dbg = .true.
#else
   logical, private, parameter          ::      GaussianBlur_dbg = .false.
#endif

   integer(kind=int64), private, parameter               ::      BADF00D = int(z'BADF00D', kind=int64)
        real(kind=real64),public,parameter                  ::      GAUSSIANBLUR_IGNORE = transfer( (BADF00D+ishft(BADF00D,32_int64)),1.0d0 )

   !---

   type, public     ::      GaussianBlur
      private
      real(kind=real64)               ::      sigma           !   width of blur in pixel units
      integer                         ::      m               !   size of blur region
      real(kind=real64), dimension(:, :), pointer    ::      g   !   [-m:m,-m:m]
   end type GaussianBlur

   !---

   interface GaussianBlur_ctor
      module procedure GaussianBlur_null
      module procedure GaussianBlur_ctor0
   end interface

   interface delete
      module procedure delete0
   end interface

   interface report
      module procedure report0
   end interface

   interface blur
      module procedure blur0
      module procedure blur1
   end interface

contains
!---^^^^^^^^

   function GaussianBlur_null() result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(GaussianBlur)           ::      this
      this%sigma = 0.0d0
      this%m = 0
      nullify (this%g)
      return
   end function GaussianBlur_null

   function GaussianBlur_ctor0(sigma) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(GaussianBlur)           ::      this
      real(kind=real64), intent(in) ::      sigma
      integer             ::      ix, iy
      real(kind=real64)   ::      i2s2, dd
      this%sigma = abs(sigma)
      this%m = ceiling(3*sigma)
      allocate (this%g(-this%m:this%m, -this%m:this%m))
      if (this%m == 0) then
         this%g = 1.0d0      !   no blur
      else
         i2s2 = 1/(2*this%sigma*this%sigma)
         do iy = -this%m, this%m
            do ix = -this%m, this%m
               dd = real(ix*ix + iy*iy, kind=real64)
               dd = dd*i2s2
               if (dd > 9.0d0) then
                  this%g(ix, iy) = 0.0d0
               else
                  this%g(ix, iy) = exp(-dd)
               end if
            end do
         end do
      end if
      return
   end function GaussianBlur_ctor0

   !---

   subroutine delete0(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^
      type(GaussianBlur), intent(inout)    ::      this
      if (this%m == 0) return
      deallocate (this%g)
      this = GaussianBlur_null()
      return
   end subroutine delete0

   !---

   subroutine report0(this, u, o)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(GaussianBlur), intent(in)    ::      this
      integer, intent(in), optional     ::      u, o
      integer     ::      uu, oo
      uu = 6; if (present(u)) uu = u
      oo = 0; if (present(o)) oo = o
      write (unit=uu, fmt='(a,f16.4,a)') repeat(" ", oo)//"GaussianBlur [sigma=", this%sigma, "]"
      return
   end subroutine report0

   !---

   function blur0(this, img, ix, iy) result(gimg)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return gaussian blur of 2d img centred on ix,iy
      type(GaussianBlur), intent(in)    ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::  img
      integer, intent(in)              ::      ix, iy
      real(kind=real64)               ::  gimg
      integer                 ::      nx, ny, kx, ky
      real(kind=real64)       ::      ff, gg, gsum
      nx = size(img, dim=1)
      ny = size(img, dim=2)
      gimg = 0.0d0; gsum = 0.0d0
      do ky = max(0, iy - this%m), min(ny - 1, iy + this%m)
         do kx = max(0, ix - this%m), min(nx - 1, ix + this%m)
            ff = img(kx, ky)
            if (ff == GAUSSIANBLUR_IGNORE) cycle
            gg = this%g(kx - ix, ky - iy)
            gsum = gsum + gg
            gimg = gimg + ff*gg
         end do
      end do
      if (gsum > 0) then
         gimg = gimg/gsum
      else
         gimg = GAUSSIANBLUR_IGNORE
      end if
      return
   end function blur0

   function blur1(this, img, ix, iy, pbc) result(gimg)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return gaussian blur of 2d img centred on ix,iy
      !*      periodic boundary conditions
      type(GaussianBlur), intent(in)    ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::  img
      integer, intent(in)              ::      ix, iy
      logical, intent(in)              ::      pbc
      real(kind=real64)               ::  gimg
      integer                 ::      nx, ny, kx, ky, jx, jy
      real(kind=real64)       ::      ff, gg, gsum

      if (.not. pbc) then
         gimg = blur0(this, img, ix, iy)
         return
      end if

      nx = size(img, dim=1)
      ny = size(img, dim=2)
      gimg = 0.0d0; gsum = 0.0d0
      do jy = -this%m, this%m
         ky = mod(iy + jy + ny, ny)
         do jx = -this%m, this%m
            kx = mod(ix + jx + nx, nx)
            ff = img(kx, ky)
            if (ff == GAUSSIANBLUR_IGNORE) cycle
            gg = this%g(jx, jy)
            gsum = gsum + gg
            gimg = gimg + ff*gg
         end do
      end do
      if (gsum > 0) then
         gimg = gimg/gsum
      else
         gimg = GAUSSIANBLUR_IGNORE
      end if
      return
   end function blur1

end module Lib_GaussianBlurs

