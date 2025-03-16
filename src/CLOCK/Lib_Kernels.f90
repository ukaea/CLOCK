
module Lib_kernels
!---^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    AnalysePngHistogram from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      very simple interface for dealing with kernels in 2 or 3d

   use iso_fortran_env
   implicit none
   private

   !---

   public      ::      Kernel_ctor
   public      ::      GaussianKernel_ctor
   public      ::      dGaussianKernel_ctor
   public      ::      report
   public      ::      delete

   public      ::      convolve
   public      ::      get
   public      ::      getNxNyNz

   !---

   type, public     ::      kernel
      private
      integer                                         ::      nx, ny, nz            !   size of kernel -nx:nx
      real(kind=real64), dimension(:, :, :), pointer      ::      k                   !   (-nx:nx,-ny:ny,-nz:nz)
   end type

   !---

   interface Kernel_ctor
      module procedure Kernel_null
      module procedure Kernel_ctor0
   end interface

   interface delete
      module procedure delete0
   end interface

   interface report
      module procedure report0
   end interface

   interface GaussianKernel_ctor
      module procedure GaussianKernel_ctor0
   end interface

   interface dGaussianKernel_ctor
      module procedure dGaussianKernel_ctor0
   end interface

   interface convolve
      module procedure convolve0
      module procedure convolve2
      module procedure convolve0a
      module procedure convolve0_2d
      module procedure convolve2a
      module procedure convolve2b
   end interface

   interface get
      module procedure get0
   end interface

contains
!---^^^^^^^^

   function Kernel_null() result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(Kernel)            ::      this
      this%nx = 0
      this%ny = 0
      this%nz = 0
      nullify (this%k)
      return
   end function Kernel_null

   function Kernel_ctor0(nx, ny, nz) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      integer, intent(in)      ::      nx, ny, nz
      type(Kernel)            ::      this
      this%nx = nx
      this%ny = ny
      this%nz = nz
      allocate (this%k(-this%nx:this%nx, -this%ny:this%ny, -this%nz:this%nz))
      this%k = 0.0d0
      this%k(0, 0, 0) = 1.0d0
      return
   end function Kernel_ctor0

   !---

   subroutine delete0(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^
      type(Kernel), intent(inout)   ::      this
      if (this%nx + this%ny == 0) return
      deallocate (this%k)
      this = Kernel_null()
      return
   end subroutine delete0

   subroutine report0(this, u, o)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(Kernel), intent(in)                 ::      this
      integer, intent(in), optional             ::      u, o
      integer             ::      uu, oo
      uu = 6; if (present(u)) uu = u
      oo = 0; if (present(o)) oo = o
      write (unit=uu, fmt='(4(a,i4))') repeat(" ", oo)//"Kernel (nx,ny,nz=", this%nx, ",", this%ny, ",", this%nz, ")"
      return
   end subroutine report0

   !---

   function GaussianKernel_ctor0(nx, ny, nz, sig) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      k(r) = exp( - r^2/2 sig^2 )
      integer, intent(in)              ::      nx, ny, nz
      real(kind=real64), intent(in)    ::      sig
      type(Kernel)                    ::      this
      integer                 ::      ix, iy, iz
      real(kind=real64)       ::      ww, yy, i2sig2
      this = Kernel_ctor0(nx, ny, nz)
      i2sig2 = sig*sig

      if (i2sig2 > 0) then
         i2sig2 = 0.5d0/i2sig2
      else
         return      !   already set k(0,0,0) = 1
      end if

      yy = 0.0d0
      do iz = -this%nz, this%nz
         do iy = -this%ny, this%ny
            do ix = -this%nx, this%nx
               ww = ix*ix + iy*iy + iz*iz
               ww = exp(-ww*i2sig2)
               yy = yy + ww
               this%k(ix, iy, iz) = ww
            end do
         end do
      end do

      yy = 1.0d0/yy
      this%k(-this%nx:this%nx, -this%ny:this%ny, -this%nz:this%nz) = this%k(-this%nx:this%nx, -this%ny:this%ny, -this%nz:this%nz)*yy

      return
   end function GaussianKernel_ctor0

   function dGaussianKernel_ctor0(nx, ny, nz, sig) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      k(r) = exp( - r^2/2 sig^2 )
      !*      note: constructs derivative k'(r)
      integer, intent(in)              ::      nx, ny, nz
      real(kind=real64), intent(in)    ::      sig
      type(Kernel)                    ::      this
      integer                 ::      ix, iy, iz
      real(kind=real64)       ::      rr, ww, yy, i2sig2
      this = Kernel_ctor0(nx, ny, nz)
      i2sig2 = sig*sig

      if (i2sig2 > 0) then
         i2sig2 = 0.5d0/i2sig2
      else
         return      !   already set k(0,0,0) = 1
      end if

      yy = 0.0d0
      do iz = -this%nz, this%nz
         do iy = -this%ny, this%ny
            do ix = -this%nx, this%nx
               rr = ix*ix + iy*iy + iz*iz
               ww = exp(-rr*i2sig2)
               ww = -2*sqrt(rr)*i2sig2*ww
               yy = yy + ww
               this%k(ix, iy, iz) = ww
            end do
         end do
      end do

      yy = 1.0d0/yy
      this%k(-this%nx:this%nx, -this%ny:this%ny, -this%nz:this%nz) = this%k(-this%nx:this%nx, -this%ny:this%ny, -this%nz:this%nz)*yy

      return
   end function dGaussianKernel_ctor0

   !---

   subroutine convolve0(this, img_in, img_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      performs a convolution with kernel
      !*          img_out(r) = int img_in(r') k(r'-r) dr' / int  k(r'-r) dr'
      !*      removes dead pixels
      type(Kernel), intent(in)                             ::      this
      real(kind=real64), dimension(:, :, :), intent(in)       ::      img_in
      real(kind=real64), dimension(:, :, :), intent(out)      ::      img_out

      integer                 ::      mx, my, mz
      integer                 ::      ix, iy, iz
      integer                 ::      jx, jy, jz, kx, ky, kz
      real(kind=real64)       ::      ksum, gg, hh

      mx = size(img_in, dim=1)
      my = size(img_in, dim=2)
      mz = size(img_in, dim=3)
      img_out(1:mx, 1:my, 1:mz) = 0.0d0
      do iz = 1, mz
         do iy = 1, my
            do ix = 1, mx
               if (img_in(ix, iy, iz) == 0) cycle

               ksum = 0.0d0; hh = 0.0d0
               do jz = max(1, iz - this%nz), min(mz, iz + this%nz)
                  kz = jz - iz
                  do jy = max(1, iy - this%ny), min(my, iy + this%ny)
                     ky = jy - iy
                     do jx = max(1, ix - this%nx), min(mx, ix + this%nx)
                        kx = jx - ix
                        gg = img_in(jx, jy, jz)
                        if (gg > 0) ksum = ksum + this%k(kx, ky, kz)
                        hh = hh + gg*this%k(kx, ky, kz)
                     end do
                  end do
               end do
               if (ksum > 0) img_out(ix, iy, iz) = hh/ksum
            end do
         end do
      end do
      return
   end subroutine convolve0

   subroutine convolve0_2d(this, img_in, img_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      performs a convolution with kernel
      !*          img_out(r) = int img_in(r') k(r'-r) dr' / int  k(r'-r) dr'
      !*      removes dead pixels
      type(Kernel), intent(in)                           ::      this
      real(kind=real64), dimension(:, :), intent(in)       ::      img_in
      real(kind=real64), dimension(:, :), intent(out)      ::      img_out

      integer                 ::      mx, my
      integer                 ::      ix, iy
      integer                 ::      jx, jy, kx, ky
      real(kind=real64)       ::      ksum, gg, hh

      mx = size(img_in, dim=1)
      my = size(img_in, dim=2)
!            img_out(1:mx,1:my) = 0.0d0
      do iy = 1, my
         do ix = 1, mx
            ksum = 0.0d0; hh = 0.0d0
            if (img_in(ix, iy) /= 0) then
               do jy = max(1, iy - this%ny), min(my, iy + this%ny)
                  ky = jy - iy
                  do jx = max(1, ix - this%nx), min(mx, ix + this%nx)
                     kx = jx - ix
                     gg = img_in(jx, jy)
                     if (gg > 0) then
                        ksum = ksum + this%k(kx, ky, 0)
                        hh = hh + gg*this%k(kx, ky, 0)
                     end if
                  end do
               end do
            end if
            if (ksum /= 0) ksum = 1/ksum
            img_out(ix, iy) = hh*ksum
         end do
      end do
      return
   end subroutine convolve0_2d

   subroutine convolve2(this, f, img_in, img_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      performs a convolution with kernel with extra multiplicative function
      !*          img_out(r) = int img_in(r') f(r') k(r'-r) dr'    / int f(r') k(r'-r) dr'
      !*      removes dead pixels
      type(Kernel), intent(in)                             ::      this
      real(kind=real64), dimension(:, :, :), intent(in)       ::      img_in, f
      real(kind=real64), dimension(:, :, :), intent(out)      ::      img_out

      integer                 ::      mx, my, mz
      integer                 ::      ix, iy, iz
      integer                 ::      jx, jy, jz, kx, ky, kz
      real(kind=real64)       ::      ksum, gg, hh

      mx = size(img_in, dim=1)
      my = size(img_in, dim=2)
      mz = size(img_in, dim=3)
      img_out(1:mx, 1:my, 1:mz) = 0.0d0

      do iz = 1, mz
         do iy = 1, my
            do ix = 1, mx
               if (img_in(ix, iy, iz) == 0) cycle

               ksum = 0.0d0; hh = 0.0d0
               do jz = max(1, iz - this%nz), min(mz, iz + this%nz)
                  kz = jz - iz
                  do jy = max(1, iy - this%ny), min(my, iy + this%ny)
                     ky = jy - iy
                     do jx = max(1, ix - this%nx), min(mx, ix + this%nx)
                        kx = jx - ix
                        gg = img_in(jx, jy, jz)
                        if (gg > 0) ksum = ksum + this%k(kx, ky, kz)*f(jx, jy, jz)
                        hh = hh + gg*this%k(kx, ky, kz)*f(jx, jy, jz)
                     end do
                  end do
               end do
               if (ksum > 0) img_out(ix, iy, iz) = hh/ksum
            end do
         end do
      end do
      return
   end subroutine convolve2

   !---

   subroutine convolve0a(this, img_in, ix, iy, iz, img_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      performs a convolution with kernel
      !*          img_out(r) = int img_in(r') k(r'-r) dr'  / int f(r') k(r'-r) dr'
      !*      removes dead pixels
      type(Kernel), intent(in)                             ::      this
      real(kind=real64), dimension(:, :, :), intent(in)       ::      img_in
      integer, intent(in)                                  ::      ix, iy, iz
      real(kind=real64), intent(out)                       ::      img_out

      integer                 ::      mx, my, mz
      integer                 ::      jx, jy, jz, kx, ky, kz
      real(kind=real64)       ::      ksum, gg, hh

      mx = size(img_in, dim=1)
      my = size(img_in, dim=2)
      mz = size(img_in, dim=3)
      img_out = 0.0d0
      if (img_in(ix, iy, iz) == 0) return

      ksum = 0.0d0; hh = 0.0d0
      do jz = max(1, iz - this%nz), min(mz, iz + this%nz)
         kz = jz - iz
         do jy = max(1, iy - this%ny), min(my, iy + this%ny)
            ky = jy - iy
            do jx = max(1, ix - this%nx), min(mx, ix + this%nx)
               kx = jx - ix
               gg = img_in(jx, jy, jz)
               if (gg > 0) ksum = ksum + this%k(kx, ky, kz)
               hh = hh + gg*this%k(kx, ky, kz)
            end do
         end do
      end do
      if (ksum > 0) img_out = hh/ksum
      return
   end subroutine convolve0a

   subroutine convolve2a(this, f, img_in, ix, iy, iz, img_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      performs a convolution with kernel with extra multiplicative function
      !*          img_out(r) = int img_in(r') f(r') k(r'-r) dr'     / int f(r') k(r'-r) dr'
      !*      removes dead pixels
      type(Kernel), intent(in)                             ::      this
      real(kind=real64), dimension(:, :, :), intent(in)       ::      img_in, f
      integer, intent(in)                                  ::      ix, iy, iz
      real(kind=real64), intent(out)                       ::      img_out

      integer                 ::      mx, my, mz
      integer                 ::      jx, jy, jz, kx, ky, kz
      real(kind=real64)       ::      ksum, gg, hh

      mx = size(img_in, dim=1)
      my = size(img_in, dim=2)
      mz = size(img_in, dim=3)
      img_out = 0.0d0
      if (img_in(ix, iy, iz) == 0) return
      ksum = 0.0d0; hh = 0.0d0

      do jz = max(1, iz - this%nz), min(mz, iz + this%nz)
         kz = jz - iz
         do jy = max(1, iy - this%ny), min(my, iy + this%ny)
            ky = jy - iy
            do jx = max(1, ix - this%nx), min(mx, ix + this%nx)
               kx = jx - ix
               gg = img_in(jx, jy, jz)
               if (gg > 0) ksum = ksum + this%k(kx, ky, kz)*f(jx, jy, jz)
               hh = hh + gg*this%k(kx, ky, kz)*f(jx, jy, jz)
            end do
         end do
      end do
      if (ksum > 0) img_out = hh/ksum
      return
   end subroutine convolve2a

   !---

   subroutine convolve2b(this, f, img_in, ix, iy, iz, img_out, img2_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      performs a convolution with kernel with extra multiplicative function
      !*          img_out(r) = int img_in(r') f(r') k(r'-r) dr'    / int f(r') k(r'-r) dr'
      !*         and           int img_in(r')^2 f(r') k(r'-r) dr'  / int f(r') k(r'-r) dr'
      !*      removes dead pixels
      type(Kernel), intent(in)                             ::      this
      real(kind=real64), dimension(:, :, :), intent(in)       ::      img_in, f
      integer, intent(in)                                  ::      ix, iy, iz
      real(kind=real64), intent(out)                       ::      img_out, img2_out

      integer                 ::      mx, my, mz
      integer                 ::      jx, jy, jz, kx, ky, kz
      real(kind=real64)       ::      ksum, gg, hh, h2

      mx = size(img_in, dim=1)
      my = size(img_in, dim=2)
      mz = size(img_in, dim=3)
      img_out = 0.0d0
      img2_out = 0.0d0
      if (img_in(ix, iy, iz) == 0) return

      ksum = 0.0d0; hh = 0.0d0; h2 = 0.0d0

      do jz = max(1, iz - this%nz), min(mz, iz + this%nz)
         kz = jz - iz
         do jy = max(1, iy - this%ny), min(my, iy + this%ny)
            ky = jy - iy
            do jx = max(1, ix - this%nx), min(mx, ix + this%nx)
               kx = jx - ix
               gg = img_in(jx, jy, jz)
               if (gg > 0) ksum = ksum + this%k(kx, ky, kz)*f(jx, jy, jz)
               hh = hh + gg*this%k(kx, ky, kz)*f(jx, jy, jz)
               h2 = h2 + gg*gg*this%k(kx, ky, kz)*f(jx, jy, jz)
            end do
         end do
      end do
      if (ksum > 0) then
         img_out = hh/ksum
         img2_out = h2/ksum
      end if
      return
   end subroutine convolve2b

   pure function get0(this, ix, iy, iz) result(k)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(Kernel), intent(in)         ::      this
      integer, intent(in)              ::      ix, iy, iz
      real(kind=real64)               ::      k
      k = this%k(ix, iy, iz)
      return
   end function get0

   pure subroutine getNxNyNz(this, nx, ny, nz)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(Kernel), intent(in)         ::      this
      integer, intent(out)             ::      nx, ny, nz
      nx = this%nx
      ny = this%ny
      nz = this%nz
      return
   end subroutine getNxNyNz

end module Lib_kernels

