
module Lib_LaplacianOfGaussians
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_LaplacianOfGaussians from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      computes the scale invariant laplacian of gaussians
!*
!*      write the discrete Gaussian
!*      with width t as
!*          g(x,y;t) = 1/Z exp[ -(x^2 + y^2)/(2 t^2) ]
!*      where sum_(x,y) g(x,y;t) = 1 defines the normalising constant Z
!*
!*      define a convolution with this Gaussian on a function f(x,y)
!*          C(x,y;t) = g(x,y;t) * f(x,y)
!*                   = 1/Z_(x,y) sum_{x',y'} g(x',y';t) f( x-x',y-y' )
!*      note that Z is a function of x,y here- the image may have finite extent, and so only f( x-x',y-y' ) in range can be sampled
!*
!*      the scale-invariant Laplacian is then
!*          L(x,y,t) = t^2 ( d^2/dx^2 C + d^2/dy^2 C )
!*      which can be computed using finite difference coefficients.
!*
!*      bright blobs radius t are local minima of L, dark blobs radius t are local maxima.
!*
!        use iso_c_binding
   use iso_fortran_env
!        use Lib_SimpleProgressBar
!        use Lib_FFTW3f
   use OMP_LIB
!        use Lib_RidlerCalvard
!        use Lib_png
   implicit none

   !include "fftw3.f03"
   private

   integer, public, parameter        ::      LIB_LoG_DARK = -1
   integer, public, parameter        ::      LIB_LoG_BRIGHT = +1
   logical, public, parameter        ::      LIB_LoG_DOG = .false.

#ifdef DEBUG
   logical, private, parameter         ::      LIB_LoG_DBG = .true.
#else
   logical, private, parameter          ::      LIB_LoG_DBG = .false.
#endif

   public          ::      scaleInvariantLaplacianOfGaussians
   public          ::      blobFinder
   public          ::      gaussianBlurImage
   public      ::      oneShotLoG
   public      ::      laplacianOfGaussian
   !public      ::      ninePointStencil

contains
!---^^^^^^^^

   subroutine blobFinder(f_in, t, nBlob, datBlob, bright, dark, quiet)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find blobs in image
      !*      test radius range t_min:t_max with nt steps
      !*      return output datBlob(1:6,1:nBlob) = (/ x,y,r,i,d,LoG /)      with d = +1 for bright and -1 for dark. note x = [0:Nx-1]

      real(kind=real64), dimension(0:, 0:), intent(in)               ::      f_in
      real(kind=real64), dimension(:), intent(in)                   ::      t
      integer, intent(out)                                         ::      nBlob
      real(kind=real64), dimension(:, :), allocatable, intent(out)    ::      datBlob
      logical, intent(in)                                          ::      bright, dark
      logical, intent(in), optional                                 ::      quiet

      real(kind=real64), dimension(:, :, :), pointer      ::      LL, f_blur
      integer             ::      Nx, Ny, Nt, NbMax
      integer             ::      ix, iy, it, sqrtNt
      real(kind=real64), dimension(:, :), pointer                  ::      bd, bd_tmp
      real(kind=real64), dimension(-1:1, -1:1, -1:1)     ::      L27
      real(kind=real64)   ::      L0, L1, L2
      real(kind=real64)   ::      ff, f0, xx, yy, tt
      logical             ::      qq

      !---    find scale-invariant stack of Laplacians of Gaussians
      Nx = size(f_in, dim=1)
      Ny = size(f_in, dim=2)
      Nt = size(t)

      sqrtNt = ceiling(sqrt(real(Nt)))
      qq = .false.; if (present(quiet)) qq = quiet

      allocate (LL(0:Nx - 1, 0:Ny - 1, 1:Nt))
      allocate (f_blur(0:Nx - 1, 0:Ny - 1, 1:Nt))

      call scaleInvariantLaplacianOfGaussians(f_in, t, LL, f_blur, qq)

      NbMax = 100
      allocate (bd(6, NbMax))
      nullify (bd_tmp)
      nBlob = 0

      do it = 2, Nt - 1

         do iy = 1, Ny - 2
            do ix = 1, Nx - 2

               ff = f_blur(ix, iy, it)

               if (dark .or. bright) then

                  L27(-1:1, -1:1, -1:1) = LL(ix - 1:ix + 1, iy - 1:iy + 1, it - 1:it + 1)

                  L0 = L27(0, 0, 0)

                  if (dark) then
                     L27(0, 0, 0) = -huge(1.0)
                     L1 = maxval(L27)
                     if ((L0 > L1)) then

                        !  SIFT algorithm

                        L27(0, 0, 0) = L0
                        call extremum(L27, f0, xx, yy, tt)

                        if (any(abs((/xx, yy, tt/)) > 0.5)) cycle                                              !   actually closer to the next pixel.
                        if (((ix - 1)*(ix - (Nx - 2))*(iy - 1)*(iy - (Ny - 2)) == 0) .and. (abs(f0) < 0.1)) cycle        !   on the edge and not a very strong signal. probably dodgy.
                        if ((abs(f0) < 0.01)) cycle                                                         !   very weak signal. probably dodgy.

                        nBlob = nBlob + 1
                        bd(1, nBlob) = ix + xx
                        bd(2, nBlob) = iy + yy
                        if (tt > 0) then
                           bd(3, nBlob) = t(it)*(1 - tt) + t(it + 1)*tt
                        else
                           bd(3, nBlob) = t(it)*(1 - tt) + t(it - 1)*tt
                        end if

                        bd(4, nBlob) = ff

                        bd(5, nBlob) = LIB_LoG_DARK

                        bd(6, nBlob) = f0

!

                     end if
                  end if
                  if (bright) then
                     L27(0, 0, 0) = +huge(1.0)
                     L2 = minval(L27)

                     if ((L0 < L2)) then

                        L27(0, 0, 0) = L0
                        call extremum(L27, f0, xx, yy, tt)
                        if (any(abs((/xx, yy, tt/)) > 0.5)) cycle                                               !   actually closer to the next pixel.
                        if (((ix - 1)*(ix - (Nx - 2))*(iy - 1)*(iy - (Ny - 2)) == 0) .and. (abs(f0) < 0.1)) cycle        !   on the edge and not a very strong signal. probably dodgy.
                        if ((abs(f0) < 0.01)) cycle                                                         !   very weak signal. probably dodgy.

                        if (LIB_LoG_DBG) &
                           print *, "bright blob ", ix, iy, it, L0, L2, f0, xx, yy, tt, ff ! , thresh(it)

                        nBlob = nBlob + 1
                        bd(1, nBlob) = ix + xx
                        bd(2, nBlob) = iy + yy
                        if (tt > 0) then
                           bd(3, nBlob) = t(it)*(1 - tt) + t(it + 1)*tt
                        else
                           bd(3, nBlob) = t(it)*(1 - tt) + t(it - 1)*tt
                        end if
                        bd(4, nBlob) = ff
                        bd(5, nBlob) = LIB_LoG_BRIGHT
                        bd(6, nBlob) = f0
                        !if (nBlob==10) stop
                     end if
                  end if
               end if

               if (nBlob == NbMax) then
                  allocate (bd_tmp(6, NbMax*2))
                  bd_tmp(1:6, 1:nBlob) = bd(1:6, 1:nBlob)
                  deallocate (bd)
                  bd => bd_tmp
                  NbMax = NbMax*2
                  !stop
               end if

            end do
         end do
      end do

      allocate (datBlob(6, nBlob))
      datBlob(1:6, 1:nBlob) = bd(1:6, 1:nBlob)
      deallocate (bd)

      if (LIB_LoG_DBG) print *, "blobFinder info - nBlob = ", nBlob

      return
   end subroutine blobFinder

   subroutine gaussianBlurImage(f_in, t, f_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      produce a Gaussian blur of input image, width t
      real(kind=real64), dimension(0:, 0:), intent(in)       ::      f_in
      real(kind=real64), dimension(0:, 0:), intent(inout)    ::      f_out
      real(kind=real64), intent(in)                        ::      t

      real(kind=real64), dimension(:), allocatable      ::      kernel, f_stripe, w_stripe

      integer             ::      Nx, Ny
      integer             ::      ix, iy, jj, kk
      integer             ::      Nk
      real(kind=real64)   ::      i2s2, ww, wf, ws!,wx

      Nx = size(f_in, dim=1)
      Ny = size(f_in, dim=2)

      !---    compute the unnormalised kernel
      Nk = max(5, ceiling(t*5))             !   range of pixels to search is +/- Nk
      allocate (kernel(0:Nk))
      i2s2 = 1/(2*t*t)
      kernel(0) = 1.0d0
      do ix = 1, Nk
         kernel(ix) = exp(-ix*ix*i2s2)
      end do

      !---    compute the output blurred image. First do x strips.
      allocate (f_stripe(0:Ny - 1))
      allocate (w_stripe(0:Nx - 1))

!$OMP PARALLEL  PRIVATE(ix,iy,wf,ws,jj,kk,ww,f_stripe,w_stripe ) SHARED( f_in,f_out,Nx,Ny,Nk,kernel)
!$OMP DO
      do ix = 0, Nx - 1
         wf = 0.0d0; ws = 0.0d0
         do jj = max(0, ix - Nk), min(Nx - 1, ix + Nk)
            kk = abs(jj - ix)
            ww = kernel(kk)
            wf = wf + ww*f_in(jj, 0)
            ws = ws + ww
         end do
         f_out(ix, 0) = wf
         w_stripe(ix) = 1/ws
      end do
!$OMP END DO

!$OMP DO
      do iy = 1, Ny - 1
         do ix = 0, Nx - 1
            wf = 0.0d0; ws = 0.0d0
            do jj = max(0, ix - Nk), min(Nx - 1, ix + Nk)
               kk = abs(jj - ix)
               ww = kernel(kk)
               wf = wf + ww*f_in(jj, iy)
               ws = ws + ww
            end do
            f_out(ix, iy) = wf
         end do
      end do
!$OMP END DO

      !   at this point, f_out has blurring in the x-direction only. weight stores the kernel weighting from this op.

      !---    now do y strips
!$OMP DO
      do ix = 0, Nx - 1
         !   make a copy of this stripe
         f_stripe(0:Ny - 1) = f_out(ix, 0:Ny - 1)
         do iy = 0, Ny - 1
            wf = 0.0d0; ws = 0.0d0
            do jj = max(0, iy - Nk), min(Ny - 1, iy + Nk)
               kk = abs(jj - iy)
               ww = kernel(kk)
               wf = wf + ww*f_stripe(jj)
               ws = ws + ww
            end do
            f_out(ix, iy) = wf*w_stripe(ix)/(ws)                !   note: ws /= 0 because there is always at least one pixel contributing. Several really.

         end do
      end do
!$OMP END DO

!$OMP END PARALLEL
      !---    now have a normalised Gaussian blur function

      return
   end subroutine gaussianBlurImage

   subroutine oneShotLoG(f_in, t, L_out, f_blur)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      compute the Laplacian of Gaussian exactly in one shot. Slow but accurate.
      real(kind=real64), dimension(0:, 0:), intent(in)       ::      f_in
      real(kind=real64), dimension(0:, 0:), intent(inout)      ::      L_out
      real(kind=real64), intent(in)                        ::      t

      real(kind=real64), dimension(0:, 0:), intent(inout)      ::      f_blur
      integer             ::      Nx, Ny, Nk
      integer             ::      ix, iy, jx, jy

      real(kind=real64), dimension(:, :), allocatable        ::      kernel
      real(kind=real64)   ::      is2, kk, ff, sumk, sumkf

      Nx = size(f_in, dim=1)
      Ny = size(f_in, dim=2)

      !---    compute the Laplacian of Gaussian kernel
      Nk = max(3, ceiling(t*3))             !   range of pixels to search is +/- Nk
      is2 = 1/(t*t)
      allocate (kernel(-Nk:Nk, -Nk:Nk)); kernel = 0.0d0

      do iy = -Nk, Nk
         do ix = -Nk, Nk
            kk = ix*ix + iy*iy
            kk = kk*is2
            if (kk < 9d0) then
               kk = exp(-kk/2)
               kernel(ix, iy) = (ix*ix + iy*iy - 2*t*t)*is2*is2*kk
            end if
         end do
      end do

      !---    now compute LoG at each point
      L_out = 0

      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            sumk = 0
            sumkf = 0
            do jy = max(0, iy - Nk), min(Ny - 1, iy + Nk)
               do jx = max(0, ix - Nk), min(Nx - 1, ix + Nk)
                  kk = kernel(jx - ix, jy - iy)
                  ff = f_in(jx, jy)
                  sumk = sumk + kk
                  sumkf = sumkf + kk*ff
               end do
            end do
            if (sumk > 0) L_out(ix, iy) = sumkf/sumk
         end do
      end do

      return

      return
   end subroutine oneShotLoG

   subroutine laplacianOfGaussian(f_in, t, L_out, f_blur)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*       produce a Laplacian of Gaussian blur of input image, width t
      real(kind=real64), dimension(0:, 0:), intent(in)       ::      f_in
      real(kind=real64), dimension(0:, 0:), intent(inout)    ::      L_out
      real(kind=real64), intent(in)                        ::      t

      real(kind=real64), dimension(0:, 0:), intent(inout)    ::      f_blur
      integer             ::      Nx, Ny
      integer             ::      ix, iy, jx, jy, kx, ky

      Nx = size(f_in, dim=1)
      Ny = size(f_in, dim=2)

      if (t > 0) then
         call gaussianBlurImage(f_in, t, f_blur)
      else
         f_blur = f_in
      end if

      do iy = 0, Ny - 1
         jy = iy - 2; ky = 0
         if (iy == 0) then
            jy = 0; ky = -2
         else if (iy == 1) then
            jy = 0; ky = -1
         else if (iy == Ny - 2) then
            jy = Ny - 5; ky = 1
         else if (iy == Ny - 1) then
            jy = Ny - 5; ky = 2
         end if
         do ix = 0, Nx - 1
            jx = ix - 2; kx = 0
            if (ix == 0) then
               jx = 0; kx = -2
            else if (ix == 1) then
               jx = 0; kx = -1
            else if (ix == Nx - 2) then
               jx = Nx - 5; kx = 1
            else if (ix == Nx - 1) then
               jx = Nx - 5; kx = 2
            end if
            L_out(ix, iy) = twentyFivePointLaplacian(f_blur(jx:jx + 4, jy:jy + 4), kx, ky)
         end do
      end do
      return

      !---    compute the Laplacian at the corners
    L_out(0, 0) = 4*f_blur(0, 0) - 5*(f_blur(1, 0) + f_blur(0, 1)) + 4*(f_blur(2, 0) + f_blur(0, 2)) - (f_blur(3, 0) + f_blur(0, 3))
            L_out(Nx-1,   0) = 4*f_blur(Nx-1,   0) - 5*(f_blur(Nx-2,   0)+f_blur(Nx-1,   1)) + 4*(f_blur(Nx-3,   0)+f_blur(Nx-1,   2)) - (f_blur(Nx-4,   0)+f_blur(Nx-1,   3))
            L_out(   0,Ny-1) = 4*f_blur(   0,Ny-1) - 5*(f_blur(   1,Ny-1)+f_blur(   0,Ny-2)) + 4*(f_blur(   2,Ny-1)+f_blur(   0,Ny-3)) - (f_blur(   3,Ny-1)+f_blur(   0,Ny-4))
            L_out(Nx-1,Ny-1) = 4*f_blur(Nx-1,Ny-1) - 5*(f_blur(Nx-2,Ny-1)+f_blur(Nx-1,Ny-2)) + 4*(f_blur(Nx-3,Ny-1)+f_blur(Nx-1,Ny-3)) - (f_blur(Nx-4,Ny-1)+f_blur(Nx-1,Ny-4))

      !---    compute the Laplacian on the edges
!$OMP PARALLEL DO SHARED(Nx,Ny,f_blur,L_out) PRIVATE (ix)
      do ix = 1, Nx - 2
         L_out(ix, 0) = f_blur(ix - 1, 0) + f_blur(ix + 1, 0) - 5*f_blur(ix, 1) + 4*f_blur(ix, 2) - f_blur(ix, 3)
         L_out(  ix,Ny-1) = f_blur(ix-1,Ny-1) + f_blur(ix+1,Ny-1) - 5*f_blur(  ix,Ny-2) + 4*f_blur(  ix,Ny-3) - f_blur(  ix,Ny-4)
      end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SHARED(Nx,Ny,f_blur,L_out) PRIVATE (iy)
      do iy = 1, Ny - 2
         L_out(0, iy) = f_blur(0, iy - 1) + f_blur(0, iy + 1) - 5*f_blur(1, iy) + 4*f_blur(2, iy) - f_blur(3, iy)
         L_out(Nx-1,  iy) = f_blur(Nx-1,iy-1) + f_blur(Nx-1,iy+1) - 5*f_blur(Nx-2,  iy) + 4*f_blur(Nx-3,  iy) - f_blur(Nx-4,  iy)
      end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SHARED(Nx,Ny,f_blur,L_out) PRIVATE (iy)
      do iy = 1, Ny - 2
         L_out(1, iy) = f_blur(0, iy) - 4*f_blur(1, iy) + f_blur(2, iy) + f_blur(1, iy - 1) + f_blur(1, iy + 1)
L_out(Nx - 2, iy) = f_blur(Nx - 3, iy) - 4*f_blur(Nx - 2, iy) + f_blur(Nx - 1, iy) + f_blur(Nx - 2, iy - 1) + f_blur(Nx - 2, iy + 1)
      end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SHARED(Nx,Ny,f_blur,L_out) PRIVATE (ix)
      do ix = 2, Ny - 3
         L_out(ix, 1) = f_blur(ix - 1, 1) - 4*f_blur(ix, 1) + f_blur(ix + 1, 1) + f_blur(ix, 0) + f_blur(ix, 2)
L_out(ix, Ny - 2) = f_blur(ix - 1, Ny - 2) - 4*f_blur(ix, Ny - 2) + f_blur(ix + 1, Ny - 2) + f_blur(ix, Ny - 3) + f_blur(ix, Ny - 1)
      end do
!$OMP END PARALLEL DO

      !---    compute the Laplacian in the interior

!$OMP PARALLEL DO SHARED(Nx,Ny,f_blur,L_out) PRIVATE (ix,iy)
      do iy = 2, Ny - 3
         do ix = 2, Nx - 3
     L_out(ix, iy) = (-f_blur(ix - 2, iy) + 16*f_blur(ix - 1, iy) - 60*f_blur(ix, iy) + 16*f_blur(ix + 1, iy) - f_blur(ix + 2, iy) &
                             - f_blur(ix, iy - 2) + 16*f_blur(ix, iy - 1) + 16*f_blur(ix, iy + 1) - f_blur(ix, iy + 2))/12

         end do
      end do
!$OMP END PARALLEL DO

      return
   end subroutine laplacianOfGaussian

   subroutine scaleInvariantLaplacianOfGaussians(f_in, t, L_out, f_blur, quiet)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      compute scale-invariant stack of laplacians.
      !*      note t(:) should be ordered in increasing size.
      real(kind=real64), dimension(0:, 0:), intent(in)       ::      f_in
      real(kind=real64), dimension(:), intent(in)           ::      t
      real(kind=real64), dimension(0:, 0:, :), intent(inout)  ::      L_out
      real(kind=real64), dimension(0:, 0:, :), intent(inout)  ::      f_blur
      logical, intent(in), optional                         ::      quiet

      real(kind=real64), dimension(:, :), allocatable        ::      tmp

      integer             ::      Nx, Ny, Nt, ix, iy, it
      real(kind=real64)   ::      tt
      logical             ::      qq

      Nx = size(f_in, dim=1)
      Ny = size(f_in, dim=2)
      Nt = size(t)
      allocate (tmp(0:Nx - 1, 0:Ny - 1))
      tt = 0
      tmp = f_in
      qq = .false.; if (present(quiet)) qq = quiet

   if (.not. qq) print *,"Lib_LaplacianOfGaussians::scaleInvariantLaplacianOfGaussians info - ",Nt," blurs on image size ",Nx,"x",Ny
            if (.not. qq) print *,"Lib_LaplacianOfGaussians::scaleInvariantLaplacianOfGaussians info - minmaxavg ",minval(f_in),maxval(f_in),sum(f_in)/(Nx*Ny)
      do it = 1, Nt

         tt = sqrt(t(it)*t(it) - tt*tt)        !   additional gaussian blur to get to total blur t.
                if (LIB_LoG_DBG) write(*,fmt='(a,i6,a,i6,a,f10.3,a,f10.3)') "Lib_LaplacianOfGaussians::scaleInvariantLaplacianOfGaussians info step ",it,"/",Nt," - add blur ",tt," total ",t(it)
         if (LIB_LoG_DoG) then
            call gaussianBlurImage(tmp, tt, f_blur(:, :, it))
         else
            call laplacianOfGaussian(tmp, tt, L_out(:, :, it), f_blur(:, :, it))
            L_out(:, :, it) = L_out(:, :, it)*t(it)*t(it)
         end if

         tt = t(it)
         tmp(0:Nx - 1, 0:Ny - 1) = f_blur(0:Nx - 1, 0:Ny - 1, it)

      end do

      if (LIB_LoG_DoG) then
         do iy = 0, Ny - 1
            do ix = 0, Nx - 1
               L_out(ix, iy, 1) = (-3*f_blur(ix, iy, 1) + 4*f_blur(ix, iy, 2) - f_blur(ix, iy, 3))/2
               do it = 2, Nt - 1
                  L_out(ix, iy, it) = (-f_blur(ix, iy, it - 1) + f_blur(ix, iy, it + 1))/2
               end do
               L_out(ix, iy, Nt) = (3*f_blur(ix, iy, Nt) - 4*f_blur(ix, iy, Nt - 1) + f_blur(ix, iy, Nt - 2))/2
            end do
         end do
      end if
      deallocate (tmp)
      return

   end subroutine scaleInvariantLaplacianOfGaussians

   subroutine twentySevenPointStencil(f, f0, dfdx, dfdy, dfdz, d2fdx2, d2fdy2, d2fdz2, d2fdxdy, d2fdydz, d2fdzdx)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(-1:1, -1:1, -1:1), intent(in)      ::      f
            real(kind=real64),intent(out)                               ::      f0, dfdx,dfdy,dfdz, d2fdx2,d2fdy2,d2fdz2, d2fdxdy, d2fdydz, d2fdzdx

      !                                                                   kernel shapes
      !                                                                   [mmm] [0mm] [pmm]    [m0m] [00m] [p0m]    [mpm] [0pm] [ppm]
      !                                                                   [mm0] [0m0] [pm0]    [m00] [000] [p00]    [mp0] [0p0] [pp0]
      !                                                                   [mmp] [0mp] [pmp]    [m0p] [00p] [p0p]    [mpp] [0pp] [ppp]

      !     kernels
      real(kind=real64), dimension(-1:1, -1:1, -1:1), parameter       ::      kernel_0 = reshape((/ &
                                                                                                 -2, 1, -2, 1, 4, 1, -2, 1, -2, &
                                                                                                 1, 4, 1, 4, 7, 4, 1, 4, 1, &
                                                                                -2, 1, -2, 1, 4, 1, -2, 1, -2/), (/3, 3, 3/))/27.0d0

      real(kind=real64), dimension(-1:1, -1:1, -1:1), parameter       ::      kernel_x = reshape((/ &
                                                                                                 -1, 0, 1, -1, 0, 1, -1, 0, 1, &
                                                                                                 -1, 0, 1, -1, 0, 1, -1, 0, 1, &
                                                                                 -1, 0, 1, -1, 0, 1, -1, 0, 1/), (/3, 3, 3/))/18.0d0
      real(kind=real64), dimension(-1:1, -1:1, -1:1), parameter       ::      kernel_y = reshape((/ &
                                                                                                 -1, -1, -1, 0, 0, 0, 1, 1, 1, &
                                                                                                 -1, -1, -1, 0, 0, 0, 1, 1, 1, &
                                                                                 -1, -1, -1, 0, 0, 0, 1, 1, 1/), (/3, 3, 3/))/18.0d0
      real(kind=real64), dimension(-1:1, -1:1, -1:1), parameter       ::      kernel_z = reshape((/ &
                                                                                               -1, -1, -1, -1, -1, -1, -1, -1, -1, &
                                                                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                                                                                    1, 1, 1, 1, 1, 1, 1, 1, 1/), (/3, 3, 3/))/18.0d0

      real(kind=real64), dimension(-1:1, -1:1, -1:1), parameter       ::      kernel_xx = reshape((/ &
                                                                                                  1, -2, 1, 1, -2, 1, 1, -2, 1, &
                                                                                                  1, -2, 1, 1, -2, 1, 1, -2, 1, &
                                                                                  1, -2, 1, 1, -2, 1, 1, -2, 1/), (/3, 3, 3/))/9.0d0
      real(kind=real64), dimension(-1:1, -1:1, -1:1), parameter       ::      kernel_yy = reshape((/ &
                                                                                                  1, 1, 1, -2, -2, -2, 1, 1, 1, &
                                                                                                  1, 1, 1, -2, -2, -2, 1, 1, 1, &
                                                                                  1, 1, 1, -2, -2, -2, 1, 1, 1/), (/3, 3, 3/))/9.0d0
      real(kind=real64), dimension(-1:1, -1:1, -1:1), parameter       ::      kernel_zz = reshape((/ &
                                                                                                  1, 1, 1, 1, 1, 1, 1, 1, 1, &
                                                                                               -2, -2, -2, -2, -2, -2, -2, -2, -2, &
                                                                                     1, 1, 1, 1, 1, 1, 1, 1, 1/), (/3, 3, 3/))/9.0d0

      real(kind=real64), dimension(-1:1, -1:1, -1:1), parameter       ::      kernel_xy = reshape((/ &
                                                                                                  1, 0, -1, 0, 0, 0, -1, 0, 1, &
                                                                                                  1, 0, -1, 0, 0, 0, -1, 0, 1, &
                                                                                  1, 0, -1, 0, 0, 0, -1, 0, 1/), (/3, 3, 3/))/12.0d0
      real(kind=real64), dimension(-1:1, -1:1, -1:1), parameter       ::      kernel_yz = reshape((/ &
                                                                                                  1, 1, 1, 0, 0, 0, -1, -1, -1, &
                                                                                                  0, 0, 0, 0, 0, 0, 0, 0, 0, &
                                                                                 -1, -1, -1, 0, 0, 0, 1, 1, 1/), (/3, 3, 3/))/12.0d0
      real(kind=real64), dimension(-1:1, -1:1, -1:1), parameter       ::      kernel_zx = reshape((/ &
                                                                                                  1, 0, -1, 1, 0, -1, 1, 0, -1, &
                                                                                                  0, 0, 0, 0, 0, 0, 0, 0, 0, &
                                                                                 -1, 0, 1, -1, 0, 1, -1, 0, 1/), (/3, 3, 3/))/12.0d0

      !                                                                   kernel shapes
      !                                                                   [mmm] [0mm] [pmm]    [m0m] [00m] [p0m]    [mpm] [0pm] [ppm]
      !                                                                   [mm0] [0m0] [pm0]    [m00] [000] [p00]    [mp0] [0p0] [pp0]
      !                                                                   [mmp] [0mp] [pmp]    [m0p] [00p] [p0p]    [mpp] [0pp] [ppp]

      f0 = sum(kernel_0*f)
      dfdx = sum(kernel_x*f)
      dfdy = sum(kernel_y*f)
      dfdz = sum(kernel_z*f)
      d2fdx2 = sum(kernel_xx*f)
      d2fdy2 = sum(kernel_yy*f)
      d2fdz2 = sum(kernel_zz*f)
      d2fdxdy = sum(kernel_xy*f)
      d2fdydz = sum(kernel_yz*f)
      d2fdzdx = sum(kernel_zx*f)

      return
   end subroutine twentySevenPointStencil

   subroutine extremum(f, f0, x, y, z)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find the extremum point in this 3x3x3 block of data

      real(kind=real64), dimension(-1:1, -1:1, -1:1), intent(in)      ::      f
      real(kind=real64), intent(out)                               ::      f0, x, y, z
 real(kind=real64)                                           ::      dfdx,dfdy,dfdz, d2fdx2,d2fdy2,d2fdz2, d2fdxdy, d2fdydz, d2fdzdx

      real(kind=real64)           ::      dd

      call twentySevenPointStencil(f, f0, dfdx, dfdy, dfdz, d2fdx2, d2fdy2, d2fdz2, d2fdxdy, d2fdydz, d2fdzdx)

            dd = (  2*d2fdxdy*d2fdydz*d2fdzdx + d2fdx2*d2fdy2*d2fdz2 - d2fdx2*d2fdydz*d2fdydz - d2fdy2*d2fdzdx*d2fdzdx - d2fdz2*d2fdxdy*d2fdxdy )        !   determinant of 3x3 matrix
      if (abs(dd) < 1.0d-16) then
         x = 0; y = 0; z = 0
         return
      end if

   x = (-d2fdydz*d2fdydz + d2fdy2*d2fdz2)*dfdx + (-d2fdxdy*d2fdz2 + d2fdzdx*d2fdydz)*dfdy + (-d2fdzdx*d2fdy2 + d2fdxdy*d2fdydz)*dfdz
   y = (-d2fdxdy*d2fdz2 + d2fdzdx*d2fdydz)*dfdx + (-d2fdzdx*d2fdzdx + d2fdx2*d2fdz2)*dfdy + (-d2fdydz*d2fdx2 + d2fdxdy*d2fdzdx)*dfdz
   z = (-d2fdzdx*d2fdy2 + d2fdxdy*d2fdydz)*dfdx + (-d2fdydz*d2fdx2 + d2fdxdy*d2fdzdx)*dfdy + (-d2fdxdy*d2fdxdy + d2fdx2*d2fdy2)*dfdz

      !print *,dd
      dd = 1/dd
      x = x*dd
      y = y*dd
      z = z*dd

      f0 = f0 + dfdx*x + dfdy*y + dfdz*z &
           + (d2fdx2*x*x + d2fdy2*y*y + d2fdz2*z*z + 2*d2fdxdy*x*y + 2*d2fdydz*y*z + 2*d2fdzdx*z*x)/2

      return
   end subroutine extremum

   pure function twentyFivePointLaplacian(f, ix, iy) result(L)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given the 5x5 block of values f
      !*      find the laplacian
      real(kind=real64), dimension(-2:2, -2:2), intent(in)       ::      f
      integer, intent(in)                                      ::      ix, iy
      real(kind=real64)                                       ::      L

      integer         ::      jx, jy

      real(kind=real64), dimension(5, 5), parameter          ::      kernel00 = reshape((/ &
                                                                                        4, 1, 0, 1, 4, &
                                                                                        1, -2, -3, -2, 1, &
                                                                                        0, -3, -4, -3, 0, &
                                                                                        1, -2, -3, -2, 1, &
                                                                                        4, 1, 0, 1, 4/), (/5, 5/))/35.0d0
      real(kind=real64), dimension(5, 5), parameter          ::      kernel10 = reshape((/ &
                                                                                        -3, 14, 0, -10, 19, &
                                                                                        -3, 11, -6, -19, 7, &
                                                                                        -3, 10, -8, -22, 3, &
                                                                                        -3, 11, -6, -19, 7, &
                                                                                        -3, 14, 0, -10, 19/), (/5, 5/))/70.0d0
      real(kind=real64), dimension(5, 5), parameter          ::      kernel11 = reshape((/ &
                                                                                        -14, 9, -3, -15, 8, &
                                                                                        9, 26, 10, -4, 19, &
                                                                                        -3, 10, -8, -22, 3, &
                                                                                        -15, -4, -22, -34, -5, &
                                                                                        8, 19, 3, -5, 30/), (/5, 5/))/70.0d0
      real(kind=real64), dimension(5, 5), parameter          ::      kernel20 = reshape((/ &
                                                                                        -7, 13, 0, -11, 15, &
                                                                                        -4, 13, -3, -17, 6, &
                                                                                        -3, 13, -4, -19, 3, &
                                                                                        -4, 13, -3, -17, 6, &
                                                                                        -7, 13, 0, -11, 15/), (/5, 5/))/35.0d0
      real(kind=real64), dimension(5, 5), parameter          ::      kernel21 = reshape((/ &
                                                                                        -25, 21, -3, -27, 19, &
                                                                                        4, 41, 10, -19, 24, &
                                                                                        -6, 26, -8, -38, 6, &
                                                                                        -20, 11, -22, -49, 0, &
                                                                                        -3, 31, 3, -17, 41/), (/5, 5/))/70.0d0
      real(kind=real64), dimension(5, 5), parameter          ::      kernel22 = reshape((/ &
                                                                                        -18, 8, -3, -16, 4, &
                                                                                        8, 28, 13, -2, 18, &
                                                                                        -3, 13, -4, -19, 3, &
                                                                                        -16, -2, -19, -32, -6, &
                                                                                        4, 18, 3, -6, 26/), (/5, 5/))/35.0d0
      real(kind=real64), dimension(5, 5)        ::      kernel

      select case (ix)
      case (-2)
         select case (iy)
         case (-2)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel22(6 - jx, 6 - jy); end do; end do
         case (-1)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel21(6 - jx, 6 - jy); end do; end do
         case (0)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel20(6 - jx, jy); end do; end do
         case (1)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel21(6 - jx, jy); end do; end do
         case (2)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel22(6 - jx, jy); end do; end do
         end select
      case (-1)
         select case (iy)
         case (-2)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel21(6 - jy, 6 - jx); end do; end do
         case (-1)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel11(6 - jx, 6 - jy); end do; end do
         case (0)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel10(6 - jy, jx); end do; end do
         case (1)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel11(6 - jx, jy); end do; end do
         case (2)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel21(6 - jy, jx); end do; end do
         end select
      case (0)
         select case (iy)
         case (-2)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel20(6 - jy, jx); end do; end do
         case (-1)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel10(6 - jy, jx); end do; end do
         case (0)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel00(jx, jy); end do; end do
         case (1)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel10(jy, jx); end do; end do
         case (2)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel20(jy, jx); end do; end do
         end select
      case (1)
         select case (iy)
         case (-2)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel21(6 - jy, jx); end do; end do
         case (-1)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel11(6 - jx, jy); end do; end do
         case (0)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel10(jx, jy); end do; end do
         case (1)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel11(jx, jy); end do; end do
         case (2)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel21(jy, jx); end do; end do
         end select
      case (2)
         select case (iy)
         case (-2)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel22(jx, 6 - jy); end do; end do
         case (-1)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel21(jx, 6 - jy); end do; end do
         case (0)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel20(jx, jy); end do; end do
         case (1)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel21(jx, jy); end do; end do
         case (2)
            do jy = 1, 5; do jx = 1, 5; kernel(jx, jy) = kernel22(jx, jx); end do; end do
         end select
      end select

      L = sum(kernel*f)
      return
   end function twentyFivePointLaplacian

end module Lib_LaplacianOfGaussians
