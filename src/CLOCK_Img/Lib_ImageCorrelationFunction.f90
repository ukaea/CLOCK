
module Lib_ImageCorrelationFunction
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_ImageCorrelationFunction from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      compute the correlation function
!*      c =                 sum (f-<f>)(g-<g>)
!*          -------------------------------------------
!*          sqrt( sum (f-<f>)^2 ) sqrt( sum (g-<g>)^2 )

   use iso_fortran_env

   use Lib_png
   implicit none
   private

   public      ::      imageCorrelationFunction
   public      ::      shrinkImage

   public      ::      optimalImageCorrelationFunction

   integer(kind=int64), private, parameter      ::      BADF00D = int(z'BADF00D', kind=int64)
   real(kind=real64), public, parameter         ::      ICF_DEAD_PIXEL = transfer((BADF00D + ishft(BADF00D, 32_int64)), 1.0d0)

contains
!---^^^^^^^^

   pure function imageCorrelationFunction(f, g, mask) result(c)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(:, :), intent(in)         ::      f, g
      logical, dimension(:, :), intent(in), optional          ::      mask
      real(kind=real64)                                   ::      c

      real(kind=real64)           ::      sumf, sumff, sumg, sumgg, sumfg
      real(kind=real64)           ::      barf, barff, barg, bargg, barfg
      integer                     ::      Nx, Ny
      integer                     ::      ix, iy, nPix
      real(kind=real64)           ::      ff, gg
      real(kind=real64)           ::      dd
      sumf = 0.0d0
      sumff = 0.0d0
      sumg = 0.0d0
      sumgg = 0.0d0
      sumfg = 0.0d0

      Nx = size(f, dim=1)
      Ny = size(f, dim=2)
      nPix = 0

      do iy = 1, Ny
         do ix = 1, Nx
            if (present(mask)) then
               if (.not. mask(ix, iy)) cycle
            end if
            nPix = nPix + 1

            ff = f(ix, iy)
            gg = g(ix, iy)

            sumf = sumf + ff
            sumff = sumff + ff*ff
            sumg = sumg + gg
            sumgg = sumgg + gg*gg
            sumfg = sumfg + ff*gg
         end do
      end do

      c = 0.0d0
      if (nPix == 0) return

      dd = 1.0d0/nPix
      barf = sumf*dd
      barff = sumff*dd
      barg = sumg*dd
      bargg = sumgg*dd
      barfg = sumfg*dd

      dd = (barff - barf*barf)*(bargg - barg*barg)
      if (dd <= 0) return

      dd = 1/sqrt(dd)
      c = (barfg - barf*barg)*dd

      return
   end function imageCorrelationFunction

   subroutine optimalImageCorrelationFunction(f, g, c, alpha, delta)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find the image correlation between f and g
      !*      and optimise wrt strain a and displacement d of image f
      real(kind=real64), dimension(:, :), intent(in)             ::      f
      real(kind=real64), dimension(:, :), intent(in)             ::      g
      real(kind=real64), intent(out)                           ::      c
      real(kind=real64), intent(out), optional                           ::      alpha
      real(kind=real64), dimension(2), intent(out), optional              ::      delta

      real(kind=real64), dimension(0:size(f, dim=1) + 1, 0:size(f, dim=2) + 1)    ::      f_tmp, g_tmp

      integer             ::      Nx, Ny       !   size of original image
      integer             ::      Mx, My       !   current scaled size of image
      integer             ::      ix, iy
      real(kind=real64), dimension(3)          ::      dc
      real(kind=real64), dimension(3, 3)        ::      d2c
      real(kind=real64), dimension(2)          ::      delta0

      real(kind=real64)   ::      dd, idd, aa, dx, dy, sx, sy, da, alpha0
      real(kind=real64)   ::      ff, gg, c0, cdenom
      real(kind=real64)   ::      barf, barg, varf, varg

      logical, dimension(size(f, dim=1), size(f, dim=2))      ::      mask
      integer         ::      NDIVa = 2, NDIVd = 2
      real(kind=real64), dimension(:, :, :), allocatable         ::      cmat
      real(kind=real64), dimension(:, :, :, :), allocatable       ::      dcmat
      real(kind=real64), dimension(:, :, :, :, :), allocatable     ::      d2cmat
      integer, dimension(3)                                ::      ii
      integer         ::      ka, kx, ky

      Nx = size(f, dim=1)
      Ny = size(f, dim=2)

      !---    compute the mean and variance of the input images
      barf = 0.0d0
      barg = 0.0d0
      varf = 0.0d0
      varg = 0.0d0
      do iy = 1, Ny
         do ix = 1, Nx
            ff = f(ix, iy)
            gg = g(ix, iy)
            barf = barf + ff        !   compute sum
            varf = varf + ff*ff     !   compute sum of squares.
            barg = barg + gg
            varg = varg + gg*gg
         end do
      end do
      !   normalise sums to averages
      idd = 1.0d0/(Nx*Ny)
      barf = barf*idd
      barg = barg*idd
      varf = varf*idd
      varg = varg*idd
      !   convert mean square to variance
      varf = varf - barf*barf
      varg = varg - barg*barg

      !   image correlation function has a denominator sqrt( varf * varg )
            print *,"Lib_ImageCorrelationFunction::optimalImageCorrelationFunction info - barf,barg,stdev f,stdev g ",barf,barg,sqrt(varf),sqrt(varg)
      cdenom = varf*varg
      mask = .true.
      if (cdenom <= 0) then
         !   no variation in the images
         if (present(alpha)) alpha = 1.0d0
         if (present(delta)) delta = 0.0d0
         c = 1.0d0
         return
      end if
      cdenom = 1/sqrt(cdenom)

      if (present(alpha)) then
         if (present(delta)) then
            NDIVa = 2
            NDIVd = 2
         else
            NDIVa = 8
            NDIVd = 0
         end if
      else if (present(delta)) then
         NDIVa = 0
         NDIVd = 6
      else
         c = imageCorrelationFunction(f(1:Nx, 1:Ny), g(1:Nx, 1:Ny))
         return
      end if

      !---    start with a tiny image, and successively quadruple in size until the correlation
      Mx = 4; My = 4
      alpha0 = 1.0d0       !   scale factor ( note x -> x + alpha0(x-x0) + delta0x
      delta0 = 0.0d0       !   offset factor in original pixels

      allocate (cmat(-NDIVa:NDIVa, -NDIVd:NDIVd, -NDIVd:NDIVd))
      allocate (dcmat(3, -NDIVa:NDIVa, -NDIVd:NDIVd, -NDIVd:NDIVd))
      allocate (d2cmat(3, 3, -NDIVa:NDIVa, -NDIVd:NDIVd, -NDIVd:NDIVd))
      da = 0.0d0
      dx = 0.0d0
      dy = 0.0d0
      c0 = imageCorrelationFunction(f(1:Nx, 1:Ny), g(1:Nx, 1:Ny))
      do

         call shrinkImage(g, Mx, My, 1/sqrt(alpha0), -delta0/2, g_tmp)

         if (NDIVa > 0) da = sqrt(1.0d0/(NDIVa*NDIVa*Mx*My))
         if (NDIVd > 0) dx = 1.0d0*Nx/(NDIVd*Mx)
         if (NDIVd > 0) dy = 1.0d0*Ny/(NDIVd*My)
         do ka = -NDIVa, NDIVa
            do kx = -NDIVd, NDIVd
               do ky = -NDIVd, NDIVd
                  call shrinkImage(f, Mx, My, sqrt(alpha0)*(1 + ka*da), delta0/2 + (/kx*dx, ky*dy/), f_tmp)
                  cmat(ka, kx, ky) = imageCorrelationFunction(f_tmp(1:Mx, 1:My), g_tmp(1:Mx, 1:My))
               end do
            end do
         end do
         print *, "scale ", Mx, My, " transform ", alpha0, delta0, " image correlation ", cmat(0, 0, 0)

         if (c0 > 0) then
            ii = maxloc(cmat)
         else
            ii = minloc(cmat)
         end if
         ii(1) = ii(1) - NDIVa - 1
         ii(2:3) = ii(2:3) - NDIVd - 1
         ka = ii(1); kx = ii(2); ky = ii(3)

         ! print *,"max ",ka,kx,ky," = ",ka*da,kx*dx,ky*dy,cmat(ka,kx,ky)
         call firstDerivs(cmat, dcmat)
         call firstDerivs(dcmat(1, :, :, :), d2cmat(:, 1, :, :, :))
         call firstDerivs(dcmat(2, :, :, :), d2cmat(:, 2, :, :, :))
         call firstDerivs(dcmat(3, :, :, :), d2cmat(:, 3, :, :, :))

         dc = dcmat(:, ka, kx, ky)
         d2c = d2cmat(:, :, ka, kx, ky)

         !---    check: can we find the necessary offset d2c b = dc
         dd = determinant3Mat(d2c)
         if (abs(dd) > 1.0d-16) then
            !   yes, can do matrix inverse no problem
            !   find matrix of cofactors
            !   find solution b = cof df / dd
            dc = dc/dd
            aa = (d2c(2, 2)*d2c(3, 3) - d2c(2, 3)*d2c(3, 2))*dc(1) &
                 + (d2c(1, 3)*d2c(3, 2) - d2c(1, 2)*d2c(3, 3))*dc(2) &
                 + (d2c(1, 2)*d2c(2, 3) - d2c(1, 3)*d2c(2, 2))*dc(3)
            sx = (d2c(2, 3)*d2c(3, 1) - d2c(2, 1)*d2c(3, 3))*dc(1) &
                 + (d2c(1, 1)*d2c(3, 3) - d2c(1, 3)*d2c(3, 1))*dc(2) &
                 + (d2c(1, 3)*d2c(2, 1) - d2c(1, 1)*d2c(2, 3))*dc(3)
            sy = (d2c(2, 1)*d2c(3, 2) - d2c(2, 2)*d2c(3, 1))*dc(1) &
                 + (d2c(1, 2)*d2c(3, 1) - d2c(1, 1)*d2c(3, 2))*dc(2) &
                 + (d2c(1, 1)*d2c(2, 2) - d2c(1, 2)*d2c(2, 1))*dc(3)

            if ((abs(aa) < 1) .and. (abs(sx) < 1) .and. (abs(sy) < 1)) then
               aa = (ka - aa)*da
               sx = (kx - sx)*dx
               sy = (ky - sy)*dy
            else
               aa = ka*da
               sx = kx*dx
               sy = ky*dy
            end if
         else
            !   matrix iss singluar??
            aa = ka*da
            sx = kx*dx
            sy = ky*dy
         end if

         alpha0 = alpha0*(1 + aa)
         delta0(1) = delta0(1) + sx
         delta0(2) = delta0(2) + sy
         call shrinkImage(f, Mx, My, sqrt(alpha0), delta0/2, f_tmp)
         call shrinkImage(g, Mx, My, 1/sqrt(alpha0), -delta0/2, g_tmp)
         c0 = imageCorrelationFunction(f_tmp(1:Mx, 1:My), g_tmp(1:Mx, 1:My), mask)

         if ((Mx == Nx) .and. (My == Ny)) exit
         Mx = min(Nx, Mx*2); My = min(Ny, My*2)

         !  stop
      end do

      call shrinkImage(f, Nx, Ny, sqrt(alpha0), delta0/2, f_tmp)
      call shrinkImage(g, Nx, Ny, 1/sqrt(alpha0), -delta0/2, g_tmp)
      c = imageCorrelationFunction(f_tmp(1:Nx, 1:Ny), g_tmp(1:Nx, 1:Ny))

      if (present(alpha)) alpha = alpha0
      if (present(delta)) delta = delta0

      return
   end subroutine optimalImageCorrelationFunction

   pure function determinant3Mat(M) result(d)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns the determinant of M
      real(kind=real64), dimension(3, 3), intent(in)      ::      M
      real(kind=real64)                                ::      d
      real(kind=real64), dimension(9)       ::      dd
      dd(1) = M(1, 1)*(M(2, 2)*M(3, 3) - M(2, 3)*M(3, 2))
      dd(2) = M(1, 2)*(M(2, 3)*M(3, 1) - M(2, 1)*M(3, 3))
      dd(3) = M(1, 3)*(M(2, 1)*M(3, 2) - M(2, 2)*M(3, 1))
      dd(4) = M(2, 1)*(M(3, 2)*M(1, 3) - M(3, 3)*M(1, 2))
      dd(5) = M(2, 2)*(M(3, 3)*M(1, 1) - M(3, 1)*M(1, 3))
      dd(6) = M(2, 3)*(M(3, 1)*M(1, 2) - M(3, 2)*M(1, 1))
      dd(7) = M(3, 1)*(M(1, 2)*M(2, 3) - M(1, 3)*M(2, 2))
      dd(8) = M(3, 2)*(M(1, 3)*M(2, 1) - M(1, 1)*M(2, 3))
      dd(9) = M(3, 3)*(M(1, 1)*M(2, 2) - M(1, 2)*M(2, 1))
      d = (1.0d0/3.0d0)*sum(dd)
      return
   end function determinant3Mat

   subroutine findHist(img_in, hist)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find simple histogram of intensities
      real(kind=real64), dimension(:, :), intent(in)         ::      img_in
      real(kind=real64), dimension(0:255), intent(out)      ::      hist

      integer         ::      ix, iy
      integer         ::      hh, nn
      real(kind=real64)   ::      gg, hsum, gbar, g2bar

      !---    find a histogram of the intensity
      hist(0:255) = 0.0d0
      gbar = 0.d0; g2bar = 0.0d0; nn = 0
      do iy = 1, size(img_in, dim=2)
         do ix = 1, size(img_in, dim=1)

            !---    add to global histogram
            gg = img_in(ix, iy)
            hh = min(255, int(gg*256))
            hist(hh) = hist(hh) + 1.0d0
            if (gg > 0) then
               gbar = gbar + gg
               g2bar = g2bar + gg*gg
               nn = nn + 1
            end if
         end do
      end do

      gbar = gbar/max(1, nn)
      g2bar = g2bar/max(1, nn)
      print *, "findHist info - average intensity ", gbar, " in ", nn, " pixels, noise ", sqrt((g2bar - gbar*gbar)*nn/(nn - 1.5d0))

      hsum = sum(hist(1:255))
      if (hsum > 0) hist(0:255) = hist(0:255)/hsum

      return
   end subroutine findHist

   !---

   subroutine findRidlerAndCalvardThresh(hist, maxPix, bb, tt, ff)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(0:255), intent(in)           ::      hist
      real(kind=real64), intent(in)                            ::      maxPix
      real(kind=real64), intent(out)                           ::      bb, tt, ff

      integer             ::      ii, hh
      real(kind=real64)   ::      nf, nb
      real(kind=real64), dimension(3)      ::      btf

      btf(1:3) = (/0.0d0, 1.0d0, 1.0d0/)
      do hh = 254, 2, -1

         !---    find average background level below hh
         bb = 0.0d0; nb = 0.0d0
         do ii = 1, hh - 1
            bb = bb + ii*hist(ii)
            nb = nb + hist(ii)
         end do
         if (nb > 0) bb = bb*0.003921568627451d0/nb

         !---    find average foreground level above hh
         ff = 0.0d0; nf = 0.0d0
         do ii = hh + 1, 255
            ff = ff + ii*hist(ii)
            nf = nf + hist(ii)
         end do
         if (nf > 0) ff = ff*0.003921568627451d0/nf

         !---    try to place the threshold at hh
         tt = hh*0.003921568627451d0       !       hh/255

         !---    does this fit R&C criterion
         if ((nf <= maxPix*(nf + nb)) .and. (2*tt > (ff + bb))) then
            btf(1:3) = (/bb, tt, ff/)
         end if

      end do

      !---    store background/threshold/foreground
      bb = btf(1)
      tt = btf(2)
      ff = btf(3)
      print *, "findRidlerAndCalvardThresh info - back,thresh,fore ", bb, tt, ff
      return
   end subroutine findRidlerAndCalvardThresh

   subroutine shrinkImage(img_in, Mx, My, alpha, delta, img_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      scale img_in by alpha
      !*      translate by delta
      !*      then reduce to Mx,My pixels
      !*      note: also provides a 1 px border
      real(kind=real64), dimension(:, :), intent(in)         ::      img_in
      integer, intent(in)                                  ::      Mx, My
      real(kind=real64), intent(in)                        ::      alpha
      real(kind=real64), dimension(2), intent(in)           ::      delta
      real(kind=real64), dimension(0:, 0:), intent(out)      ::      img_out

      integer             ::      Nx, Ny       !   original image size
      integer             ::      ix, iy       !   pixel in original image
      real(kind=real64)   ::      xx, yy       !   position of pixel after scale/translate in original image
      integer             ::      jx, jy       !   pixel in output image
      real(kind=real64)   ::      x0, y0       !   centre of original image
      real(kind=real64)   ::      sx, sy       !   scaling of x- y- axes
      real(kind=real64), dimension(-1:Mx + 2, -1:My + 2)      ::      img_tmp       !   note 2 px border
      integer, dimension(-1:Mx + 2, -1:My + 2)                ::      n_hits        !   number of contributors to pixel in output image
      integer             ::      kx, ky
      logical             ::      ok

      n_hits = 0
      Nx = size(img_in, dim=1)
      Ny = size(img_in, dim=2)
      x0 = (Nx - 1)*0.5d0       !   -1 because I take the centre of pixel(1,1) to be at real space location (0.5,0.5)
      y0 = (Ny - 1)*0.5d0
      sx = real(Mx, kind=real64)/Nx
      sy = real(My, kind=real64)/Ny
      img_tmp = 0

      do iy = 1, Ny
         yy = y0 + alpha*(iy - 0.5d0 - y0) + delta(2)
         jy = floor(yy*sy) + 1
         if ((jy + 1)*(My + 2 - jy) < 0) cycle      !   outside range -1:My+2

         do ix = 1, Nx
            !   find position of pixel after scale/translate
            xx = x0 + alpha*(ix - 0.5d0 - x0) + delta(1)
            !   find position of pixel in output image
            jx = floor(xx*sx) + 1
            if ((jx + 1)*(Mx + 2 - jx) < 0) cycle

            !   add to pixel in output image
            n_hits(jx, jy) = n_hits(jx, jy) + 1
            img_tmp(jx, jy) = img_tmp(jx, jy) + img_in(ix, iy)
         end do
      end do

      !---    now average each output pixel that can be averaged
      ok = .true.
      do jy = -1, My + 2
         do jx = -1, Mx + 2
            if (n_hits(jx, jy) > 0) then
               img_tmp(jx, jy) = img_tmp(jx, jy)/n_hits(jx, jy)
            else
               img_tmp(jx, jy) = ICF_DEAD_PIXEL
               ok = .false.
            end if

         end do
      end do
      img_out(0:Mx + 1, 0:My + 1) = img_tmp(0:Mx + 1, 0:My + 1)

      if (ok) return

      do iy = 1, 20     !   try to fill in up to 20 dead pixels in a row
         ok = .true.
         do jy = 0, My + 1
            do jx = 0, Mx + 1
               if (dead(img_out(jx, jy))) then
                  ix = 0; sx = 0.0d0; ok = .false.
                  do ky = jy - 1, jy + 1
                     do kx = jx - 1, jx + 1
                        if (live(img_tmp(kx, ky))) then
                           ix = ix + 1
                           sx = sx + img_tmp(kx, ky)
                        end if
                     end do
                  end do
                  if (ix > 0) img_out(jx, jy) = sx/ix

               end if
            end do
         end do
         if (ok) exit

         img_tmp(0:Mx + 1, 0:My + 1) = img_out(0:Mx + 1, 0:My + 1)

      end do

      return
   end subroutine shrinkImage

   elemental function dead(f) result(is)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), intent(in)            ::      f
      logical                                 ::      is
      is = (f == ICF_DEAD_PIXEL)
      return
   end function dead

   elemental function live(f) result(is)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), intent(in)            ::      f
      logical                                 ::      is
      is = (f /= ICF_DEAD_PIXEL)
      return
   end function live

   pure subroutine firstDerivs(f, df)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(:, :, :), intent(in)       ::      f
      real(kind=real64), dimension(:, :, :, :), intent(out)    ::      df

      integer         ::      Nx, Ny, Nz
      integer         ::      ix, iy, iz

      Nx = size(f, dim=1)
      Ny = size(f, dim=2)
      Nz = size(f, dim=3)
      df = 0.0d0
      do iz = 1, Nz
         do iy = 1, Ny
            do ix = 1, Nx
               if (Nx > 3) then
                  if (ix == 1) then
                     df(1, ix, iy, iz) = (-3*f(1, iy, iz) + 4*f(2, iy, iz) - f(3, iy, iz))/2
                  else if (ix == Nx) then
                     df(1, ix, iy, iz) = (3*f(Nx, iy, iz) - 4*f(Nx - 1, iy, iz) + f(Nx - 2, iy, iz))/2
                  else
                     df(1, ix, iy, iz) = (f(ix + 1, iy, iz) - f(ix - 1, iy, iz))/2
                  end if
               end if

               if (Ny > 3) then
                  if (iy == 1) then
                     df(2, ix, iy, iz) = (-3*f(ix, 1, iz) + 4*f(ix, 2, iz) - f(ix, 3, iz))/2
                  else if (iy == Ny) then
                     df(2, ix, iy, iz) = (3*f(ix, Ny, iz) - 4*f(ix, Ny - 1, iz) + f(ix, Ny - 2, iz))/2
                  else
                     df(2, ix, iy, iz) = (f(ix, iy + 1, iz) - f(ix, iy - 1, iz))/2
                  end if
               end if

               if (Nz > 3) then
                  if (iz == 1) then
                     df(3, ix, iy, iz) = (-3*f(ix, iy, 1) + 4*f(ix, iy, 2) - f(ix, iy, 3))/2
                  else if (iz == Nz) then
                     df(3, ix, iy, iz) = (3*f(ix, iy, Nz) - 4*f(ix, iy, Nz - 1) + f(ix, iy, Nz - 2))/2
                  else
                     df(3, ix, iy, iz) = (f(ix, iy, iz + 1) - f(ix, iy, iz - 1))/2
                  end if
               end if
            end do
         end do
      end do

      return
   end subroutine firstDerivs

end module Lib_ImageCorrelationFunction
