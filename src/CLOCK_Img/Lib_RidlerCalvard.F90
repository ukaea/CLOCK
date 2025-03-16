
module Lib_RidlerCalvard
!---^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_RidlerCalvard from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      provide functionality for the Ridler and Calvard thresholding criterion:
!*          I_thresh = ( I_back + I_fore )/2
!*      where I_back,I_fore is average background/foreground intensity
!*
!*
   use iso_fortran_env
   implicit none

   !external    ::      DSYSV,SSYSV

   !---    function calls provided

   public      ::      estimateMaxpixFromDarkPixels
   public      ::      findHist
   public      ::      findRidlerAndCalvardThresh
   public      ::      findSignalToNoise
   public      ::      findImageIntensityFeatures
   public      ::      findImageIntensityStdDev

#ifdef DEBUG
   logical, private, parameter         ::      LIB_RANDC_DBG = .true.
#else
   logical, private, parameter         ::      LIB_RANDC_DBG = .false.
#endif

   !---    define a 64bit floating constant "RIDLERCALVARD_IGNORE" with a really unlikely bit pattern
   integer(kind=int64), private, parameter               ::      BADF00D = int(z'BADF00D', kind=int64)
        real(kind=real64),public,parameter                  ::      RIDLERCALVARD_IGNORE = transfer( (BADF00D+ishft(BADF00D,32_int64)),1.0d0 )

   interface findRidlerAndCalvardThresh
      module procedure findRidlerAndCalvardThresh0
      module procedure findRidlerAndCalvardThresh1
   end interface

   interface estimateMaxpixFromDarkPixels
      module procedure estimateMaxpixFromDarkPixels0
      module procedure estimateMaxpixFromDarkPixels1
   end interface

contains
!---^^^^^^^^

   subroutine findImageIntensityFeatures(img, bb, tt, ff, bstd, fbar, snr, fot)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given an input image, compute the avg background, threshold and foreground levels
      !*      using the Ridler & Calvard criterion
      !*      Compute a signal:noise ratio from the threshold and std dev in background
      !*      Compute the fraction of pixels over the threshold level.
      real(kind=real64), dimension(:, :), intent(in)         ::      img         !Image array in.
      real(kind=real64), intent(out)                       ::      bb, tt, ff    !Average background, threshold and foreground intensities
      real(kind=real64), intent(out)                       ::      bstd, fbar, snr !background std dev, XXX?, signal to noise ratio
      real(kind=real64), intent(out)                       ::      fot !fraction of pixels over threshold

      real(kind=real64), dimension(0:255)          ::      hist !histogram of image intensities
      integer             ::      ii
      real(kind=real64)   ::      maxPix

      integer             ::      Nx, Ny !image x & y dimensions
      integer             ::      ix, iy, nb, nn !pixel indicies, number of background and not-ignored pixels respectively.
      real(kind=real64)   ::      gg, bbar, b2bar  !

      real(kind=real64)   ::      img_min, img_max, img_range, img_irange !image min, max, range, inverse range respectively.

      img_min = minval(img)
      img_max = maxval(img)
      img_range = img_max - img_min
      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      if (LIB_RANDC_DBG) then
         print *, "Lib_RidlerCalvard::findImageIntensityFeatures info - img size ", Nx, Ny, " minval ", img_min, " maxval ", img_max
      end if

      img_irange = 1/max(1.0d-8, img_range)

      call findHist((img - img_min)*img_irange, hist)
      call estimateMaxpixFromDarkPixels1(hist, maxPix)
      call findRidlerAndCalvardThresh1(hist, maxPix, bb, tt, ff)

      bbar = 0.0d0
      fbar = 0.0d0
      b2bar = 0.0d0
      nb = 0
      nn = 0
      do iy = 1, size(img, dim=2)
         do ix = 1, size(img, dim=1)
            if (img(ix, iy) == RIDLERCALVARD_IGNORE) cycle
            gg = (img(ix, iy) - img_min)*img_irange
            fbar = fbar + gg
            nn = nn + 1
            if (img(ix, iy) >= tt) cycle
            bbar = bbar + gg
            b2bar = b2bar + gg*gg
            nb = nb + 1
         end do
      end do

      if (nn < 3) then
         !   error - insufficient pixels
         fbar = 0.0d0

      else
         fbar = fbar/nn
      end if

      if (nb < 3) then
         !   error - insufficient pixels under threshold
         snr = 0.0d0
         bbar = 0.0d0
         gg = 0.001d0
         !return
      else

         bbar = bbar/nb
         gg = (b2bar - bbar*bbar*nb)/(nb - 1)          !   this is  (<b^2> - <b>^2) * nb/(nb-1)

      end if

      if (gg < 1.0d-8*tt) then
         !   variance appears to be zero, but the threshold is above background? Not a realistic scene!
         snr = 0.0d0
      else

         bstd = sqrt(gg)
         snr = (tt - bbar)/bstd

      end if

      fot = 0.0d0
      do ii = 0, 255
         if (ii > tt*255) fot = fot + hist(ii)
      end do

      bb = bbar

      !---    rescale back to original intensity range
      bb = bb*img_range + img_min
      tt = tt*img_range + img_min
      ff = ff*img_range + img_min
      bstd = bstd*img_range + img_min
      fbar = fbar*img_range + img_min

      return
   end subroutine findImageIntensityFeatures

   function findSignalToNoise(img, thresh) result(snr)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(:, :), intent(in)         ::      img
      real(kind=real64), intent(in)                        ::      thresh
      real(kind=real64)                                   ::      snr

      integer             ::      ix, iy, nb
      real(kind=real64)   ::      gg, bbar, b2bar, bstd

      bbar = 0.0d0
      b2bar = 0.0d0
      nb = 0
      do iy = 1, size(img, dim=2)
         do ix = 1, size(img, dim=1)
            gg = img(ix, iy)
            if (gg == RIDLERCALVARD_IGNORE) cycle
            if (gg >= thresh) cycle
            bbar = bbar + gg
            b2bar = b2bar + gg*gg
            nb = nb + 1
         end do
      end do

      if (LIB_RANDC_DBG) then
         print *, "Lib_RidlerCalvard::findSignalToNoise info - nb,bbar,b2bar ", nb, bbar/max(1, nb), b2bar/max(1, nb)
      end if

      if (nb < 3) then
         !   error - insufficient pixels under threshold
         snr = 0.0d0
         return
      end if

      bbar = bbar/nb
      gg = (b2bar - bbar*bbar*nb)/(nb - 1)          !   this is  (<b^2> - <b>^2) * nb/(nb-1)

      if (gg < 1.0d-8*thresh) then
         !   variance appears to be zero, but the threshold is above background? Not a realistic scene!
         snr = 0.0d0
         return
      end if

      bstd = sqrt(gg)

      snr = (thresh - bbar)/bstd

      return
   end function findSignalToNoise

   subroutine estimateMaxpixFromDarkPixels0(img, maxpix)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      use the dark pixels to estimate the noise level, and from this get a guess for the proportion of bright pixels
      real(kind=real64), dimension(:, :), intent(in)         ::      img
      real(kind=real64), intent(out)                       ::      maxpix

      real(kind=real64), dimension(0:255)                  ::      hist

      !---    generate a histogram
      call findHist(img, hist)
      call estimateMaxpixFromDarkPixels1(hist, maxpix)

      return
   end subroutine estimateMaxpixFromDarkPixels0

   subroutine estimateMaxpixFromDarkPixels1(hist, maxpix)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      use the dark pixels to estimate the noise level, and from this get a guess for the proportion of bright pixels
      real(kind=real64), dimension(0:255), intent(in)       ::      hist
      real(kind=real64), intent(out)                       ::      maxpix

      real(kind=real64)                   ::      xx, yy, hmin, hh, ss, y0, y1, xold
      integer                             ::      ii, mm, loop
      integer, dimension(3)                ::      ipiv
      real(kind=real64), dimension(3)      ::      bb
      real(kind=real64), dimension(3, 3)    ::      aa

      integer, parameter                   ::      MAXLOOPS = 100

      real(kind=real64), dimension(3*66)   ::      work

      ss = 0; xx = 0; xold = huge(1.0)
      do loop = 1, MAXLOOPS

         if (loop == 1) then
            !---    find maximum occupation bin ( excluding 0 and 1 )
            mm = maxloc(hist(1:254), dim=1) + 16
            mm = min(mm, 254)
         else
            mm = int((xx + 2*ss))      !   add 2 sigma to make sure...
         end if

         mm = max(1, min(255, mm))

         !---    fit to gaussian below this level
         !       hist = h0 exp( - (x-x0)^2/(2s^2) )
         !       write ln( hist ) = b0 x^2 + b1 x + b2
         !       solve with linear least squares
         !       then b0 = - x^2/(2s^2)
         !           b1 =  (x0/s^2)
         !           b2 = ln(h0) -x0^2/2s^2

   if (LIB_RANDC_DBG) print *, "Lib_RidlerCalvard::estimateMaxpixFromDarkPixels dbg -  fitting pixel intensities ", loop, xx, ss, mm

         hmin = 1.0d0/(4096*4096)      !   less than one pixel at this level is the minimum possible
         aa = 0.0d0
         bb = 0.0d0
         do ii = 1, mm

            xx = ii
            if (hist(ii) < hmin) cycle
            yy = log(hist(ii))

            aa(1, 1) = aa(1, 1) + xx*xx*xx*xx
            aa(2, 1) = aa(2, 1) + xx*xx*xx
            aa(3, 1) = aa(3, 1) + xx*xx
            aa(2, 2) = aa(2, 2) + xx*xx
            aa(3, 2) = aa(3, 2) + xx
            aa(3, 3) = aa(3, 3) + 1

            bb(1) = bb(1) + yy*xx*xx
            bb(2) = bb(2) + yy*xx
            bb(3) = bb(3) + yy

         end do

         call DSYSV("L", 3, 1, aa, 3, ipiv, bb, 3, work, size(work), ii)

         !---    compute estimated standard deviation
         if (bb(1) >= 0) then
            maxPix = 0.5d0
            if (LIB_RANDC_DBG) print *, "Lib_RidlerCalvard::estimateMaxpixFromDarkPixels dbg - maxpix ", maxpix
            return
         end if
         ss = sqrt(-1/(2*bb(1)))
         xx = bb(2)*ss*ss
         hh = exp(bb(3) + xx*xx/(2*ss*ss))
      if (LIB_RANDC_DBG) print *, "Lib_RidlerCalvard::estimateMaxpixFromDarkPixels dbg - loop ", loop, " mu,sigma,I0 = ", xx, ss, hh
         if (abs(xx - xold) < 1.0d0/(256*256)) exit
         xold = xx
      end do

      !---    compute total intensity found in and expected in histogram
      y0 = 0.0d0
      y1 = 0.0d0
      do ii = 1, 254
         yy = hh*exp(-(ii - xx)*(ii - xx)/(2*ss*ss))
         if (hist(ii) > yy) then
            y0 = y0 + hist(ii)
            y1 = y1 + yy
         end if
      end do

      if (LIB_RANDC_DBG) print *, "Lib_RidlerCalvard::estimateMaxpixFromDarkPixels info - y0,y1 = ", y0, y1

      !---    maxpixf

      maxpix = 2*(y0 - y1)/y0

      maxpix = max(0.001d0, min(0.999d0, maxpix))
      if (LIB_RANDC_DBG) print *, "Lib_RidlerCalvard::estimateMaxpixFromDarkPixels dbg - maxpix ", maxpix

      return
   end subroutine estimateMaxpixFromDarkPixels1

   subroutine findHist(img_in, hist)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find simple histogram of intensities
      real(kind=real64), dimension(:, :), intent(in)         ::      img_in
      real(kind=real64), dimension(0:255), intent(out)      ::      hist

      integer             ::      ix, iy
      integer             ::      hh, nn
      real(kind=real64)   ::      gg, hsum, gbar, g2bar

      !---    find a histogram of the intensity
      hist(0:255) = 0.0d0
      gbar = 0.d0; g2bar = 0.0d0; nn = 0
      do iy = 1, size(img_in, dim=2)
         do ix = 1, size(img_in, dim=1)
            !---    add to global histogram
            gg = img_in(ix, iy)
            if (gg == RIDLERCALVARD_IGNORE) cycle
            hh = max(0, min(255, int(gg*256)))
            hist(hh) = hist(hh) + 1.0d0
            gbar = gbar + gg
            g2bar = g2bar + gg*gg
            nn = nn + 1
         end do
      end do

      if (nn > 2) then
         gbar = gbar/nn
         g2bar = g2bar/nn
      end if

      hsum = sum(hist(0:255))
      if (hsum > 0) hist(0:255) = hist(0:255)/hsum

      if (LIB_RANDC_DBG) then
         do hh = 0, 255
            print *, "Lib_RidlerCalvard::findHist info - ", hh, hist(hh)
         end do
      end if

      return
   end subroutine findHist

   !---

   subroutine findRidlerAndCalvardThresh0(img_in, maxPix, bb, tt, ff)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(:, :), intent(in)         ::      img_in
      real(kind=real64), intent(in)                        ::      maxPix
      real(kind=real64), intent(out)                       ::      bb, tt, ff

      real(kind=real64), dimension(0:255)          ::      hist
      call findHist(img_in, hist)
      call findRidlerAndCalvardThresh1(hist, maxPix, bb, tt, ff)

      return
   end subroutine findRidlerAndCalvardThresh0

   subroutine findRidlerAndCalvardThresh1(hist, maxPix, bb, tt, ff)
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
      if (btf(1) == 0) then
         bb = 0.0d0
         tt = 0.5d0
         ff = 1.0d0
      else
         bb = btf(1)
         tt = btf(2)
         ff = btf(3)
      end if

      return
   end subroutine findRidlerAndCalvardThresh1

   subroutine findImageIntensityStdDev(img, stdDevOut)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      Finds the standard deviation of the pixel intensities of an image.
      real(kind=real64), dimension(:, :), intent(in)         ::      img         !Image array in.
      real(kind=real64), intent(out)                       ::      stdDevOut    !standard deviation of all pixel intensities

      integer             ::      ii
      real(kind=real64)   ::      maxPix

      integer             ::      Nx, Ny !image x & y dimensions
      integer             ::      ix, iy, nb, nn !pixel indicies, number of background and not-ignored pixels respectively.
      real(kind=real64)   ::      gg, bbar, b2bar  !

      real(kind=real64)   ::      img_min, img_max, img_range, img_irange !image min, max, range, inverse range respectively.

      img_min = minval(img)
      img_max = maxval(img)
      img_range = img_max - img_min
      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      img_irange = 1/max(1.0d-8, img_range)

      bbar = 0.0d0
      b2bar = 0.0d0
      nb = 0
      nn = 0
      do iy = 1, size(img, dim=2)
         do ix = 1, size(img, dim=1)
            gg = (img(ix, iy) - img_min)*img_irange
            nn = nn + 1
            bbar = bbar + gg
            b2bar = b2bar + gg*gg
            nb = nb + 1
         end do
      end do

      if (nb < 3) then
         !   error - insufficient pixels under threshold
         bbar = 0.0d0
         gg = 0.001d0
         !return
      else

         bbar = bbar/nb
         gg = (b2bar - bbar*bbar*nb)/(nb - 1)          !   this is  (<b^2> - <b>^2) * nb/(nb-1)

      end if

      stdDevOut = sqrt(gg)

      !---    rescale back to original intensity range

      stdDevOut = stdDevOut*img_range + img_min

      return
   end subroutine findImageIntensityStdDev

end module Lib_RidlerCalvard
