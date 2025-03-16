
module analyseSpots
!---^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    analyseSpots from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      This module fits 2d gaussians to spots in an input image
!*

   use Lib_Gaussian2d
   use Lib_Maxima2d
   use Lib_RidlerCalvard
   use Lib_ColourScale
   use Lib_Png
   use Lib_Quicksort
   use Lib_MultipleGaussians2d
   use Lib_SimpleProgressBar

   use iso_fortran_env
   implicit none
   private

   !---    static fields
   real(kind=real64), public                ::      F_THRESH = 0.0d0    !intensity threshold of spot detection
   !real(kind=real64),public                ::      T_THRESH = 0.0d0    !unused currently
   real(kind=real64), public                ::      D_THRESH = 0.8d0    !Peak diameter merging threshold
   !real(kind=real64),public                ::      Q_THRESH = 1.0d0    !unused currently
   integer, public                ::      ROIXMAX = 200       !maximum permitted egion of interest (roi) dimension
   real(kind=real64), public                ::      ROIRMAX = 2.0d0     !maximum permitted padding region
   real(kind=real64), public                ::      MAXPIX = 0.0d0      !maximum fraction of bright pixels
   real(kind=real64), public                ::      TOND_THRESH = 0.0d0 !unused currently
   integer, public                          ::      DBG_ROIX = -1       !print roi location debug (x coord)
   integer, public                          ::      DBG_ROIY = -1       !print roi location debug (x coord)
   logical, public                          ::      OPROI = .false.     !debug toggle
   logical, public                          ::      OPRANDC = .false.   !debug toggle
   logical, public                          ::      FITANDCOUNT = .true.!Toggle spot fitting

   character(len=2), public                 ::      LENGTHUNIT = "px"     !length unit (default pixels)

   !---

   public          ::      fitSpots

   !---

   interface fitSpots
      module procedure fitSpots0
   end interface

contains
!---^^^^^^^^

   subroutine fitSpots0(img, fitted_g2d, f_auto, maxpix_auto, roi_auto, ss, ofile)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(0:, 0:), intent(inout)   ::      img

      !---    new code
      type(MultipleGaussians2d), intent(out)           ::      fitted_g2d      !   object for storing parameters for multiple 2D gaussians
      logical, intent(in)                              ::      f_auto          !   should the threshold intensity be determined automatically from the image or fixed?
      logical, intent(in)                              ::      maxpix_auto     !   should the maximum fraction of bright pixels be determined automatically from the image or fixed?
      logical, intent(in)                              ::      roi_auto        !   should the padding of the ROIs be determined automatically not to percolate, or be fixed?
      real(kind=real64), intent(out)                   ::      ss              !   background standard deviation
      character(len=*), intent(in), optional            ::      ofile           !   spot output file

      integer         ::      nx, ny
      integer, dimension(:, :), allocatable              ::      indx            !indcies of detected maxima
      integer, dimension(:), allocatable                ::      left, right, top, bottom !indices of the fourcorners of the bounding box of a roi
      integer                                         ::      xmax, ymax       !   maximum x and y dimension of a given roi
      real(kind=real64), dimension(:, :), allocatable    ::      img_roi         !   Pixles contained in subject roi
      integer, dimension(:, :), allocatable              ::      indx_roi        !   index each point as belonging to a single maximum indx(i,j) = 1:nMax

      real(kind=real64)                   ::      bb, tt, ff    !Background, threshold and foreground intensity
      real(kind=real64)                   ::      dd, xx       !debug roi intensity, inverse number of rois
      logical                             ::      ok          !is a given roi ok (boolean)
      integer                             ::      roi, nroi, nSpot, nReducedMax, nRawMax !roi index, nmuber of: roi, spots, reduced maxima and raw maxima respectively.
      integer                             ::      ix, iy, ng    !debug roi x and y index. number of spots detected respectively.
      real(kind=real64)                   ::      detectFminOnS !F_THRESH, intensity threshold of spot detection
      real(kind=real64), dimension(:, :, :), allocatable    ::      img_dbg !array to hold colour values for roi IDs to output roi colourmap
      real(kind=real64), dimension(:), allocatable    ::      shuffle !temp array of random numbers for shuffeling the debug roi indicies.
      integer, dimension(:), allocatable    ::      indx_dbg !debug roi indicies

      type(MultipleGaussians2d)           ::      mg2d        !object for storing parameters for multiple 2D gaussians
      logical             ::      halfRoi, quarterRoi          !Booleans to control roi reduction if they are found to be oversized.

      integer             ::      nReducedMaxSum, nRawMaxSum

      nx = size(img, dim=1)
      ny = size(img, dim=2)

      !---    filter and smooth the input image
      !       determine the background noise level
      if (maxpix_auto) call estimateMaxpixFromDarkPixels(img, MAXPIX)

      print *, ""
      print *, "finding the noise level"
      print *, "^^^^^^^^^^^^^^^^^^^^^^^"
      call findSigma(img, MAXPIX, ss, bb, tt, ff)

      if (f_auto) then

         F_THRESH = 3!                   !   auto detects set signal to noise threshold to a sensible default value.

      end if
      detectFminOnS = F_THRESH

      print *, "pix over R&C threshold    ", count(img >= tt)/real(nx*ny, kind=real64)
      print *, "estimated signal to noise ", (tt - bb)/ss
      print *, "using detection level     ", F_THRESH, " sigma"

      print *, ""

      if (OPROI .or. OPRANDC) then
         allocate (img_dbg(3, 0:Nx - 1, 0:Ny - 1))
      end if

      !---    find the regions of interest
      allocate (indx(0:nx - 1, 0:ny - 1))
!

      print *, ""
      print *, "finding the raw maxima"
      print *, "^^^^^^^^^^^^^^^^^^^^^^"
      call rawMaxima2(img, bb + detectFminOnS*ss, indx, nRawMax, gte=.true.)
      fitted_g2d = MultipleGaussians2d_ctor(nRawMax)
      call setN(fitted_g2d, 0)

      print *, "raw maxima in image ", nRawMax

!
!
      print *, ""
      print *, "finding the regions of interest"
      print *, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
      if (roi_auto) then
         call nonPercolatingROI(img, bb, ss, detectFminOnS, ROIXMAX, ROIRMAX, indx, nroi, ok)
         if (.not. ok) then
           print *, "AnalyseSpots::fitSpots0 warning - can't find a suitable non-percolating threshold - increasing detection limit"
            detectFminOnS = F_THRESH*0.75d0
            call nonPercolatingROI(img, bb, ss, detectFminOnS, ROIXMAX, ROIRMAX, indx, nroi, ok)
            if (.not. ok) then
           print *, "AnalyseSpots::fitSpots0 warning - can't find a suitable non-percolating threshold - increasing detection limit"
               detectFminOnS = F_THRESH
               call nonPercolatingROI(img, bb, ss, detectFminOnS, ROIXMAX, ROIRMAX, indx, nroi, ok)
               if (.not. ok) then
                  print *, "AnalyseSpots::fitSpots0 error - can't find any suitable non-percolating threshold"
                  ROIRMAX = 2.0d0
                  print *, "reverting to R = ", ROIRMAX
                  call paddedROI(img, bb + ss*detectFminOnS, ROIRMAX, indx, nroi)
               end if
            end if
         end if

      else
         call paddedROI(img, bb + ss*detectFminOnS, ROIRMAX, indx, nroi)
      end if
      print *, "number of roi ", nroi
      print *, "detect threshold ", (bb + ss*detectFminOnS), int((bb + ss*detectFminOnS)*255)

      if (OPROI) then

         allocate (shuffle(nroi))
         call random_number(shuffle)
         allocate (indx_dbg(nroi))
         do ix = 1, nroi
            indx_dbg(ix) = ix
         end do
         call quicksort(shuffle, indx_dbg)
         deallocate (shuffle)
         xx = 1/real(max(1, nroi))
         img_dbg = 0

         do iy = 0, Ny - 1
            do ix = 0, Nx - 1
               if (indx(ix, iy) > 0) then
                  dd = indx_dbg(indx(ix, iy))*xx
                  img_dbg(1:3, ix, iy) = getRGB_double(dd)
               end if
            end do
         end do
         deallocate (indx_dbg)

         print *, "writing roi to """//trim(ofile)//".roi.png"""
         call write_rgb_png(trim(ofile)//".roi.png", img_dbg)
      end if

      if (OPRANDC) then
         img_dbg = 0
         do iy = 0, Ny - 1
            do ix = 0, Nx - 1
               dd = img(ix, iy)
               if (dd > ff) then
                  img_dbg(1:3, ix, iy) = getRGB_double(COLOURSCALE_MAGMA, ff)
               else if (dd > tt) then
                  img_dbg(1:3, ix, iy) = getRGB_double(COLOURSCALE_MAGMA, tt)
               else if (dd > bb) then
                  img_dbg(1:3, ix, iy) = getRGB_double(COLOURSCALE_MAGMA, bb)
               end if

            end do
         end do
         print *, "writing R&C to """//trim(ofile)//".randc.png"""
         call write_rgb_png(trim(ofile)//".randc.png", img_dbg)
      end if

      allocate (left(nroi))
      allocate (right(nroi))
      allocate (top(nroi))
      allocate (bottom(nroi))
      call dimensionsOfROI(indx, nroi, left, right, top, bottom, xmax, ymax)
      allocate (img_roi(xmax, ymax))
      allocate (indx_roi(xmax, ymax))
      print *, "AnalyseSpots::fitSpots0 info - largest roi ", xmax, "x", ymax

      !---    output locations of regions of interest
      if (OPROI .and. (DBG_ROIX > -1)) then
         do roi = 1, nroi

            if (indx(DBG_ROIX, DBG_ROIY) == roi) then
     print *, "roi*", roi, "x", left(roi), ":", right(roi), " y", bottom(roi), ":", top(roi), " px count ", count(indx(:, :) == roi)
            else
     print *, "roi ", roi, "x", left(roi), ":", right(roi), " y", bottom(roi), ":", top(roi), " px count ", count(indx(:, :) == roi)
            end if

         end do
      end if

      !---    find the spots
      print *, ""
      print *, "finding the spots"
      print *, "^^^^^^^^^^^^^^^^^"
      indx_roi = 0
      call cleanImage(img)
      DETECTFMIN = ss*detectFminOnS         !   convert intensity threshold back to intensity units.
      ng = 0
      nReducedMaxSum = 0
      nRawMaxSum = 0

      do roi = 1, nroi

         call progressBar(roi, nroi)

         call setLib_Maxima2d_dbg(.false.)
         call setLib_MultipleGaussians2d_dbg(.false.)
         if (DBG_ROIX > -1) then
            if (indx(DBG_ROIX, DBG_ROIY) == roi) then
               call setLib_Maxima2d_dbg(.true.)
               call setLib_MultipleGaussians2d_dbg(.true.)
            end if
         end if

         call extractROI(img, roi, indx, left, right, top, bottom, xMax, yMax, img_roi)
         if ((right(roi) < left(roi)) .or. (top(roi) < bottom(roi))) cycle

         call findMaxima(img_roi(1:xMax, 1:yMax), bb, bb + detectFminOnS*ss, indx_roi(1:xMax, 1:yMax), nRawMax, gte=.true.)

         nReducedMax = nRawMax
         if (nReducedMax > 10000) then
            print *, "AnalyseSpots::fitSpots0 WARNING - maximum fitting finds ", nReducedMax, " maxima in roi ", roi
            print *, "    this is very likely to be an error! Trying to continue with reduced resolution."
         else if (nReducedMax > 0) then
            nSpot = nReducedMax
            call mergeMaxima(img_roi(1:xMax, 1:yMax), bb, D_THRESH, indx_roi(1:xMax, 1:yMax), nSpot)
         else if (nReducedMax == 0) then
            cycle               !   no work to do.
         end if

         if ((right(roi) - left(roi) > ROIXMAX*8) .or. (top(roi) - bottom(roi) > ROIXMAX*8)) then
            print *, "AnalyseSpots::fitSpots0 warning - oversize roi x 8"
     print *, "roi ", roi, "x", left(roi), ":", right(roi), " y", bottom(roi), ":", top(roi), " px count ", count(indx(:, :) == roi)
            halfRoi = .true.; quarterRoi = .true.
            call shrinkRoi(img_roi, indx_roi, xMax, yMax, scaling=8)

         else if ((right(roi) - left(roi) > ROIXMAX*4) .or. (top(roi) - bottom(roi) > ROIXMAX*4)) then
            print *, "AnalyseSpots::fitSpots0 warning - oversize roi x 4"
     print *, "roi ", roi, "x", left(roi), ":", right(roi), " y", bottom(roi), ":", top(roi), " px count ", count(indx(:, :) == roi)
            halfRoi = .false.; quarterRoi = .true.
            call shrinkRoi(img_roi, indx_roi, xMax, yMax, scaling=4)

         else if ((right(roi) - left(roi) > ROIXMAX*2) .or. (top(roi) - bottom(roi) > ROIXMAX*2)) then
            print *, "AnalyseSpots::fitSpots0 warning - oversize roi x 2"
     print *, "roi ", roi, "x", left(roi), ":", right(roi), " y", bottom(roi), ":", top(roi), " px count ", count(indx(:, :) == roi)
            halfRoi = .true.; quarterRoi = .false.
            call shrinkRoi(img_roi, indx_roi, xMax, yMax, scaling=2)

         else
            halfRoi = .false.; quarterRoi = .false.
         end if

         if (nReducedMax > 10000) then
            !   haven't yet reduced the maxima here. Try again.
            call findMaxima(img_roi(1:xMax, 1:yMax), bb, bb + detectFminOnS*ss, indx_roi(1:xMax, 1:yMax), nRawMax)
if (nRawMax == 0) call findMaxima(img_roi(1:xMax, 1:yMax), bb, bb + detectFminOnS*ss, indx_roi(1:xMax, 1:yMax), nRawMax, gte=.true.)
            nReducedMax = nRawMax
            if (nReducedMax > 10000) then
          print *, "AnalyseSpots::fitSpots0 ERROR - maximum fitting on shrunk roi still finds ", nReducedMax, " maxima in roi ", roi
               stop
            else
          print *, "AnalyseSpots::fitSpots0 warning - maximum fitting on shrunk roi now finds ", nReducedMax, " maxima in roi ", roi
            end if
            !nSpot = nReducedMax
            call mergeMaxima(img_roi(1:xMax, 1:yMax), bb, D_THRESH, indx_roi(1:xMax, 1:yMax), nReducedMax)

         end if

         if (OPROI) then
            !   DEBUGGING- output regions of interest as .png
            xx = 1/real(max(1, nSpot))
            img_dbg = 0
            allocate (shuffle(nSpot))
            allocate (indx_dbg(nSpot))
            call random_number(shuffle)
            do ix = 1, nSpot
               indx_dbg(ix) = ix
            end do
            call quicksort(shuffle, indx_dbg)
            deallocate (shuffle)
            print *, "roi ", roi, xmax, ymax, minval(indx_roi(1:xMax, 1:yMax)), maxval(indx_roi(1:xMax, 1:yMax)), nSpot
            img_dbg = 0
            do iy = 1, yMax
               do ix = 1, xMax
                  if (indx_roi(ix, iy) > 0) then
                     dd = indx_dbg(indx_roi(ix, iy))*xx
                     img_dbg(1:3, ix - 1, iy - 1) = shadeColour(getRGB_double(COLOURSCALE_BGR, dd), img_roi(ix, iy), ff)
                  else if (indx_roi(ix, iy) == 0) then
                     img_dbg(1:3, ix - 1, iy - 1) = getRGB_double(COLOURSCALE_GREYSCALE, 0.25d0)
                  end if
               end do
            end do
            deallocate (indx_dbg)
            call write_rgb_png(trim(ofile)//".roidbg.png", img_dbg(1:3, 0:xMax - 1, 0:yMax - 1))
            !stop
         end if

         !---    new fitting code
         if (FITANDCOUNT) then
            mg2d = MultipleGaussians2d_ctor()
            call fit(mg2d, img_roi(1:xMax, 1:yMax), indx_roi(1:xMax, 1:yMax), bb, 1.0d0)
            call findT(mg2d, img_roi(1:xMax, 1:yMax), bb, 1.0d0, ss)
            call findQ(mg2d, img_roi(1:xMax, 1:yMax), bb, 1.0d0, ss)
            !nSpot = getN( mg2d,ignore=.true. )

            !---    new line 190523: check last minute ignore
            call setIgnore(mg2d)
            nSpot = getN(mg2d, ignore=.true.)

            if (OPROI .and. (DBG_ROIX > -1)) then
                        write(*,fmt='(10(a,i6))') "roi ",roi,"/",nroi," ",left(roi),":",right(roi),",",bottom(roi),":",top(roi)," spots = ",nSpot," from ",nReducedMax," peaks from ",nRawMax," maxima"
            end if

            nReducedMaxSum = nReducedMaxSum + nReducedMax
            nRawMaxSum = nRawMaxSum + nRawMax

            if (halfRoi) call xscale(mg2d, 2.0d0)
            if (quarterRoi) call xscale(mg2d, 4.0d0)
            call translate(mg2d, left(roi) - 1, bottom(roi) - 1)

            call add(fitted_g2d, mg2d)
            call delete(mg2d)
            ng = ng + nSpot
            call setN(fitted_g2d, ng)
         end if

      end do

      print *, ""
      print *, "analyseSpots::fitSpots0 info - "
      print *, "total raw max    ", nRawMaxSum
      print *, "total reduce max ", nReducedMaxSum
      print *, "total spot count ", ng
      print *, ""

      if (OPROI .or. OPRANDC) then
         deallocate (img_dbg)
      end if

      deallocate (left)
      deallocate (right)
      deallocate (top)
      deallocate (bottom)
      deallocate (img_roi)
      deallocate (indx_roi)
      deallocate (indx)

      return
   end subroutine fitSpots0

   !---

   subroutine cleanImage(img)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      set pixels with intensity zero to LIB_G2D_IGNORE
      !*      set pixels with intensity 1 with all neighbours intensity 1 to LIB_G2D_IGNORE
      real(kind=real64), dimension(0:, 0:), intent(inout)        ::      img !array holding image

      integer         ::      Nx, Ny   !x and y dimensions of img
      integer         ::      ix, iy   !x abd y pixel indices
      logical, dimension(:, :), allocatable        ::      tmp !temp array of booleans, which pixels to ignore.
      real(kind=real64)   ::      ff  !intensity of pixel of interest
      Nx = size(img, dim=1)
      Ny = size(img, dim=2)
      allocate (tmp(0:Nx - 1, 0:Ny - 1))
      tmp = .false.

      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            if (img(ix, iy) == 0) then
               tmp(ix, iy) = .true.
            else if (img(ix, iy) == 1) then !find pixels with intensity 1 with all neighbours intensity 1

               ff = img(max(0, ix - 1), iy)
               ff = min(ff, img(min(Nx - 1, ix + 1), iy))
               ff = min(ff, img(ix, max(0, iy - 1)))
               ff = min(ff, img(ix, min(Ny - 1, iy + 1)))
               tmp(ix, iy) = (ff == 1)
            end if
         end do
      end do

      where (tmp)
         img = LIB_G2D_IGNORE
      end where

      return
   end subroutine cleanImage

   subroutine shrinkRoi(img_roi, indx_roi, xMax, yMax, scaling)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(0:, 0:), intent(inout)            ::      img_roi     !array holding roi intensity values
      integer, dimension(0:, 0:), intent(inout)                      ::      indx_roi    !index each point as belonging to a single maximum indx(i,j) = 1:nMax
      integer, intent(inout)                                       ::      xMax, yMax   !maximum x and y dimensions of subject roi
      integer, intent(in)                                          ::      scaling     !factor to shrink roi by

      real(kind=real64), dimension(:), allocatable      ::      imgBarRoi       !average intensity within patch for each peak
      integer, dimension(:), allocatable                ::      nCountIndxRoi   !number of pixels within patch for each peak

      integer             ::      ix, iy, jx, jy     !x and y indices of the origional and scaled roi respectively
      integer             ::      nPxTot, ii       !total number of pixels within patch and scaled roi pixel of interest
      integer             ::      nReducedMax     !maximum possible index for shrunk roi based on origional.

      nReducedMax = maxval(indx_roi)
      allocate (imgBarRoi(0:nReducedMax))
      allocate (nCountIndxRoi(0:nReducedMax))
      do iy = 0, int((yMax - 1)/scaling)
         do ix = 0, int((xMax - 1)/scaling)
            nCountIndxRoi = 0                   !   number of pixels within patch for each peak
            imgBarRoi = 0                       !   average intensity within patch for each peak
            nPxTot = 0                          !   total number of pixels within patch
            do jy = iy*scaling, min(iy*scaling + scaling, yMax) - 1
               do jx = ix*scaling, min(ix*scaling + scaling, xMax) - 1
                  if (img_roi(jx, jy) == LIB_G2D_IGNORE) cycle
                  nPxTot = nPxTot + 1
                  ii = indx_roi(jx, jy)
                  nCountIndxRoi(ii) = nCountIndxRoi(ii) + 1
                  imgBarRoi(ii) = imgBarRoi(ii) + img_roi(jx, jy)
               end do
            end do
            indx_roi(ix, iy) = 0                 !   note: this is safe because the next patch will not look at (jx,jy)
            img_roi(ix, iy) = LIB_G2D_IGNORE
            do ii = 1, nReducedMax
               if (nCountIndxRoi(ii)*2 > (nPxTot - nCountIndxRoi(0))) then !   most pixels within roi. Lets assume is a roi patch
                  indx_roi(ix, iy) = ii
                  img_roi(ix, iy) = imgBarRoi(ii)/nCountIndxRoi(ii)      !   safe as nCountIndxRoi > 0
               end if
            end do
         end do
      end do
      deallocate (imgBarRoi)
      deallocate (nCountIndxRoi)

      xMax = ceiling(real(xMax)/scaling)
      yMax = ceiling(real(yMax)/scaling)

      return
   end subroutine shrinkRoi

end module analyseSpots

