
program AnalysePngHistogram
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
!*    A simple test program to find whether a .png file appears to be
!*    a dark background with a bright signal, or a bright background with a dark signal
!*
!*      This is done by analysing the histogram of pixel intensities as a sum of Gaussians.
!*
!*      to determine the probability that a fit is correct:
!*          p_(unnormalised) = Exp[ (aic_0 - aic_j)/2 ]
!*      where aic_0 is the best aic value we find and aic_j is the value for the jth trial set of Gaussians.
!*
!*
!*      We imagine that the background is one Gaussian contribution to the histogram
!*      and that the signal is the other Gaussians.
!*      ( note: the background may not be > 50% of the signal, rather I will say that the probability that peak j is background is proportional to the area under j )
!*
!*      We then say that a secondary peak gives a weighting to dark signal given by its area left of the mean of peak j and a weighting to bright signal
!*      given by its area to the right of peak j.
!*      This gives a probability for each fit to the data that each part of the signal is bright or dark,
!*      multiplied together with the probability that the total fit is correct, and that peak j is the background.
!*      Sum these together to give a total probability that the signal is bright or dark.
!*

   use Lib_CommandLineArguments
   use Lib_FitGaussians1d
   use Lib_Png
   use iso_fortran_env
   implicit none

   character(len=*), parameter      ::      VERSION = "1.0.0"

   !---    input parameters
   type(CommandLineArguments)      ::      cla
   character(len=256)              ::      filename = ""
   integer                         ::      nBins = 256         !   number of bins in histogram
   integer                         ::      mMax = 10           !   max number of Gaussians to try

   !---    key physical properties
   real(kind=real64), dimension(:, :), allocatable    ::      img !array holding the image
   integer                                         ::      Nx, Ny   !Dimensions of img
   real(kind=real64), dimension(:), allocatable      ::      hist    !histogram of image intensities
   type(FitGaussians1d), dimension(:), allocatable   ::      fg1d    !Object holding fit information
   real(kind=real64), dimension(:), allocatable      ::      aic     !Akaike infomation criterion
   integer, parameter                               ::      DARK_BACK = -1 !integers used to describe image background signal
   integer, parameter                               ::      MIXED_BACK = 0
   integer, parameter                               ::      BRIGHT_BACK = 1

   !---    dummies
   integer                         ::      ix, iy, mm, jj, kk
   integer                         ::      hh
   real(kind=real64)               ::      ff, aic0, pp, qq, xlo, xhi
   real(kind=real64)               ::      areasum

   real(kind=real64), dimension(:, :), allocatable        ::      fmusigax
   real(kind=real64), dimension(-1:1)                   ::      ppsum
   !--Integers--
   !ix: x address of image
   !iy: y address of image
   !mm: fit index
   !jj: background gaussian index
   !kk: non background gaussian index
   !hh: histogram instenity bin value

   !--Reals--
   !ff: pixel intensity value of img at (ix, iy)
   !aic0: lowest AIC value
   !pp: probability that given fit is "correct"
   !qq: probability that given gaussian is the background
   !xlo: integral from 0 to mean of background for a given gaussian
   !xhi, integral from meanof background to max intensity for a given gaussian
   !areasum: sum of areas

   !--Real arrays--
   !fmusigax: array holding gaussian parameters
   !ppsum: probability weighted area sum

   !---    read command line arguments
   cla = CommandLineArguments_ctor(30)

   call setProgramDescription(cla, "AnalysePngHistogram.exe" &
        //"\n    A simple test program to find whether a .png file appears to be  \n    a dark background with a bright signal, or a bright background with a dark signal" )
   call setProgramVersion(cla, VERSION)

   call get(cla, "f", filename, LIB_CLA_REQUIRED, "            input filename")
   call get(cla, "n", nBins, LIB_CLA_OPTIONAL, "            number of histogram bins")
   call get(cla, "m", mMax, LIB_CLA_OPTIONAL, "            number of gaussians to try")

   call report(cla)
   if (.not. allRequiredArgumentsSet(cla)) stop
   if (hasHelpArgument(cla)) stop
   call delete(cla)

   !---    welcome message
   print *, "AnalysePngHistogram.exe v."//VERSION
   print *, "^^^^^^^^^^^^^^^^^^^^^^^^^^"//repeat("^", len(VERSION))
   print *, ""
   print *, "   filename                """//trim(filename)//""""
   print *, "   intensity bins          ", nBins
   print *, "   max gaussians           ", mMax
   print *, ""

   !---    read input .png
   call readPng(filename, img)
   Nx = size(img, dim=1)
   Ny = size(img, dim=2)
   print *, "AnalysePngHistogram.exe info - image size ", Nx, "x", Ny
   print *, ""

   !---    construct a histogram
   allocate (hist(0:nBins - 1))
   hist = 0.0d0
   do iy = 1, Ny
      do ix = 1, Nx
         ff = img(ix, iy)
         hh = floor(nBins*ff)
         hh = max(0, min(nBins - 1, hh))
         hist(hh) = hist(hh) + 1
      end do
   end do
   hist = hist/(Nx*Ny)

   !---    it is generally a good idea to exclude the intensity 0 and 1 boxes, as these could contain
   !       salt and pepper noise. Unless, of course, the data is computer generated, and all black/white?
   !       Not sure there's a completely assumption-free thing to do here. I'll see whether the majority of the image is black/white...?
   if (hist(0) + hist(nBins - 1) < 0.5d0) then
      hist(0) = max(0.0d0, min(1.0d0, hist(1) - (-3*hist(1) + 4*hist(2) - hist(3))/2))
      hist(nBins - 1) = max(0.0d0, min(1.0d0, hist(nBins - 2) + (3*hist(nBins - 2) - 4*hist(nBins - 3) + hist(nBins - 4))/2))
   end if

   !---    fit multiple guassians to the histogram
   allocate (fg1d(mMax))
   call fit(hist, mMax, fg1d)

   !---    find AIC criterion "best" fit and report it.
   allocate (aic(mMax))
   do mm = 1, mMax
      aic(mm) = AICc(fg1d(mm))
   end do
   mm = minloc(aic, dim=1)                     !    "best" fit
   aic0 = aic(mm)
   print *, ""
   print *, "AnalysePngHistogram.exe info - best multiple gaussian result"
   call report(fg1d(mm), verbose=.true.)
   print *, ""

   !---    is it dark signal on bright or vice versa?
   allocate (fmusigax(5, mMax))                              !   intensity, mean, width, area, centre-of-mass as an array

   ppsum = 0                                               !   probability of DARK/MIXED/BRIGHT background level

   do mm = 1, mMax                              !   for each fit in the stack

      pp = exp((aic0 - aic(mm))/2)          !   this is the probability that this fit is "correct"

      if (getM(fg1d(mm)) == 1) then
         !   single gaussian fit = all background.
         ppsum(MIXED_BACK) = ppsum(MIXED_BACK) + pp

      else
         !   multiple gaussian fit = need to work out which bits are bright signal vs dark signal

         !---    extract all the data about the gaussians for this trial fit
         areasum = 0                             !   will normalise the weighting of each gaussian according to its area
         do jj = 1, getM(fg1d(mm))
            fmusigax(1:5, jj) = getFmuSigAx(fg1d(mm), jj)
            areasum = areasum + fmusigax(4, jj)
         end do

         !---    now lets assume gaussian jj is the background, and the others are foreground
         do jj = 1, getM(fg1d(mm))

            qq = fmusigax(4, jj)/max(1.0d-16, areasum)        !   probability that jj is the background

            do kk = 1, getM(fg1d(mm))
               if (kk == jj) cycle           !   consider the other gaussians in the set

               xlo = intgrl(fg1d(mm), kk, fmusigax(2, jj))      !   integral from 0: mean of background
               xlo = xlo/max(1.0d-16, fmusigax(4, kk))       !   fraction of gaussian below mean of background
               xhi = 1 - xlo

               ppsum(BRIGHT_BACK) = ppsum(BRIGHT_BACK) + qq*pp*xlo
               ppsum(DARK_BACK) = ppsum(DARK_BACK) + qq*pp*xhi

            end do
         end do

      end if
   end do

   pp = sum(ppsum)
   ppsum = ppsum/pp
   print *, ""

   !---    output the solution
   print *, "p(dark,mixed,bright) = ", ppsum
   if (ppsum(DARK_BACK) > max(ppsum(MIXED_BACK), ppsum(BRIGHT_BACK))) then
      write (*, fmt='(a,f6.2,a,f10.3)') """"//trim(filename)//""" Dark background ", ppsum(DARK_BACK)*100, "% AIC (best) ", aic0
   else if (ppsum(BRIGHT_BACK) > max(ppsum(MIXED_BACK), ppsum(DARK_BACK))) then
      write (*, fmt='(a,f6.2,a,f10.3)') """"//trim(filename)//""" Bright background ", ppsum(BRIGHT_BACK)*100, "% AIC (best) ", aic0
   else
      write (*, fmt='(a,f6.2,a,f10.3)') """"//trim(filename)//""" Mixed background ", ppsum(MIXED_BACK)*100, "% AIC (best) ", aic0
   end if

   !---    bye bye
   do mm = 1, mMax
      call delete(fg1d(mm))
   end do
   deallocate (fg1d)
   print *, ""
   print *, "done"
   print *, ""

end program AnalysePngHistogram
