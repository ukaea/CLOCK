
program linearFeatures
!---^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    linearFeatures from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      A simple test program to compute radon transformation
!*
!*      The image size is img(1:Nx,1:Ny)
!*      From this we construct a buffered image with an extra pixel all round the outside, ie buffered_img(0:Nx+1,0:Ny+1)
!*      Then I can ask for the value at any interpolated point from 0:Nx+1

   use Lib_CommandLineArguments
   use Lib_Png
   use Lib_Filenames
   use Lib_LinearFeatureDetect
   use Lib_Radon
   use Lib_FFTW3f
   use iso_fortran_env
   implicit none

   character(len=*), parameter      ::      VERSION = "1.0.0"

   !---    input parameters
   type(CommandLineArguments)      ::      cla             !object to hold command line arguments (CLA)
   character(len=256)              ::      filename = ""   !See get(cla...) calls below
   character(len=256)              ::      outfile = ""
   logical                         ::      negative = .false.
   integer                         ::      n = 0      !   block size
   logical                         ::      reduce = .true.
   logical                         ::      linesOnly = .false.
   logical                         ::      greybg = .true.
   integer                         ::      nLoops = 1    !   number of AIC min loops to try

   logical                         ::      lineFile = .false.  !print CAD style text file containing vertices and line data
   character(len=256)              ::      lineFilename = ""   !What to call the line file
   character(len=2)                ::      LENGTHUNIT = "px"     !length unit (default pixels)
   integer                         ::      nVerts = 0    !number of vertices

   !---    key physical properties

   real(kind=real64), dimension(:, :), allocatable    ::      img     !image array
   integer                                         ::      Nx, Ny   !x & y dimesnions of img/pixels
   type(LinearFeatureDetect)                       ::      lfd     !object to hold image information & detected features
   integer         ::      ii                                      !loop counter

   logical                 ::  ok1, ok2                             !booleans controlling vertex and line reduction respectively.
   !---jph hack start---
   real(kind=real64), dimension(:), allocatable         ::      cprime, rpf
   real(kind=real64), dimension(:, :), allocatable         ::      fftout
   character(len=256)              ::      tempFilename
   !---jph hack end-----

   !---    read command line arguments
   cla = CommandLineArguments_ctor(35)

   call setProgramDescription(cla, "linearFeatures" &
                              //"\n    A simple test program to compute radon transformations to detect linear features. \n")
   call setProgramVersion(cla, VERSION)

   call get(cla, "f", filename, LIB_CLA_REQUIRED, "              input filename")
   outfile = trim(removeSuffix(filename))//".recon.png"
   call get(cla, "o", outfile, LIB_CLA_OPTIONAL, "            output filename")
   call get(cla, "negative", negative, LIB_CLA_OPTIONAL, "     invert image")
   call get(cla, "n", n, LIB_CLA_OPTIONAL, "            block size (px)")
   call get(cla, "r", reduce, LIB_CLA_OPTIONAL, "            reduce complexity by combining vertices")
   call get(cla, "g", LIB_LFD_GRID, LIB_CLA_OPTIONAL, "            draw gridlines in output image")
   call get(cla, "b", greybg, LIB_CLA_OPTIONAL, "            use optimal greyscale background in output image")
   call get(cla, "l", linesOnly, LIB_CLA_OPTIONAL, "            skinny lines only in output image")
   call get(cla, "lf", lineFile, LIB_CLA_OPTIONAL, "            Output text file with lines/vertices")
   lineFilename = trim(removeSuffix(filename))//".lines.txt"
   call get(cla, "ln", lineFilename, LIB_CLA_OPTIONAL, "            name of lines/vertices file")

   call get(cla, "nLoops", nLoops, LIB_CLA_OPTIONAL, "       maximum number of AIC minimisation loops to try")
   if (hasArgument(cla, "nLoops")) reduce = (nLoops > 0)
   call report(cla)
   if (.not. allRequiredArgumentsSet(cla)) stop
   if (hasHelpArgument(cla)) stop
   call delete(cla)

   !---    welcome message
   print *, "linearFeatures.exe v."//VERSION
   print *, "^^^^^^^^^^^^^^^^^^^^^"//repeat("^", len_trim(VERSION))
   print *, ""
   print *, "   filename                """//trim(filename)//""""
   print *, "   output file             """//trim(outfile)//""""
   print *, "   negative                ", negative
   print *, "   block size              ", n
   print *, "   reduce complexity       ", reduce
   print *, "   draw gridlines          ", LIB_LFD_GRID
   print *, "   draw optimal greyscale  ", greybg
   print *, "   lines only              ", linesonly
   print *, "   Write line/vertex data to file  ", lineFile
   print *, "   Line/vertex file name   """//trim(lineFilename)//""""
   print *, ""

   !---    read input .png
   print *, "reading input .png "
   call readPng(filename, img, negative)
   Nx = size(img, dim=1)
   Ny = size(img, dim=2)
   print *, "image size ", Nx, "x", Ny
   if (n == 0) n = max(Nx, Ny)
   print *, ""
   allocate (fftout(((2*(Nx/2 + 1) - 1)), Ny))
   allocate (rpf(0:100))
   allocate (cprime(0:100))
   call FFT2d(img, fftout)
   call radialPowerSpectrum(img, (2*3.141592654d0/max(Nx, Ny)), 3.141592654d0, cprime, rpf)
   tempFilename = trim(outfile)//'.fft.png'
   call write_greyscale_png(tempFilename, fftout, negative=.false.)
   lfd = LinearFeatureDetect_ctor(img, n)
   call report(lfd)

   if (reduce) then
      do ii = 1, nLoops
         call reduceVertices(lfd, ok1)
         call reduceLines(lfd, ok2)
         call report(lfd)
         if (.not. (ok1 .or. ok2)) exit
      end do
   end if

   !---    dbg double output before and after minimsation - remove conditional for final version
   if (reduce) then
      print *, ""
      print *, "linearFeatures.exe info - writing reconstructed image to """//trim(outfile)//""""
      LIB_LFD_GREYBG = greybg
      call reconstructImage(lfd, thin=linesOnly)
      call reconstructImage(lfd, img)
      if (linesOnly) img = max(0.0d0, min(1.0d0, img))
      if (LIB_LFD_GREYBG) then
         call writePng(outfile, img, normalise=.true., negative=negative)
      else
         call writePng(outfile, img, negative=negative)
      end if
   end if
   !---

   if (lineFile) then !write the line/vertex data in an eisier to use format.
      nVerts = getNoofVerts(lfd, nVerts)
      open (unit=521, file=lineFilename, action="write")
      call printVertsInfo(lfd, 521, LENGTHUNIT)

      close (unit=521)
   end if

   do ii = 0, 100
      print *, "rpf,cprime", rpf(ii), cprime(ii)
   end do

   print *, ""
   print *, "done"
   print *, ""

end program linearFeatures
