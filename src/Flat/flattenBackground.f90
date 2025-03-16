
program flattenBackground
!---^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    flattenBackground from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      utility program designed to help feature detection in greyscale images
!*      by removing background fluctuations and performing a couple of simple cleaning tasks

   use iso_fortran_env
   use Lib_CommandLineArguments
   use Lib_LocalFilter2d
   use Lib_RidlerCalvard
   use Lib_LowPassFilter3d
   use Lib_MaxLikelihoodFilter
   use Lib_Png
   use Lib_FFTW3f
   use Lib_DataDoubler
   use Lib_Filenames
   implicit none

   !---    constants
   real(kind=real64), parameter                     ::      PI = 3.141592654d0

   !---    command line arguments (CLA)
   character(len=256)                  ::      filename = "" !see "get(cla, ...)" calls below
   character(len=256)                  ::      outfile = ""
   logical                             ::      zeroIsUnset = .true.
   logical                             ::      sandp = .false.
   logical                             ::      negative = .false.
   logical                             ::      half = .false.
   logical                             ::      flat = .true.
   logical                             ::      retone = .false.
   logical                             ::      dbg = .false.

   logical                             ::      quarter = .false.
   logical                             ::      rescale = .false.
   real(kind=real64)                   ::      lambda = LIB_CLA_NODEFAULT_R
   integer                             ::      vRC_in = LIB_CLA_NODEFAULT_I

   integer                             ::      nGB = 20            !   number of gaussian blurs used
   integer                             ::      nLoops = 2
   integer                             ::      power = 1
   integer                             ::      hpf_coeff = 0
   type(CommandLineArguments)          ::      cla
   character(len=8)                    ::      VERSION = "1.0"

   !---    key physical quantities and arrays
   real(kind=real64), dimension(:, :), allocatable    ::      img             !   image input
   integer                                         ::      Nx, Ny           !   size of image in pixels
   integer                                         ::      Mx, My           !   size of scaled working image
   real(kind=real64), dimension(:, :), allocatable    ::      img_out         !   image output
   real(kind=real64)                               ::      stddev          !   image standard deviation
   real(kind=real64)                               ::      sigma           !   feature size
   real(kind=real64), dimension(:, :), allocatable    ::      img_work        !   working image (scaled and filtered)
   real(kind=real64), dimension(:, :), allocatable    ::      img_back        !   background image
   real(kind=real64), dimension(:, :), allocatable    ::      img_fore        !   foreground image

   real(kind=real64), dimension(:, :), allocatable    ::      bg_node, sd_node, fg_node         !   background level on coarse nodes (background, background std.dev and fore ground respectively)

   real(kind=real64), dimension(:, :), allocatable    ::      logistic        !   logistic function switch between 1=background and 0=foreground

   real(kind=real64), dimension(:, :, :), allocatable  ::      blur_stack      !Array to hold several images that are increasingly blured

   integer                             ::      vRC                        !See CLA
   real(kind=real64)                   ::      b0, f0                      !after filtered background and foreground intensity

   !---    dummy
   integer                             ::      nn                          !Number of zero intensity pixels
   real(kind=real64)                   ::      area                        !iamge area/pixel squared
   real(kind=real64)                   ::      bb, th, ff, bstd, fbar, snr, fot      !   Ridler-Calvard values, back,thresh,fore, back std dev,<f>,signal-to-noise, frac over thresh
   real(kind=real64)                   ::      q_min, q_max                     !   FT values, min max q values
   real(kind=real64)                   ::      rpfsum, qrpfbar, q2rpfbar         !   radial power function variables
   real(kind=real64), dimension(:), allocatable      ::      rpf                 !   radial power function array

   integer                             ::      Cx, Cy                       !   coarse grid for computing bg level and sd dev
   integer                             ::      ix, iy, jx, jy, kx, ky           !x and y indicies for image pixels, sub image pixels and WHAT IS K???
   real(kind=real64)                   ::      ax, ay                       !distance to background node centre why is this a vector???
   real(kind=real64)                   ::      dsigma, ww, wsum              !spreading sigma, lengthscale weighted average background level and its sum.
   integer                             ::      loop                        !radial power function loop counter
   real(kind=real64)                   ::      img_min, img0avg, img_irange !min, averge and range of image intensity

   real(kind=real64), dimension(2)      ::      gg                         !Varible in logistic function???

   real(kind=real64), dimension(:, :), allocatable      ::      x            !Temporary varibles for removing fourier coefficiants???
   real(kind=real64), dimension(:), allocatable        ::      dat
   real(kind=real64), dimension(:, :, :), allocatable    ::      qdat

   !---

   cla = CommandLineArguments_ctor(50)
   call setProgramDescription(cla, "flattenBackground.exe " &
            //"\n    utility program designed to help feature detection in greyscale images \n    by removing background fluctuations and performing a couple of simple cleaning tasks" )
   call setProgramVersion(cla, VERSION)

   call get(cla, "f", filename, LIB_CLA_REQUIRED, "                    filename")
   outfile = trim(removeSuffix(filename))//".flat.png"
   call get(cla, "o", outfile, LIB_CLA_OPTIONAL, "                  output filename")
        call get( cla,"zero",zeroIsUnset,LIB_CLA_OPTIONAL,"               interpret pixels with intensity zero as 'unset' rather than just very dark" )     
   call get(cla, "sandp", sandp, LIB_CLA_OPTIONAL, "              apply salt and pepper filter")
   call get(cla, "lambda", lambda, LIB_CLA_OPTIONAL, "             set range for max likelihood filter")
   call get(cla, "negative", negative, LIB_CLA_OPTIONAL, "           invert image ")
   call get(cla, "half", half, LIB_CLA_OPTIONAL, "               compute background using 1/2 size image (default T if Nx>2000)")
   call get(cla, "quarter", quarter, LIB_CLA_OPTIONAL, "            compute background using 1/4 size image (default T if Nx>4000)")
   call get(cla, "g", nGB, LIB_CLA_OPTIONAL, "                  number of gaussian blurs to use for lengthscale estimate")
   call get(cla, "n", nLoops, LIB_CLA_OPTIONAL, "                  number of iterations")
   call get(cla, "p", power, LIB_CLA_OPTIONAL, "                  power law for gradient contributions")
   call get(cla, "flat", flat, LIB_CLA_OPTIONAL, "               flatten image using scale invariant background subtraction")
   call get(cla, "retone", retone, LIB_CLA_OPTIONAL, "             retone image to preserve foreground/background intensity")
   rescale = hasArgument(cla, "half") .or. hasArgument(cla, "quarter")
   call get(cla, "hpf", hpf_coeff, LIB_CLA_OPTIONAL, "                high pass filter - fourier coefficients to remove")
   call get(cla, "vRC", vRC_in, LIB_CLA_OPTIONAL, "                minimum variable-Ridler-Calvard length scale")
   call get(cla, "dbg", dbg, LIB_CLA_OPTIONAL, "                output intermediate images for debugging")

   if (getSuffix(outfile) /= "png") outfile = trim(outfile)//".png"

   call report(cla)
   if (.not. allRequiredArgumentsSet(cla)) stop
   if (hasHelpArgument(cla)) stop
   call delete(cla)

   if (.not. flat) nLoops = 0

   !---    welcome
   print *, "flattenBackground.exe"
   print *, "^^^^^^^^^^^^^^^^^^^^^"
   print *, ""

   print *, "       filename                """//trim(filename)//""""
   print *, "       output file             """//trim(outfile)//""""
   print *, "       treat zero as unset     ", zeroIsUnset
   print *, "       salt and pepper         ", sandp
   if (lambda /= LIB_CLA_NODEFAULT_R) then
      print *, "       max likelihood lambda   ", lambda
   else
      print *, "       no max likelihood filter"
   end if
   print *, "       number Gaussian blurs   ", nGB
   print *, "       high pass filter coeffs ", hpf_coeff
   print *, "       flatten background      ", flat
   print *, "       retone                  ", retone
   print *, "       variable Ridler-Calvard ", vRC_in
   print *, "       power law               ", power
   print *, "       flattening loops        ", nLoops
   print *, ""

   !---    load in input file
   call read_greyscale_png(filename, img)
   Nx = size(img, dim=1)
   Ny = size(img, dim=2)
   if (.not. rescale) then
      if (min(Nx, Ny) > 8000) then
         half = .true.; quarter = .true.
      else if (min(Nx, Ny) > 4000) then
         quarter = .true.
      else if (min(Nx, Ny) > 2000) then
         half = .true.
      end if
   end if
   if (half) then
      if (quarter) then
         print *, "       compute background      1/8"
      else
         print *, "       compute background      1/2  "
      end if
   else if (quarter) then
      print *, "       compute background      1/4  "
   end if
   print *, ""

   print *, "flattenBackground.exe info - read file successfully, size ", Nx, ",", Ny
   area = Nx*Ny
   allocate (img_out(0:Nx - 1, 0:Ny - 1))
   print *, ""

   !---    first step is to tidy up the file a bit
   if (zeroIsUnset) then
      if (sandp) then
         print *, ""
         print *, "Salt and pepper filter overides 'zero', not converting pixels with intensity zero to 'unset'"
      else

         print *, ""
         print *, "converting pixels with intensity zero to 'unset'"
         where (img == 0)
            img = LOCALFILTER2D_IGNORE
         end where
         nn = count(img == LOCALFILTER2D_IGNORE)
         print *, "count of dead pixels ", nn, "  (", nn*100.0/area, "%)"
         print *, ""
      end if
   end if

   if (negative) then
      print *, ""
      print *, "creating negative image"
      img = 1 - img
      print *, ""
   end if

   if (sandp) then
      print *, ""
      print *, "salt and pepper filter"
      call setLocalFilter2d_sigma(0.707d0)
      call denoise(img, 3.0d0, img_out(0:Nx - 1, 0:Ny - 1), stddev)
      img = min(1.0d0, max(0.0d0, img_out(0:Nx - 1, 0:Ny - 1)))
      nn = count(img == LOCALFILTER2D_IGNORE)
      print *, "count of dead pixels ", nn, "  (", nn*100.0/area, "%)"
      print *, ""
   end if

   if (lambda /= LIB_CLA_NODEFAULT_R) then
      print *, ""
      print *, "max likelihood filter"
      if (zeroIsUnset) then
         if (negative) then
            where (img == LOCALFILTER2D_IGNORE)
               img = 1
            end where
         else
            where (img == LOCALFILTER2D_IGNORE)
               img = 0
            end where
         end if
      end if

      if ((lambda /= LIB_CLA_NODEFAULT_R) .and. (lambda > 0)) then

         print *, ""
         print *, "log likelihood filtering file with lambda = ", lambda
         print *, ""
         call findImageIntensityFeatures(img, bb, th, ff, bstd, fbar, snr, fot)
         th = bstd/lambda
         print *, "   before filter: back, thresh, fore     ", bb, th, ff
         print *, "                  back std, <f>,snr,fot  ", bstd, fbar, snr, fot
         call maxLikelihoodFilter(img, img_out, bstd, th)
         img = min(1.0d0, max(0.0d0, img_out(0:Nx - 1, 0:Ny - 1)))

      end if

      print *, ""
   end if

   if (hpf_coeff > 0) then

      print *, ""
      print *, "high pass filter - removing ", hpf_coeff, " fourier coefficients"
      print *, ""

      allocate (x(2, nx*ny))
      allocate (dat(nx*ny))
      do iy = 1, ny
         do ix = 1, nx
            x(1:2, ix + nx*(iy - 1)) = (/ix, iy/)
            dat(ix + nx*(iy - 1)) = img(ix, iy)
         end do
      end do
      allocate (qdat(4, 0:hpf_coeff, 0:hpf_coeff))
      call generateFourierCoefficients(x, dat, real((/nx, ny/), kind=real64), (/hpf_coeff, hpf_coeff/), qdat)
      deallocate (x)
      deallocate (dat)
      img_out = interpolateFromFourierCoefficients(nx, ny, (/hpf_coeff, hpf_coeff/), qdat)
      deallocate (qdat)
      img0avg = sum(img)/(nx*ny)
      img_min = minval(img)
      img_irange = maxval(img) - img_min
      if (img_irange > 0) img_irange = 1/img_irange

      img = (img - img_out)*img_irange + (img0avg - img_min)*img_irange

      img = min(1.0d0, max(0.0d0, img))

   end if

   call findImageIntensityFeatures(img, b0, th, f0, bstd, fbar, snr, fot)
   print *, "   after filters: back, thresh, fore     ", b0, th, f0
   print *, "                  back std, <f>,snr,fot  ", bstd, fbar, snr, fot

   !stop

   if (flat) then

      !---    scale image?
      allocate (img_work(0:Nx - 1, 0:Ny - 1))
      allocate (img_back(0:Nx - 1, 0:Ny - 1))
      if (half) then
         print *, ""

         if (quarter) then
            print *, "scaling image x 1/8"
            call halfImage(Nx, Ny, img, img_work, pbc=.false.)
            Mx = int((Nx + 1)/2); My = int((Ny + 1)/2)
            call halfImage(Mx, My, img_work, img_out, pbc=.false.)
            Mx = int((Mx + 1)/2); My = int((My + 1)/2)
            call halfImage(Mx, My, img_out, img_work, pbc=.false.)
            Mx = int((Mx + 1)/2); My = int((My + 1)/2)
         else
            print *, "scaling image x 1/2"
            call halfImage(Nx, Ny, img, img_work, pbc=.false.)
            Mx = int((Nx + 1)/2); My = int((Ny + 1)/2)
         end if
      else if (quarter) then
         print *, "scaling image x 1/4"
         call halfImage(Nx, Ny, img, img_out, pbc=.false.)
         Mx = int((Nx + 1)/2); My = int((Ny + 1)/2)
         call halfImage(Mx, My, img_out, img_work, pbc=.false.)
         Mx = int((Mx + 1)/2); My = int((My + 1)/2)
      else
         img_work = img
         Mx = Nx; My = Ny
      end if
      print *, "working image size ", Mx, ",", My

      allocate (logistic(0:Nx - 1, 0:Ny - 1))
      allocate (img_fore(0:Mx - 1, 0:My - 1))
      allocate (blur_stack(0:Mx - 1, 0:My - 1, 0:nGB))
      allocate (rpf(0:100))
      img_back(0:Mx - 1, 0:My - 1) = img_work(0:Mx - 1, 0:My - 1)
      img_fore(0:Mx - 1, 0:My - 1) = img_work(0:Mx - 1, 0:My - 1)

      do loop = 1, nLoops !Main flattening loop starting here

         !---    determine feature size roughly using fourier transform

         print *, ""
         print *, "determining feature size using radial power function"
         call findImageIntensityFeatures(img_fore(0:Mx - 1, 0:My - 1), bb, th, ff, bstd, fbar, snr, fot)
         print *, "   loop ", loop
         print *, "       : back, thresh, fore     ", bb, th, ff
         print *, "         back std, <f>,snr,fot  ", bstd, fbar, snr, fot
         q_max = PI
         q_min = 2*PI/max(Mx, My)

         call radialPowerFunction(img_fore(0:Mx - 1, 0:My - 1), q_min, q_max, rpf, rpfsum, qrpfbar, q2rpfbar, snr)
         sigma = sqrt(PI)/(2*qrpfbar)
         if (half) then
            if (quarter) then
               print *, "single Gaussian sigma    = ", sigma*8
            else
               print *, "single Gaussian sigma    = ", sigma*2
            end if
         else if (quarter) then
            print *, "single Gaussian sigma    = ", sigma*4
         else
            print *, "single Gaussian sigma    = ", sigma
         end if
         print *, ""
         vRC = ceiling(sigma*2)       !   std dev assumed to vary on scale 4 x characteristic scale of image
         if (vRC_in /= LIB_CLA_NODEFAULT_I) then
            print *, "variable Ridler-Calvard pixels estimated by FFTW ", vRC, " using input value ", vRC_in
            vRC = vRC_in
         end if

         !---    the flattening proceeds by computing the probability that each pixel is part of the background
         !       which requires an estimate of the background level and std dev for each pixel
         !       compute background and std dev of background on a course node grid
         Cx = max(1, Mx/vRC)
         Cy = max(1, My/vRC)
         print *, "variable Ridler-Calvard pixels in scaled image ", vRC, " grid nodes required ", Cx, ",", Cy
         allocate (bg_node(-1:Cx, -1:Cy))      !   defined -1:Cx for convenience- I'm going to extrapolate the edges.
         allocate (sd_node(-1:Cx, -1:Cy))      !
         allocate (fg_node(-1:Cx, -1:Cy))      !   b
         jx = ceiling(Mx*1.0/Cx)           !   number of pixels per block
         jy = ceiling(My*1.0/Cy)

         do ky = 0, Cy - 1
            do kx = 0, Cx - 1
 call findImageIntensityFeatures( img_fore( kx*jx:min(kx*jx+jx-1,Mx-1) , ky*jy:min(ky*jy+jy-1,My-1) ) , bb,th,ff,bstd,fbar,snr,fot )
               bg_node(kx, ky) = bb
               sd_node(kx, ky) = bstd
               fg_node(kx, ky) = ff
            end do
         end do

         !---    finite difference extrapolation to borders.
         if (Cy > 3) then
            do kx = 0, Cx - 1
               bg_node(kx, -1) = (7*bg_node(kx, 0) - 4*bg_node(kx, 1) + bg_node(kx, 2))/4
               bg_node(kx, Cy) = (7*bg_node(kx, Cy - 1) - 4*bg_node(kx, Cy - 2) + bg_node(kx, Cy - 3))/4
               sd_node(kx, -1) = (7*sd_node(kx, 0) - 4*sd_node(kx, 1) + sd_node(kx, 2))/4
               sd_node(kx, Cy) = (7*sd_node(kx, Cy - 1) - 4*sd_node(kx, Cy - 2) + sd_node(kx, Cy - 3))/4
               fg_node(kx, -1) = (7*fg_node(kx, 0) - 4*fg_node(kx, 1) + fg_node(kx, 2))/4
               fg_node(kx, Cy) = (7*fg_node(kx, Cy - 1) - 4*fg_node(kx, Cy - 2) + fg_node(kx, Cy - 3))/4
            end do
         else
            do kx = 0, Cx - 1
               bg_node(kx, -1) = bg_node(kx, 0)
               bg_node(kx, Cy) = bg_node(kx, Cy - 1)
               sd_node(kx, -1) = sd_node(kx, 0)
               sd_node(kx, Cy) = sd_node(kx, Cy - 1)
               fg_node(kx, -1) = fg_node(kx, 0)
               fg_node(kx, Cy) = fg_node(kx, Cy - 1)
            end do
         end if
         if (Cx > 3) then
            do ky = -1, Cy
               bg_node(-1, ky) = (7*bg_node(0, ky) - 4*bg_node(1, ky) + bg_node(2, ky))/4
               bg_node(Cx, ky) = (7*bg_node(Cx - 1, ky) - 4*bg_node(Cx - 2, ky) + bg_node(Cx - 3, ky))/4
               sd_node(-1, ky) = (7*sd_node(0, ky) - 4*sd_node(1, ky) + sd_node(2, ky))/4
               sd_node(Cx, ky) = (7*sd_node(Cx - 1, ky) - 4*sd_node(Cx - 2, ky) + sd_node(Cx - 3, ky))/4
               fg_node(-1, ky) = (7*fg_node(0, ky) - 4*fg_node(1, ky) + fg_node(2, ky))/4
               fg_node(Cx, ky) = (7*fg_node(Cx - 1, ky) - 4*fg_node(Cx - 2, ky) + fg_node(Cx - 3, ky))/4
            end do
         else
            do ky = -1, Cy
               bg_node(-1, ky) = bg_node(0, ky)
               bg_node(Cx, ky) = bg_node(Cx - 1, ky)
               sd_node(-1, ky) = sd_node(0, ky)
               sd_node(Cx, ky) = sd_node(Cx - 1, ky)
               fg_node(-1, ky) = fg_node(0, ky)
               fg_node(Cx, ky) = fg_node(Cx - 1, ky)
            end do
         end if

         !*      construct logistic function using given background
         !*      and space-varying standard deviation s
         !*          1/h = 1 + exp[ (i-b)/s ] = 0 for i>>b and = 1 for i<<b
         !*      h = 0 - foreground. h = 1 - background
         do iy = 0, My - 1
            ky = int(iy/jy)                     !   coarse node block
            ay = iy/real(jy) - (ky + 0.5)         !   distance to bg node
            do ix = 0, Mx - 1
               kx = int(ix/jx)                 !   coarse node block
               ax = ix/real(jx) - (kx + 0.5)     !   distance to bg node

               gg(1:2) = (/ax*(1.5d0 - 2*ax*ax), ay*(1.5d0 - 2*ay*ay)/)

               !---    linear interpolate
               if (ax < 0) then
                  if (ay < 0) then
                                bb   =   bg_node(kx-1,ky-1)*gg(1)*gg(2) - bg_node(kx,ky-1)*(1+gg(1))*gg(2) - bg_node(kx-1,ky)*gg(1)*(1+gg(2)) + bg_node(kx,ky)*(1+gg(1))*(1+gg(2))
                                bstd =   sd_node(kx-1,ky-1)*gg(1)*gg(2) - sd_node(kx,ky-1)*(1+gg(1))*gg(2) - sd_node(kx-1,ky)*gg(1)*(1+gg(2)) + sd_node(kx,ky)*(1+gg(1))*(1+gg(2))
                                ff   =   fg_node(kx-1,ky-1)*gg(1)*gg(2) - fg_node(kx,ky-1)*(1+gg(1))*gg(2) - fg_node(kx-1,ky)*gg(1)*(1+gg(2)) + fg_node(kx,ky)*(1+gg(1))*(1+gg(2))
                  else
                                bb   = - bg_node(kx-1,ky+1)*gg(1)*gg(2) + bg_node(kx,ky+1)*(1+gg(1))*gg(2) - bg_node(kx-1,ky)*gg(1)*(1-gg(2)) + bg_node(kx,ky)*(1+gg(1))*(1-gg(2))
                                bstd = - sd_node(kx-1,ky+1)*gg(1)*gg(2) + sd_node(kx,ky+1)*(1+gg(1))*gg(2) - sd_node(kx-1,ky)*gg(1)*(1-gg(2)) + sd_node(kx,ky)*(1+gg(1))*(1-gg(2))
                                ff   = - fg_node(kx-1,ky+1)*gg(1)*gg(2) + fg_node(kx,ky+1)*(1+gg(1))*gg(2) - fg_node(kx-1,ky)*gg(1)*(1-gg(2)) + fg_node(kx,ky)*(1+gg(1))*(1-gg(2))
                  end if
               else
                  if (ay < 0) then
                                bb   = - bg_node(kx+1,ky-1)*gg(1)*gg(2) - bg_node(kx,ky-1)*(1-gg(1))*gg(2) + bg_node(kx+1,ky)*gg(1)*(1+gg(2)) + bg_node(kx,ky)*(1-gg(1))*(1+gg(2))
                                bstd = - sd_node(kx+1,ky-1)*gg(1)*gg(2) - sd_node(kx,ky-1)*(1-gg(1))*gg(2) + sd_node(kx+1,ky)*gg(1)*(1+gg(2)) + sd_node(kx,ky)*(1-gg(1))*(1+gg(2))
                                ff   = - fg_node(kx+1,ky-1)*gg(1)*gg(2) - fg_node(kx,ky-1)*(1-gg(1))*gg(2) + fg_node(kx+1,ky)*gg(1)*(1+gg(2)) + fg_node(kx,ky)*(1-gg(1))*(1+gg(2))
                  else
                                bb   =   bg_node(kx+1,ky+1)*gg(1)*gg(2) + bg_node(kx,ky+1)*(1-gg(1))*gg(2) + bg_node(kx+1,ky)*gg(1)*(1-gg(2)) + bg_node(kx,ky)*(1-gg(1))*(1-gg(2))
                                bstd =   sd_node(kx+1,ky+1)*gg(1)*gg(2) + sd_node(kx,ky+1)*(1-gg(1))*gg(2) + sd_node(kx+1,ky)*gg(1)*(1-gg(2)) + sd_node(kx,ky)*(1-gg(1))*(1-gg(2))
                                ff   =   fg_node(kx+1,ky+1)*gg(1)*gg(2) + fg_node(kx,ky+1)*(1-gg(1))*gg(2) + fg_node(kx+1,ky)*gg(1)*(1-gg(2)) + fg_node(kx,ky)*(1-gg(1))*(1-gg(2))
                  end if
               end if

               if (bstd > 0) then

                  snr = (img_fore(ix, iy) - bb)/bstd
                  if (snr > 22) then
                     logistic(ix, iy) = 0
                  else if (snr < -22) then
                     logistic(ix, iy) = 1
                  else
                     logistic(ix, iy) = 1/(1 + exp(snr))
                  end if
               else
                  if (img_fore(ix, iy) > bb) then
                     logistic(ix, iy) = 0
                  else
                     logistic(ix, iy) = 1
                  end if
               end if

            end do
         end do

         !---    create a stack of increasingly blurred images of the background

         blur_stack(0:Mx - 1, 0:My - 1, 0) = img_back(0:Mx - 1, 0:My - 1)
         ff = 0.0d0          !   previously used gaussian blur radius
         do nn = 1, nGB
            th = nn*sigma/nGB
            dsigma = sqrt(th*th - ff*ff)      !   new spreading sigma to use
            call weightedGaussianBlurImage(blur_stack(:, :, nn - 1), logistic, dsigma, blur_stack(:, :, nn))
            ff = th
         end do

         !---    for each pixel, find the weighted average derivative, reverse order, construct weighted value for background

         do iy = 0, My - 1
            do ix = 0, Mx - 1

               !---    find the lengthscale weighted average background level

               ww = abs(3*blur_stack(ix, iy, 0) - 4*blur_stack(ix, iy, 1) + blur_stack(ix, iy, 2))**power
               ff = ww*blur_stack(ix, iy, 0)
               wsum = ww
               do nn = 1, nGB - 1
                  ww = abs(blur_stack(ix, iy, nn + 1) - blur_stack(ix, iy, nn - 1))**power
                  wsum = wsum + ww
                  ff = ff + ww*blur_stack(ix, iy, nn)
               end do
               ww = abs(3*blur_stack(ix, iy, nGB) - 4*blur_stack(ix, iy, nGB - 1) + blur_stack(ix, iy, nGB - 2))**power
               ff = ff + ww*blur_stack(ix, iy, nGB)
               wsum = wsum + ww
               if (wsum > 0) ff = ff/wsum

               !---    construct the background blur image
               img_back(ix, iy) = ff

            end do
         end do

         !---    find the average background level
         ff = sum(img_back(0:Mx - 1, 0:My - 1))/(Mx*My)
         print *, "loop ", loop, " average background level ", ff

         !---    construct a "foreground" image which gradually reduces the local variation due to background,
         !       and scales intensity from 0:1. This is used to improve the estimate of bg std dev and logistic fn
         !       in the next round.
         do iy = 0, My - 1
            do ix = 0, Mx - 1

           img_fore(ix, iy) = (1 - logistic(ix, iy)/2)*(img_work(ix, iy) - img_back(ix, iy)) + (logistic(ix, iy)/2)*img_fore(ix, iy)

            end do
         end do
         img_min = minval(img_fore(0:Mx - 1, 0:My - 1))
         ff = maxval(img_fore(0:Mx - 1, 0:My - 1)) - img_min
         img_irange = 1/max(1.0d-8, ff)
         img_fore(0:Mx - 1, 0:My - 1) = img_fore(0:Mx - 1, 0:My - 1)*img_irange - (img_min*img_irange)

         if (dbg) then
       call write_greyscale_png(trim(numberFile(trim(outfile), loop))//".back.png", img_back(0:Mx - 1, 0:My - 1), negative=negative)
       call write_greyscale_png(trim(numberFile(trim(outfile), loop))//".work.png", img_work(0:Mx - 1, 0:My - 1), negative=negative)
       call write_greyscale_png(trim(numberFile(trim(outfile), loop))//".fore.png", img_fore(0:Mx - 1, 0:My - 1), negative=negative)
         end if

         deallocate (bg_node)
         deallocate (fg_node)
         deallocate (sd_node)
         print *, ""
      end do      !   loops
      !---    at this point I have a small background image, logistic fn and foreground image

      !---    create full size background and logistic
      print *, "finished loops, reconstructing final logistic, back & fore images"
      if (half) then

         print *, ""
         if (quarter) then
            print *, "upscaling background image x 8"
            Mx = int((Nx + 1)/2); My = int((Ny + 1)/2)
            Mx = int((Mx + 1)/2); My = int((My + 1)/2)
            call doubleImage(Mx, My, logistic, img_work, pbc=.false.)
            Mx = int((Nx + 1)/2); My = int((Ny + 1)/2)
            call doubleImage(Mx, My, img_work, logistic, pbc=.false.)
            call doubleImage(Nx, Ny, logistic, img_work, pbc=.false.)
            logistic = img_work

            Mx = int((Nx + 1)/2); My = int((Ny + 1)/2)
            Mx = int((Mx + 1)/2); My = int((My + 1)/2)
            call doubleImage(Mx, My, img_back, img_work, pbc=.false.)
            Mx = int((Nx + 1)/2); My = int((Ny + 1)/2)
            call doubleImage(Mx, My, img_work, img_back, pbc=.false.)
            call doubleImage(Nx, Ny, img_back, img_work, pbc=.false.)
            img_back = img_work

            Mx = int((Nx + 1)/2); My = int((Ny + 1)/2)
            Mx = int((Mx + 1)/2); My = int((My + 1)/2)
            call doubleImage(Mx, My, img_fore, img_work, pbc=.false.)
            Mx = int((Nx + 1)/2); My = int((Ny + 1)/2)
            deallocate (img_fore); allocate (img_fore(0:Nx - 1, 0:Ny - 1))
            call doubleImage(Mx, My, img_work, img_fore, pbc=.false.)
            call doubleImage(Nx, Ny, img_fore, img_work, pbc=.false.)
            img_fore = img_work

         else
            print *, "upscaling background image x 2"
            call doubleImage(Nx, Ny, logistic, img_work, pbc=.false.)
            logistic = img_work
            call doubleImage(Nx, Ny, img_back, img_work, pbc=.false.)
            img_back = img_work
            call doubleImage(Nx, Ny, img_fore, img_work, pbc=.false.)
            deallocate (img_fore); allocate (img_fore(0:Nx - 1, 0:Ny - 1))
            img_fore = img_work

         end if

      else if (quarter) then

         print *, ""
         print *, "upscaling background image x 4"
         Mx = int((Nx + 1)/2); My = int((Ny + 1)/2)
         call doubleImage(Mx, My, logistic, img_work, pbc=.false.)
         call doubleImage(Nx, Ny, img_work, logistic, pbc=.false.)

         call doubleImage(Mx, My, img_back, img_work, pbc=.false.)
         call doubleImage(Nx, Ny, img_work, img_back, pbc=.false.)

         call doubleImage(Mx, My, img_fore, img_work, pbc=.false.)
         deallocate (img_fore); allocate (img_fore(0:Nx - 1, 0:Ny - 1))
         call doubleImage(Nx, Ny, img_work, img_fore, pbc=.false.)

      end if

      if (dbg) then
         call write_greyscale_png(trim(outfile)//".bg.png", img_back, negative=negative)
         call write_greyscale_png(trim(outfile)//".fg.png", img_fore, negative=negative)
         call write_greyscale_png(trim(outfile)//".lg.png", logistic)
      end if

      !---    final subtract background from foreground
      !       logistic = 0 = foreground
      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            img_out(ix, iy) = (1 - b0/2)*(1 - logistic(ix, iy))*(img(ix + 1, iy + 1) - img_back(ix, iy)) + b0/2
         end do
      end do
      img_min = minval(img_out(0:Nx - 1, 0:Ny - 1))
      ff = maxval(img_out(0:Nx - 1, 0:Ny - 1)) - img_min
      img_irange = 1/max(1.0d-8, ff)
      img_out(0:Nx - 1, 0:Ny - 1) = img_out(0:Nx - 1, 0:Ny - 1)*img_irange - (img_min*img_irange)

      if (retone) then

         !---    last pass, this time on the full size image
         call radialPowerFunction(img, q_min, q_max, rpf, rpfsum, qrpfbar, q2rpfbar, snr)
         sigma = sqrt(PI)/(2*qrpfbar)
         print *, "single Gaussian sigma    = ", sigma
         print *, ""
         vRC = ceiling(sigma*2)       !   std dev assumed to vary on scale 4 x characteristic scale of image
         if (vRC_in /= LIB_CLA_NODEFAULT_I) then
            print *, "variable Ridler-Calvard pixels estimated by FFTW ", vRC, " using input value ", vRC_in
            vRC = vRC_in
         end if

         !---    the flattening proceeds by computing the probability that each pixel is part of the background
         !       which requires an estimate of the background level and std dev for each pixel
         !       compute background and std dev of background on a course node grid
         Cx = max(1, Nx/vRC)
         Cy = max(1, Ny/vRC)

         allocate (bg_node(-1:Cx, -1:Cy))      !   defined -1:Cx for convenience- I'm going to extrapolate the edges.
         allocate (sd_node(-1:Cx, -1:Cy))      !
         allocate (fg_node(-1:Cx, -1:Cy))      !   b
         jx = ceiling(Nx*1.0/Cx)           !   number of pixels per block
         jy = ceiling(Ny*1.0/Cy)
print *, "variable Ridler-Calvard pixels in scaled image ", vRC, " grid nodes required ", Cx, ",", Cy, " px per block ", jx, ",", jy

         do ky = 0, Cy - 1
            do kx = 0, Cx - 1
      call findImageIntensityFeatures( img_out( kx*jx:min(kx*jx+jx-1,Nx) , ky*jy:min(ky*jy+jy-1,Ny) ) , bb,th,ff,bstd,fbar,snr,fot )
               bg_node(kx, ky) = bb
               sd_node(kx, ky) = bstd
               fg_node(kx, ky) = ff

            end do
         end do

         !---    finite difference extrapolation to borders.
         if (Cy > 3) then
            do kx = 0, Cx - 1
               bg_node(kx, -1) = (7*bg_node(kx, 0) - 4*bg_node(kx, 1) + bg_node(kx, 2))/4
               bg_node(kx, Cy) = (7*bg_node(kx, Cy - 1) - 4*bg_node(kx, Cy - 2) + bg_node(kx, Cy - 3))/4
               sd_node(kx, -1) = (7*sd_node(kx, 0) - 4*sd_node(kx, 1) + sd_node(kx, 2))/4
               sd_node(kx, Cy) = (7*sd_node(kx, Cy - 1) - 4*sd_node(kx, Cy - 2) + sd_node(kx, Cy - 3))/4
               fg_node(kx, -1) = (7*fg_node(kx, 0) - 4*fg_node(kx, 1) + fg_node(kx, 2))/4
               fg_node(kx, Cy) = (7*fg_node(kx, Cy - 1) - 4*fg_node(kx, Cy - 2) + fg_node(kx, Cy - 3))/4
            end do
         else
            do kx = 0, Cx - 1
               bg_node(kx, -1) = bg_node(kx, 0)
               bg_node(kx, Cy) = bg_node(kx, Cy - 1)
               sd_node(kx, -1) = sd_node(kx, 0)
               sd_node(kx, Cy) = sd_node(kx, Cy - 1)
               fg_node(kx, -1) = fg_node(kx, 0)
               fg_node(kx, Cy) = fg_node(kx, Cy - 1)
            end do
         end if
         if (Cx > 3) then
            do ky = -1, Cy
               bg_node(-1, ky) = (7*bg_node(0, ky) - 4*bg_node(1, ky) + bg_node(2, ky))/4
               bg_node(Cx, ky) = (7*bg_node(Cx - 1, ky) - 4*bg_node(Cx - 2, ky) + bg_node(Cx - 3, ky))/4
               sd_node(-1, ky) = (7*sd_node(0, ky) - 4*sd_node(1, ky) + sd_node(2, ky))/4
               sd_node(Cx, ky) = (7*sd_node(Cx - 1, ky) - 4*sd_node(Cx - 2, ky) + sd_node(Cx - 3, ky))/4
               fg_node(-1, ky) = (7*fg_node(0, ky) - 4*fg_node(1, ky) + fg_node(2, ky))/4
               fg_node(Cx, ky) = (7*fg_node(Cx - 1, ky) - 4*fg_node(Cx - 2, ky) + fg_node(Cx - 3, ky))/4
            end do
         else
            do ky = -1, Cy
               bg_node(-1, ky) = bg_node(0, ky)
               bg_node(Cx, ky) = bg_node(Cx - 1, ky)
               sd_node(-1, ky) = sd_node(0, ky)
               sd_node(Cx, ky) = sd_node(Cx - 1, ky)
               fg_node(-1, ky) = fg_node(0, ky)
               fg_node(Cx, ky) = fg_node(Cx - 1, ky)
            end do
         end if

         do iy = 0, Ny - 1
            ky = int(iy/jy)                     !   coarse node block
            ay = iy/real(jy) - (ky + 0.5)         !   distance to bg node
            do ix = 0, Nx - 1
               kx = int(ix/jx)                 !   coarse node block
               ax = ix/real(jx) - (kx + 0.5)     !   distance to bg node
               !---    smooth interpolate
               gg(1:2) = (/ax*(1.5d0 - 2*ax*ax), ay*(1.5d0 - 2*ay*ay)/)
               if (ax < 0) then
                  if (ay < 0) then
                                bb   =   bg_node(kx-1,ky-1)*gg(1)*gg(2) - bg_node(kx,ky-1)*(1+gg(1))*gg(2) - bg_node(kx-1,ky)*gg(1)*(1+gg(2)) + bg_node(kx,ky)*(1+gg(1))*(1+gg(2))
                                bstd =   sd_node(kx-1,ky-1)*gg(1)*gg(2) - sd_node(kx,ky-1)*(1+gg(1))*gg(2) - sd_node(kx-1,ky)*gg(1)*(1+gg(2)) + sd_node(kx,ky)*(1+gg(1))*(1+gg(2))
                                ff   =   fg_node(kx-1,ky-1)*gg(1)*gg(2) - fg_node(kx,ky-1)*(1+gg(1))*gg(2) - fg_node(kx-1,ky)*gg(1)*(1+gg(2)) + fg_node(kx,ky)*(1+gg(1))*(1+gg(2))
                  else
                                bb   = - bg_node(kx-1,ky+1)*gg(1)*gg(2) + bg_node(kx,ky+1)*(1+gg(1))*gg(2) - bg_node(kx-1,ky)*gg(1)*(1-gg(2)) + bg_node(kx,ky)*(1+gg(1))*(1-gg(2))
                                bstd = - sd_node(kx-1,ky+1)*gg(1)*gg(2) + sd_node(kx,ky+1)*(1+gg(1))*gg(2) - sd_node(kx-1,ky)*gg(1)*(1-gg(2)) + sd_node(kx,ky)*(1+gg(1))*(1-gg(2))
                                ff   = - fg_node(kx-1,ky+1)*gg(1)*gg(2) + fg_node(kx,ky+1)*(1+gg(1))*gg(2) - fg_node(kx-1,ky)*gg(1)*(1-gg(2)) + fg_node(kx,ky)*(1+gg(1))*(1-gg(2))
                  end if
               else
                  if (ay < 0) then
                                bb   = - bg_node(kx+1,ky-1)*gg(1)*gg(2) - bg_node(kx,ky-1)*(1-gg(1))*gg(2) + bg_node(kx+1,ky)*gg(1)*(1+gg(2)) + bg_node(kx,ky)*(1-gg(1))*(1+gg(2))
                                bstd = - sd_node(kx+1,ky-1)*gg(1)*gg(2) - sd_node(kx,ky-1)*(1-gg(1))*gg(2) + sd_node(kx+1,ky)*gg(1)*(1+gg(2)) + sd_node(kx,ky)*(1-gg(1))*(1+gg(2))
                                ff   = - fg_node(kx+1,ky-1)*gg(1)*gg(2) - fg_node(kx,ky-1)*(1-gg(1))*gg(2) + fg_node(kx+1,ky)*gg(1)*(1+gg(2)) + fg_node(kx,ky)*(1-gg(1))*(1+gg(2))
                  else
                                bb   =   bg_node(kx+1,ky+1)*gg(1)*gg(2) + bg_node(kx,ky+1)*(1-gg(1))*gg(2) + bg_node(kx+1,ky)*gg(1)*(1-gg(2)) + bg_node(kx,ky)*(1-gg(1))*(1-gg(2))
                                bstd =   sd_node(kx+1,ky+1)*gg(1)*gg(2) + sd_node(kx,ky+1)*(1-gg(1))*gg(2) + sd_node(kx+1,ky)*gg(1)*(1-gg(2)) + sd_node(kx,ky)*(1-gg(1))*(1-gg(2))
                                ff   =   fg_node(kx+1,ky+1)*gg(1)*gg(2) + fg_node(kx,ky+1)*(1-gg(1))*gg(2) + fg_node(kx+1,ky)*gg(1)*(1-gg(2)) + fg_node(kx,ky)*(1-gg(1))*(1-gg(2))
                  end if
               end if

               img_fore(ix, iy) = ff
               img_back(ix, iy) = bb

               img_out(ix, iy) = intensityShift(img_out(ix, iy), bb, b0, ff, f0)
            end do
         end do

      end if  !   retone

   else        !   don't flat

      if (retone) then

         call findImageIntensityFeatures(img, bb, th, ff, bstd, fbar, snr, fot)
         do iy = 0, Ny - 1
            do ix = 0, Nx - 1
               img_out(ix, iy) = intensityShift(img(ix + 1, iy + 1), bb, b0, ff, f0)
            end do
         end do

      else    !   do nothing.

         img_out = img

      end if

   end if      !   flat

   !---    output result
   call findImageIntensityFeatures(img_out, bb, th, ff, bstd, fbar, snr, fot)
   print *, "   final"
   print *, "       : back, thresh, fore     ", bb, th, ff
   print *, "         back std, <f>,snr,fot  ", bstd, fbar, snr, fot

   print *, "writing file"
   call write_greyscale_png(outfile, img_out, negative=negative)

   print *, ""
   print *, "done"
   print *, ""

contains
!---^^^^^^^^

   subroutine gaussianBlurImage(f_in, t, f_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      produce a Gaussian blur of input image, width t
      real(kind=real64), dimension(0:, 0:), intent(in)       ::      f_in    !Image array in
      real(kind=real64), dimension(0:, 0:), intent(inout)    ::      f_out   !Blurred image array out
      real(kind=real64), intent(in)                        ::      t       !Gaussian width to blur

      real(kind=real64), dimension(:), allocatable      ::      kernel, f_stripe, w_stripe   !Arrays to hold kernel and the image broken into strips
      !in the x and y direction (for parallelization).

      integer             ::      Nx, Ny               !Image size in X and Y dimension/pixels
      integer             ::      ix, iy, jj, kk         !Pixel X, pixel Y, strip and kernel indicies.
      integer             ::      Nk                  !number of pixels in kernel
      real(kind=real64)   ::      i2s2, ww, wf, ws       !0.5*t^2, kernal weight, kernel result, strip weight

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

   subroutine weightedGaussianBlurImage(f_in, weight, t, f_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      produce a Gaussian blur of input image, width t
      real(kind=real64), dimension(0:, 0:), intent(in)       ::      f_in            !Image array in
      real(kind=real64), dimension(0:, 0:), intent(in)       ::      weight          !Array of weights
      real(kind=real64), dimension(0:, 0:), intent(inout)    ::      f_out           !Blurred image out
      real(kind=real64), intent(in)                        ::      t               !Gaussian width to blur

      real(kind=real64), dimension(:), allocatable      ::      kernel, f_stripe    !arrays holding kernel and a strip of the image
      real(kind=real64), dimension(:, :), allocatable    ::      kernelweight        !Array of kernel weights

      integer             ::      Nx, Ny               !Image x and y dimensions/pixles
      integer             ::      ix, iy, jj, kk         !Pixel x, pixel y, strip and kernel indices.
      integer             ::      Nk                  !number of pixels in kernel
      real(kind=real64)   ::      i2s2, ww, wf          !0.5*t^2, kernel weight, kernel result

      Nx = size(f_in, dim=1)
      Ny = size(f_in, dim=2)

      !---    compute the unnormalised kernel
      Nk = max(5, ceiling(t*3))             !   range of pixels to search is +/- Nk
      allocate (kernel(0:Nk))
      i2s2 = 1/(2*t*t)
      kernel(0) = 1.0d0
      do ix = 1, Nk
         kernel(ix) = exp(-ix*ix*i2s2)
      end do

      !---    compute the output blurred image. First do x strips.
      allocate (f_stripe(0:Ny - 1))
      allocate (kernelweight(0:Nx - 1, 0:Ny - 1))

!$OMP PARALLEL  PRIVATE(ix,iy,wf,jj,kk,ww,f_stripe  ) SHARED( f_in,weight,f_out,Nx,Ny,Nk,kernel ,kernelweight)

!   first compute the unnormalised gaussian function
      !---    x strips
!$OMP DO
      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            wf = 0.0d0
            do jj = max(0, ix - Nk), min(Nx - 1, ix + Nk)
               kk = abs(jj - ix)
               ww = kernel(kk)
               wf = wf + ww*f_in(jj, iy)*weight(jj, iy)
            end do
            f_out(ix, iy) = wf
         end do
      end do
!$OMP END DO

      !---    now do y strips
!$OMP DO
      do ix = 0, Nx - 1
         !   make a copy of this stripe
         f_stripe(0:Ny - 1) = f_out(ix, 0:Ny - 1)
         do iy = 0, Ny - 1
            wf = 0.0d0
            do jj = max(0, iy - Nk), min(Ny - 1, iy + Nk)
               kk = abs(jj - iy)
               ww = kernel(kk)
               wf = wf + ww*f_stripe(jj)
            end do
            f_out(ix, iy) = wf

         end do
      end do
!$OMP END DO

!   now have an unnormalised blur function. Repeat, this time computing the weighting
      !---    x strips
!$OMP DO
      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            wf = 0.0d0
            do jj = max(0, ix - Nk), min(Nx - 1, ix + Nk)
               kk = abs(jj - ix)
               wf = wf + kernel(kk)*weight(jj, iy)
            end do
            kernelweight(ix, iy) = wf
         end do
      end do
!$OMP END DO

      !---    now do y strips
!$OMP DO
      do ix = 0, Nx - 1
         !   make a copy of this stripe
         f_stripe(0:Ny - 1) = kernelweight(ix, 0:Ny - 1)
         do iy = 0, Ny - 1
            wf = 0.0d0
            do jj = max(0, iy - Nk), min(Ny - 1, iy + Nk)
               kk = abs(jj - iy)
               wf = wf + kernel(kk)*f_stripe(jj)
            end do
            wf = max(1.0d-16, wf)
            f_out(ix, iy) = f_out(ix, iy)/wf

         end do
      end do
!$OMP END DO

!$OMP END PARALLEL
      !---    now have a normalised Gaussian blur function

      return
   end subroutine weightedGaussianBlurImage

   elemental function intensityShift(g, b, b0, f, f0) result(g_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !      then use affine shifts to push 0->0, b->b0, t->t0, f->f0, 1->1
      !       return the new intensity level g
      !
      real(kind=real64), intent(in)            ::      g           !input pixel intensity
      real(kind=real64), intent(in)            ::      b, b0, f, f0   !background, post filter background, foreground and post filter foreground intensities
      real(kind=real64)                       ::      g_out       !output pixel intensity

      real(kind=real64)               ::      tt, tt0              !pre and post filter Threshold inensities,respectively.

      tt = (b + f)/2
      tt0 = (b0 + f0)/2

      if (g <= 0.0) then
         g_out = 0.0d0
      else if (g < b) then
         g_out = b0*g/b
      else if (g < f) then
         g_out = b0 + (g - b)*(f0 - b0)/(f - b)
      else if (g < 1) then
         g_out = f0 + (g - f)*(1.0d0 - f0)/(1.0d0 - f)
      else
         g_out = 1.0d0
      end if

      return
   end function intensityShift

end program flattenBackground

