
program countSpots
!---^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    countSpots from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
   use Lib_CommandLineArguments
   use Lib_Gaussian2d
   use Lib_MultipleGaussians2d
   use Lib_Maxima2d
   use Lib_Filenames
   use AnalyseSpots
   use Lib_ColourScale
   use Lib_Png
   use Lib_Quicksort
   use Lib_RandomSeed                          !   which might seem odd, but its so that when random colours are ascribed to output maps, they are always the same.
   use Lib_Kernels, getKernel => get
   use iso_fortran_env
   implicit none

   !---

   character(len=*), parameter      ::      VERSION = "4.3.0"

   !---        input options

   type(CommandLineArguments)      ::      cla                     !Object to hold commmand line argumments (CLA)

   character(len=256)              ::      filename = "", outfile   !input and output filename respectively
   real(kind=real64)               ::      nmperpixel = 1.0d0      !length per pixel/nm
   real(kind=real64)               ::      diamscale = 2.0d0      !diameter = nmPerPixel*(2 sigma)*diamscale
   logical                         ::      negative = .false.      !Should the image intnesites be inverted?
   logical                         ::      quiet = .false.         !Should we output .count and .spots files?

!     !---        analysis options
!         real(kind=real64)               ::      roirmax0 = 0
!         real(kind=real64)               ::      maxpix0 = 0
!         real(kind=real64)               ::      f_thresh0 = 0

   !---        output options
   logical                         ::      oppng = .false.         !Output pnf of spot locations, set using CLA

   logical                         ::      avgdiam = .false.       !Report average spot diameter in histogram, set using CLA
   integer                         ::      nbins = 20              !number of histogram bins, set using CLA
   real(kind=real64)               ::      x = 50.0d0             !value of max histogram bin, set using CLA

   !---        solution
   real(kind=real64), dimension(:, :), allocatable        ::      img !array to hold image
   integer                         ::      Nx, Ny
   type(MultipleGaussians2d)                       ::      fitted_g2d  !object to hold parameters of several gaussians.
   integer                                         ::      nTotal      !muber of detected spots
   integer, dimension(:), allocatable                ::      hist        !Array to hold histogram

   real(kind=real64)               ::      bg_sigma                    !background intensity standard deviation

   !---        dummy variables
   integer                         ::      ii, jj   !indices of detected spots and histogram bins respectively
   real(kind=real64)               ::      dd, ff, tt, aa !see below
   logical                         ::      ok      !does input file exist?
   real(kind=real64), dimension(2)  ::      work    !array to hold debug values
   real(kind=real64)               ::      test_bin, bin_width
   real(kind=real64), dimension(:), allocatable      ::      spotDiams
   character(len=256)              ::      dummy   ! string to store name of colour scale

   !real dummies
   !dd: average gaussian diameter
   !ff, fvar, fsig: gaussian intensity, variance of ff, std.dev of ff respectively, note ff is resused to hold mean(ff) later.
   !tt, tvar, tsig: t-value (t-statisitics), variance of tt, std.dev of tt respectively, note tt is resused to hold mean(tt) later.
   !aa: micrograph area excluding dead pixels
   !---    command line options

   !---    read command line arguments
   cla = CommandLineArguments_ctor(30)

   call setProgramDescription(cla, "count.exe" &
        //"\n    Count detects bright spot features (such as an inverted brightfield image of black spot radiation damage) and  \n    fits 2D gaussians to the spots and produces a size frequency histogram." )
   call setProgramVersion(cla, VERSION)

   call setCategories(cla, (/"file handling          ", &
                             "input image parameters ", &
                             "detection parameters   ", &
                             "output options         ", &
                             "debug / misc           "/))

   call get(cla, "f", filename, LIB_CLA_REQUIRED, "            input filename", 1)
   outfile = removeSuffix(filename)
   call get(cla, "o", outfile, LIB_CLA_OPTIONAL, "          output filename", 1)
   call get(cla, "png", oppng, LIB_CLA_OPTIONAL, "        output .png of spot locations", 1)
   call get(cla, "roi", oproi, LIB_CLA_OPTIONAL, "        output roi map", 1)
   call get(cla, "randc", oprandc, LIB_CLA_OPTIONAL, "      output Ridler&Calvard threshold map", 1)

   call get(cla, "s", nmperpixel, LIB_CLA_OPTIONAL, "          image resolution (nm per px)", 2)
   call get(cla, "negative", negative, LIB_CLA_OPTIONAL, "   input negative image", 2)
   call get(cla, "quiet", quiet, LIB_CLA_OPTIONAL, "   quiet mode ( don't output .spots and .count files )", 2)

   call get(cla, "t", T_THRESH, LIB_CLA_OPTIONAL, "          t* threshold", 3)
   call get(cla, "d", D_THRESH, LIB_CLA_OPTIONAL, "          intensity dip between separable maxima", 3)
   call get(cla, "e", E_THRESH, LIB_CLA_OPTIONAL, "          max eccentricity", 3)
   call get(cla, "b", F_THRESH, LIB_CLA_OPTIONAL, "          intensity threshhold (f0/sigma) ( 0 for automatic )", 3)
   call get(cla, "tond", TOND_THRESH, LIB_CLA_OPTIONAL, "       t/d threshold", 3)
   call get(cla, "q", Q_THRESH, LIB_CLA_OPTIONAL, "          q threshold", 3)
   call get(cla, "m", MAXPIX, LIB_CLA_OPTIONAL, "          maximum fraction of foreground pixels - ( 0 for automatic )", 3)
   call get(cla, "X", ROIXMAX, LIB_CLA_OPTIONAL, "          maximum roi dimension (px)", 3)
   call get(cla, "R", ROIRMAX, LIB_CLA_OPTIONAL, "          maximum roi padding (px)", 3)

   call get(cla, "x", x, LIB_CLA_OPTIONAL, "          max histogram bin size", 4)
   call get(cla, "nbins", nbins, LIB_CLA_OPTIONAL, "      number of histogram bins", 4)
   call get(cla, "scale", diamscale, LIB_CLA_OPTIONAL, "      diameter scaling", 4)
   call get(cla, "a", avgdiam, LIB_CLA_OPTIONAL, "          report avg diameter in histogram", 4)

   ii = 2; call get(cla, "dbg", work, ii, LIB_CLA_OPTIONAL, "        debug ROI in vicinty of point (x,y)", 5)
   if (hasArgument(cla, "dbg")) then
      DBG_ROIX = nint(work(1))
      DBG_ROIY = nint(work(2))
   end if
   call get(cla, "count", FITANDCOUNT, LIB_CLA_OPTIONAL, "    count spots", 5)
   dummy = "IBM"; call get(cla, "colour", dummy, LIB_CLA_OPTIONAL, "     colourscale ", 4)
   if (hasArgument(cla, "colour")) then
      COLOURSCALE_DEFAULT = getColourScale(dummy)
      if (COLOURSCALE_DEFAULT == COLOURSCALE_UNSET) then
         print *, "error: did not recognise colourscale. Options are"
         call listAvailableColourScales()
         stop
      end if

   end if

   if (hasArgument(cla, "s")) LENGTHUNIT = "nm"              !   otherwise output result in px

   call init_random_seed(12345)

   call report(cla)
   if (.not. allRequiredArgumentsSet(cla)) stop
   if (hasHelpArgument(cla)) stop
   call delete(cla)

   !---

   print *, "count v."//VERSION
   print *, "^^^^^^^^"//repeat("^", len(VERSION))
   print *, ""

   !---    check that the file exists
   inquire (file=trim(filename), exist=ok)
   if (.not. ok) then
      print *, "error: couldn't find file """//trim(filename)//""""
      stop
   end if

   !---    welcome
   print *, "   filename                """//trim(filename)//""""
   print *, "   outfile                 """//trim(outfile)//""""
   print *, "   scale nm/px             ", nmperpixel
   print *, "   diam scale              ", diamscale
   print *, "   average diameter        ", avgdiam
   print *, "   debug roi?              ", DBG_ROIX, ",", DBG_ROIY
   print *, "   t* threshold            ", T_THRESH
   print *, "   d threshold             ", D_THRESH
   print *, "   e threshold             ", E_THRESH
   print *, "   t*/d threshold          ", TOND_THRESH
   if (F_THRESH == 0) then
      print *, "   brightness threshold    auto"
   else
      print *, "   brightness threshold    ", F_THRESH
   end if
   if (MAXPIX == 0) then
      print *, "   maximum fg pix frac     auto"
   else
      print *, "   maximum fg pix frac     ", MAXPIX
   end if
   print *, "   maximum roi dimension   ", ROIXMAX
   if (ROIXMAX == 0) then
      print *, "   maximum roi padding     auto"
   else
      print *, "   maximum roi padding     ", ROIRMAX
   end if
   if (nbins > 0) then
      print *, "   histogram bins          ", nbins

   end if
   if (hasArgument(cla, "x")) then

      if (x > 0) then
         print *, "   histogram max val defined by user:", x!user defined max bin size
      elseif (x == 0) then
         print *, "   Histogram max val calculated from spot data"!get a reasonalbe value for max bin size from the data
      end if
   else
      x = 50.0d0 !default
      print *, "   histogram max val set to default:", x!
   end if

   print *, "   output png              ", oppng
   print *, "   default colourscale      """//trim(getColourScaleName(COLOURSCALE_DEFAULT))//""""
   print *, "   negative                ", negative
   print *, "   count spots?            ", FITANDCOUNT
   print *, "   output .count and .spots",.not. quiet
   print *, ""

   !---    job starts here
   print *, ""
   print *, "reading image file"
   print *, "^^^^^^^^^^^^^^^^^^"

   call readPng(filename, img, negative)
   nx = size(img, dim=1)
   ny = size(img, dim=2)

   print *, ""

   if (oproi .or. oprandc) then
      call fitSpots(img,fitted_g2d,f_auto=(F_THRESH==0),maxpix_auto=(MAXPIX==0),roi_auto=(ROIRMAX==0),ss=bg_sigma,ofile=outfile)
   else
      call fitSpots(img, fitted_g2d, f_auto=(F_THRESH == 0), maxpix_auto=(MAXPIX == 0), roi_auto=(ROIRMAX == 0), ss=bg_sigma)
   end if

   if (FITANDCOUNT) then
      nTotal = getN(fitted_g2d)

      print *, ""
      print *, "output spots"
      print *, "^^^^^^^^^^^^"

      if (x == 0) then
         if (hasArgument(cla, "nbins")) then
            print *, "Bin width is calculated automatically (-x=0), -nbins user input will be overridden"
         end if
         allocate (spotDiams(nTotal))
         spotDiams = 0.0d0
         call getScaledDiameters(fitted_g2d, spotDiams)
         x = automatic_xrange(spotDiams)
         bin_width = automatic_binwidth(spotDiams)
         call automatic_nBins(x, bin_width, nBins)
         print *, "Optimising histogram:", nBins, "bins of width", bin_width, "max bin is", x
      end if

      if (.not. quiet) then
                call outputSpots( trim(outfile)//".spots",nTotal,fitted_g2d, size(img,dim=1),size(img,dim=2),nmperpixel, count(img==0),"",diamscale,bg_sigma)

         open (unit=701, file=trim(outfile)//".count", action="write")
         call outputCount(701, nmperpixel, diamscale)
         close (unit=701)
      end if
      call outputCount(6, nmperpixel, diamscale)

      if (oppng) call outputPng(trim(outfile)//".spots.png", fitted_g2d, img, diamScale, x, negative)
   end if

   print *, ""
   print *, "done"
   print *, ""

contains
!---^^^^^^^^

   subroutine outputCount(u, nmperpixel, diamscale)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      integer, intent(in)              ::      u                   !Fortran file unit/ID
      real(kind=real64), intent(in)    ::      nmperpixel          ! length per pixel/nm
      real(kind=real64), intent(in)    ::      diamscale           !   diameter reported as nmperpixel * (2 * sigma) * diamscale

      real(kind=real64)       ::      dminor_bar, dmajor_bar, davg_bar  !mean minor, major and average gaussian diameters
      !Note the above vars are used to sum their respective comnents first.

      real(kind=real64)               ::      ff, tt, qq                !   intensity,t* and quality for individual spot

      real(kind=real64), dimension(:), allocatable      ::      fbar, f2bar !sum of intensity values(f)
      real(kind=real64), dimension(:), allocatable      ::      tbar, t2bar !sum of t values (t)
      real(kind=real64), dimension(:), allocatable      ::      qbar, q2bar !sum of quality values (q)

      real(kind=real64)               ::     fvar, tvar, qvar       !   variance of intensity,t*,quality
      real(kind=real64)               ::     fsig, tsig, qsig       !   std.dev. of intensity,t*,quality
      real(kind=real64)               ::     hpera, minhpera, maxhpera
      integer                         ::      nn      !   number in bin

      if (nbins <= 0) return

      !---    construct the histograms
      allocate (hist(0:nbins))
      allocate (fbar(0:nbins))
      allocate (f2bar(0:nbins))
      allocate (tbar(0:nbins))
      allocate (t2bar(0:nbins))
      allocate (qbar(0:nbins))
      allocate (q2bar(0:nbins))
      hist = 0
      fbar = 0.0d0; f2bar = 0.0d0
      tbar = 0.0d0; t2bar = 0.0d0
      qbar = 0.0d0; q2bar = 0.0d0
      dminor_bar = 0.0d0
      dmajor_bar = 0.0d0
      davg_bar = 0.0d0
      do ii = 1, nTotal
         dd = getSigma(fitted_g2d, ii)
         davg_bar = davg_bar + dd
         if (avgdiam) then
            dminor_bar = dminor_bar + getSigma(fitted_g2d, ii, larger=.false.)
            dmajor_bar = dmajor_bar + getSigma(fitted_g2d, ii, larger=.true.)
         else
            dminor_bar = dminor_bar + getSigma(fitted_g2d, ii, larger=.false.)
            dd = getSigma(fitted_g2d, ii, larger=.true.)
            dmajor_bar = dmajor_bar + dd
         end if
         dd = dd*nmperpixel*(2*diamscale)       !   2 for convert radius to diameter.

         ff = getF(fitted_g2d, ii)
         tt = getT(fitted_g2d, ii)
         qq = getQ(fitted_g2d, ii)
         jj = int(dd*nbins/x)
         jj = max(0, min(nbins, jj))
         hist(jj) = hist(jj) + 1
         fbar(jj) = fbar(jj) + ff
         f2bar(jj) = f2bar(jj) + ff*ff
         tbar(jj) = tbar(jj) + tt
         t2bar(jj) = t2bar(jj) + tt*tt
         qbar(jj) = qbar(jj) + qq
         q2bar(jj) = q2bar(jj) + qq*qq
      end do
      davg_bar = davg_bar*nmperpixel*(2*diamscale)/nTotal               !   2 converts radius into diameter
      dminor_bar = dminor_bar*nmperpixel*(2*diamscale)/nTotal
      dmajor_bar = dmajor_bar*nmperpixel*(2*diamscale)/nTotal

      !   output histogram
      write (unit=u, fmt='(a)') "# count v."//VERSION//" automatic sizing count"
      write (unit=u, fmt='(a)') "# "//trim(filename)
 write(unit=u,fmt='(a,2f12.3,a)') "# ",size(img,dim=1)*nmperpixel,size(img,dim=2)*nmperpixel," micrograph extent ("//LENGTHUNIT//")"

      aa = count(img > 0)*nmperpixel*nmperpixel
      write (unit=u, fmt='(a,f16.3,a)') "# ", aa, " micrograph area excluding dead pixels ("//LENGTHUNIT//"^2)"
      write (unit=u, fmt='(a,f12.3,a)') "# ", D_THRESH, " distinct spot threshold"
      write (unit=u, fmt='(a,f12.3,a)') "# ", E_THRESH, " spot maximum axis ratio"
      write (unit=u, fmt='(a,f12.3,a)') "# ", ROIRMAX, " roi padding radius"

      write (unit=u, fmt='(a,f12.3,a)') "# ", T_THRESH, " t-test threshold"
      !write(unit=u,fmt='(a,f12.3,a)') "# ",Q_THRESH," quality threshold" !unused currently
      write (unit=u, fmt='(a,f12.3,a)') "# ", DETECTFMIN, " bright pixel threshold"
      write (unit=u, fmt='(a,f12.3,a)') "# ", MAXPIX, " max pixels for threshold"
      write (unit=u, fmt='(a,f12.3,a)') "# ", diamscale, " diameter scaling: d = (2 sigma)*scale"

      write (unit=u, fmt='(a,i12,a)') "# ", nTotal, " total spot count"

      write (unit=u, fmt='(a,f12.3,a)') "# ", dmajor_bar, " average major diameter ("//LENGTHUNIT//")"
      write (unit=u, fmt='(a,f12.3,a)') "# ", dminor_bar, " average minor diameter ("//LENGTHUNIT//")"
      write (unit=u, fmt='(a,f12.3,a)') "# ", davg_bar, " (geometric) average diameter ("//LENGTHUNIT//")"
      write (unit=u, fmt='(a,i12,f12.3,a)') "# ", nTotal, davg_bar, " headline result (count,avg diam)"

      ff = sum(fbar)/max(nTotal, 1)
      fvar = sum(f2bar)/max(nTotal, 1)
      tt = sum(tbar)/max(nTotal, 1)
      tvar = sum(t2bar)/max(nTotal, 1)
      qq = sum(qbar)/max(nTotal, 1)
      qvar = sum(q2bar)/max(nTotal, 1)
      fsig = 0.0d0; tsig = 0.0d0; qsig = 0.0d0
      if (nTotal > 1) then
         fsig = sqrt(fvar - ff*ff)*sqrt(nTotal/real((nTotal + 1), kind=real64))           !   note: previous version used nTotal+1.5, which is a better unbiassing of variance
         tsig = sqrt(tvar - tt*tt)*sqrt(nTotal/real((nTotal + 1), kind=real64))           !   but everyone uses nTotal+1, so we probably should too :)
         qsig = sqrt(qvar - qq*qq)*sqrt(nTotal/real((nTotal + 1), kind=real64))
      end if
      if (nTotal > 1) fsig = sqrt(fvar - ff*ff)*sqrt((nTotal + 1)/real(nTotal, kind=real64))
      write (unit=u, fmt='(3(a,f12.3))') "# ", ff, " +/-", fsig, " global mean intens, sample std dev"
      write (unit=u, fmt='(3(a,f12.3))') "# ", tt, " +/-", tsig, " global mean t*, sample std dev"
      write (unit=u, fmt='(3(a,f12.3))') "# ", qq, " +/-", qsig, " global mean AIC quality, sample std dev"
      write (unit=u, fmt='(a,i12,a)') "# ", nbins, " number of bins"
      write (unit=u, fmt='(a,f12.3,a)') "# ", x, " maximum diameter bin ("//LENGTHUNIT//")"
                write(unit=u,fmt='(a5,a9,a8,3a16,100a14)') "# bin","diam("//LENGTHUNIT//")"," count","count min/"//LENGTHUNIT//"^2","count exp/"//LENGTHUNIT//"^2","count max/"//LENGTHUNIT//"^2"        &
         , "<intens>", " s.d. intens", " <t*>", " s.d. t*", " <quality>", " s.d. quality"
      do jj = 0, nbins
         dd = (jj + 0.5d0)*x/nbins
         ff = fbar(jj); tt = tbar(jj); qq = qbar(jj)
         fsig = 0.0d0; tsig = 0.0d0; qsig = 0.0d0

         nn = hist(jj)

         if (nn > 1) then
            ff = fbar(jj)/nn
            fvar = f2bar(jj)/nn
            tt = tbar(jj)/nn
            tvar = t2bar(jj)/nn
            qq = qbar(jj)/nn
            qvar = q2bar(jj)/nn
            fsig = sqrt(fvar - ff*ff)*sqrt(nn/real((nn + 1), kind=real64))
            tsig = sqrt(tvar - tt*tt)*sqrt(nn/real((nn + 1), kind=real64))
            qsig = sqrt(qvar - qq*qq)*sqrt(nn/real((nn + 1), kind=real64))
         end if
         if (jj < nbins) then
            write (unit=u, fmt='(i4,f9.2)', advance="no") jj, dd
         else
            write (unit=u, fmt='(a4,f9.2)', advance="no") "# > ", x
         end if

         !---    compute count per unit area, and min/max error bars
         !       note: will implement proper method here, for now use placeholder sqrt std dev.
         hpera = nn/aa
         if (nn == 0) then
            minhpera = 0
            maxhpera = 0
         else if (nn == 1) then
            minhpera = hpera
            maxhpera = hpera
         else
            minhpera = hpera*(1 - sqrt(1.0d0/nn))
            maxhpera = hpera*(1 + sqrt(1.0d0/nn))
         end if

         write (unit=u, fmt='(i8,3g16.4,100f14.5)', advance="yes") nn, minhpera, hpera, maxhpera, ff, fsig, tt, tsig, qq, qsig
      end do
      deallocate (hist)
      deallocate (fbar)
      deallocate (f2bar)
      deallocate (tbar)
      deallocate (t2bar)

      return
   end subroutine outputCount

   subroutine outputSpots(filename, ng, fitted_g2d, nx, ny, nmPerPixel, nDeadPixels, commentLine, diamScale, bg_sigma)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      output as standard .spots file

      character(len=*), intent(in)                     ::      filename        !   Filename in.
      integer, intent(in)                              ::      ng              !   number of spots fitted
      type(MultipleGaussians2d), intent(in)            ::      fitted_g2d      !   Object to hold 2D gaussians
      integer, intent(in)                              ::      nx, ny           !   width, height of input image, in pixels
      real(kind=real64), intent(in)                    ::      nmPerPixel      !   Length per pixel/nm
      integer, intent(in)                              ::      nDeadPixels     !   number of pixels out-of-scope
      character(len=*), intent(in)                     ::      commentLine     !   String to hold optional comment line
      real(kind=real64), intent(in)                    ::      diamscale       !   diameter = nmPerPixel*(2 sigma)*diamscale
      real(kind=real64), intent(in), optional           ::      bg_sigma        !   background intensity standard deviation

      real(kind=real64), parameter                     ::      PI = 3.1415926535897932384626433832795d0

      real(kind=real64)           ::      ww, hh, aa!,ss             !   width, height, area of input image, in nm
      integer                     ::      ii                      !   Indices of detected spots

      ww = nx*nmPerPixel
      hh = ny*nmPerPixel
      aa = (nx*ny - nDeadPixels)*nmPerPixel*nmPerPixel
      !ss = 1.0d0 ; if (present(bg_sigma)) ss = 1/bg_sigma     !   scale by noise level?

      open (unit=500, file=trim(filename), action="write")
      write (500, fmt='(a)') "# Clock "//trim(VERSION)
      write (500, fmt='(a)') "# "//trim(commentLine)
      write (500, fmt='(a,2f12.3,a)') "# ", ww, hh, " micrograph extent ("//LENGTHUNIT//")"

      write (500, fmt='(a2,f13.3,a)') "# ", aa, "  micrograph area ("//LENGTHUNIT//") discounting dead pixels"
      write (500, fmt='(a2,i8,a)') "# ", ng, "  number of spots"
                write (500,fmt='(a,100a12)') "# pos x ("//LENGTHUNIT//")  pos y ("//LENGTHUNIT//") diam 1 ("//LENGTHUNIT//") diam 2 ("//LENGTHUNIT//")"," angle(deg)"," intensity"," t* value"," quality"

      do ii = 1, getN(fitted_g2d)
         write (500, fmt='(100f12.3)') getX(fitted_g2d, ii)*nmPerPixel &
            , getY(fitted_g2d, ii)*nmPerPixel &
            , getSigma(fitted_g2d, ii, .true.)*nmPerPixel*(2*diamscale) &            !   2 for convert radius to diameter.
            , getSigma(fitted_g2d, ii, .false.)*nmPerPixel*(2*diamscale) &            !   2 for convert radius to diameter.
            , getTheta(fitted_g2d, ii)*180.0d0/PI &            !  convert radians to degrees
            , getF(fitted_g2d, ii) &
            , getT(fitted_g2d, ii) &
            , getQ(fitted_g2d, ii)
      end do
      close (unit=500)

      return
   end subroutine outputSpots

   subroutine outputPng(filename, fitted_g2d, img, diamScale, x, negative)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      simple dump of spots as rgb png
      character(len=*), intent(in)                     ::      filename        !   Filename of detected spots image outputted
      type(MultipleGaussians2d), intent(in)            ::      fitted_g2d      !   Object to hold 2D gaussians
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img             !   Array to hold image
      real(kind=real64), intent(in)                    ::      diamscale       !   diameter = nmPerPixel*(2 sigma)*diamscale
      real(kind=real64), intent(in)                    ::      x               !   maximum histogram bin
      logical, intent(in)                              ::      negative        !   Do image intensities need to be inverted?

      real(kind=real64), dimension(:, :, :), allocatable  ::      img_out         !   Array to hold output image

      integer             ::      Nx, Ny                                       !   Image x and y dimensions
      integer             ::      ix, iy                                       !   Image x and y indicies
      real(kind=real64)   ::      ss, s1, x0, y0, xx, yy, gg, ixx                    !   See below
      !reals
      !   ss  : Gaussian average diameter
      !   s1  : Gaussian major diameter
      !   x0  : Gaussian centre X coordinate
      !   y0  : Gaussian centre y coordinate
      !   xx  : Temp variable for calculating matrix DD, see Lib_Gaussian2d for definitions.
      !   yy  : Temp variable for calculating matrix DD, see Lib_Gaussian2d for definitions.
      !   gg  : Gaussian value at given pixel on the image
      !   ixx : Factor to convert gaussian radius to diameter and normalise to largest bin size.
      real(kind=real64), dimension(2, 2)        ::  DD  !matrix DD, see Lib_Gaussian2d for definition.
      integer             ::      ii                  !Spot index
      real(kind=real64), dimension(3)          ::  rgb_back, rgb_fore !RGB triplet for background and foregronud respectively.

      Nx = size(img, dim=1) !get image dimesions
      Ny = size(img, dim=2)
      allocate (img_out(3, 0:Nx - 1, 0:Ny - 1))
      ixx = nmperpixel*2*diamscale/x

      !Convert greyscale data to rgb
      do iy = 0, Ny - 1 !for all pixels
         do ix = 0, Nx - 1
            if (img(ix, iy) == LIB_G2D_IGNORE) then
               img_out(1:3, ix, iy) = 0 !set intensity to zero if pixel has ignore flag
            else if (negative) then
               img_out(1:3, ix, iy) = 1 - img(ix, iy)!invert intensity if negative flag is present.
            else
               img_out(1:3, ix, iy) = img(ix, iy) !if not just convert to RGB
            end if

         end do
      end do

      do ii = 1, getN(fitted_g2d)          ! for all detected spots
         if (ignore(fitted_g2d, ii)) cycle    !ignore if flagged
         x0 = getX(fitted_g2d, ii)            !Get gaussian parameters
         y0 = getY(fitted_g2d, ii)
         s1 = getSigma(fitted_g2d, ii, .true.)         !   larger radius
         ss = getSigma(fitted_g2d, ii)
         DD = getD(fitted_g2d, ii)

         rgb_fore = getRGB_double(ss*ixx)!get RGB values for overlayed spot
         do iy = max(0, nint(y0 - diamscale*s1)), min(Ny - 1, nint(y0 + diamscale*s1))
            do ix = max(0, nint(x0 - diamscale*s1)), min(Nx - 1, nint(x0 + diamscale*s1))
               xx = ix - x0
               yy = iy - y0
               gg = DD(1, 1)*xx*xx + 2*DD(2, 1)*xx*yy + DD(2, 2)*yy*yy

               gg = ss*(sqrt(2*gg) - diamscale)
               gg = exp(-2*gg*gg)
               rgb_back = img_out(1:3, ix, iy)
               img_out(1:3, ix, iy) = transparentColour(rgb_back, rgb_fore, gg)
            end do
         end do

      end do

      call write_rgb_png(filename, img_out)

      return
   end subroutine outputPng

   pure real(kind=real64) function automatic_xrange(diam)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given the input set of diameters, choose a sensible histogram range
      real(kind=real64), dimension(:), intent(in)           ::      diam !array of spot diameters

      !---    list of possible -x settings
      integer, parameter                               ::      NAUTOX = 25 !number of logarithmic bins
      real(kind=real64), dimension(NAUTOX), parameter   ::      AUTOX = (/1.0d4, 8.0d3, 5.0d3, 4.0d3, 3.0d3, 2.0d3, & !Logarithmic bins
                                                                          1.0d3, 8.0d2, 5.0d2, 4.0d2, 3.0d2, 2.0d2, &
                                                                          1.0d2, 8.0d1, 5.0d1, 4.0d1, 3.0d1, 2.0d1, &
                                                                          1.0d1, 8.0d0, 5.0d0, 4.0d0, 3.0d0, 2.0d0, 1.0d0/)
      real(kind=real64), parameter                     ::      MIN_HIST_CAPTURE = 0.95d0           !   minimum proportion of diameters I want to capture
      integer, dimension(NAUTOX)                       ::      autox_bin_count !logaritmic histogram bin count

      integer         ::      ii, jj, nDiam, minbinocc !diameter index, bin index, number of diameters, minimum bin occupancy

      !---    check for pathological cases
      nDiam = size(diam)
      if (nDiam == 0) then
         automatic_xrange = 1.0d0
         return
      else if (nDiam <= 4) then
         automatic_xrange = maxval(diam)
         return
      end if

      !---    find a histogram of the diameters
      autox_bin_count = 0
      do ii = 1, nDiam
         !   find the setting just larger than the diameter
         do jj = 2, NAUTOX
            if (diam(ii) > AUTOX(jj - 1)) then
               autox_bin_count(jj - 1) = autox_bin_count(jj - 1) + 1
               exit
            end if
         end do
      end do

      !---    now look at the largest settings found. I want to choose the largest setting
      !       which has > 1 diameter ( or 0.1%, whichever is bigger ) in it, subject to constraint that I want to <5% of diams outside range
      minbinocc = max(1, nint(nDiam*0.001d0))
      ii = floor((1 - MIN_HIST_CAPTURE)*nDiam)            !    This is the largest number of diams I'm allowed to ignore
      do jj = 2, NAUTOX

         if (autox_bin_count(jj) > minbinocc) then
            automatic_xrange = AUTOX(jj - 1)
            return
         end if

         ii = ii - autox_bin_count(jj)
         if (ii < 0) then
            !   If I look for a smaller setting, then I exclude too many diams
            automatic_xrange = AUTOX(jj)
            return
         end if

      end do

      !--- I don't think its possible to get here, but I have to return something!
      automatic_xrange = 50.0d0
      return
   end function automatic_xrange

   real(kind=real64) function automatic_binwidth(diam)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given the input set of diameters, choose a sensible histogram bin width
      !*      using the Freedman-Diaconis' rule
      !*          https://en.wikipedia.org/wiki/Histogram#Number_of_bins_and_width
      real(kind=real64), dimension(:), intent(in)           ::      diam !array of diameters

      real(kind=real64), dimension(size(diam))         ::      diam_sorted_list !array of diameters
      integer             ::      nDiam !number of diameters

      real(kind=real64)   ::      IQR         !   inter-quartile range

      !---    check for pathological cases
      nDiam = size(diam)
      if (nDiam == 0) then
         automatic_binwidth = 1.0d0
         return
      else if (nDiam <= 10) then
         automatic_binwidth = maxval(diam)/2
         return
      end if

      !---    sort the diameters into smallest to largest
      diam_sorted_list = diam
      call quicksort(diam_sorted_list)

      !---    find the interquartile range = q3 - q1
      IQR = diam_sorted_list(nint(nDiam*0.75d0)) - diam_sorted_list(nint(nDiam*0.25d0))

      automatic_binwidth = 2*IQR/(real(nDiam)**0.33333)

      return
   end function automatic_binwidth

   pure subroutine automatic_nBins(x, dx, nBins)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given the max bin size, and an ideal bin width
      !*      return a nice-looking number of bins. Adjust the bin width as appropriate
      real(kind=real64), intent(in)                    ::      x !max bin value
      real(kind=real64), intent(inout)                 ::      dx !bin width
      integer, intent(out)                             ::      nBins !number of bins

      integer                 ::      ii, nn, bestn, besti !bin edge index, current number of bins, best number of bins, best bin width index.
      real(kind=real64)       ::      dn, bestdn, nBin_FD
      !---    list of possible dx settings
      integer, parameter                               ::      NAUTODX = 26 !number of logarithmic bin edges
      real(kind=real64), dimension(NAUTODX), parameter  ::      AUTODX = (/1.0d4, 8.0d3, 5.0d3, 3.0d3, 2.0d3, & !logarithmic bin edges
                                                                           1.0d3, 8.0d2, 5.0d2, 3.0d2, 2.0d2, &
                                                                           1.0d2, 8.0d1, 5.0d1, 3.0d1, 2.0d1, &
                                                                           1.0d1, 8.0d0, 5.0d0, 3.0d0, 2.0d0, &
                                                                           1.0d0, 0.8d0, 0.5d0, 0.3d0, 0.2d0, 0.1d0/)
      nBin_FD = x/dx
      nBins = nint(nBin_FD)          !   this is the default

      !---    check pathological cases
      if (nBins <= 5) then
         dx = x/nBins
         return
      end if

      if (nBins > 20) then
         !   pretty much never want > 20 bins
         nBins = 20
         dx = x/nBins
         return
      end if

      !---    OK, at this point what I want is a nice bin width.
      !print *,"unscaled x,nBins,dx ",x,nBins,x/nBins
      bestdn = huge(1.0)
      besti = NAUTODX
      do ii = 1, NAUTODX

         nn = nint(x/AUTODX(ii))       !   how many bins would I have if I chose this bin with
         dn = (nn*AUTODX(ii) - x)
         if (abs(dn) > AUTODX(ii)*1.0d-4) cycle      !   don't get a good integer number of bins with this choice
         dn = (nn - nBin_FD)
         dn = dn*dn                      !   this is the square difference between this number of bins and my stated ideal
         if (dn <= bestdn) then
            bestdn = dn
            besti = ii
         end if
      end do

      nBins = nint(x/AUTODX(besti))

      dx = x/nBins
      return
   end subroutine automatic_nBins

   subroutine getScaledDiameters(g2dIn, diamsOut)
      !returns array of (real) spot diameters scaled to nm
      type(MultipleGaussians2d), intent(in)            ::      g2dIn      !   Object to hold 2D gaussians
      real(kind=real64), dimension(:), intent(inout)    ::      diamsOut

      integer                     ::          ii
      real(kind=real64)           ::          sig

      do ii = 1, nTotal
         sig = getSigma(fitted_g2d, ii)

         diamsOut(ii) = sig*nmperpixel*(2*diamscale)       !   2 for convert radius to diameter.

      end do

   end subroutine

end program countSpots
