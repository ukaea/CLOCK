
module Lib_Gaussian2d
!---^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_Gaussian2d from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      code which defines and fits a single 2d gaussian to input data
!*
!*      Gaussian defined internally as
!*          f(x,y) = f0 Exp[ - r . D r ]
!*      with
!*          r = ( x-x0 )
!*              ( y-y0 )
!*      and
!*          D = ( Dxx   Dxy )
!*              ( Dxy   Dyy )
!*
!*      note: has publically accessible components.
!*

   use Lib_GoldenSection                        !   golden section search for non-linear optimisation
   use iso_fortran_env
   implicit none
   private

   public          ::      Gaussian2d_ctor
   public          ::      report
   public          ::      delete

   public          ::      getDat                  !   return 6 parameters required to reconstruct gaussian in order
   public          ::      loglikelihood           !   return log likelihood assuming linear least squares fitting
   public          ::      findSigmaAndAngle       !   convert D matrix into major/minor radii

   public          ::      weight                  !   integral w = int g(r) d2r
   public          ::      overlap                 !   integral     int g1(r) g2(r) d2r / (w1 w2)
   public          ::      eccentricity            !   ratio of major:minor radii
   public          ::      add                     !   add gaussian to array of data
   public          ::      subtract                !   subtract from array of data
   public          ::      getRss                  !   residual sum of squares
   public          ::      finddrss                !   derivative of sum of squares wrt components of gaussian
   public          ::      getT                   !   compute t-statistic
   public          ::      translate               !   translate x,y centre
   public          ::      area                    !   Area under 1 std dev = pi sigma1 sigma2
   public          ::      fit                     !   fit gaussian to input array

   public          ::      computeRss              !   rss to 2 sigma range, also return number of pixels
   !---

   real(kind=real64), private, parameter         ::      PI = 3.141592653590d0
   integer(kind=int64), private, parameter       ::      BADF00D = int(z'BADF00D', kind=int64)
   real(kind=real64), public, parameter          ::      LIB_G2D_IGNORE = transfer((BADF00D + ishft(BADF00D, 32_int64)), 1.0d0)

   !---

   type, public     ::      Gaussian2d
      real(kind=real64)           ::      x0, y0               !   centre
      real(kind=real64)           ::      f0                  !   height
      real(kind=real64)           ::      Dxx, Dxy, Dyy         !   shape
   end type

   !---

   interface Gaussian2d_ctor
      module procedure Gaussian2d_null
      module procedure Gaussian2d_ctor0
      module procedure Gaussian2d_ctor1
      module procedure Gaussian2d_ctor2
   end interface

   interface delete
      module procedure delete0
   end interface

   interface report
      module procedure report0
   end interface

   interface getDat
      module procedure get0
   end interface

   interface loglikelihood
      module procedure loglikelihood0
   end interface

   interface fit
      module procedure fit0
      module procedure fit1
   end interface

   interface weight
      module procedure weight0
   end interface

   interface getT
      module procedure getT0
   end interface

   interface overlap
      module procedure overlap0
   end interface

   interface eccentricity
      module procedure eccentricity0
   end interface

   interface area
      module procedure area0
   end interface

   interface getRss
      module procedure getRss0
   end interface

   interface findDrss
      module procedure findDrss0
   end interface

   interface add
      module procedure add0
   end interface

   interface subtract
      module procedure subtract0
   end interface

   interface translate
      module procedure translate0
      module procedure translate1
      module procedure translate2
   end interface

contains
!---^^^^^^^^

   pure function Gaussian2d_null() result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns a null gaussian
      type(Gaussian2d)           ::      this
      this%f0 = 0
      this%x0 = 0
      this%y0 = 0
      this%Dxx = 0
      this%Dxy = 0
      this%Dyy = 0
      return
   end function Gaussian2d_null

   pure function Gaussian2d_ctor0(dat) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns a gaussian defined by its 6 parameters in order ( x0,y0,f0,Dxx,Dxy,Dyy )
      real(kind=real64), dimension(6), intent(in)       ::      dat
      type(Gaussian2d)           ::      this
      this%x0 = dat(1)
      this%y0 = dat(2)
      this%f0 = dat(3)
      this%Dxx = dat(4)
      this%Dxy = dat(5)
      this%Dyy = dat(6)
      return
   end function Gaussian2d_ctor0

   function Gaussian2d_ctor1(img) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns a gaussian defined by an intensity map
      real(kind=real64), dimension(:, :), intent(in)     ::      img
      type(Gaussian2d)           ::      this
      call fit(this, img)
      return
   end function Gaussian2d_ctor1

   function Gaussian2d_ctor2(img, imin, imax) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns a gaussian defined by an intensity map
      !*      with background level imin and maximum imax
      real(kind=real64), dimension(:, :), intent(in)     ::      img
      real(kind=real64), intent(in)                    ::      imin, imax
      type(Gaussian2d)           ::      this
      call fit(this, img, imin, imax)
      return
   end function Gaussian2d_ctor2

   subroutine delete0(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^
      !*      note: no dynamic memory
      type(Gaussian2d), intent(inout)    ::      this
      this = Gaussian2d_null()
      return
   end subroutine delete0

   !---

   subroutine report0(this, u, o)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      simple output. Defaults to unit 6 = screen.
      !*      optional Gaussian2d o determines left hand margin
      type(Gaussian2d), intent(in)     ::      this
      integer, intent(in), optional     ::      u, o
      integer     ::      uu, oo
      real(kind=real64)       ::      s1, s2, theta

      uu = 6; if (present(u)) uu = u
      oo = 0; if (present(o)) oo = o
      call findSigmaAndAngle(this, s1, s2, theta)

            write(unit=uu,fmt='(8(a,f10.4))') repeat(" ",oo)//"Gaussian2d [x0,y0,f0,s1,s2,theta = ",this%x0,",",this%y0,",",this%f0,",",s1,",",s2,",",theta*180/PI," deg )]"    
      return
   end subroutine report0

   !---

   pure function get0(this) result(dat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return data to construct gaussian from 6 parameters in order ( x0,y0,f0,Dxx,Dxy,Dyy )
      type(Gaussian2d), intent(in)         ::      this
      real(kind=real64), dimension(6)      ::      dat
      dat = (/this%x0, this%y0, this%f0, this%Dxx, this%Dxy, this%Dyy/)
      return
   end function get0

   !---

   function loglikelihood0(this, img) result(lambda)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given the input density function, return the log likelihood ( ignoring the noise term )
      !*      lambda = - n log (rss/n)
      type(Gaussian2d), intent(in)                     ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img
      real(kind=real64)                   ::      lambda

      integer             ::      Nx, Ny, n
      integer             ::      ix, iy
      real(kind=real64)   ::      gg, ss, xx, yy, rss

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      n = 0               !   pixel count
      rss = 0             !   residual sum of squares
      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            if (img(ix, iy) == LIB_G2D_IGNORE) cycle
            xx = (ix - this%x0)
            yy = (iy - this%y0)
            gg = this%Dxx*xx*xx + 2*this%Dxy*xx*yy + this%Dyy*yy*yy
            gg = Exp(-gg)
            ss = this%f0*gg - img(ix, iy)
            rss = rss + ss*ss
            n = n + 1
         end do
      end do

      if (n > 0) then
         lambda = -n*log(rss/n)
      else
         lambda = 0
      end if

      return
   end function loglikelihood0

   pure real(kind=real64) function getRss0(this, img)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given the input density function, return the residual sum of squares
      type(Gaussian2d), intent(in)                     ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img

      integer             ::      Nx, Ny
      integer             ::      ix, iy
      real(kind=real64)   ::      gg, ss, xx, yy

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      getRss0 = 0
      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            if (img(ix, iy) == LIB_G2D_IGNORE) cycle
            xx = (ix - this%x0)
            yy = (iy - this%y0)
            gg = this%Dxx*xx*xx + 2*this%Dxy*xx*yy + this%Dyy*yy*yy
            gg = Exp(-gg)
            ss = this%f0*gg - img(ix, iy)
            getRss0 = getRss0 + ss*ss
         end do
      end do

      return
   end function getRss0

   !---

   subroutine findDrss0(this, img, rss, drss)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return the first derivatives of the residual sum of squares with respect to change of parameters
      type(Gaussian2d), intent(in)                     ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img
      real(kind=real64), intent(out)                   ::      rss
      real(kind=real64), dimension(6), intent(out)      ::      drss

      integer             ::      Nx, Ny
      integer             ::      ix, iy
      real(kind=real64)   ::      ff, gg, ss, xx, yy

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      rss = 0
      drss = 0
      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            if (img(ix, iy) == LIB_G2D_IGNORE) cycle
            xx = (ix - this%x0)
            yy = (iy - this%y0)
            gg = this%Dxx*xx*xx + 2*this%Dxy*xx*yy + this%Dyy*yy*yy
            gg = Exp(-gg)
            ff = this%f0*gg
            ss = ff - img(ix, iy)

            rss = rss + ss*ss

            drss(1) = drss(1) + 2*ss*(2*ff*(this%Dxx*xx + this%Dxy*yy))        !    = d f/d x
            drss(2) = drss(2) + 2*ss*(2*ff*(this%Dxy*xx + this%Dyy*yy))
            drss(3) = drss(3) + 2*ss*(gg)        !    = d f/d f0
            drss(4) = drss(4) + 2*ss*(-ff*xx*xx)        !    = d f/d Dxx
            drss(5) = drss(5) + 2*ss*(-2*ff*xx*yy)
            drss(6) = drss(6) + 2*ss*(-ff*yy*yy)

         end do
      end do

      return
   end subroutine findDrss0

   !---

   subroutine fit0(this, img, imin, imax, mask, dbg)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find the best fit gaussian to the intensity map
      !*      imin is taken to be the background level
      !*      anything above imax is taken to be out-of-bounds
      !*      anything not masked is out-of-bounds
      type(Gaussian2d), intent(inout)                  ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img
      real(kind=real64), intent(in)                    ::      imin, imax
      logical, dimension(0:, 0:), intent(in), optional    ::      mask
      logical, intent(in), optional                     ::      dbg

      logical, dimension(:), allocatable                ::      hasx, hasy
      real(kind=real64), dimension(:, :), allocatable    ::      img_clean
      integer             ::      Nx, Ny
      integer             ::      ix, iy, jx, jy, kx, ky, nPx
      real(kind=real64)   ::      ff

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      !---    find the extent of the readable part of the image
      allocate (hasx(0:Nx - 1))
      allocate (hasy(0:Ny - 1))
      hasx = .false.
      hasy = .false.
      nPx = 0
      if (present(mask)) then
         if (imax > imin) then
            do iy = 0, Ny - 1
               do ix = 0, Nx - 1
                  if (.not. mask(ix, iy)) cycle
                  ff = img(ix, iy)
                  if ((ff == LIB_G2D_IGNORE) .or. (ff > imax)) cycle
                  hasx(ix) = .true.
                  hasy(iy) = .true.
                  nPx = nPx + 1
               end do
            end do
         else
            do iy = 0, Ny - 1
               do ix = 0, Nx - 1
                  if (.not. mask(ix, iy)) cycle
                  ff = img(ix, iy)
                  if ((ff == LIB_G2D_IGNORE) .or. (ff < imax)) cycle
                  hasx(ix) = .true.
                  hasy(iy) = .true.
                  nPx = nPx + 1
               end do
            end do
         end if
      else
         if (imax > imin) then
            do iy = 0, Ny - 1
               do ix = 0, Nx - 1
                  ff = img(ix, iy)
                  if ((ff == LIB_G2D_IGNORE) .or. (ff > imax)) cycle
                  hasx(ix) = .true.
                  hasy(iy) = .true.
                  nPx = nPx + 1
               end do
            end do
         else
            do iy = 0, Ny - 1
               do ix = 0, Nx - 1
                  ff = img(ix, iy)
                  if ((ff == LIB_G2D_IGNORE) .or. (ff < imax)) cycle
                  hasx(ix) = .true.
                  hasy(iy) = .true.
                  nPx = nPx + 1
               end do
            end do
         end if
      end if

      ix = Nx; jx = -1
      do kx = 0, Nx - 1
         if (hasx(kx)) then
            ix = kx; exit
         end if
      end do
      do kx = Nx - 1, 0, -1
         if (hasx(kx)) then
            jx = kx; exit
         end if
      end do

      iy = Ny; jy = -1
      do ky = 0, Ny - 1
         if (hasy(ky)) then
            iy = ky; exit
         end if
      end do
      do ky = Ny - 1, 0, -1
         if (hasy(ky)) then
            jy = ky; exit
         end if
      end do

      if ((ix > jx) .or. (iy > jy) .or. (nPx <= 4)) then    !   not enough valid points
         !---    this is the null guess: zero intensity, centred in the box
         this%x0 = Nx/2
         this%y0 = Ny/2
         this%f0 = 0.0d0
         this%Dxx = 2/(this%x0*this%x0)
         this%Dxy = 0
         this%Dyy = 2/(this%y0*this%y0)
         return
      end if

      !---    cut out the centre of the image, and set any points out of bounds
      allocate (img_clean(ix:jx, iy:jy))
      img_clean = LIB_G2D_IGNORE
      if (present(mask)) then
         if (imax > imin) then
            do ky = iy, jy
               do kx = ix, jx
                  if (.not. mask(kx, ky)) cycle
                  ff = img(kx, ky)
                  if (.not. (ff == LIB_G2D_IGNORE) .or. (ff > imax)) img_clean(kx, ky) = ff - imin
               end do
            end do
         else
            do ky = iy, jy
               do kx = ix, jx
                  if (.not. mask(kx, ky)) cycle
                  ff = img(kx, ky)
                  if (.not. (ff == LIB_G2D_IGNORE) .or. (ff < imax)) img_clean(kx, ky) = ff - imin
               end do
            end do
         end if
      else
         if (imax > imin) then
            do ky = iy, jy
               do kx = ix, jx
                  ff = img(kx, ky)
                  if (.not. (ff == LIB_G2D_IGNORE) .or. (ff > imax)) img_clean(kx, ky) = ff - imin
               end do
            end do
         else
            do ky = iy, jy
               do kx = ix, jx
                  ff = img(kx, ky)
                  if (.not. (ff == LIB_G2D_IGNORE) .or. (ff < imax)) img_clean(kx, ky) = ff - imin
               end do
            end do
         end if
      end if

      if (present(dbg)) then
         print *, "Lib_Gaussians2d::fit0 info - offset ", ix, iy
         call fit1(this, img_clean, dbg)
      else
         call fit1(this, img_clean)
      end if

      this%x0 = this%x0 + ix
      this%y0 = this%y0 + iy

      deallocate (img_clean)
      deallocate (hasx)
      deallocate (hasy)

      return
   end subroutine fit0

   subroutine fit1(this, img, dbg)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find the best fit gaussian to the intensity map
      type(Gaussian2d), intent(inout)                  ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img
      logical, intent(in), optional                     ::      dbg

      integer                         ::      Nx, Ny
      integer                         ::      ix, iy

      real(kind=real64)               ::      TT, Tx, Ty, Txx, Txy, Tyy
      real(kind=real64)               ::      ff

      real(kind=real64)               ::      rss, oldrss, rss0
      real(kind=real64), dimension(6)  ::      dat0

      integer                         ::      ii

      integer, parameter               ::      NLOOPS = 100
      real(kind=real64), parameter     ::      TOL = 1.0d-6

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      if (present(dbg)) then
         print *, "Lib_Gaussian2d::fit1 info - "
         do iy = 0, Ny - 1
            write (*, fmt='(100l7)') img(0:min(99, Nx - 1), iy) /= LIB_G2D_IGNORE
         end do

         do iy = 0, Ny - 1
            write (*, fmt='(100f7.3)') img(0:min(99, Nx - 1), iy)
         end do
         print *, ""
      end if

      !---    this is the null guess: zero intensity, centred in the box
      ff = -huge(1.0)
      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            if ((img(ix, iy) /= LIB_G2D_IGNORE) .and. (img(ix, iy) > ff)) then
               ff = img(ix, iy)
               this%f0 = ff
               this%x0 = ix
               this%y0 = iy
            end if
         end do
      end do
      this%Dxx = 8.0d0/(Nx*Nx)          !   sigma = Nx/4 -> Dxx = 1/(2 sigma sigma)
      this%Dxy = 0
      this%Dyy = 8.0d0/(Ny*Ny)

      if (present(dbg)) then
         oldrss = getRss(this, img)
         print *, "Lib_Gaussian2d::fit1 info - first guess (before)", oldrss
         call report(this)
      end if

      if (present(dbg)) then
         call fit2(this, img, rss0, dbg)
      else
         call fit2(this, img, rss0)
      end if
      dat0 = getDat(this)

      !---    step 1: find the mean and variance of the data
      TT = 0
      Tx = 0
      Ty = 0
      Txx = 0
      Txy = 0
      Tyy = 0
      ii = 0
      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            ff = img(ix, iy)
            if (ff == LIB_G2D_IGNORE) cycle
            ff = max(0.0d0, ff)
            TT = TT + ff
            Tx = Tx + ff*ix
            Ty = Ty + ff*iy
            Txx = Txx + ff*ix*ix
            Txy = Txy + ff*ix*iy
            Tyy = Tyy + ff*iy*iy
            ii = ii + 1
         end do
      end do

      if (present(dbg)) print *, "Lib_Gaussian2d::fit1 info - T,Tx,Ty,Txx,Tyy,Txy,npx ", TT, Tx, Ty, Txx, Txy, Tyy, ii

      !---    quick escape for dodgy input data
      if (TT < 1d-12) return !   there are no good pixels, or there is no bias one way or the other, or the peak is negative. The fit is obvious, and unhelpful

      !---    step 2: initial guess should be that the means and variances of the gaussian match the data
      !       I can do this assuming all the weight of the gaussian is inside the region 0:Nx and so integrating -infty:infty.
      !       this gives solvable analytic expressions.
      this%x0 = Tx/TT
      this%y0 = Ty/TT
      ff = 2*(Txx*Ty*Ty + Tyy*Tx*Tx - 2*Tx*Txy*Ty - TT*(Txx*Tyy - Txy*Txy))
      if (abs(ff) > 1d-12) then
         ff = 1/ff
         this%Dxx = TT*(Ty*Ty - TT*Tyy)*ff
         this%Dxy = -TT*(Tx*Ty - TT*Txy)*ff
         this%Dyy = TT*(Tx*Tx - TT*Txx)*ff
         ff = this%Dxx*this%Dyy - this%Dxy*this%Dxy
         if (ff > 0) then
            this%f0 = min(maxval(img), sqrt(ff)*TT/PI)       !   note: this%f0 is earlier initialised to maxval(img)
         end if
      else
         return
      end if

      if (present(dbg)) then
         call fit2(this, img, rss, dbg)
      else
         call fit2(this, img, rss)
      end if

      if (rss0 < rss) then
         !   give up "sophisticated" guess, return to "simple"
         this = Gaussian2d_ctor(dat0)
      end if

      return
   end subroutine fit1

   subroutine fit2(this, img, rss, dbg)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      improve the fit gaussian to the intensity map using Golden section search
      !*      on input, img should be cleaned up to have zero background level.
      type(Gaussian2d), intent(inout)                  ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img
      real(kind=real64), intent(out)                   ::      rss
      logical, intent(in), optional                     ::      dbg

      integer             ::      ix

      real(kind=real64)   ::      ff

      real(kind=real64)   ::      oldrss
      real(kind=real64), dimension(6)  ::      dat, drss
      type(GoldenSection) ::      gold
      logical             ::      isConverged, isWithinTimeLimit
      real(kind=real64)   ::      x1, x2, x3, y1, y2, y3
      integer             ::      ii

      integer, parameter               ::      NLOOPS = 100
      real(kind=real64), parameter     ::      TOL = 1.0d-6

      if (present(dbg)) then
         oldrss = getRss(this, img)
         print *, "Lib_Gaussian2d::fit2 info - (before)", oldrss
         call report(this)
      end if

      !---    try to fit
      do ii = 1, NLOOPS

         dat = getDat(this)
         call findDrss(this, img, rss, drss)
         if (ii == 1) oldrss = 2*rss

         if (mod(ii, 2) == 0) drss(4:6) = 0.0d0

         ff = drss(1)*drss(1) + drss(2)*drss(2) + drss(3)*drss(3) + drss(4)*drss(4) + drss(5)*drss(5) + drss(6)*drss(6)
         if (abs(ff) < 1.0d-8) exit

         !---    have the derivative direction, now need to estimate the length of the deriv vector
         x2 = -0.001*rss/ff
         y1 = huge(1.0)
         do ix = 0, 20
            y2 = getRss(Gaussian2d_ctor(dat + x2*drss), img)
            if (y2 > y1) then
               y1 = y2; x1 = x2; exit
            else
               y1 = y2
            end if
            x2 = x2*2
         end do
         x3 = 0; y3 = rss

         gold = GoldenSection_ctor(x1, x3, y1, y3)
         do
            call nextPoint(gold, x2)
            y2 = getRss(Gaussian2d_ctor(dat + x2*drss), img)
            call minimise(gold, x2, y2, isConverged, isWithinTimeLimit)
            if (isConverged) exit
            if (.not. isWithinTimeLimit) exit
         end do

         this = Gaussian2d_ctor(dat + x2*drss)
         rss = getRss(this, img)

         if (abs(rss - oldrss) < rss*TOL) exit
         oldrss = rss

      end do

      if (present(dbg)) then
         print *, "Lib_Gaussian2d::fit2 info - (after)", rss
         call report(this)
      end if

      return
   end subroutine fit2

   subroutine findSigmaAndAngle(this, s1, s2, theta, ok)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(Gaussian2d), intent(in)             ::      this
      real(kind=real64), intent(out)           ::      s1, s2, theta
      logical, intent(out), optional            ::      ok
      real(kind=real64)                       ::      lambda1, lambda2, ton2, det, zz

      if (abs(this%Dxy) < 1.0d-14) then

         if (this%Dxx >= this%Dyy) then
            s1 = 1/sqrt(2*this%Dyy)
            s2 = 1/sqrt(2*this%Dxx)
            theta = 0.0d0
         else
            s1 = 1/sqrt(2*this%Dxx)
            s2 = 1/sqrt(2*this%Dyy)
            theta = PI/2
         end if
         if (present(ok)) ok = .false.
         return

      end if

      ton2 = (this%Dxx + this%Dyy)/2                      !   trace /2
      det = this%Dxx*this%Dyy - this%Dxy*this%Dxy                !   determinant

      if (det < 0) then
         !   make a fudge - add to the diagonal sufficient to push the det positive
         zz = sqrt((this%Dxx - this%Dyy)*(this%Dxx - this%Dyy) + 4*this%Dxy*this%Dxy) - this%Dxx - this%Dyy
         ton2 = ton2 + zz
         det = det + (this%Dxx + this%Dyy)*zz + zz*zz
      end if

      zz = ton2*ton2 - det
      zz = sqrt(max(0.0d0, zz))

      lambda1 = ton2 - zz                         !
      lambda2 = ton2 + zz

      if (min(lambda1, lambda2) <= 1.0d-12) then
         s1 = 0.0
         s2 = 0.0
         theta = 0.0d0
         if (present(ok)) ok = .false.
         return
      end if

      if (lambda2 < lambda1) then                   !   if det < 0  then might need reorder
         lambda1 = ton2 + zz
         lambda2 = ton2 - zz
      end if

      theta = atan2(this%Dxy, lambda1 - this%Dyy)

      if (theta < 0) then
         theta = theta + PI
      else if (theta > PI) then
         theta = theta - PI
      end if

      s1 = 1/sqrt(2*lambda1)
      s2 = 1/sqrt(2*lambda2)

      return
   end subroutine findSigmaAndAngle

   !---

   pure real(kind=real64) function weight0(this, Nx, Ny)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      compute weight of gaussian in a rectangular box
      type(Gaussian2d), intent(in)         ::      this
      integer, intent(in), optional         ::      Nx, Ny

      integer             ::      ix, iy
      real(kind=real64)   ::      ff, xx, yy

      if (present(Nx)) then
         weight0 = 0

         do iy = 0, Ny - 1
            do ix = 0, Nx - 1
               xx = (ix - this%x0)
               yy = (iy - this%y0)
               ff = this%Dxx*xx*xx + 2*this%Dxy*xx*yy + this%Dyy*yy*yy
               ff = Exp(-ff)
               weight0 = weight0 + ff

            end do
         end do
         weight0 = weight0*this%f0
      else
         !   analytic over plane

         ff = this%Dxx*this%Dyy - this%Dxy*this%Dxy   !  note this is the determinant. sigma1 sigma2 = 1/(2 sqrt(ff))
         if (ff > 1.0d-12) then
            weight0 = PI*this%f0/sqrt(ff)
            ! weight0 = sqrt(2*PI)*this%f0/(ff**0.25)           !   a t- measure = integral f(x) d2x / sqrt(area)
         else
            weight0 = 0.0d0
         end if
      end if

      return
   end function weight0

   real(kind=real64) function getT0(this, img, sigma)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given the input density function, return the t-value
      !*          t = |<f>| / (sigma/sqrt(n))
      !*      note that on input, the background level should be zero
      type(Gaussian2d), intent(in)                             ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)           ::      img
      real(kind=real64), intent(in)                            ::      sigma       !   background noise std dev.

      integer             ::      Nx, Ny
      integer             ::      ix, iy, nn
      real(kind=real64)   ::      gg, fsum, xx, yy

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      nn = 0
      fsum = 0.0d0
      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            xx = ix - this%x0
            yy = iy - this%y0
            gg = this%Dxx*xx*xx + 2*this%Dxy*xx*yy + this%Dyy*yy*yy
            if (gg > 0.5d0) cycle               !   only search to 1 sigma
            fsum = fsum + Exp(-gg)
            nn = nn + 1
         end do
      end do

      if (nn > 0) then
         getT0 = abs(this%f0*fsum)/(sqrt(nn*1.0d0)*sigma)
      else
         getT0 = 0
      end if

      return
   end function getT0

   pure real(kind=real64) function overlap0(this, that, Nx, Ny)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      compute overlap of two gaussians in a rectangular box
      type(Gaussian2d), intent(in)         ::      this, that
      integer, intent(in)                  ::      Nx, Ny

      integer             ::      ix, iy
      real(kind=real64)   ::      f1, xx, yy, f2, w1, w2, w12

      w1 = 0
      w2 = 0
      w12 = 0
      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            xx = (ix - this%x0)
            yy = (iy - this%y0)
            f1 = this%Dxx*xx*xx + 2*this%Dxy*xx*yy + this%Dyy*yy*yy
            f1 = Exp(-f1)
            w1 = w1 + f1            !   note that this%f0 cancels out.

            xx = (ix - that%x0)
            yy = (iy - that%y0)
            f2 = that%Dxx*xx*xx + 2*that%Dxy*xx*yy + that%Dyy*yy*yy
            f2 = Exp(-f2)
            w2 = w2 + f2

            w12 = w12 + f1*f2
         end do
      end do
      if (w1*w2 > 0) then
         overlap0 = w12/(w1*w2)
      else
         overlap0 = 0
      end if

      return
   end function overlap0

   real(kind=real64) function eccentricity0(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns the ratio of the high:low widths
      type(Gaussian2d), intent(in)         ::      this

      real(kind=real64)                       ::      lambda1, lambda2, ton2, det, zz

      if (abs(this%Dxy) < 1.0d-14) then

         if (this%Dxx >= this%Dyy) then
            eccentricity0 = sqrt(abs(this%Dxx/this%Dyy))
         else
            eccentricity0 = sqrt(abs(this%Dyy/this%Dxx))
         end if
         return

      end if

      ton2 = (this%Dxx + this%Dyy)/2                      !   trace /2
      det = this%Dxx*this%Dyy - this%Dxy*this%Dxy                !   determinant

      if (det < 0) then
         !   make a fudge - add to the diagonal sufficient to push the det positive
         zz = sqrt((this%Dxx - this%Dyy)*(this%Dxx - this%Dyy) + 4*this%Dxy*this%Dxy) - this%Dxx - this%Dyy
         ton2 = ton2 + zz
         det = det + (this%Dxx + this%Dyy)*zz + zz*zz
      end if

      zz = ton2*ton2 - det
      zz = sqrt(max(0.0d0, zz))

      lambda1 = (ton2 - zz)
      lambda2 = (ton2 + zz)

      if (lambda1 > lambda2) then
         eccentricity0 = sqrt(abs(lambda1)/max(1.0d-8, lambda2))
      else
         eccentricity0 = sqrt(abs(lambda2)/max(1.0d-8, lambda1))
      end if

      return
   end function eccentricity0

   real(kind=real64) function area0(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns the area above 1 std dev = pi sigma1 sigma2
      type(Gaussian2d), intent(in)         ::      this

      real(kind=real64)       ::      det

      det = this%Dxx*this%Dyy - this%Dxy*this%Dxy
      if (det > 0) then
         area0 = PI/(2*sqrt(det))
      else
         area0 = 0.0d0
      end if

      return
   end function area0

   subroutine subtract0(img, this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      remve this gaussian from the image
      type(Gaussian2d), intent(in)                         ::  this
      real(kind=real64), dimension(0:, 0:), intent(inout)    ::  img

      integer             ::      Nx, Ny
      integer             ::      ix, iy
      real(kind=real64)   ::      ff, gg, xx, yy

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            ff = img(ix, iy)
            if (ff == LIB_G2D_IGNORE) cycle
            xx = (ix - this%x0)
            yy = (iy - this%y0)
            gg = this%Dxx*xx*xx + 2*this%Dxy*xx*yy + this%Dyy*yy*yy
            if (gg > 4.5d0) cycle           !   only consider up to 3 sigma.
            gg = Exp(-gg)
            gg = this%f0*gg
            img(ix, iy) = ff - gg
         end do
      end do

      return
   end subroutine subtract0

   subroutine add0(img, this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      add this gaussian to the image
      type(Gaussian2d), intent(in)                         ::  this
      real(kind=real64), dimension(0:, 0:), intent(inout)    ::  img

      integer             ::      Nx, Ny
      integer             ::      ix, iy
      real(kind=real64)   ::      ff, gg, xx, yy

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            ff = img(ix, iy)
            if (ff == LIB_G2D_IGNORE) cycle
            xx = (ix - this%x0)
            yy = (iy - this%y0)
            gg = this%Dxx*xx*xx + 2*this%Dxy*xx*yy + this%Dyy*yy*yy
            if (gg > 4.5d0) cycle           !   only consider up to 3 sigma.
            gg = Exp(-gg)
            img(ix, iy) = ff + this%f0*gg
         end do
      end do

      return
   end subroutine add0

   pure subroutine translate0(this, x)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      translate centre of gaussian
      type(Gaussian2d), intent(inout)                  ::      this
      real(kind=real64), dimension(2), intent(in)       ::      x
      this%x0 = this%x0 + x(1)
      this%y0 = this%y0 + x(2)
      return
   end subroutine translate0

   pure subroutine translate1(this, x, y)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      translate centre of gaussian
      type(Gaussian2d), intent(inout)                  ::      this
      real(kind=real64), intent(in)                    ::      x, y
      this%x0 = this%x0 + x
      this%y0 = this%y0 + y
      return
   end subroutine translate1

   pure subroutine translate2(this, x, y)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      translate centre of gaussian
      type(Gaussian2d), intent(inout)                  ::      this
      integer, intent(in)                              ::      x, y
      this%x0 = this%x0 + x
      this%y0 = this%y0 + y
      return
   end subroutine translate2

   pure subroutine computeRss(this, img, rss, npx)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given the input density function, return the residual sum of squares
      !*      in a 2 sigma range. Return also number of pixels used
      type(Gaussian2d), intent(in)                     ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img
      real(kind=real64), intent(out)                   ::      rss
      integer, intent(out)                             ::      npx
      integer             ::      Nx, Ny
      integer             ::      ix, iy
      real(kind=real64)   ::      gg, ss, xx, yy

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      rss = 0
      npx = 0
      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            if (img(ix, iy) == LIB_G2D_IGNORE) cycle
            xx = (ix - this%x0)
            yy = (iy - this%y0)
            gg = this%Dxx*xx*xx + 2*this%Dxy*xx*yy + this%Dyy*yy*yy
            if (gg > 2.0d0) cycle             !   only search to 2 sigma range
            gg = Exp(-gg)
            ss = this%f0*gg - img(ix, iy)
            rss = rss + ss*ss
            npx = npx + 1
         end do
      end do

      return
   end subroutine computeRss

end module Lib_Gaussian2d
