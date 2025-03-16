
module Lib_FitGaussians1d
!---^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_FitGaussians1d from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      A code to fit multiple 1d gaussians to a single array of values
!*
!*      A Gaussian is defined as
!*
!*          g(x) = (f^2) Exp[ -(x-m)^2 / s^2 ]
!*
!*      with f,m, and s being saved as internal variables.
!*      Note that this definition means the Gaussian is always positive.
!*
!*      Multiple gaussians are simply added
!*          G(x) = sum_j g_j(x)
!*
!*      Note that we actually work with an array of integer x- values 1,2,3...n
!*      in order to fit to an input data curve.
!*      Therefore the real line over which the data is defined is actually 0.5<=x<n+0.5
!*      with the knot points falling in the centre of the intervals.
!*
!*      We can fit to the input data in_dat by trying a few possible
!*      initial conditions for starting a set of Gaussians, then using a 1st derivative search to find the
!*      (local) minimum least-squares fit.
!*      Then we can look at different numbers of Gaussians, and find the one with the lowest AIC
!*

   use Lib_GoldenSection
   use iso_fortran_env
   implicit none
   private

   !---

   public      ::      FitGaussians1d_ctor         !   fg1d = FitGaussians1d_ctor( input_array, max_gaussians )
   public      ::      delete                      !   removes dynamic memory
   public      ::      report                      !   call report( fg1d [,verbose] [,u] [,o] ) simple dump to unit [u=screen]
   public      ::      fit                         !   best fit using AIC criterion
   public      ::      AICc                        !   best fit using AIC criterion
   public      ::      getM                        !   returns the number of gaussians in the fit
   public      ::      getFMuSigAx                 !   getFMuSigAx(this,j) returns array (/ f,mu,sig,Area,<x> /)
   public      ::      intgrl                      !   integrate single gaussian between 1/2 and n+1/2, or between 1/2 and x

   !---

   type, public     ::      FitGaussians1d
      private
      integer                                 ::      n               !   number of data points to fit
      integer                                 ::      m               !   number of gaussians used
      integer                                 ::      mmax            !   max number of gaussians to consider
      real(kind=real64), dimension(:), pointer  ::      in_dat          !   data to fit (1:n)
      real(kind=real64), dimension(:), pointer  ::      mu, sig, f        !   mean, inverse square width, (sqrt) max intensity (1:mmax)
   end type

   !---

   interface FitGaussians1d_ctor
      module procedure FitGaussians1d_null
      module procedure FitGaussians1d_ctor0
   end interface

   interface delete
      module procedure delete0
   end interface

   interface report
      module procedure report0
      module procedure report1
   end interface

   interface func
      module procedure func0
      module procedure func1
   end interface

   interface getM
      module procedure getM0
   end interface

   interface getFMuSigAx
      module procedure getFMuSigAx0
   end interface

   interface fit
      module procedure fit0
      module procedure fit0a
   end interface

!---

contains
!---^^^^^^^^

   function FitGaussians1d_null() result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      default null constructor - sets to zero.
      type(FitGaussians1d)           ::      this
      this%n = 0
      this%mmax = 0
      this%m = 0
      nullify (this%in_dat)
      nullify (this%mu)
      nullify (this%sig)
      nullify (this%f)
      return
   end function FitGaussians1d_null

   function FitGaussians1d_ctor0(in_dat, mmax) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      standard constructor, takes in the data to fit and the maximum number of Gaussians to try
      type(FitGaussians1d)           ::      this
      real(kind=real64), dimension(:), intent(in)       ::      in_dat
      integer, intent(in)             ::      mmax

      this%n = size(in_dat)
      allocate (this%in_dat(this%n))
      this%in_dat(1:this%n) = in_dat(1:this%n)
      this%mmax = mmax
      this%m = 0
      allocate (this%mu(this%mmax))
      allocate (this%sig(this%mmax))
      allocate (this%f(this%mmax))
      this%mu = 0
      this%sig = 1.00
      this%f = 0
      return
   end function FitGaussians1d_ctor0

   subroutine delete0(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^
      !       deallocate dynamic memory
      type(FitGaussians1d), intent(inout)    ::      this
      if (this%mmax == 0) return
      deallocate (this%in_dat)
      deallocate (this%mu)
      deallocate (this%sig)
      deallocate (this%f)
      this = FitGaussians1d_null()
      return
   end subroutine delete0

   !---

   subroutine report0(this, u, o)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      simple output to unit u ( default screen ) with margin offset o
      type(FitGaussians1d), intent(in)     ::      this
      integer, intent(in), optional         ::      u, o
      integer                             ::      uu, oo
      integer         ::      jj
      uu = 6; if (present(u)) uu = u
      oo = 0; if (present(o)) oo = o
      write (unit=uu, fmt='(4(a,i6))') repeat(" ", oo)//"FitGaussians1d[n=", this%n, ",m=", this%m, "]"
      do jj = 1, this%m
                write(unit=uu,fmt='(a,i6,5(a,f16.8))') repeat(" ",oo+4)//"Gaussian ",jj," [f,mu,sig = ",this%f(jj)*this%f(jj),",",this%mu(jj),",",this%sig(jj),", A=",intgrl(this,jj),"]"                                                                                               
      end do
      return
   end subroutine report0

   subroutine report1(this, verbose, u, o)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      output to unit u ( default screen ) with margin offset o
      !*      also returns the input array and the fitted contributions from each gaussian
      type(FitGaussians1d), intent(in)     ::      this
      logical, intent(in)                  ::      verbose
      integer, intent(in), optional         ::      u, o
      integer                             ::      uu, oo

      integer             ::      ii, jj

      uu = 6; if (present(u)) uu = u
      oo = 0; if (present(o)) oo = o
      call report0(this, uu, oo)

      if (verbose) then
         write (unit=uu, fmt='(a,a6,100a16)', advance="no") repeat(" ", oo + 4), " bin ", " input "
         do jj = 1, this%m
            write (unit=uu, fmt='(a12,i4)', advance="no") " Gaussian ", jj
         end do
         write (unit=uu, fmt='(a16)', advance="yes") " sum"
         do ii = 1, this%n
 write (unit=uu, fmt='(a,i6,100f16.8)') repeat(" ", oo + 4), ii, this%in_dat(ii), (func(this, ii, jj), jj=1, this%m), func(this, ii)
         end do
      end if
      write (unit=uu, fmt='(a)') ""
      write (unit=uu, fmt='(2(a,f16.8))') repeat(" ", oo + 4)//"RSS = ", rss(this), " AICc = ", AICc(this)
      return
   end subroutine report1

   !---

   pure integer function getM0(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      accessor for number of Gaussians
      type(FitGaussians1d), intent(in)     ::      this
      getM0 = this%m
   end function getM0

   pure function getFMuSigAx0(this, j) result(FMuSigAx)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      accessor for result of jth Gaussian.
      !*      returns an array FMuSigAx( 1:5 ) = (/ max_gaussian_height , centre_of_gaussian , width_of_gaussian , area_under_gaussian , centre_of_mass /)
      type(FitGaussians1d), intent(in)     ::      this
      integer, intent(in)                  ::      j
      real(kind=real64), dimension(5)      ::      FMuSigAx
      real(kind=real64)       ::      aa, xx
      aa = intgrl(this, j)
      xx = xintgrl(this, j)
      FMuSigAx(1:5) = (/this%f(j)*this%f(j), this%mu(j), this%sig(j), aa, xx/max(aa, 1.0d-16)/)
      return
   end function getFMuSigAx0

   !---

   !---

   pure real(kind=real64) function func0(this, i)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns value of function at data point i
      type(FitGaussians1d), intent(in)     ::      this
      integer, intent(in)                  ::      i

      real(kind=real64)   ::      dd
      integer             ::      jj

      func0 = 0
      do jj = 1, this%m
         dd = (i - this%mu(jj))/this%sig(jj)
         dd = this%f(jj)*this%f(jj)*exp(-dd*dd/2)
         func0 = func0 + dd
      end do

      return
   end function func0

   pure real(kind=real64) function func1(this, i, j)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns value of function at data point i for gaussian j
      type(FitGaussians1d), intent(in)     ::      this
      integer, intent(in)                  ::      i
      integer, intent(in)                  ::      j

      real(kind=real64)   ::      dd

      dd = (i - this%mu(j))/this%sig(j)
      func1 = this%f(j)*this%f(j)*exp(-dd*dd/2)

      return
   end function func1

   pure integer function ndof(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns number of degrees of freedom used in the fitting
      !*       = 3* number of gaussians
      type(FitGaussians1d), intent(in)     ::      this
      ndof = 3*this%m
      return
   end function ndof

   pure real(kind=real64) function rss(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns residual sum of squares
      type(FitGaussians1d), intent(in)     ::      this
      integer             ::      ii
      real(kind=real64)   ::      dd
      rss = 0.0d0
      do ii = 1, this%n
         dd = func(this, ii) - this%in_dat(ii)
         rss = rss + dd*dd
      end do
      return
   end function rss

   pure subroutine drss(this, rr, dr)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns derivative of residual sum of squares wrt changes in parameters
      type(FitGaussians1d), intent(in)                     ::      this
      real(kind=real64), intent(out)                       ::      rr              !   = rss
      real(kind=real64), dimension(3*this%m), intent(out)   ::      dr
      integer             ::      ii, jj
      real(kind=real64)   ::      dd, dx, ff
      real(kind=real64), dimension(this%m)     ::      df, dm, ds, is
      rr = 0.0d0
      dr = 0.0d0
      do jj = 1, this%m
         is(jj) = 1/this%sig(jj)
      end do
      do ii = 1, this%n
         ff = 0.0d0
         do jj = 1, this%m
            dx = (ii - this%mu(jj))*is(jj)
            dd = exp(-dx*dx/2)
            ff = ff + this%f(jj)*this%f(jj)*dd

            df(jj) = 4*dd*this%f(jj)
            dm(jj) = 2*dd*this%f(jj)*this%f(jj)*dx*is(jj)
            ds(jj) = 2*dd*this%f(jj)*this%f(jj)*dx*dx*is(jj)

         end do
         dd = ff - this%in_dat(ii)
         rr = rr + dd*dd
         do jj = 1, this%m
            dr(3*jj - 2) = dr(3*jj - 2) + df(jj)*dd
            dr(3*jj - 1) = dr(3*jj - 1) + dm(jj)*dd
            dr(3*jj - 0) = dr(3*jj - 0) + ds(jj)*dd
         end do
      end do
      return
   end subroutine drss

   !---

   pure real(kind=real64) function AICc(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns the AIC (corrected) value assuming homoscedastic errors
      type(FitGaussians1d), intent(in)     ::      this

      integer             ::      kk
      real(kind=real64)   ::      rr
      real(kind=real64)   ::      reduced_chisquare_stat
      integer             ::      chisquare_dof
      real(kind=real64), parameter   ::      TWOPI = 6.283185307d0
      if (this%n == 0) then
         AICc = 0
      else
         rr = rss(this)
         kk = ndof(this) + 1                   !   +1 for unseen noise variable
         chisquare_dof = this%n - 1
         reduced_chisquare_stat = rr/chisquare_dof       !   rss per degree of freedom
         AICc = 2*kk + this%n*log(TWOPI*reduced_chisquare_stat) + chisquare_dof                    !   standard AIC ...
         if (this%n > kk + 1) AICc = AICc + 2*kk*(kk + 1)/(this%n - kk - 1)      !   ... with correction for small sample size
      end if

      return
   end function AICc

   pure subroutine fixIntegral(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      fix the integral under the curve
      type(FitGaussians1d), intent(inout)    ::      this
      real(kind=real64)           ::      intgrlIn
      real(kind=real64)           ::      intgrlOut
      integer                     ::      jj
      if (this%m == 0) return

      intgrlIn = sum(this%in_dat)
      if (intgrlIn < 1.0d-16) return

      !---    find the integral under the curve as sum over error functions
      intgrlOut = 0.0d0
      do jj = 1, this%m
         intgrlOut = intgrlOut + intgrl(this, jj)
      end do

      !---    fix the solution to have same integral
      if (intgrlOut > 1.0d-16) this%f = this%f*sqrt(intgrlIn/intgrlOut)

      return
   end subroutine fixIntegral

   pure real(kind=real64) function intgrl(this, j, x)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find the integral under curve j
      !*      int G[x] dx
      !*      note: am using integral ranges 0.5:n+0.5 so each point has equal weight
      !*      optionally integrate from 0.5 to x
      type(FitGaussians1d), intent(in)         ::      this
      integer, intent(In)                      ::      j
      real(kind=real64), intent(in), optional   ::      x
      real(kind=real64), parameter         ::      SQRTPION2 = 1.25331413731550d0      !    = sqrt( pi/2 )
      real(kind=real64), parameter         ::      ISQRT8 = 1/sqrt(8.0d0)

      real(kind=real64)           ::      is, xx

      is = 1/this%sig(j)

      if (present(x)) then
         xx = max(0.5d0, min(this%n + 0.5d0, x))
         intgrl = this%f(j)*this%f(j)*this%sig(j)*SQRTPION2*(ERF((2*xx - 2*this%mu(j))*is*ISQRT8) &
                                                             - ERF((1 - 2*this%mu(j))*is*ISQRT8))
      else
         intgrl = this%f(j)*this%f(j)*this%sig(j)*SQRTPION2*(ERF((1 + 2*this%n - 2*this%mu(j))*is*ISQRT8) &
                                                             - ERF((1 - 2*this%mu(j))*is*ISQRT8))
      end if

      return
   end function intgrl

   pure real(kind=real64) function xintgrl(this, j)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find the first moment of the integral under curve j
      !*      int x G[x] dx
      !*      note: am using integral ranges 0.5:n+0.5 so each point has equal weight
      type(FitGaussians1d), intent(in)     ::      this
      integer, intent(In)                  ::      j
      real(kind=real64), parameter         ::      SQRTPION2 = 1.25331413731550d0      !    = sqrt( pi/2 )
      real(kind=real64), parameter         ::      ISQRT8 = 1/sqrt(8.0d0)
      real(kind=real64)           ::      is

      is = 1/this%sig(j)
      xintgrl = this%f(j)*this%f(j)*this%sig(j)*( &
                this%sig(j)*(Exp(-((1 - 2*this%mu(j))*is)**2/8) &
                             - Exp(-((1 + 2*this%n - 2*this%mu(j))*is)**2/8)) &
                + SQRTPION2*this%mu(j)*(ERF((1 + 2*this%n - 2*this%mu(j))*is*ISQRT8) &
                                        - ERF((1 - 2*this%mu(j))*is*ISQRT8)) &
                )

      return
   end function xintgrl

   !---

   subroutine initialGuess(this, m, model)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      make an initial guess for the positions of m gaussians
      !*      model = 0 - improve on previous guess if possible by adding more evenly spaced gaussians
      !*      model = 1 - find maximum value, fix narrow width
      !*      model = 2 - find mean & std dev
      !*      model = 3 - evenly separated peaks
      !*      model = 4 - evenly separated peaks, very narrow

      type(FitGaussians1d), intent(inout)      ::      this
      integer, intent(in)                      ::      m
      integer, intent(in)                      ::      model
      real(kind=real64), dimension(:), pointer  ::      in_one
      type(FitGaussians1d)                    ::      one
      integer                                 ::      ii, jj
      real(kind=real64)                       ::      xx, ss
      one = FitGaussians1d_ctor(this%in_dat, 1)
      in_one => one%in_dat

      if (model == 0) then

         !   use the previous fit in this, and add a couple of new gaussians
         do ii = 1, this%n
            in_one(ii) = in_one(ii) - func(this, ii)
         end do
         ss = this%N/max(4, 2*this%m)
         do jj = this%m + 1, m
            xx = (jj - this%m)*this%N*1.0d0/(m - this%m)

            call fitOne(one, 3, xx, ss)
            this%f(jj) = one%f(1)
            this%mu(jj) = one%mu(1)
            this%sig(jj) = one%sig(1)
            do ii = 1, this%n
               in_one(ii) = in_one(ii) - func(one, ii)
            end do
         end do
         this%m = m

      else

         !---    fit single gaussians one at a time with no prior guess, and fit them
         this%m = m
         ss = this%N/(2*this%m)
         if (model == 4) ss = this%N/(4*this%m)
         if (model == 5) ss = this%N/(8*this%m)
         if (model == 6) ss = this%N/(4*this%m)
         do jj = 1, this%m
            xx = jj*this%N*1.0d0/(this%m + 1)
            if ((model == 6) .and. (this%m == 1)) then
               xx = this%N*0.5d0
            else if (model == 6) then
               xx = (jj - 1)*this%N*1.0d0/(this%m - 1)
            end if

            call fitOne(one, model, xx, ss)
            this%f(jj) = one%f(1)
            this%mu(jj) = one%mu(1)
            this%sig(jj) = one%sig(1)
            do ii = 1, this%n
               in_one(ii) = in_one(ii) - func(one, ii)
            end do
         end do

      end if
      call delete(one)
      return
   end subroutine initialGuess

   subroutine fitOne(this, model, x, s)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      fit the first gaussian ignoring all the others.
      type(FitGaussians1d), intent(inout)      ::      this
      integer, intent(in)                      ::      model
      real(kind=real64), intent(in), optional   ::      x, s         !   hints for mean and std dev, used for models 3+

      real(kind=real64), parameter     ::      PI = 3.141592654d0
      real(kind=real64)               ::      ss, mm, vv, dd
      integer                         ::      ii

      real(kind=real64), dimension(3)  ::      qq, drdq
      real(kind=real64)               ::      x1, x2, x3, y1, y2, y3
      real(kind=real64)               ::      rr, oldrr, moddrdq2
      type(GoldenSection)             ::      gold
      logical                         ::      isConverged, isWithinTimeLimit
      integer, parameter               ::      NLOOPS = 100
      real(kind=real64), parameter     ::      GTOL = 1.0d-10
      real(kind=real64), parameter     ::      FTOL = 1.0d-6
      integer                         ::      loop

      this%m = 1

      select case (model)

      case (1)

         this%mu(1) = maxloc(this%in_dat, dim=1)    !   position of maximum
         this%sig(1) = this%N/8                      !   small width
         this%f(1) = sqrt(max(0.0d0, maxval(this%in_dat)))             !   height tbd...

      case (2)

         !---    find data sum, mean and variance
         ss = 0.0d0; mm = 0.0d0; vv = 0.0d0
         do ii = 1, this%n
            dd = this%in_dat(ii)
            ss = ss + dd
            mm = mm + ii*dd
            vv = vv + ii*ii*dd
         end do

         !---    quick escape for no contribution.
         if (ss < 1.0d-16) then
            this%mu(1) = this%n/2
            this%sig(1) = this%n/4
            this%f(1) = 0.0d0
            return
         end if

         !---    fit the gaussian sum, mean, and variance.
         !       this is much much easier to do if we imagine the gaussian tails are small.
         !       so I'll use this as an initial fit. Note that we refine the fit below.

         mm = mm/ss
         vv = vv/ss

         dd = (vv - mm*mm)
         if (dd <= 1.0d-16) then
            this%mu(1) = maxloc(this%in_dat, dim=1)    !   position of maximum
            this%sig(1) = this%N/8                      !   small width
         else
            this%mu(1) = mm
            this%sig(1) = sqrt(dd)/2       !    half width, gives a better guess
         end if
         this%f(1) = 1.0d0

      case (3:)

         !---    fixed spacing and fixed width
         this%mu(1) = x
         this%sig(1) = s
         ii = max(1, min(this%N, nint(this%mu(1))))
         this%f(1) = sqrt(max(0.0d0, this%in_dat(ii)))
      end select
      call fixIntegral(this)

      !---    refine the fit using a golden section search
      oldrr = huge(1.0)
      GOLD_TOL = GTOL
      do loop = 1, NLOOPS

         qq(1:3) = (/this%f(1), this%mu(1), this%sig(1)/)
         call drss(this, rr, drdq)                   !   find the first derivative

         if (abs(rr - oldrr) < FTOL) exit           !   not improving the estimate = done

         x3 = 0.0d0; y3 = rr
         moddrdq2 = (drdq(1)*drdq(1) + drdq(2)*drdq(2) + drdq(3)*drdq(3))
         if (moddrdq2 < 1.0d-16) exit                      !   no derivative = done
         dd = min(10.0d0, rr/moddrdq2)*0.001d0
         moddrdq2 = max(moddrdq2, 10*rr)
         do ii = 1, 10
            x1 = -dd; y1 = findRss(qq, x1, drdq)
            if (y1 > rr) exit
            dd = dd*2
         end do
         gold = GoldenSection_ctor(x1, x3, y1, y3)
         do
            call nextPoint(gold, x2)
            y2 = findRss(qq, x2, drdq)
            call minimise(gold, x2, y2, isConverged, isWithinTimeLimit)
            if (isConverged) exit
            if (.not. isWithinTimeLimit) exit
         end do

         this%f(1) = abs(qq(1) + x2*drdq(1))
         this%mu(1) = qq(2) + x2*drdq(2)
         this%sig(1) = qq(3) + x2*drdq(3)
         oldrr = rr
      end do

      return

   contains
      !---^^^^^^^^

      pure real(kind=real64) function findRss(qq, alpha, dq)
         !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         !*      quick form of finding residual sum of squares at a trial position without having to alter contents this%mu etc.
         real(kind=real64), dimension(3), intent(in)       ::      qq, dq
         real(kind=real64), intent(in)                    ::      alpha
         integer             ::      ii
         real(kind=real64)   ::      dd
         real(kind=real64)   ::      is, mm, ff
         findRss = 0.0d0
         dd = qq(1) + alpha*dq(1)
         ff = dd*dd
         mm = qq(2) + alpha*dq(2)
         dd = qq(3) + alpha*dq(3)
         is = 1/max(0.5d0, dd)
         do ii = 1, this%n
            dd = (ii - mm)*is
            dd = ff*exp(-dd*dd/2)
            dd = dd - this%in_dat(ii)
            findRss = findRss + dd*dd
         end do
         return
      end function findRss

   end subroutine fitOne

   subroutine fit0(this)
      !---^^^^^^^^^^^^^^^^^^^^^^
      !*      fit multiple gaussians - return best solution according to AIC
      type(FitGaussians1d), intent(inout)              ::      this

      real(kind=real64), dimension(3*this%mmax)       ::      qq_best
      real(kind=real64)                           ::      aic, aic_best
      real(kind=real64), dimension(this%mmax)      ::      area

      integer             ::      jj, mm, mm_best
      logical             ::      ok
      real(kind=real64)   ::      dd

      !---    make all the possible fits. Keep track of which is the best one.
      aic_best = huge(1.0)
      do mm = 1, this%mmax
         call fit2(this, mm)
         aic = AICc(this)
         print *, "fit step ", mm, " ng ", this%m, " aic ", aic
         if (aic < aic_best) then
            mm_best = this%m
            aic_best = aic
            do jj = 1, this%m
               qq_best(jj*3 - 2) = this%f(jj)
               qq_best(jj*3 - 1) = this%mu(jj)
               qq_best(jj*3) = this%sig(jj)
            end do
         else
            if (this%m > mm_best + 2) exit
         end if
      end do

      !---    recover the best solution. compute the integrated area under each curve
      this%m = mm_best
      do jj = 1, mm_best
         this%f(jj) = qq_best(jj*3 - 2)
         this%mu(jj) = qq_best(jj*3 - 1)
         this%sig(jj) = qq_best(jj*3)
         area(jj) = intgrl(this, jj)
      end do

      !---    bubble sort to get greatest area first
      do
         ok = .true.
         do jj = 1, mm_best - 1
            if (area(jj + 1) > area(jj)) then
               dd = area(jj + 1); area(jj + 1) = area(jj); area(jj) = dd
               dd = this%f(jj + 1); this%f(jj + 1) = this%f(jj); this%f(jj) = dd
               dd = this%mu(jj + 1); this%mu(jj + 1) = this%mu(jj); this%mu(jj) = dd
               dd = this%sig(jj + 1); this%sig(jj + 1) = this%sig(jj); this%sig(jj) = dd
               ok = .false.
            end if
         end do
         if (ok) exit
      end do

      return
   end subroutine fit0

   subroutine fit0a(in_dat, mmax, this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      fit gaussians - return all solutions
      real(kind=real64), dimension(:), intent(in)           ::      in_dat
      integer, intent(in)                                  ::      mmax
      type(FitGaussians1d), dimension(:), intent(out)       ::      this

      integer             ::      jj, mm, mm_best
      real(kind=real64), dimension(3*mmax)       ::      qq_best
      real(kind=real64)                           ::      aic, aic_best
      real(kind=real64), dimension(mmax, mmax)      ::      area

      logical             ::      ok
      real(kind=real64)   ::      dd

      !---    allocate necessary memory
      do mm = 1, mmax
         this(mm) = FitGaussians1d_ctor(in_dat, mmax)
      end do

      !---    find all possible solutions
      aic_best = huge(1.0)
      mm_best = 0
      do mm = 1, mmax

         if (mm_best > 0) then
            this(mm)%m = this(mm_best)%m
            this(mm)%f(:) = this(mm_best)%f(:)
            this(mm)%mu(:) = this(mm_best)%mu(:)
            this(mm)%sig(:) = this(mm_best)%sig(:)
         end if

         call fit2(this(mm), mm)
         aic = AICc(this(mm))
         print *, "Lib_FitGaussians1d::fit0a info - trial ", mm, " ng ", this(mm)%m, " aic ", aic
         if (aic < aic_best) then
            mm_best = this(mm)%m
            aic_best = aic
            do jj = 1, this(mm)%m
               qq_best(jj*3 - 2) = this(mm)%f(jj)
               qq_best(jj*3 - 1) = this(mm)%mu(jj)
               qq_best(jj*3) = this(mm)%sig(jj)
            end do
         end if

         !---    bubble sort to get greatest area first

         do jj = 1, this(mm)%m
            area(jj, mm) = intgrl(this(mm), jj)
         end do
         do
            ok = .true.
            do jj = 1, this(mm)%m - 1
               if (area(jj + 1, mm) > area(jj, mm)) then
                  dd = area(jj + 1, mm); area(jj + 1, mm) = area(jj, mm); area(jj, mm) = dd
                  dd = this(mm)%f(jj + 1); this(mm)%f(jj + 1) = this(mm)%f(jj); this(mm)%f(jj) = dd
                  dd = this(mm)%mu(jj + 1); this(mm)%mu(jj + 1) = this(mm)%mu(jj); this(mm)%mu(jj) = dd
                  dd = this(mm)%sig(jj + 1); this(mm)%sig(jj + 1) = this(mm)%sig(jj); this(mm)%sig(jj) = dd
                  ok = .false.
               end if
            end do
            if (ok) exit
         end do

      end do

      return
   end subroutine fit0a

   subroutine fit2(this, m)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      fit m gaussians
      type(FitGaussians1d), intent(inout)      ::      this
      integer, intent(in)                      ::      m

      integer             ::      model
      integer             ::      jj, mm
      real(kind=real64), dimension(3*this%mmax)       ::      qq_best
      real(kind=real64)               ::      rr, rr_best

      !---    find the best initial guess model
      rr_best = huge(1.0)
      do model = 0, 6

         call initialGuess(this, m, model)
         call fit1(this, m)

         rr = RSS(this)
         !print *," fit2 model ",model," m = ",m," rss = ",rr

         if (rr < rr_best) then
            rr_best = rr
            do jj = 1, this%m
               qq_best(jj*3 - 2) = this%f(jj)
               qq_best(jj*3 - 1) = this%mu(jj)
               qq_best(jj*3) = this%sig(jj)
            end do
         end if

      end do

      !---    recover the best solution, cutting out any unused gaussians
      mm = 0
      do jj = 1, this%m
         if (qq_best(jj*3 - 2) > 1.0d-4) then
            mm = mm + 1
            this%f(mm) = qq_best(jj*3 - 2)
            this%mu(mm) = qq_best(jj*3 - 1)
            this%sig(mm) = qq_best(jj*3)
         end if
      end do
      this%m = mm

      return
   end subroutine fit2

   subroutine fit1(this, m)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      fit m gaussians given an initial guess
      type(FitGaussians1d), intent(inout)      ::      this
      integer, intent(in)                      ::      m
      real(kind=real64), parameter             ::      PI = 3.141592654d0
      real(kind=real64)       ::      dd
      integer                 ::      ii, jj

      real(kind=real64), dimension(3*this%mmax)       ::      qq, drdq

      real(kind=real64)               ::      x1, x2, x3, y1, y2, y3
      real(kind=real64)               ::      rr, oldrr, moddrdq2
      type(GoldenSection)             ::      gold
      logical                         ::      isConverged, isWithinTimeLimit
      integer, parameter               ::      NLOOPS = 100
      integer                         ::      loop
      real(kind=real64), parameter     ::      GTOL = 1.0d-10
      real(kind=real64), parameter     ::      FTOL = 1.0d-6

      !---    refine the fit using a golden section search
      oldrr = huge(1.0)
      do loop = 1, NLOOPS

         call drss(this, rr, drdq)
         do jj = 1, this%m
            qq(jj*3 - 2) = this%f(jj)
            qq(jj*3 - 1) = this%mu(jj)
            qq(jj*3) = this%sig(jj)
         end do

         if (abs(rr - oldrr) < FTOL) exit           !   not improving the estimate = done

         x3 = 0.0d0; y3 = findRss(qq, x3, drdq)
         moddrdq2 = dot_product(drdq, drdq)
         if (moddrdq2 < 1.0d-16) exit                      !   no derivative = done
         dd = min(10.0d0, rr/moddrdq2)*0.001d0
         moddrdq2 = max(moddrdq2, 10*rr)
         do ii = 1, 10
            x1 = -dd; y1 = findRss(qq, x1, drdq)
            if (y1 > rr) exit
            dd = dd*2
         end do

         gold = GoldenSection_ctor(x1, x3, y1, y3)
         GOLD_TOL = GTOL
         do
            call nextPoint(gold, x2)
            y2 = findRss(qq, x2, drdq)
            call minimise(gold, x2, y2, isConverged, isWithinTimeLimit)
            if (isConverged) exit
            if (.not. isWithinTimeLimit) exit
         end do
         qq(1:3*this%m) = qq(1:3*this%m) + x2*drdq(1:3*this%m)

         do jj = 1, this%m
            this%f(jj) = qq(jj*3 - 2)
            this%mu(jj) = qq(jj*3 - 1)
            this%sig(jj) = qq(jj*3)
         end do

         oldrr = rr
      end do

      return

   contains
      !---^^^^^^^^

      pure real(kind=real64) function findRss(qq, alpha, dq)
         !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         !*      quickly find trial solution without changing this%mu etc
         real(kind=real64), dimension(:), intent(in)       ::      qq, dq
         real(kind=real64), intent(in)                    ::      alpha
         integer             ::      ii, jj
         real(kind=real64)   ::      dd, gg
         real(kind=real64), dimension(this%m)   ::      is, mm, ff

         do jj = 1, this%m
            dd = qq(3*jj - 2) + alpha*dq(3*jj - 2)
            ff(jj) = dd*dd
            mm(jj) = qq(3*jj - 1) + alpha*dq(3*jj - 1)
            dd = qq(3*jj) + alpha*dq(3*jj)
            is(jj) = 1/max(0.5d0, dd)
         end do
         findRss = 0.0d0
         do ii = 1, this%n
            gg = 0.0d0
            do jj = 1, this%m
               dd = (ii - mm(jj))*is(jj)
               dd = ff(jj)*exp(-dd*dd/2)
               gg = gg + dd
            end do
            dd = gg - this%in_dat(ii)
            findRss = findRss + dd*dd
         end do
         return
      end function findRss

   end subroutine fit1

end module Lib_FitGaussians1d

