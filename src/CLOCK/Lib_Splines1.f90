
module Lib_Splines1
!---^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_Splines1 from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      construct a cubic Spline fit through 1-dimensional data with n knots
!*      to interpolate ask for a real(kind=real64) number between knot 0 and knot n-1
!*      ( or if periodic a number from 0 to n )
!*      note that the lengthscale is not stored, so you need to multiply derivatives correctly
!*
!*      construct a natural Spline through these with
!*          type(Spline1)                    ::  s
!*          s = Spline1_ctor( y )
!*      or construct a periodic Spline with
!*          s = Spline1_ctor( y,periodic=.true. )            !   assumes y(0) ( if it existed ) = y(N)
!*      or define gradients at both ends with
!*          s = Spline1_ctor( y,periodic=.false.,yp0,ypn )   !   where yp0 ~ y(1) - y(0) is the gradient at knot 0.

!*      find interpolations with
!*          yy = splint( s, x )
!*          dy = dsplint( s, x )

   use iso_fortran_env
   implicit none
   private

   public          ::      Spline1_ctor         ! Spline1_ctor( dat [,periodic] [,(yp0,ypn)/zero] )
   public          ::      delete
   public          ::      report              ! report(this [,unit] [,verbose])
   public          ::      splint
   public          ::      dsplint
   public          ::      d2splint
   public          ::      derivs
   public          ::      d3splint
   public          ::      findSpline          ! findSpline( this [,(yp0,ypn)/zero] )
   public          ::      getN
   public          ::      get, set
   public          ::      isPeriodic
!         public          ::      getCubic

   !---

   type, public     ::      Spline1
      private
      integer                                 ::      n
      logical                                 ::      periodic
      real(kind=real64), dimension(:), pointer  ::      y, y2
   end type Spline1

   !---
   interface delete
      module procedure delete0
   end interface

   interface report
      module procedure report0
      module procedure report2
   end interface

   interface splint
      module procedure splint0
   end interface

   interface dsplint
      module procedure dsplint0
   end interface

   interface d2splint
      module procedure d2splint0
   end interface

   interface derivs
      module procedure derivs0
   end interface

   interface d3splint
      module procedure d3splint0
   end interface

   interface Spline1_ctor
      module procedure Spline1_ctor0
      module procedure Spline1_ctor1
      module procedure Spline1_ctor2
      module procedure Spline1_ctor3
      module procedure Spline1_null
   end interface

   interface Spline1_alloc
      module procedure Spline1_alloc0
      module procedure Spline1_alloc1
   end interface

   interface findSpline
      module procedure findSpline0
      module procedure findSpline1
      module procedure findSpline2
   end interface

   interface getN
      module procedure getN0
   end interface

   interface get
      module procedure get0
      module procedure get1
   end interface

   interface set
      module procedure set0
      module procedure set1
   end interface

   interface isPeriodic
      module procedure isPeriodic0
   end interface

contains
!---^^^^^^^^

   function Spline1_null() result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      empty ctor
      type(Spline1)                        ::      this
      this%n = 0
      this%periodic = .false.
      nullify (this%y)
      nullify (this%y2)
      return
   end function Spline1_null

   function Spline1_alloc0(y) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      allocate storage but do not interpolate
      real(kind=real64), dimension(:), intent(in)         ::      y
      type(Spline1)                        ::      this
      this%n = size(y)
      allocate (this%y(0:this%n - 1))
      allocate (this%y2(0:this%n - 1))
      this%y(0:this%n - 1) = y(1:this%n)
      this%y2(0:this%n - 1) = 0.0d0
      return
   end function Spline1_alloc0

   function Spline1_alloc1(n) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      allocate storage but do not interpolate
      integer, intent(in)                  ::      n
      type(Spline1)                        ::      this
      this%n = n
      allocate (this%y(0:this%n - 1))
      allocate (this%y2(0:this%n - 1))
      this%y(0:this%n - 1) = 0.0d0
      this%y2(0:this%n - 1) = 0.0d0
      return
   end function Spline1_alloc1

   function Spline1_ctor0(y, periodic) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      default ctor
      real(kind=real64), dimension(:), intent(in)         ::      y
      logical, intent(in), optional         ::      periodic
      type(Spline1)                        ::      this
      this = Spline1_alloc(y)
      this%periodic = .false.
      if (present(periodic)) this%periodic = periodic

      !---    original line
      call findSpline(this, zero=.not. this%periodic)

      return
   end function Spline1_ctor0

   function Spline1_ctor1(y, yp0, ypn) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(:), intent(in)  ::      y
      real(kind=real64), intent(in)            ::      yp0, ypn
      type(Spline1)                            ::      this
      this = Spline1_alloc(y)
      this%periodic = .false.
      call findSpline(this, yp0, ypn)
      return
   end function Spline1_ctor1

   function Spline1_ctor2(y, periodic, yp0, ypn) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(:), intent(in)      ::      y
      logical, intent(in)                          ::      periodic
      real(kind=real64), intent(in)                ::      yp0, ypn
      type(Spline1)                                ::      this
      this = Spline1_alloc(y)
      this%periodic = periodic
      call findSpline(this, yp0, ypn)
      return
   end function Spline1_ctor2

   function Spline1_ctor3(n, periodic) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      integer, intent(in)                  ::      n
      logical, intent(in), optional         ::      periodic
      type(Spline1)                        ::      this
      this = Spline1_alloc(n)
      this%periodic = .false.; if (present(periodic)) this%periodic = periodic
      return
   end function Spline1_ctor3

   !---

   subroutine delete0(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^
      !*      deallocate dynamic memory.
      type(Spline1), intent(inout)          ::      this
      if (this%n == 0) return
      deallocate (this%y)
      deallocate (this%y2)
      this = Spline1_null()
      return
   end subroutine delete0

   !---

   subroutine report0(this, u, o)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      simple report of Spline1
      type(Spline1), intent(in)         ::      this
      integer, intent(in), optional     ::      u
      integer, intent(in), optional     ::      o
      integer             ::      uu, oo
      oo = 0; if (present(o)) oo = o
      uu = 6; if (present(u)) uu = u
      if (this%periodic) then
         write (unit=uu, fmt='(a,i6,a)') repeat(" ", oo)//"PeriodicSpline1[knots = ", this%n, "]"
      else
         write (unit=uu, fmt='(a,i6,a)') repeat(" ", oo)//"Spline1[knots = ", this%n, "]"
      end if
      return
   end subroutine report0

   subroutine report2(this, u, verbose, o)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      verbose report of Spline1 to unit u
      type(Spline1), intent(in)     ::      this
      integer, intent(in)          ::      u
      logical, intent(in)          ::      verbose
      integer, intent(in), optional ::      o
      integer             ::      ii, oo
      oo = 0
      if (present(o)) oo = o
      call report0(this, u, oo)
      if (verbose) then
         write (unit=u, fmt='(a)') repeat(" ", oo)//"y data"
         do ii = 0, this%n - 1
            write (unit=u, fmt='(a,i6,500g12.4)') repeat(" ", oo), ii, this%y(ii), this%y2(ii)
         end do
      end if
      return
   end subroutine report2

   function splint0(this, x) result(y)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given a real point x  on the line between 0 and n-1 ( or n for periodic )
      !*      return the Spline1 interpolated value
      !*      If x is outside the range, extrapolate but assume wrong answer
      type(Spline1), intent(in)                         ::      this
      real(kind=real64), intent(in)                    ::      x
      real(kind=real64)                               ::      y
      integer                 ::      ii, jj
      real(kind=real64)       ::      aa, bb, cc, dd, xx

      if (this%periodic) then

         ii = floor(x)         !   note x = -5.4 gives i = -6
         aa = x - ii             !   ... and so b = 0.6
         bb = 1.0d0 - aa

         ii = mod(ii + 32768*this%n, this%n)                        !   this then brings i back in bounds

         jj = mod(ii + 1, this%n)
         y = aa*this%y(jj) + bb*this%y(ii) &
             + ((aa*aa - 1)*aa*this%y2(jj) + (bb*bb - 1)*bb*this%y2(ii))/6.0

      else

         if (x < 0) then
            !---    need to extrapolate. Expect the wrong answer...
            !   write y = aa + bb x + 1/2 cc x^2 + 1/6 dd x^3
            aa = this%y(0)
            bb = this%y(1) - this%y(0) - (this%y2(1) + 2*this%y2(0))/6
            cc = this%y2(0)
            dd = this%y2(1) - this%y2(0)
            xx = x
            y = aa + xx*(bb + xx*(cc/2 + xx*dd/6))

            return
         else if (x >= this%n - 1) then
            !---    need to extrapolate. Expect the wrong answer...
            !   write y = aa + bb x + 1/2 cc x^2 + 1/6 dd x^3
            aa = this%y(this%n - 1)
            bb = this%y(this%n - 1) - this%y(this%n - 2) + (2*this%y2(this%n - 1) + this%y2(this%n - 2))/6
            cc = this%y2(this%n - 1)
            dd = this%y2(this%n - 1) - this%y2(this%n - 2)

            xx = x - (this%n - 1)
            y = aa + xx*(bb + xx*(cc/2 + xx*dd/6))

            return
         else

            ii = floor(x)         !    note that ii<this%n-1
            aa = x - ii
            bb = 1.0d0 - aa

            y = aa*this%y(ii + 1) + bb*this%y(ii) &
                + ((aa*aa - 1)*aa*this%y2(ii + 1) + (bb*bb - 1)*bb*this%y2(ii))/6.0

         end if
      end if
      return
   end function splint0

   !---

   function dsplint0(this, x) result(dy)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(Spline1), intent(in)             ::      this
      real(kind=real64), intent(in)        ::      x
      real(kind=real64)                   ::      dy
      integer                 ::      ii, jj
      real(kind=real64)       ::      aa, bb, cc, dd, xx

      if (this%periodic) then

         ii = floor(x)         !   note x = -5.4 gives i = -6
         aa = x - ii             !   ... and so b = 0.6
         bb = 1.0d0 - aa

         ii = mod(ii + 32768*this%n, this%n)                        !   this then brings i back in bounds

         jj = mod(ii + 1, this%n)
         dy = this%y(jj) - this%y(ii) &
              + ((3*aa*aa - 1)*this%y2(jj) - (3*bb*bb - 1)*this%y2(ii))/6.0

      else

         if (x < 0) then
            !---    need to extrapolate. Expect the wrong answer...
            !   write y = aa + bb x + 1/2 cc x^2 + 1/6 dd x^3
            bb = this%y(1) - this%y(0) - (this%y2(1) + 2*this%y2(0))/6
            cc = this%y2(0)
            dd = this%y2(1) - this%y2(0)
            xx = x
            dy = bb + xx*(cc + xx*dd/2)

            return
         else if (x >= this%n - 1) then
            !---    need to extrapolate. Expect the wrong answer...
            !   write y = aa + bb x + 1/2 cc x^2 + 1/6 dd x^3
            bb = this%y(this%n - 1) - this%y(this%n - 2) + (2*this%y2(this%n - 1) + this%y2(this%n - 2))/6
            cc = this%y2(this%n - 1)
            dd = this%y2(this%n - 1) - this%y2(this%n - 2)

            xx = x - (this%n - 1)
            dy = bb + xx*(cc + xx*dd/2)

            return
         else

            ii = floor(x)         !    note that ii<this%n-1
            aa = x - ii
            bb = 1.0d0 - aa

            dy = this%y(ii + 1) - this%y(ii) &
                 + ((3*aa*aa - 1)*this%y2(ii + 1) - (3*bb*bb - 1)*this%y2(ii))/6.0

         end if
      end if
      return
   end function dsplint0

   !---

   pure function d2splint0(this, x) result(d2y)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(Spline1), intent(in)             ::      this
      real(kind=real64), intent(in)        ::      x
      real(kind=real64)                   ::      d2y
      integer                 ::      ii, jj
      real(kind=real64)       ::      aa, bb, cc, dd, xx

      if (this%periodic) then

         ii = floor(x)         !   note x = -5.4 gives i = -6
         aa = x - ii             !   ... and so b = 0.6
         bb = 1.0d0 - aa

         ii = mod(ii + 32768*this%n, this%n)                        !   this then brings i back in bounds

         jj = mod(ii + 1, this%n)
         d2y = aa*this%y2(jj) + bb*this%y2(ii)

      else

         if (x < 0) then
            !---    need to extrapolate. Expect the wrong answer...
            !   write y = aa + bb x + 1/2 cc x^2 + 1/6 dd x^3
            cc = this%y2(0)
            dd = this%y2(1) - this%y2(0)
            xx = x
            d2y = cc + xx*dd

            return
         else if (x >= this%n - 1) then
            !---    need to extrapolate. Expect the wrong answer...
            !   write y = aa + bb x + 1/2 cc x^2 + 1/6 dd x^3
            cc = this%y2(this%n - 1)
            dd = this%y2(this%n - 1) - this%y2(this%n - 2)

            xx = x - (this%n - 1)
            d2y = cc + xx*dd

            return
         else

            ii = floor(x)         !    note that ii<this%n-1
            aa = x - ii
            bb = 1.0d0 - aa

            d2y = aa*this%y2(ii + 1) + bb*this%y2(ii)

         end if
      end if
      return
   end function d2splint0

   subroutine derivs0(this, x, y, dy, d2y)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given a real point x  on the line between 0 and n-1 ( or n for periodic )
      !*      return the Spline1 interpolated value
      !*      If x is outside the range, extrapolate but assume wrong answer
      type(Spline1), intent(in)                        ::      this
      real(kind=real64), intent(in)                    ::      x
      real(kind=real64), intent(out)                   ::      y, dy, d2y
      integer                 ::      ii, jj
      real(kind=real64)       ::      aa, bb, cc, dd, xx

      if (this%periodic) then

         ii = floor(x)         !   note x = -5.4 gives i = -6
         aa = x - ii             !   ... and so b = 0.6
         bb = 1.0d0 - aa

         ii = mod(ii + 32768*this%n, this%n)                        !   this then brings i back in bounds
         jj = mod(ii + 1, this%n)
         y = aa*this%y(jj) + bb*this%y(ii) &
             + ((aa*aa - 1)*aa*this%y2(jj) + (bb*bb - 1)*bb*this%y2(ii))/6.0
         dy = this%y(jj) - this%y(ii) &
              + ((3*aa*aa - 1)*this%y2(jj) - (3*bb*bb - 1)*this%y2(ii))/6.0

         d2y = aa*this%y2(jj) + bb*this%y2(ii)
      else

         if (x < 0) then
            !---    need to extrapolate. Expect the wrong answer...
            !   write y = aa + bb x + 1/2 cc x^2 + 1/6 dd x^3
            aa = this%y(0)
            bb = this%y(1) - this%y(0) - (this%y2(1) + 2*this%y2(0))/6
            cc = this%y2(0)
            dd = this%y2(1) - this%y2(0)
            xx = x
            y = aa + xx*(bb + xx*(cc/2 + xx*dd/6))
            dy = bb + xx*(cc + xx*dd/2)
            d2y = cc + xx*dd
            return
         else if (x >= this%n - 1) then
            !---    need to extrapolate. Expect the wrong answer...
            !   write y = aa + bb x + 1/2 cc x^2 + 1/6 dd x^3
            aa = this%y(this%n - 1)
            bb = this%y(this%n - 1) - this%y(this%n - 2) + (2*this%y2(this%n - 1) + this%y2(this%n - 2))/6
            cc = this%y2(this%n - 1)
            dd = this%y2(this%n - 1) - this%y2(this%n - 2)

            xx = x - (this%n - 1)
            y = aa + xx*(bb + xx*(cc/2 + xx*dd/6))
            dy = bb + xx*(cc + xx*dd/2)
            d2y = cc + xx*dd
            return
         else

            ii = floor(x)         !    note that ii<this%n-1
            aa = x - ii
            bb = 1.0d0 - aa

            y = aa*this%y(ii + 1) + bb*this%y(ii) &
                + ((aa*aa - 1)*aa*this%y2(ii + 1) + (bb*bb - 1)*bb*this%y2(ii))/6.0

            dy = this%y(ii + 1) - this%y(ii) &
                 + ((3*aa*aa - 1)*this%y2(ii + 1) - (3*bb*bb - 1)*this%y2(ii))/6.0
            d2y = aa*this%y2(ii + 1) + bb*this%y2(ii)

         end if
      end if
      return
   end subroutine derivs0

   !---

   pure function d3splint0(this, x) result(d3y)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(Spline1), intent(in)             ::      this
      real(kind=real64), intent(in)        ::      x
      real(kind=real64)                   ::      d3y
      integer                 ::      ii, jj

      if (this%periodic) then

         ii = floor(x)         !   note x = -5.4 gives i = -6

         ii = mod(ii + 32768*this%n, this%n)                        !   this then brings i back in bounds

         jj = mod(ii + 1, this%n)
         d3y = this%y2(jj) - this%y2(ii)

      else

         if (x < 0) then
            !---    need to extrapolate. Expect the wrong answer...
            !   write y = aa + bb x + 1/2 cc x^2 + 1/6 dd x^3
            d3y = this%y2(1) - this%y2(0)

            return
         else if (x >= this%n - 1) then
            !---    need to extrapolate. Expect the wrong answer...
            !   write y = aa + bb x + 1/2 cc x^2 + 1/6 dd x^3
            d3y = this%y2(this%n - 1) - this%y2(this%n - 2)

            return
         else

            ii = floor(x)         !    note that ii<this%n-1

            d3y = this%y2(ii + 1) - this%y2(ii)

         end if
      end if
      return
   end function d3splint0

   !---

   subroutine findSpline0(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(Spline1), intent(inout)              ::      this

      integer                                 ::      ii
      real(kind=real64), dimension(0:this%n - 1)     ::      uu
      real(kind=real64)                       ::      pp

      if (this%periodic) then
         call findSpline2(this, zero=.false.)
         return
      end if

      this%y2(0) = 0.0
      uu(0) = 0.0

      do ii = 1, this%n - 2
         pp = 1/(this%y2(ii - 1)/2 + 2)
         this%y2(ii) = -pp/2
         uu(ii) = (3*(this%y(ii + 1) - 2*this%y(ii) + this%y(ii - 1)) - uu(ii - 1)/2)*pp
      end do

      this%y2(this%n - 1) = 0.0

      do ii = this%n - 2, 0, -1
         this%y2(ii) = this%y2(ii)*this%y2(ii + 1) + uu(ii)
      end do
      return
   end subroutine findSpline0

   subroutine findSpline1(this, yp0, ypn)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(Spline1), intent(inout)              ::      this
      real(kind=real64), intent(in)            ::      yp0
      real(kind=real64), intent(in)            ::      ypn

      integer                                     ::      ii
      real(kind=real64), dimension(0:this%n - 1)     ::      uu
      real(kind=real64)                           ::      pp, un

      if (this%periodic) then
         if (abs(yp0 - ypn) > 1.0d-8) &
            stop "Lib_Splines1::findSpline1 error - trying to set different derivatives for ends of periodic function"

         un = 3*(this%y(0) - this%y(this%n - 1) - yp0)
         pp = 4/7.0d0
         this%y2(0) = -pp/2
         uu(0) = (3*(this%y(1) - 2*this%y(0) + this%y(this%n - 1)) - un/2)*pp

         do ii = 1, this%n - 2
            pp = 1/(this%y2(ii - 1)/2 + 2)
            this%y2(ii) = -pp/2
            uu(ii) = (3*(this%y(ii + 1) - 2*this%y(ii) + this%y(ii - 1)) - uu(ii - 1)/2)*pp
         end do

         un = 3*(ypn - this%y(this%n - 1) + this%y(this%n - 2))
         this%y2(this%n - 1) = (un - uu(this%n - 2)/2)/(this%y2(this%n - 2)/2 + 1)

         do ii = this%n - 2, 0, -1
            this%y2(ii) = this%y2(ii)*this%y2(ii + 1) + uu(ii)
         end do

      else    !   not periodic

         this%y2(0) = -0.5d0
         uu(0) = 3*((this%y(1) - this%y(0)) - yp0)      !   y2-y1 = yp dx + ypp dx^2/2 so uu1 = 6 ypp

         do ii = 1, this%n - 2
            pp = 1/(this%y2(ii - 1)/2 + 2)
            this%y2(ii) = -pp/2
            uu(ii) = (3*(this%y(ii + 1) - 2*this%y(ii) + this%y(ii - 1)) - uu(ii - 1)/2)*pp
         end do

         un = 3*(ypn - this%y(this%n - 1) + this%y(this%n - 2))
         this%y2(this%n - 1) = (un - uu(this%n - 2)/2)/(this%y2(this%n - 2)/2 + 1)

         do ii = this%n - 2, 0, -1
            this%y2(ii) = this%y2(ii)*this%y2(ii + 1) + uu(ii)
         end do

      end if
      return
   end subroutine findSpline1

   subroutine findSpline2(this, zero)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      fit a Spline1 with or without zero second derivatives at each end
      type(Spline1), intent(inout)              ::      this
      logical, intent(in)                      ::      zero

      integer                                 ::      ii, jj
      real(kind=real64), dimension(0:this%n - 1) ::      uu
      real(kind=real64)                       ::      pp, un

      integer, save                            ::      lastN = 0
      type(Spline1), save                       ::      pSpline1

      if (this%periodic) then

         if (zero) then

            un = (this%y(this%n - 2) - 8*this%y(this%n - 1) + 8*this%y(0) - this%y(1))/12
            un = 3*(this%y(0) - this%y(this%n - 1) - un)
            pp = 4/7.0d0
            this%y2(0) = -pp/2
            uu(0) = (3*((this%y(1) - 2*this%y(0)) + this%y(this%n - 1)) - un/2)*pp

            do ii = 1, this%n - 2
               pp = 1/(this%y2(ii - 1)/2 + 2)
               this%y2(ii) = -pp/2
               uu(ii) = (3*(this%y(ii + 1) - 2*this%y(ii) + this%y(ii - 1)) - uu(ii - 1)/2)*pp
            end do

            un = (this%y(this%n - 2) - 8*this%y(this%n - 1) + 8*this%y(0) - this%y(1))/12
            un = 3*(un - this%y(this%n - 1) + this%y(this%n - 2))
            this%y2(this%n - 1) = (un - uu(this%n - 2)/2)/(this%y2(this%n - 2)/2 + 1)

            do ii = this%n - 2, 0, -1
               this%y2(ii) = this%y2(ii)*this%y2(ii + 1) + uu(ii)
            end do

         else

            if (this%n > lastN) then
               if (lastN > 0) call delete(pSpline1)
               lastN = this%n
               pSpline1 = Spline1_alloc(2*this%n)
               pSpline1%periodic = .false.
            end if
            do ii = 0, 2*this%n - 1
               jj = mod(ii + this%n/2, this%n)
               pSpline1%y(ii) = this%y(jj)
            end do
            call findSpline(pSpline1)
            do ii = 0, this%n - 1
               jj = ii + (this%n + 1)/2
               this%y(ii) = pSpline1%y(jj)
               this%y2(ii) = pSpline1%y2(jj)
            end do

         end if

      else

         call findSpline(this)

      end if

   end subroutine findSpline2

!-------

   pure function getN0(this) result(n)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns the number of knots
      type(Spline1), intent(in)         ::      this
      integer                         ::      n
      n = this%n
      return
   end function getN0

   pure function isPeriodic0(this) result(is)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns true if the data is designated periodic
      type(Spline1), intent(in)         ::      this
      logical                         ::      is
      is = this%periodic
      return
   end function isPeriodic0

   function get0(this) result(y)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns a pointer to the y data
      type(Spline1), intent(inout)                  ::      this
      real(kind=real64), dimension(:), pointer      ::      y
      y => this%y
      return
   end function get0

   function get1(this, i) result(f)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return stored function value at a knot point - no interpolation
      type(Spline1), intent(in)                        ::      this
      integer, intent(in)                              ::      i
      real(kind=real64)                               ::      f
      f = this%y(i)
      return
   end function get1

   pure subroutine set0(this, y)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      sets the y data, does not refit
      type(Spline1), intent(inout)                  ::      this
      real(kind=real64), dimension(0:), intent(in)  ::      y
      this%y(0:this%n - 1) = y(0:this%n - 1)
      return
   end subroutine set0

   pure subroutine set1(this, i, y)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      sets the y data at knot i, does not refit
      type(Spline1), intent(inout)                  ::      this
      integer, intent(in)                          ::      i
      real(kind=real64), intent(in)                ::      y
      this%y(i) = y
      return
   end subroutine set1

end module Lib_Splines1
