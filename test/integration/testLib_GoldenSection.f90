program testLib_GoldenSection
!---^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    testLib_GoldenSection from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
   use Lib_GoldenSection
   use Lib_ColouredTerminal
   use iso_fortran_env
   implicit none

   integer                 ::      ii, passcount
   real(kind=real64)       ::      xx

   real(kind=real64)       ::      xt, yt, x1, x3
   type(GoldenSection)     ::      gold
   logical     ::          is, ok

   !---    demonstrate presence of actual minimum

   do ii = 1, 10000
      xx = ii*0.001d0
      if (func(xx) < min(func(xx - 1.0d-3), func(xx + 1.0d-3))) then
         print *, "minimum of func found at ", xx, func(xx)
      end if
   end do
   passcount = 0 !count how many passes

   !----    guess between 1 and 1.5
   x1 = 1.0d0; x3 = 1.5d0
   print *, "guess between ", x1, x3
   gold = GoldenSection_ctor(x1, x3, func(x1), func(x3))
   call report(gold)
   do
      call nextPoint(gold, xt)
      yt = func(xt)
      call minimise(gold, xt, yt, is, ok)
      call report(gold)
      if (is) then
         passcount = passcount + 1
         exit
      end if
      if (.not. ok) then
         print *, "testLib_Gold error - minimise failed"
         exit
      end if
   end do
   print *, "golden section result ", xt, yt
   print *, ""

   !----    guess between 1.5 and 2
   x1 = 1.5d0; x3 = 2.0d0
   print *, "guess between ", x1, x3
   gold = GoldenSection_ctor(x1, x3, func(x1), func(x3))
   call report(gold)
   do
      call nextPoint(gold, xt)
      yt = func(xt)
      call minimise(gold, xt, yt, is, ok)
      call report(gold)
      if (is) then
         passcount = passcount + 1
         exit
      end if
      if (.not. ok) then
         print *, "testLib_Gold error - minimise failed"
         exit
      end if
   end do
   print *, "golden section result ", xt, yt
   print *, ""

   !----    guess between 2 and 2.5
   x1 = 2.0d0; x3 = 2.5d0
   print *, "guess between ", x1, x3
   gold = GoldenSection_ctor(x1, x3, func(x1), func(x3))
   call report(gold)
   do
      call nextPoint(gold, xt)
      yt = func(xt)
      call minimise(gold, xt, yt, is, ok)
      call report(gold)
      if (is) then
         passcount = passcount + 1
         exit
      end if
      if (.not. ok) then
         print *, "testLib_Gold error - minimise failed"
         exit
      end if
   end do
   print *, "golden section result ", xt, yt
   print *, ""

   !----    guess between 3 and 4
   x1 = 3.0d0; x3 = 4.0d0
   print *, "guess between ", x1, x3
   gold = GoldenSection_ctor(x1, x3, func(x1), func(x3))
   call report(gold)
   do
      call nextPoint(gold, xt)
      yt = func(xt)
      call minimise(gold, xt, yt, is, ok)
      call report(gold)
      if (is) then
         passcount = passcount + 1
         exit
      end if
      if (.not. ok) then
         print *, "testLib_Gold error - minimise failed"
         exit
      end if
   end do
   print *, "golden section result ", xt, yt
   print *, ""

   !----    guess between 4 and 5
   x1 = 4.0d0; x3 = 5.0d0
   print *, "guess between ", x1, x3
   gold = GoldenSection_ctor(x1, x3, func(x1), func(x3))
   call report(gold)
   do
      call nextPoint(gold, xt)
      yt = func(xt)
      call minimise(gold, xt, yt, is, ok)
      call report(gold)
      if (is) then
         passcount = passcount + 1
         exit
      end if
      if (.not. ok) then
         print *, "testLib_Gold error - minimise failed"
         exit
      end if
   end do
   print *, "golden section result ", xt, yt
   print *, ""

   !----    guess between 5 and 6
   x1 = 5.0d0; x3 = 6.0d0
   print *, "guess between ", x1, x3
   gold = GoldenSection_ctor(x1, x3, func(x1), func(x3))
   call report(gold)
   do
      call nextPoint(gold, xt)
      yt = func(xt)
      call minimise(gold, xt, yt, is, ok)
      call report(gold)
      if (is) then
         passcount = passcount + 1
         exit
      end if
      if (.not. ok) then
         print *, "testLib_Gold error - minimise failed"
         exit
      end if
   end do
   print *, "golden section result ", xt, yt
   print *, ""

   !----    guess between 1 and 3
   x1 = 1.0d0; x3 = 3.0d0
   print *, "guess between ", x1, x3
   gold = GoldenSection_ctor(x1, x3, func(x1), func(x3))
   call report(gold)
   do
      call nextPoint(gold, xt)
      yt = func(xt)
      call minimise(gold, xt, yt, is, ok)
      call report(gold)
      if (is) then
         passcount = passcount + 1
         exit
      end if
      if (.not. ok) then
         print *, "testLib_Gold error - minimise failed"
         exit
      end if
   end do
   print *, "golden section result ", xt, yt
   print *, ""

   print *, ""
   print *, "done"
   print *, ""
   print *, "Pass count is = ", passcount

   if (passcount == 7) then
      print *, colour(LIGHT_GREEN, "PASS")
   else
      print *, colour(RED, "FAIL")
   end if

contains
!---^^^^^^^^

   pure real(kind=real64) function func(x)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      real(kind=real64), intent(in)        ::      x

      if (x < 0) then
         func = huge(1.0)
      else
         func = 1.0d0 + x*(1.0d0 + x*(-0.1d0 + x*0.2d0))*sin(3*x) + cos(x)*sqrt(x)
      end if
      return
   end function func
end program testLib_GoldenSection
