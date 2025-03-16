module Lib_UtilsForTests
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_UtilsForTests from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      Contains helper functions/subroutines for testing
   use iso_fortran_env
   use Lib_ColouredTerminal
   implicit none
   private
   public          ::      haveAllSubTestsPassed
   public          ::      announceSubTest
   public          ::      announcePassOrFail

contains

   pure function haveAllSubTestsPassed(countIn, noofTestsIn) result(ok)
      !Generalised logic counting the right number of sub tests have passed
      integer, intent(in)              ::      countIn, noofTestsIn
      logical                         ::      ok
      if (countIn == noofTestsIn) then
         ok = .true.
      else
         ok = .false.
      end if
      return
   end function haveAllSubTestsPassed

   subroutine announceSubTest(libnameIn, functOrSubName, testId, noofTestsIn, condition, countInOut)
      !quick subroutine to print results of subtest
      character(len=*), intent(in)    ::          libnameIn, functOrSubName
      integer, intent(in)              ::          testId, noofTestsIn
      logical, intent(in)              ::          condition
      integer, intent(inout)          ::          countInOut

      if (condition) then
         countInOut = countInOut + 1
 print'(a,i0,a,i0,a,a,a,a,a)', "Test ", testId, " of ", noofTestsIn, " in ", trim(libnameIn), ", ", trim(functOrSubName), " SUCCESS"
      else
            print'(a,i0,a,i0,a,a,a,a,a,a)',"Test ", testId, " of ", noofTestsIn, " in ", trim(libnameIn), ", ", trim(functOrSubName)," ",colour(RED,"FAIL")
      end if

   end subroutine announceSubTest

   subroutine announcePassOrFail(ok)
      !generalisation of important PASS/FAIL logic
      logical, intent(in)      ::      ok

      if (ok) then
         print *, colour(LIGHT_GREEN, "PASS")
      else
         print *, colour(RED, "FAIL")
      end if

   end subroutine announcePassOrFail

end module Lib_UtilsForTests
