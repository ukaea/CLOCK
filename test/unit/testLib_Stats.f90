!------Contains multiple functions/subroutine tests, if these test fails then the user must rerun the failed
!------test verbosely to diagnose the problem.
program testLib_Stats
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    testLib_Stats from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
!*    Copyright (C) 2024  James Heath

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
   !   Unit tests for Lib_Stats

   use iso_fortran_env
   use Lib_ColouredTerminal !needed
   use Lib_UtilsForTests
   use Lib_Stats
   implicit none

   logical             ::      ok
   integer             ::      correctcount
   integer             ::      noofTests = 3
   character(len=32)   ::      libName = "testLib_Stats"
   logical             ::      tempCheck

   real(kind=real64), dimension(10, 10)   :: arrayA, arrayB

   correctcount = 0

   !test 1
   arrayA = 1.0d0
   arrayB = 2.0d0
   tempCheck = .false.
   if (getRSS2Dreal(arrayA, arrayB) == 100.0d0) tempCheck = .true.
   call announceSubTest(libName, "getRSS2Dreal", 1, noofTests, tempCheck, correctcount)

   !test 2
   tempCheck = .false.
   if (compareAICvals(-9.0d0, -10.0d0) .eqv. .true.) tempCheck = .true.
   call announceSubTest(libName, "compareAICvals", 2, noofTests, tempCheck, correctcount)

   !test 3
   tempCheck = .false.
   if (calcAIC(41943.0d0, 600, 4194304) == -19314284.731772382d0) tempCheck = .true. !may need a floating point tolerance on this
   !print*,"DBG calcAIC is:",calcAIC(41943.0d0, 600,4194304)
   call announceSubTest(libName, "calcAIC", 3, noofTests, tempCheck, correctcount)

   ok = haveAllSubTestsPassed(correctcount, noofTests)

   call announcePassOrFail(ok)

end program testLib_Stats
