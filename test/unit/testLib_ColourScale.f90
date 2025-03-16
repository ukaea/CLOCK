program testLib_ColourScale
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    testLib_ColourScale from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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

   use Lib_ColouredTerminal !needed
   use Lib_ColourScale
   use Lib_UtilsForTests
   implicit none

   logical             ::      ok
   integer             ::      correctcount
   integer             ::      noofTests = 6 !
   character(len=32)   ::      libName = "testLib_ColourScale"
   logical             ::      tempCheck

   correctcount = 0

   call listAvailableColourScales()

   !test 1
   tempCheck = .false.
   if (getColourScale("viridis   ") == 3) tempCheck = .true.
   call announceSubTest(libName, "getColourScale", 1, noofTests, tempCheck, correctcount)

   !test 2
   tempCheck = .false.
   if (getColourScaleName(5) == "bgr") tempCheck = .true.
   call announceSubTest(libName, "getColourScaleName", 2, noofTests, tempCheck, correctcount)

   !test 3
   tempCheck = .false. !getRGB0 interfaced through getRGB
   if (getRGB(0, 1.0d0) == 16777215) tempCheck = .true.
   !print'(a,I0)',"DBG getRGB",getRGB(0,1.0d0)
   call announceSubTest(libName, "getRGB0", 3, noofTests, tempCheck, correctcount)

   !test 4
   tempCheck = .false. !IF(ALL(abs(ARRAY1 - ARRAY2) < tolerance)))
   if (all(abs(getRGB_double(0, 1.0d0) - (/0.99609375d0, 0.99609375d0, 0.99609375d0/)) == 0d0)) tempCheck = .true. !shouldn't need a floating point tolerance for such
   !print*,"DBG getRGB_double",getRGB_double(0,1.0d0)                                                    !a trivial example
   call announceSubTest(libName, "getRGB_double", 4, noofTests, tempCheck, correctcount)

   !test 5
   tempCheck = .false. ! shadeColour( rgb,f,fc )
   if (all(abs(shadeColour((/0.5d0, 0.5d0, 0.5d0/), 0.25d0, 0.5d0) - (/0.25d0, 0.25d0, 0.25d0/)) == 0d0)) tempCheck = .true.
   !print*,"DBG shadeColour",shadeColour( (/0.5d0,0.5d0,0.5d0/),0.25d0,0.5d0 )
   call announceSubTest(libName, "shadeColour", 5, noofTests, tempCheck, correctcount)

   !test 6
   tempCheck = .false. ! transparentColour( rgb_back,rgb_fore,x )
if(all(abs(transparentColour( (/1.0d0,1.0d0,1.0d0/),(/1.0d0,1.0d0,1.0d0/),0.5d0 ) - (/1.0d0,1.0d0,1.0d0/)) == 0d0)) tempCheck=.true.
   !print*,"DBG shadeColour",transparentColour( (/1.0d0,1.0d0,1.0d0/),(/1.0d0,1.0d0,1.0d0/),0.5d0 )
   call announceSubTest(libName, "transparentColour", 6, noofTests, tempCheck, correctcount)

   ok = haveAllSubTestsPassed(correctcount, noofTests)

   call announcePassOrFail(ok)

end program testLib_ColourScale
