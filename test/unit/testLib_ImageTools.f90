!------A program that Ctest can call that can be used to test multiple functions/subroutines
!------in Lib_ImageTools. If these test fails then the user must rerun the failed
!------test verbosely to diagnose the problem.
program testLib_ImageTools
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    testLib_ImageTools from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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

   use iso_fortran_env
   use Lib_ColouredTerminal !needed
   use Lib_UtilsForTests
   use Lib_ImageTools
   implicit none

   logical             ::      ok
   integer             ::      correctcount
   integer             ::      noofTests = 9 !
   character(len=32)   ::      libName = "testLib_ImageTools"
   logical             ::      tempCheck

   real(kind=real64), dimension(32, 42)              ::      testImg1
   real(kind=real64), dimension(:, :), allocatable    ::      testBufferedImg
   integer                                         ::      intcheck1, intcheck2, intcheck3
   type(subImage)                                  ::      testSubimage1
   real(kind=real64), dimension(16, 16)              ::      testImg2
   type(subImageRGB)                               ::      testSubimageRGB1
   real(kind=real64), dimension(3, 16, 16)            ::      testImgRGB1
   real(kind=real64), dimension(0:5, 0:6)                ::      testImg3, testImg5
   type(subImageGrid)                              ::      subImgGridTest1, subImgGridTest2 !SubimageGrid_ctor
   integer                                         ::      ix, iy
   type(subImageGridRGB)                           ::      subImgGridTestRGB1, subImgGridTestRGB2
   real(kind=real64), dimension(16, 16)                ::      testImg4
   type(subImageRGB)                               ::      testSubimageRGB2

   correctcount = 0

   !test 1 - Can we find how many subimages an image can be divided into for a given path size?
   tempCheck = .false.

   testImg1 = 1
   intcheck1 = 0
   intcheck2 = 0

   call getSubImageGridDims(testImg1, 5, intcheck1, intcheck2)
   if (intcheck1 == 6 .and. intcheck2 == 8) tempCheck = .true.
   call announceSubTest(libName, "getSubImageGridDims", 1, noofTests, tempCheck, correctcount)

   !test 2 - Can we cut a buffered image?
   tempCheck = .false.

   call cutBuffered_img(testImg1, 5, 6, 10, 10, testBufferedImg, 2) !end up with 10 by 9 elements each value 1
   if (sum(testBufferedImg) == 90) tempCheck = .true. !shoudl sum to 90
   call announceSubTest(libName, "cutBuffered_img", 2, noofTests, tempCheck, correctcount)

   !test 3 - can we set up a subimage derived type?
   tempCheck = .false.
   testImg2 = 1
   testSubimage1 = SubImage_ctor(testImg2(9:, 9:), 16, 16)
   if ((sum(testSubimage1%subImage) + testSubimage1%extentX + testSubimage1%extentY) == 96) tempCheck = .true.
   call announceSubTest(libName, "SubImage_ctor", 3, noofTests, tempCheck, correctcount)

   !test 4 - can we set up a subimageRGB derived type?
   tempCheck = .false.
   testImgRGB1 = 1
   testSubimageRGB1 = SubImageRGB_ctor(testImgRGB1(:, 9:, 9:), 16, 16)
   if ((sum(testSubimageRGB1%subImage) + testSubimageRGB1%extentX + testSubimageRGB1%extentY) == 224) tempCheck = .true.
   call announceSubTest(libName, "SubImageRGB_ctor", 4, noofTests, tempCheck, correctcount)

   !test 5 -can we set find the extents a subImageGrid derived type?
   tempCheck = .false.
   allocate (subImgGridTest2%subImg(3, 3))
   call getSubImageGridExtents(3, 3, 3, 1, 0, subImgGridTest2)!test greyscale
   if ((sum(subImgGridTest2%subImg%extentX) + sum(subImgGridTest2%subImg%extentY)) == 93) tempCheck = .true.
   call announceSubTest(libName, "getSubImageGridExtents", 5, noofTests, tempCheck, correctcount)

   !

   !test 6 -can we set find the extents a subImageGrid derived type?
   tempCheck = .false.
   allocate (subImgGridTestRGB1%subImgRGB(3, 3))
   call getSubImageGridExtentsRGB(3, 3, 3, 1, 0, subImgGridTestRGB1)!test greyscale
   if ((sum(subImgGridTestRGB1%subImgRGB%extentX) + sum(subImgGridTestRGB1%subImgRGB%extentY)) == 93) tempCheck = .true.
   call announceSubTest(libName, "getSubImageGridExtentsRGB", 6, noofTests, tempCheck, correctcount)

   !test 7 -can we set up a subImageGrid derived type?
   tempCheck = .false.
   testImg3 = 1
   intcheck3 = 0
   subImgGridTest1 = SubimageGrid_ctor(2, 2, testImg3, 3)
   intcheck3 = subImgGridTest1%subImageX + subImgGridTest1%subImageY !sum to 4
   do ix = 1, 2
      do iy = 1, 2
         associate (grid => subImgGridTest1%subImg(ix, iy))!alias for clarity
            intcheck3 = intcheck3 + (sum(grid%subImage) + grid%extentX + grid%extentY)! all loops sum to 80
         end associate
      end do
   end do
   if (intcheck3 == 76) tempCheck = .true.
   call announceSubTest(libName, "SubImageRGB_ctor", 7, noofTests, tempCheck, correctcount)

   !SubImageRGB_annotation_ctor

   !test 8 - can we set up a subimageRGB derived type using a greyscale iamge for annotation?
   tempCheck = .false.
   testImg4 = 1
   testSubimageRGB2 = SubImageRGB_annotation_ctor(testImg4(9:, 9:), 16, 16)
   if ((sum(testSubimageRGB2%subImage) + testSubimageRGB2%extentX + testSubimageRGB2%extentY) == 224) tempCheck = .true.
   call announceSubTest(libName, "SubImageRGB_annotation_ctor", 8, noofTests, tempCheck, correctcount)

   !test 9 - can we set up a RGB sub image grid using a greyscale iamge for annotation?
   tempCheck = .false.
   !SubimageGridRGB_annotation_ctor
   testImg5 = 1
   intcheck3 = 0
   subImgGridTestRGB2 = SubimageGridRGB_annotation_ctor(2, 2, testImg5, 3)
   intcheck3 = subImgGridTestRGB2%subImageX + subImgGridTestRGB2%subImageY !sum to 4
   do ix = 1, 2
      do iy = 1, 2
         associate (grid => subImgGridTestRGB2%subImgRGB(ix, iy))!alias for clarity
            intcheck3 = intcheck3 + (sum(grid%subImage) + grid%extentX + grid%extentY)! all loops sum to 156
         end associate
      end do
   end do
   if (intcheck3 == 160) tempCheck = .true.
   call announceSubTest(libName, "SubimageGridRGB_annotation_ctor", 9, noofTests, tempCheck, correctcount)
   !-----------------------------------------------
   ok = haveAllSubTestsPassed(correctcount, noofTests)

   call announcePassOrFail(ok)

end program testLib_ImageTools
