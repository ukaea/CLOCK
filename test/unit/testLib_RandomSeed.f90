!   gfortran -ffree-line-length-256 src/Lib_RandomSeed.f90 src/testLib_RandomSeed.f90 -o Test/testLib_RandomSeed.exe

program testLib_RandomSeed
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    testLib_RandomSeed from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      test correct functioning of Lib_RandomSeed
!*
!*      successful operation
!*
!*           $ ./Test/testLib_RandomSeed.exe
!*            set seed to 12345
!*           seed returned        12345       12382       12419       12456       12493       12530       12567       12604
!*           first number in ran0() sequence       0.27458033
!*
!*            done
!*

   use Lib_RandomSeed
   use Lib_ColouredTerminal
   use iso_fortran_env
   implicit none

   integer                     ::      seed
   integer, dimension(8)        ::      seed_out
   real(kind=real32)           ::      rr

   character(len=256), dimension(2)            ::      output
   character(len=*), dimension(2), parameter   ::      output0 = (/"set seed to 12345                               ", &
                                                                   "first number in ran0() sequence       0.27458033"/)

   logical                     ::      ok
   integer                     ::      ii

   output(1) = "set seed to 12345"
   seed = 12345
   call init_random_seed(seed)

   call get_random_seed(seed_out)
   write (*, fmt='(a,8i12)') "seed returned ", seed_out

   rr = ran0(seed)
   write (output(2), fmt='(a,f16.8)') "first number in ran0() sequence ", rr

   !---    here is the simple test: does the output look like the stored output?
   ok = .true.
   do ii = 1, size(output)
      if (trim(cutSpaces(output(ii))) == trim(cutSpaces(output0(ii)))) then
         write (*, fmt='(a)') trim(cutSpaces(output(ii)))
      else
         ok = .false.
         write (*, fmt='(a)') trim(cutSpaces(output0(ii)))//"    "//colour(RED, trim(cutSpaces(output(ii))))
      end if
   end do

   !---    output the result "PASS" or "FAIL"
   if (ok) then
      print *, colour(LIGHT_GREEN, "PASS")
   else
      print *, colour(RED, "FAIL")
   end if

   print *, ""
   print *, "done"
   print *, ""

end program testLib_RandomSeed
