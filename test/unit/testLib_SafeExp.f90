!*      A simple test program to check working of safeExp
!*      
    program testLib_SafeExp
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    testLib_SafeExp from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
        use iso_fortran_env
        use Lib_SafeExp
        use Lib_ColouredTerminal
        implicit none



        real(kind=real64)       ::      x,calculated,expected
        logical                 ::      testok,ok



        ok = .true.

    !---    test exponential in the regular range
        calculated = safeExp(1.0d0)
        expected = exp(1.0d0)
        testok = ( abs( calculated/expected - 1.0d0 ) < 1.0d-12 )
        print *,"unit test exp( 1.0d0 ) = ",calculated," expected ",expected,testok
        ok = ok .and. testok


    !---    test exponential too low/too high
        print *,"safeExp ranges ",log( tiny(1.0d0) ),log( huge(1.0) )

        x = 2 * log( tiny(1.0d0) )
        calculated = safeExp(x)
        expected = tiny(1.0d0)
        testok = ( abs( calculated/expected - 1.0d0 ) < 1.0d-12 )
        print *,"unit test exp( ",x," ) = ",calculated," expected ",expected,testok
        ok = ok .and. testok

        x = 2 * log( huge(1.0) )
        calculated = safeExp(x)
        expected = huge(1.0)
        testok = ( abs( calculated/expected - 1.0d0 ) < 1.0d-12 )                     
        print *,"unit test exp( ",x," ) = ",calculated," expected ",expected,testok
        ok = ok .and. testok


        print *,""
        if (ok) then
            print *,colour(GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if

    end program testLib_SafeExp