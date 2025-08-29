!*  A module containing a single simple function, an exponential that won't overflow
!*      y = safeExp(x) 
!*                          = exp(x)        if x is within bounds
!*                          = 0             if x too negative
!*                          = huge(1.0)     if x too positive       ( not huge(1.0d0) so that we can find safe_exp(x1) + safe_exp(x2)
!*  
!*  Set PASS_INVALID_ARG = .true. if you want to see what happens when an invalid argument is passed to exponential - could help debugging?
!*
!*  version history
!*      0.0.1       March 2025      First working version
!*

    module Lib_SafeExp
!---^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_SafeExp from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
!*    (c) UKAEA March 2025  Daniel Mason

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
        implicit none
        private

        real(kind=real64),private,parameter             ::      MIN_ARG = log( tiny(1.0d0) )
        real(kind=real64),private,parameter             ::      MAX_ARG = log( real( huge(1.0),kind=real64 ) )
        
        public          ::      safeExp


#ifdef DEBUG        
        logical,private,parameter   ::      PASS_INVALID_ARG = .false.
        integer,private             ::      static_underflow_count = 0
        integer,private             ::      static_overflow_count = 0
#endif


    contains
!---^^^^^^^^

#ifdef DEBUG

    real(kind=real64) function safeExp(x)
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        real(kind=real64),intent(in)            ::      x
        real(kind=real64)                       ::      safeArg

        safeArg = x
        if (x<MIN_ARG) then
            static_underflow_count = static_underflow_count + 1
            if (static_underflow_count == 1) then
                print *,"Lib_SafeExp::safeExp warning - arg = ",x
            else 
                if (any( (/10,100,1000,10000/) == static_underflow_count )) print *,"Lib_SafeExp::safeExp WARNING - arg = ",x,", underflow count = ",static_underflow_count
            end if
            if (.not. PASS_INVALID_ARG) safeArg = MIN_ARG
        else if (x>MAX_ARG) then
            static_overflow_count = static_overflow_count + 1
            if (static_overflow_count == 1) then
                print *,"Lib_SafeExp::safeExp warning - arg = ",x
            else 
                if (any( (/10,100,1000,10000/) == static_overflow_count )) print *,"Lib_SafeExp::safeExp warning - arg = ",x,", overflow count = ",static_overflow_count
            end if
            if (.not. PASS_INVALID_ARG) safeArg = MAX_ARG
        else if (x/=x) then
            if (.not. PASS_INVALID_ARG) stop "Lib_SafeExp::safeExp ERROR - arg = NaN"
        end if
            
        safeExp = exp( safeArg )

        return
    end function safeExp

#else

        elemental real(kind=real64) function safeExp(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),intent(in)            ::      x
            real(kind=real64)                       ::      safeArg

            safeArg = max( MIN_ARG, min( MAX_ARG,x ) )

            safeExp = exp( safeArg )

            return
        end function safeExp

#endif


    end module Lib_SafeExp
        