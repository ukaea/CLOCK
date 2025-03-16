
module Lib_SimpleProgressBar
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_SimpleProgressBar from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      Very simple module which outputs to screen progress through a long operation
!*      usage:
!*          do i = 1,n
!*              call progressBar( i,n )
!*          end do
!*      don't try to use i<1 by mistake.
!*

   use iso_fortran_env
   private

   public      ::      progressBar

contains
!---^^^^^^^^

   subroutine progressBar(i, n)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      output 10%,20% etc as iteration counter i gets towards n
      integer, intent(in)          ::      i       !   1<=i<=n
      integer, intent(in)          ::      n       !   end point 100%

      integer             ::  non10, pct

      non10 = max(1, nint(n*0.1d0))

      pct = 0

      if (i == 1) write (*, fmt='(i3,a)', advance="no") pct, "% "

      if (mod(i, non10) == 0) then
         pct = nint(i/real(non10))*10
         if (pct == 100) then
            write (*, fmt='(i3,a)', advance="yes") pct, "% "
            return
         else
            write (*, fmt='(i3,a)', advance="no") pct, "% "
         end if
      end if

      if (i == n) then
         if (pct /= 100) then
            write (*, fmt='(i3,a)', advance="yes") 100, "% "
         else
            write (*, fmt='(a)', advance="yes")
         end if
      end if

      return
   end subroutine progressBar

end module Lib_SimpleProgressBar
