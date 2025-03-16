
module Lib_ColouredTerminal
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_ColouredTerminal from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      Persuades fortran to output coloured text on the terminal.
!*      Note: can't guarantee that this will work with all compilers and all terminals.
!*
!*
!*
!*      Basic operation:
!*
!*          print *,"this is in "//colour(RED,"red")
!*

   implicit none
   private

   integer, public, parameter        ::      DEFAULT = 0
   integer, public, parameter        ::      DARK_GREY = 1
   integer, public, parameter        ::      PEACH = 2
   integer, public, parameter        ::      LIGHT_GREEN = 3
   integer, public, parameter        ::      LIGHT_YELLOW = 4
   integer, public, parameter        ::      LIGHT_BLUE = 5
   integer, public, parameter        ::      PINK = 6
   integer, public, parameter        ::      LIGHT_AQUA = 7
   integer, public, parameter        ::      PEARL_WHITE = 8
   integer, public, parameter        ::      BLACK = 9
   integer, public, parameter        ::      RED = 10
   integer, public, parameter        ::      GREEN = 11
   integer, public, parameter        ::      YELLOW = 12
   integer, public, parameter        ::      BLUE = 13
   integer, public, parameter        ::      PURPLE = 14
   integer, public, parameter        ::      AQUA = 15

   character(len=2), dimension(15), private, parameter &
      ::      COLOUR_CHARCODES = (/ &
         "90", "91", "92", "93", "94", "95", "96", "97", &
         "30", "31", "32", "33", "34", "35", "36"/)

   character(len=1), public, parameter   ::      NL = achar(10)
   character(len=1), public, parameter   ::      BS = achar(8)

   public          ::      colour
   public          ::      cutSpaces

contains
!---^^^^^^^^

   pure function colour(c, text) result(ctext)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      produce a text string "text" decorated with terminal markup so it appears in the requested colour.
      integer, intent(in)                  ::      c
      character(len=*), intent(in)         ::      text
      character(len=len_trim(text) + 9)    ::      ctext

      if ((c - 1)*(size(COLOUR_CHARCODES) - c) >= 0) then
         ctext = achar(27)//"["//COLOUR_CHARCODES(c)//"m"//trim(text)//achar(27)//"[0m"
      else
         ctext = trim(text)
      end if

      return
   end function colour

   pure function cutSpaces(text) result(ctext)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      remove leading and repeated spaces
      character(len=*), intent(in)         ::      text
      character(len=len(text))            ::      ctext
      integer                 ::      ii, jj
      logical                 ::      lastSpace
      character(len=1)        ::      letter
      lastSpace = .true.
      jj = 0
      ctext = ""
      do ii = 1, len_trim(text)
         letter = text(ii:ii)
         if (letter == " ") then
            if (lastSpace) cycle
            jj = jj + 1
            ctext(jj:jj) = letter
            lastSpace = .true.
         else
            jj = jj + 1
            ctext(jj:jj) = letter
            lastSpace = .false.
         end if
      end do
      return
   end function cutSpaces

end module Lib_ColouredTerminal

