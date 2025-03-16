
module Lib_Filenames
!---^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_Filenames from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      short module to handle file numbering and suffixes
!*
   use iso_fortran_env
   implicit none
   private

   !---

   public          ::      getSuffix
   public          ::      removeSuffix
   public          ::      numberFile

   !--

   interface numberFile
      module procedure numberFile1
      module procedure numberFile2
   end interface

contains
!---^^^^^^^^

   function getSuffix(file_in) result(suffix)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return filename suffix - eg "directory/foo.bar.txt" becomes "txt"

      character(len=*), intent(in)     ::      file_in
      character(len=len(file_in))     ::      suffix

      integer     ::      ii
      ii = index(file_in, ".", back=.true.)

      suffix = ""
      if (ii > index(file_in, "/", back=.true.)) then                   !   last "." after last "/"
         suffix = file_in(ii + 1:)
      end if
      return
   end function getSuffix

   function removeSuffix(file_in) result(file_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return filename with suffix removed - eg "directory/foo.bar.txt" becomes "directory/foo.bar"

      character(len=*), intent(in)     ::      file_in
      character(len=len(file_in))     ::      file_out

      integer     ::      ii
      file_out = ""
      ii = index(file_in, ".", back=.true.)
      if (ii > 0) then
         file_out(1:ii - 1) = file_in(1:ii - 1)
      else
         file_out = file_in
      end if
      return
   end function removeSuffix

   function numberFile1(file_prefix, n) result(file_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return filename with number added
      !*      eg "directory/foo.bar" + [n=123] + "png" returns "directory/foo.bar.00123.png"
      !*      if suffix not supplied then it is automatically found
      !*      eg "directory/foo.bar" + [n=123]         returns "directory/foo.00123.bar"

      character(len=*), intent(in)             ::      file_prefix
      integer, intent(in)                      ::      n

      character(len=len(file_prefix) + 7)       ::      file_out
      character(len=5)                    ::      aaaa

      integer     ::      nn
      nn = mod(n + 900000, 100000)    !   ensure a number between 0:99999
      write (aaaa, fmt='(i5)') nn
      aaaa = adjustl(aaaa)
      aaaa = repeat("0", 5 - len_trim(aaaa))//trim(aaaa)

      file_out = getSuffix(file_prefix)
      if (len_trim(file_out) > 0) then
         file_out = trim(removeSuffix(file_prefix))//"."//aaaa//"."//trim(file_out)
      else
         file_out = trim(file_prefix)//"."//aaaa
      end if
      return

   end function numberFile1

   function numberFile2(file_prefix, n, file_suffix) result(file_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return filename with number added
      !*      eg "directory/foo.bar" + [n=123] + "png" returns "directory/foo.bar.00123.png"
      !*      if suffix not supplied then it is automatically found
      !*      eg "directory/foo.bar" + [n=123]         returns "directory/foo.00123.bar"

      character(len=*), intent(in)             ::      file_prefix
      integer, intent(in)                      ::      n
      character(len=*), intent(in)             ::      file_suffix

      character(len=len(file_prefix) + len(file_suffix) + 7)       ::      file_out
      character(len=5)                    ::      aaaa

      integer     ::      nn
      nn = mod(n, 100000)    !   ensure a number between 0:99999
      write (aaaa, fmt='(i5)') nn
      aaaa = adjustl(aaaa)
      aaaa = repeat("0", 5 - len_trim(aaaa))//trim(aaaa)

      if (len_trim(file_suffix) > 0) then
         file_out = trim(file_prefix)//"."//aaaa//"."//trim(file_suffix)
      else
         file_out = trim(file_prefix)//"."//aaaa
      end if

      return

   end function numberFile2

end module Lib_Filenames

