
module Lib_Png
!---^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_Png from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      Fortran side of Lib_Greyscale.c
!*      which uses Lib_Png to read/write 16 bit greyscale png.

   use iso_fortran_env
   use iso_c_binding
   implicit none
   private

   !---    interface to C side

   interface
      function read_greyscale_png_file(file_name, width, height) bind(c, name='read_greyscale_png_file') result(img)
         use, intrinsic                          ::  iso_c_binding, only: c_ptr, c_int, c_char
         use                                     ::  iso_fortran_env
         character(kind=c_char)                  ::  file_name(*)
         integer(kind=c_int)                     ::  width, height
         type(c_ptr)                             ::  img
      end function read_greyscale_png_file

      subroutine write_greyscale_png_file(file_name, width, height, img) bind(c, name='write_greyscale_png_file')
         use, intrinsic                          ::  iso_c_binding, only: c_ptr, c_int, c_char
         character(kind=c_char)                  ::  file_name(*)
         integer(kind=c_int), intent(in), value    ::  width, height
         type(c_ptr), intent(in), value          ::  img
      end subroutine write_greyscale_png_file

      subroutine write_rgb_png_file(file_name, width, height, img) bind(c, name='write_rgb_png_file')
         use, intrinsic                          ::  iso_c_binding, only: c_ptr, c_int, c_char
         character(kind=c_char)                  ::  file_name(*)
         integer(kind=c_int), intent(in), value    ::  width, height
         type(c_ptr), intent(in), value          ::  img
      end subroutine write_rgb_png_file

      subroutine destroy_storage(p) bind(c, name='destroy_png_storage')
         use, intrinsic                          ::  iso_c_binding, only: c_ptr
         type(c_ptr), intent(in), value          ::  p
      end subroutine destroy_storage

   end interface

#ifdef DEBUG
   logical, private, parameter         ::      LIB_PNG_DBG = .true.
#else
   logical, private, parameter         ::      LIB_PNG_DBG = .false.
#endif

   !---    interface to F side

   public      ::      read_greyscale_png
   public      ::      write_greyscale_png
   public      ::      write_rgb_png
   public      ::      removePngExtension
   public      ::      readPng
   public      ::      writePng
   public      ::      setClockPath

   interface readPng
      module procedure read_greyscale_png
      module procedure read_greyscale_png3d
   end interface

   interface writePng
      module procedure write_greyscale_png
      module procedure write_greyscale_png3d
   end interface

   interface setClockPath
      module procedure setClockPath_legacy
   end interface

   interface removePngExtension
      module procedure removePngExtension_legacy
   end interface

contains
!---^^^^^^^^

   subroutine read_greyscale_png(file_name, img, negative)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      character(len=*), intent(in)                                 ::      file_name
      real(kind=real64), dimension(:, :), allocatable, intent(out)    ::      img
      logical, intent(in), optional                                 ::      negative
      integer                                         ::      width, height
      type(c_ptr)                             ::      img_cptr
      real(kind=c_float), dimension(:), pointer ::      img_fptr
      character(len=len_trim(file_name) + 1)    ::      ff
      integer             ::      ii, jj
      logical             ::      neg

      ff = trim(file_name)//C_NULL_CHAR
      img_cptr = read_greyscale_png_file(ff, width, height)
      if (width*height == 0) then
         print *, "Lib_Png::read_greyscale_png() error - image read failed"
         return
      end if

      call c_f_pointer(img_cptr, img_fptr, [width*height])

      allocate (img(width, height))
      do jj = 1, height
         do ii = 1, width
            img(ii, jj) = img_fptr(ii + width*(jj - 1))
         end do
      end do

      call destroy_storage(img_cptr)

      neg = .false.; if (present(negative)) neg = negative
      if (neg) img = 1 - img

      return
   end subroutine read_greyscale_png

   subroutine write_greyscale_png(file_name, img, normalise, negative)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      character(len=*), intent(in)                                 ::      file_name
      real(kind=real64), dimension(:, :), intent(in)                 ::      img
      logical, intent(in), optional                                 ::      normalise
      logical, intent(in), optional                                 ::      negative

      type(c_ptr)                             ::      img_cptr
      real(kind=c_float), dimension(:), pointer ::      img_fptr
      integer                                 ::      width, height
      character(len=len_trim(file_name) + 5)    ::      ff
      integer             ::      ii, jj
      logical             ::      ok, norm, neg
      real(kind=real64)   ::      aa, bb, img_min, img_max
      real(kind=real64)   ::      ONE = 1.0d0 - 1.0d0/65536

      ff = trim(file_name)//".tmp"//C_NULL_CHAR
      width = size(img, dim=1)
      height = size(img, dim=2)

      norm = .false.; if (present(normalise)) norm = normalise
      neg = .false.; if (present(negative)) neg = negative
      aa = 1.0d0; bb = 0.0d0
      if (norm) then
         img_min = minval(img)
         img_max = maxval(img)
         if (img_max - img_min < 1.0d-16) then
            aa = 0.0d0; bb = 0.5d0
         else
            if (neg) then
               aa = ONE/(img_min - img_max)
               bb = -aa*img_max
            else
               aa = ONE/(img_max - img_min)
               bb = -aa*img_min
            end if
         end if
      else if (neg) then
         aa = -ONE; bb = 1.0d0
      end if

if (LIB_PNG_DBG) print *, "Lib_Png::write_greyscale_png() info - file,width,height = """//trim(file_name)//""",", width, ",", height
            if (LIB_PNG_DBG) print *,"Lib_Png::write_greyscale_png() info - min,max,avg (data)= ",minval(img),maxval(img),sum(img)/(width*height)

      allocate (img_fptr(width*height))
      do jj = 1, height
         do ii = 1, width
            img_fptr(ii + width*(jj - 1)) = real(aa*img(ii, jj) + bb, kind=c_float)
         end do
      end do

            if (LIB_PNG_DBG) print *,"Lib_Png::write_greyscale_png() info - min,max,avg (img) = ",minval(img_fptr),maxval(img_fptr),sum(img_fptr)/(width*height)
      img_cptr = c_loc(img_fptr(1))
      call write_greyscale_png_file(ff, width, height, img_cptr)
      deallocate (img_fptr)
      inquire (file=trim(file_name)//".tmp", exist=ok)
      if (.not. ok) then
         print *, "Lib_Png::write_greyscale_png() error - image write failed"
      else
         call system("mv "//trim(file_name)//".tmp "//trim(file_name))
      end if
      return
   end subroutine write_greyscale_png

   subroutine write_rgb_png(file_name, img)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      character(len=*), intent(in)                                 ::      file_name
      real(kind=real64), dimension(:, :, :), intent(in)               ::      img               !       [3,width,height]

      type(c_ptr)                             ::      img_cptr
      real(kind=c_float), dimension(:), pointer ::      img_fptr
      integer                                 ::      width, height
      character(len=len_trim(file_name) + 5)    ::      ff
      integer             ::      ii, jj
      logical             ::      ok

      ff = trim(file_name)//".tmp"//C_NULL_CHAR
      width = size(img, dim=2)
      height = size(img, dim=3)

      allocate (img_fptr(0:(3*width*height) - 1))
      do jj = 0, height - 1
         do ii = 0, width - 1
            img_fptr(3*ii + 0 + 3*width*jj) = real(img(1, ii + 1, jj + 1), kind=c_float)
            img_fptr(3*ii + 1 + 3*width*jj) = real(img(2, ii + 1, jj + 1), kind=c_float)
            img_fptr(3*ii + 2 + 3*width*jj) = real(img(3, ii + 1, jj + 1), kind=c_float)

         end do
      end do
      img_cptr = c_loc(img_fptr(0))
      call write_rgb_png_file(ff, width, height, img_cptr)
      deallocate (img_fptr)
      inquire (file=trim(file_name)//".tmp", exist=ok)
      if (.not. ok) then
         print *, "Lib_Png::write_rgb_png() error - image write failed"
      else
         call system("mv "//trim(file_name)//".tmp "//trim(file_name))
      end if
      return
   end subroutine write_rgb_png

   function removePngExtension_legacy(filename) result(prefix)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      character(len=*), intent(in)                         ::      filename
      character(len=len_trim(filename))                   ::      prefix
      integer             ::      ii

      ii = index(filename, ".png", back=.true.)
      if (ii == 0) then
         prefix = trim(filename)
      else
         prefix = filename(1:ii - 1)
      end if
      return
   end function removePngExtension_legacy

   subroutine setClockPath_legacy(path)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      character(len=*), intent(in), optional        ::      path
      if (LIB_PNG_DBG) print *, "Lib_Png::setClockPath() info - legacy call not needed "
      return
   end subroutine setClockPath_legacy

   subroutine read_greyscale_png3d(file_name, img, negative)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      character(len=*), intent(in)                                 ::      file_name
      real(kind=real64), dimension(:, :, :), allocatable, intent(out)  ::      img
      logical, intent(in), optional                                 ::      negative
      integer                                         ::      width, height, depth
      type(c_ptr)                             ::      img_cptr
      real(kind=c_float), dimension(:), pointer ::      img_fptr
      character(len=len_trim(file_name) + 1)    ::      ff
      integer             ::      ii, jj
      logical             ::      neg

      ff = trim(file_name)//C_NULL_CHAR
      img_cptr = read_greyscale_png_file(ff, width, height)
      depth = 1
      if (width*height == 0) then
         print *, "Lib_Png::read_greyscale_png() error - image read failed"
         return
      end if

      call c_f_pointer(img_cptr, img_fptr, [width*height])

      allocate (img(width, height, depth))
      do jj = 1, height
         do ii = 1, width
            img(ii, jj, 1) = img_fptr(ii + width*(jj - 1))
         end do
      end do

      call destroy_storage(img_cptr)

      neg = .false.; if (present(negative)) neg = negative
      if (neg) img = 1 - img

      return
   end subroutine read_greyscale_png3d

   subroutine write_greyscale_png3d(file_name, img, normalise, negative)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      character(len=*), intent(in)                                 ::      file_name
      real(kind=real64), dimension(:, :, :), intent(in)               ::      img
      logical, intent(in), optional                                 ::      normalise
      logical, intent(in), optional                                 ::      negative
      logical             ::      norm, neg
      norm = .false.; if (present(normalise)) norm = normalise
      neg = .false.; if (present(negative)) neg = negative
      call write_greyscale_png(file_name, img(:, :, 1), norm, neg)
      return
   end subroutine write_greyscale_png3d

end module Lib_Png
