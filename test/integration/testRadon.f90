
program test_Radon
!---^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    test_Radon from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      A simple program to test the working of the sinogram
   use Lib_Radon
   use Lib_Png
   use iso_fortran_env
   use Lib_ColouredTerminal
   implicit none

   character(len=256)                                  ::      filename
   real(kind=real64), dimension(:, :), allocatable        ::      img, img_out
   real(kind=real64), dimension(:, :), allocatable        ::      buffered_img
   real(kind=real64), dimension(:, :), allocatable        ::      sinogram

   integer         ::      pad = 2
   integer         ::      Nx, Ny
   integer         ::      nLine, globalNoofLines
   real(kind=real64), dimension(:, :), pointer            ::      line

   integer         ::      jx, jy, mx, my, kx, ky, lx, ly

   integer         ::      NPIXGRID = 100

   call get_command_argument(2, filename)
   if (len_trim(filename) > 0) read (filename, *) pad

   call get_command_argument(3, filename)
   if (len_trim(filename) > 0) read (filename, *) nPixGrid

   call get_command_argument(1, filename)
   if (len_trim(filename) == 0) stop "usage test_Radon <filename> [pad] [nPixGrid]"

   print *, "test_Radon"
   print *, "^^^^^^^^^^"
   print *, "using pad = ", pad
   print *, "file = """//trim(filename)//""""

   call readPng(filename, img)

   Nx = size(img, dim=1)
   Ny = size(img, dim=2)
   allocate (img_out(0:Nx - 1, 0:Ny - 1))
   img_out = 0
   print *, "image size ", Nx, "x", Ny
   globalNoofLines = 0

   !---    work on complete image
   call cutBuffered_img(img, 0, 0, Nx - 1, Ny - 1, buffered_img, pad)

   call writePng(trim(filename)//".bi.png", buffered_img)

   call construct_sinogram(buffered_img, sinogram, pad, pad)

   call writePng(trim(filename)//".sinogram.png", sinogram)

   call findLinesFromSinogram(sinogram, Nx, Ny, nLine, line, pad, pad)

   call reconstructImage(line(:, 1:nLine), img_out)
   call writePng(trim(filename)//".recon.png", img_out, normalise=.true.)
   deallocate (buffered_img)
   deallocate (sinogram)
   if (associated(line)) deallocate (line)

   !---    work on grid image
   kx = int(Nx/NPIXGRID)          !   number of grid blocks
   ky = int(Ny/NPIXGRID)
   mx = ceiling(Nx*1.0/kx)
   my = ceiling(Ny*1.0/ky)
   if (mx*(kx - 1) >= Nx - 2) kx = kx - 1
   if (my*(ky - 1) >= Ny - 2) ky = ky - 1
   mx = ceiling(Nx*1.0/kx)
   my = ceiling(Ny*1.0/ky)
   print *, "breaking into grid squares ", kx, "x", ky, " size(px) ", mx, "x", my

   img_out = 0
   do jx = 0, kx - 1
      img_out(jx*mx, :) = 0.2
   end do

   do jy = 0, ky - 1
      img_out(:, jy*my) = 0.2
   end do

   do jy = 0, ky - 1

      do jx = 0, kx - 1

         lx = min((jx + 1)*mx, Nx)            !   end pixel
         ly = min((jy + 1)*my, Ny)
         print *, "grid ", jx, jy, " from ", jx*mx, jy*my, " to ", lx - 1, ly - 1

         call cutBuffered_img(img, jx*mx, jy*my, lx - 1, ly - 1, buffered_img, pad)

         call construct_sinogram(buffered_img, sinogram, pad, pad)

         call findLinesFromSinogram(sinogram, lx - jx*mx, ly - jy*my, nLine, line, pad, pad)

         if (nLine > 0) then
            call reconstructImage(line(:, 1:nLine), img_out(jx*mx:lx - 1, jy*my:ly - 1))
            deallocate (line)
            globalNoofLines = globalNoofLines + nLine
         end if

         deallocate (buffered_img)
         deallocate (sinogram)

         print *, ""
      end do
   end do

   call writePng(trim(filename)//".recon_grid.png", img_out, normalise=.true.)

   print *, "total lines: ", globalNoofLines

   print *, ""
   print *, "done"
   print *, ""
   if (globalNoofLines == 6) then
      print *, colour(LIGHT_GREEN, "PASS")
   else
      print *, colour(RED, "FAIL")
   end if

   return
end program test_Radon

