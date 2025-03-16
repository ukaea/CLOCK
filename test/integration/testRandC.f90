
program testRandC
!---^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    testRandC from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      unit test for Lib_RidlerCalvard

   use iso_fortran_env
   use Lib_ColouredTerminal
   use Lib_CommandLineArguments
   use Lib_RidlerCalvard
   implicit none

   real(kind=real64), parameter     ::      PI = 3.141592654d0

   type(CommandLineArguments)      ::      cla
   logical                         ::      quiet = .false.

   real(kind=real64), dimension(:, :), allocatable        ::      img
   integer                                             ::      Nx = 1024, Ny
   logical                                             ::      ok

   integer                         ::      ix, iy
   real(kind=real64)               ::      dd, i2s2, sig
   integer                         ::      ii
   real(kind=real64), dimension(0:255)      ::      hanal, himg
   real(kind=real64)               ::      xx, dx, rr, dr, dr1, dr2, ftot

   real(kind=real64)               ::      bb, tt, ff, bstd, fbar, snr, fot

   !---    check command line args
   cla = CommandLineArguments_ctor(10)
   call setProgramDescription(cla, "testRandC.exe")
   call setProgramVersion(cla, "1.0.0")
   call get(cla, "q", quiet, LIB_CLA_OPTIONAL, " quiet mode")
   call get(cla, "Nx", Nx, LIB_CLA_OPTIONAL, " image size for test ( note: need > 512 or won't have smooth histogram )")
   call report(cla)              !   gives full output iff there is an error, or "-h" is a cla
   if (.not. allRequiredArgumentsSet(cla)) stop
   if (hasHelpArgument(cla)) stop
   call delete(cla)

   !---    set up input image.
   Ny = Nx
   allocate (img(0:Nx - 1, 0:Ny - 1))

   !---    image is dark grey, with a single bright spot in the centre,
   !       so I can compute R&C analytically

   sig = Nx/8               !   size of gaussian radius
   i2s2 = 1/(2*sig*sig)      !   1/2 sigma^2
   do iy = 0, Ny - 1
      do ix = 0, Nx - 1
         dd = (ix - Nx/2)*(ix - Nx/2) + (iy - Ny/2)*(iy - Ny/2)
         dd = exp(-dd*i2s2)
         img(ix, iy) = dd
      end do
   end do

   !---    integral img dx dy = 2 pi r^2 Erf[ Nx/(2 r sqrt 2) ] Erf[ Ny/(2 r sqrt 2) ] = Nx Ny <I>
   !       integral img drho  = 2 pi r^2 ( 1 - Exp[ - rho^2/(2 r^2) ] )                = 2 pi rho^2 <I_fore(rho)>
   !

   dd = sum(img)
   ftot = (2*PI*sig*sig)*ERF(Nx/(sqrt(8.0d0)*sig))*ERF(Ny/(sqrt(8.0d0)*sig))
   if (.not. quiet) then
      print *, "average intensity calc ", dd/(Nx*Ny), " analytic ", ftot/(Nx*Ny)
   end if
   ok = (abs(dd/ftot - 1) < (Nx*Ny)*1.0d-3)

   !---    analytic histogram - look at highest 50% intensity to avoid corners of img
   !       hist( 0 ) should contain fraction of pixels between 0:1/256
   dx = 1/256.0d0
   do ii = 1, 255
      xx = ii*dx
      !   area between xx and xx+dx should be pi rho^2 - pi (rho+drho)^2
      !   with rho = sqrt[ 2 sig^2 log 1/x ]
      !   and drho = - sig^2 dx / ( rho x )
      rr = sig*sqrt(2*log(1/xx))
      dr = -sig*sig*dx/(rr*xx)                           !   this is the first order change in radius

      dd = 1 + 2*dx*(rr*rr - sig*sig)/(rr*rr*xx)
      dr1 = 1 + sqrt(dd)
      dr2 = 1 - sqrt(dd)
      dr1 = dr1*rr*sig*sig/(rr*rr - sig*sig)
      dr2 = dr2*rr*sig*sig/(rr*rr - sig*sig)              !   this is the second order change in radius. Better.
      dr = dr2
      hanal(ii) = -PI*(2*rr*dr + dr*dr)

   end do
   hanal(0) = max(0.0d0, Nx*Ny - sum(hanal(1:)))
   call findHist(img, himg)

   dd = 0.0d0
   do ii = 64, 255
      xx = himg(ii)*(Nx*Ny)/hanal(ii) - 1
      dd = max(dd, xx*xx)
   end do
   if (.not. quiet) print *, "histogram max frac err test ", dd

   !---    here is the test:
   !       due to Gauss circle theorem we expect a bound on the error, not zero error
   ok = ok .and. (dd*(Nx*Nx*Nx) < 1.4d7)

   !---    now find R and C image properties
   call findImageIntensityFeatures(img, bb, tt, ff, bstd, fbar, snr, fot)

   if (.not. quiet) print *, "bb,tt,ff,bstd,fbar,snr,fot ", bb, tt, ff, bstd, fbar, snr, fot

   !   now extract the expected radius(squared) corresponding to fraction over thresh
   rr = Nx*Ny*fot/PI

   !   and the threshold according to this value
   dd = exp(-rr*i2s2)

   if (.not. quiet) print *, "pixels over threshold test ", abs(tt/dd - 1)
   ok = ok .and. (abs(tt/dd - 1) < 1e-2)

   !   now find the foreground level and check intensities are sensible
   !rr = sig*sqrt( 2*log(1/ff) )                            !   radius for foreground

   dd = 2*PI*sig*sig*(1 - exp(-rr*i2s2))           !   sum intensity of pixels over thresh

   xx = (Nx*Ny - PI*rr)*bb                              !   expected sum intensity pixels over backgound

   ftot = (2*PI*sig*sig)*ERF(Nx/(sqrt(8.0d0)*sig))*ERF(Ny/(sqrt(8.0d0)*sig))   !   total sum intensity pixels

   if (.not. quiet) print *, "sum intensity test ", abs((dd + xx)/ftot - 1)

   ok = ok .and. (abs((dd + xx)/ftot - 1) < 2.5e-3)

   !---    output the result "PASS" or "FAIL"
   if (ok) then
      print *, colour(LIGHT_GREEN, "PASS")
   else
      print *, colour(RED, "FAIL")
   end if

   !---    bye bye
   if (.not. quiet) then
      print *, ""
      print *, "done"
      print *, ""
   end if

end program testRandC
