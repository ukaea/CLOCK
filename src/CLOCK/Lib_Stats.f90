module Lib_Stats
   !---^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_Stats from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
   !*      Module containg various statistics functions for use in CLOCK

   use iso_fortran_env
   implicit none
   private

   public      ::      getRSS2Dreal
   public      ::      compareAICvals
   public      ::      calcAIC

contains
   !---^^^^^^^^
   !PURE FUNCTIONS

   pure logical function compareAICvals(AicOld, AicNew) result(isNewBetter)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Input two AIC values: (AicOld) from the origional model (AicNew) from a new model with reduced degrees of freedom
      !If the new model has the same AIC value or less it is an improvement and is selected.
      real(kind=real64), intent(in)                 ::      AicOld, AicNew
      if (AicNew <= AicOld) then
         isnewbetter = .true.
      else
         isnewbetter = .false.
      end if

      return
   end function compareAICvals

   pure real(kind=real64) function calcAIC(RSS, DoFs, noofObs) result(AICc)!do I need noof dimensions for model?
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Return Aikike information criterion (AIC) given the resisudal sum of squares (RSS), number of degrees of freedom
      !(DoFs) and the number of observations (noofObs). Includes a correction for a small sample size if the DoFs>noofObs
      Real(kind=real64), intent(in)                 ::      RSS
      integer, intent(in)                 ::      DoFs, noofObs
      integer                     ::      correctedDoFs !corrected for noise
      integer                     ::      chisquare_dof !corrected for size constraint.
      Real(kind=real64)           ::      reduced_chisquare_stat !RSS per degree of freedom

      if (DoFs == 0) then
         AICc = 0
      else
         correctedDoFs = DoFs + 1                   !   +1 for unseen noise variable
         chisquare_dof = noofObs - 1                   !   -1 for constraint of a fixed number of observations, 2D so should this be -2
         reduced_chisquare_stat = RSS/chisquare_dof       !   rss per chiquare degree of freedom
         !AICc = 2*correctedDoFs + noofObs*log( TWOPI * reduced_chisquare_stat ) + chisquare_dof !old DM definition
         AICc = 2*correctedDoFs + noofObs*log(reduced_chisquare_stat) !wiki definition

         if (noofObs < correctedDoFs) AICc = AICc + 2*correctedDoFs*(correctedDoFs + 1)/(noofObs - correctedDoFs - 1)      !   ... with correction for small sample size
         !was (noofObs>kk+1), wrong way round and why +1 here?
      end if
      return
   end function calcAIC

   !FUNCTIONS

   real(kind=real64) function getRSS2Dreal(array1, array2) result(RSS)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Calculates the residual sume of squares (RSS) between two 2D arrays of reals, the arrays must be the same size
      real(kind=real64), dimension(:, :), intent(in)     ::      array1, array2
      integer                                         ::      length1d1, length2d1 !size of boths arrays first dimesion
      integer                                         ::      length1d2, length2d2 !size of boths arrays first dimesion
      logical                                         ::      ok !can arrays 1 and 2 be compared
      integer                                         ::      ii, jj !array indices
      real(kind=real64)                               ::      dif !difference between each element pair

      ok = .false.
      RSS = 0.0d0

      length1d1 = size(array1, dim=1)!get all input array dimesions
      length2d1 = size(array2, dim=1)
      length1d2 = size(array1, dim=2)
      length2d2 = size(array2, dim=2)

      if ((length1d1 == length2d1) .and. (length1d2 == length2d2)) ok = .true.

      if (ok) then
         do ii = 1, length1d1
            do jj = 1, length1d2
               dif = array1(ii, jj) - array2(ii, jj)
               RSS = RSS + dif*dif
            end do
         end do
      else
         print *, "ERROR: Attempt to calcualte RSS arrays with different lengths Lib_Stats, getRSS2Dreal"
      end if

      return
   end function getRSS2Dreal

   !SUBROUTINES

end module Lib_Stats
