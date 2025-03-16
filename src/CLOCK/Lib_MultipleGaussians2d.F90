
module Lib_MultipleGaussians2d
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_MultipleGaussians2d from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      code which defines and fits a multiple 2d gaussians to input data
!*

   use Lib_GoldenSection
   use Lib_Gaussian2d
   use Lib_Quicksort
   use iso_fortran_env

   !---    for debugging
   use Lib_Png
   use Lib_Filenames

   implicit none
   private

   public          ::      MultipleGaussians2d_ctor
   public          ::      report
   public          ::      delete

   !---

   public          ::      fit
   public          ::      ignore
   public          ::      getRss
   public          ::      finddrss
   public          ::      findT
   public          ::      findQ
   public          ::      setDat
   public          ::      getDat
   public          ::      getN, setN
   public          ::      getG
   public          ::      getSigma
   public          ::      getTheta
   public          ::      getX
   public          ::      getY
   public          ::      getF
   public          ::      getQ
   public          ::      getT
   public          ::      getD
   public          ::      getW !JPH hack
   public          ::      translate
   public          ::      xscale          !   multiply length components x and D
   public          ::      add
   public          ::      setLib_MultipleGaussians2d_dbg
   public          ::      setIgnore!ignore_individual_guassian
   public          ::      ignore_individual_guassian
   ! public          ::      Cullspot    !gets rid of spot if eccentricity is greater than input

   !---

   real(kind=real64), private, parameter         ::      PI = 3.141592653590d0
   real(kind=real64), public                    ::      Q_THRESH = huge(1.0)        !   AIC per spot threshold
   real(kind=real64), public                    ::      T_THRESH = 0.0d0            !   t* threshold
   real(kind=real64), public                    ::      E_THRESH = 3.0d0            !   discard any gaussian with eccentricity over this
   real(kind=real64), public                    ::      DETECTFMIN = 0.0d0          !   discard any gaussian with intensity less than this

   logical, private                 ::      LIB_MG2D_DBG = .false.

   !---

   type, public     ::      MultipleGaussians2d
      private
      integer                                     ::      n           !   number of gaussians
      type(Gaussian2d), dimension(:), pointer       ::      g
      real(kind=real64), dimension(:), pointer      ::      w           !   weight under gaussian - set -ve for ignore
      real(kind=real64), dimension(:, :), pointer    ::      sigma       !   major/minor diameters
      real(kind=real64), dimension(:), pointer      ::      theta       !   angle to x-axis
      real(kind=real64), dimension(:), pointer      ::      t           !   t-test
      real(kind=real64), dimension(:), pointer      ::      q           !   quality
   end type

   !---

   interface MultipleGaussians2d_ctor
      module procedure MultipleGaussians2d_null
      module procedure MultipleGaussians2d_ctor0
      module procedure MultipleGaussians2d_ctor1
      module procedure MultipleGaussians2d_ctor2
   end interface

   interface delete
      module procedure delete0
   end interface

   interface report
      module procedure report0
   end interface

   interface fit
      module procedure fit0
   end interface

   interface loglikelihood
      module procedure loglikelihood0
   end interface

   interface getDat
      module procedure get0
      module procedure get1
   end interface

   interface getN
      module procedure getN0
   end interface

   interface setN
      module procedure setN0
   end interface

   interface getG
      module procedure getG0
   end interface

   interface setDat
      module procedure set0
      module procedure set1
   end interface

   interface getRss
      module procedure getRss0
   end interface

   interface findDrss
      module procedure findDrss0
   end interface

   interface getX
      module procedure getX0
   end interface

   interface getY
      module procedure getY0
   end interface

   interface getF
      module procedure getF0
   end interface

   interface getTheta
      module procedure getTheta0
   end interface

   interface getD
      module procedure getD0
   end interface

   interface getSigma
      module procedure getSigma0
   end interface

   interface getQ
      module procedure getQ0
   end interface

   interface getT
      module procedure getT0
   end interface

   interface getW
      module procedure getW0
   end interface

   interface translate
      module procedure translate0
      module procedure translate1
      module procedure translate2
      module procedure translate3
      module procedure translate4
      module procedure translate5
   end interface

   interface add
      module procedure add0
      module procedure add1
   end interface

   interface subtract
      module procedure subtract0
   end interface

   interface ignore
      module procedure ignore0
      module procedure ignore1
   end interface

   ! interface Cullspot
   !         module procedure Cullspot0
   ! end interface

contains
!---^^^^^^^^

   pure function MultipleGaussians2d_null() result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      returns a null gaussian
      type(MultipleGaussians2d)           ::      this
      this%n = 0
      nullify (this%g)
      nullify (this%w)
      nullify (this%sigma)
      nullify (this%theta)
      nullify (this%t)
      nullify (this%q)
      return
   end function MultipleGaussians2d_null

   pure function MultipleGaussians2d_ctor0(n) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      allocates memory for several gaussians
      integer, intent(in)                  ::      n
      type(MultipleGaussians2d)           ::      this
      integer         ::      ii
      this%n = n
      allocate (this%g(this%n))
      allocate (this%w(this%n))
      allocate (this%sigma(2, this%n))
      allocate (this%theta(this%n))
      allocate (this%t(this%n))
      allocate (this%q(this%n))
      do ii = 1, this%n
         this%g(ii) = Gaussian2d_ctor()
      end do
      this%w = 0
      this%sigma = 0
      this%theta = 0
      this%t = 0
      this%q = 0
      return
   end function MultipleGaussians2d_ctor0

   function MultipleGaussians2d_ctor1(g) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      construct a wrapper for an input set of gaussians
      type(Gaussian2d), dimension(:), intent(in)     ::      g
      type(MultipleGaussians2d)           ::      this
      integer         ::      ii
      this = MultipleGaussians2d_ctor0(size(g))
      do ii = 1, this%n
         this%g(ii) = Gaussian2d_ctor(getDat(g(ii)))
         this%w(ii) = weight(this%g(ii))
         call findSigmaAndAngle(this%g(ii), this%sigma(1, ii), this%sigma(2, ii), this%theta(ii))
      end do
      this%t = 0
      this%q = 0
      return
   end function MultipleGaussians2d_ctor1

   function MultipleGaussians2d_ctor2(dat) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      construct a wrapper for an input set of Gaussians defined by their raw data values
      !*      note: advanced call does not set sigma or weight.
      real(kind=real64), dimension(:, :), intent(in)     ::      dat
      type(MultipleGaussians2d)           ::      this
      integer         ::      ii
      this = MultipleGaussians2d_ctor0(size(dat, dim=2))
      do ii = 1, this%n
         this%g(ii) = Gaussian2d_ctor(dat(:, ii))
!                this%w(ii) = weight( this%g(ii) )
!                call findSigmaAndAngle( this%g(ii),this%sigma(1,ii),this%sigma(2,ii),this%theta(ii) )
      end do
      this%w = 0
      this%sigma = 0
      this%theta = 0
      this%t = 0
      this%q = 0
      return
   end function MultipleGaussians2d_ctor2

   !---

   subroutine delete0(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^
      !*      deallocates dynamic memory
      type(MultipleGaussians2d), intent(inout)    ::      this
      if (this%n == 0) return
      deallocate (this%g)
      deallocate (this%w)
      deallocate (this%sigma)
      deallocate (this%theta)
      deallocate (this%t)
      deallocate (this%q)
      this = MultipleGaussians2d_null()
      return
   end subroutine delete0

   !---

   subroutine report0(this, u, o)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      simple output. Defaults to unit 6 = screen.
      !*      optional Gaussian2d o determines left hand margin
      type(MultipleGaussians2d), intent(in)     ::      this
      integer, intent(in), optional     ::      u, o
      integer     ::      uu, oo
      integer     ::      ii

      uu = 6; if (present(u)) uu = u
      oo = 0; if (present(o)) oo = o
      write (unit=uu, fmt='(8(a,i6))') repeat(" ", oo)//"MultipleGaussians2d [n = ", count(this%w > 0), "/", this%n, "]"
      do ii = 1, this%n
         if (.not. ignore(this, ii)) then
write (unit=uu, fmt='(8(a,f10.3))', advance="no") repeat(" ", oo + 4)//"[w,t,q=", this%w(ii), ",", this%t(ii), ",", this%q(ii), "] "
            call report(this%g(ii), uu, 0)
         end if
      end do
      return
   end subroutine report0

   !---

   subroutine fit0(this, img, indx, imin, imax)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      fit multiple gaussians to the image
      !*      on input indx should show which pixels are (tentatively) associated with which peak )
      !*      imin is taken to be the background level
      !*      anything above imax is taken to be out-of-bounds
      !*      indx should be ordered so that spot i is more important than spot i+1
      type(MultipleGaussians2d), intent(inout)         ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img
      integer, dimension(0:, 0:), intent(in)             ::      indx
      real(kind=real64), intent(in)                    ::      imin, imax

      integer             ::      Nx, Ny
      integer             ::      ii, nGroups, nn, nPeaks

      real(kind=real64), dimension(:, :), allocatable    ::      img_tmp, img_dbg
      logical, dimension(:, :), allocatable              ::      overlapMatrix
      integer, dimension(:, :), allocatable              ::      group
      type(Gaussian2d)                                ::      gg

      real(kind=real64), dimension(:, :, :), allocatable  ::      dat_trial
      real(kind=real64), dimension(:), allocatable      ::      aic_trial
      integer, dimension(:), allocatable                ::      nspot_trial
      integer                                         ::      trial, nTrials
      integer                                         ::      nTrialsMax
      real(kind=real64)                               ::      aa

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      !---    temporarily relax thresholds for fitting
      E_THRESH = E_THRESH*2
      DETECTFMIN = DETECTFMIN/2

      !---    allocate memory for individual gaussians and make an initial fit to individual peaks
      nPeaks = maxval(indx)
      if (LIB_MG2D_DBG) then
         print *, "Lib_MultipleGaussians2d::fit0 info - assessing ", nPeaks, " peaks in ", Nx, "x", Ny, " px"
         allocate (img_dbg(0:Nx - 1, 0:Ny - 1))
      end if

      call delete(this)
      this = MultipleGaussians2d_ctor(nPeaks)
      allocate (img_tmp(0:Nx - 1, 0:Ny - 1))
      img_tmp = img
      nn = 0
      do ii = 1, this%n

         if (LIB_MG2D_DBG) then
            print *, "fitting peak ", ii, "/", this%n, " px ", count(indx == ii)
            call fit(gg, img_tmp, imin, imax, mask=(indx == ii), dbg=.true.)
         else
            call fit(gg, img_tmp, imin, imax, mask=(indx == ii))
         end if

         if (LIB_MG2D_DBG) then
            call report(gg)
            img_dbg = imin
            where (indx /= ii)
               img_dbg = LIB_G2D_IGNORE
            end where
            call add(img_dbg, gg)
            call writePng(trim(numberFile("output", ii))//".png", img_dbg)
            img_dbg = img_tmp
            where (indx /= ii)
               img_dbg = LIB_G2D_IGNORE
            end where
            call writePng(trim(numberFile("input", ii))//".png", img_dbg)
         end if

         if (.not. ignore(gg)) then
            nn = nn + 1
            call subtract(img_tmp, gg)
            img_tmp = max(0.0d0, img_tmp)
            this%g(nn) = Gaussian2d_ctor(getDat(gg))
            aa = area(gg)

            if (aa > 0) then
               this%w(nn) = weight(gg)*count(indx == ii)/aa
            else
               this%w(nn) = -1
            end if
         else if (LIB_MG2D_DBG) then
            print *, "Lib_MultipleGaussians2d::fit0 info - ignore peak ", ii
         end if

         if (LIB_MG2D_DBG) print *, ""

      end do

      !---    allocate memory
      this%n = nn             !   note that I could have discarded a few at this stage.
      !---    restore thresholds for fitting
      E_THRESH = E_THRESH/2
      DETECTFMIN = DETECTFMIN*2

      if (this%n == 0) return !   quick escape for nothing to do.

      allocate (overlapMatrix(this%n, this%n))
      allocate (group(0:this%n, this%n))        !   group(0,i) is number in the ith group. group(j,i) is the jth gaussian in the group.

      nTrialsMax = min(1000, this%n*this%n)
      allocate (dat_trial(0:6, this%n, nTrialsMax))
      allocate (aic_trial(nTrialsMax))
      allocate (nspot_trial(nTrialsMax))
      aic_trial = huge(1.0)                   !   if this trial hasn't been considered, then assume its a terrible fit

      call cleanImg(img, imin, imax, img_tmp)

      !---    now we need to decide what the best number of gaussians is.
      !       there could be a very large number of individual gaussians to fit here.
      !       fitting all simultaneously is hard, so try a divide-and-conquer strategy by fitting groups rather than all at once
      !       for this to be successful, first need to split the long list of possible gaussians into safe groups.
      !       I define a safe group as containing members with low overlap.

      if (LIB_MG2D_DBG) call report(this)

      trial = 1           !   trial 1 is individually fitted only
      !---    store the fit
      nn = 0
      do ii = 1, this%n
         if (this%w(ii) < 0) cycle
         nn = nn + 1
         dat_trial(0, nn, trial) = this%w(ii)
         dat_trial(1:6, nn, trial) = getDat(this, ii)
      end do
      nspot_trial(trial) = nn
      aic_trial(trial) = 2*(6*nn + 1) - 2*loglikelihood0(this, img_tmp)
      nTrials = 1

      do trial = 2, nTrialsMax             !   try different numbers of spots

         if (trial == 2) then
            !   only have the individually fitted spots to consider, so assume all weights are good and try and fit them all.
         else
            !   now have the fully fitted spots, so use the last trial as the starting guess...
            nn = nSpot_trial(trial - 1) - 1           !   less the lowest weighted spot
            if (nn == 0) then
               !   no spots left to fit. Quit
               nTrials = trial - 1
               exit
            end if
            call setDat(this, dat_trial(1:6, 1:nn, trial - 1))
            this%w(1:nn) = dat_trial(0, 1:nn, trial - 1)

         end if

         call groupSpots(this, Nx, Ny, group, nGroups)
         call fit1(this, img_tmp, group(:, 1:nGroups))

         !---    store the fit
         call orderByWeight(this)
         nn = 0
         do ii = 1, this%n
            if (this%w(ii) < 0) cycle
            nn = nn + 1
            dat_trial(0, nn, trial) = this%w(ii)
            dat_trial(1:6, nn, trial) = getDat(this, ii)
         end do
         nspot_trial(trial) = nn
         aic_trial(trial) = 2*(6*nn + 1) - 2*loglikelihood0(this, img_tmp)
         nTrials = trial

         !---    quick exit? Is the aic definitely increasing?

         if (trial > 3) then
            if ((aic_trial(trial) > aic_trial(trial - 1)) .and. (aic_trial(trial - 1) > aic_trial(trial - 2))) exit   !   for the last two turns, the aic has not improved by removing a spot.
         end if

      end do

      if (LIB_MG2D_DBG) then
         do trial = 1, nTrials
            print *, "trial ", trial, " nSpots remaining ", nSpot_trial(trial), "/", nPeaks, " aic ", aic_trial(trial)
         end do
      end if

      !---    return the best solution
      trial = minloc(aic_trial(1:nTrials), dim=1)
      this%n = nSpot_trial(trial)

      do ii = 1, this%n
         call setDat(this, ii, dat_trial(1:6, ii, trial))
         this%w(ii) = weight(this%g(ii))
         call findSigmaAndAngle(this%g(ii), this%sigma(1, ii), this%sigma(2, ii), this%theta(ii))
         this%q(ii) = 0.0d0
         this%t(ii) = 0.0d0
      end do

      if (LIB_MG2D_DBG) then
         img_dbg = imin
         do ii = 1, this%n
            call add(img_dbg, this%g(ii))
         end do
         call writePng("multiple_output.png", img_dbg)
         print *, ""
      end if

      return
   end subroutine fit0

   subroutine findT(this, img, imin, imax, sigma)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      compute t-values on all gaussians
      type(MultipleGaussians2d), intent(inout)         ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img
      real(kind=real64), intent(in)                    ::      imin, imax
      real(kind=real64), intent(in)                    ::      sigma

      integer             ::      Nx, Ny
      integer             ::      ii

      real(kind=real64), dimension(:, :), allocatable    ::      img_clean, img_tmp

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      allocate (img_clean(0:Nx - 1, 0:Ny - 1))
      allocate (img_tmp(0:Nx - 1, 0:Ny - 1))

      call cleanImg(img, imin, imax, img_clean)
      call subtract(img_clean, this)
      this%t = 0.0d0
      do ii = 1, this%n
         img_tmp = img_clean
         call add(img_tmp, this%g(ii))
         this%t(ii) = getT(this%g(ii), img_tmp, sigma)

      end do

      return
   end subroutine findT

   subroutine findQ(this, img, imin, imax, sigma)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      compute q-values on all gaussians
      type(MultipleGaussians2d), intent(inout)         ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img
      real(kind=real64), intent(in)                    ::      imin, imax
      real(kind=real64), intent(in)                    ::      sigma

      integer             ::      Nx, Ny
      integer             ::      ii
      integer             ::      npx
      real(kind=real64)   ::      rss
      real(kind=real64), dimension(:, :), allocatable    ::      img_clean, img_tmp

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      allocate (img_clean(0:Nx - 1, 0:Ny - 1))
      allocate (img_tmp(0:Nx - 1, 0:Ny - 1))

      call cleanImg(img, imin, imax, img_clean)
      call subtract(img_clean, this)
      this%q = 0.0d0
      do ii = 1, this%n
         img_tmp = img_clean
         call add(img_tmp, this%g(ii))
         !---    old line 190523
         !       this%q(ii) = log( getRss( this%g(ii),img_tmp ) / (sigma*area( this%g(ii) )) )
         !---    new line 190523
         call computeRss(this%g(ii), img_tmp, rss, npx)
         this%q(ii) = 2*(6 + 1) + 2*npx*log(rss/max(sigma, sigma*npx))      !   AIC. Note degrees of freedom k = 6 (x,y,dmax,dmin,angle,intensity) + 1 for noise. log likelihood is L = - n log( <rss> )
         !        Note: scale rss by std dev of background to compare between different images.

      end do

      return
   end subroutine findQ

   subroutine groupSpots(this, Nx, Ny, group, nGroups)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*
      type(MultipleGaussians2d), intent(inout)         ::      this
      integer, intent(in)                              ::      Nx, Ny
      integer, dimension(0:, :), intent(inout)           ::      group
      integer, intent(out)                             ::      nGroups

      integer             ::      ii, jj, ng
      real(kind=real64)   ::      dd
      logical             ::      ok
      logical, dimension(this%n, this%n)        ::      overlapMatrix
      logical, dimension(this%n)               ::      alreadyGrouped

      !       first decide which are highly correlated
      call orderByWeight(this)
      overlapMatrix = .false.
      do ii = 1, this%n - 1
         do jj = ii + 1, this%n
            dd = overlap(this%g(ii), this%g(jj), Nx, Ny)
            overlapMatrix(jj, ii) = (dd > 0.25d0)              !   this parameter is not random. This is the overlap function for two identical circular peaks are separated so that their tops are flat.
            overlapMatrix(ii, jj) = overlapMatrix(jj, ii)
         end do
      end do

      !       now can choose groups of gaussians which are safe to fit together.
      do ii = 1, this%n
         alreadyGrouped(ii) = ignore(this, ii)                !   neat trick: say the ignored spots are grouped.
      end do
      group = 0
      nGroups = 0                  !   number of groups decided
      do ii = 1, this%n
         if (alreadyGrouped(ii)) cycle

         !---    haven't decided where ii should go yet, so ii must be in its own new group
         alreadyGrouped(ii) = .true.
         nGroups = nGroups + 1
         ng = 1                                  !   number in this (new) group is 1.
         group(0, nGroups) = ng
         group(ng, nGroups) = ii

         do jj = 1, this%n
            if (alreadyGrouped(jj)) cycle       !   can't add jj to ii's group

            !   now I know that jj is a potential new member of ii's group.
            !   check that it doesn't overlap any existing members of ii's group. Note that ii's group number is nGroups
            ok = .not. any(overlapMatrix(group(1:ng, nGroups), jj))

            if (ok) then
               !   I can add jj to ii's group as jj does not overlap with any member of ii's group.
               ng = ng + 1
               group(0, nGroups) = ng
               group(ng, nGroups) = jj
               alreadyGrouped(jj) = .true.
            end if
         end do
      end do

      if (LIB_MG2D_DBG) then
         print *, "Lib_MultipleGaussians2d::fit0 info - ", nGroups, " groups"
         do ii = 1, nGroups
            ng = group(0, ii)
            print *, "group ", ii, " members (", group(1:ng, ii), ")"
         end do
         print *, ""
      end if

      return
   end subroutine groupSpots

   subroutine fit1(this, img, group)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      fit subgroups of the gaussians to the image.
      !*      on input img should be cleaned up and set to LIB_G2D_IGNORE where necessary
      !*      This algorithm first subtracts the intensity due to gaussians outside the group
      !*      then fits the gaussians inside the group.

      type(MultipleGaussians2d), intent(inout)         ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img
      integer, dimension(0:, :), intent(in)              ::      group

      integer             ::      ii, jj, kk, mm
      integer             ::      Nx, Ny
      integer             ::      nGroups, nInGroup
      real(kind=real64)   ::      rss, oldrss

      integer, parameter           ::      NLOOPS = 100
      real(kind=real64), parameter ::      TOL = 1.0d-4
      type(MultipleGaussians2d)      ::      mg2d
      real(kind=real64), dimension(:, :), allocatable   ::      img_tmp

      nGroups = size(group, dim=2)
      mg2d = MultipleGaussians2d_ctor(this%n)
      allocate (img_tmp(0:size(img, dim=1) - 1, 0:size(img, dim=2) - 1))
      rss = getRss(this, img)
      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      do jj = 1, NLOOPS

         if (jj == 1) oldrss = 2*rss
         if (LIB_MG2D_DBG) print *, "Lib_MultipleGaussians2d::fit1 info - loop ", jj, " rss ", rss

         do kk = 1, nGroups                            !   k is group number
            nInGroup = group(0, kk)                    !   number of spots in group

            !---    generate a copy of the image with all peaks not in group deleted
            mm = 0                                  !   number of spots to fit in the group.
            img_tmp = img
            do ii = 1, this%n                        !   i is gaussian number
               if (ignore(this, ii)) cycle
               if (any(group(1:nInGroup, kk) == ii)) then
                  mm = mm + 1
                  mg2d%g(mm) = Gaussian2d_ctor(getDat(this%g(ii)))
               else
                  call subtract(img_tmp, this%g(ii))
               end if
            end do
            mg2d%n = mm
            img_tmp = max(0.0d0, img_tmp)

            !---    fit remainder
            call fit2(mg2d, img_tmp, rss)     !   note: rss here is the fit of the group to the image with the reset subtracted. It shouldn't be compared to the total fit.

            !---    extract
            mm = 0
            do ii = 1, this%n                    !   i is gaussian number
               if (ignore(this, ii)) cycle
               if (any(group(1:nInGroup, kk) == ii)) then
                  mm = mm + 1
                  this%g(ii) = Gaussian2d_ctor(getDat(mg2d%g(mm)))
               end if
            end do

         end do

         !---    are we done?
         if (LIB_MG2D_DBG) call report(this)
         rss = getRss(this, img)

         if (abs(rss - oldrss) < rss*TOL) exit
         !if (nGroups==1) exit                        !   can't do any better if there is only one group.
         if (rss > oldrss) exit
         oldrss = rss
         if (LIB_MG2D_DBG) print *, ""

      end do
      call delete(mg2d)

      do ii = 1, this%n
         if (ignore(this%g(ii))) then
            this%w(ii) = -1
         else
            this%w(ii) = weight(this%g(ii))
         end if
      end do

      return

   end subroutine fit1

   subroutine fit2(this, img, rss)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      fit a subgroup of the gaussians to the image.
      !*      on input img should be cleaned up and set to LIB_G2D_IGNORE where necessary
      !*      This algorithm first subtracts the intensity due to gaussians outside the group
      !*      then fits the gaussians inside the group.

      type(MultipleGaussians2d), intent(inout)         ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img
      real(kind=real64), intent(out)                   ::      rss

      integer                                         ::      Nx, Ny
      integer             ::      ii, jj

      real(kind=real64), dimension(:, :), allocatable    ::      dat, drss
      real(kind=real64)                               ::      oldrss, ff
      real(kind=real64)                               ::      x1, x2, x3, y1, y2, y3
      type(GoldenSection)                             ::      gold
      logical                                         ::      isConverged, isWithinTimeLimit
      integer, parameter           ::      NLOOPS = 100
      real(kind=real64), parameter ::      TOL = 1.0d-4
      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      !---    allocate memory and remove peaks not part of initial fitting
      allocate (dat(6, this%n))
      allocate (drss(6, this%n))

      !---    find derivatives wrt each gaussian
      do jj = 1, NLOOPS

         dat = getDat(this)
         call findDrss(this, img, rss, drss)

         if (jj == 1) oldrss = 2*rss
         if (mod(jj, 2) == 0) drss(4:6, :) = 0.0d0

         ff = 0
         do ii = 1, this%n
                    ff = ff + drss(1,ii)*drss(1,ii) + drss(2,ii)*drss(2,ii) + drss(3,ii)*drss(3,ii) + drss(4,ii)*drss(4,ii) + drss(5,ii)*drss(5,ii) + drss(6,ii)*drss(6,ii)  
         end do
         if (abs(ff) < 1.0d-8) exit

         x1 = 0.001*rss/ff
         call setDat(this, dat + x1*drss); y1 = getRss(this, img)

         x3 = -0.001*rss/ff
         call setDat(this, dat + x3*drss); y3 = getRss(this, img)

         if (y3 < y1) then
            x2 = x3; x1 = 0; y1 = rss
            do ii = 1, 20
               call setDat(this, dat + x2*drss); y2 = getRss(this, img)
               if (y2 > y1) then
                  y1 = y2; x1 = x2; exit
               else
                  y1 = y2
               end if
               x2 = x2*2

            end do
            if (x1 == 0) stop
            x3 = 0; y3 = rss
         else
            x2 = x1; x3 = 0; y3 = rss; 
            do ii = 1, 20
               call setDat(this, dat + x2*drss); y2 = getRss(this, img)
               if (y2 > y3) then
                  y3 = y2; x3 = x2; exit
               else
                  y3 = y2
               end if
               x2 = x2*2
            end do
            if (x3 == 0) stop
            x1 = 0; y1 = rss
         end if
         gold = GoldenSection_ctor(x1, x3, y1, y3)

         do
            call nextPoint(gold, x2)
            call setDat(this, dat + x2*drss)
            y2 = getRss(this, img)
            call minimise(gold, x2, y2, isConverged, isWithinTimeLimit)
            if (isConverged) exit
            if (.not. isWithinTimeLimit) exit
         end do

         call setDat(this, dat + x2*drss)
         rss = getRss(this, img)

         if (abs(rss - oldrss) < rss*TOL) exit
         oldrss = rss

      end do

      return
   end subroutine fit2

   pure subroutine cleanImg(img_in, imin, imax, img_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      tidy the input image by removing background, preserving dead pixels
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img_in
      real(kind=real64), intent(in)                    ::      imin, imax
      real(kind=real64), dimension(0:, 0:), intent(out)  ::      img_out
      if (imax > imin) then
         where (img_in == LIB_G2D_IGNORE)
            img_out = LIB_G2D_IGNORE
         elsewhere(img_in < imax)
            img_out = img_in - imin
         elsewhere
            img_out = LIB_G2D_IGNORE
         end where
      else
         where (img_in == LIB_G2D_IGNORE)
            img_out = LIB_G2D_IGNORE
         elsewhere(img_in > imax)
            img_out = img_in - imin
         elsewhere
            img_out = LIB_G2D_IGNORE
         end where
      end if
      return
   end subroutine cleanImg

   subroutine orderByWeight(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      reorder the list of gaussians by weight
      type(MultipleGaussians2d), intent(inout)     ::      this
      integer                 ::      ii, jj
      real(kind=real64), dimension(0:6, this%n)     ::      dat
      integer, dimension(this%n)                 ::       indx

      do ii = 1, this%n
         dat(0, ii) = this%w(ii)
         dat(1:6, ii) = getDat(this, ii)
         indx(ii) = ii
      end do
      call quicksort(this%w(1:this%n), indx, back=.true.)
      do ii = 1, this%n
         jj = indx(ii)
         this%w(ii) = dat(0, jj)
         this%g(ii) = Gaussian2d_ctor(dat(1:6, jj))
      end do

      return
   end subroutine orderByWeight

   pure real(kind=real64) function getRss0(this, img)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given the input density function, return the residual sum of squares
      type(MultipleGaussians2d), intent(in)                    ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)           ::      img

      integer             ::      Nx, Ny
      integer             ::      ix, iy, ii
      real(kind=real64)   ::      gg, ss, xx, yy

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      getRss0 = 0
      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            if (img(ix, iy) == LIB_G2D_IGNORE) cycle
            ss = -img(ix, iy)
            do ii = 1, this%n
               if (ignore(this, ii)) cycle
               xx = ix - this%g(ii)%x0
               yy = iy - this%g(ii)%y0
               gg = this%g(ii)%Dxx*xx*xx + 2*this%g(ii)%Dxy*xx*yy + this%g(ii)%Dyy*yy*yy
               if (gg > 4.5d0) cycle               !   only search to 3 sigma
               gg = this%g(ii)%f0*Exp(-gg)
               ss = ss + gg
            end do
            getRss0 = getRss0 + ss*ss
         end do
      end do

      return
   end function getRss0

   !---

   subroutine findDrss0(this, img, rss, drss)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return the first derivatives of the residual sum of squares with respect to change of parameters
      type(MultipleGaussians2d), intent(in)            ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img
      real(kind=real64), intent(out)                   ::      rss
      real(kind=real64), dimension(:, :), intent(out)    ::      drss

      integer             ::      Nx, Ny
      integer             ::      ix, iy, ii
      real(kind=real64)   ::      gg, ss, xx, yy
      real(kind=real64), dimension(this%n) ::      ggi

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      rss = 0
      drss = 0

      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            if (img(ix, iy) == LIB_G2D_IGNORE) cycle
            ss = -img(ix, iy)
            ggi = 0
            do ii = 1, this%n
               if (ignore(this, ii)) cycle
               xx = (ix - this%g(ii)%x0)
               yy = (iy - this%g(ii)%y0)
               gg = this%g(ii)%Dxx*xx*xx + 2*this%g(ii)%Dxy*xx*yy + this%g(ii)%Dyy*yy*yy
               if (gg > 4.5d0) cycle               !   only search to 3 sigma
               ggi(ii) = Exp(-gg)
               ss = ss + this%g(ii)%f0*ggi(ii)
            end do

            rss = rss + ss*ss

            do ii = 1, this%n
               drss(1, ii) = drss(1, ii) + 2*ss*(2*this%g(ii)%f0*ggi(ii)*(this%g(ii)%Dxx*xx + this%g(ii)%Dxy*yy))        !    = d f/d x
               drss(2, ii) = drss(2, ii) + 2*ss*(2*this%g(ii)%f0*ggi(ii)*(this%g(ii)%Dxy*xx + this%g(ii)%Dyy*yy))
               drss(3, ii) = drss(3, ii) + 2*ss*(ggi(ii))        !    = d f/d f0
               drss(4, ii) = drss(4, ii) + 2*ss*(-this%g(ii)%f0*ggi(ii)*xx*xx)        !    = d f/d Dxx
               drss(5, ii) = drss(5, ii) + 2*ss*(-2*this%g(ii)%f0*ggi(ii)*xx*yy)
               drss(6, ii) = drss(6, ii) + 2*ss*(-this%g(ii)%f0*ggi(ii)*yy*yy)
            end do

         end do
      end do

      return
   end subroutine findDrss0

   pure function get0(this) result(dat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return data to construct gaussian from 6 parameters in order ( x0,y0,f0,Dxx,Dxy,Dyy )
      type(MultipleGaussians2d), intent(in)         ::      this
      real(kind=real64), dimension(6, this%n)        ::      dat
      integer             ::      ii
      do ii = 1, this%n
         dat(:, ii) = getDat(this%g(ii))
      end do
      return
   end function get0

   pure function get1(this, i) result(dat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return data to construct gaussian from 6 parameters in order ( x0,y0,f0,Dxx,Dxy,Dyy )
      type(MultipleGaussians2d), intent(in)        ::      this
      integer, intent(in)                          ::      i
      real(kind=real64), dimension(6)              ::      dat
      dat = getDat(this%g(i))
      return
   end function get1

   pure integer function getN0(this, ignore)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(MultipleGaussians2d), intent(in)         ::     this
      logical, intent(in), optional                 ::      ignore
      getN0 = this%n
      if (present(ignore)) then
         if (ignore) getN0 = getN0 - count(this%w(1:this%n) < 0)
      end if

      return
   end function getN0

   pure subroutine setN0(this, n)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(MultipleGaussians2d), intent(inout)      ::     this
      integer, intent(in)                           ::     n
      this%n = n
      return
   end subroutine setN0

   pure function getG0(this, i) result(g)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return data to construct gaussian from 6 parameters in order ( x0,y0,f0,Dxx,Dxy,Dyy )
      type(MultipleGaussians2d), intent(in)        ::      this
      integer, intent(in)                          ::      i
      type(Gaussian2d)        ::      g
      g = Gaussian2d_ctor(getDat(this%g(i)))
      return
   end function getG0

   pure subroutine set0(this, dat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      construct gaussian from 6 parameters in order ( x0,y0,f0,Dxx,Dxy,Dyy )
      !*      note: does not set weights - use full constructor for that
      type(MultipleGaussians2d), intent(inout)                 ::      this
      real(kind=real64), dimension(:, :), intent(in)        ::      dat
      integer             ::      ii, nn

      nn = size(dat, dim=2)
      do ii = 1, nn
         this%g(ii) = Gaussian2d_ctor(dat(:, ii))
      end do
      this%w = 0
      this%t = 0
      this%q = 0
      this%sigma = 0
      this%theta = 0
      this%n = nn
      return
   end subroutine set0

   pure subroutine set1(this, i, dat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      construct gaussian from 6 parameters in order ( x0,y0,f0,Dxx,Dxy,Dyy )
      !*      note: does not set weights - use full constructor for that
      type(MultipleGaussians2d), intent(inout)                 ::      this
      integer, intent(in)                                      ::      i
      real(kind=real64), dimension(6), intent(in)               ::      dat
      this%g(i) = Gaussian2d_ctor(dat(:))
      this%w(i) = 0
      this%t(i) = 0
      this%q(i) = 0
      this%sigma(:, i) = 0
      this%theta(i) = 0
      return
   end subroutine set1

   !---    trivial inquiry functions

   pure real(kind=real64) function getX0(this, i)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(MultipleGaussians2d), intent(in)        ::      this
      integer, intent(in)                          ::      i
      getX0 = this%g(i)%x0
      return
   end function getX0

   pure real(kind=real64) function getY0(this, i)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(MultipleGaussians2d), intent(in)        ::      this
      integer, intent(in)                          ::      i
      getY0 = this%g(i)%y0
      return
   end function getY0

   pure real(kind=real64) function getF0(this, i)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(MultipleGaussians2d), intent(in)        ::      this
      integer, intent(in)                          ::      i
      getF0 = this%g(i)%f0
      return
   end function getF0

   pure real(kind=real64) function getTheta0(this, i)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(MultipleGaussians2d), intent(in)        ::      this
      integer, intent(in)                          ::      i
      getTheta0 = this%theta(i)
      return
   end function getTheta0

   pure function getD0(this, i) result(D)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(MultipleGaussians2d), intent(in)        ::      this
      real(kind=real64), dimension(2, 2)            ::      D
      integer, intent(in)                          ::      i
      D(1, 1) = this%g(i)%Dxx
      D(2, 1) = this%g(i)%Dxy
      D(1, 2) = this%g(i)%Dxy
      D(2, 2) = this%g(i)%Dyy
      return
   end function getD0

   pure real(kind=real64) function getSigma0(this, i, larger)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(MultipleGaussians2d), intent(in)        ::      this
      integer, intent(in)                          ::      i
      logical, intent(in), optional                 ::      larger
      if (present(larger)) then
         if (larger) then
            getSigma0 = this%sigma(1, i)
         else
            getSigma0 = this%sigma(2, i)
         end if
      else
         getSigma0 = sqrt(this%sigma(1, i)*this%sigma(2, i))
      end if
      return
   end function getSigma0

   pure real(kind=real64) function getT0(this, i)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(MultipleGaussians2d), intent(in)        ::      this
      integer, intent(in)                          ::      i
      getT0 = this%t(i)
      return
   end function getT0

   pure real(kind=real64) function getQ0(this, i)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(MultipleGaussians2d), intent(in)        ::      this
      integer, intent(in)                          ::      i
      getQ0 = this%q(i)
      return
   end function getQ0

   pure real(kind=real64) function getW0(this, i)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(MultipleGaussians2d), intent(in)        ::      this
      integer, intent(in)                          ::      i
      getW0 = this%w(i)
      return
   end function getW0

   pure subroutine translate0(this, i, x)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      translate centre of gaussian
      type(MultipleGaussians2d), intent(inout)                  ::      this
      integer, intent(in)                          ::      i
      real(kind=real64), dimension(2), intent(in)       ::      x
      this%g(i)%x0 = this%g(i)%x0 + x(1)
      this%g(i)%y0 = this%g(i)%y0 + x(2)
      return
   end subroutine translate0

   pure subroutine translate1(this, i, x, y)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      translate centre of gaussian
      type(MultipleGaussians2d), intent(inout)                  ::      this
      integer, intent(in)                          ::      i
      real(kind=real64), intent(in)                    ::      x, y
      this%g(i)%x0 = this%g(i)%x0 + x
      this%g(i)%y0 = this%g(i)%y0 + y
      return
   end subroutine translate1

   pure subroutine translate2(this, i, x, y)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      translate centre of gaussian
      type(MultipleGaussians2d), intent(inout)                  ::      this
      integer, intent(in)                          ::      i
      integer, intent(in)                              ::      x, y
      this%g(i)%x0 = this%g(i)%x0 + x
      this%g(i)%y0 = this%g(i)%y0 + y
      return
   end subroutine translate2

   pure subroutine translate3(this, x)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      translate centre of gaussian
      type(MultipleGaussians2d), intent(inout)         ::      this
      real(kind=real64), dimension(2), intent(in)       ::      x
      integer         ::      ii
      do ii = 1, this%n
         this%g(ii)%x0 = this%g(ii)%x0 + x(1)
         this%g(ii)%y0 = this%g(ii)%y0 + x(2)
      end do
      return
   end subroutine translate3

   pure subroutine translate4(this, x, y)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      translate centre of gaussian
      type(MultipleGaussians2d), intent(inout)         ::      this
      real(kind=real64), intent(in)                    ::      x, y
      integer         ::      ii
      do ii = 1, this%n
         this%g(ii)%x0 = this%g(ii)%x0 + x
         this%g(ii)%y0 = this%g(ii)%y0 + y
      end do
      return
   end subroutine translate4

   pure subroutine translate5(this, x, y)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      translate centre of gaussian
      type(MultipleGaussians2d), intent(inout)         ::      this
      integer, intent(in)                              ::      x, y
      integer         ::      ii
      do ii = 1, this%n
         this%g(ii)%x0 = this%g(ii)%x0 + x
         this%g(ii)%y0 = this%g(ii)%y0 + y
      end do
      return
   end subroutine translate5

   pure subroutine xscale(this, x)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      scale length components of gaussians
      type(MultipleGaussians2d), intent(inout)         ::      this
      real(kind=real64), intent(in)       ::      x
      integer         ::      ii
      do ii = 1, this%n
         this%g(ii)%x0 = this%g(ii)%x0*x
         this%g(ii)%y0 = this%g(ii)%y0*x
         this%g(ii)%Dxx = this%g(ii)%Dxx/(x*x)
         this%g(ii)%Dxy = this%g(ii)%Dxy/(x*x)
         this%g(ii)%Dyy = this%g(ii)%Dyy/(x*x)
         this%w(ii) = this%w(ii)*(x*x)
         this%sigma(1:2, ii) = this%sigma(1:2, ii)*x
      end do
      return
   end subroutine xscale

   subroutine add0(this, that)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      add the contents of that into this
      type(MultipleGaussians2d), intent(inout)     ::      this
      type(MultipleGaussians2d), intent(in)        ::      that

      integer         ::      ii, nn

      if (this%n + that%n > size(this%w)) then
                stop "Lib_MultipleGaussians2d::add0 error - attempt to add more gaussians than allocated space. Easy code fix, but must be an error??"
      end if

      nn = this%n
      do ii = 1, that%n
         if (ignore(that, ii)) cycle
         nn = nn + 1
         this%g(nn) = Gaussian2d_ctor(getDat(that%g(ii)))
         this%w(nn) = that%w(ii)
         this%t(nn) = that%t(ii)
         this%q(nn) = that%q(ii)
         this%sigma(:, nn) = that%sigma(:, ii)
         this%theta(nn) = that%theta(ii)
      end do

      this%n = nn

      return
   end subroutine add0

   subroutine subtract0(img, this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      remve this gaussian from the image
      type(MultipleGaussians2d), intent(in)                         ::  this
      real(kind=real64), dimension(0:, 0:), intent(inout)    ::  img

      integer             ::      Nx, Ny
      integer             ::      ix, iy, ii
      real(kind=real64)   ::      ff, gg, xx, yy

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            ff = img(ix, iy)
            if (ff == LIB_G2D_IGNORE) cycle
            do ii = 1, this%n
               if (ignore(this, ii)) cycle
               xx = (ix - this%g(ii)%x0)
               yy = (iy - this%g(ii)%y0)
               gg = this%g(ii)%Dxx*xx*xx + 2*this%g(ii)%Dxy*xx*yy + this%g(ii)%Dyy*yy*yy
               if (gg > 4.5d0) cycle           !   only consider up to 3 sigma.
               gg = Exp(-gg)
               ff = ff - this%g(ii)%f0*gg
            end do
            img(ix, iy) = ff - gg
         end do
      end do

      return
   end subroutine subtract0

   subroutine add1(img, this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      add this gaussian to the image
      type(MultipleGaussians2d), intent(in)                         ::  this
      real(kind=real64), dimension(0:, 0:), intent(inout)    ::  img

      integer             ::      Nx, Ny
      integer             ::      ix, iy, ii
      real(kind=real64)   ::      ff, gg, xx, yy

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            ff = img(ix, iy)
            if (ff == LIB_G2D_IGNORE) cycle
            do ii = 1, this%n
               if (ignore(this, ii)) cycle
               xx = (ix - this%g(ii)%x0)
               yy = (iy - this%g(ii)%y0)
               gg = this%g(ii)%Dxx*xx*xx + 2*this%g(ii)%Dxy*xx*yy + this%g(ii)%Dyy*yy*yy
               if (gg > 4.5d0) cycle           !   only consider up to 3 sigma.
               gg = Exp(-gg)
               ff = ff + this%g(ii)%f0*gg
            end do
            img(ix, iy) = ff
         end do
      end do

      return
   end subroutine add1

   logical function ignore0(this, q, t)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(Gaussian2d), intent(in)             ::      this
      real(kind=real64), intent(in), optional   ::      q, t
      real(kind=real64)       ::      ee

      ee = eccentricity(this)
      ignore0 = (ee > E_THRESH)
      ignore0 = ignore0 .or. (this%f0 < DETECTFMIN)
      ignore0 = ignore0 .or. (this%Dxx*this%Dyy < this%Dxy*this%Dxy)
      if (present(q)) ignore0 = ignore0 .or. (q > Q_THRESH)
      if (present(t)) ignore0 = ignore0 .or. (t < T_THRESH)

      if (LIB_MG2D_DBG) then
         print *, " eccentricity ", ee, " > e_thresh   (", E_THRESH, ")?  ", (ee > E_THRESH)
         print *, " intensity    ", this%f0, " < DETECTFMIN (", DETECTFMIN, ")?  ", (this%f0 < DETECTFMIN)
         if (present(q)) print *, " q            ", q, " > Q_THRESH ? (", Q_THRESH, ")?  ", (q > Q_THRESH)
         if (present(t)) print *, " t            ", q, " < T_THRESH ? (", T_THRESH, ")?  ", (t < T_THRESH)
        print *, " det[D]       ", this%Dxx*this%Dyy - this%Dxy*this%Dxy, " < 0?          ", (this%Dxx*this%Dyy < this%Dxy*this%Dxy)
         print *, " ignore?      ", ignore0
      end if
      return
   end function ignore0

   elemental logical function ignore1(this, i)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(MultipleGaussians2d), intent(in)        ::      this
      integer, intent(in)                          ::      i
      ignore1 = (this%w(i) < 0)
      return
   end function ignore1

   subroutine setIgnore(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(MultipleGaussians2d), intent(inout)     ::      this
      integer     ::      ii
      do ii = 1, this%n
         if (ignore0(this%g(ii), this%q(ii), this%t(ii))) this%w(ii) = -1.0d0
      end do
      return
   end subroutine setIgnore

   function loglikelihood0(this, img) result(lambda)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given the input density function, return the log likelihood ( ignoring the noise term )
      !*      lambda = - n log (rss/n)
      type(MultipleGaussians2d), intent(in)            ::      this
      real(kind=real64), dimension(0:, 0:), intent(in)   ::      img
      real(kind=real64)                 ::      lambda

      integer             ::      Nx, Ny, n
      integer             ::      ix, iy, ii
      real(kind=real64)   ::      gg, ss, xx, yy, rss

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      n = 0
      rss = 0
      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            if (img(ix, iy) == LIB_G2D_IGNORE) cycle
            ss = -img(ix, iy)
            do ii = 1, this%n
               if (ignore(this, ii)) cycle
               xx = (ix - this%g(ii)%x0)
               yy = (iy - this%g(ii)%y0)
               gg = this%g(ii)%Dxx*xx*xx + 2*this%g(ii)%Dxy*xx*yy + this%g(ii)%Dyy*yy*yy
               gg = Exp(-gg)
               ss = ss + this%g(ii)%f0*gg
            end do
            rss = rss + ss*ss
            n = n + 1
         end do
      end do

      if (n > 0) then
         lambda = -n*log(rss/n)
      else
         lambda = 0
      end if

      return
   end function loglikelihood0

   subroutine setLib_MultipleGaussians2d_dbg(dbg)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      logical, intent(in)      ::      dbg
      LIB_MG2D_DBG = dbg
      return
   end subroutine setLib_MultipleGaussians2d_dbg

   subroutine ignore_individual_guassian(this, ii)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !sets the weigth of an individual gaussian to -1 so various functions ignore it.
      type(MultipleGaussians2d), intent(inout)             ::      this
      integer, intent(in)                         ::       ii
      this%w(ii) = -1.0d0
      return
   end subroutine ignore_individual_guassian

end module Lib_MultipleGaussians2d
