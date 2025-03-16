
module Lib_Maxima2d
!---^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_Maxima2d from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      This simple module contains functions for
!*      finding the number of maxima in a 2d function
   use iso_fortran_env
   use Lib_RidlerCalvard
   implicit none
   private

   !---

   logical, private          ::      LIBMAX_DBG = .false.

   !---

   public          ::      findMaxima
   public          ::      mergeMaxima
   public          ::      paddedROI
   public          ::      dimensionsOfROI
   public          ::      nonPercolatingROI
   public          ::      extractROI
   public          ::      bxFilter
   public          ::      polyFilter
   public          ::      polyFilter2
   public          ::      rawMaxima2
   public          ::      findSigma
   public          ::      intensityShift
   public          ::      setLib_Maxima2d_dbg
!        public          ::      estimateMaxpixFromDarkPixels

!        public          ::      screendump

   !---

contains
!---^^^^^^^^

   subroutine setLib_Maxima2d_dbg(dbg)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      logical, intent(in)      ::      dbg
      LIBMAX_DBG = dbg
      return
   end subroutine setLib_Maxima2d_dbg

   subroutine extractROI(img_in, i, indx, left, right, top, bottom, xMax, yMax, img_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(:, :), intent(in)     ::      img_in             !   (1:nx,1:nz)
      integer, intent(in)                              ::      i
      integer, dimension(:, :), intent(in)               ::      indx
      integer, dimension(:), intent(in)                 ::      left, right, top, bottom               !   (1:nMax)
      integer, intent(out)                             ::      xMax, yMax
      real(kind=real64), dimension(:, :), intent(out)    ::      img_out             !   (xMax,yMax)

      integer             ::      ix, iy, ii, jx, jy

      img_out = 0.0d0
      do iy = bottom(i), top(i)
         jy = iy + 1 - bottom(i)
         do ix = left(i), right(i)
            jx = ix + 1 - left(i)
            ii = indx(ix, iy)
            if ((ii - i) == 0) img_out(jx, jy) = img_in(ix, iy)            !   only if pixel part of roi i
         end do
      end do
      xMax = right(i) - left(i) + 1
      yMax = top(i) - bottom(i) + 1

      return
   end subroutine extractROI

   !----

   subroutine nonPercolatingROI(img, b, sig, pMax, xMax, rMax, indx, nMax, ok)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      compute regions of interest based on a background level b and std dev sig
      !*      such that the largest regions are smaller than xMax.
      !*      and the threshold level not more than

      !*      combine regions of interest using NSEW mask
      !*
      !*      index each point as belonging to a single roi indx(i,j) = 1:nMax

      real(kind=real64), dimension(:, :), intent(in)     ::      img             !   (1:nx,1:nz) -
      real(kind=real64), intent(in)                    ::      b, sig
      real(kind=real64), intent(in)                    ::      pMax            !   maximum permitted multiple of std dev for "successful" roi decomposition
      real(kind=real64), intent(in)                    ::      rMax            !   maximum permitted padding region
      integer, intent(in)                              ::      xMax            !   maximum permitted roi dimension
      integer, dimension(:, :), intent(out)              ::      indx            !   (1:nx,1:nz) -
      integer, intent(out)                             ::      nMax
      logical, intent(out)                             ::      ok

      integer, parameter               ::      NR = 19 !number of padding radii to try, with RMIN=1.0do, gives max radius of 20.0d0
      real(kind=real64), parameter     ::      RMIN = 1.0d0
      real(kind=real64)               ::      rr, rlast, pp, plast, dr
      integer                         ::      ii

      !---    check for quick escape for tiny image
      if (size(img, dim=1)*size(img, dim=2) <= xMax*xMax) then
         call nonPercolatingROIr(img, b, sig, rMax, xMax, indx, nMax, pp)
         ok = .true.
         return
      end if

            write(*,fmt='(3(a,f12.5))') "Lib_Maxima2d::nonPercolatingROI info - hunting for rois with thresh/sd < ",pMax," plus padding <",rMax," giving roi size < ",real(xMax)
      dr = (rMax - RMIN)/NR
      ok = .true.
      rlast = 0.0d0
      plast = 0.0d0
      do ii = 0, NR !
         rr = RMIN + ii*DR
         call nonPercolatingROIr(img, b, sig, rr, xMax, indx, nMax, pp)

         if (ii == 0) then
            if (pp > pMax) then
               ok = .false.
               return
            end if
         else
            if (pp > pMax) then
               rr = ((pMax - plast)/(pp - plast))*(rr - rlast) + rlast
               call nonPercolatingROIr(img, b, sig, rr, xMax, indx, nMax, pp)
               exit
            end if
         end if
         rlast = rr
         plast = pp
      end do
      write (*, fmt='(2(a,f12.5))') "Lib_maxima2d::nonPercolatingROI info - thresh/sd ", pp, " padding ", rr
      return
   end subroutine nonPercolatingROI

   subroutine nonPercolatingROIr(img, b, sig, r, xMax, indx, nMax, p)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      compute regions of interest based on a background level b and std dev sig
      !*      such that the largest regions are about xMax.
      !*      pad each roi with radius r pixels
      !*      combine regions of interest using NSEW mask
      !*
      !*      index each point as belonging to a single roi indx(i,j) = 1:nMax

      real(kind=real64), dimension(:, :), intent(in)     ::      img             !   (1:nx,1:nz) -
      real(kind=real64), intent(in)                    ::      b, sig
      real(kind=real64), intent(in)                    ::      r
      integer, intent(in)                              ::      xMax
      integer, dimension(:, :), intent(out)              ::      indx            !   (1:nx,1:nz) -
      integer, intent(out)                             ::      nMax
      real(kind=real64), intent(out)                   ::      p

      integer, parameter                           ::      NLOOPS = 16
      integer                                     ::      ii
      integer, dimension(:), pointer                ::      left, right, top, bottom
      integer, dimension(NLOOPS)                   ::      xmx, ymx
      real(kind=real64), dimension(NLOOPS)         ::      pm
      real(kind=real64)                           ::      f_c

      real(kind=real64)                           ::      p_last, dp, xx

      p_last = 0.0d0                                      !   multiplier on f_c = b + p*sig. Assume last choice was zero

      dp = 1.0d0                                          !   go up in units of sig, trying to find a non-percolating threshold
      p = huge(1.0); xmx = 0; ymx = 0; pm = 0.0d0
      do ii = 1, NLOOPS                                    !   don't go on forever - if you really can't find a solution then you may have to accept defeat.

         pm(ii) = p_last + dp
         f_c = b + pm(ii)*sig

         !---    find the regions of interest based on this guess for the threshold level
         call paddedROI(img, f_c, r, indx, nMax)
         if (nMax > 1) then
            allocate (left(nMax)); allocate (right(nMax)); allocate (top(nMax)); allocate (bottom(nMax))
            call dimensionsOfROI(indx, nMax, left, right, top, bottom, xmx(ii), ymx(ii))

            deallocate (left); deallocate (right); deallocate (top); deallocate (bottom)
         else
            xmx(ii) = size(img, dim=1)
            ymx(ii) = size(img, dim=2)
         end if

         !---    now look at the largest ROI. Is it too big?
         xx = sqrt(real(xmx(ii), kind=real64)*ymx(ii))
         if (xx > xMax) then
            !   yes. Too big. Have to look at a higher threshold.
            p_last = pm(ii)
            cycle
         else
            !   no. Not too big, so somewhere between p_last and pp is the ideal threshold.
            dp = dp/2
            if (xx > 0.75d0*xMax) then
               !   now within 75% of the ideal ROI size, so stick.
               p = pm(ii)
               exit
            else if (ii == NLOOPS) then
               p = pm(ii)
            end if

         end if
      end do
      write (*, fmt='(3(a,f12.5))') "Lib_Maxima2d::nonPercolatingROIr info - thresh/sd ", p, " padding ", r, " roi size ", xx !," f_c = ",f_c,"(",int(f_c*256),")"

      return
   end subroutine nonPercolatingROIr

   subroutine paddedROI(img, f_c, r, indx, nMax)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      scan the 2d image
      !*      find the number of roi containing at least one pixel over the level f_c
      !*      pad each roi with radius r pixels
      !*      combine regions of interest using NSEW mask
      !*
      !*      index each point as belonging to a single roi indx(i,j) = 1:nMax
      !*      optionally remove

      real(kind=real64), dimension(:, :), intent(in)     ::      img             !   (1:nx,1:nz) -
      real(kind=real64), intent(in)                    ::      f_c
      real(kind=real64), intent(in)                    ::      r
      integer, dimension(:, :), intent(out)              ::      indx            !   (1:nx,1:nz) -
      integer, intent(out)                             ::      nMax

      integer                 ::      nx, ny, nr
      integer                 ::      iy, ix, jy, jx, ii, jj, ri, rj, kk, kx, ky
      real(kind=real64)       ::      rr, fbar

      integer, dimension(:), allocatable                ::      reindex, rere

      integer, dimension(:, :), allocatable              ::      pad_indx        !   copy of index, padded so dont need to check bounds

      logical, dimension(:, :), allocatable              ::      pad_mask        !   circular mask of radius r

      logical                 ::      done
      nx = size(img, dim=1)            !   note that x- is row index
      ny = size(img, dim=2)            !   and y- is column index. Otherwise it gets confusing...
      if (LIBMAX_DBG) print *, "Lib_Maxima2d::paddedROI info - pixels = ", nx*ny
      if (LIBMAX_DBG) print *, "Lib_Maxima2d::paddedROI info - f_c = ", f_c, int(f_c*256)

      !---    establish the padding region
      nr = ceiling(r)                 !   number of pixels required for padding
      allocate (pad_mask(-nr:nr, -nr:nr))
      if (LIBMAX_DBG) print *, "pad mask"
      do iy = -nr, nr
         do ix = -nr, nr
            rr = ix*ix + iy*iy
            pad_mask(ix, iy) = (rr <= r*r)
         end do
         if (LIBMAX_DBG) write (*, fmt='(100l2)') pad_mask(:, iy)
      end do

      !---    first step: each maximum pixel is identified and labelled. Note maxima are in range 2:nx-1, ignoring boundary
      !
      indx = 0
      call aboveThreshold(img(2:Nx - 1, 2:Ny - 1), f_c, indx(2:Nx - 1, 2:Ny - 1), nMax)        !   test contig over thresh method
      if (nMax == 0) then
         deallocate (pad_mask)
         return                                                     !   well, that was easy...
      end if
      if (LIBMAX_DBG) print *, "Lib_Maxima2d::paddedROI info - rawMaxima = ", nMax

      !---    second step: splat padding region.
      allocate (pad_indx(1:nx, 1:ny))
      pad_indx(1:nx, 1:ny) = indx(1:nx, 1:ny)
!

      do iy = 2, ny - 1                  !   Note maxima are in range 2:nx-1, ignoring boundary
         do ix = 2, nx - 1
            ii = indx(ix, iy)
            if (ii > 0) then        !   over threshold

               do jy = -nr, nr
                  do jx = -nr, nr
                     if (.not. pad_mask(jx, jy)) cycle

                     kx = max(1, min(nx, ix + jx))
                     ky = max(1, min(ny, iy + jy))

                     if (indx(kx, ky) == 0) pad_indx(kx, ky) = ii
                  end do
               end do

            end if
         end do
      end do
      indx(1:nx, 1:ny) = pad_indx(1:nx, 1:ny)

      deallocate (pad_indx)
      if (LIBMAX_DBG) then
         print *, "splat mask"
         call screendump(img, indx)
         print *, ""
      end if
      if (nMax == 1) then
         deallocate (pad_mask)
         return                                                     !   well, that was easy...
      end if

      !---    third step: connect regions where pixels touch with different indices
      !   construct the array reindex(1:nMax)
      !   where all pixels with reindex( indx(i,j) ) = k belong to roi k.

      allocate (reindex(0:nMax))
      do ii = 0, nMax
         reindex(ii) = ii            !   start with assumption all roi are disconnected
      end do

!

      do
         done = .true.
         do iy = 2, ny - 1
            do ix = 2, nx - 1

               ii = indx(ix, iy)
               if (ii == 0) cycle                      !   ignore - not in a roi
               ri = reindex(ii)

               !---    checking neighbours NSEW are easy - if both over threshold they are same roi
               jj = indx(ix, iy - 1)                    !   check point north
               rj = reindex(jj)
               if (rj*(rj - ri) /= 0) then               !   point is ascribed to a different roi
                  ri = min(ri, rj)
                  reindex(jj) = ri
                  done = .false.
               end if

               jj = indx(ix - 1, iy)                    !   check point west
               rj = reindex(jj)
               if (rj*(rj - ri) /= 0) then               !   point is ascribed to a different roi
                  ri = min(ri, rj)
                  reindex(jj) = ri
                  done = .false.
               end if

               jj = indx(ix + 1, iy)                    !   check point east
               rj = reindex(jj)
               if (rj*(rj - ri) /= 0) then               !   point is ascribed to a different roi
                  ri = min(ri, rj)
                  reindex(jj) = ri
                  done = .false.
               end if

               jj = indx(ix, iy + 1)                    !   check point south
               rj = reindex(jj)
               if (rj*(rj - ri) /= 0) then               !   point is ascribed to a different roi
                  ri = min(ri, rj)
                  reindex(jj) = ri
                  done = .false.
               end if

               !---    checking neighbours NE,NW,SE,SW is harder, as we need to be topologically correct.

               jj = indx(ix - 1, iy - 1)                    !   check point north west
               rj = reindex(jj)
               if (rj*(rj - ri) /= 0) then               !   point is ascribed to a different roi
                  fbar = sum(img(ix - 1:ix, iy - 1:iy))
                  if (fbar > 4*f_c) then              !   contiguous
                     ri = min(ri, rj)
                     reindex(jj) = ri
                     done = .false.
                  end if
               end if

               jj = indx(ix + 1, iy - 1)                   !   check point north east
               rj = reindex(jj)
               if (rj*(rj - ri) /= 0) then               !   point is ascribed to a different roi
                  fbar = sum(img(ix:ix + 1, iy - 1:iy))
                  if (fbar > 4*f_c) then              !   contiguous
                     ri = min(ri, rj)
                     reindex(jj) = ri
                     done = .false.
                  end if
               end if

               jj = indx(ix - 1, iy + 1)                   !   check point south west
               rj = reindex(jj)
               if (rj*(rj - ri) /= 0) then               !   point is ascribed to a different roi
                  fbar = sum(img(ix - 1:ix, iy:iy + 1))
                  if (fbar > 4*f_c) then              !   contiguous
                     ri = min(ri, rj)
                     reindex(jj) = ri
                     done = .false.
                  end if
               end if

               jj = indx(ix + 1, iy + 1)                   !   check point south east
               rj = reindex(jj)
               if (rj*(rj - ri) /= 0) then               !   point is ascribed to a different roi
                  fbar = sum(img(ix:ix + 1, iy:iy + 1))
                  if (fbar > 4*f_c) then              !   contiguous
                     ri = min(ri, rj)
                     reindex(jj) = ri
                     done = .false.
                  end if
               end if

               reindex(ii) = ri

            end do
         end do

         if (done) exit
      end do

      allocate (rere(1:nMax))
      rere = 0
      kk = 0      !   number of distinct roi
      jj = 0      !   highest number roi seen
      do ii = 1, nMax
         if (reindex(ii) > jj) then
            !   this roi is the highest label seen, so must be a new roi
            jj = reindex(ii)
            kk = kk + 1
            rere(jj) = kk           !   the roi indexed jj is in fact only the kkth roi seen
         end if
      end do

      do ii = 1, nMax
         ri = reindex(ii)
         reindex(ii) = rere(ri)
      end do

      deallocate (rere)
      nMax = kk

      !---    rewrite the index
      do iy = 1, ny
         do ix = 1, nx
            ii = indx(ix, iy)
            ri = reindex(ii)
            indx(ix, iy) = ri
         end do
      end do

      if (LIBMAX_DBG) then
         print *, "reindexed"
         call screendump(img, indx, txt=.true.)
         print *, ""
      end if
      if (LIBMAX_DBG) print *, "Lib_Maxima2d::paddedROI info - roi = ", nMax

      !--- tidy up
      deallocate (pad_mask)
      deallocate (reindex)

      return

   end subroutine paddedROI

   subroutine dimensionsOfROI(indx, nMax, left, right, top, bottom, xMax, yMax)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find the size of the ROI
      integer, dimension(:, :), intent(in)               ::      indx
      integer, intent(in)                              ::      nMax
      integer, dimension(:), intent(out)                ::      left, right, top, bottom               !   (1:nMax)
      integer, intent(out)                             ::      xMax, yMax

      integer                 ::      nx, ny
      integer                 ::      iy, ix, ii

      nx = size(indx, dim=1)            !   note that x- is row index
      ny = size(indx, dim=2)            !   and y- is column index. Otherwise it gets confusing...

      left(1:nMax) = nx
      right(1:nMax) = 1
      top(1:nMax) = 1
      bottom(1:nMax) = ny

      do iy = 1, ny
         do ix = 1, nx

            ii = indx(ix, iy)
            if (ii == 0) cycle
            left(ii) = min(left(ii), ix)
            right(ii) = max(right(ii), ix)
            top(ii) = max(top(ii), iy)
            bottom(ii) = min(bottom(ii), iy)

         end do
      end do

      xMax = 0
      yMax = 0
      do ii = 1, nMax
         xMax = max(xMax, right(ii) - left(ii) + 1)
         yMax = max(yMax, top(ii) - bottom(ii) + 1)
      end do
      if (LIBMAX_DBG) print *, "Lib_Maxima2d::dimensionsOfROI info - xmax,ymax = ", xmax, ymax
      return
   end subroutine dimensionsOfROI

   subroutine bxFilter(img, bx)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      flatten all pixels > bx over their neighbours
      real(kind=real64), dimension(:, :), intent(inout)  ::      img
      real(kind=real64), intent(in)                    ::      bx

      integer                 ::      nx, ny
      integer                 ::      iy, ix
      real(kind=real64)       ::      ff, fmax, fmin
      real(kind=real64), dimension(3, 3)        ::      img_plaque
      integer                 ::      bx_cut
      real(kind=real64), dimension(size(img, dim=1), size(img, dim=2))        ::      img_tmp

      nx = size(img, dim=1)            !   note that x- is row index
      ny = size(img, dim=2)            !   and y- is column index. Otherwise it gets confusing...
      img_tmp(1:nx, 1:ny) = img(1:nx, 1:ny)
      bx_cut = 0
      do iy = 2, ny - 1                  !   note 2: ny-1 because dont want maxima on boundary
         do ix = 2, nx - 1
            ff = img(ix, iy)

            img_plaque(1:3, 1:3) = img(ix - 1:ix + 1, iy - 1:iy + 1)
            img_plaque(2, 2) = 0.0d0

            fmax = maxval(img_plaque)       !   maximum pixel in 3x3 block centred on ix,iy
            fmin = minval(img_plaque)       !   maximum pixel in 3x3 block centred on ix,iy

            if (ff > fmax + bx) then
               !   pixel much brighter than neighbours. Most likely noise, so flatten
               ff = sum(img_plaque)/8
               img_tmp(ix, iy) = ff
               bx_cut = bx_cut + 1
            else if (ff < fmax - bx) then
               !   pixel much darker than neighbours. Most likely noise, so flatten
               ff = sum(img_plaque)/8
               img_tmp(ix, iy) = ff
               bx_cut = bx_cut + 1
            end if

         end do
      end do
      print *, "Lib_Maxima2d::bxFilter info - number of pixels flattened ", bx_cut
      img(1:nx, 1:ny) = img_tmp(1:nx, 1:ny)
      return
   end subroutine bxFilter

   subroutine polyFilter(img, bx)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      if pixel X is different by more than bx than its neighbours,
      !*      Least square fit to a quadratic polynomial fit to find the best value for X
      !*
      !*              1   2   3
      !*          4   5   6   7   8
      !*          9   10  X   11  12
      !*          13  14  15  16  17
      !*              18  19  20

      real(kind=real64), dimension(:, :), intent(inout)  ::      img
      real(kind=real64), intent(in)                    ::      bx

      integer                 ::      nx, ny
      integer                 ::      iy, ix, jx, jy, kx, ky, kk
      real(kind=real64)       ::      ff, fmax, fmin
      real(kind=real64), dimension(-2:2, -2:2)        ::      img_plaque
      integer                 ::      bx_cut
      real(kind=real64), dimension(size(img, dim=1), size(img, dim=2))        ::      img_tmp

      !---    flat weighting
      real(kind=real64), parameter     ::      C0 = -0.0530303d0
      real(kind=real64), parameter     ::      C1 = +0.0113636d0
      real(kind=real64), parameter     ::      C2 = +0.140152d0
      real(kind=real64), parameter     ::      C3 = +0.204545d0

      integer             ::      n0, n1, n2, n3
      real(kind=real64)   ::      m0, m1, m2, m3

      integer     ::      loop

      nx = size(img, dim=1)            !   note that x- is row index
      ny = size(img, dim=2)            !   and y- is column index. Otherwise it gets confusing...

      do loop = 1, 3

         img_tmp(1:nx, 1:ny) = img(1:nx, 1:ny)
         bx_cut = 0
         do iy = 1, ny
            do ix = 1, nx
               ff = img(ix, iy)

               img_plaque(-2:2, -2:2) = 0
               do jy = max(1, iy - 2), min(ny, iy + 2)
                  ky = jy - iy
                  do jx = max(1, ix - 2), min(nx, ix + 2)
                     kx = jx - ix
                     kk = kx*kx + ky*ky
                     if (kk*(8 - kk) > 0) img_plaque(kx, ky) = img(jx, jy)       !   ignore (0,0) and (2,2)
                  end do
               end do

               fmax = maxval(img_plaque(-1:1, -1:1))       !   maximum pixel in 3x3 block centred on ix,iy
               fmin = minval(img_plaque(-1:1, -1:1), mask=(img_plaque(-1:1, -1:1) > 0))       !   minimum pixel in 3x3 block centred on ix,iy

               if ((ff > fmax + bx) .or. (ff < fmin - bx)) then

                  !   Most likely noise, so flatten
                  n0 = count((/img_plaque(-1, -2), img_plaque(1, -2), img_plaque(-2, -1), img_plaque(2, -1) &
                               , img_plaque(-2, 1), img_plaque(2, 1), img_plaque(-1, 2), img_plaque(1, 2)/) > 0)
                  n1 = count((/img_plaque(0, -2), img_plaque(-2, 0), img_plaque(2, 0), img_plaque(-0, 2)/) > 0)
                  n2 = count((/img_plaque(-1, -1), img_plaque(1, -1), img_plaque(-1, 1), img_plaque(1, 1)/) > 0)
                  n3 = count((/img_plaque(0, -1), img_plaque(-1, 0), img_plaque(1, 0), img_plaque(0, 1)/) > 0)

                  m0 = sum((/img_plaque(-1, -2), img_plaque(1, -2), img_plaque(-2, -1), img_plaque(2, -1) &
                             , img_plaque(-2, 1), img_plaque(2, 1), img_plaque(-1, 2), img_plaque(1, 2)/))
                  m1 = sum((/img_plaque(0, -2), img_plaque(-2, 0), img_plaque(2, 0), img_plaque(0, 2)/))
                  m2 = sum((/img_plaque(-1, -1), img_plaque(1, -1), img_plaque(-1, 1), img_plaque(1, 1)/))
                  m3 = sum((/img_plaque(0, -1), img_plaque(-1, 0), img_plaque(1, 0), img_plaque(0, 1)/))

                  m0 = m0*8/max(n0, 1)
                  m1 = m1*4/max(n1, 1)
                  m2 = m2*4/max(n2, 1)
                  m3 = m3*4/max(n3, 1)

                  ff = C0*m0 + C1*m1 + C2*m2 + C3*m3
                  img_tmp(ix, iy) = ff
                  !                        print *,"Lib_Maxima2d::polyFilter info - change pixel ",ix,iy," from ",img(ix,iy)," to ",ff
                  bx_cut = bx_cut + 1

               end if

            end do
         end do
         print *, "Lib_Maxima2d::polyFilter info - number of pixels flattened ", bx_cut
         img(1:nx, 1:ny) = img_tmp(1:nx, 1:ny)

      end do
      return
   end subroutine polyFilter

   subroutine polyFilter2(img, b, sig)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      Least square fit to a quadratic polynomial fit to find the best value for X
      !*
      !*              1   2   3
      !*          4   5   6   7   8
      !*          9   10  X   11  12
      !*          13  14  15  16  17
      !*              18  19  20

      real(kind=real64), dimension(:, :), intent(inout)  ::      img
      real(kind=real64), intent(in)                    ::      b, sig

      integer                 ::      nx, ny
      integer                 ::      iy, ix, jx, jy, kx, ky, kk
      real(kind=real64)       ::      ff, gg, pp, zz
      real(kind=real64), dimension(-2:2, -2:2)        ::      img_plaque
      integer                 ::      bx_cut
      real(kind=real64), dimension(size(img, dim=1), size(img, dim=2))        ::      img_tmp

      !---    flat weighting
      real(kind=real64), parameter     ::      C0 = -0.0530303d0
      real(kind=real64), parameter     ::      C1 = +0.0113636d0
      real(kind=real64), parameter     ::      C2 = +0.140152d0
      real(kind=real64), parameter     ::      C3 = +0.204545d0

      integer             ::      n0, n1, n2, n3
      real(kind=real64)   ::      m0, m1, m2, m3

      integer     ::      loop

      nx = size(img, dim=1)            !   note that x- is row index
      ny = size(img, dim=2)            !   and y- is column index. Otherwise it gets confusing...

      zz = 1/sqrt(2*3.1415926535897932384626433832795d0*sig*sig)

      do loop = 1, 3

         img_tmp(1:nx, 1:ny) = img(1:nx, 1:ny)
         bx_cut = 0
         do iy = 1, ny
            do ix = 1, nx
               ff = img(ix, iy)

               img_plaque(-2:2, -2:2) = 0
               do jy = max(1, iy - 2), min(ny, iy + 2)
                  ky = jy - iy
                  do jx = max(1, ix - 2), min(nx, ix + 2)
                     kx = jx - ix
                     kk = kx*kx + ky*ky
                     if (kk*(8 - kk) > 0) img_plaque(kx, ky) = img(jx, jy)       !   ignore (0,0) and (2,2)
                  end do
               end do

               !---    find the expected level of this point
               n0 = count((/img_plaque(-1, -2), img_plaque(1, -2), img_plaque(-2, -1), img_plaque(2, -1) &
                            , img_plaque(-2, 1), img_plaque(2, 1), img_plaque(-1, 2), img_plaque(1, 2)/) > 0)
               n1 = count((/img_plaque(0, -2), img_plaque(-2, 0), img_plaque(2, 0), img_plaque(-0, 2)/) > 0)
               n2 = count((/img_plaque(-1, -1), img_plaque(1, -1), img_plaque(-1, 1), img_plaque(1, 1)/) > 0)
               n3 = count((/img_plaque(0, -1), img_plaque(-1, 0), img_plaque(1, 0), img_plaque(0, 1)/) > 0)

               m0 = sum((/img_plaque(-1, -2), img_plaque(1, -2), img_plaque(-2, -1), img_plaque(2, -1) &
                          , img_plaque(-2, 1), img_plaque(2, 1), img_plaque(-1, 2), img_plaque(1, 2)/))
               m1 = sum((/img_plaque(0, -2), img_plaque(-2, 0), img_plaque(2, 0), img_plaque(0, 2)/))
               m2 = sum((/img_plaque(-1, -1), img_plaque(1, -1), img_plaque(-1, 1), img_plaque(1, 1)/))
               m3 = sum((/img_plaque(0, -1), img_plaque(-1, 0), img_plaque(1, 0), img_plaque(0, 1)/))

               m0 = m0*8/max(n0, 1)
               m1 = m1*4/max(n1, 1)
               m2 = m2*4/max(n2, 1)
               m3 = m3*4/max(n3, 1)

               gg = C0*m0 + C1*m1 + C2*m2 + C3*m3

               !---    how probable is ff?
               pp = exp(-(gg - ff)*(gg - ff)/(2*sig*sig))

               !---    construct a new filtered point
               gg = ff*pp + gg*(1 - pp)

               if ((ff - b)*(gg - b) < 0) then
                  !   this is strange - the modification wants to make a dark pixel bright or vice versa. Cap it.
                  gg = b
               end if

               img_tmp(ix, iy) = gg

            end do
         end do
         print *, "Lib_Maxima2d::polyFilter2 info - number of pixels flattened ", bx_cut
         img(1:nx, 1:ny) = img_tmp(1:nx, 1:ny)

      end do
      return
   end subroutine polyFilter2

   subroutine aboveThreshold(img, f_c, indx, nMax)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find all pixels which are above threshold

      real(kind=real64), dimension(:, :), intent(in)     ::      img
      real(kind=real64), intent(in)                    ::      f_c
      integer, dimension(:, :), intent(out)              ::      indx
      integer, intent(out)                             ::      nMax

      integer                 ::      nx, ny
      integer                 ::      iy, ix
      real(kind=real64)       ::      ff

      nx = size(img, dim=1)            !   note that x- is row index
      ny = size(img, dim=2)            !   and y- is column index. Otherwise it gets confusing...

      !---    first step: each maxmimum pixel is identified and labelled.
      nMax = 0
      indx = 0
      do iy = 1, ny
         do ix = 1, nx
            ff = img(ix, iy)
            if (ff >= f_c) then
               nMax = nMax + 1
               indx(ix, iy) = nMax
            end if
         end do
      end do

      if (LIBMAX_DBG) then
        print *, "Lib_Maxima2d::aboveThreshold info - number of pixels above threshold ", nMax, " ", real(nMax, kind=real64)/(nx*ny)
         call screendump(img, indx)
         print *, ""
      end if

      return
   end subroutine aboveThreshold

   subroutine rawMaxima2(img, f_c, indx, nMax, gte)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find all pixels which are a maximum -
      !*          1   2   3
      !*          4   5   6
      !*          7   8   9
      !*      ie those points where f(5)>=f(1:9) and f(5)>=f_c
      !*      mark each with an index point indx(i,j) = 1:nmax

      real(kind=real64), dimension(:, :), intent(in)     ::      img
      real(kind=real64), intent(in)                    ::      f_c
      integer, dimension(:, :), intent(out)              ::      indx
      integer, intent(out)                             ::      nMax
      logical, intent(in)                              ::      gte
      integer                 ::      nx, ny
      integer                 ::      iy, ix
      real(kind=real64)       ::      ff, fmax
      real(kind=real64), dimension(3, 3)        ::      img_plaque

      real(kind=real64)       ::      DELTA

      real(kind=real64), dimension(:, :), allocatable      ::      img_buffer

      nx = size(img, dim=1)            !   note that x- is row index
      ny = size(img, dim=2)            !   and y- is column index. Otherwise it gets confusing...

      DELTA = 1.0d0/(256*256)
      if (gte) DELTA = 0.0d0

      !---    allow maxima on the boundaries
      allocate (img_buffer(0:nx + 1, 0:ny + 1))
      img_buffer = 0
      img_buffer(1:nx, 1:ny) = img(1:nx, 1:ny)

      !---    first step: each maximum pixel is identified and labelled.
      nMax = 0
      indx = 0
      do iy = 1, ny
         do ix = 1, nx
            ff = img_buffer(ix, iy)
            if (ff < f_c) cycle

            img_plaque(1:3, 1:3) = img_buffer(ix - 1:ix + 1, iy - 1:iy + 1)
            img_plaque(2, 2) = 0.0d0

            fmax = min(1.0d0 - 2*DELTA, maxval(img_plaque))      !   maximum pixel in 3x3 block centred on ix,iy

            if (ff >= fmax + DELTA) then
               !   max pixel is ix,iy
               if (fmax > f_c) then
                  !   there is also a neighbour over threshold
                  nMax = nMax + 1
                  indx(ix, iy) = nMax
               end if
            end if

         end do
      end do

      if (LIBMAX_DBG) then
         print *, "raw max"
         call screendump(img, indx)
         print *, ""

      end if

      return
   end subroutine rawMaxima2

   subroutine findMaxima(img, b, f_c, indx, nMax, gte)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      scan the 2d image
      !*      find the number of maxima over the level f_c
      !*      b is the background level used to ascribe weights correctly
      !*      index each point as belonging to a single maximum indx(i,j) = 1:nMax

      real(kind=real64), dimension(:, :), intent(in)     ::      img
      real(kind=real64), intent(in)                    ::      f_c, b
      integer, dimension(:, :), intent(out)              ::      indx
      integer, intent(out)                             ::      nMax
      logical, intent(in), optional                     ::      gte

      integer                 ::      nx, ny, nflood
      integer                 ::      iy, ix, jy, jx, ii, jj, kx, ky
      real(kind=real64)       ::      ff, fmax
      logical                 ::      done, ok
      real(kind=real64), dimension(:), allocatable      ::      ww
      integer, dimension(:), allocatable                ::      reindex
      integer, dimension(size(indx, dim=1), size(indx, dim=2))              ::      indx_tmp

      nx = size(img, dim=1)            !   note that x- is row index
      ny = size(img, dim=2)            !   and y- is column index. Otherwise it gets confusing...

      if (LIBMAX_DBG) print *, "f_c = ", f_c, int(f_c*255)

      !---    first step: each maxmimum pixel is identified and labelled.
      if (present(gte)) then
         call rawMaxima2(img, f_c, indx, nMax, gte)
      else
         call rawMaxima2(img, f_c, indx, nMax, .false.)
      end if

      !---    second step: pair up any maxima within a short range -
      !           X Y Y
      !         X X Y Y Y
      !         X X 0 Y Y
      !         X X Y Y Y
      !           X Y Y

      do nflood = 1, max(nx, ny)
         done = .true.
         do iy = 3, ny - 2          !   note -3 because checking to iy+2,ix+2, no maxima on boundary
            do ix = 1, nx - 2
               ii = indx(ix, iy)
               if (ii == 0) cycle
               do jy = -2, 2
                  do jx = 0, 2
                     jj = indx(ix + jx, iy + jy)
                     if (jj*(jj - ii) /= 0) then
                        ii = min(ii, jj)
                        indx(ix + jx, iy + jy) = ii
                        indx(ix, iy) = ii
                        done = .false.
                     end if
                  end do
               end do
            end do
         end do
         if (done) exit
      end do

      if (LIBMAX_DBG) then
         print *, "close maxima joined"
         call screendump(img, indx, txt=.true.)
         print *, ""
      end if

      if (maxval(indx) == 0) return     !   no maxima

      !---    third step: flood fill
      do nflood = 1, max(nx, ny)/2
         done = .true.
         indx_tmp = indx
         do iy = 1, ny
            do ix = 1, nx

               ii = indx(ix, iy)
               if (ii /= 0) cycle                                !   only continue if this point not ascribed yet
               ff = img(ix, iy)

               if (ff*(1 - ff) > 0) then                           !   dont try to ascribe dead pixels
                  fmax = 0.0d0                                !   find highest value neighbour which is ascribed

                  do jy = max(1, iy - 1), min(ny, iy + 1)
                     do jx = max(1, ix - 1), min(nx, ix + 1)

!                                    if ( (jy-iy)*(jy-iy) + (jx-ix)*(jx-ix) == 2 ) cycle     !   ignore (1,1) diagonal connection
                        jj = indx(jx, jy)
                        if (jj == 0) cycle

                        if (img(jx, jy) >= fmax) then
                           fmax = img(jx, jy)
                           ii = jj
                           kx = jy; ky = jx
                        end if

                     end do
                  end do

                  if (ff < f_c) indx_tmp(ix, iy) = ii                        !   if this point is under threshold, just set to index of max neighbour
                  if (fmax > ff) indx_tmp(ix, iy) = ii                       !   if neighbour is brighter, set to index of brighter neighbour
                  ok = ((iy - 1)*(ny - iy)*(ix - 1)*(nx - ix) == 0)         !   ok if on border - if a pixel _is_ a maximum on the border just quietly drop it
                  if (ok) indx_tmp(ix, iy) = ii
                  ok = ok .or. (indx_tmp(ix, iy) > 0)                        !   ok if set to region
                  done = done .and. ok

               end if
            end do
         end do

         indx = indx_tmp

         if (done) exit
      end do

      if (LIBMAX_DBG) then
         print *, "flood fill"
         call screendump(img, indx, txt=.true.)
         print *, ""
      end if

      !---    fourth step: find the weights of each region w_i = sum f
      !       note that some regions can have weight < 0
      allocate (ww(0:nMax)); ww = 0.0d0       !   allow 0: for easier algorithm where there are a lot of dead pixels
      allocate (reindex(0:nMax)); reindex = 0
      do iy = 1, ny
         do ix = 1, nx
            ii = indx(ix, iy)
            ff = (img(ix, iy) - b)
            ww(ii) = ww(ii) + ff
         end do
      end do

      !---    final step: reindex sequentially so that 1 is the highest weight region
      nMax = 0
      do
         ii = maxloc(ww(1:), dim=1)
         if (ww(ii) <= -huge(1.0)/2) exit        !   algorithm done - no positive weight regions remain.
         nMax = nMax + 1                         !   have found a new region
         reindex(ii) = nMax                      !   all pixels indexed ii should read "nmax"
         ww(ii) = -huge(1.0)
      end do
      do iy = 1, ny
         do ix = 1, nx
            ii = indx(ix, iy)
            indx(ix, iy) = reindex(ii)
         end do
      end do

      if (LIBMAX_DBG) then
         print *, "ordered"
         call screendump(img, indx, txt=.true.)
!!                stop
      end if
      deallocate (ww)
      deallocate (reindex)
      return

   end subroutine findMaxima

   subroutine mergeMaxima(img, b, d, indx, nMax)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      merge any two peaks where there exists a path which does not drop below
      !*      a fraction d of the linear interpolated value
      !*      ie you find a point x with distance r1 from peak f1
      !*      and distance r2 from peak f2
      !*      then x should have value fx = (f1 r2 + f2 r1)/(r1+r2)
      !*      merge if f < d(fx - b)

      real(kind=real64), dimension(:, :), intent(in)     ::      img
      integer, dimension(:, :), intent(inout)            ::      indx
      real(kind=real64), intent(in)                    ::      b, d
      integer, intent(inout)                           ::      nMax

      integer                 ::      nx, ny
      integer                 ::      iy, ix, jy, jx
      integer                 ::      ii, jj, kk, nn

      real(kind=real64), dimension(0:nMax)         ::      pp          !   height for each peak
      real(kind=real64), dimension(2, 0:nMax)       ::      xx          !   position of each peak
      real(kind=real64), dimension(0:nMax)         ::      ww          !   weights of maxima w_i = sum f
      integer, dimension(0:nMax)                   ::      reindex
      integer, dimension(:, :), allocatable              ::      connex
      real(kind=real64), dimension(:, :), allocatable    ::      saddle          !   height of saddle
      real(kind=real64), dimension(:, :), allocatable    ::      saddleFrac      !   height of saddle as fraction of the expected linear interpolation
      real(kind=real64), dimension(:, :, :), allocatable  ::      saddlex          !   position of each saddle
      logical                 ::      ok, done
      real(kind=real64)       ::      fi, fj, ri, rj, fx, dx, dy
      integer, dimension(2)    ::      maxSaddle

      allocate (connex(0:nMax, 0:nMax))
      allocate (saddle(nMax, 0:nMax))
      allocate (saddleFrac(nMax, nMax))
      allocate (saddlex(2, nMax, 0:nMax))

      nx = size(img, dim=1)
      ny = size(img, dim=2)

      !---    find weights and peak heights. also find connectivity

      do
         done = .true.

         ww = 0.0d0
         pp = 0.0d0
         xx = 0.0d0
         connex = 0
         saddle = 0.0d0
         do iy = 1, ny
            do ix = 1, nx

               ii = indx(ix, iy)
               if (ii == 0) cycle            !   dead pixel

               fi = img(ix, iy)
               ww(ii) = ww(ii) + (fi - b)
               if (fi > pp(ii)) then
                  xx(1:2, ii) = (/ix, iy/)
                  pp(ii) = fi
               end if

               !---    check if this point is on an interface
               nn = connex(0, ii)
               do jy = max(1, iy - 1), min(ny, iy + 1)
                  do jx = max(1, ix - 1), min(nx, ix + 1)
                     jj = indx(jx, jy)
                     if (jj > ii) then
                        !   yes, this point is on the interface between i and j. But is it the saddle point?
                        fj = img(jx, jy)
                        ok = .false.
                        do kk = 1, nn
                           if (connex(kk, ii) == jj) then
                              ok = .true.
                              saddle(kk, ii) = max(saddle(kk, ii), (fi + fj)/2)
                              saddlex(1:2, kk, ii) = 0.50d0*(/(ix + jx), (iy + jy)/)
                           end if
                        end do
                        if (.not. ok) then
                           !   this is a new interface not seen before
                           nn = nn + 1
                           connex(0, ii) = nn
                           connex(nn, ii) = jj
                           saddle(kk, ii) = max(saddle(kk, ii), (fi + fj)/2)
                           saddlex(1:2, kk, ii) = 0.50d0*(/(ix + jx), (iy + jy)/)
                        end if
                     end if
                  end do
               end do
            end do
         end do

         saddleFrac = 0.0d0
         do ii = 1, nMax
            nn = connex(0, ii)                               !   for each peak i, look at saddles
            fi = pp(ii)
            do kk = 1, nn
               jj = connex(kk, ii)
               fj = pp(jj)

               dx = saddlex(1, kk, ii) - xx(1, ii)            !   find distance from peak i to saddle ...
               dy = saddlex(2, kk, ii) - xx(2, ii)
               ri = sqrt(real(dx*dx + dy*dy, kind=real64))

               dx = saddlex(1, kk, ii) - xx(1, jj)            !   ... and distance from peak j to saddle
               dy = saddlex(2, kk, ii) - xx(2, jj)
               rj = sqrt(real(dx*dx + dy*dy, kind=real64))

               fx = (fi*rj + fj*ri)/(ri + rj) - b            !   expected height above background assuming linear model

               saddleFrac(kk, ii) = (saddle(kk, ii) - b)/fx

            end do
         end do

         !---     code revision: connect the highest fraction peak
         fx = maxval(saddleFrac)   !   this is the highest saddle, so the best to connect
         if (fx > d) then
            !   the highest saddle is over the threshold for recombination
            maxSaddle = maxloc(saddleFrac)
            ii = maxSaddle(2)
            jj = connex(maxSaddle(1), ii)
            !---    peak j is connected to peak i
            if (LIBMAX_DBG) print *, "connect peaks ", maxSaddle, " ", ii, " and ", jj, " max saddle frac ", fx!,saddle(kk,ii)-b
            ww(ii) = ww(ii) + ww(jj)
            ww(jj) = -huge(1.0)
            done = .false.
            where (indx == jj)
               indx = ii
            end where
         end if

         if (LIBMAX_DBG) print *, ""

         !---    final step: reindex sequentially so that 1 is the highest weight region
         nMax = 0; reindex = 0
         do
            ii = maxloc(ww(1:), dim=1)
            if (ww(ii) == 0) then
               ww(ii) = -huge(1.0)
               cycle
            end if
            if (ww(ii) <= -huge(1.0)/2) exit       !   algorithm done
            nMax = nMax + 1             !   have found a new region
            reindex(ii) = nMax          !   all pixels indexed ii should read "nmax"
            ww(ii) = -huge(1.0)

         end do
         do iy = 1, ny
            do ix = 1, nx
               ii = indx(ix, iy)
               indx(ix, iy) = reindex(ii)
            end do
         end do

         if (done) exit
      end do

      if (LIBMAX_DBG) then
         print *, "merged"
         call screendump(img, indx, txt=.true.)
         print *, ""
         print *, ""
         !stop
      end if
      return
   end subroutine mergeMaxima

   subroutine screendump(img, indx, txt)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(:, :), intent(in)     ::      img
      integer, dimension(:, :), intent(in)               ::      indx
      logical, intent(in), optional                     ::      txt
      integer         ::      iy, ix, nx, ny
      logical         ::      asTxt
            character(len=2),dimension(0:16),parameter      ::      letter = (/ "  ","@@","##","**","++","==","--","::",".."," @"," #"," *"," +"," ="," -"," :"," ." /)

      nx = size(img, dim=1)        !   note nx is number of rows
      ny = size(img, dim=2)        !   and ny is number of columns

      asTxt = .false.
      if (present(txt)) asTxt = txt

      if (.not. asTxt) then
         write (*, fmt='(a3)', advance="no") ""
         do ix = 1, nx
            write (*, fmt='(i3)', advance="no") ix
         end do
         write (*, fmt='(a)', advance="yes") ""
      else
         write (*, fmt='(a)') repeat(".", (2*nx + 2))
      end if

      do iy = 1, ny
         if (asTxt) then
            write (*, fmt='(a)', advance="no") "."
            do ix = 1, nx
               if (img(ix, iy) == 0) then
                  write (*, fmt='(a2)', advance="no") "  "
               else if (indx(ix, iy) < size(letter)) then
                  write (*, fmt='(a2)', advance="no") letter(indx(ix, iy))
               else if (indx(ix, iy) < 100) then
                  write (*, fmt='(i2)', advance="no") indx(ix, iy)
               else
                  write (*, fmt='(a2)', advance="no") "XX"
               end if
            end do
            write (*, fmt='(a)', advance="no") "."
         else
            write (*, fmt='(i3)', advance="no") iy
            do ix = 1, nx
               if (indx(ix, iy) == 0) then
                  write (*, fmt='(a3)', advance="no") "   "
               else
                  write (*, fmt='(i3)', advance="no") indx(ix, iy)
               end if
            end do
         end if
         write (*, fmt='(a)', advance="yes") ""
      end do
      if (asTxt) write (*, fmt='(a)') repeat(".", (2*nx + 2))
      return
   end subroutine screendump

   subroutine findSigma(img, maxPix, s, b, t, f)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      simple initialisation of background standard deviation using Ridler & Calvard criterion
      real(kind=real64), dimension(:, :), intent(in)         ::      img
      real(kind=real64), intent(in)                        ::      maxPix
      real(kind=real64), intent(out)                       ::      s, b, t, f

      integer             ::      ix, iy, nx, ny
      real(kind=real64)   ::      gg, ksum
      real(kind=real64)   ::      s2

      nx = size(img, dim=1)
      ny = size(img, dim=2)

      call findRidlerAndCalvardThresh(img, maxPix, b, t, f)

      ksum = 0.0d0; s = 0.0d0; s2 = 0.0d0
      do iy = 1, ny
         do ix = 1, nx

            gg = img(ix, iy)
            if (gg*(t - gg) > 0) then
               s = s + gg
               s2 = s2 + gg*gg
               ksum = ksum + 1.0d0
            end if

         end do
      end do

      if (ksum > 0) then
         s = s/ksum
         s2 = s2/ksum
         s = sqrt(max(0.0d0, s2 - s*s))

      end if

      print *, "Lib_Maxima2d::findSigma info - b = ", b, "~", int(b*255), " in ", count(img < t), " px"
      print *, "Lib_Maxima2d::findSigma info - t = ", t, "~", int(t*255)
      print *, "Lib_Maxima2d::findSigma info - f = ", f, "~", int(f*255), " in ", count(img >= t), " px"
      print *, "Lib_Maxima2d::findSigma info - s = ", s, "~", int(s*255)
      return
   end subroutine findSigma

   pure subroutine intensityShift0(g, b, b0, t, t0, f, f0, g_out)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !      then use affine shifts to push 0->0, b->b0, t->t0, f->f0, 1->1
      !       return the new intensity level g
      !
      real(kind=real64), intent(in)            ::      g
      real(kind=real64), intent(in)            ::      b, b0, t, t0, f, f0
      real(kind=real64), intent(out)           ::      g_out

      if (g <= 0.0) then
         g_out = 0.0d0
      else if (g < b) then
         g_out = b0*g/b
      else if (g < t) then
         g_out = b0 + (g - b)*(t0 - b0)/(t - b)
      else if (g < f) then
         g_out = t0 + (g - t)*(f0 - t0)/(f - t)
      else if (g < 1) then
         g_out = f0 + (g - f)*(1.0 - f0)/(1 - f)
      else
         g_out = 1.0d0
      end if
!
      return
   end subroutine intensityShift0

   subroutine intensityShift(img, maxPix, b0, t0, f0)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !       shift intensity of img_in so that the levels are those stored
      real(kind=real64), dimension(:, :), intent(inout)  ::      img
      real(kind=real64), intent(in)                    ::      maxPix
      real(kind=real64), intent(in)                    ::      b0, t0, f0

      real(kind=real64), dimension(0:255)              ::      hist
      real(kind=real64)                               ::      bb, tt, ff, gg
      integer             ::      ix, iy, loop
      integer             ::      nx, ny

      nx = size(img, dim=1)
      ny = size(img, dim=2)

      do loop = 1, 4

         !---    find new intensity histogram
         call findHist(img, hist)
         call findRidlerAndCalvardThresh(hist, maxPix, bb, tt, ff)
         print *, "Lib_Maxima2d::intensityShift"
         write (*, fmt='(a,3f12.4)') "old back, thresh, fore ", b0, t0, f0
         write (*, fmt='(a,3f12.4)') "new back, thresh, fore ", bb, tt, ff

         !---    shift intensities
         do iy = 1, ny
            do ix = 1, nx
               call intensityShift0(img(ix, iy), bb, b0, tt, t0, ff, f0, gg)
               img(ix, iy) = gg
            end do
         end do

      end do

      return
   end subroutine intensityShift

end module Lib_Maxima2d
