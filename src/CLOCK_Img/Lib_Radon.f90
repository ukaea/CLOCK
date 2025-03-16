
module Lib_Radon
!---^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_Radon from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      A module to work with Radon transformations of 2d image data
!*      Basic operation is shown in testRadon.f90
!*
!*          call cutBuffered_img( img, left,bottom,right,top , buffered_img [,px_pad] )                             !   cuts a rectangle out of img(0:Nx-1,0:Ny-1) and allocates buffered_img. Left should be in range [0:Nx-1]
!*
!*          call construct_sinogram(buffered_img,sinogram [,dt_pad,px_pad])                                         !   allocates and generates the sinogram
!*
!*          call findLinesFromSinogram( sinogram, (right+1-left),(top+1-bottom), nLine,line [,dt_pad,px_pad] )      !   allocates and constructs lines with ends 0:Mx,0:My, where Mx = right+1-left
!*
!*          call reconstructImage( line(:,1:nLine),img_out )                                                        !   reconstructs the image
!*

!*

   use iso_fortran_env
   use Lib_png
   use Lib_RidlerCalvard
   use Lib_ColourScale

   implicit none
   private

   public              ::      cutBuffered_img
   public              ::      construct_sinogram
   public              ::      findLinesFromSinogram
   public              ::      reconstructImage
   public              ::      reconstructImageRGB !check if there are improvements in reconstructImage that can be
   !ported to this. Remove debug functions when integrated.

   real(kind=real64), parameter                     ::      PI = 3.141592654d0

   interface reconstructImage
      module procedure reconstructImage0
      module procedure reconstructImage1
   end interface

   interface findLinesFromSinogram
      module procedure findLinesFromSinogram0
      module procedure findLinesFromSinogram1
      module procedure findLinesFromSinogram2
   end interface

contains
!---^^^^^^^^

   subroutine cutBuffered_img(img, x1, y1, x2, y2, buffered_img, bimg_pad)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      cut out a section of img from x1:x2 , y1:y2
      !*      and add a buffer of pad pixels all round the edge
      real(kind=real64), dimension(0:, 0:), intent(in)               ::      img                 !   input image
      integer, intent(in)                                          ::      x1, x2, y1, y2         !   extent
      real(kind=real64), dimension(:, :), allocatable, intent(out)    ::      buffered_img        !   cropped image including padding
      integer, intent(in), optional                                 ::      bimg_pad                 !   padding

      integer                         ::      Nx, Ny   !   input image x&y dimensions
      integer                         ::      Mx, My   !   cropped image dimensions
      integer                         ::      ix, iy   !   pixel index
      integer                         ::      px_pad  !   pixel padding

      !---    size of input image
      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      !---    size of cropped image
      Mx = x2 - x1 + 1
      My = y2 - y1 + 1

      !---    size of padding
      px_pad = 2; if (present(bimg_pad)) px_pad = bimg_pad

      !---    set border of buffered image to zero
      allocate (buffered_img(-px_pad:Mx - 1 + px_pad, -px_pad:My - 1 + px_pad))
      buffered_img(:-1, :) = 0
      buffered_img(Mx:, :) = 0
      buffered_img(:, :-1) = 0
      buffered_img(:, My:) = 0

      !---    copy pixels from input image
      do iy = max(-px_pad, -y1), min(My - 1 + px_pad, Ny - 1 - y1)
         do ix = max(-px_pad, -x1), min(Mx - 1 + px_pad, Nx - 1 - x1)
            buffered_img(ix, iy) = img(x1 + ix, y1 + iy)
         end do
      end do

      return
   end subroutine cutBuffered_img

   subroutine construct_sinogram(img, sinogram, dtheta_pad, bimg_pad)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given an input image    (0:Nx-1,0:Ny-1)
      !*      construct a sinogram    (-rr:rr,-dtheta_pad:nTheta-1+dtheta_pad)        - using an extra dtheta_pad buffer of +2 dtheta intervals allows me to use a 25 point stencil
      !*

      real(kind=real64), dimension(0:, 0:), intent(in)               ::      img             !   image array in
      real(kind=real64), dimension(:, :), allocatable, intent(out)    ::      sinogram        !   sinogram array out
      integer, intent(in), optional                                 ::      dtheta_pad      !   padding for additional orientations in sinogram
      integer, intent(in), optional                                 ::      bimg_pad        !   padding used in buffered image - only needed for debugging
      integer                 ::      Nx, Ny                       !Buffered iamge x&y dimensions

      integer                 ::      nTheta                      !number of angles to ratoate image
      real(kind=real64)       ::      aa, bb, theta, dtheta, ww       !See below
      real(kind=real64)       ::      cost, sint, qx, qy
      integer                 ::      dt_pad

      real(kind=real64), dimension(2)  ::      q1, q2   !   pixel corners for linear interpolation.

      logical                 ::      ok              !   does the radon beam intercept the image?
      integer                 ::      nn              !   Rotation index
      integer                 ::      ix, ndivx        !   integrate x from 0:ndivx
      integer                 ::      iy, rr           !   find integral at y positions -rr:rr

      !---    find the size of the image
      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      !---    find the size of the sinogram. I want to rotate the image though 180�
      !       at each angular step, I will integrate across the image.
      !       I will therefore allocate the sinogram to just contain the rotating image,
      !       with a theta division to be 1 px height over the length
      !               _______
      !              |      /|
      !              |   r / |
      !              |    /  |
      !           2b |   x   |
      !              |       |
      !              |       |
      !              |_______|
      !                  2a
      aa = Nx*0.5d0                                   !   box half side
      bb = Ny*0.5d0
      rr = ceiling(sqrt(aa*aa + bb*bb))           !   radius of circumcircle which just contains the rotated image.
      if (present(bimg_pad)) then
         theta = ceiling(sqrt((aa - bimg_pad)*(aa - bimg_pad) + (bb - bimg_pad)*(bb - bimg_pad)))        !   gives exact dtheta
         theta = atan2(2.0d0, theta)                                                                !   independent of buffered image padding.
      else
         theta = atan2(2.0, real(rr))                   !   1 px over length of box = 2 px over half length of box
      end if
      nTheta = ceiling(PI/theta)
      dtheta = PI/nTheta
      dt_pad = 2; if (present(dtheta_pad)) dt_pad = dtheta_pad
      allocate (sinogram(-rr:rr, -dt_pad:nTheta - 1 + dt_pad))
      !print *,"Lib_Radon::construct_sinogram rr,dt_pad = ",rr,dt_pad," size sino ",size(sinogram,dim=1),size(sinogram,dim=2)," ntheta,dtheta ",ntheta,dtheta
      sinogram = 0
      !---

      !---    compute the sinogram
      do nn = -dt_pad, nTheta - 1 + dt_pad
         theta = nn*dtheta
         cost = cos(theta); sint = sin(theta)
         do iy = -rr, rr

            call boxIntercepts(aa, bb, real(iy, kind=real64), cost, sint, ok, q1, q2)

            if (ok) then
               ndivx = max(1, ceiling(norm2(q2 - q1)))         !   number of pixels = length of vector

               q2 = (q2 - q1)/ndivx                             !   (vector) step length
               q1 = q1 + (/aa, bb/)                            !  remember that my image is defined (0:Nx-1), but boxIntercepts returns (-a:a)
               ww = 0.0d0

               do ix = 0, ndivx

                  qx = q1(1) + q2(1)*ix
                  qy = q1(2) + q2(2)*ix

                  ww = ww + linint(img, qx, qy)

               end do
               sinogram(iy, nn) = ww/(ndivx + 1)
            end if
         end do

      end do
      !---

      return
   end subroutine construct_sinogram

   subroutine findLinesFromSinogram0(sinogram, Nx, Ny, nLine, line)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find spots in the sinogram. From these reconstruct lines in the image
      !*      line(1:6) = (/x1,y1,x2-x1,y2-y1,r,i/)
      !*      Note that the extent of the sinogram is determined by the longest diagonal of the original image

      real(kind=real64), dimension(0:, 0:), intent(in)               ::      sinogram    !   Input sinogram (-rr:rr,0:nTheta-1)
      integer, intent(in)                                          ::      Nx, Ny       !   original image dimensions ie img(0:Nx-1,0:Ny-1)
      integer, intent(out)                                         ::      nLine       !   number of detected lines
      real(kind=real64), dimension(:, :), pointer, intent(out)        ::      line        !   pointer to array of line objects
      call findLinesFromSinogram2(sinogram, Nx, Ny, nLine, line, r=int((size(sinogram, dim=1) - 1)/2), px_pad=2, dt_pad=2)
      return
   end subroutine findLinesFromSinogram0

   subroutine findLinesFromSinogram1(sinogram, Nx, Ny, nLine, line, px_pad, dt_pad)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find spots in the sinogram. From these reconstruct lines in the image
      !*      line(1:6) = (/x1,y1,x2-x1,y2-y1,r,i/)
      !*      Note that the extent of the sinogram is determined by the longest diagonal of the original image
      real(kind=real64), dimension(0:, 0:), intent(in)               ::      sinogram    !   Input sinogram (-rr:rr,0:nTheta-1)
      integer, intent(in)                                          ::      Nx, Ny       !   original image dimensions ie img(0:Nx-1,0:Ny-1)
      integer, intent(out)                                         ::      nLine       !   number of detected lines
      real(kind=real64), dimension(:, :), pointer, intent(out)        ::      line        !   pointer to array of line objects
      integer, intent(in)                                          ::      px_pad      !   padding used on buffered image buffered_img(-pad:Nx-1+pad,-pad:Ny-1+pad). Going to assume pad>=2.
      integer, intent(in)                                          ::      dt_pad      !   padding used on sinogram Going to assume pad>=2.
      call findLinesFromSinogram2(sinogram, Nx, Ny, nLine, line, int((size(sinogram, dim=1) - 1)/2), px_pad, dt_pad)
      return
   end subroutine findLinesFromSinogram1

   subroutine findLinesFromSinogram2(sinogram, Nx, Ny, nLine, line, r, px_pad, dt_pad)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find spots in the sinogram. From these reconstruct lines in the image
      !*      line(1:6) = (/x1,y1,x2-x1,y2-y1,r,i/)
      !*      Note that the extent of the sinogram is determined by the longest diagonal of the original image

      !*      so need to input original image extent (0:Nx-1,0:Ny-1) to place the lines at the borders of the image
      integer, intent(in)                                          ::      r           !   half width of sinogram
      integer, intent(in)                                          ::      dt_pad      !   padding used on sinogram Going to assume pad>=2.
      real(kind=real64), dimension(-r:, -dt_pad:), intent(in)        ::      sinogram    !   Input sinogram (-rr:rr,0:nTheta-1)
      integer, intent(in)                                          ::      Nx, Ny       !   original image dimensions ie img(0:Nx-1,0:Ny-1)
      integer, intent(out)                                         ::      nLine       !   number of detected lines
      real(kind=real64), dimension(:, :), pointer, intent(out)        ::      line        !   pointer to array of line objects
      integer, intent(in)                                          ::      px_pad      !   padding used on buffered image buffered_img(-pad:Nx-1+pad,-pad:Ny-1+pad). Going to assume pad>=2.

      integer             ::      nTheta                      !   number of angular divisions between 0:Pi in sinogram
      real(kind=real64)   ::      dtheta, theta, tmax, tmin    !change in angle per rotation, angle (both /rad), maximum gaussian blur, max & min number of gaussian blurs
      real(kind=real64)   ::      qx, qy, ww, mm, cc            !Qx &qy: detected line x & y postion from blob detection
      !ww: distance of blob from centre of rotated image
      !mm & cc: gradient and constant from y=mx+c
      logical             ::      ymxc                        !is the line in the form y=mx+c, else x=my+c

      !---    blob detection
      integer                                         ::      nt              !   number of gaussian blurs rounded to nearest integer
      real(kind=real64), dimension(:), allocatable      ::      tt              !   array of gaussian blurs
      integer                                         ::      nBlob           !   number of bright blobs in sinogram

      integer         ::      ii, jj, nBlobMax
      real(kind=real64), dimension(:, :), pointer    ::  bb, bb_tmp       !Pointer to object containing blob info, smaller temporary one.
      real(kind=real64), dimension(-2:2, -2:2)      ::  ff              !Sub section of sinogram
      logical                                     ::  ismax           !is the detected blob a maximum
      real(kind=real64)                           ::  f0, fx, fy, fwhm   !intensity, x & y coords, FWHM of blob maxima
      real(kind=real64), dimension(2)  ::      p1, p2                   !two element arrays to hold points

      logical         ::      ok

      real(kind=real64)           ::      randc_maxpx, randc_b, randc_t, randc_f

      nLine = 0
      nullify (line)

      nTheta = size(sinogram, dim=2) - 2*dt_pad
      dtheta = PI/nTheta

      !print *,"Lib_Radon::findLinesFromSinogram r,dt_pad,px_pad = ",r,dt_pad,px_pad," size sino ",size(sinogram,dim=1),size(sinogram,dim=2)," ntheta,dtheta ",ntheta,dtheta

      !---    blob detection
      !   first job is to decide whether a blob is significant: look at the Ridler& Calvard threshold criterion
      randc_maxpx = 0.25d0                             !   maximum proportion of pixels I'm prepared to consider "foreground"
      call findRidlerAndCalvardThresh(sinogram, randc_maxpx, randc_b, randc_t, randc_f)

      !   make an initial guess of 100 spots. If needed, we can reallocate this later.
      nBlobMax = 100
      allocate (bb(4, nBlobMax))
      nBlob = 0

      do jj = 0, nTheta - 1                                  !   note that I am only looking in the region 0<=theta<Pi, even though the sinogram is defined outside this range.

         do ii = 2 - r, r - 2                                 !   note: leaving a pad of pixels so I can use my 25 point stencil.

            !---    find a blob at this point
            ff = sinogram(ii - 2:ii + 2, jj - 2:jj + 2)
            call twentyfivePointStencil(ff, isMax, fx, fy, fwhm, f0)

            !---    do we accept the blob as a maximum
            isMax = isMax .and. (fwhm < r/2) &                   !   check width of maximum is not insane, ie > 1/4 the image
                    .and. (f0 > randc_t)                          !   check maximum is greater than R&C threshold

            !---    store the blob
            if (isMax) then
               if (nBlob == nBlobMax) then
                  !   need to reallocate the blob storage
                  allocate (bb_tmp(4, nBlobMax*2))
                  bb_tmp(1:4, 1:nBlobMax) = bb(1:4, 1:nBlobMax)
                  deallocate (bb)
                  bb => bb_tmp
                  nBlobMax = nBlobMax*2
               end if
               nBlob = nBlob + 1
               bb(1:4, nBlob) = (/ii + fx, jj + fy, fwhm/2, f0/)

               !print *,"Lib_Radon::findLinesFromSinogram info - blob ",nBlob," at ",ii+r,jj+dt_pad,fx,fy,fwhm ,f0

            end if

         end do
      end do

      !---    reconstruct the lines
      if (nBlob == 0) return                   !   no blobs = no lines... OK... (?)

      allocate (line(6, nBlob))

      nLine = 0
      do ii = 1, nBlob
         ww = bb(1, ii)                       !   distance of blob from centre of rotated image
         theta = bb(2, ii)*dtheta

         qx = Nx*0.5d0 + ww*sin(-theta)      !   remember I have found this blob at [0,y] by rotating + theta
         qy = Ny*0.5d0 + ww*cos(-theta)      !   so the real point on the line must be found by rotating -theta.

         !print *,"Lib_Radon::findLinesFromSinogram info - blob ",ii," at ",ww,theta,qx,qy

         !---    convert to y = mx + c       or   x = my + c
         ymxc = (2*theta <= PI) .or. (2*theta > 3*PI)

         if (ymxc) then       !  y = m x + c
            mm = tan(theta)
            cc = qy - mm*qx

            call ymxcInterceptBox(mm, cc, Nx, Ny, p1, p2, ok)
            if (ok) then
               nLine = nLine + 1
               line(1:2, nLine) = p1(1:2)
               line(3:4, nLine) = p2(1:2) - p1(1:2)
               line(5, nLine) = bb(3, ii)                  !   line width
               line(6, nLine) = bb(4, ii)                  !   line intensity
            end if
         else
            mm = tan(PI/2 - theta)
            cc = qx - mm*qy
            call xmycInterceptBox(mm, cc, Nx, Ny, p1, p2, ok)
            if (ok) then
               nLine = nLine + 1
               line(1:2, nLine) = p1(1:2)
               line(3:4, nLine) = p2(1:2) - p1(1:2)
               line(5, nLine) = bb(3, ii)                  !   line width
               line(6, nLine) = bb(4, ii)                  !   line intensity
            end if
         end if

         !if (ok) print *,"Lib_Radon::findLinesFromSinogram info - line ",nLine,line(:,nLine)

      end do
      ! print *,"Lib_Radon::findLinesFromSinogram info - lines ",nLine

      line(:, nLine + 1:) = 0

      deallocate (bb)

      return
   end subroutine findLinesFromSinogram2

   subroutine reconstructImage0(line, img)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given a set of lines line(1:6) = (/x1,y1,x2-x1,y2-y1,r,i/)
      !*      add to the image
      !*      The line positions are given in real units 0:Nx but the img is in integer pixels.
      !*      We assume that the img is set at positions 0.5,1.5,... Nx-0.5
      !           0     1     2     3
      !       0    _________________
      !           |     |     |     |
      !           |  x  |  x  |  x  |
      !       1   |_____|_____|_____|        pixels are squares, with the value defined at point x in the centre
      !           |     |     |     |
      !           |  x  |  x  |  x  |
      !       2   |_____|_____|_____|
      !           |     |     |     |
      !           |  x  |  x  |  x  |
      !       3   |_____|_____|_____|

      real(kind=real64), dimension(:, :), intent(in)         ::      line    !   array of lines
      real(kind=real64), dimension(0:, 0:), intent(inout)    ::      img     !   array of image

      integer             ::      Nx, Ny, nLine                         !   x & y dimensions of image/pixels, number of lines
      integer             ::      ix, iy                               !   Image x & y indices
      real(kind=real64)   ::      dd, ee, lambda                        !   see reconstructImage1, ee = reconstructed pixel intensity
      real(kind=real64), dimension(:), allocatable   ::      d1, d2   !   Prefactors to avoid duplicate calculations
      integer             ::      ii                                  !   line index

      nLine = size(line, dim=2)
      if (nLine == 0) return

      !---    size of image to draw
      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      allocate (d1(nLine))
      allocate (d2(nLine))

      do ii = 1, nLine
         d1(ii) = (line(3, ii)*line(3, ii) + line(4, ii)*line(4, ii))      !   square length of line
         d2(ii) = d1(ii)*2*line(5, ii)*line(5, ii)
         if (d1(ii) > 1.0d-16) then
            d1(ii) = 1/d1(ii)
            d2(ii) = 1/d2(ii)
         end if
      end do

      !---    construct the image
      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            ee = img(ix, iy)
            do ii = 1, nLine

               dd = line(3, ii)*(line(2, ii) - (iy + 0.5d0)) - line(4, ii)*(line(1, ii) - (ix + 0.5d0))            !   note pixel positions are at (1/2,1/2) points
               dd = dd*dd*d2(ii)
               if (dd > 4.5d0) cycle           !   only drawing out to 3 sigma range

               lambda = line(3, ii)*(ix - line(1, ii)) + line(4, ii)*(iy - line(2, ii))
               lambda = lambda*d1(ii)
               if ((line(5, ii) + lambda)*(1 + line(5, ii) - lambda) < 0) cycle

               ee = ee + line(6, ii)*exp(-dd)
            end do
            img(ix, iy) = ee
         end do
      end do

      return
   end subroutine reconstructImage0

   subroutine reconstructImage1(line, img)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given a single line(1:6) = (/x1,y1,x2-x1,y2-y1,r,i/)
      !*      add to the image
      !*      The line positions are given in real units 0:Nx but the img is in integer pixels.
      !*      We assume that the img is set at positions 0.5,1.5,... Nx-0.5
      !           0     1     2     3
      !       0    _________________
      !           |     |     |     |
      !           |  x  |  x  |  x  |
      !       1   |_____|_____|_____|        pixels are squares, with the value defined at point x in the centre
      !           |     |     |     |
      !           |  x  |  x  |  x  |
      !       2   |_____|_____|_____|
      !           |     |     |     |
      !           |  x  |  x  |  x  |
      !       3   |_____|_____|_____|

      real(kind=real64), dimension(6), intent(in)         ::      line    !   array of lines
      real(kind=real64), dimension(0:, 0:), intent(inout)    ::      img     !   array of image

      integer             ::      Nx, Ny, nLine                         !   x & y dimensions of image/pixels, number of lines
      integer             ::      ix, iy                               !   Image x & y indices
      real(kind=real64)   ::      dd, ee, lambda                        !   see reconstructImage1, ee = reconstructed pixel intensity
      real(kind=real64)   ::      d1, d2   !   Prefactors to avoid duplicate calculations
      integer             ::      ii                                  !   line index

      !nLine = size(line,dim=2)
      if (nLine == 0) return

      !---    size of image to draw
      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      d1 = (line(3)*line(3) + line(4)*line(4))      !   square length of line
      d2 = d1*2*line(5)*line(5)
      if (d1 > 1.0d-16) then
         d1 = 1/d1
         d2 = 1/d2
      end if

      !---    construct the image
      do iy = 0, Ny - 1
         do ix = 0, Nx - 1
            ee = img(ix, iy)
            dd = line(3)*(line(2) - (iy + 0.5d0)) - line(4)*(line(1) - (ix + 0.5d0))            !   note pixel positions are at (1/2,1/2) points
            dd = dd*dd*d2
            if (dd > 4.5d0) cycle           !   only drawing out to 3 sigma range

            lambda = line(3)*(ix - line(1)) + line(4)*(iy - line(2))
            lambda = lambda*d1
            if ((line(5) + lambda)*(1 + line(5) - lambda) < 0) cycle

            ee = ee + line(6)*exp(-dd)
            img(ix, iy) = ee
         end do
      end do

      return
   end subroutine reconstructImage1

   pure real(kind=real64) function linint(img, x, y)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find a linear interpolation of the buffered image at a position x (0:Nx-1)
      !*      The image is stored as pixels, but here I am returning the interpolated value at a real value 0:Nx
      !*      I imagine that the pixel value is defined at points 0.5,1.5,2.5 ... Nx-0.5
      !           0     1     2     3
      !       0    _________________
      !           |     |     |     |
      !           |  x  |  x  |  x  |
      !       1   |_____|_____|_____|        pixels are squares, with the value defined at point x in the centre
      !           |     | o  .|     |        The value at o (1.3,1.3) requires interpolating pixels 00 10 01 11
      !           |  x  |  x  |  x  |        The value at . (1.7,1.3) requires interpolating pixels 10 20 11 21
      !       2   |_____|_____|_____|
      !           |     |     |     |
      !           |  x  |  x  |  x  |
      !       3   |_____|_____|_____|

      real(kind=real64), dimension(0:, 0:), intent(in)               ::      img
      real(kind=real64), intent(in)                                ::      x, y     !x&y coordinate of point to interpolate
      integer                     ::          Nx, Ny
      integer                     ::          ix, iy   !buffered images x & y pixel indices
      real(kind=real64)           ::          ss, tt

      ix = floor(x - 0.499999d0)                     !   this is the pixel to the top left of point (x,y).
      iy = floor(y - 0.499999d0)                     !   eg o returns (00) and . returns (10)

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)

      if ((x*(Nx - x) < 0) .or. (y*(Ny - y) < 0)) then
         linint = 0
         return
      end if

      ix = max(0, min(Nx - 2, ix))        !   this places the top left pixel firmly inside the box, but
      iy = max(0, min(Ny - 2, iy))        !   means that an linear extrapolation is done for the last half of the last pixel

      ss = (x - 0.5d0) - ix             !   fraction of way across pixel. eg o returns (0.8,0.8) and . returns (0.2,0.8)
      tt = (y - 0.5d0) - iy

      linint = img(ix, iy)*(1 - ss)*(1 - tt) + img(ix + 1, iy)*ss*(1 - tt) &
               + img(ix, iy + 1)*(1 - ss)*tt + img(ix + 1, iy + 1)*ss*tt

      return
   end function linint

   pure function interceptxaxis(p1y, p2y) result(lambda)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given the line L: p = p1 + (p2-p1) lambda
      !*      at what point lambda does p(2) = 0 ie the line cross the x-axis?
      !*      note: don't need the vector direction, only the y-values as input
      !*      note: if the line is parallel to the x-axis, the intercept is +/- infinity.
      !*      I don't need to discriminate between +/- so just return + infinity.
      real(kind=real64), intent(in)        ::      p1y, p2y !y coord of the two points
      real(kind=real64)                   ::      lambda  !see above

      if (abs(p1y - p2y) < 1.0d-16) then !parallel to x axis
         lambda = huge(1.0)
      else
         lambda = p1y/(p1y - p2y)
      end if

      return
   end function interceptxaxis

   pure subroutine boxIntercepts(a, b, y, cost, sint, ok, q1, q2)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given a box defined from (-a:a,-b:b)
      !*      rotate the box by -theta  note this is passed as cos(+theta) and sin(+theta)
      !*      then find which of the sides intercept the y-axis at y
      !*                                                      p3
      !*                                                        _
      !*          p4_______p3                                  / \_
      !*           |       |                             p4   /    \_
      !*           |       |              /   ^               \------\ p2
      !*           |       |             /  -theta              \_   /
      !*           |_______|            /____ |                   \_/
      !*          p1       p2                                       p1
      !*
      !*      then back-transform to find the line across.
      !*            _______
      !*        q1 |\      |
      !*           |  \    |
      !*           |    \  |
      !*           |______\|q2
      !*
      !*      this is equivalent to finding the line across at angle +theta
      !*           _________
      !*           \        |
      !*             \    theta
      !*               \    |
      !*                 \  v
      !*
      !*      return ok if the line crosses
      real(kind=real64), intent(in)            ::      a, b                 !box extents
      real(kind=real64), intent(in)            ::      y                   !y coordinate to query
      real(kind=real64), intent(in)            ::      cost, sint           !cos(theta), sin(theta)

      real(kind=real64), dimension(2), intent(out)           ::      q1, q2  !two points of interest
      logical, intent(out)                     ::      ok                  !does the box intercept the y axis

      real(kind=real64)       ::      ppy, qqy, lambda  !rotated cordinates, see "interceptxaxis" for lambda

      logical                 ::      intercept       !does a given side of the box intercept the y axis

      q1 = 0.0d0
      q2 = 0.0d0
      ok = .false.

      !---    test p1 (-a,-b) to p2 (a,-b)
      ppy = sint*a - cost*b - y
      qqy = -sint*a - cost*b - y
      lambda = interceptxaxis(ppy, qqy)
      intercept = lambda*(1 - lambda) >= 0          !   true if 0<=lambda<=1
      if (intercept) then
         q1 = (/a*(2*lambda - 1), -b/); ok = .true.
      end if

      !---    test p2 (a,-b) to p3 (a,b)
      ppy = qqy
      qqy = -sint*a + cost*b - y
      lambda = interceptxaxis(ppy, qqy)
      intercept = lambda*(1 - lambda) >= 0          !   true if 0<=lambda<=1
      if (intercept) then
         if (ok) then
            q2 = (/a, b*(2*lambda - 1)/)
         else
            q1 = (/a, b*(2*lambda - 1)/); ok = .true.
         end if
      end if

      !---    test p3 (a,b) to p4 (-a,b)
      ppy = qqy
      qqy = sint*a + cost*b - y
      lambda = interceptxaxis(ppy, qqy)
      intercept = lambda*(1 - lambda) >= 0          !   true if 0<=lambda<=1
      if (intercept) then
         if (ok) then
            q2 = (/a*(1 - 2*lambda), b/)
         else
            q1 = (/a*(1 - 2*lambda), b/); ok = .true.
         end if
      end if

      !---    test p4 (-a,b) to p1 (-a,-b)
      ppy = qqy
      qqy = sint*a - cost*b - y
      lambda = interceptxaxis(ppy, qqy)
      intercept = lambda*(1 - lambda) >= 0          !   true if 0<=lambda<=1
      if (intercept) then
         q2 = (/-a, b*(1 - 2*lambda)/)
      end if

      return
   end subroutine boxIntercepts

   subroutine ymxcInterceptBox(m, c, Nx, Ny, p1, p2, ok)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !       given the line y = m x + c
      !       and the (real space) box 0:Nx:0:Ny
      !       find the points p1,p2 where the line intercepts the box
      !       or return ok = .false. if the line does not
      real(kind=real64), intent(in)        ::      m, c
      integer, intent(in)                  ::      Nx, Ny
      real(kind=real64), dimension(2), intent(out)  ::      p1, p2
      logical, intent(out)                 ::      ok

      real(kind=real64)       ::      qx, qy !vars for storing test results to prevent duplicate calculations.

      ok = .false.

      !---    test 1: side (0,0):(Nx,0)
      qx = -c/m
      if (qx*(Nx - qx) >= 0) then
         p1(1:2) = (/qx, 0.0d0/)
         ok = .true.
      end if
      !---    test 2: side (0,0):(0,Ny)
      qy = c
      if (qy*(Ny - qy) >= 0) then
         if (ok) then
            p2(1:2) = (/0.0d0, qy/)
            return
         else
            p1(1:2) = (/0.0d0, qy/)
            ok = .true.
         end if
      end if
      !---    test 3: side (Nx,0):(Nx,Ny)
      qy = m*Nx + c
      if (qy*(Ny - qy) >= 0) then
         if (ok) then
            p2(1:2) = (/Nx*0.999999999d0, qy/)
            return
         else
            p1(1:2) = (/Nx*0.999999999d0, qy/)
            ok = .true.
         end if
      end if
      !---    test 4: side (0,Ny):(Nx,Ny)
      qx = (Ny - c)/m
      if (qx*(Nx - qx) >= 0) then
         if (ok) then
            p2(1:2) = (/qx, Ny*0.999999999d0/)
            return
         else
            ok = .true.
            p1(1:2) = (/qx, Ny*0.999999999d0/)
         end if
      end if

      ok = .false.

      return
   end subroutine ymxcInterceptBox

   subroutine xmycInterceptBox(m, c, Nx, Ny, p1, p2, ok)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !       given the line x = m y + c
      !       and the box 0:Nx, 0:Ny
      !       find the points p1,p2 where the line intercepts the box
      !       or return ok = .false. if the line does not
      real(kind=real64), intent(in)        ::      m, c
      integer, intent(in)                  ::      Nx, Ny
      real(kind=real64), dimension(2), intent(out)  ::      p1, p2
      logical, intent(out)                 ::      ok

      real(kind=real64)       ::      qx, qy !vars for storing test results to prevent duplicate calculations.

      ok = .false.

      !---    test 1: side (0,0):(Nx,0)
      qx = c
      if (qx*(Nx - qx) >= 0) then
         p1(1:2) = (/qx, 1.0d0/)
         ok = .true.
      end if
      !---    test 2: side (0,0):(0,Ny)
      qy = -c/m
      if (qy*(Ny - qy) >= 0) then
         if (ok) then
            p2(1:2) = (/0.0d0, qy/)
            return
         else
            p1(1:2) = (/0.0d0, qy/)
            ok = .true.
         end if
      end if
      !---    test 3: side (Nx,0):(Nx,Ny)
      qy = (Nx - c)/m
      if (qy*(Ny - qy) >= 0) then
         if (ok) then
            p2(1:2) = (/Nx*0.999999999d0, qy/)
            return
         else
            p1(1:2) = (/Nx*0.999999999d0, qy/)
            ok = .true.
         end if
      end if
      !---    test 4: side (0,Ny):(Nx,Ny)
      qx = m*Ny + c
      if (qx*(Nx - qx) >= 0) then
         if (ok) then
            p2(1:2) = (/qx, 0.999999999d0*Ny/)
            return
         else
            ok = .true.
            p1(1:2) = (/qx, 0.999999999d0*Ny/)
         end if
      end if

      ok = .false.
      return
   end subroutine xmycInterceptBox

   subroutine twentyfivePointStencil(f, isMax, x, y, fwhm, fx, dbg_in)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      for the 25 points given, determine if there is a maximum near (0,0)
      !*      and if so, what is its high point and full width half max
      real(kind=real64), dimension(-2:2, -2:2), intent(in)       ::      f           !array of 9 pixels in.
      logical, intent(out)                                     ::      isMax       !does this blob have a maxima
      real(kind=real64), intent(out)                           ::      x, y, fx, fwhm ! x & y postion, intensity and FWHM of maxima

      logical, intent(in), optional ::  dbg_in !Boolean to print dd and f(0,0)
      real(kind=real64)           ::      f0, dx, dy, dxx, dxy, dyy
      real(kind=real64)           ::      tron2, dd, l1, l2, disc
      !reals
      !f0: value of quadratic fit at central pixel.
      !dx, dy: derivatives
      !dxx, dxy, dyy: partial derivatives
      !tron2: (trace of f)/2
      !dd: determinatnt
      !l1, l2: eigenvalues
      !disc: discriminant

      logical        ::      dbg     !Local version of dbg_in

      isMax = .false.
      x = 0
      y = 0
      fx = 0
      fwhm = 0

      dbg = .false.
      if (present(dbg_in)) dbg = dbg_in

!
      if (dbg) then
         print *, "twentyfivePointStencil"
         write (*, fmt='(5f8.4)') f(-2, :)
         write (*, fmt='(5f8.4)') f(-1, :)
         write (*, fmt='(5f8.4)') f(0, :)
         write (*, fmt='(5f8.4)') f(1, :)
         write (*, fmt='(5f8.4)') f(2, :)

                dd = maxval( (/maxval( f(:,-2) ),maxval( f(:,-1) ),maxval( (/f(-2,0),f(-1,0),f(1,0),f(2,0)/) ),maxval( f(:,1) ),maxval( f(:,2) )/) )  !   max excluding 0,0
         print *, "f(0,0),max ", f(0, 0), dd
      end if

          dxx = (-2*f(0,0) - 2*f(0,-1) - 2*f(0,-2) - 2*f(0,1) - 2*f(0,2) - f(-1,0) + 2*f(-2,0) - f(-1,-1) - f(-1,-2) + 2*f(-2,-1) + 2*f(-2,-2) -       &
             f(-1,1) - f(-1,2) + 2*f(-2,1) + 2*f(-2,2) - f(1,0) + 2*f(2,0) - f(1,-1) - f(1,-2) + 2*f(2,-1) + 2*f(2,-2) - f(1,1) - f(1,2) + 2*f(2,1) +  &
             2*f(2, 2))/35

      dxy = (f(-1, -1) + 2*f(-1, -2) + 2*f(-2, -1) + 4*f(-2, -2) - f(-1, 1) - 2*f(-1, 2) - 2*f(-2, 1) - 4*f(-2, 2) - f(1, -1) - &
             2*f(1, -2) - 2*f(2, -1) - 4*f(2, -2) + f(1, 1) + 2*f(1, 2) + 2*f(2, 1) + 4*f(2, 2))/100

          dyy = (-2*f(0,0) - f(0,-1) + 2*f(0,-2) - f(0,1) + 2*f(0,2) - 2*f(-1,0) - 2*f(-2,0) - f(-1,-1) + 2*f(-1,-2) - f(-2,-1) + 2*f(-2,-2) -         &
             f(-1,1) + 2*f(-1,2) - f(-2,1) + 2*f(-2,2) - 2*f(1,0) - 2*f(2,0) - f(1,-1) + 2*f(1,-2) - f(2,-1) + 2*f(2,-2) - f(1,1) + 2*f(1,2) -         &
             f(2, 1) + 2*f(2, 2))/35

      !---    do we have a maximum? check the 2nd derivatives. I expect 2 -ve eigenvalues for a max.
      tron2 = (dxx + dyy)/2                             !   trace/2
      dd = dxx*dyy - dxy*dxy                          !   determinant
      disc = sqrt(max(0.0d0, tron2*tron2 - dd))      !   discriminant
      l1 = tron2 + disc
      l2 = tron2 - disc

      if (dbg) print *, "dxx,dxy,dyy ", dxx, dxy, dyy, " eivals ", l1, l2

      if (max(l1, l2) < -1.0d-8) then
         !   yes, there is a maximum, but is it near (0,0)?
         if (abs(dd) > 1.0d-16) then

                  dx = (-f(-1,0) - 2*f(-2,0) - f(-1,-1) - f(-1,-2) - 2*f(-2,-1) - 2*f(-2,-2) - f(-1,1) - f(-1,2) - 2*f(-2,1) - 2*f(-2,2) + f(1,0) + 2*f(2,0) + &
                  f(1, -1) + f(1, -2) + 2*f(2, -1) + 2*f(2, -2) + f(1, 1) + f(1, 2) + 2*(f(2, 1) + f(2, 2)))/50

                  dy = (-f(0,-1) - 2*f(0,-2) + f(0,1) + 2*f(0,2) - f(-1,-1) - 2*f(-1,-2) - f(-2,-1) - 2*f(-2,-2) + f(-1,1) + 2*f(-1,2) + f(-2,1) + 2*f(-2,2) -  &
                  f(1, -1) - 2*f(1, -2) - f(2, -1) - 2*f(2, -2) + f(1, 1) + 2*f(1, 2) + f(2, 1) + 2*f(2, 2))/50

            dd = 1/dd
            x = (dy*dxy - dx*dyy)*dd
            y = (dx*dxy - dy*dxx)*dd

            if (dbg) then

               print *, "dx,dy ", dx, dy
!

               print *, "x,y ", x, y

                        f0 = (27*f(0,0) + 22*f(0,-1) + 7*f(0,-2) + 22*f(0,1) + 7*f(0,2) + 22*f(-1,0) + 7*f(-2,0) + 17*f(-1,-1) + 2*f(-1,-2) + 2*f(-2,-1) -     &
13*f(-2,-2) + 17*f(-1,1) + 2*f(-1,2) + 2*f(-2,1) - 13*f(-2,2) + 22*f(1,0) + 7*f(2,0) + 17*f(1,-1) + 2*f(1,-2) + 2*f(2,-1) -        &
                     13*f(2, -2) + 17*f(1, 1) + 2*(f(1, 2) + f(2, 1)) - 13*f(2, 2))/175

               fx = f0 + (2*dx*dy*dxy - dy*dy*dxx - dx*dx*dyy)*dd/2
               print *, "f0,fx,f[x,y]", f0, fx, f0 + x*dx + y*dy + (x*x*dxx + 2*x*y*dxy + y*y*dyy)/2
            end if

            if (max(abs(x), abs(y)) <= 0.5d0) then
               !   yes, there is a maximum, and its closer to (0,0) than any other point
               f0 = (27*f(0, 0) + 22*f(0, -1) + 7*f(0, -2) + 22*f(0, 1) + 7*f(0, 2) + 22*f(-1, 0) + &
                     7*f(-2, 0) + 17*f(-1, -1) + 2*f(-1, -2) + 2*f(-2, -1) - 13*f(-2, -2) + 17*f(-1, 1) + &
                     2*f(-1, 2) + 2*f(-2, 1) - 13*f(-2, 2) + 22*f(1, 0) + 7*f(2, 0) + 17*f(1, -1) + &
                     2*f(1, -2) + 2*f(2, -1) - 13*f(2, -2) + 17*f(1, 1) + 2*(f(1, 2) + f(2, 1)) - 13*f(2, 2)) &
                    /175

               fx = f0 + (2*dx*dy*dxy - dy*dy*dxx - dx*dx*dyy)*dd/2
               if (dbg) print *, "fx,f0 ", fx, f0

               !   now to get the full width half maximum, I note that a 2d gaussian has
               !   2nd derivatives daa = -f0/sa^2 , dbb = -f0/sb^2.
               !   so the geometric average sqrt( sa sb ) = ( daa dbb )**(-1/4)
               !   BUT - for this special case do I actually just want the y- width?
               isMax = (fx > 0)

               l1 = abs(l1); l2 = abs(l2)

               isMax = isMax .and. (max(l1, l2) < 5*min(l1, l2))                     !   eccentricity check - is there a factor 5 betwen the eigenvalues? there shouldnt be
               if (isMax) fwhm = sqrt(fx/sqrt(l1*l2))

               !if (isMax) print *," twentyfivePointStencil ",fx,l1,l2,fwhm
               !
            end if

         end if
      end if

      return
   end subroutine twentyfivePointStencil

   subroutine reconstructImageRGB(line, img, colourVals, imgRGB)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      given a set of lines line(1:6) = (/x1,y1,x2-x1,y2-y1,r,i/) and
      !*      an array that contains values you want to colour them by
      !*      add to the image
      real(kind=real64), dimension(:, :), intent(in)         ::      line    !array of lines
      real(kind=real64), dimension(:, :), intent(inout)      ::      img     !array of image
      real(kind=real64), dimension(:, 0:, 0:), intent(inout)      ::      imgRGB     !array of image
      real(kind=real64), dimension(:), intent(in)         ::      colourVals    !array of lines

      integer             ::      Nx, Ny, nLine                         !x & y dimensions of image/pixels, number of lines
      integer             ::      ix, iy, maxCol                               !Image x & y indices, largest colour value
      real(kind=real64)   ::      dd, ee, lambda                        !see reconstructImage1, ee = reconstructed pixel intensity
      real(kind=real64), dimension(size(line, dim=2))   ::      d1, d2   !Prefactors to avoid duplicate calculations
      integer             ::      ii                                  !line index
      real(kind=real64), dimension(3)          ::  rgb_back, rgb_fore, rgb_temp !RGB triplet for background and foregronud respectively.

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)
      nLine = size(line, dim=2)

      if (Nline == 0) return

      maxCol = maxval(colourVals)
      !print *,"DBG reconstructImageRGB  (Nx,Ny) is ", Nx, Ny

      if (maxCol > 0) then
         maxCol = maxCol
         !print *,"DBG reconstructImageRGB  (Nx,Ny) is ", Nx, Ny
      else
         maxCol = 1
      end if

      do iy = 1, (Ny) !for all pixels
         !print *,"DBG reconstructImageRGB  iy is ", iy
         do ix = 1, (Nx)
            !print *,"DBG reconstructImageRGB  ix is ", ix
            !print *,"DGB reconstructImageRGB 2, size of imgRGB is: ",size(imgRGB)
            imgRGB(1:3, (ix - 1), (iy - 1)) = img(ix, iy) !convert to RGB
            !print *,"DGB reconstructImageRGB 3, size of imgRGB is: ",size(imgRGB)
         end do
      end do

      do ii = 1, nLine
         d1(ii) = (line(3, ii)*line(3, ii) + line(4, ii)*line(4, ii))
         !d2(ii) =  d1(ii)*2*line(5,ii)*line(5,ii)
         d2(ii) = d1(ii)*2
         if (d1(ii) > 1.0d-16) then
            d1(ii) = 1/d1(ii)
            d2(ii) = 1/d2(ii)
         end if
      end do
      !^
      !---    construct the image
      !print *,"DBG Nx,Ny=  ", Nx, Ny
      do iy = 1, Ny

         do ix = 1, Nx
            ee = img(ix, iy)

            rgb_back = imgRGB(1:3, (ix - 1), (iy - 1))
            rgb_temp = (/0.0d0, 0.0d0, 0.0d0/)
            do ii = 1, nLine

               !print *, "debug RGB: ", ii
               dd = line(3, ii)*(line(2, ii) - iy) - line(4, ii)*(line(1, ii) - ix)
               dd = dd*dd*d2(ii)
               if (dd > 4.5d0) cycle           !   only drawing out to 3 sigma range
               !print *,"DBG past cycle 1  ", ix, iy

               lambda = line(3, ii)*(ix - line(1, ii)) + line(4, ii)*(iy - line(2, ii))
               lambda = lambda*d1(ii)
               if ((line(5, ii) + lambda)*(1 + line(5, ii) - lambda) < 0) cycle
               !print *,"DBG past cycle 2  ", ix, iy
               rgb_fore = getRGB_double(colourVals(ii)/maxCol) !get RGB values for overlayed line
               ee = ee + line(6, ii)*exp(-dd)
               !rgb_back = imgRGB(1:3,(ix-1),(iy-1))
               rgb_temp = rgb_temp + transparentColour(rgb_back, rgb_fore, ee)
               !imgRGB(1:3,(ix-1),(iy-1))=transparentColour( rgb_back,rgb_fore,ee)
               !print *,"DBG part 3  ", ix, iy

            end do

            !img(ix,iy) = ee

            imgRGB(1:3, (ix - 1), (iy - 1)) = rgb_temp

         end do

      end do
      return
   end subroutine reconstructImageRGB

end module Lib_Radon
