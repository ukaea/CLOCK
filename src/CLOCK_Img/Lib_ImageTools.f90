module Lib_ImageTools
   !---^^^^^^^^^^^^^^^^^^
   !*-----------------------------------------------------------------------------------------------------------------------------------
   !*    Lib_ImageTools from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
   !*      Generic 2D image array manipulation

   use iso_fortran_env
   implicit none
   private

   public              ::      cutBuffered_img !this will conflict with Lib_radon, need to sort before adding to cmakelists
   public              ::      getSubImageGridDims
   public              ::      SubImage_ctor
   public              ::      SubImageRGB_ctor
   public              ::      SubimageGrid_ctor
   public              ::      getSubImageGridExtents
   public              ::      getSubImageGridExtentsRGB
   public              ::      SubImageRGB_annotation_ctor
   public              ::      SubimageGridRGB_annotation_ctor
   !sub_image_CTOR

   !DERIVED_TYPES_HERE

   type, public         ::      subImage
      !---    defines a sub image with the bottom left corner indices stored for easy access
      real(kind=real64), dimension(:, :), allocatable    ::      subImage !pointer to array containnig sub image
      integer                                     ::      extentX, extentY     !index of the bottom right pixel
      !of a sub image w.r.t. the origional image.
   end type

   type, public         ::      subImageRGB
      !---    defines a sub image with the bottom left corner indices stored for easy access
      real(kind=real64), dimension(:, :, :), pointer      ::      subImage            !pointer to array containnig sub image
      integer                                         ::      extentX, extentY     !index of the bottom right pixel
      !of a sub image w.r.t. the origional image.
   end type

   type, public         ::         subImageGrid
      !---    defines a grid of sub images divided from a larger greyscale base image
      integer                                     ::      subImageX, subImageY !number of pixels in the X and Y directions
      !integer, dimension(:),pointer               ::      extentX,extentY     !index of the bottom right pixel of a sub image.
      !real(kind=real64),dimension(:,:),pointer    ::      subImage !pointer to array containnig sub image
      type(subImage), dimension(:, :), pointer       ::      subImg

   end type

   type, public         ::         subImageGridRGB
      !---    defines a grid of sub images divided from a larger RGB base image
      integer                                         ::      subImageX, subImageY !number of pixels in the X and Y directions
      type(subImageRGB), dimension(:, :), pointer        ::      subImgRGB

   end type

   !interface   DEFINE_INTERFACE_HERE
   !module procedure    DEFINE_INTERFACE_HERE0
   !end interface

contains
   !---^^^^^^^^

   function SubImage_ctor(subimageIn, extentX, extentY) result(SubImageOut)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Allocates the sub image and updates extent fields (x&y indices of the bottom right pixel of
      ! the sub image on the original image).
      real(kind=real64), dimension(:, :), intent(in)     ::      subimageIn !a subsection of a larger image
      integer, intent(in)                             ::      extentX, extentY
      type(subImage)                                  ::      SubImageOut

      integer                                         ::      dimX, dimY !size of the subimage

      dimX = size(subimageIn, dim=1)!get dimesnions of input inmage
      dimY = size(subimageIn, dim=2)

      allocate (SubImageOut%subImage(dimX, dimY))!allocate space to hold the image
      SubImageOut%subImage = subimageIn!put the subsection of the image in place

      SubImageOut%extentX = extentX! set extent fields
      SubImageOut%extentY = extentY
   end function SubImage_ctor

   function SubImageRGB_ctor(subimgeIn, extentX, extentY) result(SubImageOut)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Allocates the sub image and updates extent fields (x&y indices of the bottom right pixel of
      ! the sub image on the original image).
      real(kind=real64), dimension(:, :, :), intent(in)       ::      subimgeIn !a subsectino of a larger image
      integer, intent(in)                                 ::      extentX, extentY
      type(subImageRGB)                                   ::      SubImageOut

      integer                                             ::      dimX, dimY, dimZ !size of the subimage

      dimX = size(subimgeIn, dim=2)!get dimesnions of input inmage
      dimY = size(subimgeIn, dim=3)
      dimZ = 3 !RGB triplet

      allocate (SubImageOut%subImage(dimZ, dimX, dimY))!allocate space to hold the image
      SubImageOut%subImage = subimgeIn!put the subsection of the image in place

      SubImageOut%extentX = extentX! set extent fields
      SubImageOut%extentY = extentY
   end function SubImageRGB_ctor

   function SubImageRGB_annotation_ctor(subimgeIn, extentX, extentY) result(SubImageOut)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Allocates the sub image and updates extent fields (x&y indices of the bottom right pixel of
      ! the sub image on the original image). input image should be a greyscale image that you want
      !to annotate, output subimage is RGB with call channels set to the same value as the grayscale.
      real(kind=real64), dimension(:, :), intent(in)       ::      subimgeIn !a subsectino of a larger image
      integer, intent(in)                                 ::      extentX, extentY
      type(subImageRGB)                                   ::      SubImageOut

      integer                                             ::      dimX, dimY, dimZ !size of the subimage

      dimX = size(subimgeIn, dim=1)!get dimesnions of input inmage
      dimY = size(subimgeIn, dim=2)
      dimZ = 3 !RGB triplet

      allocate (SubImageOut%subImage(dimZ, dimX, dimY))!allocate space to hold the image
      SubImageOut%subImage(1, :, :) = subimgeIn!put the subsection of the image in place, channel 1
      SubImageOut%subImage(2, :, :) = subimgeIn! channel 2
      SubImageOut%subImage(3, :, :) = subimgeIn! channel 3

      SubImageOut%extentX = extentX! set extent fields
      SubImageOut%extentY = extentY
   end function SubImageRGB_annotation_ctor

   function SubimageGrid_ctor(noofSubimagesX, noofSubimagesY, imgIn, patchSizeIn) result(SubimageGridOut)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Divides imgIn into noofSubimagesX by noofSubimagesY subimages, all subimages are square except
      !for those on the bottom and right edges which contain any remaining pixels. Subimage_ctor is called
      !to allocate the subimage array and update size fields.
      integer, intent(in)     ::      noofSubimagesX
      integer, intent(in)     ::      noofSubimagesY
      real(kind=real64), dimension(0:, 0:), intent(in)     ::      imgIn    !starts at (0,0) as image read is a C function.
      integer, intent(in)      ::      patchSizeIn
      type(SubimageGrid)      ::      SubimageGridOut

      integer                 ::      ix, iy !sub image indices
      integer                 ::      remainX, remainY !remainder pixels in X & Y after division by patchSizeIn
      integer                                         ::      prevX, prevY !previous X&Y subimage extent indices
      integer                                         ::      x2, y2 !extents of current subimages indices

      !allocate subImages
      allocate (SubimageGridOut%subImg(noofSubimagesX, noofSubimagesY))
      !set size fields
      SubimageGridOut%subImageX = noofSubimagesX
      SubimageGridOut%subImageY = noofSubimagesY
      !get remainders
      remainX = size(imgIn, dim=1) - (noofSubimagesX*patchSizeIn)
      remainY = size(imgIn, dim=2) - (noofSubimagesY*patchSizeIn)
      !get extents
      call getSubImageGridExtents(noofSubimagesX, noofSubimagesY, patchSizeIn, remainX, remainY, SubimageGridOut)
      !get sub images
      prevX = 0 !initial minimum extents start at (0,0)
      prevY = 0
      do ix = 1, noofSubimagesX
         do iy = 1, noofSubimagesY
            x2 = SubimageGridOut%subImg(ix, iy)%extentX
            y2 = SubimageGridOut%subImg(ix, iy)%extentY
            SubimageGridOut%subImg(ix, iy) = SubImage_ctor(imgIn(prevX:x2, prevY:y2), x2, y2)
            prevY = y2 + 1 !use the previous Y max extent +1 to get the next min extent
         end do
         prevX = x2 + 1 !same for X
         prevY = 0 !set 0 for new column
      end do

   end function SubimageGrid_ctor

   function SubimageGridRGB_annotation_ctor(noofSubimagesX, noofSubimagesY, imgIn, patchSizeIn) result(SubimageGridOut)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Sets up a RGB subimage grid with all channels set to the a greyscale values in subimages from a greyscale iamge. The subimages
      !are set up with SubImageRGB_annotation_ctor
      integer, intent(in)     ::      noofSubimagesX
      integer, intent(in)     ::      noofSubimagesY
      real(kind=real64), dimension(0:, 0:), intent(in)     ::      imgIn    !starts at (0,0) as image read is a C function.
      integer, intent(in)      ::      patchSizeIn
      type(subImageGridRGB)      ::      SubimageGridOut

      integer                 ::      ix, iy !sub image indices
      integer                 ::      remainX, remainY !remainder pixels in X & Y after division by patchSizeIn
      integer                                         ::      prevX, prevY !previous X&Y subimage extent indices
      integer                                         ::      x2, y2 !extents of current subimages indices

      !allocate subImages
      allocate (SubimageGridOut%subImgRGB(noofSubimagesX, noofSubimagesY))
      !set size fields
      SubimageGridOut%subImageX = noofSubimagesX
      SubimageGridOut%subImageY = noofSubimagesY
      !get remainders
      remainX = size(imgIn, dim=1) - (noofSubimagesX*patchSizeIn)
      remainY = size(imgIn, dim=2) - (noofSubimagesY*patchSizeIn)
      !get extents
      call getSubImageGridExtentsRGB(noofSubimagesX, noofSubimagesY, patchSizeIn, remainX, remainY, SubimageGridOut)
      !get dummy RGB images
      prevX = 0 !initial minimum extents start at (0,0)
      prevY = 0
      do ix = 1, noofSubimagesX
         do iy = 1, noofSubimagesY
            x2 = SubimageGridOut%subImgRGB(ix, iy)%extentX
            y2 = SubimageGridOut%subImgRGB(ix, iy)%extentY
            !NEED DUMY 3 channel image array
            !SubimageGridOut%subImgRGB(ix,iy)=SubImageRGB_ctor(imgIn(prevX:x2,prevY:y2), x2, y2)
                    !!!testSubimageRGB2=SubImageRGB_annotation_ctor(testImg4(9:,9:),16,16)
            SubimageGridOut%subImgRGB(ix, iy) = SubImageRGB_annotation_ctor(imgIn(prevX:x2, prevY:y2), x2, y2)
            prevY = y2 + 1 !use the previous Y max extent +1 to get the next min extent
         end do
         prevX = x2 + 1 !same for X
         prevY = 0 !set 0 for new column
      end do

   end function SubimageGridRGB_annotation_ctor

   !function FUNCTION_NAME_HERE(VARS_IN) result(RESULT_OUT)
   !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   !DESCRIBE_FUNCTION
   !TYPE,intent(in)                 ::      VARS_IN !DEFINITION
   !TYPE                          ::      RESULT_OUT !DEFINITION
   !GET_RESULT
   !return
   !end function get0

   !SUBROUTINES
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

   subroutine getSubImageGridDims(imageIn, patchSizeIn, noofSubImgX, noofSubImgY)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Calculates how many NxN sub images will fit on a greyscale image where N is given by
      ! patchSizeIn. Outputs noofSubImgX and noofSubImgX for the horizontal and vertical number
      ! of sub images respectively. If the division is not exact the remaining pixels are added to the
      !sub images at the bottom and right side of the original image.
      real(kind=real64), dimension(:, :), intent(in)            ::      imageIn
      integer, intent(in)                         ::      patchSizeIn
      integer, intent(out)                        ::      noofSubImgX
      integer, intent(out)                        ::      noofSubImgY

      integer                                     ::      noofPixelsX, noofPixelsY !dimesions of input image

      noofPixelsX = size(imageIn, dim=1)!get input image dimensions
      noofPixelsY = size(imageIn, dim=2)

      noofSubImgX = max(1, floor(noofPixelsX/real(patchSizeIn)))       !Number of pathes in x
      noofSubImgY = max(1, floor(noofPixelsY/real(patchSizeIn)))       !Number of pathes in y

   end subroutine getSubImageGridDims

   subroutine getSubImageGridExtents(noofSubimagesX, noofSubimagesY, patchSizeIn, remainX, remainY, SubimageGridIn)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !Calculates the extents of each subimage on updates its extent fields. Fields can be written to
      !greyscale or RGB sub image grids. If a RGB and a greyscale grid is given simultaneously, the
      !greyscale grid will be selected by default and a warning displayed. If no grids are present then a
      !warning will be displayed. Must be called after allocationg subimageGrid%subimage(noofSubimagesX,noofSubimagesY)
      type(SubimageGrid), intent(inout)         ::      SubimageGridIn
      integer, intent(in)          ::      noofSubimagesX, noofSubimagesY
      integer, intent(in)          ::      remainX, remainY
      integer, intent(in)          ::      patchSizeIn

      integer, dimension(:, :), allocatable    ::      tempExtentHolderX, tempExtentHolderY
      integer                             ::      ix, iy           !subimage x & y indices

      allocate (tempExtentHolderX(noofSubimagesX, noofSubimagesX))!allocate and initialise to 0
      allocate (tempExtentHolderY(noofSubimagesX, noofSubimagesX))
      tempExtentHolderX = 0
      tempExtentHolderY = 0

      !get extents
      do ix = 1, noofSubimagesX - 1 !for all subimages except bottom and right side
         !find bottom right corner index (w.r.t. origional image)
         do iy = 1, noofSubimagesY - 1
            tempExtentHolderX(ix, iy) = patchSizeIn*ix - 1 !-1 to all extents as Images read in CLOCK start at index 0,0
            tempExtentHolderY(ix, iy) = patchSizeIn*iy - 1
         end do
      end do
      !for the right side
      do iy = 1, noofSubimagesY - 1
         tempExtentHolderX(ix, iy) = patchSizeIn*noofSubimagesX + remainX - 1
         tempExtentHolderY(ix, iy) = patchSizeIn*iy - 1
      end do
      !for the bottom
      do ix = 1, noofSubimagesX - 1
         tempExtentHolderX(ix, iy) = patchSizeIn*ix - 1
         tempExtentHolderY(ix, iy) = patchSizeIn*noofSubimagesY + remainY - 1
      end do
      !for bottom right corner
      tempExtentHolderX(ix, iy) = patchSizeIn*noofSubimagesX + remainX - 1
      tempExtentHolderY(ix, iy) = patchSizeIn*noofSubimagesY + remainY - 1

      do ix = 1, noofSubimagesX !for all subimages
         do iy = 1, noofSubimagesY
            SubimageGridIn%subImg(ix, iy)%extentX = tempExtentHolderX(ix, iy)
            SubimageGridIn%subImg(ix, iy)%extentY = tempExtentHolderY(ix, iy)
         end do
      end do

   end subroutine getSubImageGridExtents

   subroutine getSubImageGridExtentsRGB(noofSubimagesX, noofSubimagesY, patchSizeIn, remainX, remainY, SubimageGridRGBIn)
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Calculates the extents of each subimage on a RGB subgrid and updates its extent fields.
!Must be called after allocationg subimageGrid%subimage(noofSubimagesX,noofSubimagesY)
      type(SubimageGridRGB), intent(inout)      ::      SubimageGridRGBIn
      integer, intent(in)          ::      noofSubimagesX, noofSubimagesY
      integer, intent(in)          ::      remainX, remainY
      integer, intent(in)          ::      patchSizeIn

      integer, dimension(:, :), allocatable    ::      tempExtentHolderX, tempExtentHolderY
      integer                             ::      ix, iy           !subimage x & y indices

      allocate (tempExtentHolderX(noofSubimagesX, noofSubimagesX))!allocate and initialise to 0
      allocate (tempExtentHolderY(noofSubimagesX, noofSubimagesX))
      tempExtentHolderX = 0
      tempExtentHolderY = 0

      !get extents
      do ix = 1, noofSubimagesX - 1 !for all subimages except bottom and right side
         !find bottom right corner index (w.r.t. origional image)
         do iy = 1, noofSubimagesY - 1
            tempExtentHolderX(ix, iy) = patchSizeIn*ix - 1 !-1 to all extents as Images read in CLOCK start at index 0,0
            tempExtentHolderY(ix, iy) = patchSizeIn*iy - 1
         end do
      end do
      !for the right side
      do iy = 1, noofSubimagesY - 1
         tempExtentHolderX(ix, iy) = patchSizeIn*noofSubimagesX + remainX - 1
         tempExtentHolderY(ix, iy) = patchSizeIn*iy - 1
      end do
      !for the bottom
      do ix = 1, noofSubimagesX - 1
         tempExtentHolderX(ix, iy) = patchSizeIn*ix - 1
         tempExtentHolderY(ix, iy) = patchSizeIn*noofSubimagesY + remainY - 1
      end do
      !for bottom right corner
      tempExtentHolderX(ix, iy) = patchSizeIn*noofSubimagesX + remainX - 1
      tempExtentHolderY(ix, iy) = patchSizeIn*noofSubimagesY + remainY - 1

      do ix = 1, noofSubimagesX !for all RGB subimages
         do iy = 1, noofSubimagesY
            SubimageGridRGBIn%subImgRGB(ix, iy)%extentX = tempExtentHolderX(ix, iy)
            SubimageGridRGBIn%subImgRGB(ix, iy)%extentY = tempExtentHolderY(ix, iy)
         end do
      end do
   end subroutine getSubImageGridExtentsRGB

end module Lib_ImageTools
