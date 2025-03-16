
module Lib_ColourScale
!---^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_ColourScale from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      simple module for converting a value (0:1) into a colour scale 24-bit rgb
   use iso_fortran_env
   implicit none
   private

   !---

   public      ::      getRGB, getRGB_double
   public      ::      listAvailableColourScales
   public      ::      getColourScale
   public      ::      getColourScaleName
   public      ::      shadeColour
   public      ::      transparentColour

   !---
   integer, public, parameter        ::      NCOLOURSCALES = 10
   integer, public, parameter        ::      COLOURSCALE_UNSET = -1
   integer, public, parameter        ::      COLOURSCALE_GREYSCALE = 0
   integer, public, parameter        ::      COLOURSCALE_REDBLUE = 1
   integer, public, parameter        ::      COLOURSCALE_BLUERED = 2
   integer, public, parameter        ::      COLOURSCALE_VIRIDIS = 3
   integer, public, parameter        ::      COLOURSCALE_MAGMA = 4
   integer, public, parameter        ::      COLOURSCALE_BGR = 5
   integer, public, parameter        ::      COLOURSCALE_RAINBOW = 6
   integer, public, parameter        ::      COLOURSCALE_HOT = 7
   integer, public, parameter        ::      COLOURSCALE_FIRE = 8
   integer, public, parameter        ::      COLOURSCALE_IBM = 9

   integer, public                  ::      COLOURSCALE_DEFAULT = COLOURSCALE_IBM

   integer, private, parameter       ::      NKNOTS = 9

   integer, dimension(3, 0:NKNOTS - 1, 0:NCOLOURSCALES - 1), private, parameter   ::      COLOURSCALE_KNOTS = reshape((/ &
                            000,000,000 , 031,031,031 , 063,063,063 , 095,095,095 , 127,127,127 , 159,159,159 , 191,191,191 , 223,223,223 , 255,255,255 ,           &            !  COLOURSCALE_GREYSCALE
                            255,000,000 , 255,063,063 , 255,127,127 , 255,191,191 , 255,255,255 , 191,191,255 , 127,127,255 , 063,063,255 , 000,000,255 ,           &            !  COLOURSCALE_REDBLUE  
                            000,000,255 , 063,063,255 , 127,127,255 , 191,191,255 , 255,255,255 , 255,191,191 , 255,127,127 , 255,063,063 , 255,000,000 ,           &            !  COLOURSCALE_BLUERED  
                            068,021,084 , 062,053,112 , 057,086,140 , 044,118,139 , 031,150,139 , 073,173,112 , 115,208,085 , 184,219,061 , 253,231,037 ,           &            !  COLOURSCALE_VIRIDIS  
                            000,000,000 , 042,011,062 , 085,023,125 , 134,044,123 , 184,055,121 , 218,096,109 , 252,138,098 , 252,195,144 , 253,253,191 ,           &            !  COLOURSCALE_MAGMA    
                            000,000,127 , 000,000,255 , 015,111,239 , 031,223,223 , 000,255,000 , 255,255,000 , 255,127,000 , 255,000,000 , 127,000,000 ,           &            !  COLOURSCALE_BGR      
                            128,000,000 , 255,000,000 , 255,128,000 , 255,255,000 , 000,255,000 , 000,191,191 , 000,000,255 , 191,000,191 , 255,255,255 ,           &            !  COLOURSCALE_RAINBOW  
                            000,000,000 , 096,000,000 , 191,000,000 , 255,000,000 , 255,127,000 , 255,223,000 , 255,255,063 , 255,255,160 , 255,255,255 ,           &            !  COLOURSCALE_HOT      
                            255,000,000 , 255,063,000 , 255,127,000 , 255,191,000 , 255,255,000 , 255,255,063 , 255,255,127 , 255,255,191 , 255,255,255 ,           &            !  COLOURSCALE_FIRE     
                            100,143,255 , 110,119,248 , 120,094,240 , 170,066,184 , 220,038,127 , 237,068,063 , 254,097,000 , 255,137,000 , 255,176,000             &            !  COLOURSCALE_IBM      
                                                                                                   /), (/3, NKNOTS, NCOLOURSCALES/))

   character(len=10), dimension(0:NCOLOURSCALES - 1), public, parameter   ::      COLOURSCALE_NAME = (/ &
                                                                             "greyscale ", &
                                                                             "redblue   ", &
                                                                             "bluered   ", &
                                                                             "viridis   ", &
                                                                             "magma     ", &
                                                                             "bgr       ", &
                                                                             "rainbow   ", &
                                                                             "hot       ", &
                                                                             "fire      ", &
                                                                             "IBM       " &
                                                                             /)

   !---

   interface getRGB
      module procedure getRGB0
      module procedure getRGB1
   end interface

   interface getRGB_double
      module procedure getRGB_double0
      module procedure getRGB_double1
   end interface

   !---

contains
!---^^^^^^^^

   subroutine listAvailableColourScales()
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      integer     ::      ii
      print *, "Lib_ColourScale::listAvailableColourScales() info"
      do ii = 0, NCOLOURSCALES - 1
         write (*, fmt='(a)') """"//trim(COLOURSCALE_NAME(ii))//""""
      end do
      return
   end subroutine listAvailableColourScales

   function getColourScale(cs) result(colourscale)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      convert character string to integer
      character(len=*), intent(in)     ::      cs
      integer                         ::     colourscale
      integer         ::      ii
      colourscale = COLOURSCALE_UNSET
      do ii = 0, NCOLOURSCALES - 1
         if (trim(cs) == trim(COLOURSCALE_NAME(ii))) then
            colourscale = ii
            return
         end if
      end do

      return
   end function getColourScale

   function getColourScaleName(cs) result(colourscale)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      convert integer to character string
      integer, intent(in)    ::      cs
      character(len=10)     ::      colourscale
      colourscale = "unset"
      if (cs*(NCOLOURSCALES - 1 - cs) >= 0) colourscale = COLOURSCALE_NAME(cs)

      return
   end function getColourScaleName

   function getRGB0(colourscale, x) result(rgb)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      integer, intent(in)                      ::      colourscale
      real(kind=real64), intent(in)            ::      x
      integer                                 ::      rgb

      integer         ::      rr, gg, bb
      integer         ::      ii
      real(kind=real64)       ::      xx

      ii = floor(x*NKNOTS)
      if (ii < 0) then
         rr = COLOURSCALE_KNOTS(1, 0, colourscale)
         gg = COLOURSCALE_KNOTS(2, 0, colourscale)
         bb = COLOURSCALE_KNOTS(3, 0, colourscale)
      else if (ii >= NKNOTS - 1) then
         rr = COLOURSCALE_KNOTS(1, NKNOTS - 1, colourscale)
         gg = COLOURSCALE_KNOTS(2, NKNOTS - 1, colourscale)
         bb = COLOURSCALE_KNOTS(3, NKNOTS - 1, colourscale)
      else
         xx = x*NKNOTS - ii
         rr = int(COLOURSCALE_KNOTS(1, ii, colourscale)*(1 - xx) + COLOURSCALE_KNOTS(1, ii + 1, colourscale)*xx)
         gg = int(COLOURSCALE_KNOTS(2, ii, colourscale)*(1 - xx) + COLOURSCALE_KNOTS(2, ii + 1, colourscale)*xx)
         bb = int(COLOURSCALE_KNOTS(3, ii, colourscale)*(1 - xx) + COLOURSCALE_KNOTS(3, ii + 1, colourscale)*xx)
      end if

      rgb = ishft(rr, 16) + ishft(gg, 8) + bb
      return
   end function getRGB0

   function getRGB_double0(colourscale, x) result(rgb)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      integer, intent(in)                      ::      colourscale
      real(kind=real64), intent(in)            ::      x
      real(kind=real64), dimension(3)          ::      rgb

      integer                 ::      ii
      real(kind=real64)       ::      xx

      ii = floor(x*NKNOTS)
      if (ii < 0) then
         rgb(1:3) = COLOURSCALE_KNOTS(1:3, 0, colourscale)
      else if (ii >= NKNOTS - 1) then
         rgb(1:3) = COLOURSCALE_KNOTS(1:3, NKNOTS - 1, colourscale)
      else
         xx = x*NKNOTS - ii
         rgb(1:3) = COLOURSCALE_KNOTS(1:3, ii, colourscale)*(1 - xx) + COLOURSCALE_KNOTS(1:3, ii + 1, colourscale)*xx
      end if
      rgb(1:3) = rgb(1:3)/256.0

      return
   end function getRGB_double0

   function getRGB1(x) result(rgb)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), intent(in)            ::      x
      integer                                 ::      rgb
      rgb = getRGB(COLOURSCALE_DEFAULT, x)
      return
   end function getRGB1

   function getRGB_double1(x) result(rgb)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), intent(in)            ::      x
      real(kind=real64), dimension(3)          ::      rgb
      rgb = getRGB_double0(COLOURSCALE_DEFAULT, x)
      return
   end function getRGB_double1

   pure function shadeColour(rgb, f, fc) result(rgb_tone)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      change the tint/shade of a colour by mixing with intensity 0<f<1
      !*      so that f=0 returns black, f=fc returns rgb and f=1 returns white
      real(kind=real64), dimension(3), intent(in)       ::      rgb
      real(kind=real64), intent(in)                    ::      f, fc
      real(kind=real64), dimension(3)                  ::      rgb_tone

      real(kind=real64)   ::      xx

      if (f < 0.5) then
         rgb_tone = rgb*f/fc
      else
         xx = (f - fc)/(1 - fc)
         rgb_tone = xx + (1 - xx)*rgb
      end if

      return
   end function shadeColour

   pure function transparentColour(rgb_back, rgb_fore, x) result(rgb)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(3), intent(in)       ::      rgb_back, rgb_fore
      real(kind=real64), intent(in)                    ::      x
      real(kind=real64), dimension(3)                  ::      rgb

      rgb = rgb_back*(1 - x) + rgb_fore*x
      return
   end function transparentColour

end module Lib_ColourScale
