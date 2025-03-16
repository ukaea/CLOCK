module Lib_LinearFeatureDetect
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_LinearFeatureDetect from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      A module to find local linear features in a 2d image
!*      using the Radon transformation to find lines
!*      and AIC to construct a graph out of line segments
!*
!*      _________________________
!*      |     |     |     |     |         An image size Nx x Ny eg 1024 x 1024
!*      |  A  |  B  |  C  |  D  |         is broken up into Mx x My patches - here shown 4x3
!*      |_____|_____|_____|_____|
!*      |     |     |     |     |         Vertices are stored in each patch using a linked-cell list
!*      |  E  |  F  |  G  |  H  |         Note that vertices can _only_ be on the left side or the bottom side of a patch
!*      |_____|_____|_____|_____|         A line in F can therefore only connect vertices in F to F,B,G.
!*      |     |     |     |     |         Lines are stored as edges of a graph. Each line has a patch in which it is drawn,
!*      |  I  |  J  |  K  |  L  |         plus a "from" vertex ( in the block ) and a "to" vertex
!*      |_____|_____|_____|_____|         plus an intensity and a gaussian radius.

   use Lib_Radon
   use Lib_Png         !   for debugging
   use Lib_Quicksort
   use iso_fortran_env
   implicit none
   private

   public      ::      Line_ctor
   public      ::      LinkCell_ctor
   public      ::      LinearFeatureDetect_ctor
   public      ::      delete
   public      ::      report
   public      ::      clone

   public      ::      reconstructImage
   public      ::      reduceVertices
   public      ::      reduceLines
   public      ::      getNoofVerts
   public      ::      printVertsInfo
   public      ::      getLineDat
   public      ::      getNoofLines
   public      ::      getNoofActiveLines
   public      ::      getActiveLineDat
   public      ::      getActiveLineIDs
   public      ::      graph_ctor
   public      ::      getVertDestination
   public      ::      getLength
   public      ::      getContinousLineLengths
   public      ::      dbg_print_multi_seg_lines
   public      ::      dbg_get_line_array_for_image
   public      ::      makeLinesImage
   public      ::      makeLinesImagePatches
   public      ::      getSmallestLId
   public      ::      whatSubGraph
   public      ::      finalStitchCellsRGB
   public      ::      getLinesInCellRGB
   public      ::      outputPngLinesRGBsubgraph
   public      ::      dbg_print_line_subgraphs
   public      ::      dbg_report_RGB_patch_size
   public      ::      dbg_report_zerosum_patch
   public      ::      cellLineRecon

   integer, private, parameter       ::      LIB_LFD_CUT = 0     !   indicating a vertex has been cut from the list.
   logical, public                  ::      LIB_LFD_GRID = .false. !Draw gridlines in output image.
   logical, public                  ::      LIB_LFD_GREYBG = .false.!.true. !use optimal greyscale background in output image.

   type, public     ::      Line
      !---    defines a line between two vertices which can be drawn to reconstruct the image
      private
      integer                     ::      id, subgraph
      integer                     ::      from, to         !   from and to vertices
      real(kind=real64)           ::      f               !   intensity
      real(kind=real64)           ::      r               !   width of line
   end type

   type, public     ::      Vertex
      !---    defines a point in space connected to one or more lines
      integer                             ::      id
      integer                             ::      ix, iy               !   in which link cell is it found?
      real(kind=real64), dimension(2)      ::      p                   !   position of the vertices. Stored including offset, ie (1:Nx,1:Ny)
      integer                             ::      nConnex             !   number of lines connected
      integer                             ::      nConnex0            !   max number of lines space allocated
      integer, dimension(:), pointer        ::      connex              !   (1:nConnex) lines connected
      real(kind=real64)                   ::      r                   !   mean radius
   end type

   type, public     ::      LinkCell
      !---    defines a subsection of the image which contains a number of vertices and lines
      private
      integer                     ::      nV                                      !   number of vertices stored
      integer                     ::      nV0                                     !   max storage space available for vertices
      integer, dimension(:), pointer                        ::  id                  !   (1:nV) which vertices are linked to this cell
      integer                     ::      w, h                                     !   width, height of link cell block
      integer                     ::      jx, jy                                   !   offset of each link cell block
      real(kind=real64)           ::      bimg_bar                                !   average intensity level in block
      real(kind=real64), dimension(:, :), pointer            ::  buffered_img        !   (0:w+1,0:h+1) a section of the image, with 1 px border
      real(kind=real64), dimension(:, :), pointer            ::  recon_img           !   (1:w,1:h) a section of the reconstruction during optimisation
      real(kind=real64)                                   ::  rss                 !   residual sum of squares
      real(kind=real64), dimension(:, :, :), pointer            ::  recon_img_RGB     !   (3,1:w,1:h) a section of the reconstruction where lines are coloured by an RGB value.
   end type

   type, public     ::      LinearFeatureDetect
      private
      integer                                             ::      Nx, Ny           !   the size of the image

      integer                                             ::      Mx, My           !   number of patches
      integer, dimension(:), pointer                        ::      jx, jy           !   (0:Mx) jx(ix) is position of left side of block ix. right side is jx(ix+1)-1.

      type(LinkCell), dimension(:, :), pointer               ::      c               !   (0:Mx-1,0:My-1)
      integer                                             ::      nVertex         !   current number of vertices
      integer                                             ::      nLines          !   current number of lines
      integer                                             ::      nVertex0, nLine0 !   allocated number of vertices, lines
      type(Vertex), dimension(:), pointer                   ::      V               !   (0:nVertex0) array holding vertex data
      type(Line), dimension(:), pointer                     ::      L               !   (0:nLine0) array  holding line data
   end type

   type, public     ::      graph
      private

      integer, dimension(:), allocatable                     ::      start, end           !   start and end vertices
      real(kind=real64), dimension(:), allocatable           ::      length              !   total length of the connected line segments
      logical, dimension(:), allocatable                     ::      between, multi, printThis             !   is a given line segment part of another graph (no ends)
      !   May need bools for used (dont look at this one agian)
      !   and loop (this is part of a loop)
   end type

   interface Line_ctor
      module procedure Line_null
      module procedure Line_ctor0
   end interface

   interface Vertex_ctor
      module procedure Vertex_null
      module procedure Vertex_ctor0
   end interface

   interface LinkCell_ctor
      module procedure LinkCell_null
      module procedure LinkCell_ctor0
   end interface

   interface LinearFeatureDetect_ctor
      module procedure LinearFeatureDetect_null
      module procedure LinearFeatureDetect_ctor0
   end interface

   interface delete
      module procedure delete0
      module procedure delete1
      module procedure delete2
      module procedure delete3
   end interface

   interface report
      module procedure report0
      module procedure report1
      module procedure report2
      module procedure report2a
      module procedure report3
   end interface

   interface clone
      module procedure clone0
   end interface

   interface reconstructImage
      module procedure reconstructImage0
      module procedure reconstructImage1
      module procedure reconstructImage2
      module procedure reconstructImage3
   end interface

   interface addVertex
      module procedure addVertex0
      module procedure addVertex1
   end interface

   interface graph_ctor
      module procedure graph_ctor0
   end interface

contains
!---^^^^^^^^

   function Line_null() result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(Line)           ::      this! imput line
      this%from = LIB_LFD_CUT
      this%to = LIB_LFD_CUT
      this%f = 0
      this%r = 0
      this%id = LIB_LFD_CUT
      this%subgraph = LIB_LFD_CUT
      return
   end function Line_null

   function Line_ctor0(from, to, r, f, id, sub) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      integer, intent(in)              ::      from, to !vertex IDs
      real(kind=real64), intent(in)    ::      r, f     !line width & intensity
      integer, intent(in)              ::      id, sub      !Line ID.
      type(Line)                 ::      this         !Line object
      this%from = from
      this%to = to
      this%r = r
      this%f = f
      this%id = id
      this%subgraph = sub
      return
   end function Line_ctor0

   !---

   function Vertex_null() result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(Vertex)           ::      this !vertex object
      this%p = 0
      this%id = LIB_LFD_CUT
      this%ix = -1
      this%iy = -1
      this%nConnex = 0
      this%nConnex0 = 0
      nullify (this%connex)
      this%r = 0
      return
   end function Vertex_null

   function Vertex_ctor0(p, id, ix, iy) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      real(kind=real64), dimension(2), intent(in)   ::      p       !exact position of vertex
      integer, intent(in)                          ::      id      !Pixel ID
      integer, intent(in)                          ::      ix, iy   !x & y link cell corrdinates of vertex
      type(Vertex)                 ::      this                   !vertex object
      this%p = p
      this%id = id
      this%ix = ix
      this%iy = iy
      this%nConnex = 0
      this%nConnex0 = 6
      allocate (this%connex(1:this%nConnex0))
      this%r = 1.0d0
      return
   end function Vertex_ctor0

   !---

   function LinkCell_null() result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(LinkCell)           ::      this !link cell object
      this%nV0 = 0
      this%w = 0
      this%h = 0
      this%nV = 0
      nullify (this%id)
      nullify (this%buffered_img)
      nullify (this%recon_img)
      nullify (this%recon_img_RGB)
      this%rss = 0
      return
   end function LinkCell_null

   function LinkCell_ctor0(nV0, jx, jy, img, buffSizeIn) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      allocates memory, using the image chunk img(w,h) as input
      integer, intent(in)                              ::      nV0     !Max number of vertices stored
      integer, intent(in)                              ::      jx, jy   !offset of each link cell block
      !real(kind=real64),dimension(0:,0:),intent(in)     ::      img     !image chunk array in
      real(kind=real64), dimension(-(buffSizeIn - 1):, -(buffSizeIn - 1):), intent(in)     ::      img     !image chunk array in
      integer, intent(in)                                  ::      buffSizeIn !number of pixels to use for the buffer
      type(LinkCell)                 ::      this                     !link cell object
      integer         ::      ix, iy                                   !x & y pixel indices for img

      this = LinkCell_null()
      this%nV0 = nV0
      this%w = size(img, dim=1) - (2*buffSizeIn) !-2 for making the buffered image in the lfd ctor.
      this%h = size(img, dim=2) - (2*buffSizeIn) !-10 for  5 px buffer
      this%jx = jx
      this%jy = jy
      !print*,"DBG in the link cell ctor jx,jy,size(img),dim1,dim2",jx,jy,size(img),size(img,dim=1),size(img,dim=2)

      allocate (this%id(this%nV0))
      this%id = LIB_LFD_CUT

      allocate (this%buffered_img((1 - buffSizeIn):(this%w + buffSizeIn), (1 - buffSizeIn):(this%h + buffSizeIn)))

            this%buffered_img((1-buffSizeIn):this%w+buffSizeIn,(1-buffSizeIn):this%h+buffSizeIn) = img((1-buffSizeIn):this%w+buffSizeIn,(1-buffSizeIn):this%h+buffSizeIn)

      allocate (this%recon_img(this%w, this%h))
      allocate (this%recon_img_RGB(1:3, 0:(this%w - 1), 0:(this%h - 1)))
      this%recon_img = 0
      this%recon_img_RGB = 0
      this%rss = 0
      return
   end function LinkCell_ctor0

   !---

   function LinearFeatureDetect_null() result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(LinearFeatureDetect)           ::      this    !LinearFeatureDetect object
      this%Nx = 0
      this%Ny = 0
      this%Mx = 0
      this%My = 0
      this%nVertex = 0
      this%nVertex0 = 0
      this%nLines = 0
      this%nLine0 = 0
      nullify (this%L)
      nullify (this%jx)
      nullify (this%jy)
      nullify (this%c)
      nullify (this%V)
      return
   end function LinearFeatureDetect_null

   function LinearFeatureDetect_ctor0(img, n) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      allocates memory using img as input, with target block size n
      real(kind=real64), dimension(:, :), intent(in)      ::      img    !Image (full or chunk) array input
      integer, intent(in)                  ::      n
      type(LinearFeatureDetect)           ::      this                !Linear feature display
      integer             ::      Kx, Ky                               !Image chunk x & y dimensions
      integer             ::      ix, iy, ii                            !Image chunk x & y indices, vertex AND line index
      integer, parameter   ::      nV0 = 100                           !Max number of vertices stored
      real(kind=real64)   ::      img_bar                             !Average image intensity.
      real(kind=real64), dimension(:, :), allocatable       ::     bigBufferedImg    !origional image with one pixel
      !buffer around edges for lininterp
      integer             ::      buffer_size = 2                         !noof pixels added in each direction to make buffered image

      this%Nx = size(img, dim=1)
      this%Ny = size(img, dim=2)

      this%Mx = max(1, floor(this%Nx/real(n)))                         !Number of pathes in x
      this%My = max(1, floor(this%Ny/real(n)))                         !Number of pathes in y

      Kx = nint(this%Nx/real(this%Mx))
      Ky = nint(this%Ny/real(this%My))

      print *, "Lib_LinearFeatureDetect::LinearFeatureDetect_ctor0 info - image size       ", this%Nx, "x", this%Ny, "px"
      print *, "Lib_LinearFeatureDetect::LinearFeatureDetect_ctor0 info - block size (in)  ", n, "x", n, "px"
      print *, "Lib_LinearFeatureDetect::LinearFeatureDetect_ctor0 info - block count      ", this%Mx, "x", this%My
      print *, "Lib_LinearFeatureDetect::LinearFeatureDetect_ctor0 info - block size       ", Kx, "x", Ky, "px"

      allocate (this%jx(0:this%Mx))
      this%jx(0) = 1
      do ix = 1, this%Mx - 1
         this%jx(ix) = min(this%Nx, this%jx(ix - 1) + Kx)
      end do
      this%jx(this%Mx) = this%Nx + 1

      allocate (this%jy(0:this%My))
      this%jy(0) = 1
      do iy = 1, this%My - 1
         this%jy(iy) = min(this%Ny, this%jy(iy - 1) + Ky)
      end do
      this%jy(this%My) = this%Ny + 1

      !new big buffed image bit)
      allocate (bigBufferedImg((1 - buffer_size):this%Nx + buffer_size, (1 - buffer_size):this%Ny + buffer_size))
      bigBufferedImg = 0.0d0
      bigBufferedImg(1:this%Nx, 1:this%Ny) = img(1:this%Nx, 1:this%Ny)

      allocate (this%c(0:this%Mx - 1, 0:this%My - 1))
      do iy = 0, this%My - 1
         do ix = 0, this%Mx - 1
            this%c(ix, iy) = LinkCell_ctor0(nV0, this%jx(ix), this%jy(iy), &
                    bigBufferedImg( ((this%jx(ix))-buffer_size+1):(this%jx(ix+1)+buffer_size) , ((this%jy(iy))-buffer_size+1):(this%jy(iy+1)+buffer_size) ),buffer_size )
         end do
      end do
      img_bar = sum(img(:, :))/(this%Nx*this%Ny)
      this%c(:, :)%bimg_bar = img_bar
      print *, "average image level ", img_bar

      this%nVertex0 = this%Mx*this%My*nV0
      this%nVertex = 0

      allocate (this%V(0:this%nVertex0))
      do ii = 0, this%nVertex0
         this%V(ii) = Vertex_ctor()
      end do

      this%nLine0 = this%nVertex0
      this%nLines = 0
      allocate (this%L(0:this%nLine0))
      do ii = 0, this%nLine0
         this%L(ii) = Line_ctor()
      end do

      call initialLineDetect(this, buffer_size)
      return
   end function LinearFeatureDetect_ctor0

!-------

   subroutine delete0(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^
      !*      no dynamic memory to deallocate
      type(Line), intent(inout)    ::      this
      this = Line_null()
      return
   end subroutine delete0

   subroutine delete1(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^
      !*      deallocate dynamic memory
      type(LinkCell), intent(inout)    ::      this
      if (this%nV0 == 0) return
      deallocate (this%id)
      deallocate (this%buffered_img)
      deallocate (this%recon_img)
      deallocate (this%recon_img_RGB)
      this = LinkCell_null()
      return
   end subroutine delete1

   subroutine delete2(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^
      !*      deallocate dynamic memory
      type(LinearFeatureDetect), intent(inout)    ::      this !LinearFeatureDetect object
      integer         ::      ix, iy, ii    !link cell indecies, line AND vertex index
      if (this%Mx == 0) return
      do iy = 0, this%My - 1
         do ix = 0, this%Mx - 1
            call delete(this%c(ix, iy))
         end do
      end do
      do ii = 1, this%nLine0
         call delete(this%L(ii))
      end do
      deallocate (this%L)
      deallocate (this%c)
      do ii = 1, this%nVertex0
         call delete(this%V(ii))
      end do
      deallocate (this%V)
      deallocate (this%jx)
      deallocate (this%jy)
      this = LinearFeatureDetect_null()
      return
   end subroutine delete2

   subroutine delete3(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^
      !*      deallocate dynamic memory
      type(Vertex), intent(inout)    ::      this !vertex object
      if (this%nConnex0 == 0) return
      deallocate (this%connex)
      this = Vertex_null()
      return
   end subroutine delete3

   !---

   subroutine clone0(this, that)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      deep copy with allocate this = that
      type(Vertex), intent(inout)  ::      this    !vertex object to store copy
      type(Vertex), intent(in)     ::      that    !vertex object to be copied
      call delete(this)
      this%id = that%id
      this%ix = that%ix
      this%iy = that%iy
      this%p = that%p
      this%nConnex = that%nConnex
      this%nConnex0 = that%nConnex0
      allocate (this%connex(1:this%nConnex0))
      this%connex(1:this%nConnex0) = that%connex(1:this%nConnex0)
      this%r = that%r
      return
   end subroutine clone0

   !---

   subroutine report0(this, u, o)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(Line), intent(in)          ::      this !line object in
      integer, intent(in), optional     ::      u, o !file unit ID & numberof whitespace printed
      integer     ::      uu, oo
      uu = 6; if (present(u)) uu = u
      oo = 0; if (present(o)) oo = o
            write(unit=uu,fmt='(3(a,i8),3(a,f10.4))') repeat(" ",oo)//"Line [lid = ",this%id," from,to = ",this%from,",",this%to," , f,r = ",this%f,",",this%r,"]"
      return
   end subroutine report0

   subroutine report3(this, u, o)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(Vertex), intent(in)         ::      this    !vertex object in
      integer, intent(in), optional     ::      u, o     !file unit ID & number of whitespace printed
      integer     ::      uu, oo
      uu = 6; if (present(u)) uu = u     !internal storeage of u & o.
      oo = 0; if (present(o)) oo = o
            write(unit=uu,fmt='(a,i8,a,2f16.6,a,f16.6,a,i6,a,i6,a,100i6)') repeat(" ",oo)//"Vertex [id = ",this%id," p,r ",this%p(:),",",this%r," ix,iy ",this%ix,",",this%iy," connex ",this%connex(1:this%nConnex)
      return
   end subroutine report3

   subroutine report1(this, u, o)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(LinkCell), intent(in)          ::      this !LinkCell object in
      integer, intent(in), optional     ::      u, o     !file unit ID & number of whitespace printed
      integer     ::      uu, oo                       !internal storage of u & o.
      integer     ::      nV, ii                       !number of vertices, LinkCell index
      uu = 6; if (present(u)) uu = u
      oo = 0; if (present(o)) oo = o
      nV = 0
      do ii = 1, this%nV
         if (this%id(ii) /= LIB_LFD_CUT) nV = nV + 1
      end do

      write (unit=uu, fmt='(a,i8,3(a,f10.4))') repeat(" ", oo)//"LinkCell [vertices = ", nV, ", rss = ", this%rss, "]"
      return
   end subroutine report1

   subroutine report2(this, u, o)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(LinearFeatureDetect), intent(in)          ::      this  !LinearFeatureDetect object in
      integer, intent(in), optional     ::      u, o                 !file unit ID & number of whitespace printed
      integer     ::      uu, oo                                   !internal storage of u & o.
      integer     ::      ix, iy                                   !x & y patch indices
      real(kind=real64)   ::      rss                             !Residual sum of squares

      integer             ::      ii, nL, v1, v2                     !Line idex, line counter, vertices IDs
      real(kind=real64)   ::      lineLen                         !Line length
      real(kind=real64), dimension(2)  ::      pp                  !Displacment(v1,v2)

      uu = 6; if (present(u)) uu = u
      oo = 0; if (present(o)) oo = o
      rss = 0
      do iy = 0, this%My - 1
         do ix = 0, this%Mx - 1
            rss = rss + this%c(ix, iy)%rss
         end do
      end do

      nL = 0
      lineLen = 0
      do ii = 1, this%nLines
         if (this%L(ii)%id /= LIB_LFD_CUT) then
            nL = nL + 1
            v1 = this%L(ii)%from
            v2 = this%L(ii)%to
            pp = this%V(v2)%p - this%V(v1)%p
            lineLen = lineLen + norm2(pp)
            !print *,"report2 line",nL,norm2(pp)
            print *, "report2 line, length,vid1,vid2, xy1, xy2 ", ii, norm2(pp), v1, v2, this%V(v1)%p, this%V(v2)%p
         end if
      end do
      write (unit=uu, fmt='(2(a,i8),3(a,f16.4))') repeat(" ", oo)//"LinearFeatureDetect [vertices,lines = ", this%nVertex &
         , ",", nL, " , line length = ", lineLen, " rss = ", rss, "]"

      return
   end subroutine report2

   subroutine report2a(this, verbose, u, o)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      type(LinearFeatureDetect), intent(in)        ::      this
      logical, intent(in)                          ::      verbose
      integer, intent(in), optional     ::      u, o
      integer     ::      uu, oo
      integer     ::      ix, iy
      uu = 6; if (present(u)) uu = u
      oo = 0; if (present(o)) oo = o
      call report(this, uu, oo)
      if (verbose) then
         do iy = 0, this%My - 1
            do ix = 0, this%Mx - 1
               call report(this%c(ix, iy), uu, oo + 4)
            end do
         end do
      end if
      return
   end subroutine report2a

!-------

   pure integer function degreesOfFreedom(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      find the total number of degrees of freedom represented
      !*      in the link cell.
      !*      TBH not sure this is right. What does the connectivity add?
      type(LinkCell), intent(in)       ::      this                              !     LinkCell object in
      degreesOfFreedom = 2*count(this%id(1:this%nV) /= LIB_LFD_CUT) &         !     vertex positions
                         + 1                                                      !     background level
      return
   end function degreesOfFreedom

   pure real(kind=real64) function AICc(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      compute the AIC for cell (ix,iy)
      type(LinkCell), intent(in)       ::      this        !LinkCell object in
      integer             ::      kk
      real(kind=real64)   ::      reduced_chisquare_stat
      integer             ::      chisquare_dof           !chisquare degrees of freedom
      real(kind=real64), parameter   ::      TWOPI = 6.283185307d0

      kk = degreesOfFreedom(this)

      if (kk == 0) then
         AICc = 0
      else
         kk = kk + 1                   !   +1 for unseen noise variable
         chisquare_dof = this%w*this%h - 1
         reduced_chisquare_stat = this%rss/chisquare_dof       !   rss per degree of freedom
         AICc = 2*kk + this%w*this%h*log(TWOPI*reduced_chisquare_stat) + chisquare_dof                    !   standard AIC ...
         if (this%w*this%h > kk + 1) AICc = AICc + 2*kk*(kk + 1)/(this%w*this%h - kk - 1)      !   ... with correction for small sample size
      end if

      return
   end function AICc

   !pure real(kind=real64) function AICpatch( this,ix,iy )
   real(kind=real64) function AICpatch(this, ix, iy)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      compute the AIC for a small patch of link cells centred on (ix,iy)
      type(LinearFeatureDetect), intent(in)        ::      this    !LinearFeatureDetect object in
      integer, intent(in)                          ::      ix, iy   !patch x & y coordinates in
      integer             ::      jx, jy                           !3 by 3 Kernal of patches centered on patch(ix,iy)
      integer             ::      kk, ii, vid   !number of degrees of freedom, vertex index, vertex ID
      real(kind=real64)   ::      reduced_chisquare_stat
      integer             ::      chisquare_dof, area  !chisquare degrees of freedom, are/pixels
      real(kind=real64), parameter   ::      TWOPI = 6.283185307d0
      integer             ::      jj!jph hack

      kk = 0
      reduced_chisquare_stat = 0
      area = 0
      do jy = max(0, iy - 1), min(this%My - 1, iy + 1)
         do jx = max(0, ix - 1), min(this%Mx - 1, ix + 1)
            do ii = 1, this%c(jx, jy)%nV
               vid = this%c(jx, jy)%id(ii)

               if (vid /= LIB_LFD_CUT) then
                  kk = kk + 2                         !   2 for vertex location
                  kk = kk + this%V(vid)%nConnex       !   1/2 * 2 for line width intensity
               end if
            end do

            kk = kk + 1                                 !   +1 for background
            reduced_chisquare_stat = reduced_chisquare_stat + this%c(jx, jy)%rss
            area = area + this%c(jx, jy)%w*this%c(jx, jy)%h
         end do
      end do
      if (kk == 0) then
         AICpatch = 0
      else
         kk = kk + 1                   !   +1 for unseen noise variable
         chisquare_dof = area - 1
         reduced_chisquare_stat = reduced_chisquare_stat/chisquare_dof       !   rss per degree of freedom
         AICpatch = 2*kk + area*log(TWOPI*reduced_chisquare_stat) + chisquare_dof                    !   standard AIC ...
         if (area > kk + 1) AICpatch = AICpatch + 2*kk*(kk + 1)/(area - kk - 1)      !   ... with correction for small sample size
      end if

      return
   end function AICpatch

   !---

   pure function isConnected(this, vidi, vidj) result(lid)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      are vertices i and j connected? return the line number or LIB_LFD_CUT
      type(LinearFeatureDetect), intent(in)            ::      this !LinearFeatureDetect in
      integer, intent(in)                              ::      vidi, vidj !ith and jth vertex ID
      integer                                         ::      lid !line ID

      integer             ::      ii                              !Connextioin index

      do ii = 1, this%V(vidi)%nConnex
         lid = this%V(vidi)%connex(ii)
         if (this%L(lid)%from == vidi) then
            if (this%L(lid)%to == vidj) return
         else if (this%L(lid)%to == vidi) then
            if (this%L(lid)%from == vidj) return
         end if
      end do
      lid = LIB_LFD_CUT
      return
   end function isConnected

   subroutine reduceVertices(this, ok)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      attempt to reduce AIC by removing vertices
      type(LinearFeatureDetect), intent(inout)         ::      this    !LinearFeatureDetect object in
      logical, intent(out)                             ::      ok      !does the operation reduce AIC?
      integer                     ::      ix, iy                       !x & y patch indices
      logical                     ::      ok_tmp                      !Temporary holder for "ok"
      ok = .false.
      do iy = 0, this%My - 1
         do ix = 0, this%Mx - 1
            call reduceVerticesInLinkCell(this, ix, iy, ok_tmp)
            ok = ok .or. ok_tmp
         end do
      end do
      return
   end subroutine reduceVertices

   subroutine reduceVerticesInLinkCell(this, ix, iy, ok)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      attempt to reduce AIC by removing one of the vertices in link cell ix,iy
      type(LinearFeatureDetect), intent(inout)         ::      this    !LinearFeatureDetect in
      integer, intent(in)                              ::      ix, iy   !link cell x & y indices
      logical, intent(out)                             ::      ok      !does this operation reduce AIC

      real(kind=real64), dimension(:), allocatable      ::      ir      !array of inverse line widths
      real(kind=real64), dimension(:, :), allocatable    ::      d2      !scaled distance between vertices
      integer                     ::      ii, jj, kk, lid, nli, nlj, vidi, vidj, vidk, mm, lidik, lidjk

      !integers
      !ii, jj, kk: vertex, cell vertices, line indices (NB kk also used to store number of vertices)
      !lid, nli, nlj: line ID, number of lines connecting the ith and jth vertex.
      !vidi, vidj, vidk: IDs of the ith, jth and kth vertices
      !mm: used as line and cell index.
      !lidik, lidjk: ID of lines between vertices i&k and j&k
      type(LinkCell), pointer      ::      cc !pointer to LinkCell object.
      real(kind=real64)           ::      rbar, rss0, aic0, aic1, dd
      real(kind=real64), dimension(2)  ::      pp, oldpi, oldpj !displacement bretween vertices, storage for the origional ith and jth vertex.
      !reals
      !rbar: compute average width of vertex
      !rss0: residula sum of squares
      !aic0, aic1: Akaike information criterion (AIC) of patch before and after operation.
      !dd: var for storing scaled distance between vertices
      integer                     ::      trial           !trial counter
      integer, dimension(2)        ::      indx            !var to hold indices of smallest value in d2

      integer                             ::      nVid    !number of vertices in region of interest.
      integer, dimension(:), allocatable    ::      vid     !array to hold vertex IDs in region of interest.
      integer                             ::      jx, jy   !x & y cell indices.

      !---jph stuff---
      integer, dimension(4)        ::      boundingCells_i, boundingCells_j    !array to hold cells that bound the
      !lines of vertices i&j in format
      !(/ix_min,iy_min,ix_max,iy_max/)
      !Note: after findin boundingCells i & j, the min and max etents of both are stored in i for further use.

      cc => this%c(ix, iy)
      if (cc%nV == 0) return

      !---    find the vertices in this and the neighbouring cells
      nVid = 0
      do jy = max(0, iy - 1), min(this%My - 1, iy + 1)
         do jx = max(0, ix - 1), min(this%Mx - 1, ix + 1)
            nVid = nVid + this%c(jx, jy)%nV
         end do
      end do
      allocate (vid(nVid))
      nVid = cc%nV
      vid(1:nVid) = cc%id(1:nVid)
      do jy = max(0, iy - 1), min(this%My - 1, iy + 1)
         do jx = max(0, ix - 1), min(this%Mx - 1, ix + 1)
            if ((jx == ix) .and. (jy == iy)) cycle        !   have added ix,iy cell already
            kk = this%c(jx, jy)%nV
            vid(nVid + 1:nVid + kk) = this%c(jx, jy)%id(1:kk)
            nVid = nVid + kk
         end do
      end do

      !---    compute average width of vertex points
      allocate (ir(nVid))

      do ii = 1, nVid                         !   loop through all vertices in cells

         rbar = 0.0d0
         vidi = vid(ii)                          !   this is the id on the long list
         nli = this%V(vidi)%nConnex               !   number of lines connecting to this vertex
         do kk = 1, nli                            !   loop through all the lines connecting to the vertex
            lid = this%V(vidi)%connex(kk)       !   this is the id of the line on the long list
            rbar = rbar + this%L(lid)%r
         end do
         if (rbar > 0) then
            ir(ii) = nli/rbar
         else
            ir(ii) = huge(1.0)              !   making the (inverse) width huge means this vertex won't be near any others.
         end if

      end do

      !---    now compute scaled distance between vertices
      allocate (d2(nVid, cc%nV))
      d2 = huge(1.0)

      do jj = 1, cc%nV                !   only consider "from" link cell (ix,iy)
         vidj = vid(jj)
         do ii = jj + 1, nVid
            vidi = vid(ii)
            pp(1:2) = this%V(vidj)%p(1:2) - this%V(vidi)%p(1:2)
            dd = (pp(1)*pp(1) + pp(2)*pp(2))*ir(ii)*ir(jj)
            d2(ii, jj) = dd
         end do
      end do

      !---    try combining vertices
      ok = .false.
      do trial = 1, cc%nV

         !---    where should the combined vertex go?    weight with inverse radius
         indx = minloc(d2)

         ii = indx(1); jj = indx(2)
         vidi = vid(ii); vidj = vid(jj)
         if (d2(ii, jj) > 4) exit       !   not going to combine anything over 2 sigma.

         if ((vidj == LIB_LFD_CUT) .or. (vidi == LIB_LFD_CUT) .or. (vidi == vidj)) exit   !   don't think this can actually happen

         oldpi = this%V(vidi)%p
         oldpj = this%V(vidj)%p

         pp(1:2) = (oldpi*ir(ii) + oldpj*ir(jj))/(ir(ii) + ir(jj))
         !pp(1:2) = ( oldpi + oldpj )/2
         call whatCellIsVetex(this, oldpi, jx, jy)
         print *, "DBG reduceVerticesInLinkCell old pi in cell:", oldpi, jx, jy
         call whatCellIsVetex(this, pp, jx, jy)
         print *, "DBG reduceVerticesInLinkCell new pi in cell:", pp, jx, jy
         !---    what was aic before? what is it with shift?

         rss0 = sum(this%c(max(0, ix - 1):min(this%Mx - 1, ix + 1), max(0, iy - 1):min(this%My - 1, iy + 1))%rss)       !   debug only

         !---    if aic1 > aic0, then linking vertices didn't help. Undo
         call getCellsBoundingVertexLines(this, vidi, boundingCells_i) !get bounding cells of all lines connected
         call getCellsBoundingVertexLines(this, vidj, boundingCells_j) !to vertices i & j
         boundingCells_i(1) = min(boundingCells_i(1), boundingCells_j(1))   !min max extents of both sets
         boundingCells_i(2) = min(boundingCells_i(2), boundingCells_j(2))   !of bounding cells.
         boundingCells_i(3) = max(boundingCells_i(3), boundingCells_j(3))
         boundingCells_i(4) = max(boundingCells_i(4), boundingCells_j(4))

         call boundingCellsAIC(this, boundingCells_i, aic0) !get AIC for bounding cells before removing vertex.
         call vertexShiftReconBoundingCells(this, vidi, oldpi, pp, boundingCells_i)  !could calculate bounding cells per line
         call vertexShiftReconBoundingCells(this, vidj, oldpj, pp, boundingCells_i)  !may reduce computation
         call balanceIntensityBoundingCells(this, boundingCells_i)
         call findRssBoundingCells(this, boundingCells_i)

         call boundingCellsAIC(this, boundingCells_i, aic1) !This should count dofs, is the changed vertex marked for deletion.
         aic1 = aic1 - 4.0d0 !removed one vertex
         if (aic1 > aic0) then !if aic1 > aic0, then linking vertices didn't help. Undo
            call vertexShiftReconBoundingCells(this, vidi, pp, oldpi, boundingCells_i)
            call vertexShiftReconBoundingCells(this, vidj, pp, oldpj, boundingCells_i)
            call balanceIntensityBoundingCells(this, boundingCells_i)
            call findRssBoundingCells(this, boundingCells_i)
            exit
         end if

         !---    OK, If I get here, I know that combining vertices ii and jj on link cell ix,iy
         !       reduces the information content AIC and so is a Good Thing.
         !       so remove vertex jj, and update all the lines as appropriate.
         !   change vertex connection on the lines from vidj to vidi

         nlj = this%V(vidj)%nConnex

         do mm = 1, nlj                               !   check the connections from j. These must all now link to i instead.

            lidjk = this%V(vidj)%connex(mm)           !   this line connects j to k     !        x k
            if (this%L(lidjk)%from == vidj) then                                        !       /
               vidk = this%L(lidjk)%to                                                 !      /
            else if (this%L(lidjk)%to == vidj) then                                     !   j x     x i
               vidk = this%L(lidjk)%from                                               !
            end if

            if (vidk == vidi) then                                                      !        x i            x i
               !   there was already a line between i and j                            !       /
               !   ... which I'm now cutting                                           !      /        =>
               call cutLine(this, lidjk)                                                !   j x            j x

            end if
         end do

         nli = this%V(vidi)%nConnex
         nlj = this%V(vidj)%nConnex

         do mm = 1, nlj                               !   check the connections from j. These must all now link to i instead.

            lidjk = this%V(vidj)%connex(mm)           !   this line connects j to k     !        x k
            if (this%L(lidjk)%from == vidj) then                                        !       /
               vidk = this%L(lidjk)%to                                                 !      /
            else if (this%L(lidjk)%to == vidj) then                                     !   j x     x i
               vidk = this%L(lidjk)%from                                               !
            end if

            lidik = isConnected(this, vidi, vidk)       !    does i already connect to k?
            !
            if (lidik /= LIB_LFD_CUT) then                                          !            x k             x k
               !   need to remove line lidjk                                       !           / \      =>       \
               this%L(lidik)%r = (this%L(lidik)%r + this%L(lidjk)%r)/2           !          /   \               \
               this%L(lidik)%f = (this%L(lidik)%f + this%L(lidjk)%f)/2           !       j x     x i     j x     x i
               call cutLine(this, lidjk)

            else
               !   need to change lidjk to now connect ik instead                   !
               if (this%L(lidjk)%from == vidj) then                                 !            x k             x k
                  this%L(lidjk)%from = vidi                                        !           /        =>       \
               else if (this%L(lidjk)%to == vidj) then                              !          /                   \
                  this%L(lidjk)%to = vidi                                          !       j x     x i     j x     x i
               end if

               call addConnection(this%V(vidi), this%V(vidk), this%L(lidjk))

            end if

         end do
         call cutVertex(this, vidj)

         !   combining was a success, and vertex jj in the link cell should be marked for removal.
         cc%id(jj) = LIB_LFD_CUT
         d2(jj, :) = huge(1.0)
         d2(:, jj) = huge(1.0)

         !   update position of vertex i
         this%V(vidi)%p = pp(1:2)
         do jj = ii + 1, cc%nV
            if (cc%id(jj) /= LIB_LFD_CUT) then
               pp(1:2) = this%V(cc%id(jj))%p - this%V(vidi)%p
               dd = (pp(1)*pp(1) + pp(2)*pp(2))*ir(ii)*ir(jj)
               d2(jj, ii) = dd          !   fill in lower triangle
            end if
         end do

         call whatCellIsVetex(this, this%V(vidi)%p, jx, jy)
         print *, "DBG reduceVerticesInLinkCell cell indices, act cell indices:", this%V(vidi)%ix, this%V(vidi)%iy, jx, jy
         print *, "DBG reduceVerticesInLinkCell verts indices, cells merged:", vidj, vidj
         !---

         ok = .true.

      end do

      !---    finally tidy up link cell
      if (ok) then
         mm = 0
         do ii = 1, cc%nV
            if (cc%id(ii) == LIB_LFD_CUT) then
               mm = mm + 1
               cc%id(mm) = cc%id(ii)
            end if
         end do
         cc%nV = mm
      end if

      print *, "DBG reduceVerticesInLinkCell, got here veticies merged in cell:", ix, iy
      return
   end subroutine reduceVerticesInLinkCell

   subroutine reduceLines(this, ok)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      attempt to reduce AIC by removing lines
      type(LinearFeatureDetect), intent(inout)         ::      this !LinearFeatureDetect object

      logical, intent(out)                             ::      ok  !do we want to cut a given line?

      integer                     ::      ii, lidij        !line index, line ID
      real(kind=real64)           ::      rss0, aic0, aic1  !residual sum of squares, AIC before and after line reduction.
      integer                     ::      ix, iy           !centre of effect x & y coordinates/pixles
      type(Vertex), pointer        ::      vi, vj           !pointers to vertex i and j
      integer, dimension(4)        ::      boundingCells   !array to hold cells that bound a given line in format
      !(/ix_min,iy_min,ix_max,iy_max/)

      !---    find all lines, and check whether two have the same vertices
      ok = .false.
      do ii = 1, this%nLines

         lidij = this%L(ii)%id
         if (lidij == LIB_LFD_CUT) cycle         !   ignoring this line already

         vi => this%V(this%L(lidij)%from)
         vj => this%V(this%L(lidij)%to)

         !---    find centre of effect of the lines
         ix = nint((vi%ix + vj%ix)*0.5d0)
         iy = nint((vj%iy + vj%iy)*0.5d0)

         rss0 = sum(this%c(max(0, ix - 1):min(this%Mx - 1, ix + 1), max(0, iy - 1):min(this%My - 1, iy + 1))%rss)       !   debug only

         call getCellsBoundingLine(this, lidij, boundingCells)
         call boundingCellsAIC(this, boundingCells, aic0)

         call reconImageSingleLineAdjacentCells(this, lidij, .true., boundingCells) !DBG test
         call balanceIntensityBoundingCells(this, boundingCells)
         call findRssBoundingCells(this, boundingCells)
         call boundingCellsAIC(this, boundingCells, aic1)
         !aic1=aic1-4.0d0

         if ((aic1) > aic0) then
            call reconImageSingleLineAdjacentCells(this, lidij, .false., boundingCells) !DBG test
            call balanceIntensityBoundingCells(this, boundingCells)
            call findRssBoundingCells(this, boundingCells)
            cycle
         end if

         !---    if we get this far, it means that cutting the line was a good idea.
         call cutLine(this, lidij)

         ok = .true.

      end do

      return
   end subroutine reduceLines

   !---

   subroutine initialLineDetect(this, pad)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !
      type(LinearFeatureDetect), intent(inout)         ::      this        !LinearFeatureDetect object
      integer, intent(in)                              ::      pad         !how many padding pixels to use for the sinogram

      real(kind=real64), dimension(:, :), allocatable    ::      sinogram, butteredSinogram    !array containing sinogram
      real(kind=real64), dimension(:, :), pointer        ::      blockLine   !array pointers to line objects
      integer             ::      ix, iy           !x & y patch indices
      integer             ::      ii, nn           !line index, number of detected lines
      integer             ::      v1, v2           !vertex IDs 1 & 2
      real(kind=real64), dimension(2)  ::      pp  !vetex coordinates
      type(LinkCell), pointer          ::      cc  !pointer to LinkCell object.

      character(len=256)              ::      lineFilename = "sinogram.combi.png"  !JPH HACK for printing sinogram
      character(len=256)              ::      tempFilename! =""
      character(len=8) :: fmt !JPH HACK format descriptor
      character(len=8)                    ::      fileIter = ""
      real(kind=real64), dimension(:, :), allocatable    ::      butterflyMask  !array to hold butterfly mask weights
      !integer         ::      pad = 2

      this%nLines = 0
      this%nVertex = 0
      this%c(:, :)%nV = 0
      fmt = '(I0,I0)'
      do iy = 0, this%My - 1
         do ix = 0, this%Mx - 1

            butteredSinogram = 0.0d0
            cc => this%c(ix, iy)
            call construct_sinogram(cc%buffered_img, sinogram, pad, pad)
            cc%bimg_bar = sum(cc%buffered_img(1:cc%w, 1:cc%h))/(cc%w*cc%h)
            write (fileIter, fmt) ix, iy
            tempFilename = trim(fileIter)//'.sinogram.combi.png' !JPH HACK
            tempFilename = trim(fileIter)//'.baseimg.combi.png' !JPH HACK
            tempFilename = trim(fileIter)//'.buffimg.combi.png' !JPH HACK
            call findLinesFromSinogram(sinogram, cc%w, cc%h, nn, blockLine) !from pre ovelap analyusis

            if (nn > 0) then

               do ii = 1, nn
                  pp = (/this%jx(ix) - 1, this%jy(iy) - 1/) + blockLine(1:2, ii)
                  call addVertex(this, pp, v1)
                  pp = pp + blockLine(3:4, ii)
                  call addVertex(this, pp, v2)
                  call addLine(this, v1, v2, blockLine(5, ii), blockLine(6, ii))
                  call addConnection(this%V(v1), this%V(v2), this%L(this%nLines))
               end do

               deallocate (blockLine)
            end if
            deallocate (sinogram)

            call reconstructImage(this, ix, iy, thin=.false.)

         end do
      end do

      return
   end subroutine initialLineDetect

   !---

   subroutine addVertex0(this, p, vid)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      add a new vertex v at point p
      type(LinearFeatureDetect), intent(inout)             ::      this    !LinearFeatureDetect object

      real(kind=real64), dimension(2), intent(in)           ::      p       !Coordinates of vertex to be added
      integer, intent(out)                                 ::      vid     !vertex ID
      type(Vertex), dimension(:), pointer                   ::      V_tmp   !temporary pointer to vertex object

      integer     ::  ix, iy, ii    !patch x & y, vertex indices

      !---    ensure sufficient memory to store vertex location
      vid = this%nVertex + 1
      this%nVertex = vid

      if (this%nVertex > this%nVertex0) then
         allocate (V_tmp(0:this%nVertex0*2))
         do ii = 0, this%nVertex0*2
            V_tmp(ii) = vertex_null()
         end do
         do ii = 0, this%nVertex0
            call clone(V_tmp(ii), this%V(ii))
         end do
         deallocate (this%V)
         this%V => V_tmp
         this%nVertex0 = this%nVertex0*2

      end if

      this%V(vid) = Vertex_ctor(p, vid, 0, 0)

      !---    now find which link cell is being referenced, and add in correct place

      do ix = 0, this%Mx - 1
         if (this%jx(ix + 1) > p(1)) then
            do iy = 0, this%My - 1
               if (this%jy(iy + 1) > p(2)) then
                  !   link cell is ix,iy
                  call addVertex1(this%c(ix, iy), vid)
                  this%V(vid)%ix = ix
                  this%V(vid)%iy = iy

                  return
               end if
            end do
         end if
      end do

      !---    if we get here, there's a problem
      print *, "Lib_LinearFeatureDetect::addVertex0 error - p = ", p, " Nx,Ny ", this%Nx, this%Ny
      stop "addVertex0 error - trying to place vertex out of bounds"

      return
   end subroutine addVertex0

   subroutine addVertex1(this, vid)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      add a new vertex vid to a link cell at point p (0:w,0:h)
      type(LinkCell), intent(inout)                        ::      this !LinkCell object
      integer, intent(in)                                  ::      vid  !Vertex ID
      integer, dimension(:), pointer    ::      id_tmp                   !temporary pointer to maximum allowed number of vertices.

      !---    ensure sufficient memory on link cell to store vertex
      this%nV = this%nV + 1

      if (this%nV > this%nV0) then
         allocate (id_tmp(2*this%nV0))
         id_tmp(1:this%nV0) = this%id(1:this%nV0)
         deallocate (this%id)
         this%id => id_tmp
         this%nV0 = 2*this%nV0
      end if

      this%id(this%nV) = vid

      return
   end subroutine addVertex1

   subroutine addLine(this, v1, v2, r, f)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      add a new line
      type(LinearFeatureDetect), intent(inout)        ::      this !LinearFeatureDetect object
      integer, intent(in)                  ::      v1, v2           !Vertices 1 & 2 IDs
      real(kind=real64), intent(in)        ::      r, f             !line thickness & intensity
      type(Line), dimension(:), pointer     ::      L_tmp           !temporary pointer to number of lines
      integer                             ::      ii              !line indices

      !---    ensure sufficient memory to store line
      this%nLines = this%nLines + 1

      if (this%nLines > this%nLine0) then
         allocate (L_tmp(0:this%nLine0*2))
         do ii = 0, this%nLine0
            L_tmp(ii) = this%L(ii)
         end do
         do ii = this%nLine0 + 1, this%nLine0*2
            L_tmp(ii) = Line_ctor()
         end do
         deallocate (this%L)
         this%L => L_tmp
         this%nLine0 = this%nLine0*2
      end if

      this%L(this%nLines) = Line_ctor(v1, v2, r, f, this%nLines, 0)

      return
   end subroutine addLine

   subroutine addConnection(vi, vj, lij)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      exisiting line lij links vertices vi and vj

      type(Vertex), intent(inout)          ::      vi, vj       !vertex i & j
      type(Line), intent(inout)            ::      lij         !line connecting vertex i & j

      integer, dimension(:), pointer        ::      connex_tmp  !temporary pointer to lines connexted to a given vertex

      integer             ::      kk              !Connection index
      integer             ::      nCi, nCj, nC0     !number of connections to vertices i & j, maximum number of lines allowed
      integer             ::      lid             !line ID
      logical             ::      oki, okj         !Is given line connected to vertices i & j

      !---    find where the connection should go on i and j
      nCi = vi%nConnex + 1; oki = .false.
      do kk = 1, vi%nConnex
         lid = vi%connex(kk)
         if (lid == lij%id) then
            !   i already connects to line . I hope j does too!
            nCi = kk
            oki = .true.
            exit
         end if
      end do

      nCj = vj%nConnex + 1; okj = .false.
      do kk = 1, vj%nConnex
         lid = vj%connex(kk)
         if (lid == lij%id) then
            !   j already connects to line. I hope i does too!
            nCj = kk
            okj = .true.
            exit
         end if
      end do

      !---    ensure sufficient memory to store connection on vertex i and j
      nC0 = vi%nConnex0
      if (nCi > nC0) then
         allocate (connex_tmp(1:max(6, nC0*2)))
         connex_tmp(1:nC0) = vi%connex(1:nC0)
         if (nC0 > 0) deallocate (vi%connex)
         vi%connex => connex_tmp
         vi%nConnex0 = max(6, nC0*2)
      end if

      nC0 = vj%nConnex0
      if (nCj > nC0) then
         allocate (connex_tmp(1:max(6, nC0*2)))
         connex_tmp(1:nC0) = vj%connex(1:nC0)
         if (nC0 > 0) deallocate (vj%connex)
         vj%connex => connex_tmp
         vj%nConnex0 = max(6, nC0*2)
      end if

      !---
      if (oki .and. okj) then
         !   this line is already present.
         return
      else if (.not. (oki .or. okj)) then
         !   this line is not present. Add it

         vi%nConnex = nCi
         vi%connex(nCi) = lij%id
         vj%nConnex = nCj
         vj%connex(nCj) = lij%id
         lij%from = vi%id
         lij%to = vj%id
      else if (oki) then
         !   this line connects i and now needs to be added to j
         vj%nConnex = nCj
         vj%connex(nCj) = lij%id
         lij%from = vi%id
         lij%to = vj%id
      else if (okj) then
         !   this line connects j and now needs to be added to i
         vi%nConnex = nCi
         vi%connex(nCi) = lij%id
         lij%from = vi%id
         lij%to = vj%id
      end if

      return
   end subroutine addConnection

   !---

   subroutine cutVertex(this, vid)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      cut vertex vid. Note that this does not delete any connected lines or remove from link cell- that's done separately
      type(LinearFeatureDetect), intent(inout)             ::      this    !LinearFeatureDetect object
      integer, intent(in)                                  ::      vid     !vertex ID

      !---    remove from the list
      call delete(this%V(vid))
      this%nVertex = this%nVertex - 1

      return
   end subroutine cutVertex

   subroutine cutLine(this, lid)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      cuts all traces of a line from this and the vertices
      type(LinearFeatureDetect), intent(inout)        ::      this !LinearFeatureDetect object
      integer, intent(in)                  ::          lid         !line ID

      type(Vertex), pointer    ::  vv
      integer         ::      vid
      integer         ::      ii

      vid = this%L(lid)%from
      if (vid /= LIB_LFD_CUT) then
         !   need to remove the line from the list of connections on i
         vv => this%V(vid)
         do ii = 1, vv%nConnex
            if (vv%connex(ii) == lid) then
               !   this is the one
               vv%connex(ii) = vv%connex(vv%nConnex)
               vv%nConnex = vv%nConnex - 1
               if (vv%nConnex == 0) call cutVertex(this, vv%id)
               exit
            end if
         end do
      end if

      vid = this%L(lid)%to
      if (vid /= LIB_LFD_CUT) then
         !   need to remove the line from the list of connections on i
         vv => this%V(vid)
         do ii = 1, vv%nConnex
            if (vv%connex(ii) == lid) then
               !   this is the one
               vv%connex(ii) = vv%connex(vv%nConnex)
               vv%nConnex = vv%nConnex - 1
               if (vv%nConnex == 0) call cutVertex(this, vv%id)
               exit
            end if
         end do
      end if

      this%L(lid) = Line_null()

      return
   end subroutine cutLine

   !---

   subroutine singleLineRecon(this, lid, cut)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      add or remove single line.
      type(LinearFeatureDetect), intent(inout)         ::      this    !LinearFeatureDetect object
      integer, intent(in)                              ::      lid     !line ID
      logical, intent(in)                              ::      cut     !is this line to be cut?
      real(kind=real64), dimension(6)      ::      linedat   !array holding data for a single line
      integer         ::      ix, iy                       !patch coordinates
      type(Vertex), pointer        ::      vi, vj           !pointers to vertices i & j

      if (this%L(lid)%id == LIB_LFD_CUT) return

      vi => this%V(this%L(lid)%from)
      vj => this%V(this%L(lid)%to)

      if (min(vi%ix, vj%ix) == -1) return        !   this vertex has been cut. Odd that the line hasn't ???

      linedat(1:2) = vi%p(1:2)

      linedat(3:4) = vj%p(1:2) - linedat(1:2)

      !---    remove cell offset

      linedat(5) = this%L(lid)%r
      if (cut) then
         linedat(6) = -this%L(lid)%f         ! the -ve makes it cut
      else
         linedat(6) = +this%L(lid)%f         ! the +ve makes it add
      end if

      do iy = min(vi%iy, vj%iy), max(vi%iy, vj%iy)
         do ix = min(vi%ix, vj%ix), max(vi%ix, vj%ix)
            linedat(1:2) = linedat(1:2) - (/this%jx(ix) - 1, this%jy(iy) - 1/)
            !print*,"DBG singleLineRecon called for cell(ix,iy):",ix,iy
            call reconstructImage(linedat, this%c(ix, iy)%recon_img)
            linedat(1:2) = linedat(1:2) + (/this%jx(ix) - 1, this%jy(iy) - 1/)
         end do
      end do
      !print*,"DBG singleLineRecon, got here modififed, cut =",cut
      return
   end subroutine singleLineRecon

   subroutine vertexShiftRecon(this, vid, pcut, padd)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      remove single vertex vid from pbefore and all its lines, and replace at after.
      !*      note: does not update the positions of the vertices
      type(LinearFeatureDetect), intent(inout)         ::      this    !LinearFeatureDetect object
      integer, intent(in)                              ::      vid     !vertex ID
      real(kind=real64), dimension(2), intent(in)       ::      pcut, padd   !updated vertex coordinates

      integer         ::      ii, lid  !connection index, line ID

      !   move vertex to cut position and remove lines from image
      this%V(vid)%p(1:2) = pcut(1:2)
      do ii = 1, this%V(vid)%nConnex                       !   loop through all the lines connecting to the vertex
         lid = this%V(vid)%connex(ii)
         call singleLineRecon(this, lid, cut=.true.)
      end do

      !   move vertex to add position and add lines to image
      this%V(vid)%p(1:2) = padd(1:2)
      do ii = 1, this%V(vid)%nConnex                   !   loop through all the lines connecting to the vertex
         lid = this%V(vid)%connex(ii)
         call singleLineRecon(this, lid, cut=.false.)
      end do

      !   put the vertex back
      this%V(vid)%p(1:2) = pcut(1:2)

      return
   end subroutine vertexShiftRecon

   subroutine balanceIntensity(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      set the average intensity level
      type(LinkCell), intent(inout)            ::      this        !LinkCell object
      real(kind=real64)               ::      img_bar             !average patch intensity
      img_bar = sum(this%recon_img)/(this%w*this%h)
      this%recon_img = this%recon_img + this%bimg_bar - img_bar

      return
   end subroutine balanceIntensity

   subroutine balanceIntensityAndFindRss_patch(this, ix, iy)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      set the average intensity level across a patch local to ix,iy and find rss
      type(LinearFeatureDetect), intent(inout) ::      this    !LinearFeatureDetect object
      integer, intent(in)                      ::      ix, iy   !x & y pixel coordinates

      integer                 ::      kx, ky   !x & y Patch indices
      type(LinkCell), pointer  ::      cc  !pointer to LinkCell object

      do ky = max(0, iy - 1), min(this%My - 1, iy + 1)
         do kx = max(0, ix - 1), min(this%Mx - 1, ix + 1)
            cc => this%c(kx, ky)
            call balanceIntensity(cc)
            call findRss(cc)
         end do
      end do

      return
   end subroutine balanceIntensityAndFindRss_patch

   subroutine reconstructImage2(this, thin)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      stitch together image from scratch, optionally using thin lines
      type(LinearFeatureDetect), intent(inout)         ::      this    !LinearFeatureDetect object
      logical, intent(in)         ::      thin !use thin lines to reconstruct the image?

      integer         ::      ix, iy           !patch indices

      do iy = 0, this%My - 1
         do ix = 0, this%Mx - 1
            call reconstructImage(this, ix, iy, thin)
         end do
      end do

      return
   end subroutine reconstructImage2

   subroutine reconstructImage0(this, img)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      stitch together image from all the link cells
      type(LinearFeatureDetect), intent(inout)         ::      this    !LinearFeatureDetect object
      real(kind=real64), dimension(:, :), intent(out)    ::      img     !array holding reconstructed image

      integer         ::      ix, iy       !patch indices
      integer         ::      jx, jy       !patch pixel coordinates
      type(LinkCell), pointer  ::      cc  !poniter to LinkCell object

      !---    copy blocks of image back into original
      do iy = 0, this%My - 1
         do ix = 0, this%Mx - 1
            cc => this%c(ix, iy)
            jx = this%jx(ix)
            jy = this%jy(iy)
            img(jx:jx + cc%w - 1, jy:jy + cc%h - 1) = cc%recon_img(1:cc%w, 1:cc%h)
         end do
      end do

      if (LIB_LFD_GRID) then
         !---    draw the grid
         do ix = 0, this%Mx - 1
            do iy = 1, this%Ny
               img(this%jx(ix), iy) = 0.25
               img(this%jx(ix + 1) - 1, iy) = 0.25
            end do
         end do

         do iy = 0, this%My - 1
            do ix = 1, this%Nx
               img(ix, this%jy(iy)) = 0.25
               img(ix, this%jy(iy + 1) - 1) = 0.25
            end do
         end do
      end if

      return
   end subroutine reconstructImage0

   subroutine reconstructImage1(this, ix, iy, thin)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      redraw image from scratch in link cell ix,iy using all lines
      type(LinearFeatureDetect), intent(inout)         ::      this    !LinearFeatureDetect object
      integer, intent(in)                              ::      ix, iy   !x & y patch indices
      logical, intent(in)                              ::      thin    !reconstruct image with thin lines?

      integer         ::      ii, nn           !Line indices, line counter
      integer         ::      jx, jy           !link cell indices
      type(Vertex), pointer    ::      vi, vj   !pointers to vertices i & j

      real(kind=real64), dimension(:, :), allocatable      ::      lines ! line(1:6) = (/x1,y1,x2-x1,y2-y1,r,i/)
      type(LinkCell), pointer  ::      cc  !pointer to LinkCell object
      logical         ::      ok          !Are vertices i and j in the same LinkCell

      cc => this%c(ix, iy)
      if (LIB_LFD_GREYBG) then
         cc%recon_img = cc%bimg_bar
      else
         cc%recon_img = 0
      end if
      nn = 0
      do ii = 1, this%nLines
         vi => this%V(this%L(ii)%from)
         vj => this%V(this%L(ii)%to)
         do jy = min(vi%iy, vj%iy), max(vi%iy, vj%iy)
            do jx = min(vi%ix, vj%ix), max(vi%ix, vj%ix)
               if ((jx == ix) .and. (jy == iy)) nn = nn + 1
            end do
         end do

      end do

      if (nn == 0) return

      allocate (lines(6, nn))

      print *, ""
            print *,"Lib_LinearFeatureDetect::reconstructImage1 info - link cell ",ix,iy," [",this%jx(ix),",",this%jy(iy),":",this%jx(ix+1)-1,",",this%jy(iy+1)-1,"]"

      nn = 0
      do ii = 1, this%nLines
         vi => this%V(this%L(ii)%from)
         vj => this%V(this%L(ii)%to)
         ok = .false.
         do jy = min(vi%iy, vj%iy), max(vi%iy, vj%iy)
            do jx = min(vi%ix, vj%ix), max(vi%ix, vj%ix)
               if ((jx == ix) .and. (jy == iy)) ok = .true.
            end do
         end do
         if (ok) then

            nn = nn + 1

            lines(1:2, nn) = vi%p(1:2)

            lines(3:4, nn) = vj%p(1:2) - lines(1:2, nn)

            lines(1, nn) = lines(1, nn) - (this%jx(ix) - 1)
            lines(2, nn) = lines(2, nn) - (this%jy(iy) - 1)

            lines(5, nn) = this%L(ii)%r
            lines(6, nn) = this%L(ii)%f

            !print *,"line ", nn," in ",ix,iy
            !call report(this%L(ii))
            !call report(vi)
            !call report(vj)

         end if

      end do

      if (thin) lines(5, :) = 1.0d0
      call reconstructImage(lines, cc%recon_img)
      if (LIB_LFD_GREYBG) then
         call balanceIntensity(cc)
         call findRss(cc)
      end if
      deallocate (lines)

      return
   end subroutine reconstructImage1

   subroutine findRss(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^
      !*      compute rss for this link cell
      type(LinkCell), intent(inout)        ::      this    !LinkCell object
      integer                 ::      kx, ky               !x & y Link cell pixel indices
      real(kind=real64)       ::      dd                  !difference between reconstructed image and input image.
      this%rss = 0
      do ky = 1, this%h
         do kx = 1, this%w
            dd = this%recon_img(kx, ky) - this%buffered_img(kx, ky)
            this%rss = this%rss + dd*dd
         end do
      end do
      return
   end subroutine findRss

   subroutine printVertsInfo(this, u, LENGTHUNIT)
      !---^^^^^^^^^^^^^^^^^^^^^^^^
      !*      compute rss for this link cell
      type(LinearFeatureDetect), intent(in)        ::      this    !Linear feature detect object in
      integer, intent(in)    ::      u !file unit ID
      character(len=2), intent(in)                ::      LENGTHUNIT
      integer                 ::  i, countVerts, countLines !vertex/line index, counters for active vertices and lines respectively.
      character(len=2)                ::      placeholder = "na"     !length unit (default pixels)
      countVerts = 0

      do i = 0, this%nVertex !count how many vertices neglecting nulled ones.

         if (this%V(i)%ix .ne. -1) then
            countVerts = countVerts + 1
         end if
      end do

      write (unit=u, fmt='(I6,A,A)') countVerts, " Vertices, Base unit: ", LENGTHUNIT !vertices header
      write (unit=u, fmt='(A,A,A,A)') "vertex ID, vertex x coordinate/", LENGTHUNIT, ", vertex y coordinate/", LENGTHUNIT

      do i = 0, this%nVertex !write vertex data
         if (this%V(i)%ix .ne. -1) then
            write (unit=u, fmt='(I6,A,E24.9,A,E24.9)') i, ", ", this%V(i)%p(1), ", ", this%V(i)%p(2)
         end if
      end do

      countLines = 0

      do i = 0, this%nLines !count how many line neglecting lines containing nulled vertices.

         if ((this%V(this%L(i)%from)%ix .ne. -1) .or. (this%V(this%L(i)%to)%ix .ne. -1)) then
            countLines = countLines + 1
         end if
      end do

      write (unit=u, fmt='(I6,A)') countLines, " Lines" !line headers
            write(unit=u,fmt='(A,A,A,A,A)')"line ID, vertex 1 ID, vertex 2 ID, line type, intensity, thickness/",LENGTHUNIT,", length/",LENGTHUNIT,", Detection confidence, fit quality (AIC)"
      do i = 0, this%nLines !write line data

         if ((this%V(this%L(i)%from)%ix .ne. -1) .or. (this%V(this%L(i)%to)%ix .ne. -1)) then
            !write line data here
            !“line ID”, “vertex 1 ID”, “vertex 2 ID”, “line type”, “intensity”, “thickness/< unit type >”, “length/< unit type >”, “Detection confidence”, “fit quality (AIC)”
       write(unit=u,fmt='(I6,A,I6,A,I6,A,A,A,E24.9,A,E24.9,A,A,A,A,A,A)')i,", ",this%L(i)%from,", ",this%L(i)%to,", ",placeholder, &
               ", ", this%L(i)%f, ", ", this%L(i)%r, ", ", placeholder, ", ", placeholder, ", ", placeholder

         end if
      end do

      return
   end subroutine printVertsInfo

   integer function getNoofVerts(this, noof_verts)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      Returns number of vertices in a given LFD object
      type(LinearFeatureDetect), intent(in)        ::      this !LinearFeatureDetect object
      integer, intent(out)                  ::      noof_verts           !number of vertices in given LFD object

      noof_verts = this%nVertex
      return
   end function getNoofVerts

   pure function getLineDat(this) result(dat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return array of line data, dat(1:6) = (/x1,y1,x2-x1,y2-y1,r,i/)
      type(LinearFeatureDetect), intent(in)         ::      this
      real(kind=real64), dimension(6, this%nLines)      ::      dat
      integer             ::      ii

      do ii = 1, this%nLines
         dat(1, ii) = this%V(this%L(ii)%from)%p(1)
         dat(2, ii) = this%V(this%L(ii)%from)%p(2)

         dat(3, ii) = this%V(this%L(ii)%from)%p(1) - this%V(this%L(ii)%to)%p(1)
         dat(4, ii) = this%V(this%L(ii)%from)%p(2) - this%V(this%L(ii)%to)%p(2)

         dat(5, ii) = this%L(ii)%r
         dat(6, ii) = this%L(ii)%f

      end do
      return
   end function getLineDat

   integer function getNoofLines(this, noof_lines)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      Returns number of vertices in a given LFD object
      type(LinearFeatureDetect), intent(in)        ::      this !LinearFeatureDetect object
      integer, intent(out)                  ::      noof_lines           !number of vertices in given LFD object

      noof_lines = this%nLines
      return
   end function getNoofLines

   integer function getNoofActiveLines(this, noof_lines)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      Returns number of vertices in a given LFD object
      type(LinearFeatureDetect), intent(in)        ::      this !LinearFeatureDetect object
      integer, intent(out)                  ::      noof_lines           !number of vertices in given LFD object
      integer             ::      ii

      noof_lines = 0
      do ii = 1, this%nLines

         if (this%L(ii)%id /= LIB_LFD_CUT) then
            noof_lines = noof_lines + 1
         end if
      end do
      getNoofActiveLines = noof_lines
      return
   end function getNoofActiveLines

   pure function getActiveLineDat(this, nooflines) result(dat)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return array of line data, dat(1:6) = (/x1,y1,x2-x1,y2-y1,r,i/)
      type(LinearFeatureDetect), intent(in)         ::      this
      integer, intent(in)                             ::   nooflines
      real(kind=real64), dimension(6, nooflines)      ::      dat
      integer             ::      ii, lineCount

      lineCount = 0
      do ii = 1, this%nLine0
         if (this%L(ii)%id /= LIB_LFD_CUT) then
            lineCount = lineCount + 1
            dat(1, lineCount) = this%V(this%L(ii)%from)%p(1)
            dat(2, lineCount) = this%V(this%L(ii)%from)%p(2)

            dat(3, lineCount) = this%V(this%L(ii)%to)%p(1) - this%V(this%L(ii)%from)%p(1)
            dat(4, lineCount) = this%V(this%L(ii)%to)%p(2) - this%V(this%L(ii)%from)%p(2)

            dat(5, lineCount) = this%L(ii)%r
            dat(6, lineCount) = this%L(ii)%f
         end if

      end do
      return
   end function getActiveLineDat

   pure function getActiveLineIDs(this, nooflines) result(lids)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      return array of line IDs
      type(LinearFeatureDetect), intent(in)         ::      this
      integer, intent(in)                             ::   nooflines
      integer, dimension(nooflines)      ::      lids
      integer             ::      ii, lineCount

      lineCount = 0
      do ii = 1, this%nLine0

         if (this%L(ii)%id /= LIB_LFD_CUT) then
            lineCount = lineCount + 1
            lids(lineCount) = this%L(ii)%id
         end if
         !lineCount=lineCount+1
         !lids(lineCount) = this%L(ii)%id

      end do
      return

   end function getActiveLineIDs

   function graph_ctor0(noof_active_lines) result(this)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      integer, intent(in)                          ::      noof_active_lines

      type(graph)                 ::      this                   !graph object
      integer                     ::      ii

      allocate (this%start(1:noof_active_lines))
      allocate (this%end(1:noof_active_lines))
      allocate (this%between(1:noof_active_lines))
      allocate (this%multi(1:noof_active_lines))
      allocate (this%printThis(1:noof_active_lines))
      allocate (this%length(1:noof_active_lines))
      do ii = 1, noof_active_lines
         this%start(ii) = 0
         this%end(ii) = 0
         this%between = .false.
         this%length = 0.0d0
         this%multi = .false.
         this%printThis = .false.

      end do
      return
   end function graph_ctor0

   pure function getVertDestination(this, lineId, vertStartId) result(vertEndId)
      !takes in a vertex ID for a given line and finds the other associated vertex accesed from a
      !linear feature detect object.
      type(LinearFeatureDetect), intent(in)            ::      this
      integer, intent(in)                              ::   lineId, vertStartId
      integer                                         ::      vertEndId
      if (vertStartId == (this%L(lineId)%from)) then
         vertEndId = this%L(lineId)%to
      else
         vertEndId = this%L(lineId)%from
      end if
   end function getVertDestination

   pure function getLength(lfdIn, vertexIn1, vertexIn2) result(lengthOut)
      !Finds the distance between to vertices accessed from a linear feature detect object
      type(LinearFeatureDetect), intent(in)            ::      lfdIn
      integer, intent(in)                              ::      vertexIn1, vertexIn2
      real(kind=real64)                               ::      lengthOut
      real(kind=real64), dimension(2)                  ::      pp

      pp = lfdIn%V(vertexIn1)%p - lfdIn%V(vertexIn2)%p
      lengthOut = norm2(pp)

   end function getLength

   subroutine iterateConnectivity(vertDirectionIn, graphIn, lineIdIn, lfdIn, gIdx)
      !takes in a line known to be an end point/knot and iterates through the lines
      !connected to it's other end to get the length
      type(LinearFeatureDetect), intent(in)            ::      lfdIn
      integer, intent(in)                             ::      vertDirectionIn, lineIdIn, gIdx
      type(graph), intent(inout)                       ::      graphIn
      !logical,dimension(:),intent(inout)              ::      lineOfInterest

      integer                                         ::      nextLine, vertBegin, vertEnd
      real(kind=real64)                               ::      lengthTemp

      logical                                         ::      endLine = .false.
      !get length of ini segment from graphIn(i)
      lengthtemp = lengthtemp + getLength(lfdIn, graphIn%start(gIdx), graphIn%end(gIdx))
      !set begin vertex
      vertEnd = graphIn%end(gIdx)
      !set nextline to current for loop
      nextLine = lineIdIn

      !while there is a next line
      endLine = .false.
      do while (endLine .eqv. .false.)
         !set begin vertex=end vertex
         vertBegin = vertEnd
         !from the begin vertex find end vertex
         vertEnd = getVertDestination(lfdIn, lineIdIn, vertBegin)
         !mark this graph has more than one line segment.
         graphIn%multi(gIdx) = .true.
         graphIn%printThis(gIdx) = .true.
         !is there a next line?
         if (lfdIn%V(vertEnd)%nConnex == 2) then
            !then get the next line
            if (lfdIn%V(vertEnd)%connex(1) == nextLine) then
               nextline = lfdIn%V(vertEnd)%connex(2)
            else
               nextline = lfdIn%V(vertEnd)%connex(1)
            end if
         else
            endLine = .true.
         end if
         lengthtemp = lengthtemp + getLength(lfdIn, vertbegin, vertend)
      end do
      graphIn%length(gIdx) = lengthTemp
   end subroutine iterateConnectivity

   subroutine getContinousLineLengths(noofLines, lineIdsIn, graphIn, lfdIn)
      !iterates through all active lines finds all connected line segments and returns
      !the continous length of each set of lines (including single lines). These are
      !stored in a "graph" data structure. Currently a node with three or more lines
      !connectd to it (reffered to as knot) is considered a stopping point.
      type(LinearFeatureDetect), intent(in)            ::      lfdIn
      integer, intent(in)                             ::      noofLines !active lines
      integer, dimension(noofLines)                    ::      lineIdsIn
      type(graph), intent(inout)                       ::      graphIn

      integer                                         ::      ii

      do ii = 1, noofLines
         if ((lfdIn%V(lfdIn%L(lineIdsIn(ii))%to)%nConnex == 2) .and. (lfdIn%V(lfdIn%L(lineIdsIn(ii))%from)%nConnex == 2)) then
            graphIn%between(ii) = .true.
            cycle
         end if

         if (lfdIn%V(lfdIn%L(lineIdsIn(ii))%from)%nConnex == 1) then
            graphIn%start(ii) = lfdIn%V(lfdIn%L(lineIdsIn(ii))%from)%id
            if (lfdIn%V(lfdIn%L(lineIdsIn(ii))%to)%nConnex == 1) then

               !single line segment
               graphIn%end(ii) = lfdIn%V(lfdIn%L(lineIdsIn(ii))%to)%id
               graphIn%length(ii) = getLength(lfdIn, graphIn%start(ii), graphIn%end(ii))
            elseif (lfdIn%V(lfdIn%L(lineIdsIn(ii))%to)%nConnex == 2) then
               !there is another line connected to vert 2
               graphIn%end(ii) = lfdIn%V(lfdIn%L(lineIdsIn(ii))%to)%id
               call iterateConnectivity(graphIn%end(ii), graphIn, lineIdsIn(ii), lfdIn, ii)
            else
               !must be a knot
               graphIn%end(ii) = lfdIn%V(lfdIn%L(lineIdsIn(ii))%to)%id
               graphIn%length(ii) = getLength(lfdIn, graphIn%start(ii), graphIn%end(ii))
            end if
         elseif (lfdIn%V(lfdIn%L(lineIdsIn(ii))%from)%nConnex == 2) then
            !there is connectivty to another line, find the rest
            if (lfdIn%V(lfdIn%L(lineIdsIn(ii))%to)%nConnex == 1) then
               !vert 2 is the start point
               graphIn%start(ii) = lfdIn%V(lfdIn%L(lineIdsIn(ii))%to)%id
               graphIn%end(ii) = lfdIn%V(lfdIn%L(lineIdsIn(ii))%from)%id
               call iterateConnectivity(graphIn%end(ii), graphIn, lineIdsIn(ii), lfdIn, ii)
            else
               !vert 2 must be a knot
               graphIn%start(ii) = lfdIn%V(lfdIn%L(lineIdsIn(ii))%to)%id
               graphIn%end(ii) = lfdIn%V(lfdIn%L(lineIdsIn(ii))%from)%id
               call iterateConnectivity(graphIn%end(ii), graphIn, lineIdsIn(ii), lfdIn, ii)
            end if

         else !vert 1 must be a knot
            graphIn%start(ii) = lfdIn%V(lfdIn%L(lineIdsIn(ii))%from)%id
            if (lfdIn%V(lfdIn%L(lineIdsIn(ii))%to)%nConnex == 1) then
               !single line segment connected to knot
               graphIn%end(ii) = lfdIn%V(lfdIn%L(lineIdsIn(ii))%to)%id
               graphIn%length(ii) = getLength(lfdIn, graphIn%start(ii), graphIn%end(ii))
            elseif (lfdIn%V(lfdIn%L(lineIdsIn(ii))%to)%nConnex == 2) then
               !there is another line connected to vert 2
               graphIn%end(ii) = lfdIn%V(lfdIn%L(lineIdsIn(ii))%to)%id
               call iterateConnectivity(graphIn%end(ii), graphIn, lineIdsIn(ii), lfdIn, ii)
            else
               !must be two knots joined by single line
               graphIn%end(ii) = lfdIn%V(lfdIn%L(lineIdsIn(ii))%to)%id
               graphIn%length(ii) = getLength(lfdIn, graphIn%start(ii), graphIn%end(ii))
            end if
         end if

      end do
   end subroutine getContinousLineLengths

   subroutine dbg_print_multi_seg_lines(nooflines, graphIn)
      !loops through all graphs and prints any composed of more than one line.
      integer, intent(in)             ::      nooflines
      type(graph), intent(in)          ::      graphIn
      integer                         ::      ii
      print *, "dbg_print_multi_seg_lines starting"
      do ii = 1, nooflines
         if (graphIn%multi(ii) .eqv. .true.) then
            print *, "DEBUG, idx: ", ii, ", length is ", graphIn%length(ii)
         end if
      end do

   end subroutine dbg_print_multi_seg_lines

   subroutine dbg_get_line_array_for_image(nooflines, graphIn, lfdIn, lidsIn, Nx, Ny, imgOrig)
      integer, intent(in)             ::      nooflines
      type(graph), intent(inout)          ::      graphIn
      type(LinearFeatureDetect), intent(in)    ::  lfdIn
      integer, dimension(:), intent(in)         ::  lidsIn
      real(kind=real64), dimension(:, :), allocatable   ::      lines_to_print
      real(kind=real64), dimension(:, :), intent(in)       ::      imgOrig !origional image
      real(kind=real64), dimension(:, :), allocatable        ::      img !array to hold image
      integer                             ::      ii, count, currentLid, nconnexLine, Nx, Ny
      count = 0
      nconnexLine = 0

      do ii = 1, nooflines
         nconnexLine = 0
         nconnexLine = nconnexLine + lfdIn%V(lfdIn%L(lidsIn(ii))%to)%nConnex
         nconnexLine = nconnexLine + lfdIn%V(lfdIn%L(lidsIn(ii))%from)%nConnex
         if (nconnexLine > 2) then
            graphIn%printThis(ii) = .true.
            count = count + 1
         end if
         !end if
      end do
      allocate (lines_to_print(6, count))
      count = 0
      do ii = 1, nooflines
         if (graphIn%printThis(ii)) then
            count = count + 1
            !line(1:6) = (/x1,y1,x2-x1,y2-y1,r,i/)
            lines_to_print(1, count) = lfdIn%V(lfdIn%L(lidsIn(ii))%to)%p(1)
            lines_to_print(2, count) = lfdIn%V(lfdIn%L(lidsIn(ii))%to)%p(2)
            lines_to_print(3, count) = lfdIn%V(lfdIn%L(lidsIn(ii))%from)%p(1) - lines_to_print(count, 1)
            lines_to_print(4, count) = lfdIn%V(lfdIn%L(lidsIn(ii))%from)%p(2) - lines_to_print(count, 2)
            lines_to_print(5, count) = lfdIn%L(lidsIn(ii))%r
            lines_to_print(6, count) = lfdIn%L(lidsIn(ii))%f
                    print *,"DBG: the lines to print are:",lines_to_print(count,1),lines_to_print(count,2),lines_to_print(count,3),lines_to_print(count,4),lines_to_print(count,5),lines_to_print(count,6)
         end if
      end do

   end subroutine dbg_get_line_array_for_image

   subroutine makeLinesImage(line, img)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !copy of reconstructImage0 from lib_radon, think there is a name clash here
      !*      given a set of lines line(1:6) = (/x1,y1,x2-x1,y2-y1,r,i/)
      !*      add to the image
      real(kind=real64), dimension(:, :), intent(in)         ::      line    !array of lines
      real(kind=real64), dimension(:, :), intent(inout)      ::      img     !array of image

      integer             ::      Nx, Ny, nLine                         !x & y dimensions of image/pixels, number of lines
      integer             ::      ix, iy                               !Image x & y indices
      real(kind=real64)   ::      dd, ee, lambda                        !see reconstructImage1, ee = reconstructed pixel intensity
      real(kind=real64), dimension(size(line, dim=2))   ::      d1, d2   !Prefactors to avoid duplicate calculations
      integer             ::      ii                                  !line index

      Nx = size(img, dim=1)
      Ny = size(img, dim=2)
      nLine = size(line, dim=2)
      print *, "DBG, printing lines noof: ", nLine
      if (Nline == 0) return

      do ii = 1, nLine
         d1(ii) = (line(3, ii)*line(3, ii) + line(4, ii)*line(4, ii))
         d2(ii) = d1(ii)*2*line(5, ii)*line(5, ii)
         if (d1(ii) > 1.0d-16) then
            d1(ii) = 1/d1(ii)
            d2(ii) = 1/d2(ii)
         end if
      end do

      !---    construct the image
      do iy = 1, Ny
         do ix = 1, Nx
            ee = 0.0d0! img(ix,iy)
            do ii = 1, nLine

               dd = line(3, ii)*(line(2, ii) - iy) - line(4, ii)*(line(1, ii) - ix)
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
   end subroutine makeLinesImage

   subroutine makeLinesImagePatches(lfdIn, graphIn, lidsIn, noofLinesInterest, Nxin, Nyin)
      !reconstructs an image patch by patch selecting lines via connectivty
      integer, intent(in)             ::      noofLinesInterest, Nxin, Nyin
      type(graph), intent(inout)          ::      graphIn
      type(LinearFeatureDetect), intent(inout)    ::  lfdIn
      integer, dimension(:), intent(in)         ::  lidsIn
      real(kind=real64), dimension(:, :), allocatable     ::      img     !array of image

      integer         ::      ix, iy           !patch indices

      do iy = 0, lfdIn%My - 1
         do ix = 0, lfdIn%Mx - 1
            call reconstructImage3(lfdIn, ix, iy, graphIn, lidsIn, noofLinesInterest) !calls reconstructImage3
         end do
      end do

      allocate (img(Nxin, Nyin))
      call reconstructImage(lfdin, thin=.false.)
      call reconstructImage0(lfdIn, img)
      !write image
      call writePng("DBG_linerecon.png", img, .true.) !negative = negative

   end subroutine makeLinesImagePatches

   subroutine reconstructImage3(this, ix, iy, graphIn, lidsIn, noofLinesInterest)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      redraw image from scratch in link cell ix,iy using all lines selcetd from graph
      type(LinearFeatureDetect), intent(inout)         ::      this    !LinearFeatureDetect object
      integer, intent(in)                              ::      ix, iy   !x & y patch indices
      !logical,intent(in)                              ::      thin    !reconstruct image with thin lines?
      type(graph), intent(in)                          ::      graphIn
      integer, dimension(:), intent(in)                 ::      lidsIn
      integer, intent(in)                              ::      noofLinesInterest

      integer         ::      ii, nn           !Line indices, line counter
      integer         ::      jx, jy           !link cell indices
      type(Vertex), pointer    ::      vi, vj   !pointers to vertices i & j

      real(kind=real64), dimension(:, :), allocatable      ::      lines ! line(1:6) = (/x1,y1,x2-x1,y2-y1,r,i/)
      type(LinkCell), pointer  ::      cc  !pointer to LinkCell object
      logical         ::      ok          !Are vertices i and j in the same LinkCell

      cc => this%c(ix, iy)
      if (LIB_LFD_GREYBG) then
         cc%recon_img = cc%bimg_bar
      else
         cc%recon_img = 0
      end if
      nn = 0
      do ii = 1, noofLinesInterest
         vi => this%V(this%L(lidsIn(ii))%from)
         vj => this%V(this%L(lidsIn(ii))%to)
         do jy = min(vi%iy, vj%iy), max(vi%iy, vj%iy)
            do jx = min(vi%ix, vj%ix), max(vi%ix, vj%ix)
               if ((jx == ix) .and. (jy == iy)) then
                  if (graphIn%printThis(ii) .eqv. .true.) nn = nn + 1
               end if
            end do
         end do

      end do
      if (nn == 0) return

      allocate (lines(6, nn))

      nn = 0
      do ii = 1, noofLinesInterest
         vi => this%V(this%L(lidsIn(ii))%from)
         vj => this%V(this%L(lidsIn(ii))%to)
         ok = .false.
         do jy = min(vi%iy, vj%iy), max(vi%iy, vj%iy)
            do jx = min(vi%ix, vj%ix), max(vi%ix, vj%ix)
               if ((jx == ix) .and. (jy == iy)) then
                  if (graphIn%printThis(ii) .eqv. .true.) ok = .true.
               end if
            end do
         end do
         if (ok) then

            nn = nn + 1

            lines(1:2, nn) = vi%p(1:2)

            lines(3:4, nn) = vj%p(1:2) - lines(1:2, nn)

            lines(1, nn) = lines(1, nn) - (this%jx(ix) - 1)
            lines(2, nn) = lines(2, nn) - (this%jy(iy) - 1)

            lines(5, nn) = this%L(ii)%r
            lines(6, nn) = this%L(ii)%f

            !print *,"line ", nn," in ",ix,iy
            !call report(this%L(ii))
            !call report(vi)
            !call report(vj)

         end if

      end do

      !if (thin) lines(5,:) = 1.0d0
      call reconstructImage(lines, cc%recon_img)

      if (LIB_LFD_GREYBG) then
         call balanceIntensity(cc)
         call findRss(cc)
      end if
      deallocate (lines)

      return
   end subroutine reconstructImage3

   subroutine whatSubGraph(lfdIn, lidsIn, noofLineIn)
      !iterates throug all lines indicated by array of line IDs in (lidsIn)
      type(LinearFeatureDetect), intent(inout)         ::      lfdIn    !LinearFeatureDetect object
      integer, dimension(:), intent(in)                 ::      lidsIn  !array of line IDs that are of interest.
      integer, intent(in)                              ::      noofLineIn

      integer                                         ::      ii, lowestLid, currentVertNconnex, tempLId, iterCount, changeCount
      logical                                         ::      changesMade

      do ii = 1, noofLineIn
         lfdIn%L(lidsIn(ii))%subgraph = lidsIn(ii) !so there is a starting lowest value
      end do

      changesMade = .true.
      itercount = 0

      do while (changesMade)
         iterCount = iterCount + 1
         changeCount = 0
         do ii = 1, noofLineIn
            tempLId = getSmallestLId(lfdIn%L(lidsIn(ii))%from, lfdIn)
            if (tempLId < lfdIn%L(lidsIn(ii))%subgraph) then
               lfdIn%L(lidsIn(ii))%subgraph = tempLId
               changeCount = changeCount + 1
            end if
            !look at V2
            tempLId = getSmallestLId(lfdIn%L(lidsIn(ii))%to, lfdIn)
            if (tempLId < lfdIn%L(lidsIn(ii))%subgraph) then
               lfdIn%L(lidsIn(ii))%subgraph = tempLId
               changeCount = changeCount + 1

            end if

            !what is its connectivity
            !look at each line connected, find lowest ID
            !setlowsetLid to lowest
            !do the same for V2

         end do
         if (changeCount == 0) then
            changesMade = .false.
         else
            changesMade = .true.
         end if
         print *, "DBG: whatGraph iteration no: ", iterCount, changesMade, changeCount

      end do

   end subroutine whatSubGraph

   pure function getSmallestLId(vIdIn, lfdIn) result(lIdOut)
      !for a given vertex ID line the smallest line ID connected to it,
      !all accessed from linear feature detect object.
      type(LinearFeatureDetect), intent(in)            ::      lfdIn
      integer, intent(in)                              ::      vIdIn

      integer                                         ::      lIdOut, ii
      lIdOut = lfdIn%V(vIdIn)%connex(1)
      if (lfdIn%V(vIdIn)%nConnex > 1) then
         do ii = 1, lfdIn%V(vIdIn)%nConnex
            if (lfdIn%V(vIdIn)%connex(ii) < lIdOut) then
               lIdOut = lfdIn%V(vIdIn)%connex(ii)
            end if
         end do
      end if

   end function getSmallestLId

   subroutine outputPngLinesRGBsubgraph(lfdIn, lidsIn, Nxin, Nyin, noofLinesInterest, filenameIn)
      !reconstructs an image patch by patch colouring lines by sub graph ID
      integer, intent(in)             ::       Nxin, Nyin
      type(LinearFeatureDetect), intent(inout)    ::  lfdIn
      integer, dimension(:), intent(in)         ::  lidsIn
      integer, intent(in)                              ::      noofLinesInterest
      character(len=*), intent(in)                     ::      filenameIn        !   Filename in.
      real(kind=real64), dimension(:, :, :), allocatable     ::      img     !array of image

      integer         ::      ix, iy           !patch indices
      do iy = 0, lfdIn%My - 1
         do ix = 0, lfdIn%Mx - 1
            call getLinesInCellRGB(lfdIn, ix, iy, lidsIn, noofLinesInterest) !calls reconstructImage3
         end do
      end do

      allocate (img(3, 0:Nxin - 1, 0:Nyin - 1))
      call finalStitchCellsRGB(lfdIn, img)
      call write_rgb_png(filenameIn, img)
   end subroutine outputPngLinesRGBsubgraph

   subroutine getLinesInCellRGB(this, ix, iy, lidsIn, noofLinesInterest)
      type(LinearFeatureDetect), intent(inout)         ::      this    !LinearFeatureDetect object
      integer, intent(in)                              ::      ix, iy   !x & y patch indices
      integer, dimension(:), intent(in)                 ::      lidsIn
      integer, intent(in)                              ::      noofLinesInterest

      integer         ::      ii, nn           !Line indices, line counter
      integer         ::      jx, jy           !link cell indices
      type(Vertex), pointer    ::      vi, vj   !pointers to vertices i & j

      real(kind=real64), dimension(:, :), allocatable      ::      lines ! line(1:7) = (/x1,y1,x2-x1,y2-y1,r,i/)
      real(kind=real64), dimension(:), allocatable      ::      subgraphIds
      type(LinkCell), pointer  ::      cc  !pointer to LinkCell object                 !where sgid is sub graph id
      logical         ::      ok          !Are vertices i and j in the same LinkCell

      cc => this%c(ix, iy)

      nn = 0
      do ii = 1, noofLinesInterest
         vi => this%V(this%L(lidsIn(ii))%from)
         vj => this%V(this%L(lidsIn(ii))%to)
         do jy = min(vi%iy, vj%iy), max(vi%iy, vj%iy)
            do jx = min(vi%ix, vj%ix), max(vi%ix, vj%ix)
               if ((jx == ix) .and. (jy == iy)) then
                  nn = nn + 1
               end if
            end do
         end do

      end do

      if (nn == 0) return
      allocate (lines(6, nn))
      allocate (subgraphIds(nn))

      nn = 0
      do ii = 1, noofLinesInterest
         vi => this%V(this%L(lidsIn(ii))%from)
         vj => this%V(this%L(lidsIn(ii))%to)
         ok = .false.
         do jy = min(vi%iy, vj%iy), max(vi%iy, vj%iy)
            do jx = min(vi%ix, vj%ix), max(vi%ix, vj%ix)
               if ((jx == ix) .and. (jy == iy)) then
                  ok = .true.
               end if
            end do
         end do
         if (ok) then

            nn = nn + 1

            lines(1:2, nn) = vi%p(1:2)

            lines(3:4, nn) = vj%p(1:2) - lines(1:2, nn)

            lines(1, nn) = lines(1, nn) - (this%jx(ix) - 1)
            lines(2, nn) = lines(2, nn) - (this%jy(iy) - 1)

            lines(5, nn) = this%L(lidsIn(ii))%r
            lines(6, nn) = this%L(lidsIn(ii))%f
            subgraphIds(nn) = this%L(lidsIn(ii))%subgraph

         end if

      end do

      call reconstructImageRGB(lines, cc%recon_img, subgraphIds, cc%recon_img_RGB)!imgRGB

      deallocate (lines)

      return
   end subroutine getLinesInCellRGB

   subroutine finalStitchCellsRGB(this, img)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      stitch together image from all the link cells
      type(LinearFeatureDetect), intent(inout)         ::      this    !LinearFeatureDetect object
      real(kind=real64), dimension(:, :, :), intent(out)    ::      img     !array holding reconstructed image

      integer         ::      ix, iy       !patch indices
      integer         ::      jx, jy       !patch pixel coordinates
      type(LinkCell), pointer  ::      cc  !poniter to LinkCell object

      !---    copy blocks of image back into original
      do iy = 0, this%My - 1
         do ix = 0, this%Mx - 1
            cc => this%c(ix, iy)
            jx = this%jx(ix)
            jy = this%jy(iy)
            img(:, jx:jx + cc%w - 1, jy:jy + cc%h - 1) = cc%recon_img_RGB(:, 0:cc%w - 1, 0:cc%h - 1)
         end do
      end do
      return
   end subroutine finalStitchCellsRGB

   subroutine dbg_print_line_subgraphs(nooflines, lidsIn, lfdIn)
      !loops through all graphs and prints any composed of more than one line.
      integer, intent(in)                      ::      nooflines
      integer, dimension(:), intent(in)         ::      lidsIn
      type(LinearFeatureDetect), intent(in)    ::      lfdIn    !LinearFeatureDetect object
      integer                                 ::      ii, jj, count
      integer, dimension(nooflines)            ::      holder
      print *, "dbg_print_multi_seg_lines starting"
      do ii = 1, nooflines
         holder(ii) = 0
         print *, "DBG line ID is:", (lfdIn%L(lidsIn(ii))%id), "subgraph is:", lfdIn%L(lidsIn(ii))%subgraph
      end do

      count = 0

      do ii = 1, nooflines
         do jj = 1, nooflines
            if (ii == jj) cycle
            if (lfdIn%L(lidsIn(ii))%subgraph == lfdIn%L(lidsIn(jj))%subgraph) then
               holder(lfdIn%L(lidsIn(jj))%subgraph) = holder(lfdIn%L(lidsIn(jj))%subgraph) + 1
            end if
         end do
      end do
      count = maxval(holder)
      !print*,"most common subgraph is",(maxloc(holder)),"with",count,"segments"
   end subroutine dbg_print_line_subgraphs

   subroutine dbg_report_RGB_patch_size(lfdIn, ix, iy)
      !reports the size of the RGB reconstructed image asocated wit hthe cell ix,iy
      !Also states jx,jy (patch position) for debuggging
      type(LinearFeatureDetect), intent(in)    ::      lfdIn    !LinearFeatureDetect object
      integer, intent(in)                      ::      ix, iy
    print *, "DBG combi RGB patch(11,12) jx, jy, size:", lfdIn%c(11, 12)%jx, lfdIn%c(11, 12)%jy, size(lfdIn%c(11, 12)%recon_img_RGB)
    print *, "DBG combi RGB patch(12,12) jx, jy, size:", lfdIn%c(12, 12)%jx, lfdIn%c(12, 12)%jy, size(lfdIn%c(12, 12)%recon_img_RGB)
    print *, "DBG combi RGB patch(12,11) jx, jy, size:", lfdIn%c(12, 11)%jx, lfdIn%c(12, 11)%jy, size(lfdIn%c(12, 11)%recon_img_RGB)

      print *, "DBG combi patch(11,12) jx, jy, size:", lfdIn%c(11, 12)%jx, lfdIn%c(11, 12)%jy, size(lfdIn%c(11, 12)%recon_img)
      print *, "DBG combi patch(12,12) jx, jy, size:", lfdIn%c(12, 12)%jx, lfdIn%c(12, 12)%jy, size(lfdIn%c(12, 12)%recon_img)
      print *, "DBG combi patch(12,11) jx, jy, size:", lfdIn%c(12, 11)%jx, lfdIn%c(12, 11)%jy, size(lfdIn%c(12, 11)%recon_img)
   end subroutine dbg_report_RGB_patch_size

   subroutine dbg_report_zerosum_patch(lfdin, marker)
      !checks all patches (greyscale and RGB) and prints the patch indices is all elemnts sum to zero
      !marker is an int that will also printed to mark where this is called
      type(LinearFeatureDetect), intent(in)    ::      lfdIn    !LinearFeatureDetect object
      integer, intent(in)                      ::      marker

      integer                                 ::      ii, jj, kk, ll
      real(kind=real64)                       ::      sumhold

      print *, "dbg_report_zerosum_patch call mark:", marker

      do ii = 0, lfdIn%Mx - 1
         do jj = 0, lfdIn%My - 1
            sumhold = 0.0d0
            sumhold = sum(lfdin%c(ii, jj)%recon_img)
            if (sumhold == 0.0d0) then
               print *, "zero sum greyscale patch ix,iy", ii, jj
            end if
         end do
      end do

      do ii = 0, lfdIn%Mx - 1
         do jj = 0, lfdIn%My - 1
            sumhold = 0.0d0
            sumhold = sum(lfdin%c(ii, jj)%recon_img_RGB)
            if (sumhold == 0.0d0) then
               print *, "zero sum RGB patch ix,iy", ii, jj
            end if
         end do
      end do

   end subroutine dbg_report_zerosum_patch

   subroutine applyButterflyMask(sinoIn, maskIn, butteredSinoOut)
      !applys a mxm pixel butterfly mask to a sinogram, the image is assumed to be periodic
      !for pixels close to the boundary.
      real(kind=real64), dimension(:, :), intent(in)         ::  maskIn
      real(kind=real64), dimension(-107:107, 0:118), intent(in)         ::  sinoIn
      real(kind=real64), dimension(-107:107, 0:118), intent(inout)        ::  butteredSinoOut

      integer, dimension(2, 2)                              ::  sinogramBounds
      integer                                             ::  Nx, Ny, ix, iy !sino xy dimensions and indices
      integer                                             ::  Mx, My, jx, jy, px, py !mask xy dimensions and indices
      integer             ::      maxy, nTheta                 !Half the y dimension of the sinogram/pixels
      real(kind=real64)                                   ::  maxBut, minBut

      maxy = int((size(sinoIn, dim=1) - 1)/2)
      nTheta = size(sinoIn, dim=2) - 1

      Nx = 2*maxy - 1!size(sinoIn,dim=1)
      Ny = nTheta - 1!size(sinoIn,dim=2)
      Mx = size(maskIn, dim=1)
      My = size(maskIn, dim=2)

      sinogramBounds(1, 1) = -maxy !sinogram rho min
      sinogramBounds(1, 2) = maxy !sinogram rho max
      sinogramBounds(2, 1) = 0 !sinogram theta min
      sinogramBounds(2, 2) = nTheta !sinogram theta max

      do iy = 0, nTheta - 1
         do ix = -maxy, maxy!0,2*maxy
            butteredSinoOut(ix, iy) = getMaskValue(sinoIn, maskIn, ix, iy, Mx, My, Nx, Ny, sinogramBounds)
         end do
      end do
      !a=(h-np.min(h))/(np.max(h)-np.min(h))
      maxBut = maxval(butteredSinoOut)
      minBut = minval(butteredSinoOut)
      !butteredSinoOut=(butteredSinoOut-minval(butteredSinoOut))/(maxval(butteredSinoOut)-minval(butteredSinoOut))
      butteredSinoOut = (butteredSinoOut - minBut)/(maxBut - minBut)
   end subroutine applyButterflyMask

   !pure function  getMaskValue(arrayIn,maskIn,ixIn,iyIn,MxIn,MyIn,NxIn,NyIn,sinoBoundsIn) result(maskValueOut)
   function getMaskValue(arrayIn, maskIn, ixIn, iyIn, MxIn, MyIn, NxIn, NyIn, sinoBoundsIn) result(maskValueOut)
      !Calculates value of pixel iy ix of input array for given mask
      !if the any of the mask array is not within the input array then vlues are obtained periodically.
      real(kind=real64), dimension(-4:4, -4:4), intent(in)         ::  maskIn
      real(kind=real64), dimension(0:, 0:), intent(in)         ::  arrayIn
      real(kind=real64)                                   ::  maskValueOut
      integer, dimension(2, 2), intent(in)                   ::  sinoBoundsIn !sinogram rho min, sinogram rho max; sinogram theta min, sinogram theta max
      integer, intent(in)                                  ::  ixIn, iyIn !idices of pixel of interest
      integer, intent(in)                                  ::  MxIn, MyIn !xy dimensions of the mask
      integer, intent(in)                                  ::  NxIn, NyIn  !xy dimensions of the input array
      integer                                             ::  jx, jy !mask xy indices
      integer                                             ::  tempx, tempy
      integer                                             ::  Nx, Ny   !arrayIn xy dimensions
      integer                                             ::  offset !offset for the kernel
      real(kind=real64), dimension(-4:4, -4:4)      ::  ff              !Sub section of sinogram

      offset = (MxIn - 1)/2 !kernel must be odd
      maskValueOut = 0.0d0 !arrayIn(ixIn,iyIn)

      do jx = -4, 4

         tempx = ixIn + jx !tempx=ixIn-jx
         if (tempx < sinoBoundsIn(1, 1)) then
            tempx = tempx + sinoBoundsIn(1, 2)
         else if ((ixIn + jx) > sinoBoundsIn(1, 2)) then
            tempx = ixIn + jx - sinoBoundsIn(1, 2)
         else
            tempx = ixIn + jx
         end if

         do jy = -4, 4

            tempy = iyIn + jy !tempy=iyIn-jy
            if (tempy < sinoBoundsIn(2, 1)) then
               tempy = tempy + sinoBoundsIn(2, 2)
            else if ((iyIn + jy) > sinoBoundsIn(2, 2)) then
               tempy = iyIn + jy - sinoBoundsIn(2, 2)
            else
               tempy = iyIn + jy
            end if

            maskValueOut = maskValueOut + arrayIn(tempx, tempy)*maskIn(jx, jy)

         end do
      end do
      maskValueOut = maskValueOut/sum(maskIn)

      return
   end function getMaskValue

   subroutine return9by9ButterflyMask(butterflyOut)
      real(kind=real64), dimension(:, :), allocatable, intent(out)    ::  butterflyOut !array holding the weights

      allocate (butterflyOut(9, 9))
      butterflyOut = transpose(reshape((/-10.0d0, -15.0d0, -22.0d0, -22.0d0, -22.0d0, -22.0d0, -22.0d0, -15.0d0, -10.0d0, &
                                         1.0d0, -6.0d0, -13.0d0, -22.0d0, -22.0d0, -22.0d0, -13.0d0, -6.0d0, 1.0d0, &
                                         3.0d0, 6.0d0, 4.0d0, -3.0d0, -22.0d0, -3.0d0, 4.0d0, 6.0d0, 3.0d0, &
                                         3.0d0, 11.0d0, 19.0d0, 28.0d0, 42.0d0, 28.0d0, 19.0d0, 11.0d0, 3.0d0, &
                                         3.0d0, 11.0d0, 27.0d0, 42.0d0, 42.0d0, 42.0d0, 27.0d0, 11.0d0, 3.0d0, &
                                         3.0d0, 11.0d0, 19.0d0, 28.0d0, 42.0d0, 28.0d0, 19.0d0, 11.0d0, 3.0d0, &
                                         3.0d0, 6.0d0, 4.0d0, -3.0d0, -22.0d0, -3.0d0, 4.0d0, 6.0d0, 3.0d0, &
                                         1.0d0, -6.0d0, -13.0d0, -22.0d0, -22.0d0, -22.0d0, -13.0d0, -6.0d0, 1.0d0, &
                            -10.0d0, -15.0d0, -22.0d0, -22.0d0, -22.0d0, -22.0d0, -22.0d0, -15.0d0, -10.0d0/), shape(butterflyOut)))
      !may need to transpose
      !array = transpose(reshape((/ 1, ..., 9 /), shape(array)))
   end subroutine return9by9ButterflyMask

   subroutine cellLineRecon(this, lid, cut, ix, iy)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      add or remove single line.
      type(LinearFeatureDetect), intent(inout)         ::      this    !LinearFeatureDetect object
      integer, intent(in)                              ::      lid     !line ID
      logical, intent(in)                              ::      cut     !is this line to be cut?
      integer, intent(in)                              ::      ix, iy
      real(kind=real64), dimension(6)      ::      line   !array holding data for a single line
      !patch coordinates
      type(Vertex), pointer        ::      vi, vj           !pointers to vertices i & j

      if (this%L(lid)%id == LIB_LFD_CUT) return

      vi => this%V(this%L(lid)%from)
      vj => this%V(this%L(lid)%to)

      if (min(vi%ix, vj%ix) == -1) return        !   this vertex has been cut. Odd that the line hasn't ???

      line(1:2) = vi%p(1:2)

      line(3:4) = vj%p(1:2) - line(1:2)

      !---    remove cell offset

      line(5) = this%L(lid)%r
      if (cut) then
         line(6) = -this%L(lid)%f         ! the -ve makes it cut
      else
         line(6) = +this%L(lid)%f         ! the +ve makes it add
      end if

      line(1:2) = line(1:2) - (/this%jx(ix) - 1, this%jy(iy) - 1/)
      call reconstructImage(line, this%c(ix, iy)%recon_img)

      return
   end subroutine cellLineRecon

   subroutine whatCellIsVetex(this, p, ix, iy)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      For a vertex at point p, returns what cell it should be in.
      type(LinearFeatureDetect), intent(inout)             ::      this    !LinearFeatureDetect object
      integer, intent(out)     ::  ix, iy    !output patch x & y indices

      real(kind=real64), dimension(2), intent(in)           ::      p       !Coordinates of vertex to be queried
      !integer,intent(out)                                 ::      vid     !vertex ID
      type(Vertex), dimension(:), pointer                   ::      V_tmp   !temporary pointer to vertex object

      integer     ::  xx, yy, ii    !patch x & y, vertex indices

      do xx = 0, this%Mx - 1
         if (this%jx(xx + 1) > p(1)) then
            do yy = 0, this%My - 1
               if (this%jy(yy + 1) > p(2)) then
                  !   link cell is ix,iy

                  ix = xx
                  iy = yy

                  return
               end if
            end do
         end if
      end do

      !---    if we get here, there's a problem
      print *, "Lib_LinearFeatureDetect::addVertex0 error - p = ", p, " Nx,Ny ", this%Nx, this%Ny
      stop "addVertex0 error -  vertex is out of bounds"

      return
   end subroutine whatCellIsVetex

   subroutine reconImageSingleLineAdjacentCells(lfdIn, lidIn, cut, lineBounds)
      !Takes in a line id, finds which cells bound it and redraws it in the bounding cells
      !and an additional one cell layer round it.

      !The redraw is done on a single image composed of all affected cells stitched together, the stich
      !is then broken back down into sub images and fed back to each respective link cell. Boolean "cut"
      !determins whether the line's pixel intensity contribution is subtracted (TRUE) or subtracted
      !(FALSE). linebounds is returned so the caluclation of which cells bound the line can bbe reused
      !caluclating the RSS.

      !Later it may be worth comparing the cell dimensions with the line thickness and orientation to
      !decide if is necessary to have the padding cells or not.
      integer, intent(in)                                  ::      lidIn       !line id of line_ij
      logical, intent(in)                                  ::      cut         !cut out line intensity?
      type(LinearFeatureDetect), intent(inout)             ::      lfdIn       !LinearFeatureDetect object
      integer, dimension(4), intent(in)                    ::      lineBounds  !cells bounding the line
      !in form ix_min,iy_min,ix_max,iy_max
      !type(Vertex),pointer    ::      vi,vj           !pointers to vertices i & j
      !integer,dimension(6)    ::      lineDat            !array of line data(x1,y1,x2-x1,y1-y2,r,i)
      real(kind=real64), dimension(:, :), allocatable        ::      imgStitch   !Array holding recon images stiched together

      integer                 ::      ix, iy           !link cell xy indices.

      do ix = max(0, lineBounds(1) - 1), min(lfdIn%Mx - 1, lineBounds(3) + 1)
         do iy = max(0, lineBounds(2) - 1), min(lfdIn%My - 1, lineBounds(4) + 1)
            call cellLineRecon(lfdIn, lidIn, cut, ix, iy)
         end do
      end do

   end subroutine reconImageSingleLineAdjacentCells

   subroutine getCellsBoundingLine(lfdIn, lidIn, lineBounds)
      !Takes in a line id, finds which cells bound it.
      integer, intent(in)                                  ::      lidIn       !line id of line_ij
      type(LinearFeatureDetect), intent(inout)             ::      lfdIn       !LinearFeatureDetect object
      integer, dimension(4), intent(out)                    ::      lineBounds  !cells bounding the line
      type(Vertex), pointer    ::      vi, vj           !pointers to vertices i & j

      vi => lfdIn%V(lfdIn%L(lidIn)%from)
      vj => lfdIn%V(lfdIn%L(lidIn)%to)

      lineBounds(1) = min(vi%ix, vj%ix)!get bounding cells
      lineBounds(2) = min(vi%iy, vj%iy)
      lineBounds(3) = max(vi%ix, vj%ix)
      lineBounds(4) = max(vi%iy, vj%iy)

   end subroutine getCellsBoundingLine

   subroutine getCellsBoundingVertexLines(lfdIn, vidIn, lineBounds)
      !For a given vertex finds the bounding cells for all connected lines.
      integer, intent(in)                                  ::      vidIn !vertex id
      type(LinearFeatureDetect), intent(inout)             ::      lfdIn       !LinearFeatureDetect object
      integer, dimension(4), intent(out)                    ::      lineBounds  !cells bounding the line
      type(Vertex), pointer    ::      vi, vj           !pointers to vertices i & j
      integer                                             ::      ii          !line index.

      lineBounds(1) = lfdIn%Mx - 1 !largest xy vals the minimum could be
      lineBounds(2) = lfdIn%My - 1
      lineBounds(3) = 0 !smallest xy vals the maximum could be
      lineBounds(4) = 0

      do ii = 1, lfdIn%V(vidIn)%nConnex
         if (lfdIn%L(lfdIn%V(vidIn)%connex(ii))%id == LIB_LFD_CUT) CYCLE
         vi => lfdIn%V(lfdIn%L(lfdIn%V(vidIn)%connex(ii))%from)
         vj => lfdIn%V(lfdIn%L(lfdIn%V(vidIn)%connex(ii))%to)

         lineBounds(1) = min(lineBounds(1), (min(vi%ix, vj%ix))) !get bounding cells
         lineBounds(2) = min(lineBounds(2), (min(vi%iy, vj%iy)))
         lineBounds(3) = max(lineBounds(3), (max(vi%ix, vj%ix)))
         lineBounds(4) = max(lineBounds(4), (max(vi%iy, vj%iy)))
      end do

   end subroutine getCellsBoundingVertexLines

   subroutine getCellsBoundingMultiLines(lfdIn, lidIn, lineBounds)
      !Takes in a line id, finds which cells bound it.
      integer, intent(in)                                  ::      lidIn       !line id of line_ij
      type(LinearFeatureDetect), intent(inout)             ::      lfdIn       !LinearFeatureDetect object
      integer, dimension(4), intent(out)                    ::      lineBounds  !cells bounding the line
      type(Vertex), pointer    ::      vi, vj           !pointers to vertices i & j

      vi => lfdIn%V(lfdIn%L(lidIn)%from)
      vj => lfdIn%V(lfdIn%L(lidIn)%to)

      lineBounds(1) = min(vi%ix, vj%ix)!get bounding cells
      lineBounds(2) = min(vi%iy, vj%iy)
      lineBounds(3) = max(vi%ix, vj%ix)
      lineBounds(4) = max(vi%iy, vj%iy)

   end subroutine getCellsBoundingMultiLines

   subroutine boundingCellsAIC(lfdIn, lineBounds, AICout)
      !Calcualtes AIC for a bounding box of cells defined by min_x,min_y,max_x,max_y
      type(LinearFeatureDetect), intent(inout)             ::      lfdIn       !LinearFeatureDetect object
      integer, dimension(4), intent(in)                    ::      lineBounds  !cells bounding the line
      real(kind=real64), intent(out)                      ::       AICout

      real(kind=real64), parameter   ::      TWOPI = 6.283185307d0
      integer                         ::          ix, iy   !linkcell indices
      integer                         ::          ii      !vertex index
      integer                         ::          nDoFs, vid   !number of degrees of freedom. vertex ID
      real(kind=real64)               ::          tRSS, tArea    !total RSS, area of all cells
      real(kind=real64)               ::          chisquare_dof, reduced_chisquare_stat
      !call findRss(lfdIn%c(ix,iy))
      nDoFs = 0
      tRSS = 0.0d0

      do ix = lineBounds(1), linebounds(3)
         do iy = linebounds(2), linebounds(4)
            !get RSS
            tRSS = tRSS + lfdIn%c(ix, iy)%rss
            !Count Noof DoFs

            do ii = 1, lfdIn%c(ix, iy)%nV
               vid = lfdIn%c(ix, iy)%id(ii)

               if (vid /= LIB_LFD_CUT) then
                  nDoFs = nDoFs + 2                         !   2 for vertex location
                  nDoFs = nDoFs + lfdIn%V(vid)%nConnex       !   1/2 * 2 for line width intensity

               end if
            end do
            nDoFs = nDoFs + 1 !add background DoF for each cell
            tArea = tArea + lfdIn%c(ix, iy)%w*lfdIn%c(ix, iy)%h
         end do
      end do
      chisquare_dof = tArea - 1
      reduced_chisquare_stat = tRSS/chisquare_dof       !   rss per degree of freedom
      AICout = 2*nDoFs + tArea*log(TWOPI*reduced_chisquare_stat) + chisquare_dof

   end subroutine boundingCellsAIC

   subroutine findRssBoundingCells(lfdIn, lineBounds)
      !---^^^^^^^^^^^^^^^^^^^^^^^^
      !*      compute rss for a bounding box of link cells
      type(LinearFeatureDetect), intent(inout)             ::      lfdIn       !LinearFeatureDetect object
      integer, dimension(4), intent(in)                    ::      lineBounds  !cells bounding the line
      integer                 ::      ix, iy               !x & y Link cell indices
      integer                 ::      kx, ky               !x & y Link cell pixel indices
      real(kind=real64)       ::      dd                  !difference between reconstructed image and input image.

      do ix = lineBounds(1), linebounds(3)
         do iy = linebounds(2), linebounds(4)
            lfdIn%c(ix, iy)%rss = 0
            do ky = 1, lfdIn%c(ix, iy)%h
               do kx = 1, lfdIn%c(ix, iy)%w
                  dd = lfdIn%c(ix, iy)%recon_img(kx, ky) - lfdIn%c(ix, iy)%buffered_img(kx, ky)
                  lfdIn%c(ix, iy)%rss = lfdIn%c(ix, iy)%rss + dd*dd
               end do
            end do
         end do
      end do
      return
   end subroutine findRssBoundingCells

   subroutine balanceIntensityBoundingCells(lfdIn, lineBounds)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      set the average intensity level for a bounding box of link cells
      type(LinearFeatureDetect), intent(inout)             ::      lfdIn       !LinearFeatureDetect object
      integer, dimension(4), intent(in)                    ::      lineBounds  !cells bounding the line
      integer                 ::      ix, iy               !x & y Link cell indices
      real(kind=real64)               ::      img_bar             !average patch intensity

      do ix = lineBounds(1), linebounds(3)
         do iy = linebounds(2), linebounds(4)
            img_bar = sum(lfdIn%c(ix, iy)%recon_img)/(lfdIn%c(ix, iy)%w*lfdIn%c(ix, iy)%h)
            lfdIn%c(ix, iy)%recon_img = lfdIn%c(ix, iy)%recon_img + lfdIn%c(ix, iy)%bimg_bar - img_bar
         end do
      end do

      return
   end subroutine balanceIntensityBoundingCells

   subroutine vertexShiftReconBoundingCells(this, vid, pcut, padd, lineBounds)
      !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !*      remove single vertex vid from pbefore and all its lines, and replace at after.
      !*      note: does not update the positions of the vertices. Calls bounding cells reconstruction
      type(LinearFeatureDetect), intent(inout)         ::      this    !LinearFeatureDetect object
      integer, intent(in)                              ::      vid     !vertex ID
      real(kind=real64), dimension(2), intent(in)       ::      pcut, padd   !updated vertex coordinates
      integer, dimension(4), intent(in)                    ::      lineBounds  !cells bounding the line
      integer         ::      ii, lid  !connection index, line ID

      !   move vertex to cut position and remove lines from image
      this%V(vid)%p(1:2) = pcut(1:2)
      do ii = 1, this%V(vid)%nConnex                       !   loop through all the lines connecting to the vertex
         lid = this%V(vid)%connex(ii)
         call reconImageSingleLineAdjacentCells(this, lid, .true., lineBounds)
      end do

      !   move vertex to add position and add lines to image
      this%V(vid)%p(1:2) = padd(1:2)
      do ii = 1, this%V(vid)%nConnex                   !   loop through all the lines connecting to the vertex
         lid = this%V(vid)%connex(ii)
         call reconImageSingleLineAdjacentCells(this, lid, .false., lineBounds)
      end do

      !   put the vertex back
      this%V(vid)%p(1:2) = pcut(1:2)

      return
   end subroutine vertexShiftReconBoundingCells

end module Lib_LinearFeatureDetect

