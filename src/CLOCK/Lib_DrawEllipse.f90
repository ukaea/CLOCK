
    module Lib_DrawEllipse
!---^^^^^^^^^^^^^^^^^^^^^^
!*      A simple module which draws a 2d  ellipse into a box
!*
!*      _______________
!*      |        __/ || 
!*      |     __/    /| 
!*      |   _/     _/ | 
!*      | _/     _/   | 
!*      |/    __/     | 
!*      ||__/_________|                                            
!*      
!*      Daniel Mason
!*      (c) UKAEA August 2024
!*
!*      version history
!*      v.0.0.1     Aug 2024                first version
!*      v.0.1.0     Mar 2025                Added support for getRss,findDrss,add,getT,getWeight
!*
!*      The general ellipse equation used here is  x.Dx = c
!*      and then shading is either undertaken by 
!*      using a Gaussian profile of the                     f = exp( - c )
!*      or a ring can be drawn at intensity level c0 = 1
!*      or a flat object can be drawn with intensity profile flat for c>c0
!*
!*      recall that in normal graphics output theta is counted clockwise because y=0 is the top row of pixels.
!*
!*          +-------------------->  x
!*          |\  theta
!*          | \ / 
!*          |  \
!*          |   \
!*          |    \
!*          |     \                               
!*          v
!*      
!*          y
!*
!*      A key routine is to _blend_ an ellipse with the background
!*          drawEllipse(D,p,img,mode,f,ignore,dellipse)
!*      with    p = (fractional) offset in pixels
!*              img = (inout) image buffer
!*              mode = LIB_DRAWELLIPSE_SHADE_GAUSSIAN / RING / FLAT
!*              f = (optional) intensity. 
!*              ignore = (optional) what to do with background. If true, background is set to LIB_DRAWELLIPSE_IGNORE outside of ellipse
!*              dellipse = (optional) derivative of img wrt [ x0,y0,f0,Dxx,Dxy,Dyy ] 
!*      
!*      or to _add_ an ellipse into the background
!*           add(D,p,img,mode,f,mult)
!*
!*      Another key routine 
!*          getRss( ellipse,img_in,mode )
!*      takes in all the ellipse information as a 6-vector
!*          ellipse = [ x0,y0,f0,Dxx,Dxy,Dyy ] 
!*      and finds the residual sum of squares if this ellipse is removed from the image.
!*      ( note- only the pixels within the ellipse are included, to a range LIB_DRAWELLIPSE_GAUSSMULT )
!*



        use Lib_SafeExp
        use Lib_ColourScale
        use iso_fortran_env
        implicit none
        private

        real(kind=real64),private,parameter     ::      PI = 3.14159265390d0
        integer(kind=int64),private,parameter   ::      BADF00D = int( z'BADF00D',kind=int64 )
        real(kind=real64),public,parameter      ::      LIB_DRAWELLIPSE_IGNORE = transfer( (BADF00D+ishft(BADF00D,32_int64)),1.0d0 )


        integer,public,parameter                ::      LIB_DRAWELLIPSE_SHADE_GAUSSIAN = 1
        integer,public,parameter                ::      LIB_DRAWELLIPSE_SHADE_RING     = 2
        integer,public,parameter                ::      LIB_DRAWELLIPSE_SHADE_FLAT     = 3

        integer,public,parameter                ::      LIB_DRAWELLIPSE_IOU_RECTANGLE  = 1
        integer,public,parameter                ::      LIB_DRAWELLIPSE_IOU_CIRCLE     = 2
        integer,public,parameter                ::      LIB_DRAWELLIPSE_IOU_ELLIPSE    = 3

        logical,public                          ::      LIB_DRAWELLIPSE_RSS_REMOVE_BG = .true.
        real(kind=real64),public                ::      LIB_DRAWELLIPSE_GAUSSMULT = 5           !   draw gaussian to 5 std devs
        real(kind=real64),public                ::      LIB_DRAWELLIPSE_RINGWIDTH = 1.0d0       !   draw ring with 1 pixel width
        

        public          ::      getDMatrix                  !   construct a D matrix from major/minor radii and angle
        public          ::      findBoundingBox             
        public          ::      drawEllipse
        public          ::      intersectionOverUnion
        public          ::      getRss
        public          ::      findDrss
        public          ::      add
        public          ::      getT
        public          ::      getWeight
        public          ::      do_spots_overlap
        public          ::      check_overlap_with_prexisting_spots   
        

        interface       findBoundingBox
            module procedure    findBoundingSquare1         !   find a square large enough to contain one ellipse - note overestimate for allocating image array
            module procedure    findBoundingBox2            !   find a rectangle large enough to contain two ellipse - note overestimate for allocating image array
            module procedure    findBoundingRectangle3      !   find minimum bounding rectangle containing pixels of intensity >= 0.5
        end interface

        interface       drawEllipse
            module procedure    drawEllipse_greyscale
            module procedure    drawEllipse_RGB
        end interface



        interface       getRss
            module procedure    getRss0
            module procedure    getRss1
        end interface


        interface       findDrss
            module procedure    findDrss0
            module procedure    findDrss1
        end interface        

        interface       add
            module procedure    add0
            module procedure    add1
        end interface        
 
        interface       getT
            module procedure    getT0
        end interface        

        interface       getWeight
            module procedure    getWeight0
        end interface        

    contains
!---^^^^^^^^

        subroutine drawEllipse_greyscale(D,p,img,mode,f,ignore,dellipse)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that D is in pixels
    !*      mode = LIB_DRAWELLIPSE_SHADE_GAUSSIAN etc describes the drawing mode to use
    !*      if ignore, then fills the image with LIB_DRAWELLIPSE_IGNORE outside gaussian
    !*      The ordering of the output derivatives in dellipse
    !*          [ x0,y0,f0,Dxx,Dxy,Dyy ]            
    !*      see also add, which directly adds intensity onto img instead of blending
    
            real(kind=real64),dimension(2,2),intent(in)                 ::      D
            real(kind=real64),dimension(2),intent(in)                   ::      p                           !   fractional pixel offset of centre  
            real(kind=real64),dimension(0:,0:),intent(inout)            ::      img
            integer,intent(in)                                          ::      mode
            real(kind=real64),intent(in),optional                       ::      f                           !   intensity peak (needed for absolute image intensity and for deriv)
            logical,intent(in),optional                                 ::      ignore
            real(kind=real64),dimension(:,:,:),allocatable,intent(out),optional       ::      dellipse      !   derivative of image wrt ellipse parameters


            real(kind=real64),dimension(:,:),allocatable                ::      img_blit
            real(kind=real64),dimension(:,:,:),allocatable              ::      dimg_blit
            real(kind=real64)   ::      xx
            integer         ::      i0,j0           !   pixel at centre of ellipse
            integer         ::      nx,ny           !   half width of ellipse image
            integer         ::      ii,jj
            integer         ::      Mx,My           !   size of img
            logical         ::      ignore_
            
            ignore_ = .false. ; if (present(ignore)) ignore_ = ignore

        !---    find the pixel centre of the ellipse
            i0 = floor(p(1))
            j0 = floor(p(2))

        !---    construct the image that will be added
            select case(mode)
                case (LIB_DRAWELLIPSE_SHADE_RING)
                    call drawRingProfile(D,p(1)-i0,p(2)-j0,img_blit)
                    if (present(dellipse)) stop "Lib_DrawEllipse::drawEllipse_greyscale error - have not coded deriv for ring"
                case (LIB_DRAWELLIPSE_SHADE_FLAT)
                    call drawFlatProfile(D,p(1)-i0,p(2)-j0,img_blit)
                    if (present(dellipse)) stop "Lib_DrawEllipse::drawEllipse_greyscale error - have not coded deriv for flat"
                case default
                    if (present(dellipse)) then
                        call drawGaussianProfile(D,p(1)-i0,p(2)-j0,img_blit,ignore = .true.,dellipse = dimg_blit)
                    else 
                        call drawGaussianProfile(D,p(1)-i0,p(2)-j0,img_blit,ignore = ignore_)
                    end if
            end select

            nx = ubound(img_blit,dim=1)
            ny = ubound(img_blit,dim=2)
            Mx = size(img,dim=1)
            My = size(img,dim=2)

            if (present(dellipse)) allocate(dellipse(6,0:Mx-1,0:My-1))
 

        !---    add the image
            if (present(f)) then
                if (present(dellipse)) then
                    do jj = max(0,j0-ny),min(My-1,j0+ny)
                        do ii = max(0,i0-nx),min(Mx-1,i0+nx)
                            xx = img_blit( ii-i0 , jj-j0 )
                            if (ignore_ .and. (xx==LIB_DRAWELLIPSE_IGNORE) ) cycle
                            img(ii,jj) = img(ii,jj)*(1-xx) + f*xx 
                            dellipse(1,ii,jj) = dellipse(1,ii,jj) + f*dimg_blit(1,ii-i0,jj-j0)
                            dellipse(2,ii,jj) = dellipse(2,ii,jj) + f*dimg_blit(2,ii-i0,jj-j0)
                            dellipse(3,ii,jj) = dellipse(3,ii,jj) +   dimg_blit(3,ii-i0,jj-j0)  !   note deriv 3 is wrt f 
                            dellipse(4,ii,jj) = dellipse(4,ii,jj) + f*dimg_blit(4,ii-i0,jj-j0)  
                            dellipse(5,ii,jj) = dellipse(5,ii,jj) + f*dimg_blit(5,ii-i0,jj-j0)
                            dellipse(6,ii,jj) = dellipse(6,ii,jj) + f*dimg_blit(6,ii-i0,jj-j0)
                        end do
                    end do
                else
                    do jj = max(0,j0-ny),min(My-1,j0+ny)
                        do ii = max(0,i0-nx),min(Mx-1,i0+nx)
                            xx = img_blit( ii-i0 , jj-j0 )
                            if (ignore_ .and. (xx==LIB_DRAWELLIPSE_IGNORE) ) cycle
                            img(ii,jj) = img(ii,jj)*(1-xx) + f*xx 
                        end do
                    end do
                end if
            else
                do jj = max(0,j0-ny),min(My-1,j0+ny)
                    do ii = max(0,i0-nx),min(Mx-1,i0+nx)
                        xx = img_blit( ii-i0 , jj-j0 )
                        if (ignore_ .and. (xx==LIB_DRAWELLIPSE_IGNORE) ) cycle
                        img(ii,jj) = img(ii,jj)*(1-xx) + xx 
                    end do
                end do
            end if

            return
        end subroutine drawEllipse_greyscale

        subroutine drawEllipse_RGB(D,p,img,colourscale,f,mode)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that D is in pixels
    !*      mode = LIB_DRAWELLIPSE_SHADE_GAUSSIAN etc describes the drawing mode to use
    
            real(kind=real64),dimension(2,2),intent(in)                 ::      D
            real(kind=real64),dimension(2),intent(in)                   ::      p               !   fractional pixel offset of centre  
            real(kind=real64),dimension(:,0:,0:),intent(inout)          ::      img
            integer,intent(in)                                          ::      mode
            integer,intent(in)                                          ::      colourscale
            real(kind=real64),intent(in)                                ::      f               !   colour (0:1)

            real(kind=real64),dimension(3)                              ::      rgb_back,rgb_fore       !RGB triplet for background and foregronud respectively.
            real(kind=real64),dimension(:,:),allocatable                ::      img_blit
            real(kind=real64)   ::      xx
            integer         ::      i0,j0,nx,ny,Mx,My
            integer         ::      ii,jj

        !---    find the pixel centre of the ellipse
            i0 = floor(p(1))
            j0 = floor(p(2))

        !---    construct the image that will be added
            select case(mode)
                case (LIB_DRAWELLIPSE_SHADE_RING)
                    call drawRingProfile(D,p(1)-i0,p(2)-j0,img_blit)
                case (LIB_DRAWELLIPSE_SHADE_FLAT)
                    call drawFlatProfile(D,p(1)-i0,p(2)-j0,img_blit)
                case default
                    call drawGaussianProfile(D,p(1)-i0,p(2)-j0,img_blit)
            end select

            nx = ubound(img_blit,dim=1)
            ny = ubound(img_blit,dim=2)
            Mx = size(img,dim=2)
            My = size(img,dim=3)

            rgb_fore = getRGB_double( colourscale,f )       !   get RGB values for overlayed spot 
           
        !---    add the image
            do jj = max(0,j0-ny),min(My-1,j0+ny)
                do ii = max(0,i0-nx),min(Mx-1,i0+nx)
                    rgb_back = img(1:3,ii,jj)
                    xx = img_blit( ii-i0 , jj-j0 )
                    img(1:3,ii,jj) = transparentColour( rgb_back,rgb_fore,xx )
                    
                end do
            end do


            return
        end subroutine drawEllipse_RGB

!-------


        pure function getDMatrix( major_rad, minor_rad, theta, scale_factor ) result(D)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      produce a D matrix for this module representing an ellipse with an angle theta 
    !*      between the x axis and the major radius. 
    !*      also optionally add a scaling factor, so that ellipses will be drawn at multiples of the major/minor radii
    !*      note: doesn't check for major > minor, this is left to user.
            real(kind=real64),intent(in)                ::      major_rad,minor_rad
            real(kind=real64),intent(in)                ::      theta                   !in radians
            real(kind=real64),intent(in),optional       ::      scale_factor
            real(kind=real64),dimension(2,2)            ::      D
            real(kind=real64)           ::      lambda1,lambda2,ss,cost,sint

            ss = 1.0d0 ; if (present(scale_factor)) ss = scale_factor


            lambda1 = 1/(max(0.5d0,major_rad)*ss)
            lambda2 = 1/(max(0.5d0,minor_rad)*ss)
            lambda1 = lambda1*lambda1/2
            lambda2 = lambda2*lambda2/2
 
            
            sint = sin(theta) ; cost = cos(theta) 
            D(1,1) = lambda1*cost*cost + lambda2*sint*sint 
            D(2,1) = lambda1*sint*cost - lambda2*sint*cost
            D(1,2) = D(2,1)
            D(2,2) = lambda1*sint*sint + lambda2*cost*cost 

            return
        end function getDMatrix



        subroutine eigenval2x2(D,lambda1,lambda2)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      perform an eigendecomposition of the symmetric 2x2 matrix 2
    !*      return the eigenvalues 
    !*      with lambda1 <= lambda2
            real(kind=real64),dimension(2,2),intent(in)             ::      D
            real(kind=real64),intent(out)                           ::      lambda1,lambda2

            real(kind=real64)           ::      tron2,det,disc

            tron2 = ( D(1,1) + D(2,2) )/2

            if (abs(D(1,2))<1.0d-12*tron2) then
                if (D(1,1)>D(2,2)) then
                    lambda1 = D(2,2)
                    lambda2 = D(1,1)
                else
                    lambda1 = D(1,1)
                    lambda2 = D(2,2)
                end if
                return
            end if

            det = D(1,1)*D(2,2) - D(1,2)*D(1,2)
            if (abs(det)<1.0d-12) then
                !   can't find eigendecomp
                lambda1 = D(1,1)
                lambda2 = D(2,2)
                return
            end if

            disc = sqrt(max(0.0d0,tron2*tron2 - det))       !   can't be negative, max is to avoid -0.0 error.
            lambda1 = tron2 - disc 
            lambda2 = tron2 + disc 

            return
        end subroutine eigenval2x2        

!-------

        subroutine findBoundingSquare1(D,dx,mult)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the bounding square box required to draw a single ellipse. 
    !*      note that D is in pixels
    !*      dx is half widths of box.
    !*      Assume gaussian profile
            real(kind=real64),dimension(2,2),intent(in)             ::      D
            integer,intent(out)                                     ::      dx 
            real(kind=real64),intent(in),optional                   ::      mult
 
            real(kind=real64)               ::      ss,lambda1,lambda2
            

        !---    find the eigenvalues  of D
            call eigenval2x2(D,lambda1,lambda2)
 
        !       now D = lambda1 x1 x1^T + lambda2 x2 x2^T
        !   the smaller eigenvalue will correspond to the longer axis direction, giving the max radius covered by the ellipse
             
            if (lambda1 > 0) then       
                ss = 1/sqrt(2*lambda1)         !   major radius
                if (present(mult)) then                    
                    dx = ceiling( 1 + mult*ss )
                else
                    dx = ceiling( 1 + LIB_DRAWELLIPSE_GAUSSMULT*ss )
                end if
            else
                dx = 1 
            end if
 
            return
        end subroutine findBoundingSquare1


        subroutine findBoundingBox2(D0,p0,D1,p1,dx,dy,mult)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the bounding rectangular box required to draw two ellipses
    !*      note that D is in pixels
    !*      dx,dy is half widths of box.
    !*      Assume gaussian profile
            real(kind=real64),dimension(2,2),intent(in)             ::      D0,D1
            real(kind=real64),dimension(2),intent(in)               ::      p0,p1           !   fractional pixel offset of centre 
            integer,intent(out)                                     ::      dx,dy
            real(kind=real64),intent(in),optional                   ::      mult

            integer                     ::      dx0,dx1
            real(kind=real64)           ::      mult_

            mult_ = LIB_DRAWELLIPSE_GAUSSMULT ; if (present(mult)) mult_ = LIB_DRAWELLIPSE_GAUSSMULT
 
            call findBoundingSquare1(D0,dx0,mult_)
            call findBoundingSquare1(D1,dx1,mult_)

            dx = ceiling( abs(p1(1) - p0(1)) + (dx0+dx1) )
            dy = ceiling( abs(p1(2) - p0(2)) + (dx0+dx1) )

            return
        end subroutine findBoundingBox2


        subroutine findBoundingRectangle3(img,x0,y0,x1,y1,thresh)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find the bounding rectangular box (x0,y0)-(x1,y1) containing all pixels 
    !*      with intensity >= 0.5d0 ( or thresh if provided )
    !*      on output the columns img(x0,:) and img(x1,:) have at least one pixel intensity >= 0.5
    !*      returns 0 if no pixels above 0.5
            real(kind=real64),dimension(0:,0:),intent(in)       ::      img
            integer,intent(out)                                 ::      x0,y0,x1,y1
            real(kind=real64),intent(in),optional               ::      thresh
            real(kind=real64)   ::      f0
            integer             ::      ii,jj

            x0 = huge(1)
            y0 = huge(1)
            x1 = -huge(1)
            y1 = -huge(1)
            f0 = 0.5d0 ; if (present(thresh)) f0 = thresh

            do jj = 0,size(img,dim=2)-1
                do ii = 0,size(img,dim=1)-1
                    if (img(ii,jj)>=f0) then
                        x0 = min(x0,ii)
                        y0 = min(y0,jj)
                        x1 = max(x1,ii)
                        y1 = max(y1,jj)
                    end if
                end do
            end do

            if (x1<x0) then
                !   didn't find any pixels over f0
                x0 = 0
                y0 = 0
                x1 = -1
                y1 = -1          
            end if
            return
        end subroutine findBoundingRectangle3

!-------

        real(kind=real64) function intersectionOverUnion(D0,p0,D1,p1,mode)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the intersection over union for two ellipses
            real(kind=real64),dimension(2,2),intent(in)             ::      D0,D1
            real(kind=real64),dimension(2),intent(in)               ::      p0,p1           !   fractional pixel offset of centre 
            integer,intent(in)                                      ::      mode

            integer                     ::      dx0,dx1,dy0,dy1,i0,j0,i1,j1
            integer                     ::      nx0,ny0,nx1,ny1
            integer                     ::      x00,x01,y00,y01,x10,x11,y10,y11
            real(kind=real64),dimension(:,:),allocatable            ::      img0,img1,img
            real(kind=real64)           ::      ss0,ss1,lambda1,lambda2,dd,uu,aa,xx


            intersectionOverUnion = 0

        !---    check for easy overlap, circles 
            if (mode == LIB_DRAWELLIPSE_IOU_CIRCLE) then
            !---    find the major radii
                call eigenval2x2(D0,lambda1,lambda2)
                ss0 = 1/sqrt(2*lambda1)
                call eigenval2x2(D1,lambda1,lambda2)
                ss1 = 1/sqrt(2*lambda1)

                !print *,"major radii ",ss0,ss1

            !---    too far?
                dd = norm2( p1-p0 )
                if (dd>=ss0+ss1) return

            !---    total area both circles   
                aa = PI*(ss0*ss0+ss1*ss1)       
                !print *,"area ",aa

                if (dd<1.0d-12) then
                    !   circles overap
                    xx = min(ss0,ss1)
                    uu = PI * xx*xx
                    intersectionOverUnion = uu/(aa - uu)
                    return
                end if

            
            !---    area of intersection https://mathworld.wolfram.com/Circle-CircleIntersection.html
                xx = (dd*dd + ss0*ss0 - ss1*ss1)/(2*dd*ss0) ; xx = max(-1.0d0,min(1.0d0,xx))      
                uu = ss0*ss0*acos( xx )

                xx = (dd*dd + ss1*ss1 - ss0*ss0)/(2*dd*ss1) ; xx = max(-1.0d0,min(1.0d0,xx))      
                uu = uu + ss1*ss1*acos( xx )

                xx = (-dd+ss0+ss1)*(dd+ss0-ss1)*(dd-ss0+ss1)*(dd+ss0+ss1) ; xx = max(0.0d0,xx)
                uu = uu - sqrt(xx)/2
                                
                intersectionOverUnion = uu/(aa - uu)
                return
            end if

                

        !---    check if any possibility of overlap
            
            call findBoundingSquare1(D0,dx0,LIB_DRAWELLIPSE_GAUSSMULT)
            call findBoundingSquare1(D1,dx1,LIB_DRAWELLIPSE_GAUSSMULT) 
            if (abs(p1(1)-p0(1))>dx0+dx1+1) return
            if (abs(p1(2)-p0(2))>dx0+dx1+1) return


        !---    find pixel offset of centre of each ellipse
            i0 = floor(p0(1))   
            j0 = floor(p0(2))   
            i1 = floor(p1(1))   
            j1 = floor(p1(2))   
            
        !---    construct flat ellipse images with narrow (but non-zero) antaliasing
            call drawFlatProfile(D0,p0(1)-i0,p0(2)-j0,img0,mult=1.0d0)
            call drawFlatProfile(D1,p1(1)-i1,p1(2)-j1,img1,mult=1.0d0)
            nx0 = ubound(img0,dim=1)            !   now I know that the first ellipse is contained in box i0-nx0:i0+nx0
            ny0 = ubound(img0,dim=2)
            nx1 = ubound(img1,dim=1)
            ny1 = ubound(img1,dim=2)

 
            if (mode == LIB_DRAWELLIPSE_IOU_RECTANGLE) then
                    
            !---    find the true bounding boxes for the flat ellipses taking account of the antialiasing
                call findBoundingRectangle3(img0,x00,y00,x01,y01)           !   now I know that the first ellipse is contained in box i0-nx0+x00:i0+nx0+x01
                call findBoundingRectangle3(img1,x10,y10,x11,y11)
 
            !---    total area both rectangles   
                dx0 = (i0+x01) - (i0+x00)
                dy0 = (j0+y01) - (j0+y00)
                dx1 = (i1+x11) - (i1+x10)
                dy1 = (j1+y11) - (j1+y10)         
  
                aa = dx0*dy0 + dx1*dy1



            !---    intersection
                if ( (i1-nx1+x10)>=(i0-nx0+x00) ) then
                    dx0 = max(0,(i0-nx0+x01)-(i1-nx1+x10))
                else
                    dx0 = max(0,(i1-nx1+x11)-(i0-nx0+x00))
                end if
        
                if ( (j1-ny1+y10)>=(j0-ny0+y00) ) then
                    dy0 = max(0,(j0-ny0+y01)-(j1-ny1+y10))
                else
                    dy0 = max(0,(j1-ny1+y11)-(j0-ny0+y00))
                end if

                uu = dx0*dy0

                intersectionOverUnion = uu/(aa - uu)
                return
            end if

            if (mode == LIB_DRAWELLIPSE_IOU_ELLIPSE) then

            !---    flatten the ellipse images
                where(img0>=0.5d0)
                    img0 = 1.0d0
                elsewhere
                    img0 = 0.0d0
                endwhere

                where(img1>=0.5d0)
                    img1 = 1.0d0
                elsewhere
                    img1 = 0.0d0
                endwhere

            !---    total area both ellipses
                aa = sum(img0) + sum(img1)

            !---    blit both images onto the same bitmap
                allocate(img( min(i0-nx0,i1-nx1):max(i0+nx0,i1+nx1) , min(j0-ny0,j1-ny1):max(j0+ny0,j1+ny1) ))
                img = 0
                img(i0-nx0:i0+nx0,j0-ny0:j0+ny0) = img0(-nx0:nx0,-ny0:ny0)
                img(i1-nx1:i1+nx1,j1-ny1:j1+ny1) = img(i1-nx1:i1+nx1,j1-ny1:j1+ny1) + img1(-nx1:nx1,-ny1:ny1)
                uu = count( img > 1.9999d0 )
                


                intersectionOverUnion = uu/(aa - uu)
                return
            end if

            return    
        end function intersectionOverUnion

!-------        

        subroutine drawGaussianProfile(D,x0,y0,img ,ignore, dellipse , mult)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that D is in pixels
    !*      The ordering of the optional output derivative data array is 
    !*          [ x0,y0,f0,Dxx,Dxy,Dyy ]                      
    !*      Note that f0 is not strictly needed here, so is not input.
    !*      if ignore, then pixels are set to IGNORE outside gaussian range
            real(kind=real64),dimension(2,2),intent(in)                 ::      D
            real(kind=real64),intent(in)                                ::      x0,y0           !   fractional pixel offset of centre in range (0:1)
            real(kind=real64),dimension(:,:),allocatable,intent(out)    ::      img
            logical,intent(in),optional                                 ::      ignore
            real(kind=real64),intent(in),optional                       ::      mult
            real(kind=real64),dimension(:,:,:),allocatable,intent(out),optional         ::      dellipse        !   (1:6,-nx:nx,-ny:ny)


             
            integer             ::      dx,nx,ny
            real(kind=real64)   ::      cc,xx,yy,mult_
            integer             ::      ii,jj

        !---    find a square box big enough to fit the ellipse
            mult_ = LIB_DRAWELLIPSE_GAUSSMULT ; if (present(mult)) mult_ = mult
            call findBoundingBox(D,dx,mult_)  

        !---    draw the ellipse into a square box. First find box bounds
            nx = -huge(1)
            ny = -huge(1)
            do jj = -dx,dx
                yy = jj - y0
                do ii = -dx,dx
                    xx = ii - x0
                    cc = xx*xx*D(1,1) + 2*xx*yy*D(1,2) + yy*yy*D(2,2)
                    if (cc <= mult_*mult_/2) then
                        nx = max(abs(ii),nx)                    !   illuminated pixel furthest along x-direction 
                        ny = max(abs(jj),ny)                    !   illuminated pixel furthest along y-direction                         
                    end if
                end do
            end do 

            nx = min(10000,max(0,nx))
            ny = min(10000,max(0,ny))
        !---    now draw ellipse            
            allocate(img(-nx:nx,-ny:ny)) 
            img = 0 
            if (present(ignore)) then
                if (ignore) img = LIB_DRAWELLIPSE_IGNORE
            end if

            if (present(dellipse)) then
                allocate(dellipse(6,-nx:nx,-ny:ny))
                dellipse = LIB_DRAWELLIPSE_IGNORE
                if (present(ignore)) then
                    if (.not. ignore) dellipse = 0
                end if
    
                do jj = -ny,ny
                    yy = jj - y0
                    do ii = -nx,ny
                        xx = ii - x0
                        cc = xx*xx*D(1,1) + 2*xx*yy*D(1,2) + yy*yy*D(2,2)
                        if (cc <= mult_*mult_/2) then
                            cc = safeExp(-cc)                    
                            img(ii,jj) = cc
                            dellipse(1,ii,jj) = 2*cc*( D(1,1)*xx + D(1,2)*yy )          !    = d f/d x0 
                            dellipse(2,ii,jj) = 2*cc*( D(1,2)*xx + D(2,2)*yy )          !    = d f/d y0  
                            dellipse(3,ii,jj) = cc                                      !    = d f/d f0 
                            dellipse(4,ii,jj) = - cc*xx*xx                              !    = d f/d Dxx
                            dellipse(5,ii,jj) = - cc*2*xx*yy                            !    = d f/d Dxy
                            dellipse(6,ii,jj) = - cc*yy*yy                              !    = d f/d Dyy
                        end if
                    end do
                end do     
            else            
                do jj = -ny,ny
                    yy = jj - y0
                    do ii = -nx,nx
                        xx = ii - x0
                        cc = xx*xx*D(1,1) + 2*xx*yy*D(1,2) + yy*yy*D(2,2)
                        if (cc <= mult_*mult_/2) img(ii,jj) = safeExp(-cc)                        
                    end do
                end do     
            end if

            return
        end subroutine drawGaussianProfile


        
        subroutine drawRingProfile(D,x0,y0,img,ignore,mult)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that D is in pixels
            real(kind=real64),dimension(2,2),intent(in)                 ::      D
            real(kind=real64),intent(in)                                ::      x0,y0           !   fractional pixel offset of centre in range (0:1)
            real(kind=real64),dimension(:,:),allocatable,intent(out)    ::      img
            real(kind=real64),intent(in),optional                       ::      mult
            logical,intent(in),optional                                 ::      ignore
            real(kind=real64),dimension(:,:),allocatable            ::      img_square
            
            integer             ::      dx,nx,ny
            real(kind=real64)   ::      cc,xx,yy,ss,lambda1,lambda2,mult_
           
            integer             ::      ii,jj

        !---    start with finding the profile c = x.D x
        !       and from this the distance from the ring (c-c0)^2
            

        !---    find a square box big enough to fit the ellipse
            mult_ = LIB_DRAWELLIPSE_RINGWIDTH ; if (present(mult)) mult_ = mult
            call findBoundingBox(D,dx,mult_)

            !print *,"drawRingProfile info - findBoundingBox ",dx
            allocate(img_square(-dx:dx,-dx:dx))

            img_square = 0 
            if(present(ignore)) then
                if (ignore) img_square = LIB_DRAWELLIPSE_IGNORE
            end if

        !---    draw the ellipse ring
            nx = -huge(1)
            ny = -huge(1)
            call eigenval2x2(D,lambda1,lambda2)
            ss = min( abs(lambda1),abs(lambda2) )
            ss = max(1.0d0,1/sqrt(2*ss)) / mult_              !   larger radius over width of ring
           
            do jj = -dx,dx

                yy = jj - y0
                do ii = -dx,dx
                    xx = ii - x0
                    cc = xx*xx*D(1,1) + 2*xx*yy*D(1,2) + yy*yy*D(2,2)
                    cc = ss*(sqrt(2*cc) - 1.0d0)                    
                    if (abs(cc) <= 2.0d0) then
                        img_square(ii,jj) = safeExp( - cc*cc/2  )
                        nx = max(abs(ii),nx)                !   illuminated pixel furthest along x-direction 
                        ny = max(abs(jj),ny)                     !   illuminated pixel furthest along y-direction                         
                    end if
                end do
            end do
            !print *,"drawRingProfile info - image extent ",nx,ny

            nx = max(0,nx)
            ny = max(0,ny)

            allocate(img(-nx:nx,-ny:ny))
            img = 0
            if(present(ignore)) then
                if (ignore) img = LIB_DRAWELLIPSE_IGNORE
            end if


        !---    fill in the ellipse
            do jj = -ny,ny
                do ii = -nx,nx
                    img(ii,jj) = img_square(ii,jj)
                end do
            end do
 
            return
        end subroutine drawRingProfile


        
        subroutine drawFlatProfile(D,x0,y0,img,ignore,mult)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that D is in pixels
            real(kind=real64),dimension(2,2),intent(in)                 ::      D
            real(kind=real64),intent(in)                                ::      x0,y0           !   fractional pixel offset of centre in range (0:1)
            real(kind=real64),dimension(:,:),allocatable,intent(out)    ::      img
            logical,intent(in),optional                                 ::      ignore
            real(kind=real64),intent(in),optional                       ::      mult

            real(kind=real64),dimension(:,:),allocatable            ::      img_square
            
            integer             ::      dx,nx,ny
            real(kind=real64)   ::      cc,xx,yy,ss,lambda1,lambda2,mult_
       
            integer             ::      ii,jj

        !---    start with finding the profile c = x.D x
        !       and from this the distance from the ring (sqrt(c)-c0)^2
            

        !---    find a square box big enough to fit the ellipse
            mult_ = LIB_DRAWELLIPSE_RINGWIDTH ; if (present(mult)) mult_ = mult
            call findBoundingBox(D,dx,mult_)
            !print*,"DBG| dx,D is:",dx,D

            allocate(img_square(-dx:dx,-dx:dx))
            img_square=0
            if(present(ignore)) then
                if (ignore) img_square = LIB_DRAWELLIPSE_IGNORE
            end if


        !---    draw the (top half) of the ellipse
            nx = -huge(1)
            ny = -huge(1)
            call eigenval2x2(D,lambda1,lambda2)
            ss = min( abs(lambda1),abs(lambda2) )
            ss = max(1.0d0,1/sqrt(2*ss)) / mult_              !   larger radius over width of ring
            do jj = -dx,dx

                yy = jj - y0
                do ii = -dx,dx
                    xx = ii - x0
                    cc = xx*xx*D(1,1) + 2*xx*yy*D(1,2) + yy*yy*D(2,2)
                    cc = ss*(sqrt(2*cc) - 1.0d0)                    
                    if (cc < 0) then
                        img_square(ii,jj) = 1.0d0
                    else if (cc <= 2.0d0) then
                        img_square(ii,jj) = safeExp( - cc*cc/2 )
                        nx = max(abs(ii),nx)                !   illuminated pixel furthest along x-direction 
                        ny = max(abs(jj),ny)                     !   illuminated pixel furthest along y-direction                         
                    end if
                end do
            end do
            nx = max(0,nx)
            ny = max(0,ny)
            !print *,"drawFlatProfile info - image extent ",nx,ny
            allocate(img(-nx:nx,-ny:ny))
            img = 0
            if(present(ignore)) then
                if (ignore) img = LIB_DRAWELLIPSE_IGNORE
            end if


            
        !---    fill in the ellipse
            do jj = -ny,ny
                do ii = -nx,nx
                    img(ii,jj) = img_square(ii,jj)
                end do
            end do
 
            return
        end subroutine drawFlatProfile


       
        real(kind=real64) function getRss0( ellipse,img_in,mode,mult )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the residual sum of squares
    !*          rss = sum_i ( g_i - f_i )^2
    !*      between an input image "ground truth" g_i = img(x,y) and its representation as a 2d gaussian
    !*      defined by
    !*          f(x,y) = f0 safeExp[ - (x-x0,y-y0).D (x-x0,y-y0) ]
    !*      Some pixels will be ignored, marked up with the special code LIB_RSS_IGNORE
    !*      The ordering of the input data array and output derivatives 
    !*          [ x0,y0,f0,Dxx,Dxy,Dyy ]            

          
            real(kind=real64),dimension(6),intent(in)               ::      ellipse
            real(kind=real64),dimension(0:,0:),intent(in)           ::      img_in
            integer,intent(in)                                      ::      mode
            real(kind=real64),intent(in),optional                   ::      mult
                        
            integer             ::      Nx,Ny
            integer             ::      i0,j0               !   pixel at centre of ellipse
            integer             ::      ix,iy,dx,npx
            real(kind=real64)   ::      ss,xx,yy,gg,mult_,bg
            real(kind=real64)   ::      Dxx,Dxy,Dyy
             
        
        !---    allocate an image for the gaussian
            Nx = size(img_in,dim=1)
            Ny = size(img_in,dim=2)


        !---    find the pixel centre of the ellipse
            i0 = floor(ellipse(1))
            j0 = floor(ellipse(2))

            getRss0 = 0.0d0

        !---    construct the image that will be added
            select case(mode)
                case (LIB_DRAWELLIPSE_SHADE_GAUSSIAN)                    
                    mult_ = LIB_DRAWELLIPSE_GAUSSMULT; if (present(mult)) mult_ = mult
                    Dxx = ellipse(4)
                    Dxy = ellipse(5)
                    Dyy = ellipse(6)

                !---    compute background level
                    bg = 0
                    if (LIB_DRAWELLIPSE_RSS_REMOVE_BG) then
                        npx = 0
                        do iy = 0,Ny-1
                            yy = iy - ellipse(2)
                            do ix = 0,Nx-1                    
                                if ( img_in(ix,iy) == LIB_DRAWELLIPSE_IGNORE ) cycle
                                xx = ix - ellipse(1)
                                gg = xx*xx*Dxx + 2*xx*yy*Dxy + yy*yy*Dyy
                                if ( gg < -10.0d0 ) cycle   !   gg should be positive, so this must mean a very bad gaussian
                                gg = safeExp( - gg )
                                ss = ellipse(3) * gg - img_in(ix,iy) 
                                bg = bg + ss
                                npx = npx + 1
                            end do
                        end do
                        bg = bg / max(1,npx)                    
                    end if

                    do iy = 0,Ny-1
                        yy = iy - ellipse(2)
                        do ix = 0,Nx-1
                            if ( img_in(ix,iy) == LIB_DRAWELLIPSE_IGNORE ) cycle
                            xx = ix - ellipse(1)
                            gg = xx*xx*Dxx + 2*xx*yy*Dxy + yy*yy*Dyy
                            gg = bg + ellipse(3) * safeExp( - gg )
                            if ( gg < -10.0d0 ) cycle   !   gg should be positive, so this must mean a very bad gaussian
                            ss = img_in(ix,iy) - gg
                            getRss0 = getRss0 + ss*ss
                        end do
                    end do                 
                case (LIB_DRAWELLIPSE_SHADE_RING)
                    stop "Lib_DrawEllipse::getRss0 error - not coded for ring"
                case (LIB_DRAWELLIPSE_SHADE_FLAT)
                    stop "Lib_DrawEllipse::getRss0 error - not coded for flat"                    
            end select
 
            return
        end function getRss0


        real(kind=real64) function getRss1( ellipse,img_in,mode,mult )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the residual sum of squares
    !*          rss = sum_i ( g_i - f_i )^2
    !*      between an input image "ground truth" g_i = img(x,y) and its representation as a sum over 2d gaussians
    !*      note: this is the multi-ellipse version of getRss0

          
            real(kind=real64),dimension(:,:),intent(in)             ::      ellipse
            real(kind=real64),dimension(0:,0:),intent(in)           ::      img_in
            integer,intent(in)                                      ::      mode
            real(kind=real64),intent(in),optional                   ::      mult

            integer             ::      Nx,Ny
            integer             ::      i0,j0               !   pixel at centre of ellipse
            integer             ::      ix,iy,dx,ii,npx
            real(kind=real64)   ::       ss,xx,yy,gg,mult_,bg,ff,zz
            real(kind=real64)   ::      Dxx,Dxy,Dyy
            real(kind=real64),dimension(:,:),allocatable        ::      img_ell
        
        !---    allocate an image for the gaussian
            Nx = size(img_in,dim=1)
            Ny = size(img_in,dim=2)

            getRss1 = 0.0d0

        !---    construct the image that will be added
            select case(mode)
                case (LIB_DRAWELLIPSE_SHADE_GAUSSIAN)                    
                    mult_ = LIB_DRAWELLIPSE_GAUSSMULT; if (present(mult)) mult_ = mult
                    allocate(img_ell(0:Nx-1,0:Ny-1))
                    

                !---    find the image difference
                    img_ell = img_in
                    bg = 0.0d0      !   background level: bg = sum_i (img_i - (sum g)_i )  then will be divided by n
                    npx = 0                   
                    ss = 0          !   ss = sum_i (img_i - (sum g)_i )^2
                    do iy = 0,Ny-1                        
                        do ix = 0,Nx-1
                            zz = img_ell(ix,iy)     !       will compute zz = img_i - (sum g)_i
                            if ( zz == LIB_DRAWELLIPSE_IGNORE ) cycle

                            do ii = 1,size(ellipse,dim=2)                    
                                xx  = ix - ellipse(1,ii)
                                yy  = iy - ellipse(2,ii)
                                ff  = ellipse(3,ii)
                                Dxx = ellipse(4,ii)
                                Dxy = ellipse(5,ii)
                                Dyy = ellipse(6,ii)
                                gg = xx*xx*Dxx + 2*xx*yy*Dxy + yy*yy*Dyy
                                !if ( gg < -10.0d0 ) cycle   !   gg should be positive, so this must mean a very bad gaussian                                
                                zz = zz - ff * safeExp( - gg )
                            end do
                            img_ell(ix,iy) = zz             !   put back 
                            bg = bg + zz          
                            npx = npx + 1
                            ss = ss + zz*zz
                        end do                 

                    end do
                    if (LIB_DRAWELLIPSE_RSS_REMOVE_BG) then
                        bg = bg / max(1,npx)                !   bg = sum_i (img_i - (sum g)_i ) / n
                        getRss1 = ss - npx*bg*bg            !   = sum_i (img_i - (sum g)_i - bg)^2 
                                                            !   = sum_i (img_i - (sum g)_i )^2 - 2 bg sum_i (img_i - (sum g)_i )  + sum_i bg^2
                                                            !   = sum_i (img_i - (sum g)_i )^2 - n bg^2 
                    else
                        getRss1 = ss
                    end if

                case (LIB_DRAWELLIPSE_SHADE_RING)
                    stop "Lib_DrawEllipse::getRss0 error - not coded for ring"
                case (LIB_DRAWELLIPSE_SHADE_FLAT)
                    stop "Lib_DrawEllipse::getRss0 error - not coded for flat"                    
            end select
 
            
            return
        end function getRss1


        subroutine findDrss0( ellipse,img_in,mode,rss,drss,mult )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the residual sum of squares
    !*          rss = sum_i ( g_i - f_i )^2
    !*      between an input image "ground truth" g_i = img(x,y) and its representation as a sum over 2d gaussians
    !*      defined by f_i = f(x,y) = 
    !*          f(x,y) = f0 safeExp[ - (x-x0,y-y0).D (x-x0,y-y0) ]
    !*      Some pixels will be ignored, marked up with the special code LIB_RSS_IGNORE
    !*      The ordering of the input data array and output derivatives 
    !*          [ x0,y0,f0,Dxx,Dxy,Dyy ]            
            real(kind=real64),dimension(6),intent(in)               ::      ellipse
            real(kind=real64),dimension(0:,0:),intent(in)           ::      img_in
            integer,intent(in)                                      ::      mode
            real(kind=real64),intent(out)                           ::      rss
            real(kind=real64),dimension(6),intent(out)              ::      drss            
            real(kind=real64),intent(in),optional                   ::      mult
                        
            integer             ::      Nx,Ny
            !integer             ::      i0,j0               !   pixel at centre of ellipse
            integer             ::      ix,iy,dx,npx
            real(kind=real64)   ::      ss,xx,yy,gg,mult_,bg,ff
            real(kind=real64)   ::      Dxx,Dxy,Dyy
             
        
        !---    allocate an image for the gaussian
            Nx = size(img_in,dim=1)
            Ny = size(img_in,dim=2)

 

            rss = 0.0d0
            drss = 0.0d0
            !f0=ellipse(3)

        !---    construct the image that will be added
            select case(mode)
                case (LIB_DRAWELLIPSE_SHADE_GAUSSIAN)                    
                    mult_ = LIB_DRAWELLIPSE_GAUSSMULT; if (present(mult)) mult_ = mult
                    Dxx = ellipse(4)
                    Dxy = ellipse(5)
                    Dyy = ellipse(6)

                 !---    compute background level
                     bg = 0
                     if (LIB_DRAWELLIPSE_RSS_REMOVE_BG) then
                        npx = 0
                        do iy = 0,Ny-1
                            yy = iy - ellipse(2)
                            do ix = 0,Nx-1                    
                                if ( img_in(ix,iy) == LIB_DRAWELLIPSE_IGNORE ) cycle
                                xx = ix - ellipse(1)
                                gg = xx*xx*Dxx + 2*xx*yy*Dxy + yy*yy*Dyy
                                if ( gg < -10.0d0 ) cycle   !   gg should be positive, so this must mean a very bad gaussian
                                gg = safeExp( - gg )
                                ss = ellipse(3) * gg - img_in(ix,iy) 
                                bg = bg + ss
                                npx = npx + 1
                            end do
                        end do
                        bg = bg / max(1,npx)
                    end if
                     
                    
                    do iy = 0,Ny-1
                        yy = iy - ellipse(2)
                        do ix = 0,Nx-1                    
                            if ( img_in(ix,iy) == LIB_DRAWELLIPSE_IGNORE ) cycle
                            xx = ix - ellipse(1)
                            gg = xx*xx*Dxx + 2*xx*yy*Dxy + yy*yy*Dyy                            
                            if ( gg < -10.0d0 ) cycle   !   gg should be positive, so this must mean a very bad gaussian

                            gg = safeExp( - gg )
                            
                            ff = ellipse(3)
                            ss = bg + ff * gg - img_in(ix,iy) 
                            rss = rss + ss*ss

                            drss(1) = drss(1) + 2*ss*( 2*gg*ff*( Dxx*xx + Dxy*yy ) )         !    = d ss^2 /d x0 
                            drss(2) = drss(2) + 2*ss*( 2*gg*ff*( Dxy*xx + Dyy*yy ) )         !    = d ss^2 /d y0  
                            drss(3) = drss(3) + 2*ss*( gg                          )         !    = d ss^2 /d f0 
                            drss(4) = drss(4) + 2*ss*( -   gg*ff*xx*xx             )         !    = d ss^2 /d Dxx
                            drss(5) = drss(5) + 2*ss*( - 2*gg*ff*xx*yy             )         !    = d ss^2 /d Dxy
                            drss(6) = drss(6) + 2*ss*( -   gg*ff*yy*yy             )         !    = d ss^2 /d Dyy

                            
                        end do
                    end do                 
                case (LIB_DRAWELLIPSE_SHADE_RING)
                    stop "Lib_DrawEllipse::getRss0 error - not coded for ring"
                case (LIB_DRAWELLIPSE_SHADE_FLAT)
                    stop "Lib_DrawEllipse::getRss0 error - not coded for flat"                    
            end select
 
            
            return
        end subroutine findDrss0

        subroutine findDrss1( ellipse,img_in,mode,rss,drss,mult )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the residual sum of squares
    !*          rss = sum_i ( g_i - f_i )^2
    !*      between an input image "ground truth" g_i = img(x,y) and its representation as a sum over 2d gaussians
    !*      defined by f_i = f(x,y) = 
    !*          f(x,y) = bg + f0 safeExp[ - (x-x0,y-y0).D (x-x0,y-y0) ]
    !*      Some pixels will be ignored, marked up with the special code LIB_RSS_IGNORE
    !*      The ordering of the input data array and output derivatives per gaussian is
    !*          [ x0,y0,f0,Dxx,Dxy,Dyy ]            
    !*      computes the background level and returns
            real(kind=real64),dimension(:,:),intent(in)             ::      ellipse     !   (6,nn)
            real(kind=real64),dimension(0:,0:),intent(in)           ::      img_in
            integer,intent(in)                                      ::      mode
            real(kind=real64),intent(out)                           ::      rss
            real(kind=real64),dimension(:,:),intent(out)            ::      drss        !   (6,nn)
            real(kind=real64),intent(in),optional                   ::      mult
                         
            
            integer             ::      Nx,Ny , nn,npx
            !integer             ::      i0,j0               !   pixel at centre of ellipse
            integer             ::      ix,iy,dx,ii
            real(kind=real64)   ::      ss,xx,yy,gg,mult_,ff,bg
            real(kind=real64)   ::      Dxx,Dxy,Dyy
            real(kind=real64),dimension(:,:,:),allocatable          ::      img_ell
        
        !---    allocate an image for the gaussian
            Nx = size(img_in,dim=1)
            Ny = size(img_in,dim=2)
            nn = size(ellipse,dim=2)

            rss = 0.0d0
            drss = 0.0


        !---    construct the image that will be added
            select case(mode)
                case (LIB_DRAWELLIPSE_SHADE_GAUSSIAN)    

                    mult_ = LIB_DRAWELLIPSE_GAUSSMULT; if (present(mult)) mult_ = mult
                    allocate(img_ell(0:nn,0:Nx-1,0:Ny-1))

 
                !---    compute difference between image and sum(gaussian)                   
                    img_ell(0,:,:) = -img_in 
                    img_ell(1:,:,:) = 0
                    do ii = 1,nn

                        Dxx = ellipse(4,ii)
                        Dxy = ellipse(5,ii)
                        Dyy = ellipse(6,ii)
 
                        do iy = 0,Ny-1
                            yy = iy - ellipse(2,ii)
                            do ix = 0,Nx-1                        
                                if ( img_in(ix,iy) == LIB_DRAWELLIPSE_IGNORE ) cycle
                                xx = ix - ellipse(1,ii)
                                gg = xx*xx*Dxx + 2*xx*yy*Dxy + yy*yy*Dyy
                                if ( gg < -10.0d0 ) cycle   !   gg should be positive, so this must mean a very bad gaussian
                                gg = safeExp( - gg )
                                img_ell(0,ix,iy) = img_ell(0,ix,iy) + ellipse(3,ii) * gg 
                                img_ell(ii,ix,iy) = gg
                            end do
                        end do  

                    end do

                !---    find the background level : bg = < img - sum(gaussian) > . 
                    bg = 0.0d0
                    if (LIB_DRAWELLIPSE_RSS_REMOVE_BG) then
                        npx = 0
                        do iy = 0,Ny-1
                            do ix = 0,Nx-1                        
                                if ( img_in(ix,iy) /= LIB_DRAWELLIPSE_IGNORE ) then
                                    bg = bg + img_ell(0,ix,iy)              
                                    npx = npx + 1
                                end if        
                            end do
                        end do        
                        bg = bg/max(1,npx)
                    end if

                !---    now compute rss and drss
                    do iy = 0,Ny-1
                        do ix = 0,Nx-1
                            if ( img_in(ix,iy) == LIB_DRAWELLIPSE_IGNORE ) cycle
                            ss = img_ell(0,ix,iy) - bg
                            rss = rss + ss*ss
                            do ii = 1,nn
                                xx = ix - ellipse(1,ii)
                                yy = iy - ellipse(2,ii)
                                gg = img_ell(ii,ix,iy)
                                ff = ellipse(3,ii)                            
                                Dxx = ellipse(4,ii)
                                Dxy = ellipse(5,ii)
                                Dyy = ellipse(6,ii)

                                drss(1,ii) = drss(1,ii) + 2*ss*( 2*ff * gg *( Dxx*xx + Dxy*yy ) )        !    = d f/d x  
                                drss(2,ii) = drss(2,ii) + 2*ss*( 2*ff * gg *( Dxy*xx + Dyy*yy ) )                    
                                drss(3,ii) = drss(3,ii) + 2*ss*(        gg                      )        !    = d f/d f0 
                                drss(4,ii) = drss(4,ii) + 2*ss*(  -ff * gg *xx*xx               )        !    = d f/d Dxx
                                drss(5,ii) = drss(5,ii) + 2*ss*(-2*ff * gg *xx*yy               )                    
                                drss(6,ii) = drss(6,ii) + 2*ss*(  -ff * gg *yy*yy               )
                            end do
                        end do
                    end do             

                case (LIB_DRAWELLIPSE_SHADE_RING)
                    stop "Lib_DrawEllipse::getRss0 error - not coded for ring"
                case (LIB_DRAWELLIPSE_SHADE_FLAT)
                    stop "Lib_DrawEllipse::getRss0 error - not coded for flat"                    
            end select

            return
        end subroutine findDrss1

                            
        subroutine add0(D,p,img,mode,f,mult)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      note that D is in pixels
    !*      mode = LIB_DRAWELLIPSE_SHADE_GAUSSIAN etc describes the drawing mode to use 
    !*      this routine directly adds intensity onto img            
    !*      see also: drawEllipse which blends intensity levels and offers more functionality
            real(kind=real64),dimension(2,2),intent(in)                 ::      D
            real(kind=real64),dimension(2),intent(in)                   ::      p                           !   fractional pixel offset of centre  
            real(kind=real64),dimension(0:,0:),intent(inout)            ::      img
            integer,intent(in)                                          ::      mode
            real(kind=real64),intent(in)                                ::      f                           !   intensity peak (needed for absolute image intensity and for deriv)
            real(kind=real64),intent(in),optional                       ::      mult
            real(kind=real64),dimension(:,:),allocatable                ::      img_blit
            real(kind=real64)   ::      xx
            real(kind=real64)                       ::      mult_
            integer             ::      i0,j0           !   pixel at centre of ellipse
            integer             ::      nx,ny           !   half width of ellipse image
            integer             ::      ii,jj
            integer             ::      Mx,My           !   size of img


        !---    find the pixel centre of the ellipse
            i0 = floor(p(1))
            j0 = floor(p(2))

        !---    construct the image that will be added
            select case(mode)
                case (LIB_DRAWELLIPSE_SHADE_RING)                    
                    mult_ = LIB_DRAWELLIPSE_RINGWIDTH; if (present(mult)) mult_ = mult
                    call drawRingProfile(    D,p(1)-i0,p(2)-j0,img_blit,mult=mult_ )
                case (LIB_DRAWELLIPSE_SHADE_FLAT)
                    mult_ = LIB_DRAWELLIPSE_GAUSSMULT; if (present(mult)) mult_ = mult
                    call drawFlatProfile(    D,p(1)-i0,p(2)-j0,img_blit,mult=mult_)
                case default
                    mult_ = LIB_DRAWELLIPSE_GAUSSMULT; if (present(mult)) mult_ = mult
                    call drawGaussianProfile(D,p(1)-i0,p(2)-j0,img_blit,mult=mult_)
            end select

            nx = ubound(img_blit,dim=1)
            ny = ubound(img_blit,dim=2)
            Mx = size(img,dim=1)
            My = size(img,dim=2)

        !---    add the image
            do jj = max(0,j0-ny),min(My-1,j0+ny)
                do ii = max(0,i0-nx),min(Mx-1,i0+nx)
                    xx = img_blit( ii-i0 , jj-j0 )
                    img(ii,jj) = img(ii,jj) + f*xx 
                end do
            end do

            return
        end subroutine add0       

        subroutine add1( ellipse,img,mode,mult )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      mode = LIB_DRAWELLIPSE_SHADE_GAUSSIAN etc describes the drawing mode to use 
    !*      this routine directly adds intensity onto img            
    !*      see also: drawEllipse which blends intensity levels and offers more functionality
    !*      The ordering of the input data array 
    !*          [ x0,y0,f0,Dxx,Dxy,Dyy ]                       
            real(kind=real64),dimension(6),intent(in)                   ::      ellipse
            real(kind=real64),dimension(0:,0:),intent(inout)            ::      img
            integer,intent(in)                                          ::      mode
            real(kind=real64),intent(in),optional                       ::      mult
            real(kind=real64)                       ::      mult_
            
        !---     
            mult_ = LIB_DRAWELLIPSE_GAUSSMULT; if (present(mult)) mult_ = mult

            call add0( reshape( (/ellipse(4),ellipse(5),ellipse(5),ellipse(6)/),(/2,2/) )                &
                            ,(/ ellipse(1),ellipse(2) /)                                                 &
                            ,img , mode, f=ellipse(3)  )


            return
        end subroutine add1   
        
        


        real(kind=real64) function getT0( ellipse,mode,sigma,mult )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the input density function, return the t-value
    !*          t = |<f>| / (sigma/sqrt(n))
    !*          f(x,y) = f0 safeExp[ - (x-x0,y-y0).D (x-x0,y-y0) ]    
    !*      The ordering of the input data array  
    !*          [ x0,y0,f0,Dxx,Dxy,Dyy ]            

          
            real(kind=real64),dimension(6),intent(in)               ::      ellipse
            integer,intent(in)                                      ::      mode
            real(kind=real64),intent(in)                            ::      sigma       !   background noise std dev.
            real(kind=real64),intent(in),optional                   ::      mult
                                    
            real(kind=real64),dimension(:,:),allocatable                ::      img_blit
            real(kind=real64)                       ::      weight , mult_
            integer                                 ::      nn
            integer                                 ::      i0,j0           !   pixel at centre of ellipse
            integer                                 ::      nx,ny           !   half width of ellipse image
            integer                                 ::      ii,jj
            real(kind=real64),dimension(2,2)        ::      DD
            
        !---    for a real T-value should be over 1 std dev, but can try longer range with sigma=1 to get avg intensity, say.
            mult_ = 1.0d0; if (present(mult)) mult_ = mult


        !---    find the pixel centre of the ellipse
            i0 = floor(ellipse(1))
            j0 = floor(ellipse(2))

        !---    construct the image that will be added
            DD = reshape( (/ellipse(4),ellipse(5),ellipse(5),ellipse(6)/),(/2,2/) )
            select case(mode)
                case (LIB_DRAWELLIPSE_SHADE_RING)
                    call drawRingProfile( DD ,ellipse(1)-i0,ellipse(2)-j0,img_blit,ignore = .true. , mult = mult_)      !   note mult=1, we only look out to 1 std dev range.
                case (LIB_DRAWELLIPSE_SHADE_FLAT)
                    call drawFlatProfile(DD,ellipse(1)-i0,ellipse(2)-j0,img_blit,ignore = .true. , mult = mult_)
                case default
                    call drawGaussianProfile(DD,ellipse(1)-i0,ellipse(2)-j0,img_blit,ignore = .true. , mult = mult_)
            end select

            nx = ubound(img_blit,dim=1)
            ny = ubound(img_blit,dim=2)

            weight = 0                !   calculates integral intensity. I'm not doing this analytically as I might have a ring or something 'orrible.
            do jj = -ny,ny
                do ii = -nx,nx
                    if ( img_blit(ii,jj) == LIB_DRAWELLIPSE_IGNORE ) cycle
                    weight = weight + img_blit(ii,jj)
                end do
            end do
            weight = weight * ellipse(3)

            nn = (2*nx+1)*(2*ny+1) - count( img_blit == LIB_DRAWELLIPSE_IGNORE )        !   sums lit pixels
            getT0 = 0
            if (nn*sigma>0) getT0 = abs(weight) / sqrt( nn*sigma*sigma )
              
            return
        end function getT0        


        real(kind=real64) function getWeight0( ellipse,mode,mult )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the input density function, return the integral
    !*          f(x,y) = f0 safeExp[ - (x-x0,y-y0).D (x-x0,y-y0) ]    
    !*      The ordering of the input data array  
    !*          [ x0,y0,f0,Dxx,Dxy,Dyy ]            

          
            real(kind=real64),dimension(6),intent(in)               ::      ellipse
            integer,intent(in)                                      ::      mode
            real(kind=real64),intent(in),optional                   ::      mult
                                    
            real(kind=real64),dimension(:,:),allocatable                ::      img_blit
            real(kind=real64)                       ::      mult_            
            integer                                 ::      i0,j0           !   pixel at centre of ellipse
            integer                                 ::      nx,ny           !   half width of ellipse image
            integer                                 ::      ii,jj
            real(kind=real64),dimension(2,2)        ::      DD
            
         
            mult_ = LIB_DRAWELLIPSE_GAUSSMULT; if (present(mult)) mult_ = mult


        !---    find the pixel centre of the ellipse
            i0 = floor(ellipse(1))
            j0 = floor(ellipse(2))

        !---    construct the image that will be added
            DD = reshape( (/ellipse(4),ellipse(5),ellipse(5),ellipse(6)/),(/2,2/) )
            select case(mode)
                case (LIB_DRAWELLIPSE_SHADE_RING)
                    call drawRingProfile( DD ,ellipse(1)-i0,ellipse(2)-j0,img_blit, mult = mult_ )      !   note mult=1, we only look out to 1 std dev range.
                case (LIB_DRAWELLIPSE_SHADE_FLAT)
                    call drawFlatProfile(DD,ellipse(1)-i0,ellipse(2)-j0,img_blit,mult = mult_ )
                case default
                    call drawGaussianProfile(DD,ellipse(1)-i0,ellipse(2)-j0,img_blit, mult = mult_ )
            end select

            nx = ubound(img_blit,dim=1)
            ny = ubound(img_blit,dim=2)

        !---    find integral            
            getWeight0 = 0                !   calculates integral intensity. I'm not doing this analytically as I might have a ring or something 'orrible.
            do jj = -ny,ny
                do ii = -nx,nx
                    getWeight0 = getWeight0 + img_blit(ii,jj)
                end do
            end do
            getWeight0 = getWeight0 * ellipse(3)
  
            return
        end function getWeight0
        
        subroutine do_spots_overlap(spot1_data,spot2_data,iou_in,noof_sigma,do_they_overlap)
            !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !*  takes in two arrays of spot data (x,y,dmaj,dmin,f,theta) converts them into covariance matrices,
            !*  checks for ciruclar overlap, then if over threshold chekcs for numerical ellipse overlap both 
            !*  using IoU threshold.
                    real(kind=real64), intent(in), dimension(6)         ::      spot1_data, spot2_data
                    real(kind=real64), intent(in)                       ::      iou_in
                    real(kind=real64), intent(in)                       ::      noof_sigma !how many stnadard deviation to define spot radii
                    logical,intent(out)                                 ::      do_they_overlap
                    real(kind=real64),dimension(2,2)                    ::      d1,d2 !coveriance matrices
                    real(kind=real64)                                   ::      iou_circle,iou_ellipse
        
                    !get D matrices
                    d1=getDMatrix( spot1_data(3), spot1_data(5), spot1_data(6), noof_sigma )
                    d2=getDMatrix( spot2_data(3), spot2_data(5), spot2_data(6), noof_sigma )
        
                    !circular overlap check
                    !intersectionOverUnion(D0,p0,D1,p1,mode)
                    iou_circle=intersectionOverUnion(d1,spot1_data(1:2),d2,spot2_data(1:2),LIB_DRAWELLIPSE_IOU_CIRCLE)
        
                    if(iou_circle>=iou_in)then
                        !ellitical ovelap check
                        iou_ellipse=intersectionOverUnion(d1,spot1_data(1:2),d2,spot2_data(1:2),LIB_DRAWELLIPSE_IOU_ELLIPSE)
                        if(iou_ellipse>=iou_in) then
                            do_they_overlap=.true.
                            return
                        end if
                    else
                        do_they_overlap=.false.
                    end if
                end subroutine do_spots_overlap
        
        subroutine check_overlap_with_prexisting_spots(spot_data_in,existing_spots_array_in,iou_tol,noof_sigma,do_any_overlap)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*  takes in an array of spot data (x,y,dmaj,dmin,f,theta) and a array containing,
    !*  data of multiple spots and checks if the new spot ovealps with any exsitng one  
    !*  using IoU threshold.
            real(kind=real64), intent(in), dimension(6)         ::      spot_data_in
            real(kind=real64), intent(in), dimension(:,:)       ::      existing_spots_array_in
            real(kind=real64), intent(in)                       ::      iou_tol
            real(kind=real64), intent(in)                       ::      noof_sigma !how many stnadard deviation to define spot radii
            logical,intent(out)                                 ::      do_any_overlap
            integer                                             ::      noof_spots,ii

            noof_spots = size(existing_spots_array_in,dim=2)
            do_any_overlap=.false.

            do ii=1,noof_spots
                call do_spots_overlap(spot_data_in,existing_spots_array_in(:,ii),iou_tol,noof_sigma,do_any_overlap)
                if (do_any_overlap) exit
            end do

        end subroutine check_overlap_with_prexisting_spots
 
    end module Lib_DrawEllipse