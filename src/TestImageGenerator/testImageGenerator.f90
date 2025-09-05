    
    program testImageGenerator
!---^^^^^^^^^^^^^^^^^^^^^^^^^^
!*
!*      A simple program to generate some noisy test images
!*
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    testImageGenerator from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
!*    Copyright (C) 2022  Daniel Mason

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
!*
!*      
        use iso_fortran_env
        use Lib_CommandLineArguments
        use Lib_RandomSeed
        use Lib_Perlin2dPlane
        use Lib_LogNormal
        use Lib_Png
        use Lib_RidlerCalvard
        use Lib_DrawEllipse
        implicit none
        
        
        character(len=8)            ::      VERSION = "0.0.1"
        
    !---    properties of the image determined by the command line arguments
        type(CommandLineArguments)      ::      cla
        integer                     ::          Nx = 1024,Ny = 1024                     !   image size (px)
        integer                     ::          seed = 12345                            !   random seed
        
        real(kind=real64)           ::          perlin_lengthscale = 64.0d0
        real(kind=real64)           ::          perlin_avg = 0.5d0
        real(kind=real64)           ::          perlin_height = 0.1d0
        logical                     ::          usePerlin = .false.
        
        real(kind=real64)           ::          noise_stddev = 0.0
        real(kind=real64)           ::          salt = 0.0d0
        real(kind=real64)           ::          pepper = 0.0d0
        
        integer                     ::          spotCount = 0
        integer                     ::          max_noof_spot_placements=100
        real(kind=real64)           ::          spotLogNormalMu = 10
        real(kind=real64)           ::          spotLogNormalSig = 1
        real(kind=real64)           ::          spotHeightMu  = 0.1d0
        real(kind=real64)           ::          spotHeightSig = 0.1d0
        real(kind=real64)           ::          background_ini_value = 0.0d0
        real(kind=real64)           ::          sin_wave_amplitude = 0.0d0
        real(kind=real64)           ::          max_overlap = 1.0d0
        
        logical                     ::          circles = .false.
        logical                     ::          ellipses = .false.
        logical                     ::          test_dbg_centre = .false. 

        real(kind=real64)           ::          max_eccentricity=1.0d0
        real(kind=real64)           ::          orientation=0.0d0
        real(kind=real64),dimension(2,2)    ::  dmatrix
        real(kind=real64),dimension(2)      ::  position
        
        
        character(len=256)          ::          filename = "test.png"
        character(len=256)          ::          filename_spots = "test.spots"
        real(kind=real64),parameter         ::      TWO_PI = 3.141592653590d0*2.0d0

        real(kind=real64),dimension(:,:),allocatable    ::   spot_array !x,y,diameter,intensity
        
    !---    the image itself        
        real(kind=real64),dimension(:,:),allocatable              ::      img, perlin_img
        
        type(LogNormal)             ::          logn
        real(kind=real64)           ::      bb,tt,ff,bstd,fbar,snr,fot !for Ridler Calvard
        
        
    !---    dummy        
        integer                     ::      ii,jj , nn,ll
        integer                     ::      ix,iy,kk,noof_spot_placement
        real(kind=real64)           ::      zeta,dd,rr,i2s2,xx
        logical                     ::      do_any_overlap
        real(kind=real64),dimension(6)           :: temp_spot_data
        
        
        
        
        
    !---    check command line arguments
         
        cla = CommandLineArguments_ctor(100)
        call setProgramDescription( cla, "testImageGenerator.exe" )
        call setProgramVersion( cla, VERSION )
        
        
        call get( cla,"f",filename ,LIB_CLA_OPTIONAL,"      output filename for image" )
        call get( cla,"o",filename_spots ,LIB_CLA_OPTIONAL,"      output filename for .spots file" )
        
        call get( cla,"Nx",Nx ,LIB_CLA_OPTIONAL,"     width of image" )
        Ny = Nx
        call get( cla,"Ny",Ny ,LIB_CLA_OPTIONAL,"     height of image" )
        call get( cla,"seed",seed ,LIB_CLA_OPTIONAL,"   random seed ( set to 0 to use system clock )" )  
        
        call get( cla,"length",perlin_lengthscale ,LIB_CLA_OPTIONAL," perlin procedural noise lengthscale, at least one perlin parameter must be present to activate" )  
        call get( cla,"avg",perlin_avg ,LIB_CLA_OPTIONAL,"    perlin procedural noise average level" )  
        call get( cla,"height",perlin_height ,LIB_CLA_OPTIONAL," perlin procedural noise half height" )  
        usePerlin = hasArgument(cla,"length").or. hasArgument(cla,"avg").or.  hasArgument(cla,"height")
        
        
        call get( cla,"count",spotCount ,LIB_CLA_OPTIONAL,"  count of spots" )  
        call get( cla,"Rmu",spotLogNormalMu ,LIB_CLA_OPTIONAL,"     spot radius distribution mu" )  
        call get( cla,"Rsig",spotLogNormalSig ,LIB_CLA_OPTIONAL,"    spot radius distribution sigma" )  
        call get( cla,"Imu",spotHeightMu,LIB_CLA_OPTIONAL,"    spot height above bg mean " )  
        call get( cla,"Isig",spotHeightSig,LIB_CLA_OPTIONAL,"   spot height above bg stddev" )  
        call get( cla,"c",circles,LIB_CLA_OPTIONAL,"      draw circles not gaussians" )  
        call get( cla,"el",ellipses,LIB_CLA_OPTIONAL,"      draw elliptical not circular gaussians" ) 
        call get( cla,"emax",max_eccentricity,LIB_CLA_OPTIONAL,"      If drawing elliptical gaussians, what max eccentricity allowed?" ) 
        
        call get( cla,"stddev",noise_stddev ,LIB_CLA_OPTIONAL," noise standard deviation" )  
        call get( cla,"salt",salt ,LIB_CLA_OPTIONAL,"   fraction of pixels set to white" )  
        call get( cla,"pepper",pepper ,LIB_CLA_OPTIONAL," fraction of pixels set to black" )
        call get( cla,"b",background_ini_value ,LIB_CLA_OPTIONAL,"   value to iniialise the background of the image to" )  
        call get( cla,"sinb",sin_wave_amplitude ,LIB_CLA_OPTIONAL,"   Sets the amplitude of a sin wave added to th image from top to bottom" )
        call get( cla,"centre",test_dbg_centre ,LIB_CLA_OPTIONAL,"   put a single spots in the centre of the image for testing" ) 
        call get( cla,"maxol",max_overlap,LIB_CLA_OPTIONAL,"   Maximum allowed (elliptical) IoU, checks each spot on placement ")
        
        
        call report(cla) 
        if (.not. allRequiredArgumentsSet(cla)) stop
        if (hasHelpArgument(cla)) stop 
        call delete(cla)
        
        if (seed>0) then
            call init_random_seed(seed)
        else
            call init_random_seed()
        end if
        
        
    !---    welcome message
            
        print *,"testImageGenerator.exe"
        print *,"^^^^^^^^^^^^^^^^^^^^^^"
        print *,"    image size             : ",Nx,"x",Ny
        print *,"    random seed            : ",seed
        print *,"    output image filename        : """//trim(filename)//""""
        if (hasArgument(cla,"o")) print *,"    output .spots filename        : """//trim(filename_spots)//""""
        if (usePerlin) then
            print *,"    Perlin procedural noise length,avg,height ",perlin_lengthscale,perlin_avg,perlin_height
        else
            print *,"    no Perlin procedural noise"
        end if
        print *,"    spot count             : ",spotCount
        print *,"    spot diam dist mu      : ",spotLogNormalMu
        print *,"    spot diam dist sigma   : ",spotLogNormalSig
        print *,"    spot height mean,sigma : ",spotHeightMu,spotHeightSig
        print *,"    draw circles?          : ",circles
        print *,"    noise_stddev           : ",noise_stddev 
        print *,"    salt,pepper            : ",salt,pepper
        print *,""
        
        
    !---    If user wants to change background, is background sensible.
        if (background_ini_value>1.0d0) then
            print*,"Back ground must be between 1.0 and 0.0, defaulting to 0.0"
            background_ini_value=0.0d0
        else if (background_ini_value<0.0d0) then
            print*,"Back ground must be between 1.0 and 0.0, defaulting to 0.0"
            background_ini_value=0.0d0
        else
            print*,"Back ground initial vaule set to:",background_ini_value
        end if
    !---    create the image
        allocate(img(0:Nx-1,0:Ny-1))
        allocate(perlin_img(0:Nx-1,0:Ny-1))
        img = background_ini_value
        perlin_img = background_ini_value

    !---    Add a sin wave to the image intensity from top to bottom, use if long wavelength variation is required.
     
    !---    Check if added sinwave amplitude is sensible
        if ((abs(sin_wave_amplitude)+background_ini_value)>1.0d0)then
            print*,"Sum of sin wave amplitude (-sinb) and background intensity (-b) must be less than 1."
            print*,"Setting sin wave amplitude to zero."
            sin_wave_amplitude=0.0d0
        end if
        
        if (abs(sin_wave_amplitude)>background_ini_value) then
            print*,"Sin wave amplitude (-sinb) needs to be less or equal to the background intensity (-b) to avoid negative intnesity."
            print*,"Setting Sin wave amplitude to background intensity:", background_ini_value
            sin_wave_amplitude=background_ini_value
        end if

        if (abs(sin_wave_amplitude)>0.0d0) then
            do iy =0,Ny-1
                img(:,iy)=img(:,iy)+sin_wave_amplitude*sin(TWO_PI*((iy*1.0d0)/(Ny-1)))
            end do
        end if



        
    ! !---    generate background
    !     if (usePerlin) then
    !         call perlinNoise( perlin_lengthscale , img )
    !         img = perlin_avg + img * perlin_height
    !     end if

        allocate(spot_array(6,spotCount))

    !---    debug centeal spot mode
        if(test_dbg_centre) then
            print*,"Test mode active, count will be overidden and set to 1, this will be placed in the centre of the image"
    !----------------------------------------------------------
            logn = LogNormal_ctor(spotLogNormalMu,spotLogNormalSig)
            if (circles) then
                dd = moment(logn,1)                         !   diameter
                kk = ceiling( dd*1.1 )                      !   pixel range to search
                i2s2 = 1/(2*dd*dd)                          !   1/(2 sigma^2)
                      
                                                            
                    
                xx = gaussianVariate()                  !   intensity ( normal )
                xx = spotHeightSig*xx + spotHeightMu    !   ... scaled
                
                !   find centre of spot

                ii = floor((Nx*1.0d0)/2.0d0)

                jj = floor((Ny*1.0d0)/2.0d0)
                print *,"spot at ",ii,jj," diameter ",dd," intensity ",xx
                spot_array(1,1)=ii !save spot paramater for .spots write
                spot_array(2,1)=jj !x,y,diameter,intensity
                spot_array(3,1)=dd
                spot_array(4,1)=xx
                spot_array(5,1)=dd
                spot_array(6,1)=0.0d0

                do iy = max(0,jj-kk),min(Ny-1,jj+kk)
                    do ix = max(0,ii-kk),min(Nx-1,ii+kk)
                        rr = ( (ix-ii)*(ix-ii) + (iy-jj)*(iy-jj) )*i2s2
                        if (rr > 0.5) cycle
                        img(ix,iy) = img(ix,iy) + xx  
                        
                    end do
                end do

            else
                
                dd = spotLogNormalMu!variate(logn)                !   diameter
                kk = ceiling( dd*5 )
                i2s2 = 1/(2*dd*dd)
                            
                
                xx = gaussianVariate()
                xx = spotHeightSig*xx + spotHeightMu
                
                        
                ii = floor((Nx*1.0d0)/2.0d0)

                jj = floor((Ny*1.0d0)/2.0d0)
                print *,"spot at ",ii,jj," diameter ",dd," intensity ",xx
                spot_array(1,1)=ii !save spot paramater for .spots write
                spot_array(2,1)=jj !x,y,diameter,intensity,theta
                spot_array(3,1)=dd
                spot_array(4,1)=xx
                spot_array(5,1)=dd
                spot_array(6,1)=0.0d0
                
                
                
                do iy = max(0,jj-kk),min(Ny-1,jj+kk)
                    do ix = max(0,ii-kk),min(Nx-1,ii+kk)
                        rr = ( (ix-ii)*(ix-ii) + (iy-jj)*(iy-jj) )*i2s2
                        if (rr > 12.5d0) cycle      ! now same as other gaussian functions LIB_G2D_SEARCH_RANGE**2/2!was 12.5d0, 5 sigma  
                        
                        
                        img(ix,iy) = img(ix,iy) + xx * exp( - rr )
                        
                    end do
                end do
                    
            end if
            spotCount=-1 !so spots are not written to the image after this.
        end if 
            
    !----------------------------------------------------------

        
        
    !---    add spots
        logn = LogNormal_ctor(spotLogNormalMu,spotLogNormalSig)
        if (spotCount>0) then
        
            if (circles) then
                dd = moment(logn,1)                         !   diameter
                kk = ceiling( dd*1.1 )                      !   pixel range to search
                i2s2 = 1/(2*dd*dd)                          !   1/(2 sigma^2)
                do nn = 1,spotCount                         
                                                            
                    
                    xx = gaussianVariate()                  !   intensity ( normal )
                    xx = spotHeightSig*xx + spotHeightMu    !   ... scaled
                    
                    !   find centre of spot
                    call random_number(zeta)
                    ii = floor( zeta * (Nx+2*dd)-dd )
                    call random_number(zeta)
                    jj = floor( zeta * (Ny+2*dd)-dd )
                    print *,"spot at ",ii,jj," diameter ",dd," intensity ",xx
                    spot_array(1,nn)=ii !save spot paramater for .spots write
                    spot_array(2,nn)=jj !x,y,diameter,intensity
                    spot_array(3,nn)=dd
                    spot_array(4,nn)=xx
                    spot_array(5,nn)=dd
                    spot_array(6,nn)=0.0d0
                    
                    do iy = max(0,jj-kk),min(Ny-1,jj+kk)
                        do ix = max(0,ii-kk),min(Nx-1,ii+kk)
                            rr = ( (ix-ii)*(ix-ii) + (iy-jj)*(iy-jj) )*i2s2
                            if (rr > 0.5) cycle
                            img(ix,iy) = img(ix,iy) + xx  
                            
                        end do
                    end do
                end do
            else if (ellipses) then
            !draw elliptical Gaussians
                do nn = 1,spotCount
                    dd = variate(logn)                !   diameter
                    kk = ceiling( dd*5 )
                    i2s2 = 1/(2*dd*dd)
                                
                    
                    xx = gaussianVariate()
                    xx = spotHeightSig*xx + spotHeightMu

                    if (max_eccentricity>1.0d0) then
                        max_eccentricity=max_eccentricity-1
                        call random_number(zeta)
                        max_eccentricity=max_eccentricity*zeta
                        max_eccentricity=max_eccentricity+1
                        call random_number(zeta)
                        orientation=(TWO_PI/2.0d0)*zeta
                    end if
                    
                    spot_array(5,nn)=dd/max_eccentricity
                    call random_number(zeta)
                    ii = floor( zeta * (Nx+2*dd)-dd )
                    call random_number(zeta)
                    jj = floor( zeta * (Ny+2*dd)-dd )
                    print *,"spot at ",ii,jj,"maj diameter ",dd,"min diameter ",spot_array(5,nn)," intensity ",xx
                    spot_array(1,nn)=ii !save spot paramater for .spots write
                    spot_array(2,nn)=jj !x,y,diameter,intensity
                    spot_array(3,nn)=dd
                    spot_array(4,nn)=xx
                    
                    spot_array(6,nn)=orientation

                    dmatrix=getDMatrix(dd,spot_array(5,nn),orientation,1.0d0)
                    position(1)=ii
                    position(2)=jj
                    call drawEllipse(dmatrix,position,img,mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN,f=xx)
                    
                end do   

            else
                temp_spot_data=0.0d0
                do nn = 1,spotCount
                    dd = variate(logn)                !   diameter
                    kk = ceiling( dd*5 )
                    i2s2 = 1/(2*dd*dd)
                                
                    
                    xx = gaussianVariate()
                    xx = spotHeightSig*xx + spotHeightMu
                    
                            
                    call random_number(zeta)
                    ii = floor( zeta * (Nx+2*dd)-dd )
                    call random_number(zeta)
                    jj = floor( zeta * (Ny+2*dd)-dd )
                    if(max_overlap<1.0d0) then
                        if (nn>1) then
                            do_any_overlap=.false.
                            call check_overlap_with_prexisting_spots((/(ii*1.0d0),(jj*1.0d0),dd,xx,dd,orientation/),spot_array(:(nn-1),:),max_overlap,1.0d0,do_any_overlap)
                            if(do_any_overlap)then
                                noof_spot_placement=0
                                do ll=1,max_noof_spot_placements
                                    do_any_overlap=.false.
                                    call random_number(zeta)
                                    ii = floor( zeta * (Nx+2*dd)-dd )
                                    call random_number(zeta)
                                    jj = floor( zeta * (Ny+2*dd)-dd )
                                    call check_overlap_with_prexisting_spots((/(ii*1.0d0),(jj*1.0d0),dd,xx,dd,orientation/),spot_array(:(nn-1),:),max_overlap,1.0d0,do_any_overlap)
                                    if(.not.do_any_overlap) exit 
                                end do
                                if(do_any_overlap) print*,"testImageGenerator: WARN could not place spot that didn't overlap, consider reducing IoU tolerance"
                            end if
                        end if
                    end if
                    print *,"spot at ",ii,jj," diameter ",dd," intensity ",xx
                    spot_array(1,nn)=ii !save spot paramater for .spots write
                    spot_array(2,nn)=jj !x,y,diameter,intensity
                    spot_array(3,nn)=dd
                    spot_array(4,nn)=xx
                    spot_array(5,nn)=dd
                    spot_array(6,nn)=orientation
                    
                    
                    
                    do iy = max(0,jj-kk),min(Ny-1,jj+kk)
                        do ix = max(0,ii-kk),min(Nx-1,ii+kk)
                            rr = ( (ix-ii)*(ix-ii) + (iy-jj)*(iy-jj) )*i2s2
                            if (rr > 12.5d0) cycle      !   5 sigma
                            
                            
                            img(ix,iy) = img(ix,iy) + xx * exp( - rr )
                            
                        end do
                    end do
                    
                end do   
            end if
        end if
                
        
        
        
    !---    add noise        
        if (noise_stddev>0) then
            do jj = 0,Ny-1
                img(:,jj) = img(:,jj) + gaussianVariate(Nx)*noise_stddev
            end do
        end if
        
        if (salt+pepper>0) then
            do jj = 0,Ny-1
                do ii = 0,Nx-1
                    call random_number(zeta)
                    if (zeta < salt) then
                        img(ii,jj) = 1.0d0
                    else if (zeta < salt + pepper) then
                        img(ii,jj) = 0.0d0
                    end if
                end do
            end do
        end if
    
    !---    generate background
        if (usePerlin) then
            call perlinNoise( perlin_lengthscale , perlin_img )
            perlin_img = perlin_avg + perlin_img * perlin_height
            img=img+perlin_img
        end if
    !---    output the image        
        call writePng( filename, img )


        
        if (hasArgument(cla,"o")) then
            if (spotCount==-1) spotCount=1 !correct spot count if using single central spot test case "-centre".
            !output spots here
            !outputSpots( filename,ng,spot_array_in, nx,ny,nmPerPixel, nDeadPixels,commentLine, diamScale,bg_sigma )
            call findImageIntensityFeatures( img, bb,tt,ff,bstd,fbar,snr,fot )
            call outputSpots( filename_spots,spotCount,spot_array, Nx,Ny,1.0d0, (nint(Nx*Ny*(salt+pepper))),"", 1.0d0 )

        end if
        
        print *,""
        print *,"done"
        print *,""

        contains

        subroutine outputSpots( filename_in,ng,spot_array_in, nxin,nyin,nmPerPixel, nDeadPixels,commentLine, diamScale )
            !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !*      output as standard .spots file 
                    
                    character(len=*),intent(in)                     ::      filename_in        !   Filename in.
                    integer,intent(in)                              ::      ng              !   number of spots fitted
                    real(kind=real64),dimension(:,:),intent(in)      ::      spot_array_in   !   Object to hold 2D gaussians
                    integer,intent(in)                              ::      nxin,nyin           !   width, height of input image, in pixels
                    real(kind=real64),intent(in)                    ::      nmPerPixel      !   Length per pixel/nm
                    integer,intent(in)                              ::      nDeadPixels     !   number of pixels out-of-scope
                    character(len=*),intent(in)                     ::      commentLine     !   String to hold optional comment line
                    real(kind=real64),intent(in)                    ::      diamscale       !   diameter = nmPerPixel*(2 sigma)*diamscale
                    !real(kind=real64),intent(in),optional           ::      bg_sigma        !   background intensity standard deviation
                    
                    real(kind=real64),parameter                     ::      PI = 3.1415926535897932384626433832795d0
                    real(kind=real64),parameter                     ::      Q_DUMMY = 0.0d0
                    real(kind=real64),parameter                     ::      T_DUMMY = 0.0d0
                    real(kind=real64),parameter                     ::      THETA_DUMMY = 0.0d0
                    character(len=2)                                ::      LENGTHUNIT="px"     !length unit (default pixels)

                    
                    real(kind=real64)           ::      ww,hh,aa!,ss             !   width, height, area of input image, in nm
                    integer                     ::      ii                      !   Indices of detected spots
                     
                    ww = nxin*nmPerPixel
                    hh = nyin*nmPerPixel
                    aa = (nxin*nyin - nDeadPixels)*nmPerPixel*nmPerPixel
                    !ss = 1.0d0 ; if (present(bg_sigma)) ss = 1/bg_sigma     !   scale by noise level?
                    
                    
                    open(unit=500,file=trim(filename_in),action="write")
                        write (500,fmt='(a)') "# Clock Image generator "//trim(VERSION)
                        write (500,fmt='(a)') "# "//trim(commentLine)
                        write (500,fmt='(a,2f12.3,a)') "# ",ww,hh," micrograph extent ("//LENGTHUNIT//")"
        
                        write (500,fmt='(a2,f13.3,a)') "# ",aa,"  micrograph area ("//LENGTHUNIT//") discounting dead pixels"
                        write (500,fmt='(a2,i8,a)') "# ",ng,"  number of spots"
                        write (500,fmt='(a,100a12)') "# pos x ("//LENGTHUNIT//")  pos y ("//LENGTHUNIT//") diam 1 ("//LENGTHUNIT//") diam 2 ("//LENGTHUNIT//")"," angle(deg)"," intensity"," t* value"," quality"
                        
                        do ii = 1,ng
                            write(500,fmt='(100f12.3)')  spot_array_in(1,ii)*nmPerPixel                            &
                                                             ,spot_array_in(2,ii)*nmPerPixel                       &
                                                             ,spot_array_in(3,ii)*nmPerPixel*(diamscale)   &            !   2 for convert radius to diameter.
                                                             ,spot_array_in(5,ii)*nmPerPixel*(diamscale)  &            !   2 for convert radius to diameter.
                                                             ,spot_array_in(6,ii)*180.0d0/PI                        &            !  convert radians to degrees
                                                             ,spot_array_in(4,ii)                                  &
                                                             ,T_DUMMY                                   &
                                                             ,Q_DUMMY                             
                        end do
                    close(unit=500)
        end subroutine




        
    end program testImageGenerator