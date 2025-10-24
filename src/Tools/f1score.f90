    program f1score
    !---^^^^^^^^^^^^^^^^ 
    !*-----------------------------------------------------------------------------------------------------------------------------------
    ! *    f1score from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
    ! *    Copyright (C) 2024  Daniel Mason

    ! *    This program is free software: you can redistribute it and/or modify
    ! *    it under the terms of the GNU General Public License as published by
    ! *    the Free Software Foundation, either version 3 of the License, or
    ! *    (at your option) any later version.

    ! *    This program is distributed in the hope that it will be useful,
    ! *    but WITHOUT ANY WARRANTY; without even the implied warranty of
    ! *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ! *    GNU General Public License for more details.

    ! *    You should have received a copy of the GNU General Public License
    ! *    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    ! *-----------------------------------------------------------------------------------------------------------------------------------
            use iso_fortran_env
    
            use Lib_CommandLineArguments
            use NBAX_StringTokenizers
            use Lib_Filenames
            use Lib_ColourScale
            use Lib_Png
            use Lib_SimpleProgressBar
            use Lib_DrawEllipse
    
            implicit none
            
        !---
        !   version numbering

        !   0.0.1                   first save
        !   0.1.0       Aug 2024    James adds option for different intersection-over-union algorithms
        !   0.1.1       Sep 2024    Daniel adds option to output ground truth and comparison file into one png.
        !   
            character(len=*),parameter      ::      VERSION = "0.1.1"
            
        !---        parameters
            real(kind=real64),parameter     ::      PI = 3.1415926535897932384626433832795d0
            

            integer,parameter               ::      IOU_MODE_RECTANGLE = LIB_DRAWELLIPSE_IOU_RECTANGLE
            integer,parameter               ::      IOU_MODE_CIRCLE = LIB_DRAWELLIPSE_IOU_CIRCLE
            integer,parameter               ::      IOU_MODE_ELLIPSE = LIB_DRAWELLIPSE_IOU_ELLIPSE

            integer,parameter               ::      UNMATCHED_SPOT = -1
            integer,parameter               ::      N_COL_DEF = 8
            integer,parameter               ::      COL_X  = 1          !   input file column definitions
            integer,parameter               ::      COL_Y  = 2
            integer,parameter               ::      COL_D1 = 3
            integer,parameter               ::      COL_D2 = 4
            integer,parameter               ::      COL_THETA = 5
            integer,parameter               ::      COL_I = 6
            integer,parameter               ::      COL_T = 7
            integer,parameter               ::      COL_AIC = 8
    
            integer,parameter               ::      COL_PLOT = 6        !   output column is always number 6.
            integer,parameter               ::      COL_ID   = 8        !   spot id
            
            character(len=16),dimension(N_COL_DEF),parameter    ::      COL_NAME=(/ "Position x ",          & 
                                                                                    "Position y ",          &
                                                                                    "Major diam ",          &
                                                                                    "Minor diam ",          &
                                                                                    "Theta (deg)",          &
                                                                                    "Intensity  ",          &
                                                                                    "t* value   ",          &
                                                                                    "AIC value  "           &
                                                                                     /)
          
            character(len=16),dimension(1:3),parameter          ::      IOU_NAME=(/ "rectangle",            &
                                                                                    "circle   ",            & 
                                                                                    "ellipse  "             &
                                                                                     /)
          

    
        !---        input options
        
            type(CommandLineArguments)      ::      cla                     !Object to hold commmand line argumments (CLA)
        
            character(len=256)              ::      groundtruth_filename = ""       !   filename of ground truth .spots file
            character(len=256)              ::      spot_filename = ""              !   filename of comparison .spots file
            character(len=256)              ::      image_filename = ""             !   filename of image which ground truth and comparison are supposed to refer to
            integer                         ::      iou_mode=IOU_MODE_CIRCLE        !   0 = Rectangle, 1= circle 2=Numerical ellipse.
            real(kind=real64)               ::      iou_cutoff=0.5d0                !   what value of iou required for a match
            logical                         ::      DBG = .false.
            real(kind=real64)               ::      i_thresh = 0.0d0                !   ignore spots below this intensity threshold
            real(kind=real64)               ::      q_thresh = huge(1.0d0)                !   ignore spots below this quality threshold
            logical                         ::      ignore_unit_mismatch = .false.           !   sometimes the unit missmatch warning is unhelpful, see scaling study


        !---    read from files
            real(kind=real64),dimension(:,:),allocatable    ::      groundtruth_dat,spot_dat
            character(len=4)                ::      units = "(px)", unitsg = "(px)"
            integer                         ::      nCol
            integer                         ::      nHeaderLines
            integer                         ::      groundtruth_nSpots,spot_nSpots           
            real(kind=real64),dimension(:,:,:),allocatable    ::      rgb_img       !   output image - with coloured ellipses drawn on 
            integer                         ::      Nx,Ny                           !   pixel size of image
            real(kind=real64)               ::      nmperpx = 1.0d0                 !   scale - needed if .spots units are in nm
            logical                         ::      hasScale = .false.
            real(kind=real64)               ::      diamscale = 1.0d0               !   draw ring at 1 std dev? 2 std devs?
            real(kind=real64)               ::      diamscaleg = 1.0d0              !   draw ring at 1 std dev? 2 std devs?
    
        !---        dummy variables        
            

            integer                         ::      ii,jj,trial

            
            integer                         ::      ioerr
            real(kind=real64)               ::      ww,hh,wwg,hhg !see below
            logical                         ::      ok      !does input file exist?  
            character(len=256)              ::      dummy

            real(kind=real64),dimension(:,:),allocatable    ::      img             !   greyscale image to which ground truth and comparison are supposed to refer
            real(kind=real64),dimension(:,:),allocatable    ::      distance
            integer,dimension(:),allocatable                ::      indx,indxg
            integer,dimension(2)                            ::      pair
            logical,dimension(:),allocatable                ::      good
            integer                                         ::      tp,fn,fp
            real(kind=real64)                               ::      f1, average_iou
            real(kind=real64),dimension(2,2)                ::      D_matrix, D_matrixg        
     
        !---    read command line arguments
            cla = CommandLineArguments_ctor(30)  
             
            call setProgramDescription( cla, "f1score" )
            call setProgramVersion( cla, VERSION )   
           
                        
            
        !---    filename options        
            call get( cla,"f",spot_filename ,LIB_CLA_REQUIRED,"           comparison .spots filename" )                                                        
            call get( cla,"g",groundtruth_filename ,LIB_CLA_OPTIONAL,"         ground truth .spots filename" )
            call get( cla,"b",image_filename ,LIB_CLA_OPTIONAL,"         background image filename" )
            dummy = IOU_NAME(iou_mode)
            call get( cla,"m",dummy ,LIB_CLA_OPTIONAL,"         How to calculate the IoU (circle/rectangle,ellise)" )
            if (hasArgument( cla,"m" ) ) then
                iou_mode = -1
                do ii = 1,ubound(IOU_NAME,dim=1)
                    if (trim(dummy)==trim(IOU_NAME(ii))) then
                        iou_mode = ii
                        exit
                    end if
                end do
                if (iou_mode == -1) then
                    print *,"f1score error - iou mode """//trim(dummy)//""" not recognised- use one of"
                    do ii = 1,ubound(IOU_NAME,dim=1)
                        print *,""""//trim(IOU_NAME(ii))//""""
                    end do
                    iou_mode = IOU_MODE_CIRCLE
                    print *,"setting default """//trim(IOU_NAME(iou_mode))//""""
                    print *,""
                end if
            end if
            call get( cla,"u",iou_cutoff ,LIB_CLA_OPTIONAL,"         What cutoff to use for the IoU" )         
            iou_cutoff = max(1.0d-8,iou_cutoff)                            
            call get( cla,"s",nmperpx ,LIB_CLA_OPTIONAL,"         scale of .png image (nm per px)" )
            hasScale = hasArgument( cla,"s" )       !       will need to check later if the units in the .spots file are "nm", then I'll need the scale. 
            call get( cla,"scale",diamscale ,LIB_CLA_OPTIONAL,"     diameter scaling (number of std devs to place ring), multiplier on diam in comparison .spots file" )                                                        
            call get( cla,"scaleg",diamscaleg ,LIB_CLA_OPTIONAL,"     diameter scaling (number of std devs to place ring), multiplier on diam in ground truth .spots file" )                                                                    
            call get( cla,"dbg",DBG ,LIB_CLA_OPTIONAL,"       debug output" )       
            call get( cla,"i",i_thresh ,LIB_CLA_OPTIONAL,"         intensity threshold" )   
            call get( cla,"q",q_thresh ,LIB_CLA_OPTIONAL,"         quality threshold" )
            call get( cla,"nomiswarn",ignore_unit_mismatch ,LIB_CLA_OPTIONAL,"         ignore unit mismatch" )
            ignore_unit_mismatch= hasArgument( cla,"nomiswarn" )
                                                             
            call report(cla)
            if (.not. allRequiredArgumentsSet(cla)) stop
            if (hasHelpArgument(cla)) stop
            call delete(cla)
            
        !---        
            
            
            
    
    
    
            print *,"f1score v."//VERSION
            print *,"^^^^^^^^^^"//repeat("^",len_trim(VERSION))
            print *,""

            print *,"    ground truth spot filename """//trim(groundtruth_filename)//""""
            print *,"    comparison spot filename   """//trim(spot_filename)//""""
            print *,"    IOU mode            : """//trim(IOU_NAME(iou_mode))//""""
            print*,"    IoU cutoff          : ",iou_cutoff
            print*,"    nm per pixel        : ",nmperpx
            print*,"    diameter scaling    : ",diamscale," (comparison file)"      
            print*,"    diameter scaling    : ",diamscaleg," (ground truth file)"      
            print*,"    intensity threshold : ",i_thresh         
            print*,"    quality threshold : ",q_thresh     
            print *,""
    
            !---    check that the file exists        
            inquire(file=trim(spot_filename),exist=ok)
            if (.not. ok) then
                 print *,"f1score error - couldn't find .spots file """//trim(spot_filename)//""""
                 stop
            end if
    
            if (len_trim(groundtruth_filename)/=0) then    
                inquire(file=trim(groundtruth_filename),exist=ok)
                if (.not. ok) then
                    print *,"f1score error - couldn't find ground truth .spots file """//trim(groundtruth_filename)//""""
                    stop
                end if
            end if 
            
            if (len_trim(image_filename)/=0) then    
                inquire(file=trim(image_filename),exist=ok)
                if (.not. ok) then
                    print *,"f1score error - couldn't find image file """//trim(image_filename)//""""
                    stop
                end if
                print *,"f1score info - reading image file """//trim(image_filename)//""""
                call readPng( image_filename,img )
                Nx = size(img,dim=1)
                Ny = size(img,dim=2)
                allocate(rgb_img(3,Nx,Ny))
                do jj = 1,Ny
                    do ii = 1,Nx
                        rgb_img(:,ii,jj) = img(ii,jj)
                    end do
                end do
                deallocate(img)
            end if 
            
        !---    welcome        
            print *,"   .spots filename                 """//trim(spot_filename)//""""
            print *,"   ground truth .spots filename    """//trim(groundtruth_filename)//""""
            
        !---    job starts here    
            call readSpotsFile(spot_filename,spot_dat,ww,hh,nHeaderLines,units) 
            call readSpotsFile(groundtruth_filename,groundtruth_dat,wwg,hhg,nHeaderLines,unitsg)  

            if((ww .ne. wwg).or.(hh .ne. hhg)) then 
                print*,"f1score WARNING - ground truth and detection images appear to be different dimensions"
                print*,"    (read from .spots file, are you comparing the right files?)"
            end if
            if (trim(units)/=trim(unitsg)) then
                print*,"f1score ERROR - ground truth and detection images are using different units"
                print*,"    (read from .spots file, are you comparing the right files?)"
                if(ignore_unit_mismatch) then
                    print*,"f1score WARNING - ignoring unit mismatch at your request."
                else
                    stop
                end if
            end if
            if ((len_trim(image_filename)/=0) .and. (trim(unitsg)/="(px)")) then
                if (.not. hasScale) then
                    print*,"f1score ERROR - ground truth image uses units """//trim(unitsg)//""""
                    print*,"and a comparison image is required, but scale (-s) is not defined"
                    stop
                end if                
            end if  


        !---    convert data from nm scale to pixel scale. 
        !       this is a useful if comparing intersection over union with ellipses, 
        !       as this is computed with bitmaps, but vital if also producing an output image with ground truth and comparison overlaid.
            groundtruth_dat(COL_X,:)  = groundtruth_dat(COL_X,:)  / nmperpx
            groundtruth_dat(COL_Y,:)  = groundtruth_dat(COL_Y,:)  / nmperpx
            groundtruth_dat(COL_D1,:) = groundtruth_dat(COL_D1,:) / nmperpx
            groundtruth_dat(COL_D2,:) = groundtruth_dat(COL_D2,:) / nmperpx
            spot_dat(COL_X,:)  = spot_dat(COL_X,:)  / nmperpx
            spot_dat(COL_Y,:)  = spot_dat(COL_Y,:)  / nmperpx
            spot_dat(COL_D1,:) = spot_dat(COL_D1,:) / nmperpx
            spot_dat(COL_D2,:) = spot_dat(COL_D2,:) / nmperpx





            nCol = size(groundtruth_dat,dim=1)
            groundtruth_nSpots = size(groundtruth_dat,dim=2)
            spot_nSpots = size(spot_dat,dim=2)

        !---    compute distance between spots in the files
            print *,""
            print *,"compute distance between spots in the files"
            allocate(distance(spot_nSpots,groundtruth_nSpots))

            do jj = 1,groundtruth_nSpots
                D_matrixg = getDMatrix(groundtruth_dat(COL_D1,jj )*diamscaleg/2, groundtruth_dat(COL_D2,jj )*diamscaleg/2, groundtruth_dat(COL_THETA,jj)*(PI/180))
                do ii = 1,spot_nSpots
                    D_matrix = getDMatrix(spot_dat(COL_D1,ii )*diamscale/2, spot_dat(COL_D2,ii )*diamscale/2, spot_dat(COL_THETA,ii)*(PI/180))
                    distance(ii,jj) = 1 -intersectionOverUnion( D_matrix,spot_dat(COL_X:COL_Y,ii),                 &
                                                                D_matrixg,groundtruth_dat(COL_X:COL_Y,jj),         &
                                                                mode = iou_mode)  
                end do
            end do  

            
        !---    for each spot in turn, find the best match.
            print *,""
            average_iou=0.0d0
            allocate(indx(spot_nSpots)) 
            allocate(indxg(groundtruth_nSpots)) 
            indx = UNMATCHED_SPOT
            indxg = UNMATCHED_SPOT
            if (spot_nSpots >= groundtruth_nSpots) then  
                allocate(good(groundtruth_nSpots))
                good=.false.
                do trial = 1,groundtruth_nSpots
                    pair = minloc( distance )
                    ii = pair(1)            !   comparison spot
                    jj = pair(2)            !   ground truth spot                    
                    
                    !indx(pair(2)) = pair(1)
            !        if (iou_mode==0) then !do circle IoU
                    good(jj) = distance(ii,jj) <= 1 - iou_cutoff
                    

                    if (good(jj)) then
                        indx(ii) = jj
                        indxg(jj) = ii
                        if (DBG) print *,"matched ",jj," at ",groundtruth_dat(COL_X:COL_Y,jj)," to ",ii," at ",spot_dat(COL_X:COL_Y,ii)," d ",distance(ii,jj),good(jj)
                        average_iou=average_iou+(1-distance(ii,jj)) ! was average_iou=average_iou+distance(ii,jj) ???
                    end if
                    distance(:,jj) = huge(1.0)
                    distance(ii,:) = huge(1.0)

                end do

                average_iou=average_iou/spot_nSpots
                tp = count(good)                                !   true positive counts = number of good matches
                fp = spot_nSpots-groundtruth_nSpots             !   false positive counts = number of spots with no match in ground truth
                fn = groundtruth_nSpots - tp                    !   false negative counts = number of ground truth spots without a match

            else
                !allocate(indx(spot_nSpots))
                allocate(good(spot_nSpots))
                do trial = 1,spot_nSpots
                    pair = minloc( distance )
                    ii = pair(1)            !   comparison spot
                    jj = pair(2)            !   ground truth spot                    
                    !indx(pair(1)) = pair(2)
                    good(ii) = distance(ii,jj) <= 1 - iou_cutoff
                    if (good(ii)) then
                        indx(ii) = jj
                        indxg(jj) = ii
                        if (DBG) print *,"matched ",ii," at ",spot_dat(COL_X:COL_Y,ii)," to ",jj," at ",groundtruth_dat(COL_X:COL_Y,jj)," d ",distance(ii,jj),good(ii)
                        average_iou=average_iou+(1-distance(ii,jj)) ! was average_iou=average_iou+distance(ii,jj) ???
                    end if
                    distance(:,jj) = huge(1.0)
                    distance(ii,:) = huge(1.0)

                end do

                average_iou=average_iou/groundtruth_nSpots
                tp = count(good)                                !   true positive counts = number of good matches
                fp = spot_nSpots - tp                           !   false positive counts = number of spots with no match in ground truth
                fn = groundtruth_nSpots - spot_nSpots           !   false negative counts = number of ground truth spots without a match

            end if
 
        !---    compute the f1 score
            
            f1 = tp*2.0d0 / ( 2.0d0*tp + fp + fn )
            print *,""
            write(*,fmt='(7a16)')   " ground nSpots "," spots nSpots ","true positive", "false positive", "false negative", "F1","<IoU>"
            write(*,fmt='(5i16,2f16.8)') groundtruth_nspots,spot_nspots,tp,fp,fn,f1,average_iou
   
        !---    construct the output image, if needed
            if (len_trim(image_filename)/=0) then
                image_filename = trim( removeSuffix(image_filename) )//".comparison.png"
                print *,"f1score info - constructing output image file """//trim(image_filename)//""""
                
                do ii = 1,groundtruth_nSpots
                    D_matrixg = getDMatrix( groundtruth_dat(COL_D1,ii )*diamscaleg/2, groundtruth_dat(COL_D2,ii )*diamscaleg/2, groundtruth_dat(COL_THETA,ii)*(PI/180))

                    if (indxg(ii)/=UNMATCHED_SPOT) then
                        call drawEllipse(D_matrixg,groundtruth_dat(COL_X:COL_Y,ii ),rgb_img,colourscale=COLOURSCALE_REDBLUE,f=0.70d0,mode=LIB_DRAWELLIPSE_SHADE_RING)                    
                    else
                        call drawEllipse(D_matrixg,groundtruth_dat(COL_X:COL_Y,ii ),rgb_img,colourscale=COLOURSCALE_REDBLUE,f=0.5d0,mode=LIB_DRAWELLIPSE_SHADE_RING)                    
                    end if

                end do
                
                do ii = 1,spot_nSpots
                    D_matrix = getDMatrix( spot_dat(COL_D1,ii )*diamscale/2, spot_dat(COL_D2,ii )*diamscale/2, spot_dat(COL_THETA,ii)*(PI/180))
                    if (indx(ii)/=UNMATCHED_SPOT) then
                        call drawEllipse(D_matrix,spot_dat(COL_X:COL_Y,ii ),rgb_img,colourscale=COLOURSCALE_RAINBOW,f=0.175d0,mode=LIB_DRAWELLIPSE_SHADE_RING)                    
                    else
                        call drawEllipse(D_matrix,spot_dat(COL_X:COL_Y,ii ),rgb_img,colourscale=COLOURSCALE_RAINBOW,f=0.250d0,mode=LIB_DRAWELLIPSE_SHADE_RING)                    
                    end if
                end do
                
                print *,"f1score info - writing output image file """//trim(image_filename)//""""
                print *,"f1score info - blue = ground truth (true positive), white = ground truth (false negative), red = matched, yellow = unmatched"
                
                call write_rgb_png( image_filename,rgb_img ) !also jph eddit
                print *,""
            end if 



            print *,""
            print *,"done"
            print *,""
    
    
        contains
    !---^^^^^^^^    
    
            pure real(kind=real64) function spot_distance( spot_dat,groundtruth_dat )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !       simple distance metric 
        !       https://mathworld.wolfram.com/Circle-CircleIntersection.html
                real(kind=real64),dimension(:),intent(in)       ::      spot_dat,groundtruth_dat
                real(kind=real64)           ::      spot_r2bar,groundtruth_r2bar
                real(kind=real64)           ::      dx,dy,d2,dd
                real(kind=real64)           ::      area_intersection,area_total

                real(kind=real64)           ::      spot_rbar,groundtruth_rbar
                spot_r2bar = spot_dat(COL_D1)*spot_dat(COL_D2)/4                            !   1/4 converts diameter^2 to radius^2
                groundtruth_r2bar = groundtruth_dat(COL_D1)*groundtruth_dat(COL_D2)/4
                dx = spot_dat(COL_X) - groundtruth_dat(COL_X)
                dy = spot_dat(COL_Y) - groundtruth_dat(COL_Y)
                d2 =  dx*dx + dy*dy

                if (d2 < 1.0d-8) then
                    !   spots overlap perfectly 
                    area_intersection = PI*min(spot_r2bar,groundtruth_r2bar)
                else if (min(spot_r2bar,groundtruth_r2bar) < 1.0d-8) then
                    !   one spot has zero area
                    area_intersection = 0
                else
                    dd = sqrt( d2 )
                    spot_rbar = sqrt( spot_r2bar )
                    groundtruth_rbar = sqrt( groundtruth_r2bar )
                    area_intersection = spot_r2bar       *acos( (d2 + spot_r2bar - groundtruth_r2bar)/(2*dd*spot_rbar) )                        &
                                      + groundtruth_r2bar*acos( (d2 - spot_r2bar + groundtruth_r2bar)/(2*dd*groundtruth_rbar) )                 &
                                      - sqrt( (-dd+spot_rbar+groundtruth_rbar)*(dd+spot_rbar-groundtruth_rbar)*(dd-spot_rbar+groundtruth_rbar)*(dd+spot_rbar+groundtruth_rbar) )/2
                end if

                area_total = PI*(spot_r2bar + groundtruth_r2bar) - area_intersection

                spot_distance = 1 - area_intersection/(max(1.0d-8,area_total))

                return
            end function spot_distance
 
                
                
            subroutine readSpotsFile(spot_filename,dat,ww,hh,nHeaderLines,units)  
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                character(len=*),intent(in)                                 ::      spot_filename
                real(kind=real64),dimension(:,:),intent(out),allocatable    ::      dat  
                real(kind=real64),intent(out)                               ::      ww,hh
                integer,intent(out)                                         ::      nHeaderLines
                character(len=4),intent(out)                                ::      units
                integer                                 ::      ii,jj,nCol
                integer                                 ::      nSpots,nOverThresh
                character(len=256)                      ::      dummy
                real(kind=real64),dimension(100)        ::      datline
                logical                                 ::      isOverThresh
    
        !---    read the header of the spots file
                open(unit=500,file=trim(spot_filename),action="read")
                !   try to read the header. I'm expecting two important lines:
                !   "micrograph extent (px)"
                !   "number of spots"
                    ww = -1 ; hh = -1       !   size of micrograph in nm
                    nSpots = - 1
                    nHeaderLines = -1
    
                    do ii = 1,20            !   I'm guessing max 20 lines of header junk.
                        read(unit=500,fmt='(a)',iostat=ioerr) dummy
                        dummy = adjustl(dummy)
                        if (ioerr /= 0) exit
                        if (index(dummy,"extent")/=0) then  
                            read(dummy(2:),fmt=*,iostat=ioerr) ww,hh                    
                            if (ioerr /= 0) exit
                            dummy = adjustl( dummy(index(dummy,"extent")+6:))
                            units = dummy(1:4)
                            print *,"ww,hh,units",ww,hh,units
                        else if (index(dummy,"number")/=0) then  
                            read(dummy(2:),fmt=*,iostat=ioerr) nSpots
                            if (ioerr /= 0) exit
                        else if (.not. (dummy(1:1)=="#" .or. dummy(1:1)=="$" .or. dummy(1:1)=="%" .or. dummy(1:1)=="!")) then    
                            nHeaderLines = ii-1 
                            call parse(dummy,datline,nCol)
                            exit
                        end if
                    end do    
                    if ( min(ww,hh)<=0 ) stop "spotsToPng::readSpotsFile error - could not read width or height from .spots file"        
                    if ( nSpots<=0 ) stop "spotsToPng::readSpotsFile error - could not read number of spots from .spots file"        
    
                    print *,"spotsToPng::readSpotsFile info - micrograph width, height = ",ww,hh,units
                    print *,"spotsToPng::readSpotsFile info - number of spots = ",nSpots
                    
    
                close(unit=500)
    
            !---    count number of spots over threshold
                nOverThresh = 0
                open(unit=501,file=trim(spot_filename),action="read")
                    do ii = 1,nHeaderLines
                        read(unit=501,fmt='(a)',iostat=ioerr) dummy
                    end do
                    do ii = 1,nSpots
                        read(unit=501,fmt='(a)',iostat=ioerr) dummy
                        call parse(dummy,datline,jj)

                        if (ioerr /= 0) then    
                            print *,"spotsToPng::readSpotsFile error - could not read data for spot ",ii
                            stop
                        end if
                        if (jj/=nCol) then    
                            print *,"spotsToPng::readSpotsFile error - could not parse sufficient columns for spot",ii
                            stop
                        end if
                        isOverThresh=.false.
                        isOverThresh=(datline(COL_I) >= i_thresh).and.(datline(COL_AIC) <= q_thresh)
                        !isOverThresh=(datline(COL_AIC) <= q_thresh)
                        !if (datline(COL_I) >= i_thresh) nOverThresh = nOverThresh + 1
                        if (isOverThresh) nOverThresh = nOverThresh + 1
                    end do
                close(unit=501)
                print *,"spotsToPng::readSpotsFile info - number of spots over threshold= ",nOverThresh
                

    
                    
            !---    read in the spots data
                allocate(dat(nCol,nOverThresh))
                nOverThresh = 0
                open(unit=501,file=trim(spot_filename),action="read")
                    do ii = 1,nHeaderLines
                        read(unit=501,fmt='(a)',iostat=ioerr) dummy
                    end do
                    do ii = 1,nSpots
                        read(unit=501,fmt='(a)',iostat=ioerr) dummy
                        call parse(dummy,datline,jj)
                        isOverThresh=.false.
                        isOverThresh=(datline(COL_I) >= i_thresh).and.(datline(COL_AIC) <= q_thresh)
                        !if (datline(COL_I) >= i_thresh) then
                        if (isOverThresh) then
                            nOverThresh = nOverThresh + 1
                            dat(1:nCol,nOverThresh) = datline(1:nCol)
                        end if

                        !if (ii==1) print*,"DBG|first line of read",trim(dummy),dat(1:nCol,ii)
        
                    end do
                close(unit=501)
                nSpots = nOverThresh


            end subroutine readSpotsFile
            
        end program f1score