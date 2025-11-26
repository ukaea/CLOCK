!*-----------------------------------------------------------------------------------------------------------------------------------
! *    multi_f1score from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
! *    Copyright (C) 2024  Daniel Mason & James Heath

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

program multif1score
    !---^^^^^^^^^^^^^^^^ 

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

!   0.0.1                   first save after converting origional f1score
!   
    character(len=*),parameter      ::      VERSION = "0.0.1"
    
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
    integer,parameter               ::      COL_Q = 8

    integer,parameter               ::      COL_PLOT = 6        !   output column is always number 6.
    integer,parameter               ::      COL_ID   = 8        !   spot id
    
    character(len=16),dimension(N_COL_DEF),parameter    ::      COL_NAME=(/ "Position x ",          & 
                                                                            "Position y ",          &
                                                                            "Major diam ",          &
                                                                            "Minor diam ",          &
                                                                            "Theta (deg)",          &
                                                                            "Intensity  ",          &
                                                                            "t* value   ",          &
                                                                            "Quality    "           &
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
    character(len=256)              ::      res_outfile                     !   results filename respectively
    integer                         ::      iou_mode=IOU_MODE_CIRCLE        !   0 = Rectangle, 1= circle 2=Numerical ellipse.
    real(kind=real64)               ::      iou_cutoff=0.5d0                !   what value of iou required for a match
    logical                         ::      DBG = .false.
    real(kind=real64)               ::      i_thresh = 0.0d0                !   ignore spots below this intensity threshold
    real(kind=real64)               ::      q_thresh = huge(1.0d0)          !   ignore spots above this quality threshold
    !logical                         ::      cut_top=.false., cut_bot=.false.    !cut spots that are not in the top, bottom,
    !ogical                         ::      cut_left=.false., cut_right=.false. !left and right edges of the image
    logical                         ::      ignore_unit_mismatch = .false.           !   sometimes the unit missmatch warning is unhelpful, see scaling study
!---        Define caculation matrix
    real(kind=real64)               ::      i_thresh_max = 1.0d0                !   ignore spots below this intensity threshold
    real(kind=real64)               ::      i_thresh_min = 0.0d0                !   ignore spots below this intensity threshold
    real(kind=real64)               ::      q_thresh_max = 0.0d0                !   ignore spots min this quality threshold
    real(kind=real64)               ::      q_thresh_min = huge(1.0d0)          !   ignore spots above this quality threshold
    real(kind=real64)               ::      diamscale_min = 1.0d0               !   Min noof std devs to draw ellipse - scales detected not gt.
    real(kind=real64)               ::      diamscale_max = 1.0d0              !   Max noof std devs to draw ellipse
    integer                         ::      noof_i_thresh                       !   number of intensity thresholds to trial
    integer                         ::      noof_q_thresh                       !   number of quality thresholds to trial.
    integer                         ::      noof_d_scale                        !   number of diameter scales to trail (on scales detected not gt)

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
    logical                         ::      max_recall = .false.            !   Maximise reccall instead of F1

!---        dummy variables        
    
    real(kind=real64),dimension(:,:,:),allocatable  ::      Dg,Df           !D matrix of groundtruth and detected spots respectively.

    integer                         ::      ii,jj,trial,kk,ll,mm

    
    integer                         ::      ioerr
    real(kind=real64)               ::      ww,hh,dd,ff,tt,aa,wwg,hhg !see below
    logical                         ::      ok      !does input file exist?  
    character(len=256)              ::      dummy

    real(kind=real64),dimension(:,:),allocatable    ::      img             !   greyscale image to which ground truth and comparison are supposed to refer
    real(kind=real64),dimension(:,:),allocatable    ::      distance, looped_distance
    integer,dimension(:),allocatable                ::      indx,indxg,best_indx,best_indxg
    integer,dimension(2)                            ::      pair
    logical,dimension(:),allocatable                ::      good
    integer                                         ::      tp,fn,fp
    real(kind=real64)                               ::      f1, average_iou
    real(kind=real64),dimension(2,2)                ::      D_matrix, D_matrixg
    real(kind=real64),dimension(:,:,:),allocatable  ::      array_D_matrixg
    real(kind=real64),dimension(:,:,:),allocatable  ::      average_iou_matrix
    real(kind=real64)                               ::      current_average_iou,best_average_iou
    real(kind=real64)                               ::      current_ithresh, current_qthresh, current_dscale
    real(kind=real64)                               ::      best_ithresh, best_qthresh, best_dscale             !to get highest <IoU>
    logical                                         ::      in_threshold !Keep a spot with current thresholds
    integer                                         ::      loop_noof_good,loop_noof_good_gt
    integer                                         ::      best_fn, best_fp, best_tp, best_noof_GT, best_noof_spots
    real(kind=real64)                               ::      best_f1, recall, best_recall

    integer,dimension(3)                            ::      f1_dat
    integer                                         ::      noof_valid_spots !number of detected spots that fall within the current thresholds
    !logical                                         ::      escape_edge_logic!If the user sets all the edge cuts to false don't check each edge
    real(kind=real64),dimension(:),allocatable      ::      i_thresh_array,d_thresh_array,f1_array
    integer                                         ::      noof_runs, results_indx



!---    read command line arguments
    cla = CommandLineArguments_ctor(45)  
        
    call setProgramDescription( cla, "f1score" )
    call setProgramVersion( cla, VERSION )   
    
                
    
!---    filename options        
    call get( cla,"f",spot_filename ,LIB_CLA_REQUIRED,"           comparison .spots filename" )                                                        
    call get( cla,"g",groundtruth_filename ,LIB_CLA_OPTIONAL,"         ground truth .spots filename" )
    call get( cla,"b",image_filename ,LIB_CLA_OPTIONAL,"         background image filename" )
    call get( cla,"o",res_outfile ,LIB_CLA_OPTIONAL,"         output filename" )
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
    call get( cla,"scaleg",diamscaleg ,LIB_CLA_OPTIONAL,"    diameter scaling (number of std devs to place ring), multiplier on diam in ground truth .spots file" )                                                                    
    call get( cla,"dbg",DBG ,LIB_CLA_OPTIONAL,"       debug output" )       
    call get( cla,"i",i_thresh ,LIB_CLA_OPTIONAL,"         intensity threshold" )   
    call get( cla,"q",q_thresh ,LIB_CLA_OPTIONAL,"         quality threshold" )
    call get( cla,"rec",max_recall ,LIB_CLA_OPTIONAL,"       maximise recall" )   
    
    
    call get( cla,"ni",noof_i_thresh ,LIB_CLA_OPTIONAL,"        number of intensity thresholds to trial" )
    call get( cla,"nq",noof_q_thresh ,LIB_CLA_OPTIONAL,"        number of quality thresholds to trial" )
    call get( cla,"nd",noof_d_scale ,LIB_CLA_OPTIONAL,"        number of diameter scales to trial" )

    call get( cla,"imin",i_thresh_min ,LIB_CLA_OPTIONAL,"      Minimum intensity threshold to trial" )
    call get( cla,"imax",i_thresh_max ,LIB_CLA_OPTIONAL,"      Maximum intensity threshold to trial" )
    call get( cla,"dmin",diamscale_min ,LIB_CLA_OPTIONAL,"      Minimum  diameter scale to trial" )
    call get( cla,"dmax",diamscale_max ,LIB_CLA_OPTIONAL,"      Maximum  diameter scale to trial" )
    call get( cla,"qmin",q_thresh_min ,LIB_CLA_OPTIONAL,"      Minimum of quality threshold to trial" )
    call get( cla,"qmax",q_thresh_max ,LIB_CLA_OPTIONAL,"      Maximum of quality threshold to trial" )
    call get( cla,"nomiswarn",ignore_unit_mismatch ,LIB_CLA_OPTIONAL," ignore unit mismatch" )
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

    if (hasArgument(cla,"rec")) then
        max_recall=.true.
        print*,"Maximising recall instead of F1 at users request"
    end if

    !--- multi F1 bits
    noof_d_scale=noof_d_scale+1
    noof_i_thresh=noof_i_thresh+1
    noof_q_thresh=noof_q_thresh+1
    print*,"multi F1, number of calculations to do:",noof_d_scale*noof_i_thresh*noof_q_thresh
    !allocate(average_iou_matrix(noof_d_scale,noof_i_thresh,noof_q_thresh))
    allocate(average_iou_matrix(noof_q_thresh,noof_i_thresh,noof_d_scale))
    i_thresh=0.0d0 !need to include all spots in read
    q_thresh=huge(1.0d0)
    diamscale=1.0d0
    diamscaleg=1.0d0
    i_thresh_min=0.0d0
    i_thresh_max=1.0d0

    noof_runs=max(1,noof_d_scale)*max(1,noof_i_thresh)*max(1,noof_q_thresh)
    allocate(f1_array(noof_runs))
    allocate(d_thresh_array(noof_runs))
    allocate(i_thresh_array(noof_runs))
    results_indx=0
    f1_array=-1.0d0
    d_thresh_array=-1.0d0
    i_thresh_array=-1.0d0


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

    print*,"the max quality value is:",maxval(spot_dat(COL_Q,:))

    nCol = size(groundtruth_dat,dim=1)
    groundtruth_nSpots = size(groundtruth_dat,dim=2)
    spot_nSpots = size(spot_dat,dim=2)

    allocate(indx(spot_nSpots)) !need these out of loop for final image write
    allocate(indxg(groundtruth_nSpots)) 
    allocate(best_indx(spot_nSpots)) !need these out of loop for final image write
    allocate(best_indxg(groundtruth_nSpots))
    
    !initialise storage vars
    allocate(distance(spot_nSpots,groundtruth_nSpots))
    best_average_iou=0.0d0
    current_average_iou=0.0d0
    best_f1=0.0d0
    best_recall=0.0d0

    best_fn=0
    best_tp=0
    best_fp=0

    best_dscale=0.0d0
    best_ithresh=0.0d0
    best_qthresh=0.0d0

    best_noof_GT=0
    best_noof_spots=0
    loop_noof_good=0
    loop_noof_good_gt=0

    !for all diameter scales
    do kk=1,noof_d_scale
        !get diam scale
        diamscale=diamscale_min+((diamscale_max-diamscale_min)*((kk-1.0d0)/(noof_d_scale-1.0d0))) 
    !---for all  intensity thresholds
        do ll=1,noof_i_thresh !for all intenisty thresholds
            !get instensity threshold
            i_thresh=i_thresh_min+(i_thresh_max-i_thresh_min)*((ll-1.0d0)/(noof_i_thresh-1.0d0))
    !-------for all quality thresholds
            do mm=1,noof_q_thresh !for all quality thresholds
                !get quality threshold
                q_thresh=huge(1.0d0)!q_thresh_min+(q_thresh_max-q_thresh_min)*((mm-1.0d0)/(noof_q_thresh-1.0d0))
                !for this loop, how many spots meet the thresholds
                noof_valid_spots=spot_nSpots
                !for this loop, how many spots are sucessfully matched
                loop_noof_good=spot_nSpots
                !get distance matrix
                    !for all detected spots
                do ii = 1,spot_nSpots
                    !is the spot in the intensity threshold
                    if(spot_dat(COL_I,ii)<i_thresh) then
                        noof_valid_spots=noof_valid_spots-1!discard it from count
                    else if (spot_dat(COL_Q,ii)>q_thresh) then
                        noof_valid_spots=noof_valid_spots-1!discard it from count
                    !else if(.not.(is_spot_in_boundary(cut_top,cut_bot,cut_left,cut_right,spot_dat(COL_X:COL_THETA,ii),ww,hh,diamscale))) then
                        !print*,"DBG| removed point ii, pos(x,y)",ii,spot_dat(COL_X:COL_Y,ii)
                        !noof_valid_spots=noof_valid_spots-1!discard it from count
                    end if
                end do
                !for all GT spots
                do jj = 1,groundtruth_nSpots
                    D_matrixg = getDMatrix(groundtruth_dat(COL_D1,jj )*diamscaleg/2, groundtruth_dat(COL_D2,jj )*diamscaleg/2, groundtruth_dat(COL_THETA,jj)*(PI/180))
                    !for all detected spots
                    do ii = 1,spot_nSpots
                        !is the spot in the intensity threshold
                        if(spot_dat(COL_I,ii)<i_thresh) then
                            distance(ii,jj) = huge(1.0d0)!discard it from distance calc
                            
                        else if (spot_dat(COL_Q,ii)>q_thresh) then
                            distance(ii,jj) = huge(1.0d0)!discard it from distance calc
                            
                        else !in thresholds so calc distance
                            D_matrix = getDMatrix(spot_dat(COL_D1,ii )*diamscale/2, spot_dat(COL_D2,ii )*diamscale/2, spot_dat(COL_THETA,ii)*(PI/180))
                            distance(ii,jj) = 1 -intersectionOverUnion( D_matrix,spot_dat(COL_X:COL_Y,ii),                 &
                                                                D_matrixg,groundtruth_dat(COL_X:COL_Y,jj),         &
                                                                mode = iou_mode)  

                        end if
                    end do
                end do 
                average_iou=0.0d0 !initialise
                recall=0.0d0
                !calc F1
                call get_f1_score(distance,iou_cutoff,f1_dat,f1,average_iou,indxg,indx,noof_valid_spots)
                !is new F1 <IoU> the best
                current_average_iou=average_iou
                recall=f1_dat(1)/(real((f1_dat(1)+f1_dat(3)),kind=real64) ) !  
                !if (current_average_iou>best_average_iou) then

                if (max_recall) then
                    if (recall>best_recall) then
                        best_recall=recall
                        best_average_iou=current_average_iou
                        best_f1=f1
                        best_tp=f1_dat(1)
                        best_fp=f1_dat(2)
                        best_fn=f1_dat(3)
                        best_dscale=diamscale
                        best_ithresh=i_thresh
                        best_qthresh=q_thresh
                        best_noof_spots=noof_valid_spots
                        best_indx=indx
                        best_indxg=indxg
                    end if
                else
                    if (f1>best_f1) then
                    !if (recall>best_recall) then
                        !store if so
                        best_recall=recall
                        best_average_iou=current_average_iou
                        best_f1=f1
                        best_tp=f1_dat(1)
                        best_fp=f1_dat(2)
                        best_fn=f1_dat(3)
                        best_dscale=diamscale
                        best_ithresh=i_thresh
                        best_qthresh=q_thresh
                        best_noof_spots=noof_valid_spots
                        best_indx=indx
                        best_indxg=indxg

                    end if
                end if
                results_indx=results_indx+1
                f1_array(results_indx)=f1
                d_thresh_array(results_indx)=diamscale
                i_thresh_array(results_indx)=i_thresh
                    
            end do
        end do

        write(*,fmt='(f6.2,a)',advance="no") ((real(kk)/real(noof_d_scale))*100.0d0),"% "
    end do
    write(*,fmt='(a)',advance="yes") " "

    print *,""
    print*,"noof(dscale,ithresh,qthresh):",noof_d_scale,noof_i_thresh,noof_q_thresh
    print*,"d_scale(min,max)", diamscale_min,diamscale_max,"i_thresh(min,max)",i_thresh_min,i_thresh_max,"q_thresh(min,max):",q_thresh_min,q_thresh_max
    write(*,fmt='(1a16,1f16.8,1a26,1f16.8,1a16,1f16.8,1a16,1f16.8)')   " Best <IoU> is: ",best_average_iou,", using diam scale of: ", best_dscale, ", Qthresh of: ", best_qthresh,", i_thresh of: ",best_ithresh
    write(*,fmt='(8a16,1a18)')   " ground nSpots "," spots nSpots ","true positive", "false positive", "false negative", "F1","<IoU>","Recall","spots in threshold"
    write(*,fmt='(5i16,3f16.8,1i16)') groundtruth_nspots,best_noof_spots,best_tp,best_fp,best_fn,best_f1,best_average_iou,best_recall,best_noof_spots

    !---    construct the output image, if needed
    if (len_trim(image_filename)/=0) then
        image_filename = trim( removeSuffix(image_filename) )//".multi_comparison.png"
        print *,"f1score info - constructing output image file """//trim(image_filename)//""""
        noof_valid_spots=spot_nSpots
        
        !---------------------
        print*,"DBG there are",groundtruth_nSpots,"GT spots"
        do ii = 1,groundtruth_nSpots
            D_matrixg = getDMatrix( groundtruth_dat(COL_D1,ii )*diamscaleg/2, groundtruth_dat(COL_D2,ii )*diamscaleg/2, groundtruth_dat(COL_THETA,ii)*(PI/180))

            if (best_indxg(ii)/=UNMATCHED_SPOT) then
                call drawEllipse(D_matrixg,groundtruth_dat(COL_X:COL_Y,ii ),rgb_img,colourscale=COLOURSCALE_REDBLUE,f=0.70d0,mode=LIB_DRAWELLIPSE_SHADE_RING)                    
            else
                call drawEllipse(D_matrixg,groundtruth_dat(COL_X:COL_Y,ii ),rgb_img,colourscale=COLOURSCALE_REDBLUE,f=0.5d0,mode=LIB_DRAWELLIPSE_SHADE_RING)                    
            end if

        end do
        
        do ii = 1,spot_nSpots
            D_matrix = getDMatrix( spot_dat(COL_D1,ii )*best_dscale/2, spot_dat(COL_D2,ii )*best_dscale/2, spot_dat(COL_THETA,ii)*(PI/180))
            if (spot_dat(COL_I,ii )>=best_ithresh) then
                if (best_indx(ii)/=UNMATCHED_SPOT) then
                    call drawEllipse(D_matrix,spot_dat(COL_X:COL_Y,ii ),rgb_img,colourscale=COLOURSCALE_RAINBOW,f=0.175d0,mode=LIB_DRAWELLIPSE_SHADE_RING)                    
                else
                    call drawEllipse(D_matrix,spot_dat(COL_X:COL_Y,ii ),rgb_img,colourscale=COLOURSCALE_RAINBOW,f=0.250d0,mode=LIB_DRAWELLIPSE_SHADE_RING)                    
                end if
            else
                cycle
            end if
        end do
        
        print *,"f1score info - writing output image file """//trim(image_filename)//""""
        print *,"f1score info - blue = ground truth (true positive), white = ground truth (false negative), red = matched, yellow = unmatched"
        
        call write_rgb_png( image_filename,rgb_img ) !also jph eddit
        print *,""
    end if 

    if (hasArgument(cla,"o")) then    
        print*,"Printing results to """//trim(res_outfile)//""""
        call output_multi_f1_results( trim(res_outfile)//".txt",results_indx,d_thresh_array,i_thresh_array,f1_array )

    end if 



contains
!---^^^^^^^^
    subroutine  get_f1_score(distance_matrix,iou_cutoff_in,f1_dat,f1_score_out,average_iou_out,indxg_io,indx_io,noof_spots_in_thresh)
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!       Gets F1 score for a distance matrix of i detected and j groundtruth spots
        real(kind=real64),dimension(:,:),intent(inout)      ::      distance_matrix
        real(kind=real64),intent(in)                        ::      iou_cutoff_in
        integer,intent(in)                                  ::      noof_spots_in_thresh !of the detected spots how many in the thrshold.
        integer,dimension(3),intent(out)                    ::      f1_dat !count of (true postives,false positives,false negatives)
        real(kind=real64),intent(out)                       ::      average_iou_out
        real(kind=real64),intent(out)                       ::      f1_score_out
        integer,dimension(:),intent(inout)                  ::      indxg_io,indx_io
        integer                                             ::      ii,jj,s_noof_gt_spots,s_noof_detected_spots,s_trial
        logical,dimension(:),allocatable                    ::      good
        integer,dimension(2)                                ::      pair
        integer                                             ::      noof_spots_discarded

        s_noof_gt_spots=size(distance_matrix,DIM=2)!get noof gt, detected spots from matrix dimensions
        s_noof_detected_spots=size(distance_matrix,DIM=1)

        average_iou_out=0.0d0

        indx = UNMATCHED_SPOT
        indxg = UNMATCHED_SPOT
        noof_spots_discarded=s_noof_detected_spots-noof_spots_in_thresh
        if(noof_spots_discarded>=s_noof_detected_spots)then
            f1_score_out=0.0d0
            f1_dat(1)=0
            f1_dat(2)=s_noof_detected_spots
            f1_dat(3)=s_noof_gt_spots
            return
        end if
        if (s_noof_detected_spots >= s_noof_gt_spots) then  
            allocate(good(s_noof_gt_spots))
            good=.false.
            do trial = 1,s_noof_gt_spots
                pair = minloc( distance_matrix )
                ii = pair(1)            !   comparison spot
                jj = pair(2)            !   ground truth spot                    
                
                !indx(pair(2)) = pair(1)
        !        if (iou_mode==0) then !do circle IoU
                good(jj) = distance_matrix(ii,jj) <= 1 - iou_cutoff_in
                

                if (good(jj)) then
                    indx_io(ii) = jj
                    indxg_io(jj) = ii
                    average_iou_out=average_iou_out+(1.0d0-distance_matrix(ii,jj))
                end if
                distance_matrix(:,jj) = huge(1.0)
                distance_matrix(ii,:) = huge(1.0)

            end do

            ! average_iou_out=average_iou_out/noof_spots_in_thresh
            ! f1_dat(1) = count(good)                                !   true positive counts = number of good matches
            ! f1_dat(2) = noof_spots_in_thresh-s_noof_gt_spots             !   false positive counts = number of spots with no match in ground truth
            ! f1_dat(3) = s_noof_gt_spots - f1_dat(1)                    !   false negative counts = number of ground truth spots without a match

        else
            !allocate(indx(spot_nSpots))
            allocate(good(s_noof_detected_spots))
            do trial = 1,s_noof_detected_spots
                pair = minloc( distance_matrix )
                ii = pair(1)            !   comparison spot
                jj = pair(2)            !   ground truth spot                    
                !indx(pair(1)) = pair(2)
                good(ii) = distance_matrix(ii,jj) <= 1 - iou_cutoff_in
                if (good(ii)) then
                    indx_io(ii) = jj
                    indxg_io(jj) = ii
                    average_iou_out=average_iou_out+(1.0d0-distance_matrix(ii,jj))
                end if
                distance_matrix(:,jj) = huge(1.0)
                distance_matrix(ii,:) = huge(1.0)

            end do

            ! average_iou_out=average_iou_out/s_noof_gt_spots
            ! f1_dat(1) = count(good)                                !   true positive counts = number of good matches
            ! f1_dat(2) = noof_spots_in_thresh - f1_dat(1)                           !   false positive counts = number of spots with no match in ground truth
            ! f1_dat(3) = s_noof_gt_spots - noof_spots_in_thresh           !   false negative counts = number of ground truth spots without a match

        end if

        !---    get the f1 parameters
        if (noof_spots_in_thresh>=s_noof_gt_spots) then
            !for more detected than GT spots
            average_iou_out=average_iou_out/noof_spots_in_thresh
            f1_dat(1) = count(good)                                !   true positive counts = number of good matches
            f1_dat(2) = noof_spots_in_thresh-s_noof_gt_spots             !   false positive counts = number of spots with no match in ground truth
            f1_dat(3) = s_noof_gt_spots - f1_dat(1)                    !   false negative counts = number of ground truth spots without a match
        else
            average_iou_out=average_iou_out/s_noof_gt_spots
            f1_dat(1) = count(good)                                !   true positive counts = number of good matches
            f1_dat(2) = noof_spots_in_thresh - f1_dat(1)                           !   false positive counts = number of spots with no match in ground truth
            f1_dat(3) = s_noof_gt_spots - noof_spots_in_thresh !s_noof_gt_spots - noof_spots_in_thresh           !   false negative counts = number of ground truth spots without a match
        end if
        !---    compute the f1 score
    
        f1_score_out = f1_dat(1)*2.0d0 / ( 2.0d0*f1_dat(1) + f1_dat(2) + f1_dat(3) )
        return
    end subroutine get_f1_score



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

    pure logical function is_spot_in_thresh(spot_dat,groundtruth_dat,q_thresh_in,i_thresh_in)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !       Enquires if a spot/ground truth pair meets the currednt defined thresholds of i and q.
            real(kind=real64),dimension(:),intent(in)       ::      spot_dat,groundtruth_dat
            real(kind=real64),intent(in)                    ::      q_thresh_in,i_thresh_in
            is_spot_in_thresh=.true.
            if(spot_dat(COL_I)<i_thresh_in) is_spot_in_thresh = .false.
            if(groundtruth_dat(COL_I)<i_thresh_in) is_spot_in_thresh = .false.

            if(spot_dat(COL_Q)>q_thresh_in) is_spot_in_thresh = .false.
            if(groundtruth_dat(COL_Q)>q_thresh_in) is_spot_in_thresh = .false.

            return
    end function is_spot_in_thresh




        
        
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
                isOverThresh=(datline(COL_I) >= i_thresh).and.(datline(COL_Q) <= q_thresh)
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
                isOverThresh=(datline(COL_I) >= i_thresh).and.(datline(COL_Q) <= q_thresh)
                !if (datline(COL_I) >= i_thresh) then
                if (isOverThresh) then
                    nOverThresh = nOverThresh + 1
                    dat(1:nCol,nOverThresh) = datline(1:nCol)
                end if

            end do
        close(unit=501)
        nSpots = nOverThresh


    end subroutine readSpotsFile

    subroutine output_multi_f1_results( filename,noof_res_lines,diam_scale_array_in,inten_array_in,f1_array_in )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      output as standard .spots file 
                
                character(len=*),intent(in)                     ::      filename        !   Filename in.
                integer,intent(in)                              ::      noof_res_lines              !   number of spots fitted
                real(kind=real64),dimension(:),intent(in)                    ::      diam_scale_array_in,inten_array_in,f1_array_in       !   diameter = nmPerPixel*(2 sigma)*diamscale
                

                integer                     ::      ii                      !   Indices of detected spots
                    
                
                
                ! open(unit=502,file=trim(filename),action="write")
                !     write (500,fmt='(a)') "# Clock "//trim(VERSION)
                !     write (500,fmt='(a)') "# "//trim(commentLine)
                !     write (500,fmt='(a,2f12.3,a)') "# ",ww,hh," micrograph extent ("//LENGTHUNIT//")"
    
                !     write (500,fmt='(a2,f13.3,a)') "# ",aa,"  micrograph area ("//LENGTHUNIT//") discounting dead pixels"
                !     write (500,fmt='(a2,i8,a)') "# ",ng,"  number of spots"
                !     write (500,fmt='(a,100a12)') "# pos x ("//LENGTHUNIT//")  pos y ("//LENGTHUNIT//") diam 1 ("//LENGTHUNIT//") diam 2 ("//LENGTHUNIT//")"," angle(deg)"," intensity"," t* value"," quality"
                    
                !     do ii = 1,getN( fitted_g2d )
                !         write(500,fmt='(100f12.3)')  getX(fitted_g2d,ii)*nmPerPixel                            &
                !                                          ,getY    (fitted_g2d,ii)*nmPerPixel                        &
                !                                          ,getSigma(fitted_g2d,ii,.true.)*nmPerPixel*(2*diamscale)   &            !   2 for convert radius to diameter.
                !                                          ,getSigma(fitted_g2d,ii,.false.)*nmPerPixel*(2*diamscale)  &            !   2 for convert radius to diameter.
                !                                          ,getTheta(fitted_g2d,ii) *180.0d0/PI                       &            !  convert radians to degrees
                !                                          ,getF    (fitted_g2d,ii)                                   &
                !                                          ,getT    (fitted_g2d,ii)                                   &
                !                                          ,getQ    (fitted_g2d,ii)                             
                !     end do
                ! close(unit=500)

                open(unit=502,file=trim(filename),action="write")

                    write (502,fmt='(10a20, 10a20, 10a20)') "#Diamter_threshold","intensity_threshold","F1_score"
                    
                    do ii = 1,noof_res_lines
                        write(502,fmt='(100f12.3)')  diam_scale_array_in(ii),inten_array_in(ii),f1_array_in(ii)                           
                    end do
                close(unit=502)
                
                    
                return
            end subroutine output_multi_f1_results


    
end program multif1score





