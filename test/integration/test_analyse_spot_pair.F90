program test_analyse_spot_pair
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      A simple test program which
    !*      draws two 2d gaussians
    !*      adds noise
    !*      fits using analysespots/lib_multigaussian2d
    !*      returns the intersection over union
    !*      Daniel Mason & James heath
    !*      (c) UKAEA 2025
    !*

    
    
        use Lib_Gaussian2d
        use Lib_MultipleGaussians2d
        use Lib_DrawEllipse
        use Lib_Png
        use Lib_RandomSeed        
        use Lib_CommandLineArguments
        use Lib_MaxLikelihoodFilter
        use Lib_ColouredTerminal
        use analyseSpots
        
        use iso_fortran_env
#ifdef MPI
        use mpi_f08
#endif        
        implicit none


        character(len=8),parameter  ::      VERSION = "0.1.1"
        integer,parameter           ::      Ny = 256                    !   (0:Nx-1) size of the box. Gaussian will be at (Nx/2,Ny/2)
        integer,parameter           ::      Nx = 256 
        real(kind=real64),parameter ::      PI = 3.141592654d0
        real(kind=real64),parameter         ::      TWO_PI = 3.141592653590d0*2.0d0
        character(len=256)              ::      outfile=""   !output filename 

        real(kind=real64),dimension(0:Nx-1,0:Nx-1)          ::      img,img_smooth
        logical,dimension(0:Nx-1,0:Nx-1)                    ::      mask

        real(kind=real64)           ::      major_rad = Nx/16, minor_rad = Nx/16, theta           !   size of ellipse (radii) and angle
        real(kind=real64)           ::      background = 0.1d0 ,f0 = 0.25d0 ,sig = 0.0d0           !   background intensity, peak intensity (over bg), noise level
        integer                     ::      seed = 12345
        logical                     ::      quiet = .false. , verbose = .false., test=.false., disable_filter=.false.
        integer,dimension(8)        ::      seed_out
        type(CommandLineArguments)  ::      cla
        
        real(kind=real64),dimension(2)          ::      p_in1,p_out1,p_in2,p_out2
        real(kind=real64),dimension(2,2)        ::      D_in1,D_in2

        integer                     ::      ii
        type(Gaussian2d)            ::      g2d1,g2d2,detected_g2d,detected_g2d1,detected_g2d2

        real(kind=real64)           ::      spot_speparation=0.0d0        !if this is set with CLA -d to be greater than 0, two spots will be drawn spot_speration/2 pixels away from the image centre
        real(kind=real64)           ::      spot_radii_ratio=1.0d0      !How many times larger is the right spot than the left spot, ignored if -d = 0
        real(kind=real64)           ::      spot_intensity_ratio=1.0d0  !How many times more bright is the right spot than the left spot, ignored if -d = 0 
        type(MultipleGaussians2d)                       ::      fitted_g2d  !object to hold parameters of several gaussians.
        real(kind=real64)           ::      mean_iou
        logical                     ::      failed_to_detect_2_spots=.false.   !marks if CLOCK has detected more  or less than two spots, hence has broken down.
        real(kind=real64),dimension(:),allocatable  ::  squared_distances1,squared_distances2 !squared idstances between detected points and groundtruth left and right points respectively.
        integer                     ::      best_indx_1,best_indx_2        !Indices of the detected ponits  closest to the leftmost and rigtmost GT points respectively
        integer                     ::      noof_spots_detected
        real(kind=real64)           ::      sin_wave_amplitude = 0.0d0
        integer                     ::      iy
        real(kind=real64)           ::     bg_sigma,bg_bar

    

    !---    parallelization
        integer                 ::      rank,nprocs 
#ifdef MPI
        integer                 ::      ierror
#endif        


#ifdef MPI
        call MPI_INIT(ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
#else
        rank = 0
        nprocs = 1
#endif        
    
    !---    read command line arguments
        cla = CommandLineArguments_ctor(20)  

        call setProgramDescription( cla, "test_analyse_spot_pair" )
        call setProgramVersion( cla, VERSION )   
            
        call get( cla,"s1",major_rad ,LIB_CLA_OPTIONAL,  "          major radius" )       
        call get( cla,"s2",minor_rad ,LIB_CLA_OPTIONAL,  "          minor radius" )       
        
        call get( cla,"theta",theta ,LIB_CLA_OPTIONAL,   "       angle (rad)" )        
        call get( cla,"bg",background ,LIB_CLA_OPTIONAL, "          background intensity" )        
        call get( cla,"f0",f0 ,LIB_CLA_OPTIONAL,         "          peak intensity" )        
        call get( cla,"sig",sig ,LIB_CLA_OPTIONAL,       "         noise std dev" )        
        call get( cla,"seed",seed ,LIB_CLA_OPTIONAL,     "        random seed" )        
        call get( cla,"q",quiet ,LIB_CLA_OPTIONAL,       "           quiet mode- suppresses unnecessary output" )        
        call get( cla,"verbose",verbose ,LIB_CLA_OPTIONAL,       "     verbose mode" )
        call get( cla,"o",outfile ,LIB_CLA_OPTIONAL,     "          output filename, note this is the noisy image for imageJ comparsion" )
        call get( cla,"test",test ,LIB_CLA_OPTIONAL,     "          Runs program as a test routine" )
        call get( cla,"d",spot_speparation ,LIB_CLA_OPTIONAL,     "          Draws two spots spearted by px value set here." )
        call get( cla,"r",spot_radii_ratio ,LIB_CLA_OPTIONAL,     "          If drawing two spots, how many times larger should the right one be?" )
        call get( cla,"i",spot_intensity_ratio ,LIB_CLA_OPTIONAL,     "           If drawing two spots, how many times brighter should the right one be?" )
        call get( cla,"df",disable_filter ,LIB_CLA_OPTIONAL,     "          Turns off the filtering" )
        call get( cla,"sinb",sin_wave_amplitude ,LIB_CLA_OPTIONAL,"   Sets the amplitude of a sin wave added to th image from top to bottom" )



        
    if (hasArgument(cla,"test")) then
        test=.true.
        major_rad=10.0d0
        minor_rad=10.0d0
        sig=0.0d0!0.18d0
        theta=0.0d0
        spot_speparation=30.0d0
    end if
    
    !---    enact command line args
        if (rank==0) then
            if (.not. hasArgument(cla,"seed")) then        
                call init_random_seed()        
                call get_random_seed(seed_out)
                seed = mod( seed_out(1),65535 )
            end if
        end if
#ifdef MPI
        call MPI_BCAST(seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
#endif
        call init_random_seed(seed)   
        if (rank==0) then   
            if (.not. hasArgument(cla,"theta")) then
                do ii = 1,10
                    theta = ran0(seed) 
                end do
                theta = PI * theta
            end if
        end if
#ifdef MPI
        call MPI_BCAST(theta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
#endif
                    
        quiet = quiet .or. (rank/=0)

        if (rank==0) call report(cla)
        if (hasHelpArgument(cla)) stop
        if (.not. allRequiredArgumentsSet(cla)) stop
        call delete(cla)

        if (.not. quiet) then
            print *,"test_analyse_spot_pair info - random seed ",seed
            print *,"test_analyse_spot_pair info - s1,s2       ",major_rad , minor_rad
            print *,"test_analyse_spot_pair info - theta       ",theta
        end if

    !---    draw the ellipse first without bg to determine cropped region
        if (rank==0) then
            img = 0
            LIB_DRAWELLIPSE_GAUSSMULT = 5       !   note I want a nice smooth ellipse drawn ...




            print*,"Doing double spots"
            ! do left hand spot (1)
            p_in1 = Nx*0.5d0
            p_in1(1)=p_in1(1)-spot_speparation/2.0d0
            D_in1 = getDMatrix( major_rad, minor_rad, theta )
            g2d1 = Gaussian2d_ctor( (/p_in1(1),p_in1(2),f0,D_in1(1,1),D_in1(2,1),D_in1(2,2) /))
            call report(g2d1)
            call delete(g2d1)        !    no cheating!
            call drawEllipse( getDMatrix( major_rad, minor_rad, theta ) , p_in1 , img , mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN , f=f0)      
            where (img >=  f0 * exp(-4.5d0) )       !   3 sigma
                mask = .true.
            elsewhere
                mask = .false.
            endwhere


            
            !do right spot, this has more potentail modifiers (2)
            p_in2 = Nx*0.5d0
            p_in2(1)=p_in2(1)+spot_speparation/2.0d0
            D_in2 = getDMatrix( major_rad*spot_radii_ratio, minor_rad*spot_radii_ratio, theta )
            g2d2 = Gaussian2d_ctor( (/p_in2(1),p_in2(2),f0,D_in2(1,1),D_in2(2,1),D_in2(2,2) /))
            call report(g2d2)
            call delete(g2d2)        !    no cheating!
            call drawEllipse( getDMatrix( major_rad*spot_radii_ratio, minor_rad*spot_radii_ratio, theta ) , p_in2 , img , mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN , f=f0*spot_intensity_ratio)      
            where (img >=  f0 * exp(-4.5d0) )       !   3 sigma
                mask = .true.
            elsewhere
                mask = .false.
            endwhere

            !---    draw ellipses again, 
            img = reshape( gaussianVariate( Nx*Nx ) , (/Nx,Nx/) )

            !---    Add a sin wave to the image intensity from top to bottom, use if long wavelength variation is required.
            !---    Check if added sinwave amplitude is sensible
            if ((abs(sin_wave_amplitude)+background)>1.0d0)then
                print*,"Sum of sin wave amplitude (-sinb) and background intensity (-b) must be less than 1."
                print*,"Setting sin wave amplitude to zero."
                sin_wave_amplitude=0.0d0
            end if
            
            if (abs(sin_wave_amplitude)>background) then
                print*,"Sin wave amplitude (-sinb) needs to be less or equal to the background intensity (-b) to avoid negative intnesity."
                print*,"Setting Sin wave amplitude to background intensity:", background
                sin_wave_amplitude=background
            end if

            if (abs(sin_wave_amplitude)>0.0d0) then
                do iy =0,Ny-1
                    img(:,iy)=img(:,iy)*sig+background+sin_wave_amplitude*sin((PI/2.0d0)+TWO_PI*((iy*1.0d0)/(Ny-1)))
                end do
            else
                img = background + img*sig
            end if

            !img = background + img*sig      
            call add( D_in1 , (/p_in1(1),p_in1(2)/) , img , mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN , f=f0)
            call add( D_in2 , (/p_in2(1),p_in2(2)/) , img , mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN , f=f0*spot_intensity_ratio)

        !---    ... but fitting is over 2 sigma only 
            LIB_DRAWELLIPSE_GAUSSMULT = 2

            end if 

#ifdef MPI        
        call MPI_BCAST(img,Nx*Nx,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
#endif

    !---    test draw .png. Note that LIB_DRAWELLIPSE_IGNORE << 1 so will render black
        if (.not. quiet) call writePng("test_analyse_spot_pair_input_ellipse.png",img) 

        if(.not.(hasArgument(cla,"df"))) then
            print*,"Denoising with maxliklihood filter."
            call maxLikelihoodFilter( img,img_smooth, 2.0d0 )
        else 
            print*,"No noise, skipping denoising."
            img_smooth=img
        end if
        
        
        if (.not. quiet) call writePng("test_analyse_spot_pair_smoothed_ellipse.png",img_smooth)   
            

    !---    pass to fit !!! and print
        if (rank==0) then
            mean_iou=0.0d0
            E_THRESH = 6.0d0!1d6  
            print*,"doing second part of double spot"
            call fitSpots(img_smooth,fitted_g2d,f_auto=(F_THRESH==0),maxpix_auto=(MAXPIX==0),roi_auto=(ROIRMAX==0),ss=bg_sigma,bb=bg_bar)
            !get noof spots detected
            noof_spots_detected=getN(fitted_g2d,ignore=.true.)
            print*,"noof_spots_detected is:", noof_spots_detected
            if (noof_spots_detected/=2) failed_to_detect_2_spots=.true.
            allocate(squared_distances1(noof_spots_detected))
            allocate(squared_distances2(noof_spots_detected))
            do ii=1,noof_spots_detected
                !get g2d
                detected_g2d=getg(fitted_g2d,ii)
                !print*,ii,"Detected spot pos at:",detected_g2d%x0,", ",detected_g2d%y0
                !get squared distances
                ! squared_distances1(ii)=((p_in1(1)-detected_g2d%x0)*(p_in1(1)-detected_g2d%x0))+((p_in1(2)-detected_g2d%y0)*(p_in1(2)-detected_g2d%y0))
                ! squared_distances2(ii)=((p_in2(1)-detected_g2d%x0)*(p_in2(1)-detected_g2d%x0))+((p_in2(2)-detected_g2d%y0)*(p_in2(2)-detected_g2d%y0))

                p_out1 = (/ detected_g2d%x0,detected_g2d%y0/) !try 1-Iou for "distance" metric, avoids the case of a bad spot with a perfect centroid.
                p_out2 = (/ detected_g2d%x0,detected_g2d%y0/)
                squared_distances1(ii)=1.0d0-intersectionOverUnion(D_in1,p_in1,getD(fitted_g2d,ii),p_out1,mode = LIB_DRAWELLIPSE_IOU_ELLIPSE)
                squared_distances2(ii)=1.0d0-intersectionOverUnion(D_in2,p_in2,getD(fitted_g2d,ii),p_out2,mode = LIB_DRAWELLIPSE_IOU_ELLIPSE)
            end do
            best_indx_1=minloc(squared_distances1,1) !detected with smallest square distance to point n is most likely detection of point n
            best_indx_2=minloc(squared_distances2,1)

            if(noof_spots_detected==0)then
                !no spots
                mean_iou=0.0d0
            else
                detected_g2d1=getg(fitted_g2d,best_indx_1)
                p_out1 = (/ detected_g2d1%x0,detected_g2d1%y0/)
                detected_g2d2=getg(fitted_g2d,best_indx_2)
                p_out2 = (/ detected_g2d2%x0,detected_g2d2%y0/)

                if (noof_spots_detected==1) then
                    if (squared_distances1(best_indx_1)>=squared_distances2(best_indx_2)) then
                        mean_iou=intersectionOverUnion(D_in2,p_in2,getD(fitted_g2d,best_indx_2),p_out2,mode = LIB_DRAWELLIPSE_IOU_ELLIPSE)!LIB_DRAWELLIPSE_IOU_CIRCLE
                        mean_iou=mean_iou/2.0d0
                    else 
                        mean_iou=intersectionOverUnion(D_in1,p_in1,getD(fitted_g2d,best_indx_1),p_out1,mode = LIB_DRAWELLIPSE_IOU_ELLIPSE)
                        mean_iou=mean_iou/2.0d0
                    end if
                else
                        mean_iou=intersectionOverUnion(D_in1,p_in1,getD(fitted_g2d,best_indx_1),p_out1,mode = LIB_DRAWELLIPSE_IOU_ELLIPSE)
                        mean_iou=mean_iou+intersectionOverUnion(D_in2,p_in2,getD(fitted_g2d,best_indx_2),p_out2,mode = LIB_DRAWELLIPSE_IOU_ELLIPSE)
                        mean_iou=mean_iou/2.0d0
                end if
            end if

            if(noof_spots_detected==0)then
                print*,"mean IoU is:",0.0d0
                print*,"Radii (maj1, min1, maj2, min2, mean1, mean2, mean12)",0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0
            else
            print*,"mean IoU is:",mean_iou
                print*,"Radii (maj1, min1, maj2, min2, mean1, mean2, mean12)",getsigma(fitted_g2d,best_indx_1,.true.),getsigma(fitted_g2d,best_indx_1,.false.),&
                    getsigma(fitted_g2d,best_indx_1),getsigma(fitted_g2d,best_indx_2,.true.),getsigma(fitted_g2d,best_indx_2,.false.),getsigma(fitted_g2d,best_indx_2),&
                    ((getsigma(fitted_g2d,best_indx_2)+getsigma(fitted_g2d,best_indx_1))/2.0d0)
            end if


            if(test) then
                print*,"Mean Iou should be greater than 0.85"
                if ((mean_iou>0.85d0).and..not.failed_to_detect_2_spots) then
                    print *,colour(GREEN,"PASS")
                else
                    print *,colour(RED,"FAIL")
                end if
            end if

            if (hasArgument(cla,"o")) then
                call write_greyscale_png(outfile,img)
            end if

        end if






                


        if (.not. quiet) print *,""  
        call errorExit()

    contains
!---^^^^^^^^

    

        subroutine errorExit(message)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            character(len=*),intent(in),optional             ::      message
#ifdef MPI
        integer                 ::      ierror
#endif        

            if (rank==0) then
                if (present(message)) print *,"test_analyse_spot_pair "//trim(message)
            end if
#ifdef MPI            
            call MPI_FINALIZE(ierror)
#endif          
            stop
        end subroutine errorExit
        
    


    end program test_analyse_spot_pair