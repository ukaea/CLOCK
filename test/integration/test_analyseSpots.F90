program test_analyseSpots
    !---^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      A simple test program which
    !*      draws a single 2d gaussian
    !*      adds noise
    !*      fits using analyseSpots
    !*      returns the intersection over union
    !*      note: this is slightly different to test_Lib_Gaussian2d because the ROI is tested too.
    !*      Daniel Mason
    !*      (c) UKAEA 2025
    !*
    !*      run with eg
    !*          for i in `seq 1 10` ; do ./test/test_analyseSpots -s1 10 -s2 20 -sig 0.05 -q | grep "iou" ; done
    
    
        use analyseSpots
        use Lib_DrawEllipse
        use Lib_Png
        use Lib_RandomSeed        
        use Lib_CommandLineArguments
        use Lib_MaxLikelihoodFilter
        use Lib_Gaussian2d
        use Lib_MultipleGaussians2d
        use Lib_ColouredTerminal     
    
        use iso_fortran_env
#ifdef MPI
        use mpi_f08
#endif        
        implicit none


        character(len=8),parameter  ::      VERSION = "0.0.1"
        integer,parameter           ::      Nx = 256                    !   (0:Nx-1) size of the box. Gaussian will be at (Nx/2,Ny/2)
        real(kind=real64),parameter ::      PI = 3.141592654d0

        real(kind=real64),dimension(0:Nx-1,0:Nx-1)          ::      img,img_smooth
        !logical,dimension(0:Nx-1,0:Nx-1)                    ::      mask

        logical                     ::      test=.false., failed_to_detect_single_spot=.false.
        real(kind=real64)           ::      major_rad = Nx/16, minor_rad = Nx/16, theta           !   size of ellipse (radii) and angle
        real(kind=real64)           ::      background = 0.1d0 ,f0 = 0.5d0 ,sig = 0.0d0           !   background intensity, peak intensity (over bg), noise level
        real(kind=real64)           ::      lambda = 2.0d0  
        real(kind=real64)           ::      tolerance=1e-6                                        !   floating point tolerance for test
        integer                     ::      seed = 12345
        logical                     ::      quiet = .false. , verbose = .false.
        integer,dimension(8)        ::      seed_out
        type(CommandLineArguments)  ::      cla
        
        real(kind=real64),dimension(2)          ::      p_in,p_out,work
        real(kind=real64),dimension(2,2)        ::      D_in,D_out

        integer                     ::      ii,nn
        real(kind=real64)           ::      iou , bg_sigma , bg_bar
        type(Gaussian2d)            ::      g2d
        type(MultipleGaussians2d)   ::      fitted_g2d
    
        real(kind=real64),dimension(:),allocatable          ::      iou_all

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

        call setProgramDescription( cla, "test_analyseSpots" )
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
        call get( cla,"m",MAXPIX ,LIB_CLA_OPTIONAL,"          maximum fraction of foreground pixels - ( 0 for automatic )" )                                                        
        call get( cla,"R",ROIRMAX ,LIB_CLA_OPTIONAL,"          maximum roi padding (px)" )
        call get( cla,"b",F_THRESH ,LIB_CLA_OPTIONAL,"          intensity threshhold (f0/sigma) ( 0 for automatic )" )                      
        call get( cla,"d",D_THRESH ,LIB_CLA_OPTIONAL,"          intensity dip between separable maxima" )                                            
        call get( cla,"lambda",lambda ,LIB_CLA_OPTIONAL,"     max likelihood filtering lengthscale (0 for off)" )
        call get( cla,"test",test ,LIB_CLA_OPTIONAL,     "          Runs program as a test routine" )                                            

        ii = 2; call get( cla,"dbg",work,ii ,LIB_CLA_OPTIONAL,"        debug ROI in vicinty of point (x,y)",5 )                      
        if (hasArgument(cla,"dbg")) then
            DBG_ROIX = nint( work(1) )
            DBG_ROIY = nint( work(2) )
        end if    
        
    !---    enact command line args
        if (rank==0) then
            if (.not. hasArgument(cla,"seed")) then        
                call init_random_seed()        
                call get_random_seed(seed_out)
                seed = mod( seed_out(1),65535 )
            end if
            if (hasArgument(cla,"test")) then !set up test parameters
                lambda=0.0d0


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
        OPRANDC = .not. quiet
        if (rank==0) call report(cla)
        if (hasHelpArgument(cla)) stop
        if (.not. allRequiredArgumentsSet(cla)) stop
        call delete(cla)

        if (.not. quiet) then
            print *,"test_analyseSpots info - random seed ",seed
            print *,"test_analyseSpots info - s1,s2       ",major_rad , minor_rad
            print *,"test_analyseSpots info - theta       ",theta
        end if

    !---    draw the ellipse first without bg to determine cropped region
        if (rank==0) then
            img = 0
            LIB_DRAWELLIPSE_GAUSSMULT = 5       !   note I want a nice smooth ellipse drawn ...
            p_in = Nx*0.5d0
            D_in = getDMatrix( major_rad, minor_rad, theta )
            g2d = Gaussian2d_ctor( (/p_in(1),p_in(2),f0,D_in(1,1),D_in(2,1),D_in(2,2) /))
            call report(g2d)
            call delete(g2d)        !    no cheating!
                

        !---    draw ellipse again, 
            img = reshape( gaussianVariate( Nx*Nx ) , (/Nx,Nx/) )
            img = background + img*sig      
            call add( getDMatrix( major_rad, minor_rad, theta ) , (/Nx*0.5d0,Nx*0.5d0/) , img , mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN , f=f0)        

        !---    ... but fitting is over 2 sigma only 
            LIB_DRAWELLIPSE_GAUSSMULT = 2
        end if
#ifdef MPI        
        call MPI_BCAST(img,Nx*Nx,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
#endif

    !---    test draw .png. Note that LIB_DRAWELLIPSE_IGNORE << 1 so will render black
        if (.not. quiet) call writePng("test_analyseSpots_input_ellipse.png",img)        
            
        if (lambda==0) then
            img_smooth = img        !   no smoothing
        else
            call gaussianBlurImage( img, lambda,img_smooth )

            !call maxLikelihoodFilter( img,img_smooth, lambda )  
            if (.not. quiet) call writePng("test_analyseSpots_smoothed_ellipse.png",img_smooth)   
        end if
        

            

    !---    pass to fit
        if (rank==0) then
            E_THRESH = 1d6      

            if (OPRANDC) then
                call fitSpots(img_smooth,fitted_g2d,f_auto=(F_THRESH==0),maxpix_auto=(MAXPIX==0),roi_auto=(ROIRMAX==0),ss=bg_sigma,bb=bg_bar,ofile="test_analyseSpots.png")
            else
                call fitSpots(img_smooth,fitted_g2d,f_auto=(F_THRESH==0),maxpix_auto=(MAXPIX==0),roi_auto=(ROIRMAX==0),ss=bg_sigma,bb=bg_bar)
            end if

        !---    find the iou of each spot, report the best divided by number of spots nn ( ground truth is 1 )
            nn = getN( fitted_g2d )
            allocate( iou_all( nn ))
            do ii = 1,nn
                g2d = getG( fitted_g2d,ii )
                call report(g2d)
                p_out = (/ g2d%x0,g2d%y0 /)
                call findSigmaAndAngle( g2d,major_rad, minor_rad,theta )
                D_out = getDMatrix( major_rad, minor_rad, theta )
                iou_all(ii) = intersectionOverUnion(D_in,p_in,D_out,p_out,mode = LIB_DRAWELLIPSE_IOU_CIRCLE)                
            end do
            ii = maxloc(iou_all,dim=1)
            g2d = getG( fitted_g2d,ii )
            call findSigmaAndAngle( g2d,major_rad, minor_rad,theta )
            iou = iou_all(ii)
            print *,"iou = ",iou," spot count ",nn," radii ",major_rad, minor_rad

            if (hasArgument(cla,"test")) then
                if (nn==1) then
                    print*,"single spot detected sucessfully"
                else
                    print*,"Error: dected", nn,"spots when there should only be one."
                    failed_to_detect_single_spot=.true.
                end if

                print*,"IoU should be greater than 0.99999"
                if (((iou-1.0d0)<tolerance).and..not.failed_to_detect_single_spot) then
                    print *,colour(GREEN,"PASS")
                else
                    print *,colour(RED,"FAIL")
                end if

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
                if (present(message)) print *,"test_analyseSpots "//trim(message)
            end if
#ifdef MPI            
            call MPI_FINALIZE(ierror)
#endif          
            stop
        end subroutine errorExit
        
    


            subroutine gaussianBlurImage( f_in,t , f_out )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      produce a Gaussian blur of input image, width t
            real(kind=real64),dimension(0:,0:),intent(in)       ::      f_in    !Image array in
            real(kind=real64),dimension(0:,0:),intent(inout)    ::      f_out   !Blurred image array out
            real(kind=real64),intent(in)                        ::      t       !Gaussian width to blur
            
            
            
            real(kind=real64),dimension(:),allocatable      ::      kernel ,f_stripe,w_stripe   !Arrays to hold kernel and the image broken into strips 
                                                                                                !in the x and y direction (for parallelization).
                
            
            integer             ::      Nx,Ny               !Image size in X and Y dimension/pixels
            integer             ::      ix,iy,jj,kk         !Pixel X, pixel Y, strip and kernel indicies. 
            integer             ::      Nk                  !number of pixels in kernel
            real(kind=real64)   ::      i2s2,ww,wf,ws       !0.5*t^2, kernal weight, kernel result, strip weight
            
            Nx = size(f_in,dim=1)
            Ny = size(f_in,dim=2)
                
        !---    compute the unnormalised kernel
            Nk = max(5,ceiling( t*5 ))             !   range of pixels to search is +/- Nk
            allocate(kernel( 0:Nk ))
            i2s2 = 1/(2*t*t)             
            kernel(0) = 1.0d0   
            do ix = 1,Nk                
                kernel( ix )  = exp( -ix*ix*i2s2 )                                
            end do
            
                    
        !---    compute the output blurred image. First do x strips.
            allocate(f_stripe(0:Ny-1))
            allocate(w_stripe(0:Nx-1))
            
!$OMP PARALLEL  PRIVATE(ix,iy,wf,ws,jj,kk,ww,f_stripe,w_stripe ) SHARED( f_in,f_out,Nx,Ny,Nk,kernel)
!$OMP DO             
            do ix = 0,Nx-1
                wf = 0.0d0 ; ws = 0.0d0
                do jj = max(0,ix-Nk),min(Nx-1,ix+Nk)
                    kk = abs(jj-ix)           
                    ww = kernel( kk )                                         
                    wf = wf + ww*f_in(jj,0)
                    ws = ws + ww
                end do
                f_out(ix,0) = wf
                w_stripe(ix) = 1/ws
            end do
!$OMP END DO     
        
!$OMP DO  
            do iy = 1,Ny-1
                do ix = 0,Nx-1
                    wf = 0.0d0 ; ws = 0.0d0
                    do jj = max(0,ix-Nk),min(Nx-1,ix+Nk)
                        kk = abs(jj-ix)           
                        ww = kernel( kk )                                         
                        wf = wf + ww*f_in(jj,iy)
                        ws = ws + ww
                    end do
                    f_out(ix,iy) = wf                     
                end do
            end do
!$OMP END DO            
        
    
        !   at this point, f_out has blurring in the x-direction only. weight stores the kernel weighting from this op.            
            
                                
        !---    now do y strips
!$OMP DO         
            do ix = 0,Nx-1
            !   make a copy of this stripe
                f_stripe(0:Ny-1) = f_out(ix,0:Ny-1)
                do iy = 0,Ny-1
                    wf = 0.0d0 ; ws = 0.0d0
                    do jj = max(0,iy-Nk),min(Ny-1,iy+Nk)
                        kk = abs(jj-iy)                                  
                        ww = kernel( kk )
                        wf = wf + ww*f_stripe(jj)
                        ws = ws + ww
                    end do
                    f_out(ix,iy) = wf *w_stripe(ix) /(ws)                !   note: ws /= 0 because there is always at least one pixel contributing. Several really.
                                        
                end do
            end do
!$OMP END DO            
                    
!$OMP END PARALLEL       
        !---    now have a normalised Gaussian blur function
                        
            
            return
        end subroutine gaussianBlurImage
                            

    end program test_analyseSpots