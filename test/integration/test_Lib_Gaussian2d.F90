program test_Lib_Gaussian2d
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      A simple test program which
    !*      draws a single 2d gaussian
    !*      adds noise
    !*      fits using Lib_Gaussian2d
    !*      returns the intersection over union
    !*      Daniel Mason
    !*      (c) UKAEA 2025
    !*
    !*      run with eg
    !*          for i in `seq 1 10` ; do ./test/test_Lib_Gaussian2d -s1 10 -s2 20 -sig 0.05 -q | grep "iou" ; done
    
    
        use Lib_Gaussian2d
        use Lib_DrawEllipse
        use Lib_Png
        use Lib_RandomSeed        
        use Lib_CommandLineArguments
        use Lib_MaxLikelihoodFilter
        use Lib_ColouredTerminal
        
        use iso_fortran_env
#ifdef MPI
        use mpi_f08
#endif        
        implicit none


        character(len=8),parameter  ::      VERSION = "0.1.1"
        integer,parameter           ::      Nx = 256                    !   (0:Nx-1) size of the box. Gaussian will be at (Nx/2,Ny/2)
        real(kind=real64),parameter ::      PI = 3.141592654d0
        character(len=256)              ::      outfile   !output filename 

        real(kind=real64),dimension(0:Nx-1,0:Nx-1)          ::      img,img_smooth
        logical,dimension(0:Nx-1,0:Nx-1)                    ::      mask

        real(kind=real64)           ::      major_rad = Nx/16, minor_rad = Nx/16, theta           !   size of ellipse (radii) and angle
        real(kind=real64)           ::      background = 0.1d0 ,f0 = 0.5d0 ,sig = 0.0d0           !   background intensity, peak intensity (over bg), noise level
        integer                     ::      seed = 12345
        logical                     ::      quiet = .false. , verbose = .false., test=.false.
        integer,dimension(8)        ::      seed_out
        type(CommandLineArguments)  ::      cla
        
        real(kind=real64),dimension(2)          ::      p_in,p_out
        real(kind=real64),dimension(2,2)        ::      D_in,D_out

        integer                     ::      ii
        real(kind=real64)           ::      iou
        logical                     ::      ok
        type(Gaussian2d)            ::      g2d

    

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
        cla = CommandLineArguments_ctor(15)  

        call setProgramDescription( cla, "test_Lib_Gaussian2d" )
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
        call get( cla,"o",outfile ,LIB_CLA_OPTIONAL,     "          output filename, note this is the noisy image for imageJ comparsion." )
        call get( cla,"test",test ,LIB_CLA_OPTIONAL,     "          Runs program as a test routine" )


        
    if (hasArgument(cla,"test")) then
        test=.true.
        major_rad=32.0d0
        minor_rad=10.0d0
        sig=0.18d0
        theta=0.0d0
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
            print *,"test_Lib_Gaussian2d info - random seed ",seed
            print *,"test_Lib_Gaussian2d info - s1,s2       ",major_rad , minor_rad
            print *,"test_Lib_Gaussian2d info - theta       ",theta
        end if

    !---    draw the ellipse first without bg to determine cropped region
        if (rank==0) then
            img = 0
            !LIB_DRAWELLIPSE_RSS_REMOVE_BG = .true.
            LIB_DRAWELLIPSE_GAUSSMULT = 5       !   note I want a nice smooth ellipse drawn ...
            p_in = Nx*0.5d0
            D_in = getDMatrix( major_rad, minor_rad, theta )
            g2d = Gaussian2d_ctor( (/p_in(1),p_in(2),f0,D_in(1,1),D_in(2,1),D_in(2,2) /))
            call report(g2d)
            call delete(g2d)        !    no cheating!
            call drawEllipse( getDMatrix( major_rad, minor_rad, theta ) , p_in , img , mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN , f=f0)      
            where (img >=  f0 * exp(-4.5d0) )       !   3 sigma
                mask = .true.
            elsewhere
                mask = .false.
            endwhere

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
        if (.not. quiet) call writePng("test_Lib_Gaussian2d_input_ellipse.png",img)        
            
        call maxLikelihoodFilter( img,img_smooth, 2.0d0 )
        
        if (.not. quiet) call writePng("test_Lib_Gaussian2d_smoothed_ellipse.png",img_smooth)   
            

    !---    pass to fit !!! and print
        if (rank==0) then
            E_THRESH = 1d6      
            call fit( g2d,img_smooth,imin = background,imax = 1.0d0,mask = mask,ok = ok, dbg = verbose)
            if (ok) then
                    
                call report(g2d)
                p_out = (/ g2d%x0,g2d%y0 /)
                call findSigmaAndAngle( g2d,major_rad, minor_rad,theta )
                D_out = getDMatrix( major_rad, minor_rad, theta )
                
                iou = intersectionOverUnion(D_in,p_in,D_out,p_out,mode = LIB_DRAWELLIPSE_IOU_CIRCLE)
                
            else
                print *,"no reasonable gaussian found"
                iou = 0.0d0
            end if
            print *,"iou = ",iou
            if (hasArgument(cla,"o")) then
                call write_greyscale_png(outfile,img)
            end if
            if(test) then
                print*,"Expected Iou should be greater than 0.9"
                if (iou>0.9d0) then
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
                if (present(message)) print *,"test_Lib_Gaussian2d "//trim(message)
            end if
#ifdef MPI            
            call MPI_FINALIZE(ierror)
#endif          
            stop
        end subroutine errorExit
        
    


    end program test_Lib_Gaussian2d