!
!   A simple standalone code to remove noise from an image
!   The program operates in  steps
!       1)  read in an image from a png file
!       2)  compute a smoothed version
!
!   Part of the CLOCK toolkit for automatic image characterization
!
!   Daniel Mason
!   (c) UKAEA Jan 2025
!
!   Version history
!       0.0.1       Feb 2025        First working version
!       0.0.2       Jul 2025        Moved filter into its own library
!
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    MLEFilter from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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


program MLEFilter
    !---^^^^^^^^^^^^^^^^^^ 
    !       
        use iso_fortran_env
        use Lib_Png
        use Lib_Filenames
        use Lib_LaplacianKernel        
        use Lib_CommandLineArguments
        use Lib_MaxLikelihoodFilter
        use Lib_ExtrapolateImage
#ifdef MPI
        use mpi_f08
#endif
    
        implicit none

    !---    magic numbers fixed in the code 
        character(len=*),parameter      ::      VERSION = "0.0.1"

    !---    information about the run from command line params
        type(CommandLineArguments)          ::      cla       
        character(len=256)                  ::      filename = "test"          !   input filename
        character(len=256)                  ::      outfile = ""
        real(kind=real64)                   ::      lambda = 2.0d0         !    lengthscale
        integer                             ::      kernelSize = 5          !   which kernel to use

    !---    information about the image
        integer                                         ::      Nx,Ny               !   image size in pixels
        real(kind=real64),dimension(:,:),allocatable    ::      img_in              !   (1:Nx,1:Ny)     input image 
        real(kind=real64),dimension(:,:),allocatable    ::      img_out             !   output image


    !---    dummy variables
        logical                 ::      ok
        real(kind=real64),parameter                     ::      X0 = 1.0d0 , Y0 = 1.0d0 
        real(kind=real64)       ::      stdev,curv,dd,xx,yy,zz
        real(kind=real64),parameter         ::  DX = 0.25d0

        integer                 ::      ii,jj
        real(kind=real64),dimension(5)      ::  calc_val = 0.0d0
        real(kind=real64),dimension(5)      ::  target_val = 0.0d0
    
    !---    parallelization
        integer                 ::      rank,nprocs,ierror

#ifdef MPI
        call MPI_INIT(ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
#else
        rank = 0
        nprocs = 1
        ierror = 0
#endif


!---    read command line arguments
        cla = CommandLineArguments_ctor(30)  
                
        call setProgramDescription( cla, "MLEFilter" )
        call setProgramVersion( cla, VERSION ) 
        

        !---    filename options        
        call get( cla,"f",filename ,LIB_CLA_REQUIRED,"           image filename (use test for test program)" )
        outfile = trim( removeSuffix(filename) )//".mle.png"
        call get( cla,"o",outfile ,LIB_CLA_OPTIONAL,"      output image filename" )                      
        call get( cla,"lambda",lambda ,LIB_CLA_OPTIONAL," characteristic smoothing lengthscale" )    
        call get( cla,"k",kernelSize ,LIB_CLA_OPTIONAL,"      which kernel to use for biharmonic (5 or 7)" )                    


        if (rank==0) call report(cla)
        if (hasHelpArgument(cla)) call errorExit()
        if (.not. allRequiredArgumentsSet(cla)) call errorExit()
        call delete(cla)
    

    !---    read in the input image
        if (trim(filename)=="test") then
            print *,"test program"
            Nx = 13; Ny = 13
            allocate(img_in(-Nx/2:Nx/2,-Ny/2:Ny/2))
            print *,"function"
            do jj = -Ny/2,Ny/2
                yy = Y0 + jj*DX
                do ii = -Nx/2,Nx/2
                    xx = X0 + ii*DX
                    dd = func(xx,yy)
                    img_in(ii,jj) = dd
                    !print *,ii,jj,xx,yy,func(xx,yy)
                    write(*,fmt='(f16.8)',advance="no") dd
                end do
                print *,""
            end do
            target_val(1) = 0.0d0
            print *,"laplacian"
            do jj = -Ny/2,Ny/2
                yy = Y0 + jj*DX
                do ii = -Nx/2,Nx/2
                    xx = X0 + ii*DX
                    dd = laplacianfunc(xx,yy)
                    target_val(1) = target_val(1) + dd*dd
                    write(*,fmt='(f16.8)',advance="no") dd
                end do
                print *,""
            end do
            target_val(1) = target_val(1) / (Nx*Ny)
            dd = laplacianfunc(X0,Y0)
            target_val(2) = dd
            calc_val(2) = discreteLaplacian( img_in(-2:2,-2:2),noisy=.false. )/(DX*DX)
            target_val(3) = dd*dd
            calc_val(3) = discreteLaplacian_squared( img_in(-2:2,-2:2),noisy=.false. )/(DX*DX*DX*DX)
            calc_val(1) = laplacianSquared( img_in )
            print *,"target < ( del^2(g) ) ^2 > ",target_val(1),calc_val(1), " note: this is sensitive to boundary conditions, and so likey to be inaccurate"
            print *,"target ( del^2 f (0,0) )     ",target_val(2),calc_val(2) 
            print *,"target ( del^2 f (0,0) ) ^2  ",target_val(3),calc_val(3)
            call writePng( "test.png",img_in )
            

            dd = 1.0d-4
            img_in(0,0) = img_in(0,0) + dd
            xx = laplacianSquared( img_in ) * (Nx*Ny)
            img_in(0,0) = img_in(0,0) - 2*dd
            yy = laplacianSquared( img_in ) * (Nx*Ny)
            img_in(0,0) = img_in(0,0) + dd
            zz = laplacianSquared( img_in ) * (Nx*Ny)
            target_val(4) = (xx-yy) / (DX*DX*DX*DX * 2*dd )
            calc_val(4) = 2 * discreteBiharmonic( img_in(-4:4,-4:4) ) / (DX*DX*DX*DX)
            print *,"deriv ",target_val(4),calc_val(4)

            target_val(5) = (xx-2*zz+yy) / (DX*DX*DX*DX * dd*dd )
            print *,"2nd deriv ",target_val(5)

        else
            if (rank==0) then
                inquire(file=trim(filename),exist=ok)
                if (.not. ok) call errorExit( "error - file not found """//trim(filename)//"""" )
                call readPng( filename,img_in )
            end if
#ifdef MPI
            Nx = size(img_in,dim=1)
            Ny = size(img_in,dim=2)
            call MPI_BCAST(Nx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
            call MPI_BCAST(Ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
            if (rank/=0) allocate(img_in(Nx,Ny))
            call MPI_BCAST(img_in,Nx*Ny,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
#endif
        end if
        Nx = size(img_in,dim=1)
        Ny = size(img_in,dim=2)
        allocate(img_out(Nx,Ny))
        if (rank==0) print *,"MLEFilter info - read in image with ",Nx,",",Ny," px"
        if (Nx*Ny<=0) call errorExit( "error - incorrect read?? pixel extent error")

    

        call maxLikelihoodFilter( img_in,img_out, lambda, kernelSize )
            
        if (rank==0) then
            call errorTerms( img_in,img_out , stdev,curv )
            print *,"MLEFilter info - after filter stdev = sqrt( <f-g>^2 ) = ",stdev
            print *,"                           curv = < ( del^2(g) ) ^2 > = ",curv
        end if

    !---    output the result
        if ((rank==0) .and. (trim(filename)/="test")) then
            print *,"write to """//trim(outfile)//""""
            call writePng( outfile,img_out )
        end if
    
    !---    bye bye
        call errorExit("done")  
    
        
    contains
!---^^^^^^^^

        subroutine errorTerms( f,g , stdev,curv )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      given two images f and g
    !*      compute stdev = sqrt( <f-g>^2 ) - a measure of the distance between the pixel intensities
    !*      and curv = < ( del^2(g) ) ^2 >  - a measure of the curvature of the second image
            real(kind=real64),dimension(:,:),intent(in)     ::      f               !   observed image
            real(kind=real64),dimension(:,:),intent(in)     ::      g               !   ground truth
            real(kind=real64),intent(out)                   ::      stdev,curv

    
            integer         ::      Nx,Ny
            integer         ::      ix,iy,jx,jy,kx,ky 
            integer         ::      mm
            real(kind=real64),dimension(:,:),allocatable    ::      kernel
            real(kind=real64)       ::      laplacian,dfg,sumdfg2,sumdel22
            
            
        !---    find the correct kernel      
            mm = 2                              !   kernel half-width            
            allocate(kernel(-mm:mm,-mm:mm))
            kernel = getKernel5x5(noisy=.true.)
                

            
        !---    establish size of problem and allocate memory
            Nx = size(f,dim=1)
            Ny = size(f,dim=2)
            sumdfg2 = 0
            sumdel22 = 0
            
            do iy = 1,Ny
                do ix = 1,Nx        
                    
                    dfg = f(ix,iy) - g(ix,iy)
                    sumdfg2 = sumdfg2 + dfg*dfg

                    laplacian = 0.0d0
                    do ky = -mm,mm
                        jy = ky + iy
                        if ( (jy<1).or.(jy>Ny) ) jy = iy - ky
                                                
                        do kx = -mm,mm
                            jx = kx + ix
                            if ( (jx<1).or.(jx>Nx) ) jx = ix - kx

                            laplacian = laplacian + kernel( kx,ky ) * f(jx,jy)

                        end do
                    end do
                    sumdel22 = sumdel22 + laplacian*laplacian

                end do
            end do
    
            stdev = sqrt( max(0.0d0,sumdfg2/(Nx*Ny)) )
            curv = ( max(0.0d0,sumdel22/(Nx*Ny)) )
            
            return
        end subroutine errorTerms     


        !pure 
        real(kind=real64) function laplacianSquared( g )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given two images f and g
    !*      compute < ( del^2 g ) ^2 >  - a measure of the curvature of the second image
            real(kind=real64),dimension(:,:),intent(in)     ::      g               !   ground truth
                

    
            integer         ::      Nx,Ny
            integer         ::      ix,iy,jx,jy,kx,ky 
            integer         ::      mm
            real(kind=real64),dimension(:,:),allocatable    ::      kernel
            real(kind=real64)       ::      laplacian
            real(kind=real64),dimension(:,:),allocatable    ::      g_border
            
        !---    find the correct kernel      
            mm = 2                              !   kernel half-width            
            allocate(kernel(-mm:mm,-mm:mm))
            kernel = getKernel5x5()
                

            
        !---    establish size of problem and allocate memory
            Nx = size(g,dim=1)
            Ny = size(g,dim=2)
            laplacianSquared = 0
            allocate(g_border(1-mm:Nx+mm,1-mm:Ny+mm))
            call extrapolate( g,mm,g_border,FADE_STYLE_AVG )

            do iy = 1,Ny
                do ix = 1,Nx
                        

                    laplacian = 0.0d0
                    do ky = -mm,mm
                        jy = ky + iy
                                                
                        do kx = -mm,mm
                            jx = kx + ix

                            laplacian = laplacian + kernel( kx,ky ) * g_border(jx,jy)

                        end do
                    end do
                    laplacianSquared = laplacianSquared + laplacian*laplacian

                end do
            end do
    
            laplacianSquared = laplacianSquared/(Nx*Ny)
            
            return
        end function laplacianSquared     


    
                        
            
        pure real(kind=real64) function func(x,y)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),intent(in)        ::      x,y
            real(kind=real64),parameter         ::      A = 0.1d0 , B = 0.2d0 
            func = cos( A*x*x ) * sin( B*y ) 
            return
        end function func


        pure real(kind=real64) function laplacianfunc(x,y)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),intent(in)        ::      x,y
            real(kind=real64),parameter         ::      A = 0.1d0 , B = 0.2d0 
            laplacianfunc = (- 4*A*A*x*x*cos( A*x*x ) - 2*A*sin( A*x*x ))* sin( B*y )  - B*B*cos( A*x*x ) * sin( B*y )  
            return
        end function laplacianfunc

        subroutine errorExit(message)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            character(len=*),intent(in),optional             ::      message
            if (rank==0) then
                if (present(message)) print *,"MLEFilter "//trim(message)
            end if
#ifdef MPI            
            call MPI_FINALIZE(ierror)
#endif          
            stop
        end subroutine errorExit
        
    




            
    end program MLEFilter 