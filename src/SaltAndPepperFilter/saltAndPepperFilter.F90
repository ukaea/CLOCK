program saltAndPepperFilter
!---^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    saltAndPepperFilter from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!   A simple standalone code to remove salt and pepper noise from an image
!   The program operates in  steps
!       1)  read in an image from a png file
!       2)  compute a smoothed version
!       3)  find the standard deviation of the difference between smoothed and input pixels
!       4)  where the input image is > 3 std devs away from the smoothed, use the smoothed image pixel instead
!       5)  repeat until converged
!
!   Part of the CLOCK toolkit for automatic image characterization
!
!   Daniel Mason
!   (c) UKAEA Jan 2025
!
!   Version history
!       0.0.1       Jan 2025        First working version
!  
        use iso_fortran_env
        use Lib_Png
        use Lib_Filenames
        use Lib_CommandLineArguments
#ifdef MPI         
        use mpi_f08 
#endif        
        implicit none

    !---    magic numbers fixed in the code
        real(kind=real64),parameter         ::      MAXSTDEV = 5.0d0    !   change pixels if deviating by more than this number of stdevs
        character(len=*),parameter      ::      VERSION = "0.0.1"

    !---    information about the run from command line params
        type(CommandLineArguments)          ::      cla       
        character(len=256)                  ::      filename = ""          !   input filename
        character(len=256)                  ::      outfile = ""


    !---    information about the image
        integer                             ::      Nx,Ny               !   image size in pixels
        real(kind=real64),dimension(:,:),allocatable        ::      img_in              !   (1:Nx,1:Ny)     input image 
        real(kind=real64),dimension(:,:),allocatable        ::      img_smoothed        !   intermediate smoothed image
        real(kind=real64),dimension(:,:),allocatable        ::      img_out             !   output image




    !---    dummy variables
        logical                 ::      ok
        integer                 ::      loop
        integer                 ::      nn,nlast
        real(kind=real64)       ::      mu,sig
        real(kind=real64)       ::      iN

        integer                 ::      nProcs,rank,ierror
                
                 
#ifdef MPI        
        call MPI_INIT(ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)   
#else
        nProcs=1
        rank=0
        ierror=0 !to supress warning if compiled without MPI
#endif      
        
        
!---    read command line arguments
        cla = CommandLineArguments_ctor(30)  
             
        call setProgramDescription( cla, "saltAndPepperFilter" )
        call setProgramVersion( cla, VERSION ) 
       

        !---    filename options        
        call get( cla,"f",filename ,LIB_CLA_REQUIRED,"           image filename" )
        outfile = trim( removeSuffix(filename) )//".sp.png"
        call get( cla,"o",outfile ,LIB_CLA_OPTIONAL,"      output image filename" )                      


        call report(cla)
        if (hasHelpArgument(cla)) stop
        if (.not. allRequiredArgumentsSet(cla)) stop
        call delete(cla)

        if (rank==0) print *,"running on ",nProcs," procs"

    !---    read in the input image
        inquire(file=trim(filename),exist=ok)
        if (.not. ok) call errorExit( "error - file not found """//trim(filename)//"""" )
        if (rank==0) then
            call readPng( filename,img_in )
            Nx = size(img_in,dim=1)
            Ny = size(img_in,dim=2)
        end if
#ifdef MPI
        call MPI_BCAST(Nx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
        call MPI_BCAST(Ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
#endif     
        if (rank/=0) allocate(img_in(Nx,Ny))
#ifdef MPI
        call MPI_BCAST(img_in,Nx*Ny,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
#endif     
        allocate(img_smoothed(Nx,Ny))
        allocate(img_out(Nx,Ny))
        if (rank==0) print *,"saltAndPepperFilter info - read in image with ",Nx,",",Ny," px"
        if (Nx*Ny<=0) call errorExit( "error - incorrect read?? pixel extent error")
        iN = 1.0d0/(Nx*Ny)

        img_out = img_in
        nlast = 0
        do loop = 1,2
        !---    create a smoothed version of the image
            call generateSmoothedImage( img_out,img_smoothed )
            
            
        !---    compute the mean and stdev of the difference between smoothed and input images
            call getMeanAndStdevOfDifference( img_in,img_smoothed, mu,sig )
            if (rank==0) print *,"loop ",loop," mean and stdev of difference ",mu,",",sig

        !---    flatten the most outrageous pixels        
            call flattenOutrageousPixels( img_in,img_smoothed, mu,sig ,img_out , nn )
            if (rank==0) print *,"loop ",loop," flattened ",nn," px (",100.0d0*nn*iN," %)"
            if (nn==nlast) exit
            nlast = nn
        end do

    !---    output the result
        if (rank==0) then
            print *,"write to """//trim(outfile)//""""
            !call writePng( outfile,img_smoothed )
            call writePng( outfile,img_out )
        end if

 
    !---    bye bye
        if (rank==0) print *,""
        call errorExit("done")  

 
    contains
!---^^^^^^^^

    
        subroutine errorExit(message)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            character(len=*),intent(in),optional             ::      message
            if (rank==0) then
                if (present(message)) print *,trim(message)
            end if
#ifdef MPI            
            call MPI_FINALIZE(ierror)
#endif          
            stop
        end subroutine errorExit
        


    !---    compute the mean and stdev of the difference between smoothed and input images
        subroutine getMeanAndStdevOfDifference( img_in,img_smoothed, mu,sig )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(:,:),intent(in)     ::      img_in              !   (1:Nx,1:Ny)     input image 
            real(kind=real64),dimension(:,:),intent(in)     ::      img_smoothed        !   intermediate smoothed image
            real(kind=real64),intent(out)                   ::      mu,sig              !   mean, stdev of difference between input and smoothed
            integer             ::      Nx,Ny
            integer             ::      ii,jj
            real(kind=real64)   ::      delta

            Nx = size(img_in,dim=1)
            Ny = size(img_in,dim=2)            

        !---    compute sum of difference, sum of difference squared
            mu = 0.0            !   will accumulate sum of difference
            sig = 0.0           !   will accumulate sum of difference squared
            do jj = 1,Ny
                do ii = 1,Nx
                    delta = img_in(ii,jj) - img_smoothed(ii,jj)
                    mu = mu + delta
                    sig = sig + delta*delta
                end do
            end do

        !---    convert sums into means by dividing through by pixel count. Note we have checked above Nx*Ny>0
            delta = 1.0d0/(Nx*Ny)
            mu = mu * delta
            sig = sig * delta

        !---    convert sum of squares into a stdev. 
        !       Ignore bias factor sqrt(N/N-1) because we expect N >> 1
            sig = sqrt( max(0.0d0, sig - mu*mu ) )          !   max function avoids sqrt(-0)


            return
        end subroutine getMeanAndStdevOfDifference

!---    flatten the most outrageous pixels        
        subroutine flattenOutrageousPixels( img_in,img_smoothed, mu,sig, img_out,n ) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(:,:),intent(in)     ::      img_in              !   (1:Nx,1:Ny)     input image 
            real(kind=real64),dimension(:,:),intent(in)     ::      img_smoothed        !   intermediate smoothed image
            real(kind=real64),intent(in)                    ::      mu,sig              !   mean, stdev of difference between input and smoothed
            real(kind=real64),dimension(:,:),intent(out)    ::      img_out             !   output smoothed image
            integer,intent(out)                             ::      n                   !   number of pixels flattened
            integer             ::      Nx,Ny
            integer             ::      ii,jj
            integer             ::      i1,i2,j1,j2
            real(kind=real64)   ::      delta,crit
            real(kind=real64),dimension(-1:1,-1:1)          ::      local
            real(kind=real64),dimension(:,:),allocatable    ::      img_tmp

            Nx = size(img_in,dim=1)
            Ny = size(img_in,dim=2)            
            allocate(img_tmp(Nx,Ny))

        !---    first pass- find the absolute difference between the input and smoothed image. Store in img_out
            do jj = 1,Ny
                do ii = 1,Nx
                    img_tmp(ii,jj) = abs( img_in(ii,jj) - img_smoothed(ii,jj) - mu )
                end do
            end do


        !---    second pass- identify the pixels which obey two criteria
        !       1) must be over critical difference
        !       2) must be local maximum 
            crit = max( 1.0d0/256,sig*MAXSTDEV )
            img_out = img_in
            n = 0
            do jj = 1,Ny
                j1 = - min(jj-1,1)
                j2 = min(Ny-jj,1)                
                do ii = 1,Nx
                    delta = img_tmp(ii,jj)               
                    if (abs(delta)>crit) then
                        i1 = - min(ii-1,1)
                        i2 = min(Nx-ii,1) 
                        local(i1:i2,j1:j2) = img_tmp(ii+i1:ii+i2,jj+j1:jj+j2)
                        local(0,0) = 0
                        if (delta > maxval( local(i1:i2,j1:j2) ) ) then
                          !  print *,"changing pixel at ",ii,jj," img in,smoothed,tmp ",img_in(ii,jj),img_smoothed(ii,jj),img_tmp(ii,jj)
                            img_out(ii,jj) = img_smoothed(ii,jj)
                            n = n + 1
                        end if
                    end if
                end do
            end do

            return
        end subroutine flattenOutrageousPixels



    !---    create a smoothed version of the image
        subroutine generateSmoothedImage( img_in,img_smoothed )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(:,:),intent(in)     ::      img_in              !   (1:Nx,1:Ny)     input image 
            real(kind=real64),dimension(:,:),intent(out)    ::      img_smoothed        !   intermediate smoothed image

            integer             ::      Nx,Ny
            integer             ::      ii,jj
            integer             ::      i1,i2,j1,j2
            real(kind=real64),dimension(-2:2,-2:2)  ::      ff
#ifdef MPI            
            real(kind=real64),dimension(:,:),allocatable    ::      img_tmp
#endif

            Nx = size(img_in,dim=1)
            Ny = size(img_in,dim=2)          

            img_smoothed = 0

            do jj = 1,Ny
                if (mod(jj-1,nprocs)/=rank) cycle
                j1 = - min(jj-1,2)
                j2 = min(Ny-jj,2)                
                do ii = 1,Nx
                    i1 = - min(ii-1,2)
                    i2 = min(Nx-ii,2)                
                    
                    ff = 0  !   not needed, for debug
                    ff(i1:i2,j1:j2) = img_in(ii+i1:ii+i2,jj+j1:jj+j2)




                    img_smoothed(ii,jj) = estimateCentralPixel( ff,i1,i2,j1,j2 )

                    ! !   !if ( kernelIndex( i1,i2,j1,j2 )==1) then
                    ! if ((ii==1).and.( abs(jj-276)<3 )) then
                    !       print *,ii,jj,"i1,i2,j1,j2", i1,i2,j1,j2,"indx", kernelIndex( i1,i2,j1,j2 )
                    !       write (*,fmt='(5f12.5)') ff(:,-2)
                    !       write (*,fmt='(5f12.5)') ff(:,-1)
                    !       write (*,fmt='(5f12.5)') ff(:, 0)
                    !       write (*,fmt='(5f12.5)') ff(:, 1)
                    !       write (*,fmt='(5f12.5)') ff(:, 2)
                    !       print *,img_smoothed(ii,jj)


                    !       call rotp90(ff)
                    !       write (*,fmt='(5f12.5)') ff(:,-2)
                    !       write (*,fmt='(5f12.5)') ff(:,-1)
                    !       write (*,fmt='(5f12.5)') ff(:, 0)
                    !       write (*,fmt='(5f12.5)') ff(:, 1)
                    !       write (*,fmt='(5f12.5)') ff(:, 2)
                    !       !stop
                    !   end if

                end do
            end do
        
#ifdef MPI
            allocate(img_tmp(Nx,Ny))
            img_tmp = img_smoothed
            call MPI_ALLREDUCE(img_tmp,img_smoothed,Nx*Ny,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
#endif  
        
            return
        end subroutine generateSmoothedImage

    !---    given a rectangular block of pixels from i1:i2 and j1:j2
    !       with i1 = -2,-1,0  and i2 = 0,1,2
    !       compute an estimate for the central pixel (0,0)
    !       based on minimising the merit function
    !           S = (p(x,y) - f(x,y))^2
    !       where p is a polynomial function
    !           p(x,y) = p0 + p_x x + p_y y + ( p_xx x^2 + 2 p_xy x y + p_yy y^2 )/2 
    !       and f(x,y) is the values observed in a 5x5 block (excluding corners)
    !       The kernels have been computed with mathematica

        real(kind=real64) function estimateCentralPixel( f,i1,i2,j1,j2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)                                          ::      i1,i2,j1,j2
            real(kind=real64),dimension(-2:2,-2:2),intent(in)           ::      f


        !--- kernel for best fit central px given all neighbours
            real(kind=real64),dimension(-2:2,-2:2),parameter    ::      kernel00 = reshape( (/      &           !     X X X
                    0   ,   -14 ,   3   ,   -14 ,   0   ,                                           &           !   X X X X X   
                    -14 ,   37  ,   54  ,   37  ,   -14 ,                                           &           !   X X 0 X X
                    3   ,   54  ,   0   ,   54  ,   3   ,                                           &           !   X X X X X
                    -14 ,   37  ,   54  ,   37  ,   -14 ,                                           &           !     X X X
                    0   ,   -14 ,   3   ,   -14 ,   0               /)  ,  (/5,5/) ) / 264.0d0

        !--- kernel for best fit central px bounded on one side
            real(kind=real64),dimension(-2:2,-1:2),parameter    ::      kernel01 = reshape( (/      &           !   X X X X X   
                    -127,   128 ,   213 ,   128 ,   -127,                                           &           !   X X 0 X X
                    35  ,   290 ,   0   ,   290 ,   35  ,                                           &           !   X X X X X
                    -41 ,   214 ,   299 ,   214 ,   -41 ,                                           &           !     X X X
                    0   ,   -100,   -15 ,  -100 ,   0               /)  ,  (/5,4/) ) / 1295.0d0

        !--- kernel for best fit central px severely bounded on one side
            real(kind=real64),dimension(-2:2,0:2),parameter     ::      kernel02 = reshape( (/      &           !   X X 0 X X
                    13  ,   58  ,   0   ,   58  ,   13  ,                                           &           !   X X X X X
                    -30 ,   15  ,   30  ,   15 ,   -30  ,                                           &           !     X X X
                    0   ,   -5  ,   10  ,  -5  ,   0               /)  ,  (/5,3/) ) / 142.0d0

        !--- kernel for best fit central px bounded on one corner
            real(kind=real64),dimension(-1:2,-1:2),parameter    ::      kernel11 = reshape( (/      &           !   X X X X  
                    -11 ,   97  ,   80  ,   -62 ,                                                   &           !   X 0 X X
                    97  ,   0   ,   162 ,   7   ,                                                   &           !   X X X X
                    80  ,   162 ,   119 ,   -49 ,                                                   &           !   X X X
                    -62 ,   7   ,   -49 ,   0       /)  ,  (/4,4/) ) / 578.0d0

        !--- kernel for best fit central px severely bounded on one side and corner
            real(kind=real64),dimension(-1:2,0:2),parameter     ::      kernel12 = reshape( (/      &           !   X 0 X X
                    39  ,   0   ,   47  ,   6   ,                                                   &           !   X X X X
                    -3  ,   19  ,   11  ,   -27 ,                                                   &           !   X X X
                    -15 ,   10  ,   5   ,   0        /)  ,  (/4,3/) ) / 92.0d0


        !--- kernel for best fit central px severely bounded on one side and corner
            real(kind=real64),dimension(0:2,-1:2),parameter     ::      kernel21 = reshape( (/      &           !   X X X  
                    39  ,   -3  ,   -15 ,                                                           &           !   0 X X
                    0   ,   19  ,   10  ,                                                           &           !   X X X
                    47  ,   11  ,   5   ,                                                           &           !   X X 
                    6   ,   -27 ,   0                 /)  ,  (/3,4/) ) / 92.0d0


        !--- kernel for best fit central px severely bounded on both sides
            real(kind=real64),dimension(0:2,0:2),parameter      ::      kernel22 = reshape( (/      &           !   0 X X
                    0   ,   3   ,   -1  ,                                                           &           !   X X X
                    3   ,   -4  ,   1   ,                                                           &           !   X X
                    -1  ,   1   ,   0                 /)  ,  (/3,3/) ) / 2.0d0

                    
            real(kind=real64),dimension(-2:2,-2:2)              ::      ff
            integer             ::      indx,kk

        !---    make a unique indexing for each possible input case 
            indx =  kernelIndex( i1,i2,j1,j2 )

!      0                1               2           3               6
!         X X X           X X X            X X           X X X           X X  
!       X X X X X         X X X X          X X X       X X X X         X X X     
!       X X 0 X X         X 0 X X          0 X X       X X 0 X         X X 0   
!       X X X X X         X X X X          X X X       X X X X         X X X  
!         X X X           X X X            X X           X X X           X X  
!       
!      9                10              11          12              15
!          
!       X X X X X         X X X X          X X X       X X X X         X X X     
!       X X 0 X X         X 0 X X          0 X X       X X 0 X         X X 0   
!       X X X X X         X X X X          X X X       X X X X         X X X  
!         X X X           X X X            X X           X X X           X X  
!       
!      18               19              20          21              24
!       
!       
!       X X 0 X X         X 0 X X          0 X X       X X 0 X         X X 0  
!       X X X X X         X X X X          X X X       X X X X         X X X  
!         X X X           X X X            X X           X X X           X X  
!       
!      27               28              29          30              33
!         X X X           X X X            X X           X X X           X X  
!       X X X X X         X X X X          X X X       X X X X         X X X     
!       X X 0 X X         X 0 X X          0 X X       X X 0 X         X X 0  
!       X X X X X         X X X X          X X X       X X X X         X X X  
!          
!       
!      54               55              56          57              60
!         X X X           X X X            X X           X X X           X X  
!       X X X X X         X X X X          X X X       X X X X         X X X     
!       X X 0 X X         X 0 X X          0 X X       X X 0 X         X X 0  
!      

        !---    rotate the input to standard form 
            ff = 0
            select case(indx)
                case (0)
                    ff(-2:2,-2:2) = f(-2:2,-2:2)
                    kk = 00
                case (1)
                    ff(-1:2,-2:2) = f(-1:2,-2:2)
                    call rotp90(ff)
                    kk = 01
                case (2)
                    ff( 0:2,-2:2) = f( 0:2,-2:2)
                    call rotp90(ff)
                    kk = 02
                case (3)
                    ff(-2:1,-2:2) = f(-2:1,-2:2)
                    call rotm90(ff)
                    kk = 01                    
                case (6)
                    ff(-2:0,-2:2) = f(-2:0,-2:2)
                    call rotm90(ff)
                    kk = 02                

                case (9)
                    ff(-2:2,-1:2) = f(-2:2,-1:2)
                    kk = 01
                case (10)
                    ff(-1:2,-1:2) = f(-1:2,-1:2)
                    kk = 11
                case (11)
                    ff( 0:2,-1:2) = f( 0:2,-1:2)
                    kk = 21
                case (12)
                    ff(-2:1,-1:2) = f(-2:1,-1:2)
                !     write (*,fmt='(5f12.5)') ff(:,-2)
                ! write (*,fmt='(5f12.5)') ff(:,-1)
                ! write (*,fmt='(5f12.5)') ff(:, 0)
                ! write (*,fmt='(5f12.5)') ff(:, 1)
                ! write (*,fmt='(5f12.5)') ff(:, 2)
                    call xflip(ff)
                    kk = 11
                case (15)
                    ff(-2:0,-1:2) = f(-2:0,-1:2)
                    call rotm90(ff)
                    kk = 12

                case (18)
                    ff(-2:2, 0:2) = f(-2:2, 0:2)
                    kk = 02
                case (19)
                    ff(-1:2, 0:2) = f(-1:2, 0:2)
                    kk = 12
                case (20)
                    ff( 0:2, 0:2) = f( 0:2, 0:2)
                    kk = 22
                case (21)
                    ff(-2:1, 0:2) = f(-2:1, 0:2)
                    call xflip(ff)
                    kk = 12
                case (24)
                    ff(-2:0, 0:2) = f(-2:0, 0:2)
                    call xyflip(ff)
                    kk = 22

                case (27)
                    ff(-2:2,-2:1) = f(-2:2,-2:1)
                    call yflip(ff)
                    kk = 01
                case (28)
                    ff(-1:2,-2:1) = f(-1:2,-2:1)
                    call rotp90(ff)
                    kk = 11                
                case (29)
                    ff( 0:2,-2:1) = f( 0:2,-2:1)
                    call rotp90(ff)
                    kk = 12
                case (30)
                    ff(-2:1,-2:1) = f(-2:1,-2:1)
                    call xyflip(ff)
                    kk = 11
                case (33)
                    ff(-2:0,-2:1) = f(-2:0,-2:1)
                    call xyflip(ff)
                    kk = 21

                case (54)
                    ff(-2:2,-2:0) = f(-2:2,-2:0)
                    call yflip(ff)
                    kk = 02
                case (55)
                    ff(-1:2,-2:0) = f(-1:2,-2:0)
                    call yflip(ff)
                    kk = 12
                case (56)
                    ff( 0:2,-2:0) = f( 0:2,-2:0)
                    call yflip(ff)
                    kk = 22                    
                case (57)
                    ff(-2:1,-2:0) = f(-2:1,-2:0)
                    call xyflip(ff)
                    kk = 12
                case (60)
                    ff(-2:0,-2:0) = f(-2:0,-2:0)
                    call xyflip(ff)
                    kk = 22
 

                case default
                    print *,i1,i2,j1,j2,indx
                    call errorExit("estimateCentralPixel error - incorrect code")

            end select



        !---    now find the kernel compute
            select case (kk)
                case(00)
                    estimateCentralPixel = sum( kernel00 * ff )
                case(01)
                    estimateCentralPixel = sum( kernel01 * ff(-2:2,-1:2) )
                case(02)
                    estimateCentralPixel = sum( kernel02 * ff(-2:2,0:2) )
                case(11)
                    estimateCentralPixel = sum( kernel11 * ff(-1:2,-1:2) )
                case(12)
                    estimateCentralPixel = sum( kernel12 * ff(-1:2,0:2) )
                case(21)
                    estimateCentralPixel = sum( kernel21 * ff(0:2,-1:2) )
                case(22)
                    estimateCentralPixel = sum( kernel22 * ff(0:2,0:2) )
                    
                case default
                    print *,kk
                    call errorExit("estimateCentralPixel error - incorrect kernel")

            end select 

            
! ! 
!             if ( indx == 1 ) then
!                 print *,"i1,i2,j1,j2", i1,i2,j1,j2,"indx", indx
!                 write (*,fmt='(5f12.5)') ff(:,-2)
!                 write (*,fmt='(5f12.5)') ff(:,-1)
!                 write (*,fmt='(5f12.5)') ff(:, 0)
!                 write (*,fmt='(5f12.5)') ff(:, 1)
!                 write (*,fmt='(5f12.5)') ff(:, 2)
!                 print *,estimateCentralPixel
!               !  stop
!             end if
! 

            return

        end function estimateCentralPixel

    !---    find a unique code for the position of the kernel
        pure integer function kernelIndex( i1,i2,j1,j2 )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            integer,intent(in)          ::      i1,i2,j1,j2
            kernelIndex = (i1+2) + 3*(2-i2) + 9*(j1+2) + 27*(2-j2)
            return
        end function kernelIndex

        pure subroutine xflip(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(-2:2,-2:2),intent(inout)        ::  f
            real(kind=real64)           ::      dd
            dd = f(-1,-2) ; f(-1,-2) = f( 1,-2) ; f( 1,-2) = dd
            dd = f(-1,-1) ; f(-1,-1) = f( 1,-1) ; f( 1,-1) = dd
            dd = f(-1, 0) ; f(-1, 0) = f( 1, 0) ; f( 1, 0) = dd
            dd = f(-1, 1) ; f(-1, 1) = f( 1, 1) ; f( 1, 1) = dd
            dd = f(-1, 2) ; f(-1, 2) = f( 1, 2) ; f( 1, 2) = dd
            dd = f(-2,-1) ; f(-2,-1) = f( 2,-1) ; f( 2,-1) = dd
            dd = f(-2, 0) ; f(-2, 0) = f( 2, 0) ; f( 2, 0) = dd
            dd = f(-2, 1) ; f(-2, 1) = f( 2, 1) ; f( 2, 1) = dd    
            return
        end subroutine xflip           
                    
        pure subroutine yflip(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(-2:2,-2:2),intent(inout)        ::  f
            real(kind=real64)           ::      dd
            dd = f(-1,-2) ; f(-1,-2) = f(-1, 2) ; f(-1, 2) = dd
            dd = f( 0,-2) ; f( 0,-2) = f( 0, 2) ; f( 0, 2) = dd
            dd = f( 1,-2) ; f( 1,-2) = f( 1, 2) ; f( 1, 2) = dd
            dd = f(-2,-1) ; f(-2,-1) = f(-2, 1) ; f(-2, 1) = dd
            dd = f(-1,-1) ; f(-1,-1) = f(-1, 1) ; f(-1, 1) = dd
            dd = f( 0,-1) ; f( 0,-1) = f( 0, 1) ; f( 0, 1) = dd
            dd = f( 1,-1) ; f( 1,-1) = f( 1, 1) ; f( 1, 1) = dd
            dd = f( 2,-1) ; f( 2,-1) = f( 2, 1) ; f( 2, 1) = dd
            return
        end subroutine yflip      

        pure subroutine xyflip(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(-2:2,-2:2),intent(inout)        ::  f
            real(kind=real64)           ::      dd
            dd = f(-1,-2) ; f(-1,-2) = f( 1, 2) ; f( 1, 2) = dd
            dd = f( 0,-2) ; f( 0,-2) = f( 0, 2) ; f( 0, 2) = dd
            dd = f( 1,-2) ; f( 1,-2) = f(-1, 2) ; f(-1, 2) = dd
            dd = f(-2,-1) ; f(-2,-1) = f( 2, 1) ; f( 2, 1) = dd
            dd = f(-1,-1) ; f(-1,-1) = f( 1, 1) ; f( 1, 1) = dd
            dd = f( 0,-1) ; f( 0,-1) = f( 0, 1) ; f( 0, 1) = dd
            dd = f( 1,-1) ; f( 1,-1) = f(-1, 1) ; f(-1, 1) = dd
            dd = f( 2,-1) ; f( 2,-1) = f(-2, 1) ; f(-2, 1) = dd
            dd = f(-2, 0) ; f(-2, 0) = f( 2, 0) ; f( 2, 0) = dd
            dd = f(-1, 0) ; f(-1, 0) = f( 1, 0) ; f( 1, 0) = dd
            return
        end subroutine xyflip      

        pure subroutine rotp90(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(-2:2,-2:2),intent(inout)        ::  f
            real(kind=real64),dimension(-2:2,-2:2)      ::      f_tmp
            f_tmp = f
            f(-1,-2) = f_tmp(-2, 1)
            f( 0,-2) = f_tmp(-2, 0)
            f( 1,-2) = f_tmp(-2,-1)

            f(-2,-1) = f_tmp(-1, 2)
            f(-1,-1) = f_tmp(-1, 1)
            f( 0,-1) = f_tmp(-1, 0)
            f( 1,-1) = f_tmp(-1,-1)
            f( 2,-1) = f_tmp(-1,-2)

            f(-2, 0) = f_tmp( 0, 2)
            f(-1, 0) = f_tmp( 0, 1)
            f( 1, 0) = f_tmp( 0,-1)
            f( 2, 0) = f_tmp( 0,-2)

            f(-2, 1) = f_tmp( 1, 2)
            f(-1, 1) = f_tmp( 1, 1)
            f( 0, 1) = f_tmp( 1, 0)
            f( 1, 1) = f_tmp( 1,-1)
            f( 2, 1) = f_tmp( 1,-2)

            f(-1, 2) = f_tmp( 2, 1)
            f( 0, 2) = f_tmp( 2, 0)
            f( 1, 2) = f_tmp( 2,-1)
            
            return
        end subroutine rotp90      
        
        pure subroutine rotm90(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(-2:2,-2:2),intent(inout)        ::  f
            real(kind=real64),dimension(-2:2,-2:2)      ::      f_tmp
            f_tmp = f
            f(-2, 1) = f_tmp(-1,-2)
            f(-2, 0) = f_tmp( 0,-2)
            f(-2,-1) = f_tmp( 1,-2)
            f(-1, 2) = f_tmp(-2,-1)
            f(-1, 1) = f_tmp(-1,-1)
            f(-1, 0) = f_tmp( 0,-1)
            f(-1,-1) = f_tmp( 1,-1)
            f(-1,-2) = f_tmp( 2,-1)
            f( 0, 2) = f_tmp(-2, 0)
            f( 0, 1) = f_tmp(-1, 0)
            f( 0,-1) = f_tmp( 1, 0)
            f( 0,-2) = f_tmp( 2, 0)
            f( 1, 2) = f_tmp(-2, 1)
            f( 1, 1) = f_tmp(-1, 1)
            f( 1, 0) = f_tmp( 0, 1)
            f( 1,-1) = f_tmp( 1, 1)
            f( 1,-2) = f_tmp( 2, 1)
            f( 2, 1) = f_tmp(-1, 2)
            f( 2, 0) = f_tmp( 0, 2)
            f( 2,-1) = f_tmp( 1, 2)
            
            return
        end subroutine rotm90      
        

    end program saltAndPepperFilter 

        
        
        
        