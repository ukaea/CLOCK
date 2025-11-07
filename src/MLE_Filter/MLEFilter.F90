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

    

        subroutine maxLikelihoodFilter( f,g, t, kernelSize )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      compute a maximum likelihood filtering with characteristic lengthscale given by t
    !*      The probability of finding scene f given that the true solution is g
    !*      is  P = p1 p2
    !*      where   p1 = Exp[ - (f-g)^2/(2 s^2 ) ]
    !*              p2 = Exp[ - ( t^2 del^2(g) )^2/(2 s^2) ]
    !*      so the log likelihood is proportional to
    !*          lambda = - (f-g)^2 - ( t^2 del^2(g) )^2     
    !*      which we can maximise wrt g with
    !*          d lambda / dg = 0
    !*      which is a set of linear equations in f,g,t
    !*          d lambda / dg_i = (-1/s^2)( g_i - f_i + t^4 B_ij g_j )
    !*          d2 lambda / dg_i dg_j = (-1/s^2) ( delta_ij + t^4 B_ij )
    !*      where B is the biharmonic kernel 
    !*          B = del^2 del^2 
    !*      This can be solved by assuming initial g = <f>
    !*      then solving
    !*          ( <f> - f_i )  +  ( delta_ij + t^4 B_ij ) dg_j = 0                
    !*              dldg                     AA
    !*      optionally use a 5x5 biharmonic kernel (poor noise control) kernelSize=5
    !*      or a 7x7 biharmonic kernel (good noise control but slower)  kernelSize=7      
    !*  

    
            real(kind=real64),dimension(:,:),intent(in)     ::      f               !   observed image
            real(kind=real64),dimension(:,:),intent(out)    ::      g               !   ground truth
            real(kind=real64),intent(in)                    ::      t
            integer,intent(in),optional                     ::      kernelSize

        !--     properties of the kernel
            integer                                         ::      mm                      !   kernel half-width
            integer                                         ::      bandwidth               !   number of kernel entries
            real(kind=real64),dimension(:,:),allocatable    ::      kernel                  !   2d biharmonic kernel
            logical,dimension(:,:),allocatable              ::      useKernel0,useKernel    !   non-zero entries

        !---    key properties of the max likelihood estimator
            integer                                         ::      Nx,Ny               !   image size
            integer                                         ::      my_size,my_centre   !   size of sparse matrix A stored, central pixel
            real(kind=real64),dimension(:,:),allocatable    ::      my_AA               !   (1:bandwidth,0:my_size)    -   filter kernel 
            integer,dimension(:,:),allocatable              ::      my_indx             !   (0:bandwidth,0:my_size)    -   index of non-zero rows/cols                
            real(kind=real64),dimension(:),allocatable      ::      dldg,dg             !   ( <f> - f_i ) and solution vector

            real(kind=real64)       ::       eps , fbar  
            integer                 ::      ii,jj,ix,iy,jx,jy,kx,ky , mi,my_ii
            real(kind=real64),dimension(:,:),allocatable    ::      kappa1          !   2d laplacian kernel

 
            
        !---    find the correct biharmonic kernel. 
            mm = 4                              !   kernel half-width
            if (present(kernelSize)) then
                if (kernelSize == 5) mm = 2
                if (kernelSize == 9) mm = 4
            end if
            allocate(kernel(-mm:mm,-mm:mm))
            bandwidth = (2*mm+1)**2             !   number of elements in biharmonic kernel

 
            if (mm==2) then
                allocate(kappa1(-1:1,-1:1))
                !kappa1(-1:1,-1:1) = reshape( (/ 0,1,0, 1,-4,1, 0,1,0 /) , (/3,3/) )
                kappa1 = (t**2) * getKernel3x3( noisy = .false. ) 
                kernel = 0
                do ky = -1,1
                    do kx = -1,1
                        do jy = -1,1
                            do jx = -1,1
                                kernel( jx+kx,jy+ky ) = kernel( jx+kx,jy+ky ) + kappa1(kx,ky)*kappa1(jx,jy)
                            end do
                        end do
                    end do
                end do
            else if (mm==4) then
                allocate(kappa1(-2:2,-2:2))
                kappa1 = (t**2) * getKernel5x5( noisy = .false. )
                kernel = 0
                do ky = -mm/2,mm/2
                    do kx = -mm/2,mm/2
                        do jy = -mm/2,mm/2
                            do jx = -mm/2,mm/2
                                kernel( jx+kx,jy+ky ) = kernel( jx+kx,jy+ky ) + kappa1(kx,ky)*kappa1(jx,jy)
                            end do
                        end do
                    end do
                end do
            end if
        


        !---    find which kernel elements need to be used in first pass of constructing A
            allocate(useKernel0(-mm:mm,-mm:mm))
            allocate(useKernel(-mm:mm,-mm:mm))
            do ky = -mm,mm
                do kx = -mm,mm
                    useKernel0(kx,ky) = abs(kernel(kx,ky))>0
                end do
            end do
  

            
        !---    establish size of problem and allocate memory
            Nx = size(f,dim=1)
            Ny = size(f,dim=2)
            allocate(dg(0:Nx*Ny-1))
            my_size = (Nx*Ny)/nProcs  
            allocate(dldg(0:Nx*Ny-1))
            allocate(my_AA(bandwidth,0:my_size))
            allocate(my_indx(0:bandwidth,0:my_size))
            
            
 
                 

        !---    compute RHS dldg = ( <f> - f_i ) 
            fbar = sum(f)/(Nx*Ny)      
            dldg = fbar - pack(f,.true.)  



        !---    compute matrix
            my_AA = 0
            my_indx = 0
                 
        !---    compute the MLE filter kernel as a sparse matrix
        !       and (minus) first derivative. Assume that "ground truth" g=f for this calc          
            do iy = 1,Ny
                do ix = 1,Nx                
                    ii = ix-1 + (iy-1)*Nx     

                    if (mod(ii,nprocs)/=rank) cycle

                !---    find which parts of the Biharmonic kernel are useable here
                    mi = 0
                    useKernel = useKernel0
                    do ky = -mm,mm
                        jy = ky + iy
                        if ( (jy<1).or.(jy>Ny) ) useKernel(:,ky) = .false.                                
                        do kx = -mm,mm
                            jx = kx + ix
                            if ( (jx<1).or.(jx>Nx) ) useKernel(kx,:) = .false.        
                        end do
                    end do


                !---    compute the matrix A  
                    my_ii = ii/nProcs
                    my_centre = 0                    
                    do ky = -mm,mm
                        do kx = -mm,mm
                            if (useKernel(kx,ky)) then
                                jy = ky + iy                                            
                                jx = kx + ix                            
                                jj = jx-1 + (jy-1)*Nx  
                                mi = mi + 1
                                my_AA(mi,my_ii) = kernel( kx,ky )                                                               
                                my_indx(mi,my_ii) = jj
                                if (ii==jj) my_centre = mi
                            end if
                        end do
                    end do

             !---    correct central element with kronecker delta term, and ensure sum_k A_ki = 0   at edges of image.
                    my_indx(0,my_ii) = mi 
                    eps = sum( my_AA(1:mi,my_ii) )
                    my_AA(my_centre,my_ii) = 1.0d0  + my_AA(my_centre,my_ii) - eps
                     
                end do
            end do


            
            !---    extract solution dldg + A dg = 0
                eps = 1e-10 
                dg = 0
                call mpi_conjgrad( my_AA,my_indx, dg,dldg , eps )     

                    
            !---    pack input observed image into vector                
                g(1:Nx,1:Ny) = fbar - reshape(dg,(/Nx,Ny/))

     
 
            return
        end subroutine maxLikelihoodFilter     


        subroutine mpi_conjgrad( my_A,my_indx, x,b , eps )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      On input A_ki is the kth non-zero column of row i of the sparse symmetric matrix A
    !*      indx_0i is the number of non-zero columns of row i
    !*      indx_ki is the column index of the kth non-zero column
    !*      where t2he data is stored st the ith row of A is row rank + i*nProc 
    !*      x is the initial guess, b is the target
    !*      eps is the tolerance for the solution - want |Ax-b|^2 < N eps where N is the number of rows
    !*      on output
    !*      x is the solution, eps is the error |Ax-b|^2/N
    
            real(kind=real64),dimension(:,0:),intent(in)        ::      my_A
            integer,dimension(0:,:),intent(in)                  ::      my_indx
            real(kind=real64),dimension(0:),intent(in)          ::      b
            real(kind=real64),dimension(0:),intent(inout)       ::      x
            real(kind=real64),intent(inout)                     ::      eps
            
            real(kind=real64),dimension(:),allocatable      ::      my_rr,pp,my_qq,tmp
            real(kind=real64)                               ::      r2,aa,bb
                   
            integer                     ::      ii,kk,NN,maxSteps
                      

            !print *,"rank,nprocs ",rank,nprocs

            NN = size(x)
            allocate(pp(0:NN-1))
            allocate(my_rr(0:NN-1))
            allocate(my_qq(0:NN-1))   
            my_rr = 0           !   remains zero where i%nProcs/=rank
            my_qq = 0             
#ifdef MPI                     
            allocate(tmp(0:NN-1))           
            tmp = 0 
#else
            allocate(tmp(0:1)) !   just to supress unused var warning.           
            tmp = 0 
#endif
            
            call mpi_sparseMatVec( my_A,my_indx,x, my_rr )       

                 
#ifdef MPI
            do ii = rank,NN-1,nprocs
                my_rr(ii) = b(ii) - my_rr(ii)
            end do      
            call MPI_ALLREDUCE( my_rr,pp,NN,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror )                  
#else
            my_rr = b - my_rr                                            
            pp = my_rr
            ii=0
#endif                       
            maxSteps =  ceiling(sqrt(1.0d0*NN))
            
            do kk = 1,maxSteps
            
            !---    r2 = r.r
                r2 = mpi_vecsq(my_rr)                       

                if (r2 <= NN*eps) exit
                !if (rank==0) print *,"mpi_conjgrad r2 ",kk,r2
                       
            !---    q = A p
                call mpi_sparseMatVec( my_A,my_indx,pp, my_qq )                   
                 
            !---    a = p.q
                aa = mpi_dotprod(pp,my_qq)                 
                
            !---    a = r2/a
                if (abs(aa)<1.0d-12) exit
                    !   weird result: p.Ap = 0 can only be true if p is a zero eigenmode of A, or if p=0
                    !   if this is true, the best solution is probably to try again with a different guess vector
                aa = r2/aa
                                
            !---    x = x + a p
                x = x + aa*pp                            
               
            !---    r = r - a q
                my_rr = my_rr - aa*my_qq                                        
                
            !---    aa = r.r/r2
                bb = mpi_vecsq(my_rr)                      !   
                bb = bb/r2          !   note r2>0
                            
             !---   p = r + ap
#ifdef MPI
                do ii = rank,NN-1,nprocs
                    tmp(ii) = my_rr(ii) + bb*pp(ii)
                end do     
                call MPI_ALLREDUCE( tmp,pp,NN,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror )              
#else
                pp = my_rr + bb*pp 
#endif            
                 
            end do
            
        !---    at this point we should have the solution
            call mpi_sparseMatVec( my_A,my_indx,x, my_rr )
            my_rr = b - my_rr                                            
            r2 = mpi_vecsq(my_rr)          
            eps = r2/NN

            
            deallocate(my_rr)
            deallocate(pp)
            deallocate(my_qq)           
          
            
            return
        end subroutine mpi_conjgrad           
            
            
        subroutine mpi_sparseMatVec( my_A,my_indx,x, my_y )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute sparse matrix multiply Ax = y
    !*      On input my_A_ki is the kth non-zero column of the ith non-zero row of the sparse symmetric matrix A
    !*      indx_0i is the number of non-zero columns of row i
    !*      indx_ki is the column index of the kth non-zero column
    !*      note: mpi version only computes for known A rows
    
            real(kind=real64),dimension(:,0:),intent(in)         ::      my_A
            integer,dimension(0:,0:),intent(in)                  ::      my_indx
            real(kind=real64),dimension(0:),intent(in)           ::      x
            real(kind=real64),dimension(0:),intent(out)          ::      my_y
            real(kind=real64)           ::      yi
            integer                     ::      ii,jj,kk,NN  , my_ii    
            NN = size(x)
            my_y = 0.0d0
            do ii = rank,NN-1,nprocs 
                my_ii = ii/nprocs
                yi = 0                                  !   compute y_i = A_ij x_j
                do kk = 1,my_indx(0,my_ii)              !   for each non-zero column at row i
                    jj = my_indx(kk,my_ii)              !   kth non-zero column is indexed j
                    yi = yi + my_A(kk,my_ii)*x(jj)      !   sum A_ij x_j
                end do
                my_y(ii) = yi
            end do
            return
        end subroutine mpi_sparseMatVec

        function mpi_dotprod( a,b ) result(dp)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the parallel dot product a.b assuming different ranks have different parts of a,b set
            real(kind=real64),dimension(0:),intent(in)          ::      a,b
            real(kind=real64)                                   ::      dp
            integer                 ::      ii,ierror
            real(kind=real64)       ::      dp_tmp
            dp_tmp = 0.0d0
            do ii = rank,size(a)-1,nprocs
                dp_tmp = dp_tmp + a(ii)*b(ii)
            end do
#ifdef MPI
            call MPI_ALLREDUCE( dp_tmp,dp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror )
#else
            dp = dp_tmp
            ierror=0
#endif
            return
        end function mpi_dotprod

        function mpi_vecsq( a ) result(dp)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the parallel dot product a.a assuming different ranks have different parts of a  set
            real(kind=real64),dimension(0:),intent(in)          ::      a 
            real(kind=real64)                                   ::      dp
            integer                                 ::      ii,ierror
            real(kind=real64)                       ::      dp_tmp
            dp_tmp = 0.0d0
            do ii = rank,size(a)-1,nprocs
                dp_tmp = dp_tmp + a(ii)*a(ii)
            end do
#ifdef MPI
            call MPI_ALLREDUCE( dp_tmp,dp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror )
#else
            dp = dp_tmp
            ierror=0
#endif
            return
        end function mpi_vecsq

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

        pure real(kind=real64) function numericallaplacianfunc(x,y)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),intent(in)    ::      x,y
            real(kind=real64)               ::      d2x,d2y
            real(kind=real64),parameter     ::      DELTA = 1.0d-4

            d2x = ( -func(x-2*DELTA,y) + 16*func(x-DELTA,y) -30*func(x,y) + 16*func(x+DELTA,y) - func(x+2*DELTA,y) )/12
            d2y = ( -func(x,y-2*DELTA) + 16*func(x,y-DELTA) -30*func(x,y) + 16*func(x,y+DELTA) - func(x,y+2*DELTA) )/12
            
             
            numericallaplacianfunc = ( d2x + d2y  )/(DELTA*DELTA)
            return
        end function numericallaplacianfunc        



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





