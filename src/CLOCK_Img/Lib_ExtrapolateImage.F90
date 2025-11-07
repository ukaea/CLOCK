!   A simple module to extrapolate a 2d image out past its boundaries using 1st and 2nd derivative kernels to estimate the function at the boundary
!
!   Daniel Mason
!   (c) UKAEA Feb 2025
!
!   Version history
!   0.0.1       Feb 2025    First working version
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_ExtrapolateImage from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
    module Lib_ExtrapolateImage
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
        use iso_fortran_env
#ifdef MPI
        use mpi_f08
#endif
        implicit none
        private
		
    !---    define the 5x5 kernels used to find derivatives
        real(kind=real64),dimension(-2:2,-2:2),private,parameter        ::      G0 = reshape(  (/       &
                         -13,   2,   7,   2, -13,           &
                           2,  17,  22,  17,   2,           &
                           7,  22,  27,  22,   7,           &
                           2,  17,  22,  17,   2,           &
                         -13,   2,   7,   2, -13            /) , (/5,5/) ) / 175.0d0

        real(kind=real64),dimension(-2:2,-2:2),private,parameter        ::      Gx = reshape(  (/       &
                          -2,  -1,   0,   1,   2,           &
                          -2,  -1,   0,   1,   2,           &
                          -2,  -1,   0,   1,   2,           &
                          -2,  -1,   0,   1,   2,           &
                          -2,  -1,   0,   1,   2            /) , (/5,5/) ) / 150.0d0

        real(kind=real64),dimension(-2:2,-2:2),private,parameter        ::      Gxx = reshape( (/       &
                           2,  -1,  -2,  -1,   2,           &
                           2,  -1,  -2,  -1,   2,           &
                           2,  -1,  -2,  -1,   2,           &
                           2,  -1,  -2,  -1,   2,           &
                           2,  -1,  -2,  -1,   2            /) , (/5,5/) ) / 35.0d0

        real(kind=real64),dimension(-2:2,-2:2),private,parameter        ::      Gy = reshape(  (/       &
                          -2,  -2,  -2,  -2,  -2,           &
                          -1,  -1,  -1,  -1,  -1,           &
                           0,   0,   0,   0,   0,           &
                           1,   1,   1,   1,   1,           &
                           2,   2,   2,   2,   2            /) , (/5,5/) ) / 150.0d0

        real(kind=real64),dimension(-2:2,-2:2),private,parameter        ::      Gyy = reshape( (/       &
                           2,   2,   2,   2,   2,           &
                          -1,  -1,  -1,  -1,  -1,           &
                          -2,  -2,  -2,  -2,  -2,           &
                          -1,  -1,  -1,  -1,  -1,           &
                           2,   2,   2,   2,   2            /) , (/5,5/) ) / 35.0d0
 
        real(kind=real64),dimension(-2:2,-2:2),private,parameter        ::      Gxy = reshape( (/       &
                          -4,  -2,   0,   2,   4,           &
                          -2,  -1,   0,   1,   2,           &
                           0,   0,   0,   0,   0,           &
                           2,   1,   0,  -1,  -2,           &
                           4,   2,   0,  -2,  -4            /) , (/5,5/) ) / 100.0d0

        integer,public,parameter                    ::      FADE_STYLE_NONE = 0             !   use polynomial extrapolation only
        integer,public,parameter                    ::      FADE_STYLE_ZERO = 1             !   fade to zero at long range
        integer,public,parameter                    ::      FADE_STYLE_AVG = 2              !   fade to average at long range
        
        integer,private     ::      nProcs = 1, rank = 0 
        
        public              ::      extrapolate

		
    contains
!---^^^^^^^^

        subroutine extrapolate( img_in, border, img_out, fade )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(0:,0:),intent(in)                   ::      img_in
            integer,intent(in)                                              ::      border
            real(kind=real64),dimension(-border:,-border:),intent(out)      ::      img_out
            integer,intent(in)                                              ::      fade

            integer             ::      Nx,Ny
            integer             ::      ii,jj
            real(kind=real64),dimension(-2:2,-2:2)      ::      ff,gg

            if (border<=0) return       !   nothing to do

            do ii = -2,2
                write(*,fmt='(100f9.3)') G0(ii,:),Gx(ii,:),Gy(ii,:),Gxx(ii,:),Gyy(ii,:),Gxy(ii,:)
            end do



        !---    determine size of problem
            Nx = size(img_in,dim=1)  
            Ny = size(img_in,dim=2)  

            if (min(Nx,Ny)<5) call errorExit("error - can't extrapolate a small image < 5x5 px")
            if ( (size(img_out,dim=1)<Nx+2*border).or.(size(img_out,dim=2)<Ny+2*border) ) call errorExit("error - insufficient boundary for img_out")
            img_out(0:Nx-1,0:Ny-1) = img_in(0:Nx-1,0:Ny-1)

            do ii = 1,border                            !   pixel row/col to modify
                
            !---    north west
                ff = img_in( 0:4,0:4 )
                do jj = -ii,1
                    gg = kernel( jj - 2 , -2-ii )
                    img_out( jj,0-ii ) = img_in( 0,0) + extrap( ff,gg )
                end do
                do jj = 1-ii,1
                    gg = kernel( -2-ii, jj-2 )
                    img_out( 0-ii,jj ) = img_in( 0,0) + extrap( ff,gg )
                end do                

            !---    north
                gg = kernel( 0,-2-ii )
                do jj = 2,Nx-3
                    ff = img_in( jj-2:jj+2,0:4 )
                    img_out( jj,0-ii ) = img_in( jj,0) + extrap( ff,gg )
                end do

            !---    north east
                ff = img_in( Nx-5:Nx-1,0:4 )
                do jj = Nx-2,Nx-1+ii
                    gg = kernel( jj + 3 - Nx, -2-ii )
                    img_out( jj,0-ii ) = img_in( Nx-1,0) + extrap( ff,gg )
                end do
                do jj = 1-ii,1
                    gg = kernel( 2+ii , jj-2 )
                    img_out( Nx-1+ii,jj ) = img_in( Nx-1,0) + extrap( ff,gg )
                end do

            !---    east
                gg = kernel( 2+ii,0 )
                do jj = 2,Ny-3
                    ff = img_in( Nx-5:Nx-1,jj-2:jj+2 )
                    img_out( Nx-1+ii,jj ) = img_in( Nx-1,jj) + extrap( ff,gg )
                end do

            !---    south west
                ff = img_in( 0:4,Ny-5:Ny-1 )
                do jj = -ii,1
                    gg = kernel( jj - 2 , 2+ii )
                    img_out( jj,Ny-1+ii ) = img_in( 0,Ny-1) + extrap( ff,gg )
                end do
                do jj = Ny-2,Ny-2+ii
                    gg = kernel( -2-ii, jj + 3 - Ny )
                    img_out( 0-ii,jj ) = img_in( 0,Ny-1) + extrap( ff,gg )
                end do        

            !---    south
                gg = kernel( 0,2+ii )
                do jj = 2,Nx-3
                    ff = img_in( jj-2:jj+2,Ny-5:Ny-1 )
                    img_out( jj,Ny-1+ii ) = img_in( jj,Ny-1) + extrap( ff,gg )
                end do

            !---    south east
                ff = img_in( Nx-5:Nx-1,Ny-5:Ny-1 )
                do jj = Nx-2,Nx-1+ii
                    gg = kernel( jj + 3 - Nx, 2+ii )
                    img_out( jj,Ny-1+ii ) = img_in( Nx-1,Ny-1) +  extrap( ff,gg )
                end do
                do jj = Ny-2,Ny-2+ii
                    gg = kernel( 2+ii, jj + 3 - Ny )
                    img_out( Nx-1+ii,jj ) = img_in( Nx-1,Ny-1) + extrap( ff,gg )
                end do                

            !---    west
                gg = kernel( -2-ii,0 )
                do jj = 2,Ny-3
                    ff = img_in( 0:4,jj-2:jj+2 )
                    img_out( 0-ii,jj ) = img_in( 0,jj) + extrap( ff,gg )
                end do

            end do


            call addFade( img_out,border,fade )
            return
        end subroutine extrapolate

        subroutine addFade( img,border,fade )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      add the requested fade to black/avg value
            integer,intent(in)                                              ::      border
            real(kind=real64),dimension(-border:,-border:),intent(inout)    ::      img 
            integer,intent(in)                                              ::      fade

            integer             ::      Nx,Ny
            integer             ::      ii,jj
            real(kind=real64)   ::      sigma,dd,level,ff


            if (fade == FADE_STYLE_NONE) return


        !---    determine size of problem
            Nx = size(img,dim=1) - 2*border
            Ny = size(img,dim=2) - 2*border
            sigma = 5.0d0 

        !---    determine fade level
            level = 0.0d0
            if (fade == FADE_STYLE_AVG) level = sum(img(0:Nx-1,0:Ny-1))/(Nx*Ny)

            do ii = 1,border                            !   pixel row/col to modify
                dd = exp( - ii*ii/(2*sigma*sigma) )     !   fraction of original to take

            !---    east
                do jj = -border,Ny-1+border         !   pixel col
                    ff = img( Nx-1+ii,jj )          !   original pixel
                    ff = ff * dd + level*(1-dd)     !   new pixel
                    img( Nx-1+ii,jj ) = ff             !   reset
                end do

            !---    north
                do jj = -border,Nx-1+border         !   pixel col
                    ff = img( jj,0-ii )             !   original pixel
                    ff = ff * dd + level*(1-dd)     !   new pixel
                    img( jj,0-ii ) = ff             !   reset
                end do

            !---    south
                do jj = -border,Nx-1+border         !   pixel col
                    ff = img( jj,Ny-1+ii )          !   original pixel
                    ff = ff * dd + level*(1-dd)     !   new pixel
                    img( jj,Ny-1+ii ) = ff          !   reset
                end do

            !---    west
                do jj = -border,Ny-1+border         !   pixel col
                    ff = img( 0-ii,jj )             !   original pixel
                    ff = ff * dd + level*(1-dd)     !   new pixel
                    img( 0-ii,jj ) = ff             !   reset
                end do
            
            end do

            return
        end subroutine addFade


        pure real(kind=real64) function extrap( f,g )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(-2:2,-2:2),intent(in)       ::      f,g
            extrap = sum( f*g )
            extrap = tanh( extrap ) / 2
            return
        end function extrap

        pure function kernel( dx,dy ) result(g)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      Find a kernel which will give an estimate for the value at dx,dy
            integer,intent(in)                          ::      dx,dy
            real(kind=real64),dimension(-2:2,-2:2)      ::      g
            !g = G0 + Gx*dx + Gy*dy + ( Gxx*dx*dx + Gyy*dy*dy + 2*Gxy*dx*dy )/2
            g =  Gx*dx + Gy*dy + ( Gxx*dx*dx + Gyy*dy*dy + 2*Gxy*dx*dy )/2
            return
        end function kernel



        subroutine errorExit(message)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      exit with a message
            character(len=*),intent(in),optional             ::      message
#ifdef MPI            
            integer             ::      ierror
            call MPI_BARRIER(MPI_COMM_WORLD,ierror)
            if (rank==0) print *,""
            call MPI_BARRIER(MPI_COMM_WORLD,ierror)
            nProcs=1
#else
            print *,""
            nProcs=1
#endif                    
            if ((rank==0) .and. present(message)) print *,"Lib_ExtrapolateImage::"//trim(message)
#ifdef MPI            
            call MPI_FINALIZE(ierror)
#endif        
            stop
        end subroutine errorExit
        
		
!-------

        

        
    end module Lib_ExtrapolateImage