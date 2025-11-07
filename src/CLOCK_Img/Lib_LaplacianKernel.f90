
    module Lib_LaplacianKernel
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_LaplacianKernel from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!---^^^^^^^^^^^^^^^^^^^^^^^^^^ 
!*      simple code to compute the 3x3x3 discrete Laplacian operator kernel in 3d
!*      or the 3x3 kernel in 2d
!*
!*      We need to compute quantities like
!*          L₁ = sum_i  sum_{j ∈ N_i}   kappa_ij f_j                 (1)
!*          L₂ = sum_i  sum_{j ∈ N_i} | kappa_ij f_j |               (2)
!*          L₃ = sum_i  sum_{j ∈ N_i} ( kappa_ij f_j )²              (3)
!*
!*      with i,j being labels for points in (2|3)d space, ie i=(ix,iy) or i=(ix,iy,iz)
!*      
!*      Note that 
!*          d L₁ / d f_i     =   sum_{j ∈ N_i} kappa_ij              
!*          d L₂ / d f_i     =   sum_{j ∈ N_i} kappa_ij sign( L1_j )
!*          d L₃ / d f_i     =   sum_{j ∈ N_i} 2 kappa_ij L1_j 
!*

!*      the call discreteLaplacian_deriv(i,j [,k]) returns the 3x3 kernel elements
!*
!*      
!*      We may also need to find laplacians of concentrations, where sum_a f_aj = f_j
!*          G_ai                = sum_{j ∈ N_i} ( kappa_ij f_aj/f_j )
!*          F_i                 = (f_i/2) sum_a ( G_ia/f_i )^2       
!*          d G_ai / d f_bj     = (kappa_ij/f_j) ( delta_ab - f_aj/f_j )
!*          d F_i / d f_bj      = ( sum_a G_ai d G_ai / d f_bj  - delta_ij F_j ) / f_i
!*          d² G_ai / d f_bj d f_ck     = (delta_jk kappa_ij / rho_j^2) ( 2 f_aj/f_j - delta_ab - delta_ac )
!*          d² F_i / d f_bj d f_ck      = ( sum_a (1/f_i) d G_ai / d f_bj  d G_ai / d f_ck + sum_a G_ai  d² G_ai / d f_bj d f_ck - (1/f_i) ( delta_ij  d F_i / d f_bj +  delta_ik d F_i / d f_ck )/f_i
!*      

!*      Daniel Mason
!*      (c) UKAEA June 2024
!*
!*      Version history
!*          0.0.1           Jun 2024        First working version
!*          0.0.2           Nov 2024        updated kernels for 5x5,5x5x5. removal of buggy "square" lapalacian kernels
!*          0.1.0           Feb 2025        Two kernels are provided for 5x5, 5x5x5. Default assumes noisy background, can switch to high accuracy with noisy=.false.


        use iso_fortran_env
        implicit none
        private

        include "Lib_LaplacianKernel_data.h"

        integer,public,parameter        ::      LIB_LK_INDX_INVSUM = 0

 

        public          ::      discreteLaplacian               !   note: returns ∑_j k_ij f_j = ∑_jx=-1^+1 ∑_jy=-1^+1  kernel(jx,jy) f(ix+jx,iy+jy)
        public          ::      discreteLaplacian_deriv         !   note: returns k_ij
        public          ::      discreteBiharmonic         !   note: returns the vector ∑_j ∑_l k_ij k_jl f_l = ∑_jx=-2^+2 ∑_jy=-2^+2  sq_kernel(jx,jy) f(ix+jx,iy+jy)
        public          ::      discreteLaplacian_Squared       !   note: returns ( ∑_j k_ij f_j )²

        public          ::      getGai
        public          ::      getFi      
        public          ::      getdGaidfbj
        public          ::      getdFidfbj
        public          ::      getd2Gaidfbjdfcj
        public          ::      getd2Fidfbjdfck

        public          ::      getKernel3x3
        public          ::      getKernel3x3x3
        public          ::      getKernel5x5
        public          ::      getKernel5x5x5
        public          ::      getSQ_Kernel5x5                 !   note: the squared kernel _derived_ from 3x3  k(2)_ij = ∑_l k_il k_lj is of size 5x5.
        public          ::      getSQ_Kernel5x5x5
        public          ::      getSQ_Kernel7x7
        public          ::      getSQ_Kernel7x7x7
        public          ::      getSQ_Kernel9x9
        public          ::      getSQ_Kernel9x9x9



        interface       discreteLaplacian
            module procedure    discreteLaplacian2d
            module procedure    discreteLaplacian3d
        end interface

        interface       discreteBiharmonic
            module procedure    discreteBiharmonic2d
            module procedure    discreteBiharmonic3d
        end interface
 
        interface       discreteLaplacian_Squared
            module procedure    discreteLaplacian_Squared2d
            module procedure    discreteLaplacian_Squared3d
        end interface

        interface       discreteLaplacian_deriv
            module procedure    discreteLaplacian_deriv2d
            module procedure    discreteLaplacian_deriv3d
        end interface
 
 

    !---    advanced laplacian kernels

        interface       getGai
            module procedure    getGai3d
        end interface

        interface       getFi
            module procedure    getFi3d
        end interface
         
        interface       getdGaidfbj
            module procedure    getdGaidfbj3d
        end interface
            
        interface       getdFidfbj
            module procedure    getdFidfbj3d
        end interface
            
        interface       getd2Gaidfbjdfcj
            module procedure    getd2Gaidfbjdfcj3d
        end interface

        interface       getd2Fidfbjdfck
            module procedure    getd2Fidfbjdfck3d
        end interface
     
      
 
        
    contains
!---^^^^^^^^

        pure real(kind=real64) function discreteLaplacian2d( f,noisy )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete laplacian in 2d assuming unit rectangular mesh
    !*      optionally switch to high accuracy solver using noisy = .false.
            real(kind=real64),dimension(:,:),intent(in)         ::      f
            logical,intent(in),optional                         ::      noisy
            if (present(noisy)) then
                if (size(f,dim=1)==3) then
                    discreteLaplacian2d = discreteLaplacian_33_2d( f,noisy )
                else
                    discreteLaplacian2d = discreteLaplacian_55_2d( f,noisy )
                end if
            else
                if (size(f,dim=1)==3) then
                    discreteLaplacian2d = discreteLaplacian_33_2d( f,.true. )
                else
                    discreteLaplacian2d = discreteLaplacian_55_2d( f,.true. )
                end if
            end if
            return
        end function discreteLaplacian2d

        pure real(kind=real64) function discreteLaplacian3d( f,noisy )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete laplacian in 3d assuming unit cubic mesh
    !*      optionally switch to high accuracy solver using noisy = .false.        
            real(kind=real64),dimension(:,:,:),intent(in)           ::      f
            logical,intent(in),optional                             ::      noisy
            if (present(noisy)) then
                if (size(f,dim=1)==3) then
                    discreteLaplacian3d = discreteLaplacian_333_3d( f,noisy )
                else
                    discreteLaplacian3d = discreteLaplacian_555_3d( f,noisy )
                end if
            else
                if (size(f,dim=1)==3) then
                    discreteLaplacian3d = discreteLaplacian_333_3d( f,.true. )
                else
                    discreteLaplacian3d = discreteLaplacian_555_3d( f,.true. )
                end if
            end if
            ! if (size(f,dim=1)==3) then
            !     discreteLaplacian3d = discreteLaplacian_333_3d( f )
            ! else if (present(noisy)) then
            !     discreteLaplacian3d = discreteLaplacian_555_3d( f,noisy )
            ! else 
            !     discreteLaplacian3d = discreteLaplacian_555_3d( f,.true. )                
            ! end if            
            return
        end function discreteLaplacian3d


        pure real(kind=real64) function discreteLaplacian_33_2d( f,noisy  )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete laplacian in 2d assuming unit rectangular mesh        
            real(kind=real64),dimension(-1:,-1:),intent(in)         ::      f
            logical,intent(in)                                      ::      noisy
            if (noisy) then
                discreteLaplacian_33_2d = sum( NOISY_KERNEL_33_2D(-1:1,-1:1)*f(-1:1,-1:1) )
            else
                discreteLaplacian_33_2d = sum( KERNEL_33_2D(-1:1,-1:1)*f(-1:1,-1:1) )
            end if
            return
        end function discreteLaplacian_33_2d

        pure real(kind=real64) function discreteLaplacian_333_3d( f,noisy  )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete laplacian in 3d assuming unit cubic mesh
            real(kind=real64),dimension(-1:,-1:,-1:),intent(in)         ::      f
            logical,intent(in)                                          ::      noisy
            !discreteLaplacian_333_3d = sum( KERNEL_333_3D*f )
            if (noisy) then
                discreteLaplacian_333_3d = sum( NOISY_KERNEL_333_3D(-1:1,-1:1,-1:1)*f(-1:1,-1:1,-1:1) )
            else
                discreteLaplacian_333_3d = sum( KERNEL_333_3D(-1:1,-1:1,-1:1)*f(-1:1,-1:1,-1:1) )
            end if
            return
        end function discreteLaplacian_333_3d



        pure real(kind=real64) function discreteLaplacian_55_2d( f,noisy )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete laplacian in 2d assuming unit rectangular mesh
    !*      optionally switch to high accuracy solver using noisy = .false.        
            real(kind=real64),dimension(-2:,-2:),intent(in)         ::      f
            logical,intent(in)                                      ::      noisy
            if (noisy) then
                discreteLaplacian_55_2d = sum( NOISY_KERNEL_55_2D(-2:2,-2:2)*f(-2:2,-2:2) )
            else
                discreteLaplacian_55_2d = sum( KERNEL_55_2D(-2:2,-2:2)*f(-2:2,-2:2) )
            end if
            return
        end function discreteLaplacian_55_2d

        pure real(kind=real64) function discreteLaplacian_555_3d( f,noisy )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete laplacian in 3d assuming unit cubic mesh
    !*      optionally switch to high accuracy solver using noisy = .false.        
            real(kind=real64),dimension(-2:,-2:,-2:),intent(in)         ::      f
            logical,intent(in)                                          ::      noisy
            if (noisy) then
                discreteLaplacian_555_3d = sum( NOISY_KERNEL_555_3D(-2:2,-2:2,-2:2)*f(-2:2,-2:2,-2:2) )
            else
                discreteLaplacian_555_3d = sum( KERNEL_555_3D(-2:2,-2:2,-2:2)*f(-2:2,-2:2,-2:2) )
            end if
            return
        end function discreteLaplacian_555_3d






        pure real(kind=real64) function discreteBiharmonic2d( f ) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete Biharmonic in 2d assuming unit rectangular mesh
            real(kind=real64),dimension(:,:),intent(in)         ::      f

            if (size(f,dim=1)==5) then
                discreteBiharmonic2d = discreteBiharmonic_55_2d( f )
            else if (size(f,dim=1)==7) then
                discreteBiharmonic2d = discreteBiharmonic_77_2d( f )
            else
                discreteBiharmonic2d = discreteBiharmonic_99_2d( f )
            end if
            return
        end function discreteBiharmonic2d

        pure real(kind=real64) function discreteBiharmonic3d( f ) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete Square laplacian in 3d assuming unit cubic mesh
            real(kind=real64),dimension(:,:,:),intent(in)         ::      f

            if (size(f,dim=1)==5) then
                discreteBiharmonic3d = discreteBiharmonic_555_3d( f )
            else if (size(f,dim=1)==7) then
                discreteBiharmonic3d = discreteBiharmonic_777_3d( f )
            else  
                discreteBiharmonic3d = discreteBiharmonic_999_3d( f )                
            end if            
            return
        end function discreteBiharmonic3d


        pure real(kind=real64) function discreteBiharmonic_55_2d( f )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete Square laplacian in 2d assuming unit rectangular mesh
            real(kind=real64),dimension(-2:,-2:),intent(in)         ::      f
            discreteBiharmonic_55_2d = sum( SQ_KERNEL_55_2D(-2:2,-2:2)*f(-2:2,-2:2) )
            return
        end function discreteBiharmonic_55_2d


        pure real(kind=real64) function discreteBiharmonic_77_2d( f )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete Square laplacian in 2d assuming unit rectangular mesh
            real(kind=real64),dimension(-3:,-3:),intent(in)         ::      f
            discreteBiharmonic_77_2d = sum( SQ_KERNEL_77_2D(-3:3,-3:3)*f(-3:3,-3:3) )
            return
        end function discreteBiharmonic_77_2d        

        pure real(kind=real64) function discreteBiharmonic_555_3d( f ) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete Square laplacian in 3d assuming unit cubic mesh
            real(kind=real64),dimension(-2:,-2:,-2:),intent(in)         ::      f
 
            discreteBiharmonic_555_3d = sum( SQ_KERNEL_555_3D(-2:2,-2:2,-2:2)*f(-2:2,-2:2,-2:2) )
            return
        end function discreteBiharmonic_555_3d



        pure real(kind=real64) function discreteBiharmonic_99_2d( f ) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete Square laplacian in 2d assuming unit rectangular mesh
            real(kind=real64),dimension(-4:,-4:),intent(in)         ::      f

            ! integer         ::      ii,jj
            ! discreteBiharmonic_55_2d = 0
            ! do jj = -4,4
            !     do ii = -4,4
            !         discreteBiharmonic_55_2d = discreteBiharmonic_55_2d + SQ_KERNEL_55_2D(ii,jj)*f(ii,jj)
            !     end do
            ! end do
            discreteBiharmonic_99_2d = sum( SQ_KERNEL_99_2D(-4:4,-4:4)*f(-4:4,-4:4) )
            return
        end function discreteBiharmonic_99_2d

        pure real(kind=real64) function discreteBiharmonic_777_3d( f ) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete Square laplacian in 3d assuming unit cubic mesh
            real(kind=real64),dimension(-3:,-3:,-3:),intent(in)         ::      f
 
            discreteBiharmonic_777_3d = sum( SQ_KERNEL_777_3D(-3:3,-3:3,-3:3)*f(-3:3,-3:3,-3:3) )
            return
        end function discreteBiharmonic_777_3d

        pure real(kind=real64) function discreteBiharmonic_999_3d( f ) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete Square laplacian in 3d assuming unit cubic mesh
            real(kind=real64),dimension(-4:,-4:,-4:),intent(in)         ::      f
 
            discreteBiharmonic_999_3d = sum( SQ_KERNEL_999_3D(-4:4,-4:4,-4:4)*f(-4:4,-4:4,-4:4) )
            return
        end function discreteBiharmonic_999_3d






        pure real(kind=real64) function discreteLaplacian_squared2d( f,noisy )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete laplacian in 2d assuming unit cubic mesh
            real(kind=real64),dimension(-1:,-1:),intent(in)         ::      f
            logical,intent(in),optional                         ::      noisy
            if (present(noisy)) then
                if (size(f,dim=1)==3) then
                    discreteLaplacian_squared2d = discreteLaplacian_33_2d( f,noisy ) 
                else 
                    discreteLaplacian_squared2d = discreteLaplacian_55_2d( f,noisy )             
                end if
            else
                if (size(f,dim=1)==3) then
                    discreteLaplacian_squared2d = discreteLaplacian_33_2d( f,.true. ) 
                else 
                    discreteLaplacian_squared2d = discreteLaplacian_55_2d( f,.true. )             
                end if
            end if
            ! if (size(f,dim=1)==3) then
            !     discreteLaplacian_squared2d = discreteLaplacian_33_2d( f ) 
            ! else if (present(noisy)) then
            !     discreteLaplacian_squared2d = discreteLaplacian_55_2d( f,noisy ) 
            ! else
            !     discreteLaplacian_squared2d = discreteLaplacian_55_2d( f,.true. ) 
            ! end if
            discreteLaplacian_squared2d = discreteLaplacian_squared2d*discreteLaplacian_squared2d
            return
        end function discreteLaplacian_squared2d

        pure real(kind=real64) function discreteLaplacian_squared3d( f,noisy )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the discrete laplacian in 3d assuming unit cubic mesh
            real(kind=real64),dimension(-1:,-1:,-1:),intent(in)         ::      f
            logical,intent(in),optional                         ::      noisy
            if (present(noisy)) then
                if (size(f,dim=1)==3) then
                    discreteLaplacian_squared3d = discreteLaplacian_333_3d( f,noisy ) 
                else 
                    discreteLaplacian_squared3d = discreteLaplacian_555_3d( f,noisy )             
                end if
            else
                if (size(f,dim=1)==3) then
                    discreteLaplacian_squared3d = discreteLaplacian_333_3d( f,.true. ) 
                else 
                    discreteLaplacian_squared3d = discreteLaplacian_555_3d( f,.true. )             
                end if
            end if
            discreteLaplacian_squared3d = discreteLaplacian_squared3d*discreteLaplacian_squared3d
            return
        end function discreteLaplacian_squared3d


        pure real(kind=real64) function discreteLaplacian_deriv2d( i,j )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the derivative of the discrete laplacian in 2d assuming unit rectangular mesh
            integer,intent(in)                  ::      i,j
            discreteLaplacian_deriv2d = KERNEL_33_2D(i,j)
            return
        end function discreteLaplacian_deriv2d

        pure real(kind=real64) function discreteLaplacian_deriv3d( i,j,k )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the derivative of the discrete laplacian in 3d assuming unit cubic mesh
            integer,intent(in)                  ::      i,j,k
            discreteLaplacian_deriv3d = KERNEL_333_3D(i,j,k)
            return
        end function discreteLaplacian_deriv3d

 

  

        pure real(kind=real64) function getGai3d( f,a ) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute G_ai = sum_{j ∈ N_i} ( kappa_ij f_aj/f_j )
    !*      With many of these routines I need the normalising factors 1 / sum_(a>=1) f_aj
    !*      so on input f(0,j) = 1 / sum_(a>=1) f_aj, with f(0,j)=0 if sum_a f_aj = 0 to avoid conditional
            real(kind=real64),dimension(0:,-1:,-1:,-1:),intent(in)          ::      f
            integer,intent(in)                                              ::      a
            
            integer             ::      jx,jy,jz
            getGai3d = 0
            do jz = -1,1
                do jy = -1,1
                    do jx = -1,1
                        getGai3d = getGai3d + KERNEL_333_3D(jx,jy,jz) * f(a,jx,jy,jz) * f(LIB_LK_INDX_INVSUM,jx,jy,jz)
                    end do
                end do
            end do
            

            return
        end function getGai3d


            

        pure real(kind=real64) function getFi3d( f ) 
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute F_i                 = 1/(2 f_i) sum_a ( G_ia )^2       
    !*      With many of these routines I need the normalising factors 1 / sum_(a>=1) f_aj
    !*      so on input f(0,j) = 1 / sum_(a>=1) f_aj, with f(0,j)=0 if sum_a f_aj = 0 to avoid conditional
            real(kind=real64),dimension(0:,-1:,-1:,-1:),intent(in)          ::      f            
            integer             ::      Na
            integer             ::      aa
            real(kind=real64)   ::      Gai
            
            Na = size(f,dim=1)-1
            getFi3d = 0
            do aa = 1,Na
                Gai = getGai3d( f,aa ) 
                getFi3d = getFi3d + Gai*Gai
            end do
            getFi3d = getFi3d * f(LIB_LK_INDX_INVSUM,0,0,0) / 2
            return
        end function getFi3d
     


        pure function getdGaidfbj3d( f,a ) result(dGaidfbj)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      d G_ai / d f_bj     = (kappa_ij/f_j) ( delta_ab - f_aj/f_j ) 
    !*      With many of these routines I need the normalising factors 1 / sum_(a>=1) f_aj
    !*      so on input f(0,j) = 1 / sum_(a>=1) f_aj, with f(0,j)=0 if sum_a f_aj = 0 to avoid conditional
            real(kind=real64),dimension(0:,-1:,-1:,-1:),intent(in)          ::      f
            integer,intent(in)                                              ::      a
            real(kind=real64),dimension(size(f,dim=1)-1,-1:1,-1:1,-1:1)     ::      dGaidfbj
            integer             ::      Na
            integer             ::      bb
            integer             ::      jx,jy,jz
            Na = size(f,dim=1)-1
                         
            do jz = -1,1
                do jy = -1,1
                    do jx = -1,1
                        do bb = 1,Na
                            dGaidfbj(bb,jx,jy,jz) = - f(a,jx,jy,jz)*f(LIB_LK_INDX_INVSUM,jx,jy,jz)                            
                        end do
                        dGaidfbj(a,jx,jy,jz) = dGaidfbj(a,jx,jy,jz) + 1
                        do bb = 1,Na
                            dGaidfbj(bb,jx,jy,jz) = dGaidfbj(bb,jx,jy,jz) * KERNEL_333_3D(jx,jy,jz) * f(LIB_LK_INDX_INVSUM,jx,jy,jz)                    
                        end do
                    end do
                end do
            end do


            return
        end function getdGaidfbj3d
     


        pure function getdFidfbj3d( f ) result(dFifbj)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      d F_i / d f_bj      = ( sum_a G_ai d G_ai / d f_bj  - delta_ij F_j ) / f_i
    !*      With many of these routines I need the normalising factors 1 / sum_(a>=1) f_aj
    !*      so on input f(0,j) = 1 / sum_(a>=1) f_aj, with f(0,j)=0 if sum_a f_aj = 0 to avoid conditional
            real(kind=real64),dimension(0:,-1:,-1:,-1:),intent(in)              ::      f
            real(kind=real64),dimension(size(f,dim=1)-1,-1:1,-1:1,-1:1)         ::      dFifbj
            integer             ::      Na
            integer             ::      aa
            real(kind=real64),dimension(size(f,dim=1)-1,-1:1,-1:1,-1:1)        ::      dGaidfbj
            real(kind=real64)   ::      Gai,Fi,f_j

            Na = size(f,dim=1)-1
            dFifbj = 0
            Fi = 0
            f_j = f(LIB_LK_INDX_INVSUM,0,0,0)
            do aa = 1,Na
                Gai = getGai3d( f,aa ) 
                dGaidfbj = getdGaidfbj3d( f,aa )  
                dFifbj(:,:,:,:) = dFifbj(:,:,:,:) + Gai*dGaidfbj(:,:,:,:)  
                Fi = Fi + Gai*Gai
            end do
            Fi = Fi * f_j/2
             
            dFifbj(:,0,0,0) = dFifbj(:,0,0,0) - Fi
            
            dFifbj = dFifbj*f_j

            return
        end function getdFidfbj3d
     



        pure function getd2Gaidfbjdfcj3d( f,a ) result(d2Gaidfbjdfcj)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      d² G_ai / d f_bj d f_ck     = (delta_jk kappa_ij / rho_j^2) ( 2 f_aj/f_j - delta_ab - delta_ac )
    !*      With many of these routines I need the normalising factors 1 / sum_(a>=1) f_aj
    !*      so on input f(0,j) = 1 / sum_(a>=1) f_aj, with f(0,j)=0 if sum_a f_aj = 0 to avoid conditional
    !*      note that the solution includes delta_jk, so in fact this function returns
    !*      d² G_ai / d f_bj d f_cj
            real(kind=real64),dimension(0:,-1:,-1:,-1:),intent(in)          ::      f
            integer,intent(in)                                              ::      a
            real(kind=real64),dimension(size(f,dim=1)-1,size(f,dim=1)-1,-1:1,-1:1,-1:1)        ::      d2Gaidfbjdfcj

            integer             ::      Na
            integer             ::      bb
            integer             ::      jx,jy,jz
            real(kind=real64)   ::      f_j
            real(kind=real64),dimension(size(f,dim=1)-1,size(f,dim=1)-1)    ::      d2G
            Na = size(f,dim=1)-1
            
            do jz = -1,1
                do jy = -1,1
                    do jx = -1,1
                        f_j = f(LIB_LK_INDX_INVSUM,jx,jy,jz) 
                        d2G(:,:) = 2*f(a,jx,jy,jz)*f_j

                        do bb = 1,Na
                            d2G(a,bb) = d2G(a,bb) - 1
                            d2G(bb,a) = d2G(bb,a) - 1
                        end do
                        d2G(:,:) = d2G(:,:) * KERNEL_333_3D(jx,jy,jz) * f_j*f_j


                        d2Gaidfbjdfcj(:,:,jx,jy,jz) = d2G(:,:)
                    end do
                end do
            end do
 
            return
        end function getd2Gaidfbjdfcj3d
     

        pure function getd2Fidfbjdfck3d( f ) result(d2Fidfbjdfck)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      d² F_i / d f_bj d f_ck      = ( sum_a d G_ai / d f_bj  d G_ai / d f_ck + sum_a G_ai  d² G_ai / d f_bj d f_ck -( delta_ij  d F_i / d f_bj +  delta_ik d F_i / d f_ck )/f_i
    !*      With many of these routines I need the normalising factors 1 / sum_(a>=1) f_aj
    !*      so on input f(0,j) = 1 / sum_(a>=1) f_aj, with f(0,j)=0 if sum_a f_aj = 0 to avoid conditional
            real(kind=real64),dimension(0:,-1:,-1:,-1:),intent(in)                      ::      f
            real(kind=real64),dimension(size(f,dim=1)-1,-1:1,-1:1,-1:1,     &
                                        size(f,dim=1)-1,-1:1,-1:1,-1:1)                 ::      d2Fidfbjdfck

            integer             ::      Na
            integer             ::      aa,bb,cc 
            real(kind=real64),dimension(size(f,dim=1)-1,-1:1,-1:1,-1:1)                 ::      dGaidfbj
            real(kind=real64),dimension(size(f,dim=1)-1,size(f,dim=1)-1,-1:1,-1:1,-1:1) ::      d2Gaidfbjdfcj
            real(kind=real64),dimension(size(f,dim=1)-1,-1:1,-1:1,-1:1)                 ::      dFifbj
            real(kind=real64)   ::      Gai 
            integer             ::      jx,jy,jz
            integer             ::      kx,ky,kz
            
            Na = size(f,dim=1)-1            
            d2Fidfbjdfck = 0             
            dFifbj(:,:,:,:) = getdFidfbj3d( f )

            do aa = 1,Na
                Gai = getGai3d(f,aa)
                d2Gaidfbjdfcj = getd2Gaidfbjdfcj3d( f,aa )
                dGaidfbj(:,:,:,:) = getdGaidfbj3d( f,aa )
                

                do kz = -1,1
                    do ky = -1,1
                        do kx = -1,1
                            do cc = 1,Na                                
                                do bb = 1,Na                            
                                    d2Fidfbjdfck(bb,0,0,0,cc,kx,ky,kz)    = d2Fidfbjdfck(bb,0,0,0,cc,kx,ky,kz)    - dFifbj(cc,kx,ky,kz)
                                    d2Fidfbjdfck(bb,kx,ky,kz,cc,0,0,0)    = d2Fidfbjdfck(bb,kx,ky,kz,cc,0,0,0)    - dFifbj(bb,kx,ky,kz)
                                    d2Fidfbjdfck(bb,kx,ky,kz,cc,kx,ky,kz) = d2Fidfbjdfck(bb,kx,ky,kz,cc,kx,ky,kz) + Gai*d2Gaidfbjdfcj(bb,cc,kx,ky,kz)
                                end do
                            end do
                            do cc = 1,Na
                                do jz = -1,1
                                    do jy = -1,1
                                        do jx = -1,1                                        
                                            do bb = 1,Na                
                                                d2Fidfbjdfck(bb,jx,jy,jz,cc,kx,ky,kz) = d2Fidfbjdfck(bb,jx,jy,jz,cc,kx,ky,kz) + dGaidfbj(bb,jx,jy,jz)*dGaidfbj(cc,kx,ky,kz)
                                            end do
                                        end do
                                    end do
                                end do
                            end do                                                  
                        end do
                    end do
                end do
            end do
            d2Fidfbjdfck = d2Fidfbjdfck * f(LIB_LK_INDX_INVSUM,0,0,0)
              
            return
        end function getd2Fidfbjdfck3d


        pure function getKernel3x3(noisy) result(kernel)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(3,3)        ::      kernel
            logical,intent(in),optional             ::      noisy
            if (present(noisy)) then
                if (noisy) then
                    kernel = NOISY_KERNEL_33_2D
                else
                    kernel = KERNEL_33_2D
                end if
            else
                kernel = KERNEL_33_2D
            end if

            return
        end function getKernel3x3  
     
        pure function getKernel3x3x3(noisy) result(kernel)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(3,3,3)      ::      kernel
            logical,intent(in),optional             ::      noisy
            if (present(noisy)) then
                if (noisy) then
                    kernel = NOISY_KERNEL_333_3D
                else
                    kernel = KERNEL_333_3D
                end if
            else
                kernel = KERNEL_333_3D
            end if
        
            return
        end function getKernel3x3x3
     
        pure function getKernel5x5(noisy) result(kernel)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(5,5)        ::      kernel
            logical,intent(in),optional             ::      noisy
            if (present(noisy)) then
                if (noisy) then
                    kernel = NOISY_KERNEL_55_2D
                else
                    kernel = KERNEL_55_2D
                end if
            else
                kernel = KERNEL_55_2D
            end if
            return
        end function getKernel5x5  
     
        pure function getKernel5x5x5(noisy) result(kernel)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(5,5,5)      ::      kernel
            logical,intent(in),optional             ::      noisy
            if (present(noisy)) then
                if (noisy) then
                    kernel = NOISY_KERNEL_555_3D
                else
                    kernel = KERNEL_555_3D
                end if
            else
                kernel = KERNEL_555_3D
            end if
            
            return
        end function getKernel5x5x5  
     
        pure function getSQ_Kernel5x5() result(kernel)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(5,5)        ::      kernel
            kernel = SQ_KERNEL_55_2D
            return
        end function getSQ_Kernel5x5  
     
        pure function getSQ_Kernel5x5x5() result(kernel)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(5,5,5)      ::      kernel
            kernel = SQ_KERNEL_555_3D
            return
        end function getSQ_Kernel5x5x5
     
        pure function getSQ_Kernel7x7() result(kernel)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(7,7)        ::      kernel
            kernel = SQ_KERNEL_77_2D
            return
        end function getSQ_Kernel7x7  
     
        pure function getSQ_Kernel7x7x7() result(kernel)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(7,7,7)      ::      kernel
            kernel = SQ_KERNEL_777_3D
            return
        end function getSQ_Kernel7x7x7  


        pure function getSQ_Kernel9x9() result(kernel)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(9,9)        ::      kernel
            kernel = SQ_KERNEL_99_2D
            return
        end function getSQ_Kernel9x9
     
        pure function getSQ_Kernel9x9x9() result(kernel)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(9,9,9)      ::      kernel
            kernel = SQ_KERNEL_999_3D
            return
        end function getSQ_Kernel9x9x9

     
    end module Lib_LaplacianKernel


