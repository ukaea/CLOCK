!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_LaplacianKernel_data from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*      supplemental data file containing the kernels

        real(kind=real64),dimension(-1:1,-1:1,-1:1),public,parameter           ::      POLYF_333_3D = reshape(   (/                                      &
                                        -2, 1,-2,                                                   &
                                         1, 4, 1,                                                   &
                                        -2, 1,-2,                                                   &
                                         1, 4, 1,                                                   &
                                         4, 7, 4,                                                   &
                                         1, 4, 1,                                                   &
                                        -2, 1,-2,                                                   &
                                         1, 4, 1,                                                   &
                                        -2, 1,-2    /),(/3,3,3/) ) / 27.0d0

        real(kind=real64),dimension(-1:1,-1:1,-1:1),public,parameter           ::      NOISY_KERNEL_333_3D = reshape(   (/                                      &
                                         1, 0, 1,                                                   &
                                         0,-1, 0,                                                   &
                                         1, 0, 1,                                                   &
                                         0,-1, 0,                                                   &
                                        -1,-2,-1,                                                   &
                                         0,-1, 0,                                                   &
                                         1, 0, 1,                                                   &
                                         0,-1, 0,                                                   &
                                         1, 0, 1    /),(/3,3,3/) ) / 3.0d0
     

          real(kind=real64),dimension(-1:1,-1:1,-1:1),public,parameter           ::      KERNEL_333_3D = reshape(   (/                                      &
                                        -1,  2, -1,                                                   &
                                         2, -1,  2,                                                   &
                                        -1,  2, -1,                                                   &
                                         2, -1,  2,                                                   &
                                        -1,-10, -1,                                                   &
                                         2, -1,  2,                                                   &
                                        -1,  2, -1,                                                   &
                                         2, -1,  2,                                                   &
                                        -1,  2, -1    /),(/3,3,3/) ) / 3.0d0


 

        real(kind=real64),dimension(-1:1,-1:1),private,parameter                ::      KERNEL_33_2D = reshape(   (/                                      &
                            1,2,1,  2,-12,2,    1,2,1    /),(/3,3/) ) / 4.0d0


        real(kind=real64),dimension(-1:1,-1:1),private,parameter                ::      NOISY_KERNEL_33_2D = reshape(   (/                                      &
                            24,-11,24,     -11,-52,-11,     24,-11,24     /),(/3,3/) ) / 37.0d0



    !---    these magic numbers are defined by setting sum k_ij g_ij = nabla^2 f
    !       with f = f0 + f" x^2/2 + f""" x^4/24                            ( note, I don't need to worry about f' terms for a symmetric kernel )
    !       and then finding the minimum value for s = sum kij kij                                        
        real(kind=real64),private,parameter     ::      k00 = -3420.0d0/2940
        real(kind=real64),private,parameter     ::      k01 =  -926.0d0/2940
        real(kind=real64),private,parameter     ::      k11 =  1448.0d0/2940
        real(kind=real64),private,parameter     ::      k02 = -1039.0d0/2940
        real(kind=real64),private,parameter     ::      k12 =   975.0d0/2940
        real(kind=real64),private,parameter     ::      k22 =  -578.0d0/2940
          

        real(kind=real64),private,parameter     ::      noisy_k00 = -4d0/35
        real(kind=real64),private,parameter     ::      noisy_k01 = -3d0/35
        real(kind=real64),private,parameter     ::      noisy_k11 = -2d0/35
        real(kind=real64),private,parameter     ::      noisy_k02 = 0
        real(kind=real64),private,parameter     ::      noisy_k12 = 1d0/35
        real(kind=real64),private,parameter     ::      noisy_k22 = 4d0/35
                            

        real(kind=real64),dimension(-2:2,-2:2),private,parameter                ::      KERNEL_55_2D = reshape( (/           &
                                                k22,k12,k02,k12,k22,                                                         &
                                                k12,k11,k01,k11,k12,                                                         &
                                                k02,k01,k00,k01,k02,                                                         &
                                                k12,k11,k01,k11,k12,                                                         &
                                                k22,k12,k02,k12,k22                               /) , (/5,5/) ) 

        real(kind=real64),dimension(-2:2,-2:2),private,parameter                ::      NOISY_KERNEL_55_2D = reshape( (/     &
                                                noisy_k22,noisy_k12,noisy_k02,noisy_k12,noisy_k22,                           &
                                                noisy_k12,noisy_k11,noisy_k01,noisy_k11,noisy_k12,                           &
                                                noisy_k02,noisy_k01,noisy_k00,noisy_k01,noisy_k02,                           &
                                                noisy_k12,noisy_k11,noisy_k01,noisy_k11,noisy_k12,                           &
                                                noisy_k22,noisy_k12,noisy_k02,noisy_k12,noisy_k22                               /) , (/5,5/) ) 

                                    
     !---     values required for a high accuracy kernel
         real(kind=real64),private,parameter     ::      k000 = -39d0/98
         real(kind=real64),private,parameter     ::      k001 = -779d0/3675
         real(kind=real64),private,parameter     ::      k011 = -251d0/7350
         real(kind=real64),private,parameter     ::      k111 = 166d0/1225
         real(kind=real64),private,parameter     ::      k002 = -2509d0/14700
         real(kind=real64),private,parameter     ::      k012 = -17d0/980
         real(kind=real64),private,parameter     ::      k112 =  1879d0/14700
         real(kind=real64),private,parameter     ::      k022 = -272d0/3675
         real(kind=real64),private,parameter     ::      k122 = 7d0/150
         real(kind=real64),private,parameter     ::      k222 = -529d0/4900


     !--- values required for a kernel expected to work on slightly noisy data.
     !    note: these are slightly better performing for smoothly varying data than those below, but exaggerate effect of corners
     !    returns (fxx+fyy+fzz) + 0.369 (fxxxx+fyyyy+fzzzz) + 2 (fxxyy+fyyzz+fzzxx)
     !    real(kind=real64),private,parameter     ::      noisy_k000 = -6d0/175
     !    real(kind=real64),private,parameter     ::      noisy_k001 = -5d0/175
     !    real(kind=real64),private,parameter     ::      noisy_k011 = -4d0/175 
     !    real(kind=real64),private,parameter     ::      noisy_k111 = -3d0/175
     !    real(kind=real64),private,parameter     ::      noisy_k002 = -2d0/175
     !    real(kind=real64),private,parameter     ::      noisy_k012 = -1d0/175
     !    real(kind=real64),private,parameter     ::      noisy_k112 = 0
     !    real(kind=real64),private,parameter     ::      noisy_k022 = 2d0/175
     !    real(kind=real64),private,parameter     ::      noisy_k122 = 3d0/175
     !    real(kind=real64),private,parameter     ::      noisy_k222 = 6d0/175


     !--- values required for a kernel expected to work on slightly noisy data
     !     returns (fxx+fyy+fzz) + 0.363 (fxxxx+fyyyy+fzzzz) + 1.585 (fxxyy+fyyzz+fzzxx)
          real(kind=real64),private,parameter     ::      noisy_k000 = -218d0/4825
          real(kind=real64),private,parameter     ::      noisy_k001 = -179d0/4825
          real(kind=real64),private,parameter     ::      noisy_k011 = -140d0/4825
          real(kind=real64),private,parameter     ::      noisy_k111 = -101d0/4825
          real(kind=real64),private,parameter     ::      noisy_k002 = -62d0/4825
          real(kind=real64),private,parameter     ::      noisy_k012 = -23d0/4825
          real(kind=real64),private,parameter     ::      noisy_k112 = 16d0/4825
          real(kind=real64),private,parameter     ::      noisy_k022 = 94d0/4825
          real(kind=real64),private,parameter     ::      noisy_k122 = 133d0/4825
          real(kind=real64),private,parameter     ::      noisy_k222 = 0

         real(kind=real64),dimension(-2:2,-2:2,-2:2),public,parameter                ::      KERNEL_555_3D = reshape( (/                                    &
                                                 k222,k122,k022,k122,k222,                           &    
                                                 k122,k112,k012,k112,k122,                           &
                                                 k022,k012,k002,k012,k022,                           &
                                                 k122,k112,k012,k112,k122,                           &
                                                 k222,k122,k022,k122,k222,                           &
                                                 k122,k112,k012,k112,k122,                           &
                                                 k112,k111,k011,k111,k112,                           &
                                                 k012,k011,k001,k011,k012,                           &
                                                 k112,k111,k011,k111,k112,                           &
                                                 k122,k112,k012,k112,k122,                           &
                                                 k022,k012,k002,k012,k022,                           &
                                                 k012,k011,k001,k011,k012,                           &
                                                 k002,k001,k000,k001,k002,                           &
                                                 k012,k011,k001,k011,k012,                           &
                                                 k022,k012,k002,k012,k022,                           &
                                                 k122,k112,k012,k112,k122,                           &
                                                 k112,k111,k011,k111,k112,                           &
                                                 k012,k011,k001,k011,k012,                           &
                                                 k112,k111,k011,k111,k112,                           &
                                                 k122,k112,k012,k112,k122,                           &
                                                 k222,k122,k022,k122,k222,                           &
                                                 k122,k112,k012,k112,k122,                           &
                                                 k022,k012,k002,k012,k022,                           &
                                                 k122,k112,k012,k112,k122,                           &
                                                 k222,k122,k022,k122,k222       /) , (/5,5,5/) )  


         real(kind=real64),dimension(-2:2,-2:2,-2:2),public,parameter                ::      NOISY_KERNEL_555_3D = reshape( (/     &  
                                                 noisy_k222,noisy_k122,noisy_k022,noisy_k122,noisy_k222,                           &    
                                                 noisy_k122,noisy_k112,noisy_k012,noisy_k112,noisy_k122,                           &
                                                 noisy_k022,noisy_k012,noisy_k002,noisy_k012,noisy_k022,                           &
                                                 noisy_k122,noisy_k112,noisy_k012,noisy_k112,noisy_k122,                           &
                                                 noisy_k222,noisy_k122,noisy_k022,noisy_k122,noisy_k222,                           &
                                                 noisy_k122,noisy_k112,noisy_k012,noisy_k112,noisy_k122,                           &
                                                 noisy_k112,noisy_k111,noisy_k011,noisy_k111,noisy_k112,                           &
                                                 noisy_k012,noisy_k011,noisy_k001,noisy_k011,noisy_k012,                           &
                                                 noisy_k112,noisy_k111,noisy_k011,noisy_k111,noisy_k112,                           &
                                                 noisy_k122,noisy_k112,noisy_k012,noisy_k112,noisy_k122,                           &
                                                 noisy_k022,noisy_k012,noisy_k002,noisy_k012,noisy_k022,                           &
                                                 noisy_k012,noisy_k011,noisy_k001,noisy_k011,noisy_k012,                           &
                                                 noisy_k002,noisy_k001,noisy_k000,noisy_k001,noisy_k002,                           &
                                                 noisy_k012,noisy_k011,noisy_k001,noisy_k011,noisy_k012,                           &
                                                 noisy_k022,noisy_k012,noisy_k002,noisy_k012,noisy_k022,                           &
                                                 noisy_k122,noisy_k112,noisy_k012,noisy_k112,noisy_k122,                           &
                                                 noisy_k112,noisy_k111,noisy_k011,noisy_k111,noisy_k112,                           &
                                                 noisy_k012,noisy_k011,noisy_k001,noisy_k011,noisy_k012,                           &
                                                 noisy_k112,noisy_k111,noisy_k011,noisy_k111,noisy_k112,                           &
                                                 noisy_k122,noisy_k112,noisy_k012,noisy_k112,noisy_k122,                           &
                                                 noisy_k222,noisy_k122,noisy_k022,noisy_k122,noisy_k222,                           &
                                                 noisy_k122,noisy_k112,noisy_k012,noisy_k112,noisy_k122,                           &
                                                 noisy_k022,noisy_k012,noisy_k002,noisy_k012,noisy_k022,                           &
                                                 noisy_k122,noisy_k112,noisy_k012,noisy_k112,noisy_k122,                           &
                                                 noisy_k222,noisy_k122,noisy_k022,noisy_k122,noisy_k222       /) , (/5,5,5/) )                                                   
 

     !--- values required for a biharmonic kernel 555 
     !    warning- poor noise control
          real(kind=real64),private,parameter     ::      b_555_000 = 1002d0/1225
          real(kind=real64),private,parameter     ::      b_555_001 =  472d0/1225
          real(kind=real64),private,parameter     ::      b_555_011 =  -48d0/1225
          real(kind=real64),private,parameter     ::      b_555_111 = -558d0/1225
          real(kind=real64),private,parameter     ::      b_555_002 =  597d0/1225
          real(kind=real64),private,parameter     ::      b_555_012 =  107d0/1225
          real(kind=real64),private,parameter     ::      b_555_112 = -373d0/1225
          real(kind=real64),private,parameter     ::      b_555_022 =  352d0/1225
          real(kind=real64),private,parameter     ::      b_555_122 =  -98d0/1225
          real(kind=real64),private,parameter     ::      b_555_222 =  267d0/1225



         real(kind=real64),dimension(-2:2,-2:2,-2:2),public,parameter                ::      SQ_KERNEL_555_3D = reshape( (/                                    &
                                                 b_555_222,b_555_122,b_555_022,b_555_122,b_555_222,                           &    
                                                 b_555_122,b_555_112,b_555_012,b_555_112,b_555_122,                           &
                                                 b_555_022,b_555_012,b_555_002,b_555_012,b_555_022,                           &
                                                 b_555_122,b_555_112,b_555_012,b_555_112,b_555_122,                           &
                                                 b_555_222,b_555_122,b_555_022,b_555_122,b_555_222,                           &
                                                 b_555_122,b_555_112,b_555_012,b_555_112,b_555_122,                           &
                                                 b_555_112,b_555_111,b_555_011,b_555_111,b_555_112,                           &
                                                 b_555_012,b_555_011,b_555_001,b_555_011,b_555_012,                           &
                                                 b_555_112,b_555_111,b_555_011,b_555_111,b_555_112,                           &
                                                 b_555_122,b_555_112,b_555_012,b_555_112,b_555_122,                           &
                                                 b_555_022,b_555_012,b_555_002,b_555_012,b_555_022,                           &
                                                 b_555_012,b_555_011,b_555_001,b_555_011,b_555_012,                           &
                                                 b_555_002,b_555_001,b_555_000,b_555_001,b_555_002,                           &
                                                 b_555_012,b_555_011,b_555_001,b_555_011,b_555_012,                           &
                                                 b_555_022,b_555_012,b_555_002,b_555_012,b_555_022,                           &
                                                 b_555_122,b_555_112,b_555_012,b_555_112,b_555_122,                           &
                                                 b_555_112,b_555_111,b_555_011,b_555_111,b_555_112,                           &
                                                 b_555_012,b_555_011,b_555_001,b_555_011,b_555_012,                           &
                                                 b_555_112,b_555_111,b_555_011,b_555_111,b_555_112,                           &
                                                 b_555_122,b_555_112,b_555_012,b_555_112,b_555_122,                           &
                                                 b_555_222,b_555_122,b_555_022,b_555_122,b_555_222,                           &
                                                 b_555_122,b_555_112,b_555_012,b_555_112,b_555_122,                           &
                                                 b_555_022,b_555_012,b_555_002,b_555_012,b_555_022,                           &
                                                 b_555_122,b_555_112,b_555_012,b_555_112,b_555_122,                           &
                                                 b_555_222,b_555_122,b_555_022,b_555_122,b_555_222       /) , (/5,5,5/) )  



     !--- values required for a biharmonic kernel 777 
     !    good noise control
          real(kind=real64),private,parameter     ::      b_777_000 =  2796d0/67914
          real(kind=real64),private,parameter     ::      b_777_001 =  2078d0/67914
          real(kind=real64),private,parameter     ::      b_777_011 =  1371d0/67914
          real(kind=real64),private,parameter     ::      b_777_111 =   675d0/67914
          real(kind=real64),private,parameter     ::      b_777_002 =   806d0/67914
          real(kind=real64),private,parameter     ::      b_777_012 =   132d0/67914
          real(kind=real64),private,parameter     ::      b_777_112 =  -531d0/67914
          real(kind=real64),private,parameter     ::      b_777_022 = -1008d0/67914
          real(kind=real64),private,parameter     ::      b_777_122 = -1638d0/67914
          real(kind=real64),private,parameter     ::      b_777_003 =  1626d0/67914
          real(kind=real64),private,parameter     ::      b_777_013 =  1007d0/67914
          real(kind=real64),private,parameter     ::      b_777_113 =   399d0/67914
          real(kind=real64),private,parameter     ::      b_777_222 = -2646d0/67914
          real(kind=real64),private,parameter     ::      b_777_023 =    32d0/67914
          real(kind=real64),private,parameter     ::      b_777_123 =  -543d0/67914
          real(kind=real64),private,parameter     ::      b_777_223 = -1386d0/67914
          real(kind=real64),private,parameter     ::      b_777_033 =  1347d0/67914
          real(kind=real64),private,parameter     ::      b_777_133 =   827d0/67914
          real(kind=real64),private,parameter     ::      b_777_233 =   149d0/67914
          real(kind=real64),private,parameter     ::      b_777_333 =  1959d0/67914
          


         real(kind=real64),dimension(-3:3,-3:3,-3:3),public,parameter                ::      SQ_KERNEL_777_3D = reshape( (/       &
                                                 b_777_333,b_777_233,b_777_133,b_777_033,b_777_133,b_777_233,b_777_333,           &
                                                 b_777_233,b_777_223,b_777_123,b_777_023,b_777_123,b_777_223,b_777_233,           &
                                                 b_777_133,b_777_123,b_777_113,b_777_013,b_777_113,b_777_123,b_777_133,           &
                                                 b_777_033,b_777_023,b_777_013,b_777_003,b_777_013,b_777_023,b_777_033,           &
                                                 b_777_133,b_777_123,b_777_113,b_777_013,b_777_113,b_777_123,b_777_133,           &
                                                 b_777_233,b_777_223,b_777_123,b_777_023,b_777_123,b_777_223,b_777_233,           &
                                                 b_777_333,b_777_233,b_777_133,b_777_033,b_777_133,b_777_233,b_777_333,           &
                                                 b_777_233,b_777_223,b_777_123,b_777_023,b_777_123,b_777_223,b_777_233,           &      
                                                 b_777_223,b_777_222,b_777_122,b_777_022,b_777_122,b_777_222,b_777_223,           &      
                                                 b_777_123,b_777_122,b_777_112,b_777_012,b_777_112,b_777_122,b_777_123,           &      
                                                 b_777_023,b_777_022,b_777_012,b_777_002,b_777_012,b_777_022,b_777_023,           &                      
                                                 b_777_123,b_777_122,b_777_112,b_777_012,b_777_112,b_777_122,b_777_123,           &      
                                                 b_777_223,b_777_222,b_777_122,b_777_022,b_777_122,b_777_222,b_777_223,           &      
                                                 b_777_233,b_777_223,b_777_123,b_777_023,b_777_123,b_777_223,b_777_233,           &      
                                                 b_777_133,b_777_123,b_777_113,b_777_013,b_777_113,b_777_123,b_777_133,           &     
                                                 b_777_123,b_777_122,b_777_112,b_777_012,b_777_112,b_777_122,b_777_123,           &     
                                                 b_777_113,b_777_112,b_777_111,b_777_011,b_777_111,b_777_112,b_777_113,           &     
                                                 b_777_013,b_777_012,b_777_011,b_777_001,b_777_011,b_777_012,b_777_013,           &      
                                                 b_777_113,b_777_112,b_777_111,b_777_011,b_777_111,b_777_112,b_777_113,           &     
                                                 b_777_123,b_777_122,b_777_112,b_777_012,b_777_112,b_777_122,b_777_123,           &     
                                                 b_777_133,b_777_123,b_777_113,b_777_013,b_777_113,b_777_123,b_777_133,           &     
                                                 b_777_033,b_777_023,b_777_013,b_777_003,b_777_013,b_777_023,b_777_033,           &     
                                                 b_777_023,b_777_022,b_777_012,b_777_002,b_777_012,b_777_022,b_777_023,           &     
                                                 b_777_013,b_777_012,b_777_011,b_777_001,b_777_011,b_777_012,b_777_013,           &     
                                                 b_777_003,b_777_002,b_777_001,b_777_000,b_777_001,b_777_002,b_777_003,           &     
                                                 b_777_013,b_777_012,b_777_011,b_777_001,b_777_011,b_777_012,b_777_013,           &     
                                                 b_777_023,b_777_022,b_777_012,b_777_002,b_777_012,b_777_022,b_777_023,           &     
                                                 b_777_033,b_777_023,b_777_013,b_777_003,b_777_013,b_777_023,b_777_033,           &     
                                                 b_777_133,b_777_123,b_777_113,b_777_013,b_777_113,b_777_123,b_777_133,           &     
                                                 b_777_123,b_777_122,b_777_112,b_777_012,b_777_112,b_777_122,b_777_123,           &     
                                                 b_777_113,b_777_112,b_777_111,b_777_011,b_777_111,b_777_112,b_777_113,           &     
                                                 b_777_013,b_777_012,b_777_011,b_777_001,b_777_011,b_777_012,b_777_013,           &      
                                                 b_777_113,b_777_112,b_777_111,b_777_011,b_777_111,b_777_112,b_777_113,           &     
                                                 b_777_123,b_777_122,b_777_112,b_777_012,b_777_112,b_777_122,b_777_123,           &     
                                                 b_777_133,b_777_123,b_777_113,b_777_013,b_777_113,b_777_123,b_777_133,           &     
                                                 b_777_233,b_777_223,b_777_123,b_777_023,b_777_123,b_777_223,b_777_233,           &      
                                                 b_777_223,b_777_222,b_777_122,b_777_022,b_777_122,b_777_222,b_777_223,           &      
                                                 b_777_123,b_777_122,b_777_112,b_777_012,b_777_112,b_777_122,b_777_123,           &      
                                                 b_777_023,b_777_022,b_777_012,b_777_002,b_777_012,b_777_022,b_777_023,           &                      
                                                 b_777_123,b_777_122,b_777_112,b_777_012,b_777_112,b_777_122,b_777_123,           &      
                                                 b_777_223,b_777_222,b_777_122,b_777_022,b_777_122,b_777_222,b_777_223,           &      
                                                 b_777_233,b_777_223,b_777_123,b_777_023,b_777_123,b_777_223,b_777_233,           &      
                                                 b_777_333,b_777_233,b_777_133,b_777_033,b_777_133,b_777_233,b_777_333,           &
                                                 b_777_233,b_777_223,b_777_123,b_777_023,b_777_123,b_777_223,b_777_233,           &
                                                 b_777_133,b_777_123,b_777_113,b_777_013,b_777_113,b_777_123,b_777_133,           &
                                                 b_777_033,b_777_023,b_777_013,b_777_003,b_777_013,b_777_023,b_777_033,           &
                                                 b_777_133,b_777_123,b_777_113,b_777_013,b_777_113,b_777_123,b_777_133,           &
                                                 b_777_233,b_777_223,b_777_123,b_777_023,b_777_123,b_777_223,b_777_233,           &
                                                 b_777_333,b_777_233,b_777_133,b_777_033,b_777_133,b_777_233,b_777_333            &
                                                                                                                        /) , (/7,7,7/) )  


     !--- values required for a biharmonic kernel 77  
     !    good noise control
          real(kind=real64),private,parameter     ::      b_77_000 =  1688d0/9702
          real(kind=real64),private,parameter     ::      b_77_001 =  1014d0/9702
          real(kind=real64),private,parameter     ::      b_77_011 =   351d0/9702
          real(kind=real64),private,parameter     ::      b_77_002 =  -126d0/9702
          real(kind=real64),private,parameter     ::      b_77_012 =  -756d0/9702
          real(kind=real64),private,parameter     ::      b_77_022 = -1764d0/9702
          real(kind=real64),private,parameter     ::      b_77_003 =   914d0/9702
          real(kind=real64),private,parameter     ::      b_77_013 =   339d0/9702
          real(kind=real64),private,parameter     ::      b_77_023 =  -504d0/9702
          real(kind=real64),private,parameter     ::      b_77_033 =  1031d0/9702
          


         real(kind=real64),dimension(-3:3,-3:3),public,parameter                ::      SQ_KERNEL_77_2D = reshape( (/   &
                                                 b_77_033,b_77_023,b_77_013,b_77_003,b_77_013,b_77_023,b_77_033,        &     
                                                 b_77_023,b_77_022,b_77_012,b_77_002,b_77_012,b_77_022,b_77_023,        &     
                                                 b_77_013,b_77_012,b_77_011,b_77_001,b_77_011,b_77_012,b_77_013,        &     
                                                 b_77_003,b_77_002,b_77_001,b_77_000,b_77_001,b_77_002,b_77_003,        &     
                                                 b_77_013,b_77_012,b_77_011,b_77_001,b_77_011,b_77_012,b_77_013,        &     
                                                 b_77_023,b_77_022,b_77_012,b_77_002,b_77_012,b_77_022,b_77_023,        &     
                                                 b_77_033,b_77_023,b_77_013,b_77_003,b_77_013,b_77_023,b_77_033         /) , (/7,7/) )  


     !--- values required for a biharmonic kernel 55
     !    warning- poor noise control
          real(kind=real64),private,parameter     ::      b_55_000 =  628d0/245
          real(kind=real64),private,parameter     ::      b_55_001 =  118d0/245
          real(kind=real64),private,parameter     ::      b_55_011 = -382d0/245
          real(kind=real64),private,parameter     ::      b_55_002 =  303d0/245
          real(kind=real64),private,parameter     ::      b_55_012 = -167d0/245
          real(kind=real64),private,parameter     ::      b_55_022 =  138d0/245


        real(kind=real64),dimension(-2:2,-2:2),private,parameter                ::      SQ_KERNEL_55_2D = reshape( (/   &        
                              b_55_022,b_55_012,b_55_002,b_55_012,b_55_022,               &
                              b_55_012,b_55_011,b_55_001,b_55_011,b_55_012,               &
                              b_55_002,b_55_001,b_55_000,b_55_001,b_55_002,               &
                              b_55_012,b_55_011,b_55_001,b_55_011,b_55_012,               &
                              b_55_022,b_55_012,b_55_002,b_55_012,b_55_022     /),(/5,5/) )



        real(kind=real64),dimension(-4:4,-4:4),private,parameter                ::      SQ_KERNEL_99_2D = reshape( (/   &        
   0.013061224489796d0,   0.006530612244898d0,   0.000816326530612d0,   0.006530612244898d0,   0.027755102040816d0,   0.006530612244898d0,   0.000816326530612d0,   0.006530612244898d0,   0.013061224489796d0,   &
   0.006530612244898d0,  -0.011428571428571d0,  -0.022857142857143d0,  -0.016326530612245d0,   0.006530612244898d0,  -0.016326530612245d0,  -0.022857142857143d0,  -0.011428571428571d0,   0.006530612244898d0,   &
   0.000816326530612d0,  -0.022857142857143d0,  -0.032653061224490d0,  -0.019591836734694d0,   0.005714285714286d0,  -0.019591836734694d0,  -0.032653061224490d0,  -0.022857142857143d0,   0.000816326530612d0,   &
   0.006530612244898d0,  -0.016326530612245d0,  -0.019591836734694d0,   0.006530612244898d0,   0.045714285714286d0,   0.006530612244898d0,  -0.019591836734694d0,  -0.016326530612245d0,   0.006530612244898d0,   &
   0.027755102040816d0,   0.006530612244898d0,   0.005714285714286d0,   0.045714285714286d0,   0.114285714285714d0,   0.045714285714286d0,   0.005714285714286d0,   0.006530612244898d0,   0.027755102040816d0,   &
   0.006530612244898d0,  -0.016326530612245d0,  -0.019591836734694d0,   0.006530612244898d0,   0.045714285714286d0,   0.006530612244898d0,  -0.019591836734694d0,  -0.016326530612245d0,   0.006530612244898d0,   &
   0.000816326530612d0,  -0.022857142857143d0,  -0.032653061224490d0,  -0.019591836734694d0,   0.005714285714286d0,  -0.019591836734694d0,  -0.032653061224490d0,  -0.022857142857143d0,   0.000816326530612d0,   &
   0.006530612244898d0,  -0.011428571428571d0,  -0.022857142857143d0,  -0.016326530612245d0,   0.006530612244898d0,  -0.016326530612245d0,  -0.022857142857143d0,  -0.011428571428571d0,   0.006530612244898d0,   &
   0.013061224489796d0,   0.006530612244898d0,   0.000816326530612d0,   0.006530612244898d0,   0.027755102040816d0,   0.006530612244898d0,   0.000816326530612d0,   0.006530612244898d0,   0.013061224489796d0    &
                                                                                /),(/9,9/) )                                                                                



        real(kind=real64),dimension(-4:4,-4:4,-4:4),private,parameter            ::      SQ_KERNEL_999_3D = reshape( (/   &        
   0.000000000000000d0,   0.000000000000000d0,   0.000759816370909d0,   0.001074026148353d0,   0.001899175816800d0,   0.001074026148353d0,   0.000759816370909d0,   0.000000000000000d0,   0.000000000000000d0,        &
   0.000000000000000d0,   0.001519632741819d0,   0.001256839109775d0,   0.001386045262960d0,   0.000179892077640d0,   0.001386045262960d0,   0.001256839109775d0,   0.001519632741819d0,   0.000000000000000d0,        &
   0.000759816370909d0,   0.001256839109775d0,   0.000244495154232d0,   0.000331090767537d0,   0.000538086928508d0,   0.000331090767537d0,   0.000244495154232d0,   0.001256839109775d0,   0.000759816370909d0,        &
   0.001074026148353d0,   0.001386045262960d0,   0.000331090767537d0,   0.001212682219657d0,   0.002387221133453d0,   0.001212682219657d0,   0.000331090767537d0,   0.001386045262960d0,   0.001074026148353d0,        &
   0.001899175816800d0,   0.000179892077640d0,   0.000538086928508d0,   0.002387221133453d0,   0.007896695213294d0,   0.002387221133453d0,   0.000538086928508d0,   0.000179892077640d0,   0.001899175816800d0,        &
   0.001074026148353d0,   0.001386045262960d0,   0.000331090767537d0,   0.001212682219657d0,   0.002387221133453d0,   0.001212682219657d0,   0.000331090767537d0,   0.001386045262960d0,   0.001074026148353d0,        &
   0.000759816370909d0,   0.001256839109775d0,   0.000244495154232d0,   0.000331090767537d0,   0.000538086928508d0,   0.000331090767537d0,   0.000244495154232d0,   0.001256839109775d0,   0.000759816370909d0,        &
   0.000000000000000d0,   0.001519632741819d0,   0.001256839109775d0,   0.001386045262960d0,   0.000179892077640d0,   0.001386045262960d0,   0.001256839109775d0,   0.001519632741819d0,   0.000000000000000d0,        &
   0.000000000000000d0,   0.000000000000000d0,   0.000759816370909d0,   0.001074026148353d0,   0.001899175816800d0,   0.001074026148353d0,   0.000759816370909d0,   0.000000000000000d0,   0.000000000000000d0,        &
   0.000000000000000d0,   0.001519632741819d0,   0.001256839109775d0,   0.001386045262960d0,   0.000179892077640d0,   0.001386045262960d0,   0.001256839109775d0,   0.001519632741819d0,   0.000000000000000d0,        &
   0.001519632741819d0,   0.000548438884265d0,  -0.001528395393165d0,  -0.001930016913206d0,  -0.000309871405944d0,  -0.001930016913206d0,  -0.001528395393165d0,   0.000548438884265d0,   0.001519632741819d0,        &
   0.001256839109775d0,  -0.001528395393165d0,  -0.004481151171844d0,  -0.004736814411125d0,  -0.002072796585143d0,  -0.004736814411125d0,  -0.004481151171844d0,  -0.001528395393165d0,   0.001256839109775d0,        &
   0.001386045262960d0,  -0.001930016913206d0,  -0.004736814411125d0,  -0.003563306397487d0,   0.000536583532444d0,  -0.003563306397487d0,  -0.004736814411125d0,  -0.001930016913206d0,   0.001386045262960d0,        &
   0.000179892077640d0,  -0.000309871405944d0,  -0.002072796585143d0,   0.000536583532444d0,   0.002224167091734d0,   0.000536583532444d0,  -0.002072796585143d0,  -0.000309871405944d0,   0.000179892077640d0,        &
   0.001386045262960d0,  -0.001930016913206d0,  -0.004736814411125d0,  -0.003563306397487d0,   0.000536583532444d0,  -0.003563306397487d0,  -0.004736814411125d0,  -0.001930016913206d0,   0.001386045262960d0,        &
   0.001256839109775d0,  -0.001528395393165d0,  -0.004481151171844d0,  -0.004736814411125d0,  -0.002072796585143d0,  -0.004736814411125d0,  -0.004481151171844d0,  -0.001528395393165d0,   0.001256839109775d0,        &
   0.001519632741819d0,   0.000548438884265d0,  -0.001528395393165d0,  -0.001930016913206d0,  -0.000309871405944d0,  -0.001930016913206d0,  -0.001528395393165d0,   0.000548438884265d0,   0.001519632741819d0,        &
   0.000000000000000d0,   0.001519632741819d0,   0.001256839109775d0,   0.001386045262960d0,   0.000179892077640d0,   0.001386045262960d0,   0.001256839109775d0,   0.001519632741819d0,   0.000000000000000d0,        &
   0.000759816370909d0,   0.001256839109775d0,   0.000244495154232d0,   0.000331090767537d0,   0.000538086928508d0,   0.000331090767537d0,   0.000244495154232d0,   0.001256839109775d0,   0.000759816370909d0,        &
   0.001256839109775d0,  -0.001528395393165d0,  -0.004481151171844d0,  -0.004736814411125d0,  -0.002072796585143d0,  -0.004736814411125d0,  -0.004481151171844d0,  -0.001528395393165d0,   0.001256839109775d0,        &
   0.000244495154232d0,  -0.004481151171844d0,  -0.007640473569760d0,  -0.006896550242960d0,  -0.003284877446374d0,  -0.006896550242960d0,  -0.007640473569760d0,  -0.004481151171844d0,   0.000244495154232d0,        &
   0.000331090767537d0,  -0.004736814411125d0,  -0.006896550242960d0,  -0.003554887379527d0,   0.002103379956509d0,  -0.003554887379527d0,  -0.006896550242960d0,  -0.004736814411125d0,   0.000331090767537d0,        &
   0.000538086928508d0,  -0.002072796585143d0,  -0.003284877446374d0,   0.002103379956509d0,   0.007415823243577d0,   0.002103379956509d0,  -0.003284877446374d0,  -0.002072796585143d0,   0.000538086928508d0,        &
   0.000331090767537d0,  -0.004736814411125d0,  -0.006896550242960d0,  -0.003554887379527d0,   0.002103379956509d0,  -0.003554887379527d0,  -0.006896550242960d0,  -0.004736814411125d0,   0.000331090767537d0,        &
   0.000244495154232d0,  -0.004481151171844d0,  -0.007640473569760d0,  -0.006896550242960d0,  -0.003284877446374d0,  -0.006896550242960d0,  -0.007640473569760d0,  -0.004481151171844d0,   0.000244495154232d0,        &
   0.001256839109775d0,  -0.001528395393165d0,  -0.004481151171844d0,  -0.004736814411125d0,  -0.002072796585143d0,  -0.004736814411125d0,  -0.004481151171844d0,  -0.001528395393165d0,   0.001256839109775d0,        &
   0.000759816370909d0,   0.001256839109775d0,   0.000244495154232d0,   0.000331090767537d0,   0.000538086928508d0,   0.000331090767537d0,   0.000244495154232d0,   0.001256839109775d0,   0.000759816370909d0,        &
   0.001074026148353d0,   0.001386045262960d0,   0.000331090767537d0,   0.001212682219657d0,   0.002387221133453d0,   0.001212682219657d0,   0.000331090767537d0,   0.001386045262960d0,   0.001074026148353d0,        &
   0.001386045262960d0,  -0.001930016913206d0,  -0.004736814411125d0,  -0.003563306397487d0,   0.000536583532444d0,  -0.003563306397487d0,  -0.004736814411125d0,  -0.001930016913206d0,   0.001386045262960d0,        &
   0.000331090767537d0,  -0.004736814411125d0,  -0.006896550242960d0,  -0.003554887379527d0,   0.002103379956509d0,  -0.003554887379527d0,  -0.006896550242960d0,  -0.004736814411125d0,   0.000331090767537d0,        &
   0.001212682219657d0,  -0.003563306397487d0,  -0.003554887379527d0,   0.003998689897715d0,   0.013151880587398d0,   0.003998689897715d0,  -0.003554887379527d0,  -0.003563306397487d0,   0.001212682219657d0,        &
   0.002387221133453d0,   0.000536583532444d0,   0.002103379956509d0,   0.013151880587398d0,   0.023580122956321d0,   0.013151880587398d0,   0.002103379956509d0,   0.000536583532444d0,   0.002387221133453d0,        &
   0.001212682219657d0,  -0.003563306397487d0,  -0.003554887379527d0,   0.003998689897715d0,   0.013151880587398d0,   0.003998689897715d0,  -0.003554887379527d0,  -0.003563306397487d0,   0.001212682219657d0,        &
   0.000331090767537d0,  -0.004736814411125d0,  -0.006896550242960d0,  -0.003554887379527d0,   0.002103379956509d0,  -0.003554887379527d0,  -0.006896550242960d0,  -0.004736814411125d0,   0.000331090767537d0,        &
   0.001386045262960d0,  -0.001930016913206d0,  -0.004736814411125d0,  -0.003563306397487d0,   0.000536583532444d0,  -0.003563306397487d0,  -0.004736814411125d0,  -0.001930016913206d0,   0.001386045262960d0,        &
   0.001074026148353d0,   0.001386045262960d0,   0.000331090767537d0,   0.001212682219657d0,   0.002387221133453d0,   0.001212682219657d0,   0.000331090767537d0,   0.001386045262960d0,   0.001074026148353d0,        &
   0.001899175816800d0,   0.000179892077640d0,   0.000538086928508d0,   0.002387221133453d0,   0.007896695213294d0,   0.002387221133453d0,   0.000538086928508d0,   0.000179892077640d0,   0.001899175816800d0,        &
   0.000179892077640d0,  -0.000309871405944d0,  -0.002072796585143d0,   0.000536583532444d0,   0.002224167091734d0,   0.000536583532444d0,  -0.002072796585143d0,  -0.000309871405944d0,   0.000179892077640d0,        &
   0.000538086928508d0,  -0.002072796585143d0,  -0.003284877446374d0,   0.002103379956509d0,   0.007415823243577d0,   0.002103379956509d0,  -0.003284877446374d0,  -0.002072796585143d0,   0.000538086928508d0,        &
   0.002387221133453d0,   0.000536583532444d0,   0.002103379956509d0,   0.013151880587398d0,   0.023580122956321d0,   0.013151880587398d0,   0.002103379956509d0,   0.000536583532444d0,   0.002387221133453d0,        &
   0.007896695213294d0,   0.002224167091734d0,   0.007415823243577d0,   0.023580122956321d0,   0.048497409326425d0,   0.023580122956321d0,   0.007415823243577d0,   0.002224167091734d0,   0.007896695213294d0,        &
   0.002387221133453d0,   0.000536583532444d0,   0.002103379956509d0,   0.013151880587398d0,   0.023580122956321d0,   0.013151880587398d0,   0.002103379956509d0,   0.000536583532444d0,   0.002387221133453d0,        &
   0.000538086928508d0,  -0.002072796585143d0,  -0.003284877446374d0,   0.002103379956509d0,   0.007415823243577d0,   0.002103379956509d0,  -0.003284877446374d0,  -0.002072796585143d0,   0.000538086928508d0,        &
   0.000179892077640d0,  -0.000309871405944d0,  -0.002072796585143d0,   0.000536583532444d0,   0.002224167091734d0,   0.000536583532444d0,  -0.002072796585143d0,  -0.000309871405944d0,   0.000179892077640d0,        &
   0.001899175816800d0,   0.000179892077640d0,   0.000538086928508d0,   0.002387221133453d0,   0.007896695213294d0,   0.002387221133453d0,   0.000538086928508d0,   0.000179892077640d0,   0.001899175816800d0,        &
   0.001074026148353d0,   0.001386045262960d0,   0.000331090767537d0,   0.001212682219657d0,   0.002387221133453d0,   0.001212682219657d0,   0.000331090767537d0,   0.001386045262960d0,   0.001074026148353d0,        &
   0.001386045262960d0,  -0.001930016913206d0,  -0.004736814411125d0,  -0.003563306397487d0,   0.000536583532444d0,  -0.003563306397487d0,  -0.004736814411125d0,  -0.001930016913206d0,   0.001386045262960d0,        &
   0.000331090767537d0,  -0.004736814411125d0,  -0.006896550242960d0,  -0.003554887379527d0,   0.002103379956509d0,  -0.003554887379527d0,  -0.006896550242960d0,  -0.004736814411125d0,   0.000331090767537d0,        &
   0.001212682219657d0,  -0.003563306397487d0,  -0.003554887379527d0,   0.003998689897715d0,   0.013151880587398d0,   0.003998689897715d0,  -0.003554887379527d0,  -0.003563306397487d0,   0.001212682219657d0,        &
   0.002387221133453d0,   0.000536583532444d0,   0.002103379956509d0,   0.013151880587398d0,   0.023580122956321d0,   0.013151880587398d0,   0.002103379956509d0,   0.000536583532444d0,   0.002387221133453d0,        &
   0.001212682219657d0,  -0.003563306397487d0,  -0.003554887379527d0,   0.003998689897715d0,   0.013151880587398d0,   0.003998689897715d0,  -0.003554887379527d0,  -0.003563306397487d0,   0.001212682219657d0,        &
   0.000331090767537d0,  -0.004736814411125d0,  -0.006896550242960d0,  -0.003554887379527d0,   0.002103379956509d0,  -0.003554887379527d0,  -0.006896550242960d0,  -0.004736814411125d0,   0.000331090767537d0,        &
   0.001386045262960d0,  -0.001930016913206d0,  -0.004736814411125d0,  -0.003563306397487d0,   0.000536583532444d0,  -0.003563306397487d0,  -0.004736814411125d0,  -0.001930016913206d0,   0.001386045262960d0,        &
   0.001074026148353d0,   0.001386045262960d0,   0.000331090767537d0,   0.001212682219657d0,   0.002387221133453d0,   0.001212682219657d0,   0.000331090767537d0,   0.001386045262960d0,   0.001074026148353d0,        &
   0.000759816370909d0,   0.001256839109775d0,   0.000244495154232d0,   0.000331090767537d0,   0.000538086928508d0,   0.000331090767537d0,   0.000244495154232d0,   0.001256839109775d0,   0.000759816370909d0,        &
   0.001256839109775d0,  -0.001528395393165d0,  -0.004481151171844d0,  -0.004736814411125d0,  -0.002072796585143d0,  -0.004736814411125d0,  -0.004481151171844d0,  -0.001528395393165d0,   0.001256839109775d0,        &
   0.000244495154232d0,  -0.004481151171844d0,  -0.007640473569760d0,  -0.006896550242960d0,  -0.003284877446374d0,  -0.006896550242960d0,  -0.007640473569760d0,  -0.004481151171844d0,   0.000244495154232d0,        &
   0.000331090767537d0,  -0.004736814411125d0,  -0.006896550242960d0,  -0.003554887379527d0,   0.002103379956509d0,  -0.003554887379527d0,  -0.006896550242960d0,  -0.004736814411125d0,   0.000331090767537d0,        &
   0.000538086928508d0,  -0.002072796585143d0,  -0.003284877446374d0,   0.002103379956509d0,   0.007415823243577d0,   0.002103379956509d0,  -0.003284877446374d0,  -0.002072796585143d0,   0.000538086928508d0,        &
   0.000331090767537d0,  -0.004736814411125d0,  -0.006896550242960d0,  -0.003554887379527d0,   0.002103379956509d0,  -0.003554887379527d0,  -0.006896550242960d0,  -0.004736814411125d0,   0.000331090767537d0,        &
   0.000244495154232d0,  -0.004481151171844d0,  -0.007640473569760d0,  -0.006896550242960d0,  -0.003284877446374d0,  -0.006896550242960d0,  -0.007640473569760d0,  -0.004481151171844d0,   0.000244495154232d0,        &
   0.001256839109775d0,  -0.001528395393165d0,  -0.004481151171844d0,  -0.004736814411125d0,  -0.002072796585143d0,  -0.004736814411125d0,  -0.004481151171844d0,  -0.001528395393165d0,   0.001256839109775d0,        &
   0.000759816370909d0,   0.001256839109775d0,   0.000244495154232d0,   0.000331090767537d0,   0.000538086928508d0,   0.000331090767537d0,   0.000244495154232d0,   0.001256839109775d0,   0.000759816370909d0,        &
   0.000000000000000d0,   0.001519632741819d0,   0.001256839109775d0,   0.001386045262960d0,   0.000179892077640d0,   0.001386045262960d0,   0.001256839109775d0,   0.001519632741819d0,   0.000000000000000d0,        &
   0.001519632741819d0,   0.000548438884265d0,  -0.001528395393165d0,  -0.001930016913206d0,  -0.000309871405944d0,  -0.001930016913206d0,  -0.001528395393165d0,   0.000548438884265d0,   0.001519632741819d0,        &
   0.001256839109775d0,  -0.001528395393165d0,  -0.004481151171844d0,  -0.004736814411125d0,  -0.002072796585143d0,  -0.004736814411125d0,  -0.004481151171844d0,  -0.001528395393165d0,   0.001256839109775d0,        &
   0.001386045262960d0,  -0.001930016913206d0,  -0.004736814411125d0,  -0.003563306397487d0,   0.000536583532444d0,  -0.003563306397487d0,  -0.004736814411125d0,  -0.001930016913206d0,   0.001386045262960d0,        &
   0.000179892077640d0,  -0.000309871405944d0,  -0.002072796585143d0,   0.000536583532444d0,   0.002224167091734d0,   0.000536583532444d0,  -0.002072796585143d0,  -0.000309871405944d0,   0.000179892077640d0,        &
   0.001386045262960d0,  -0.001930016913206d0,  -0.004736814411125d0,  -0.003563306397487d0,   0.000536583532444d0,  -0.003563306397487d0,  -0.004736814411125d0,  -0.001930016913206d0,   0.001386045262960d0,        &
   0.001256839109775d0,  -0.001528395393165d0,  -0.004481151171844d0,  -0.004736814411125d0,  -0.002072796585143d0,  -0.004736814411125d0,  -0.004481151171844d0,  -0.001528395393165d0,   0.001256839109775d0,        &
   0.001519632741819d0,   0.000548438884265d0,  -0.001528395393165d0,  -0.001930016913206d0,  -0.000309871405944d0,  -0.001930016913206d0,  -0.001528395393165d0,   0.000548438884265d0,   0.001519632741819d0,        &
   0.000000000000000d0,   0.001519632741819d0,   0.001256839109775d0,   0.001386045262960d0,   0.000179892077640d0,   0.001386045262960d0,   0.001256839109775d0,   0.001519632741819d0,   0.000000000000000d0,        &
   0.000000000000000d0,   0.000000000000000d0,   0.000759816370909d0,   0.001074026148353d0,   0.001899175816800d0,   0.001074026148353d0,   0.000759816370909d0,   0.000000000000000d0,   0.000000000000000d0,        &
   0.000000000000000d0,   0.001519632741819d0,   0.001256839109775d0,   0.001386045262960d0,   0.000179892077640d0,   0.001386045262960d0,   0.001256839109775d0,   0.001519632741819d0,   0.000000000000000d0,        &
   0.000759816370909d0,   0.001256839109775d0,   0.000244495154232d0,   0.000331090767537d0,   0.000538086928508d0,   0.000331090767537d0,   0.000244495154232d0,   0.001256839109775d0,   0.000759816370909d0,        &
   0.001074026148353d0,   0.001386045262960d0,   0.000331090767537d0,   0.001212682219657d0,   0.002387221133453d0,   0.001212682219657d0,   0.000331090767537d0,   0.001386045262960d0,   0.001074026148353d0,        &
   0.001899175816800d0,   0.000179892077640d0,   0.000538086928508d0,   0.002387221133453d0,   0.007896695213294d0,   0.002387221133453d0,   0.000538086928508d0,   0.000179892077640d0,   0.001899175816800d0,        &
   0.001074026148353d0,   0.001386045262960d0,   0.000331090767537d0,   0.001212682219657d0,   0.002387221133453d0,   0.001212682219657d0,   0.000331090767537d0,   0.001386045262960d0,   0.001074026148353d0,        &
   0.000759816370909d0,   0.001256839109775d0,   0.000244495154232d0,   0.000331090767537d0,   0.000538086928508d0,   0.000331090767537d0,   0.000244495154232d0,   0.001256839109775d0,   0.000759816370909d0,        &
   0.000000000000000d0,   0.001519632741819d0,   0.001256839109775d0,   0.001386045262960d0,   0.000179892077640d0,   0.001386045262960d0,   0.001256839109775d0,   0.001519632741819d0,   0.000000000000000d0,        &
   0.000000000000000d0,   0.000000000000000d0,   0.000759816370909d0,   0.001074026148353d0,   0.001899175816800d0,   0.001074026148353d0,   0.000759816370909d0,   0.000000000000000d0,   0.000000000000000d0         &
                                                                                /),(/9,9,9/) )                                                                                


!----------------------- JUNK KERNELS TO BE DELETED
!     3x3x3 testing
     !     error    0.0000000000000000        9.0070491243193062E-003   5.0061681478386078E-004
     !     error    8.9185788867613547E-005   1.5828088062424387E-003 -0.82418201953960035
     !     error    1.7837157773522709E-004  -2.2326173480203753E-002  -3.4799854014059060
     !     error    2.6755736660284061E-004   1.4031576896916805E-005 -0.99844137617698292
     !     error    3.5674315547045419E-004   1.2889236723164636E-002  0.43173298089151424
     !     error    4.4592894433806765E-004  -9.8938532657566455E-002  -11.990065845580226
     !     error    5.3511473320568123E-004   2.7465058763624643E-002   2.050811409440093802
     !     error    6.2430052207329480E-004  -3.2593195190576403E-003  -1.3620443437361707
     !     error    7.1348631094090838E-004   3.0570680184837863E-002   2.3957830094203314
     !     error    8.0267209980852184E-004  -4.4654887431434455E-002  -5.9602529976565801
     !     error    8.9185788867613531E-004  -2.3313974083405405E-002  -3.5897100292115698

     !   real(kind=real64),dimension(-1:1,-1:1,-1:1),public,parameter           ::      NOISY_KERNEL_333_3D = reshape(   (/        &
     !                                     35,  2, 35,                                                                             &   
     !                                      2,-37,  2,                                                                             &
     !                                     35,  2, 35,                                                                             &
     !                                      2,-37,  2,                                                                             &
     !                                    -37,-82,-37,                                                                             &
     !                                      2,-37,  2,                                                                             & 
     !                                     35,  2, 35,                                                                             &   
     !                                      2,-37,  2,                                                                             &
     !                                     35,  2, 35       /),(/3,3,3/) ) / 111.0d0


     !    real(kind=real64),dimension(-1:1,-1:1,-1:1),private,parameter           ::      KERNEL_333_3D = reshape(   (/                                      &
     !                        2,3,2,  3,6,3,  2,3,2       ,       3,6,3,  6,-88,6,    3,6,3       ,       2,3,2,  3,6,3,  2,3,2    /),(/3,3,3/) ) / 26.0d0

        ! real(kind=real64),dimension(-1:1,-1:1,-1:1),private,parameter           ::      KERNEL_333_3D = reshape(   (/                                      &
        !                             -1,  2, -1,                         &
        !                              2, -1,  2,                         &
        !                             -1,  2, -1,                         &
        !                              2, -1,  2,                         &
        !                             -1,-10, -1,                         &
        !                              2, -1,  2,                         &
        !                             -1,  2, -1,                         &
        !                              2, -1,  2,                         &
        !                             -1,  2, -1          /) , (/3,3,3/) ) / 3.0d0
                  
     !   real(kind=real64),private,parameter     ::      noisy_k00 = -0.295771d0
     !   real(kind=real64),private,parameter     ::      noisy_k01 = -0.166276d0
     !   real(kind=real64),private,parameter     ::      noisy_k11 = -0.0745175d0
     !   real(kind=real64),private,parameter     ::      noisy_k02 = 0.133127d0
     !   real(kind=real64),private,parameter     ::      noisy_k12 =  0.111679d0
     !   real(kind=real64),private,parameter     ::      noisy_k22 =  -0.041748d0

     
!        real(kind=real64),private,parameter     ::      k000 = -5850.0d0/14700 
!        real(kind=real64),private,parameter     ::      k001 = -3116.0d0/14700 
!        real(kind=real64),private,parameter     ::      k002 = -2509.0d0/14700 
!        real(kind=real64),private,parameter     ::      k011 =  -502.0d0/14700 
!        real(kind=real64),private,parameter     ::      k012 =  -255.0d0/14700 
!        real(kind=real64),private,parameter     ::      k022 = -1088.0d0/14700 
!        real(kind=real64),private,parameter     ::      k111 =  1992.0d0/14700 
!        real(kind=real64),private,parameter     ::      k112 =  1879.0d0/14700 
!        real(kind=real64),private,parameter     ::      k122 =   686.0d0/14700 
!        real(kind=real64),private,parameter     ::      k222 = -1587.0d0/14700     
!         
!         real(kind=real64),private,parameter     ::      k000 = -39943.0d0/78450 
!         real(kind=real64),private,parameter     ::      k001 = -10517.0d0/39225 
!         real(kind=real64),private,parameter     ::      k002 = -6753.0d0/52300 
!         real(kind=real64),private,parameter     ::      k011 =  -769.0d0/15690 
!         real(kind=real64),private,parameter     ::      k012 =  3799.0d0/156900
!         real(kind=real64),private,parameter     ::      k022 = -1306.0d0/13075
!         real(kind=real64),private,parameter     ::      k111 =  5812d0/39225
!         real(kind=real64),private,parameter     ::      k112 =  8139d0/52300
!         real(kind=real64),private,parameter     ::      k122 =  -2687d0/78450
!         real(kind=real64),private,parameter     ::      k222 = 0.0d0     
!
!          real(kind=real64),private,parameter     ::      k000 = -54d0/81 
!          real(kind=real64),private,parameter     ::      k001 = -15d0/81 
!          real(kind=real64),private,parameter     ::      k002 = -12d0/81
!          real(kind=real64),private,parameter     ::      k011 =   2d0/81
!          real(kind=real64),private,parameter     ::      k012 =   1d0/81
!          real(kind=real64),private,parameter     ::      k022 =  -4d0/81
!          real(kind=real64),private,parameter     ::      k111 =   8d0/81
!          real(kind=real64),private,parameter     ::      k112 =   6d0/81
!          real(kind=real64),private,parameter     ::      k122 =   1d0/81
!          real(kind=real64),private,parameter     ::      k222 =  -2d0/81


     !--- values required for a kernel expected to work on noisy data
     !     real(kind=real64),private,parameter     ::      noisy_k000 = -94179d0/603950
     !     real(kind=real64),private,parameter     ::      noisy_k001 = -18733d0/181185
     !     real(kind=real64),private,parameter     ::      noisy_k011 = -106223d0/1811850
     !     real(kind=real64),private,parameter     ::      noisy_k111 = -6536d0/301975
     !     real(kind=real64),private,parameter     ::      noisy_k002 = -2593d0/3623700
     !     real(kind=real64),private,parameter     ::      noisy_k012 = 25007d0/1207900
     !     real(kind=real64),private,parameter     ::      noisy_k112 = 24887d0/724740
     !     real(kind=real64),private,parameter     ::      noisy_k022 = 27172d0/905925
     !     real(kind=real64),private,parameter     ::      noisy_k122 = 143d0/7050
     !     real(kind=real64),private,parameter     ::      noisy_k222 = -77077d0/1207900


     
!        real(kind=real64),dimension(-2:2,-2:2,-2:2),private,parameter            ::      SQ_KERNEL_333_3D = reshape( (/   &        
!   0.005917159763314d0,  0.017751479289941d0,   0.025147928994083d0,   0.017751479289941d0,   0.005917159763314d0,          &
!   0.017751479289941d0,  0.062130177514793d0,   0.088757396449704d0,   0.062130177514793d0,   0.017751479289941d0,          &
!   0.025147928994083d0,  0.088757396449704d0,   0.130177514792899d0,   0.088757396449704d0,   0.025147928994083d0,          &
!   0.017751479289941d0,  0.062130177514793d0,   0.088757396449704d0,   0.062130177514793d0,   0.017751479289941d0,          &
!   0.005917159763314d0,  0.017751479289941d0,   0.025147928994083d0,   0.017751479289941d0,   0.005917159763314d0,          &
!   0.017751479289941d0,  0.062130177514793d0,   0.088757396449704d0,   0.062130177514793d0,   0.017751479289941d0,          &
!   0.062130177514793d0, -0.360946745562130d0,  -0.550295857988166d0,  -0.360946745562130d0,   0.062130177514793d0,          &
!   0.088757396449704d0, -0.550295857988166d0,  -1.278106508875740d0,  -0.550295857988166d0,   0.088757396449704d0,          &
!   0.062130177514793d0, -0.360946745562130d0,  -0.550295857988166d0,  -0.360946745562130d0,   0.062130177514793d0,          &
!   0.017751479289941d0,  0.062130177514793d0,   0.088757396449704d0,   0.062130177514793d0,   0.017751479289941d0,          &
!   0.025147928994083d0,  0.088757396449704d0,   0.130177514792899d0,   0.088757396449704d0,   0.025147928994083d0,          &
!   0.088757396449704d0, -0.550295857988166d0,  -1.278106508875740d0,  -0.550295857988166d0,   0.088757396449704d0,          &
!   0.130177514792899d0, -1.278106508875740d0,  11.982248520710066d0,  -1.278106508875740d0,   0.130177514792899d0,          &
!   0.088757396449704d0, -0.550295857988166d0,  -1.278106508875740d0,  -0.550295857988166d0,   0.088757396449704d0,          &
!   0.025147928994083d0,  0.088757396449704d0,   0.130177514792899d0,   0.088757396449704d0,   0.025147928994083d0,          &
!   0.017751479289941d0,  0.062130177514793d0,   0.088757396449704d0,   0.062130177514793d0,   0.017751479289941d0,          &
!   0.062130177514793d0, -0.360946745562130d0,  -0.550295857988166d0,  -0.360946745562130d0,   0.062130177514793d0,          &
!   0.088757396449704d0, -0.550295857988166d0,  -1.278106508875740d0,  -0.550295857988166d0,   0.088757396449704d0,          &
!   0.062130177514793d0, -0.360946745562130d0,  -0.550295857988166d0,  -0.360946745562130d0,   0.062130177514793d0,          &
!   0.017751479289941d0,  0.062130177514793d0,   0.088757396449704d0,   0.062130177514793d0,   0.017751479289941d0,          &
!   0.005917159763314d0,  0.017751479289941d0,   0.025147928994083d0,   0.017751479289941d0,   0.005917159763314d0,          &
!   0.017751479289941d0,  0.062130177514793d0,   0.088757396449704d0,   0.062130177514793d0,   0.017751479289941d0,          &
!   0.025147928994083d0,  0.088757396449704d0,   0.130177514792899d0,   0.088757396449704d0,   0.025147928994083d0,          &
!   0.017751479289941d0,  0.062130177514793d0,   0.088757396449704d0,   0.062130177514793d0,   0.017751479289941d0,          &
!   0.005917159763314d0,  0.017751479289941d0,   0.025147928994083d0,   0.017751479289941d0,   0.005917159763314d0           &
!                                                                                /),(/5,5,5/) )

