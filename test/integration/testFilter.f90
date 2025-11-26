!------This is a test program to check the main filter routines (salt & pepper filter, max likelyhoodfilter ) in Flat are working correctly.
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    testFilter from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
!*     version history
!*      0.1.0       Jul 2025        Changed old maxlikelihood to new 
program testFilter

   use iso_fortran_env
   use Lib_ColouredTerminal !needed
   use Lib_UtilsForTests
   use Lib_CommandLineArguments
   use Lib_LocalFilter2d
   use Lib_RidlerCalvard
   use Lib_LowPassFilter3d
   use Lib_MaxLikelihoodFilter
   use Lib_Png
   use Lib_FFTW3f
   use Lib_DataDoubler
   use Lib_Filenames
   implicit none
   

   !**----testing vars----**
   logical             ::      ok
   integer             ::      correctcount
   integer             ::      noofTests=2 !
   character(len=32)   ::      libName="testFilter"
   logical             ::      tempCheck
   real(kind=real64),dimension(0:255)          ::      hist

   
   !**----From flat---**
!---    constants
   real(kind=real64),parameter                     ::      PI = 3.141592654d0         
   character(len=8)                    ::      VERSION = "0.1.0"
   
!---    command line arguments (CLA)
   character(len=256)                  ::      filename = "" !see "get(cla, ...)" calls below
   character(len=256)                  ::      outfile = ""
   logical                             ::      zeroIsUnset = .true.
   logical                             ::      sandp = .false.
   logical                             ::      negative = .false.
   logical                             ::      half = .false.
   logical                             ::      flat = .true.
   logical                             ::      retone = .false.
   logical                             ::      dbg = .false.
   
   logical                             ::      quarter = .false.
   logical                             ::      rescale = .false.
   real(kind=real64)                   ::      lambda = LIB_CLA_NODEFAULT_R
   integer                             ::      vRC_in = LIB_CLA_NODEFAULT_I
   
   integer                             ::      nGB = 20            !   number of gaussian blurs used
   integer                             ::      nLoops = 2
   integer                             ::      power = 1
   integer                             ::      hpf_coeff = 0
   type(CommandLineArguments)          ::      cla                
   !character(len=8)                    ::      VERSION = "1.0"
   
   
   
!---    key physical quantities and arrays
   real(kind=real64),dimension(:,:),allocatable    ::      img             !   image input
   integer                                         ::      Nx,Ny           !   size of image in pixels
   real(kind=real64),dimension(:,:),allocatable    ::      img_out         !   image output
   real(kind=real64)                               ::      stddev          !   image standard deviation
   
   
   
!---    dummy 
   integer                             ::      nn                          !Number of zero intensity pixels
   real(kind=real64)                   ::      area                        !iamge area/pixel squared
   real(kind=real64)                   ::      bb,th,ff,bstd,fbar,snr,fot      !   Ridler-Calvard values, back,thresh,fore, back std dev,<f>,signal-to-noise, frac over thresh

         

   !---more testing set up---

   correctcount=0
   !-----------------------------------------------

   !----main prog start----

   cla = CommandLineArguments_ctor(50)
   call setProgramDescription( cla, "flattenBackground.exe "    &
      //"\n    utility program designed to help feature detection in greyscale images \n    by removing background fluctuations and performing a couple of simple cleaning tasks" )
   call setProgramVersion( cla, VERSION )
   
   call get( cla,"f",filename ,LIB_CLA_REQUIRED,"                    filename" )  
   outfile = trim(removeSuffix(filename))//".flat.png"
   call get( cla,"o",outfile ,LIB_CLA_OPTIONAL,"                  output filename" )     
   call get( cla,"zero",zeroIsUnset,LIB_CLA_OPTIONAL,"               interpret pixels with intensity zero as 'unset' rather than just very dark" )     
   call get( cla,"sandp",sandp,LIB_CLA_OPTIONAL,"              apply salt and pepper filter" )     
   call get( cla,"lambda",lambda ,LIB_CLA_OPTIONAL,"             set range for max likelihood filter" )  
   call get( cla,"negative",negative,LIB_CLA_OPTIONAL,"           invert image " )     
   call get( cla,"half",half,LIB_CLA_OPTIONAL,"               compute background using 1/2 size image (default T if Nx>2000)" )     
   call get( cla,"quarter",quarter,LIB_CLA_OPTIONAL,"            compute background using 1/4 size image (default T if Nx>4000)" )     
   call get( cla,"g",nGB,LIB_CLA_OPTIONAL,"                  number of gaussian blurs to use for lengthscale estimate" )     
   call get( cla,"n",nLoops,LIB_CLA_OPTIONAL,"                  number of iterations" )     
   call get( cla,"p",power,LIB_CLA_OPTIONAL,"                  power law for gradient contributions" )     
   call get( cla,"flat",flat,LIB_CLA_OPTIONAL,"               flatten image using scale invariant background subtraction" )     
   call get( cla,"retone",retone,LIB_CLA_OPTIONAL,"             retone image to preserve foreground/background intensity" )     
   rescale = hasArgument(cla,"half") .or. hasArgument(cla,"quarter")
   call get( cla,"hpf",hpf_coeff,LIB_CLA_OPTIONAL,"                high pass filter - fourier coefficients to remove" )     
   call get( cla,"vRC",vRC_in,LIB_CLA_OPTIONAL,"                minimum variable-Ridler-Calvard length scale" )     
   call get( cla,"dbg",dbg,LIB_CLA_OPTIONAL,"                output intermediate images for debugging" )     
   
   if (getSuffix(outfile)/="png") outfile = trim(outfile)//".png"
   
   
   call report(cla)
   if (.not. allRequiredArgumentsSet(cla)) stop
   if (hasHelpArgument(cla)) stop
   call delete(cla)
   

   if (.not. flat) nLoops = 0

   !---    load in input file
   call read_greyscale_png( filename ,img )        
   Nx = size(img,dim=1)
   Ny = size(img,dim=2)
   if (.not. rescale) then
      if (min(Nx,Ny)>8000) then
         half = .true.; quarter = .true.
      else if (min(Nx,Ny)>4000) then
         quarter = .true.
      else if (min(Nx,Ny)>2000) then
         half = .true.
      end if
   end if

   area = Nx*Ny
   allocate(img_out(0:Nx-1,0:Ny-1))

   if (zeroIsUnset) then
      if(sandp) then
         print *,""
         print *,"Salt and pepper filter overides 'zero', not converting pixels with intensity zero to 'unset'"
      else

         print *,""
         print *,"converting pixels with intensity zero to 'unset'"
         where (img==0) 
               img = LOCALFILTER2D_IGNORE
         end where
         nn = count(img==LOCALFILTER2D_IGNORE)
         print *,"count of dead pixels ",nn,"  (",nn*100.0/area,"%)"
         print *,""
      end if
   end if


   if (negative) then
      img = 1 - img
   end if

   if (sandp) then
      print *,""
      print *,"salt and pepper filter"
      call setLocalFilter2d_sigma(0.707d0)
      call denoise( img,3.0d0,img_out(0:Nx-1,0:Ny-1),stddev )  
      img = min(1.0d0,max(0.0d0,img_out(0:Nx-1,0:Ny-1)))
      nn = count(img==LOCALFILTER2D_IGNORE)
      print *,"count of dead pixels ",nn,"  (",nn*100.0/area,"%)"
      print *,"" 
   end if

   !test 1
   tempCheck=.false.

   if((maxval(img_out)<1).and.(minval(img_out)>0)) tempCheck=.true.
   call announceSubTest(libName,"salt and pepper filter",1,noofTests,tempCheck,correctcount)


   if (lambda/=LIB_CLA_NODEFAULT_R) then
      print *,""
      print *,"max likelihood filter"
      if (zeroIsUnset) then
         if (negative) then
               where (img==LOCALFILTER2D_IGNORE) 
                  img = 1                  
               end where
         else
               where (img==LOCALFILTER2D_IGNORE) 
                  img = 0                  
               end where
         end if
      end if
      
         
      if ( (lambda /= LIB_CLA_NODEFAULT_R).and.(lambda>0) ) then
         
         print *,""
         print *,"log likelihood filtering file with lambda = ",lambda
         print *,""

         call findHist(img,hist)

         call findImageIntensityFeatures( img, bb,th,ff,bstd,fbar,snr,fot )

      !---    old maxlikelihood 
      !    th = bstd/lambda
      !    call maxLikelihoodFilter( img,img_out, bstd,th )
      !---    new maxlikelihood 16/07/25
         call maxLikelihoodFilter( img,img_out, lambda )

         img = min(1.0d0,max(0.0d0,img_out(0:Nx-1,0:Ny-1)))
         call findImageIntensityFeatures( img, bb,th,ff,bstd,fbar,snr,fot )

         !test 2
         tempCheck=.false.
         if(snr>4) tempCheck=.true.
         call announceSubTest(libName,"max likelihood filter",2,noofTests,tempCheck,correctcount)
         
      end if   

      print *,"" 
   end if

   !----main prog end------
   !call write_greyscale_png( outfile,img_out,negative=negative )  


   ok=haveAllSubTestsPassed(correctcount,noofTests)

   call announcePassOrFail(ok)


end program testFilter