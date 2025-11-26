program imageUtilities

    !---^^^^^^^^^^^^^^^^ 
    !*      James Heath & Daniel Mason
    !*      (c) UKAEA October 2024
    !---^^^^^^^^^^^^^^^^ 
    !*-----------------------------------------------------------------------------------------------------------------------------------
    ! *    imageUtilities from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
    ! *    Copyright (C) 2024  James Heath

    ! *    This program is free software: you can redistribute it and/or modify
    ! *    it under the terms of the GNU General Public License as published by
    ! *    the Free Software Foundation, either version 3 of the License, or
    ! *    (at your option) any later version.

    ! *    This program is distributed in the hope that it will be useful,
    ! *    but WITHOUT ANY WARRANTY; without even the implied warranty of
    ! *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ! *    GNU General Public License for more details.

    ! *    You should have received a copy of the GNU General Public License
    ! *    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    ! *-----------------------------------------------------------------------------------------------------------------------------------
    !   Program for little utilities that don't warrant their own program but need a home.
    use iso_fortran_env
    
    use Lib_CommandLineArguments
    use NBAX_StringTokenizers
    use Lib_Filenames
    use Lib_Png
    use Lib_RandomSeed
    use Lib_Maxima2d
    use Lib_RidlerCalvard
    use Lib_FFTW3f

    
    implicit none
                
    !---
    !   version numbering
    
    !   0.0.1                   first save
    
            character(len=*),parameter      ::      VERSION = "0.0.1"
    !---        parameters
            type(CommandLineArguments)      ::      cla                     !Object to hold commmand line argumments (CLA)
        
            character(len=256)              ::      image_filename = ""
            character(len=256)              ::      output_filename = ""
       !---        dummy variables  
            real(kind=real64),dimension(:,:),allocatable    ::      img, img_out
            logical                                         ::      ok      !does input file exist? 
            integer                                         ::      Nx,Ny, jj, ii
            real(kind=real64)                               ::          noise_stddev = 0.0
            logical                                         ::      bias=.false.
            real(kind=real64)           ::      zeta
            !real(kind=real64)                               ::      fbar, fstddev ! image intenisty avearge and standard devaition
            !real(kind=real64)                               ::      ff, f2bar !var to store pixel values in

            integer                     ::          seed = 12345
            real(kind=real64)           ::          salt = 0.0d0
            real(kind=real64)           ::          pepper = 0.0d0
            logical                     ::          write_img=.false.
            logical                     ::          get_fft_power_spec=.false.
            real(kind=real64),dimension(:),allocatable         ::      cprime,rpf 
            real(kind=real64),dimension(:,:),allocatable         ::      fftout
            character(len=256)              ::      tempFilename
        !---    read command line arguments
            cla = CommandLineArguments_ctor(30)  
             
            call setProgramDescription( cla, "imageUtilities" )
            call setProgramVersion( cla, VERSION ) 
           
                        
            
        !---    filename options        
            call get( cla,"f",image_filename ,LIB_CLA_REQUIRED,"           image filename" )
            call get( cla,"o",output_filename ,LIB_CLA_OPTIONAL,"      output image filename" )  
            call get( cla,"seed",seed ,LIB_CLA_OPTIONAL,"   random seed ( set to 0 to use system clock )" )
            call get( cla,"stddev",noise_stddev ,LIB_CLA_OPTIONAL," noise standard deviation" )   
            call get( cla,"salt",salt ,LIB_CLA_OPTIONAL,"   fraction of pixels set to white" )  
            call get( cla,"pepper",pepper ,LIB_CLA_OPTIONAL," fraction of pixels set to black" ) 
            call get( cla,"bias",bias ,LIB_CLA_OPTIONAL," Estimate if background standard deviation calcualtion is biased" )
            call get( cla,"fftps",get_fft_power_spec,LIB_CLA_OPTIONAL,"Calculate FFT power spectrum of input image.")
    
            call report(cla)
            if (.not. allRequiredArgumentsSet(cla)) stop
            if (hasHelpArgument(cla)) stop
            call delete(cla)
    
            if (len_trim(image_filename)/=0) then    
                inquire(file=trim(image_filename),exist=ok)
                if (.not. ok) then
                    print *,"imageUtilities  error - couldn't find image file """//trim(image_filename)//""""
                    stop
                end if
                print *,"imageUtilities info - reading image file """//trim(image_filename)//""""
                call readPng( image_filename,img )
                Nx = size(img,dim=1)
                Ny = size(img,dim=2)
                allocate(img_out(Nx,Ny))
            end if

            if (seed>0) then
                call init_random_seed(seed)
            else
                call init_random_seed()
            end if
            
            img_out=img !just incase we do nothing

        !---    add noise        
        if (noise_stddev>0) then
            write_img=.true.
            do jj = 1,Ny
                img_out(:,jj) = img(:,jj) + gaussianVariate(Nx)*noise_stddev
            end do
        end if

        if (salt+pepper>0) then
            write_img=.true.
            do jj = 1,Ny
                do ii = 1,Nx
                    call random_number(zeta)
                    if (zeta < salt) then
                        img_out(ii,jj) = 1.0d0
                    else if (zeta < salt + pepper) then
                        img_out(ii,jj) = 0.0d0
                    end if
                end do
            end do
        end if

        if(bias) then
            print*,"Calculating bias in background standard deviation caluclation."
            call estimate_bias_in_bg_stddev_calc(img)
        end if

        if (.not.hasArgument(cla,"o")) then
            output_filename= trim(removeSuffix(image_filename))//".out.png"
        end if

        !if(hasArgument(cla,"fftps")) then
        if (get_fft_power_spec) then
            allocate(fftout(0:((2*(Nx/2+1)-1)),0:Ny)) 
            allocate( rpf(0:100) )
            allocate( cprime(0:100) )
            call FFT2d( img,fftout )
            !radialPowerSpectrum2d( in,q_min,q_max, cprime,rpf )
            !call radialPowerSpectrum( img,0.0d0,(Nx/2.0d0), cprime,rpf )
            call radialPowerSpectrum( img,(2*3.141592654d0/max(Nx,Ny)),3.141592654d0, cprime,rpf )
            tempFilename=trim(image_filename)//'.fft.png'
            call write_greyscale_png(tempFilename,fftout,negative=.false.)
        end if

        
        if(write_img) then
            print *,"imageUtilities info - writing image file """//trim(output_filename)//""""
            call write_greyscale_png( output_filename,img_out )
        end if


        contains

        subroutine estimate_bias_in_bg_stddev_calc(img_in)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*  Esitmates bias in the background intneisty standard deviation by integrating x^2p(x) with respect to x in the interval
        !*  0 to t where x is intensity, p(x) is the probability of a given intensity bin (from histogram) and t is the Ridler
        !*  Calvard threshold
            real(kind=real64),dimension(:,:),intent(in)         ::      img_in
            real(kind=real64),dimension(0:255)                  ::      hist
            real(kind=real64),dimension(0:255)                  ::      intensity_bins
            real(kind=real64)                                   ::      bias_estimate
            real(kind=real64)                                   ::      s,b,t,f !ridler clavard image intneisty parameters, background standard devaition, mean background, threshold and mean foreground.
            integer                                             ::      ii,nx,ny,ix,iy
            real(kind=real64)                                   ::      max_fraction_bright_pixels !estimate
            real(kind=real64)                                   ::      mu_0t, var_0t, img_f
            integer                                             ::      noof_px_0t

            nx = size(img_in,dim=1)
            ny = size(img_in,dim=2)
            call findHist(img_in,hist)
            call estimateMaxpixFromDarkPixels(img_in, max_fraction_bright_pixels)
            call findSigma(img_in,max_fraction_bright_pixels,s,b,t,f)

            bias_estimate=0.0d0

            do ii=0,255
                intensity_bins(ii)=ii/255+0.50d0
            end do

            mu_0t=0.0d0
            noof_px_0t=0
            var_0t=0.0d0
            
            do ii=0,(nint(t)) !evalulate integral
                bias_estimate=bias_estimate+(intensity_bins(ii)*intensity_bins(ii)*hist(ii))
            end do

            do ix=1,nx ! go through image 1st to find average vallue of pixels with intensity between 0 and t
                do iy=1,ny
                    img_f=img_in(ix,iy)
                    if (img_f < 0 ) cycle
                    if (img_f > t) cycle
                    mu_0t=mu_0t+img_f
                    noof_px_0t=noof_px_0t+1
                end do
            end do
            if (noof_px_0t>0) then !escape if no pixels over threshold
                mu_0t=mu_0t/noof_px_0t
            else 
                mu_0t=0
            end if
            

            do ix=1,nx ! go through image 2nd time to get variance of pixels between 0 and t
                do iy=1,ny
                    img_f=img_in(ix,iy)
                    if (img_f < 0) cycle
                    if (img_f > t) cycle
                    var_0t=var_0t+((img_f-mu_0t)*(img_f-mu_0t))

                end do
            end do
            if (noof_px_0t>0) then !escape if no pixels over threshold
                var_0t=var_0t/noof_px_0t
            else 
                var_0t=0
            end if


            print*,"image utilities, bg variance =",(s*s)
            print*,"image utilities, variance of pixels with intneisty in range(0,t)=",var_0t
            print*,"image utilities, bg stddev bias estimator =",bias_estimate

            end subroutine estimate_bias_in_bg_stddev_calc

    
        end program imageUtilities
    
    
    
    