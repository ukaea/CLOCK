program scaleImage
    !simple program for interpolating an image over more/less pixels
    !---^^^^^^^^^^^^^^^^ 
    !*-----------------------------------------------------------------------------------------------------------------------------------
    ! *    scaleImage from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
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
            use iso_fortran_env
    
            use Lib_CommandLineArguments
            use NBAX_StringTokenizers
            use Lib_Filenames
            use Lib_ColourScale
            use Lib_Png
            use Lib_SimpleProgressBar
            use Lib_DrawEllipse
            use Lib_DataDoubler
    
            implicit none

        !   version numbering

        !   0.0.1                   first save
        !---        input options
            character(len=*),parameter      ::      VERSION = "0.0.1"
            type(CommandLineArguments)      ::      cla                     !Object to hold commmand line argumments (CLA)
        
            character(len=256)              ::      input_filename = ""         !   filename of image to down scale
            character(len=256)              ::      output_filename = "scaled_image.png"        !   filename of output image
            logical                         ::      shrink=.true.               !   If ture downscale, otherwise up scale
            integer                         ::      noof_scales=0               !   How many times to half or double the image size
            integer                         ::      ii,jj                       !   image indices
            integer                         ::      Nx,Ny,Mx,My                 !   XY Dimesions of the input and scaled image respectively
            logical                         ::      ok                          !   For checks
            real(kind=real64),dimension(:,:),allocatable    ::      img_in      !   greyscale image to scale
            real(kind=real64),dimension(:,:),allocatable    ::      img_out     !   Scaled greyscale image to output
            real(kind=real64),dimension(:,:),allocatable    ::      img_tmp     !   temp array to work on            
            
        !---    read command line arguments
            cla = CommandLineArguments_ctor(30)  
             
            call setProgramDescription( cla, "scaleImage" )
            call setProgramVersion( cla, VERSION )
            
            
        !---        input options
            call get( cla,"f",input_filename ,LIB_CLA_REQUIRED,"           Input image filename" )                                                        
            call get( cla,"o",output_filename ,LIB_CLA_OPTIONAL,"         Output scaled image filename" )
            call get( cla,"n",noof_scales ,LIB_CLA_OPTIONAL,"         How many times to half/double" )
            call get( cla,"s",shrink ,LIB_CLA_OPTIONAL,"         Shink if true, double if false" )

            call report(cla)
            if (.not. allRequiredArgumentsSet(cla)) stop
            if (hasHelpArgument(cla)) stop
            call delete(cla)

                        !---    check that the file exists        
            inquire(file=trim(input_filename),exist=ok)
            if (.not. ok) then
                 print *,"scaleImage error - couldn't find image """//trim(input_filename)//""""
                 stop
            end if

            if(noof_scales<=0) then
                 print *,"scaleImage error -n<=0, no scaling to be done"
                 stop
            end if
            print *,"scaleImage info - reading image file """//trim(input_filename)//""""
            call readPng( input_filename,img_in )
            Nx = size(img_in,dim=1)
            Ny = size(img_in,dim=2)

            allocate(img_tmp(0:Nx-1,0:Ny-1))

            Mx=Nx
            My=Ny
            img_tmp=img_in
            if(shrink) then
                do ii=1,noof_scales
                    Mx = int((Mx+1)/2) ; My = int((My+1)/2) !get new dowscaled image dims
                    allocate(img_out(0:Mx-1,0:My-1))    !allocate to hold the downscaled image
                    call halfImage( (size(img_tmp,dim=1)),(size(img_tmp,dim=2)),img_tmp,img_out, pbc=.false. )!down scale it
                    deallocate(img_tmp)!dealloacte tmp so we can ...
                    allocate(img_tmp(0:Mx-1,0:My-1)) !...reallocate it for the smaller dimensions
                    img_tmp=img_out !switch so the argument order is preserved and...
                    deallocate(img_out) !... out deallocated so it can be reallocated in the next loop
                end do
                allocate(img_out(0:Mx-1,0:My-1))
                img_out=img_tmp
            else
                do ii=1,noof_scales
                    Mx = int((Mx)*2) ; My = int((My)*2) !get new upscaled image dims
                    allocate(img_out(0:Mx-1,0:My-1))
                    call doubleImage( 2*(size(img_tmp,dim=1)),2*(size(img_tmp,dim=2)),img_tmp,img_out, pbc=.false. )
                    deallocate(img_tmp)
                    allocate(img_tmp(0:Mx-1,0:My-1))
                    img_tmp=img_out
                    deallocate(img_out)
                end do
                allocate(img_out(0:Mx-1,0:My-1))
                img_out=img_tmp
            end if

              
            print*,"Wrinting image"
            call writePng(output_filename, img_out)


end program scaleImage