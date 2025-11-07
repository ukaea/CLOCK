program imageCorrelation
    !---^^^^^^^^^^^^^^^^ 
    !*-----------------------------------------------------------------------------------------------------------------------------------
    !*    ImageCorrelation from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
    !*      James Heath & Daniel Mason
    !*      (c) UKAEA October 2024

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
            use iso_fortran_env
    
            use Lib_CommandLineArguments
            use NBAX_StringTokenizers
            use Lib_Filenames
            use Lib_Png
            use Lib_ImageCorrelationFunction


    
            implicit none
            
    !---
    !   version numbering

    !   0.0.1                   first save

            character(len=*),parameter      ::      VERSION = "0.0.1"
    !---        parameters
            type(CommandLineArguments)      ::      cla                     !Object to hold commmand line argumments (CLA)
        
            character(len=256)              ::      denoised_image_filename = ""
            character(len=256)              ::      gt_image_filename = ""   
            character(len=256)              ::      diff_image_filename = ""
            
    !---        dummy variables  
            real(kind=real64),dimension(:,:),allocatable    ::      img_gt,img_denoised
            real(kind=real64),dimension(:,:,:),allocatable    ::      img_rgb
            logical                                         ::      ok      !does input file exist? 
            integer                                         ::      Nx,Ny,Nx_gt,Ny_gt
            real(kind=real64)                               ::      image_correlation_result,max_error=0.0d0
        !---    read command line arguments
            cla = CommandLineArguments_ctor(45)  
             
            call setProgramDescription( cla, "imageCorrelation" )
            call setProgramVersion( cla, VERSION ) 
           
                        
            
        !---    filename options        
            call get( cla,"f",denoised_image_filename ,LIB_CLA_REQUIRED,"           comparison image filename" )                                                        
            call get( cla,"g",gt_image_filename ,LIB_CLA_REQUIRED,"         ground truth image* filename" ) 
            call get( cla,"diff",diff_image_filename ,LIB_CLA_OPTIONAL,"         image diff filename" ) 
            call get( cla,"max",max_error ,LIB_CLA_OPTIONAL,"         Manually set the max image diff value" )

            call report(cla)
            if (.not. allRequiredArgumentsSet(cla)) stop
            if (hasHelpArgument(cla)) stop
            call delete(cla)
            
            
            if (len_trim(denoised_image_filename)/=0) then    
                inquire(file=trim(denoised_image_filename),exist=ok)
                if (.not. ok) then
                    print *,"imageCorrelation  error - couldn't find image file """//trim(denoised_image_filename)//""""
                    stop
                end if
                print *,"imageCorrelation info - reading image file """//trim(denoised_image_filename)//""""
                call readPng( denoised_image_filename,img_denoised )
                Nx = size(img_denoised,dim=1)
                Ny = size(img_denoised,dim=2)
            end if

            if (len_trim(gt_image_filename)/=0) then    
                inquire(file=trim(gt_image_filename),exist=ok)
                if (.not. ok) then
                    print *,"imageCorrelation error - couldn't find image file """//trim(gt_image_filename)//""""
                    stop
                end if
                print *,"imageCorrelation info - reading image file """//trim(gt_image_filename)//""""
                call readPng( gt_image_filename,img_gt )
                Nx_gt = size(img_gt,dim=1)
                Ny_gt = size(img_gt,dim=2)
            end if

            if((Nx.ne.Nx_gt).or.(Ny.ne.Ny_gt)) then
                print *,"imageCorrelation error - The images that you are trying to compare must have the same dimensions" 
                print*,"yours do not, have you got the right images?"
                stop
            end if

            image_correlation_result=imageCorrelationFunction( img_denoised,img_gt)

            print*,"image correlation is:",image_correlation_result

            if (len_trim(diff_image_filename)/=0) then
                print*,"Calculating image diffs"
                ! if (max_error>0) then
                !     call get_image_diff(img_denoised,img_gt,img_rgb,max_error)  
                ! else
                !     call get_image_diff(img_denoised,img_gt,img_rgb,max_error)
                ! end if
                call get_image_diff(img_denoised,img_gt,img_rgb,max_error)
                call write_rgb_png(diff_image_filename,img_rgb)  
            end if

end program imageCorrelation




            