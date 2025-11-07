!------This is just a test for the maximum likelhood filter (MLEFilter.f90)
!------Attempts to denoise an image with and without added shot noise and checks  
!------image correlation of ground truth (GT) and filtered fall within tolerance. 
program testMLE_Filter

        use iso_fortran_env
        use Lib_ColouredTerminal !needed
        use Lib_UtilsForTests
        use Lib_Png
        !use Lib_HERE
        implicit none
        
        logical             ::      ok
        integer             ::      correctcount
        integer             ::      noofTests=2 !
        character(len=32)   ::      libName="testMLE_Filter"
        logical             ::      tempCheck
        real(kind=real64),dimension(:,:), allocatable   ::     gt_img, noisy_img, filtered_noisy_img, filtered_gt_img
        real(kind=real64)   ::      dic_test1,dic_test2 

        !to make noisless image: ./build/bin/testImageGenerator -f data/shot_test_clean.png -Nx 256 -Ny 256 -seed 123 -count 50 -Rmu 8 -Rsig 0 -Isig 0.0 -Imu 0.25 -b 0.3 -sinb 0.1
        !to make noisy image:./build/bin/testImageGenerator -f data/shot_test_noisy.png -Nx 256 -Ny 256 -seed 123 -count 50 -Rmu 8 -Rsig 0 -Isig 0.0 -Imu 0.25 -b 0.3 -sinb 0.1 -stddev 0.05
    
        correctcount=0

        
        !call system("./../../src/MLE_Filter/MLEFilter -f ../../../data/shot_test_noisy.png")
        
    

        !test 1, if a noisless image f(x) is denoised with the MLE filter, is DIC(GT,f(x)) <1e-6
        tempCheck=.false.
        dic_test1=0.0d0
        call system("./../../src/MLE_Filter/MLEFilter -f ../../../data/shot_test_clean.png")
        call readPng( "../../../data/shot_test_clean.png",gt_img )
        call readPng( "../../../data/shot_test_clean.mle.png",filtered_gt_img )
        dic_test1=1.0d0-dic(gt_img,filtered_gt_img)
        print*,"DIC val is",dic_test1,"should be less than 1e-3."
        if(dic_test1<1.0e-3) tempCheck=.true.
        call announceSubTest(libName,"MLE Filter on noisless image",1,noofTests,tempCheck,correctcount)

        !test 2, if a noisy image f(x) is denoised with the MLE filter, is DIC(GT,f(x)) <1e-2
        tempCheck=.false.
        dic_test2=0.0d0
        call system("./../../src/MLE_Filter/MLEFilter -f ../../../data/shot_test_noisy.png")
        call readPng( "../../../data/shot_test_noisy.mle.png",filtered_noisy_img )
        dic_test2=1.0d0-dic(gt_img,filtered_noisy_img)
        print*,"DIC val is",dic_test2,"should be less than 1e-2."
        if(dic_test2<1.0e-2) tempCheck=.true.
        call announceSubTest(libName,"MLE Filter on noisy image",1,noofTests,tempCheck,correctcount)
        !-----------------------------------------------
        ok=haveAllSubTestsPassed(correctcount,noofTests)

        call announcePassOrFail(ok)

        contains

                pure real(kind=real64) function dic( f,g )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute dic = sum (f-<f>)(g-<g>) / sqrt( sum(f-<f>)^2 sum(g-<g>)^2 )
            real(kind=real64),dimension(:,:),intent(in)       ::      f,g
            integer             ::      nn,nx,ny,ii,jj
            real(kind=real64)   ::      fsum,gsum,f2sum,g2sum,fgsum
            real(kind=real64)   ::      fbar,gbar
            

            nx = size(f,dim=1)
            ny = size(f,dim=2)
            nn = nx*ny
            fsum = 0
            gsum = 0
            f2sum = 0
            g2sum = 0
            fgsum = 0
            do jj = 1,ny
                do ii = 1,nx
                    fsum    = fsum  + f(ii,jj)
                    f2sum   = f2sum + f(ii,jj)*f(ii,jj)
                    gsum    = gsum  + g(ii,jj)
                    g2sum   = g2sum + g(ii,jj)*g(ii,jj)
                    fgsum   = fgsum + f(ii,jj)*g(ii,jj)
                end do
            end do
            fbar = fsum/nn
            gbar = gsum/nn

           ! print *,nn,fsum,gsum,fbar,gbar,f2sum,g2sum,fgsum

            dic = (f2sum - nn*fbar*fbar)*(g2sum - nn*gbar*gbar)
            if (dic <= 0) then      !   <0 for -0 error
                dic = 0
                return
            end if

            dic = (fgsum - nn*fbar*gbar) / sqrt(dic)

            return
        end function dic


    end program testMLE_Filter