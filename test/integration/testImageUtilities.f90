!------This is  a test for iamgeUtilities.f90 that Ctest can call that can be used to test
!------multiple functions/subroutines, if these test fails then the user must rerun the failed 
!------test verbosely to diagnose the problem.
program testLib_boilerPlate
    !   Describe test
    !   Any other details

        use iso_fortran_env
        use Lib_ColouredTerminal !needed
        use Lib_UtilsForTests
        use Lib_Png
        implicit none
        
        logical             ::      ok
        integer             ::      correctcount
        integer             ::      noofTests=2 !
        character(len=32)   ::      libName="testIamgeUtilities"
        logical             ::      tempCheck

        real(kind=real64)   ::      max_pix,min_pix
        real(kind=real64),dimension(:,:), allocatable   ::  fft_read, noisy_read
        real(kind=real64)   ::      tolerance=1e-6
    
        correctcount=0
        
    
        !run the command & get the outputs
        call system("./../../src/Tools/imageUtilities -f ../../../data/test_grey.png -o ../../../data/test_grey_imut.png -seed 123 -stddev 0.2 -salt 1e-3 -pepper 1e-3 -bias -fftps")
        call readPng("../../../data/test_grey_imut.png", noisy_read)
        call readPng("../../../data/test_grey.png.fft.png", fft_read)
        
        !test 1, has salt & pepper noise been added
        tempCheck=.false.
        max_pix=maxval(noisy_read)
        min_pix=minval(noisy_read)
        print*,"S&P max val:",max_pix,"min val:",min_pix
        if(((max_pix-min_pix)-1.0d0)<tolerance) tempCheck=.true.
        call announceSubTest(libName,"Salt & pepper addition",1,noofTests,tempCheck,correctcount)

        !test 2, is the FFT of the flat image correct (single pixel white pixel in center)
        tempCheck=.false.
        max_pix=0.0d0
        max_pix=fft_read(513,513)
        print*,"val off fft at (512,512) is:",max_pix
        if((max_pix-1.0d0)<tolerance) tempCheck=.true.
        call announceSubTest(libName,"FFT calc",2,noofTests,tempCheck,correctcount)

        !-----------------------------------------------
        ok=haveAllSubTestsPassed(correctcount,noofTests)

        call announcePassOrFail(ok)


    end program testLib_boilerPlate