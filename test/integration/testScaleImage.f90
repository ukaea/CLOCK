!------This is just a test for the scaleImage program that Ctest can call that can be used to test
!------multiple functions/subroutines, if these test fails then the user must rerun the failed 
!------test verbosely to diagnose the problem.
program testScaleImage
    !   Describe test
    !   Any other details

        use iso_fortran_env
        use Lib_ColouredTerminal !needed
        use Lib_UtilsForTests
        use Lib_Png
        
        !use Lib_HERE
        implicit none
        
        logical             ::      ok
        integer             ::      correctcount
        integer             ::      noofTests=2 !
        character(len=32)   ::      libName="testScaleImage"
        logical             ::      tempCheck
        real(kind=real64),dimension(2,2)    ::      make_me_big
        real(kind=real64),dimension(4,4)    ::      make_me_little
        real(kind=real64),dimension(:,:),allocatable    ::      made_you_little 
        real(kind=real64),dimension(:,:),allocatable    ::      made_you_big
        real(kind=real64)   ::      tolerance=5e-4 !large as the test images are tiny
        real(kind=real64)   ::      value_check=huge(1.0d0)
    
        correctcount=0

        !set up test data
        make_me_big=reshape((/1.0d0,0.0d0,1.0d0,0.0d0/),(/2,2/))
        make_me_little=reshape((/   1.0d0,1.0d0,0.0d0,0.0d0, &
                                    1.0d0,1.0d0,0.0d0,0.0d0, &
                                    0.0d0,0.0d0,1.0d0,1.0d0, &
                                    0.0d0,0.0d0,1.0d0,1.0d0/),(/4,4/))

        ! call writePng("../../../data/make_me_big.png", make_me_big)
        ! call writePng("../../../data/make_me_little.png", make_me_little)
        call writePng("make_me_big.png", make_me_big)
        call writePng("make_me_little.png", make_me_little)
    
        !call system("./../../src/MLE_Filter/MLEFilter -f ../../../data/shot_test_clean.png")

        !test 1, can we downscale an image
        tempCheck=.false.
        !call system("./../../src/Tools/scaleImage -f ../../../data/make_me_little.png -n 1 -o ../../../data/made_you_little.png")
        call system("./../../src/Tools/scaleImage -f make_me_little.png -n 1 -o made_you_little.png")
        !call readPng("../../../data/made_you_little.png",made_you_little)
        call readPng("made_you_little.png",made_you_little)
        print*,"orig large"
        print*,make_me_little
        print*,"downscaled small"
        print*,made_you_little

        print*,"sum of diff",abs(sum(made_you_little-make_me_big))
        value_check=abs(sum(made_you_little-make_me_big))
        if(value_check<tolerance) tempCheck=.true.
        call announceSubTest(libName,"test downscaling",1,noofTests,tempCheck,correctcount)

        !test 2, can we downscale an image
        tempCheck=.false.
        !call system("./../../src/Tools/scaleImage -f ../../../data/make_me_big.png -nos -n 1 -o ../../../data/made_you_big.png")
        !call readPng("../../../data/made_you_big.png",made_you_big)
        call system("./../../src/Tools/scaleImage -f make_me_big.png -nos -n 1 -o made_you_big.png")
        call readPng("made_you_big.png",made_you_big)
        print*,"orig small"
        print*,make_me_big
        print*,"upscaled large "
        print*,made_you_big

        print*,"sum of diff",abs(sum(made_you_big-make_me_little))
        value_check=abs(sum(made_you_big-make_me_little))
        if(value_check<tolerance) tempCheck=.true.
        call announceSubTest(libName,"test upscaling",2,noofTests,tempCheck,correctcount)
        !-----------------------------------------------
        ok=haveAllSubTestsPassed(correctcount,noofTests)

        call announcePassOrFail(ok)


    end program testScaleImage