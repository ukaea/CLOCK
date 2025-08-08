!------Integration tests for the imageCorrelation program that Ctest can call that can be used to test
!------multiple functions/subroutines, if these test fails then the user must rerun the failed 
!------test verbosely to diagnose the problem.
program testImageCorrelation
    !   Describe test
    !   Any other details

        use iso_fortran_env
        use Lib_ColouredTerminal !needed
        use Lib_UtilsForTests
        use Lib_ImageCorrelationFunction
        use pipes_module
        use Lib_Png
        implicit none
        
        logical             ::      ok
        integer             ::      correctcount
        integer             ::      noofTests=3 !
        character(len=32)   ::      libName="testImageCorrelation"
        logical             ::      tempCheck
        character(len=:),allocatable :: res
        real(kind=real64)           ::  num
    
        correctcount=0
        
    

        !test 1, compare image with itself, should return 1.0
        tempCheck=.false.
        res = get_command_as_string('./../../src/ImageCorrelation/imageCorrelation -f ../../../data/fe6cr_n_irrad_1.8dpa_fib_dam_BF.png -g ../../../data/fe6cr_n_irrad_1.8dpa_fib_dam_BF.png | grep " is:" | awk '//"'"//'{print $4}'//"'")
        num=str2real(res)
        print*,"expected: 1.0, result:",num
        if(num==1.0d0) tempCheck=.true.
        call announceSubTest(libName,"compare identical",1,noofTests,tempCheck,correctcount)

        !test 2, compare darkfield image with bright field image of same dimensions, should return <0 due to small negative corrlation
        tempCheck=.false.
        res = get_command_as_string('./../../src/ImageCorrelation/imageCorrelation -f ../../../data/fe6cr_n_irrad_1.8dpa_fib_dam_BF.png -g ../../../data/fe6cr_n_irrad_1.8dpa_fib_dam_DF.png | grep " is:" | awk '//"'"//'{print $4}'//"'")
        num=str2real(res)
        print*,"expected: <0.0, result:",num
        if(num<0.0d0) tempCheck=.true.
        call announceSubTest(libName,"compare differnet BF/DF",2,noofTests,tempCheck,correctcount)

        !test 3, attempt to compare images with different dimension, check for correct error message
        tempCheck=.false.
        !res = get_command_as_string('./../../src/ImageCorrelation/imageCorrelation -f ../../../data/test_cameraman.png -g ../../../data/test_cat.png | grep " is:" | awk '//"'"//'{print $1}'//"'")
        res = get_command_as_string('./../../src/ImageCorrelation/imageCorrelation -f ../../../data/test_cameraman.png -g ../../../data/test_cat.png| grep " you got"')
        print*,"expected:   yours do not, have you got the right images?"
        print*,"got:",res
        if(lle(res,' yours do not, have you got the right images?')) tempCheck=.true.
        call announceSubTest(libName,"compare differnet BF/DF",3,noofTests,tempCheck,correctcount)

        !test 4, check diff print works "-diff" TODO

        !test 5, check if diff print  "-max" TODO
        !-----------------------------------------------
        ok=haveAllSubTestsPassed(correctcount,noofTests)

        call announcePassOrFail(ok)


    end program testImageCorrelation