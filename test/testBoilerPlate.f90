!------This is just a template for making a program that Ctest can call that can be used to test
!------multiple functions/subroutines, if these test fails then the user must rerun the failed 
!------test verbosely to diagnose the problem.
program testLib_boilerPlate
    !   Describe test
    !   Any other details

        use iso_fortran_env
        use Lib_ColouredTerminal !needed
        use Lib_UtilsForTests
        !use Lib_HERE
        implicit none
        
        logical             ::      ok
        integer             ::      correctcount
        integer             ::      noofTests=0 !
        character(len=32)   ::      libName="testLib_boilerPlate"
        logical             ::      tempCheck
    
        correctcount=0
        
    

        !test 1
        tempCheck=.false.
        !if(TEST_FUNCTION_WITH_CONDITION_HERE) tempCheck=.true.
        call announceSubTest(libName,"functOrSubName",1,noofTests,tempCheck,correctcount)
        !-----------------------------------------------
        ok=haveAllSubTestsPassed(correctcount,noofTests)

        call announcePassOrFail(ok)


    end program testLib_boilerPlate