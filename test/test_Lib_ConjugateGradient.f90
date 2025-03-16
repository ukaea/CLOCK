!------This is just a template for making a program that Ctest can call that can be used to test
!------multiple functions/subroutines, if these test fails then the user must rerun the failed 
!------test verbosely to diagnose the problem.
program testLib_ConjugateGradient
     
    !   Tests the conjugate graident module can return know results within a floating point tolerance.
    !   Any other details

        use Lib_ColouredTerminal !needed
        use Lib_UtilsForTests
        use Lib_ConjugateGradient
        use OMP_LIB

        implicit none
        
        logical             ::      ok
        integer             ::      correctcount
        integer             ::      noofTests=1 !
        character(len=32)   ::      libName="testLib_ConjugateGradient"
        logical             ::      tempCheck

        !conjgrad( A,indx, x,b , eps )
        real(kind=real64),dimension(5,5)         ::      A !solving  A x = b, where A is sparse and symmetric and real
        integer,dimension(0:,:)                  ::      indx
        real(kind=real64),dimension(5)           ::      b
        real(kind=real64),dimension(5)        ::      x
        real(kind=real64),intent(inout)                     ::      eps
        
    
        correctcount=0
        
    

        !test 1
        tempCheck=.false.
        A=reshape((/4.0d0, 1.0d0, 2.0d0, 0.5d0, 2.0d0&
                    1.0d0, 0.5d0, 0.0d0, 0.0d0, 0.0d0&
                    2.0d0, 0.0d0, 3.0d0, 0.0d0, 0.0d0&
                    0.5d0, 0.0d0, 0.0d0, 0.625d0 0.0d0&
                    0.2d0, 0.0d0, 0.0d0, 0.0d0, 16.0d0/),shape(A))
        
        



        !if(TEST_FUNCTION_WITH_CONDITION_HERE) tempCheck=.true.
        call announceSubTest(libName,"functOrSubName",1,noofTests,tempCheck,correctcount)
 
        ok=haveAllSubTestsPassed(correctcount,noofTests)

        call announcePassOrFail(ok)


    end program testLib_ConjugateGradient
         