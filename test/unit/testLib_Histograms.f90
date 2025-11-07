!------This is a test suite for Lib_Histograms that Ctest can call that can be used to test
!------multiple functions/subroutines, if these test fails then the user must rerun the failed 
!------test verbosely to diagnose the problem.
program testLib_Histograms

        use iso_fortran_env
        use Lib_ColouredTerminal !needed
        use Lib_UtilsForTests
        use Lib_Histograms
        implicit none
        
        logical             ::      ok
        integer             ::      correctcount
        integer             ::      noofTests=3 !
        character(len=32)   ::      libName="testLib_Histograms"
        logical             ::      tempCheck

        real(kind=real64),dimension(100)        ::      testData,x
        real(kind=real64)                       ::      f,m,s!g(x) = (f^2) Exp[ -(x-m)^2 / s^2 ]
        real(kind=real64)                       ::      tolerance !floating point tolerance.
        real(kind=real64)                       ::      dx !bin widths with a nice value
        integer                                 ::      ii, nBinsOut

    
        correctcount=0
        

        tolerance=1.0e-8 !set floating point tolerance

        f=5.0d0 !generate some normally distirbuted data to do some tests on.
        m=10.0d0
        s=11.0d0
        do ii =1,100
            x(ii)=ii/5.0d0
            testData(ii)=(f*f)*exp((-(x(ii)-m)*(x(ii)-m))/(s*s))
        end do

        !test 1, Can we suggest a sensible range of bin values for normally distirbuted data?
        tempCheck=.false.

        if(automatic_xrange(testData)==30.0d0) tempCheck=.true.
        call announceSubTest(libName,"automatic_xrange",1,noofTests,tempCheck,correctcount)

        !test 2, Can we suggest a suitable bin width for a histogram of normally distirbuted data?
        tempCheck=.false.
        if((automatic_binwidth(testData)- 3.5036187498017330d0)<tolerance) tempCheck=.true.
        call announceSubTest(libName,"automatic_binwidth",2,noofTests,tempCheck,correctcount)

        !test 3, can we suggest a sensible number of bins with a "nice" bin width for a histogram of normally distirbuted data?
        tempCheck=.false. 
        dx=3.5036187498017330d0
        nBinsOut=0.0d0
        call automatic_nBins(30.0d0,dx,nBinsOut)
        if((dx==3.0d0).and.(nBinsOut==10)) tempCheck=.true.
        call announceSubTest(libName,"automatic_nBins",3,noofTests,tempCheck,correctcount)
        !-----------------------------------------------
        ok=haveAllSubTestsPassed(correctcount,noofTests)

        call announcePassOrFail(ok)


    end program testLib_Histograms