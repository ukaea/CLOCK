    program test_Lib_Heteroscedastic_Fit
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!TODO: test needs to check values.
        use Lib_Heteroscedastic_Fit
        use Lib_RandomSeed
        use iso_fortran_env
        implicit none


        integer             ::      N = 100         !   number of points
        real(kind=real64),dimension(:),allocatable      ::      x,y         !   points
        real(kind=real64)   ::      a,b,c,d         !   distribution actually used to generate random points
   
        integer             ::      ii  
        real(kind=real64)   ::      sig,zeta,xx,yy,eta
        real(kind=real64)   ::      aa,bb,cc,dd
        real(kind=real64)   ::      tolerance=1e-6

        call init_random_seed(12345)
        a = 1 ; b = 0.5
        c = 0.5 ; d = 0.4
        allocate(x(N))
        allocate(y(N))
        call random_number(x) ; x = x * 10
      
      
        do ii = 1,N
            xx = x(ii)
            sig = c*xx + d
            zeta = gaussianVariate() 
            
            call random_number(eta)
            if (eta<0.1) then
                zeta = zeta * 10
            end if
            yy = a*xx + b + zeta * (c*xx + d)

            
                y(ii) = yy
                x(ii) = xx
                !print *,ii,x(ii),(a*x(ii) + b),(c*x(ii) + d),y(ii)                    
   
           
        end do
        
        LIB_HET_FIT_DBG = .true.
        print *,""
        call Heteroscedastic_Fit( x,y , aa,bb,cc,dd , sig_min=1.0d-6)
        print *,"fit result y = a x + b , sigma = c x + d "
        print *,"           a = ",a,aa
        print *,"           b = ",b,bb
        print *,"           c = ",c,cc
        print *,"           d = ",d,dd



        print *,""
        print *,"PASS"
        print *,""


    end program test_Lib_Heteroscedastic_Fit