
    module Lib_Heteroscedastic_Fit
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_Heteroscedastic_Fit from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
!*      Daniel Mason
!*      (c) UKAEA June 2024

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
!*      A simple module to fit data {xi,yi}
!*      to a line y = a x + b
!*      where the errors are _not_ constant but are a linear function of x
!*      ie sigma = c x + d
!*
!*      then the probability of observing (xi,yi) is
!*          pi = Exp[ - (yi-y)^2/(2 s^2) ] /  Sqrt[ 2 pi s^2 ]
!*      where y = a xi + b, s = c xi + d
!*        
!*      For c = 0, this reduces to the simple linear least squares fit, which can be fitted analytically 
!*      but for c/=0, we must make a guess for a,b,c,d, then find the log-likelihood
!*          Lambda = sum_i Log[ pi ]      
!*      and maximise wrt {a,b,c,d}
!*      
!*      Fortunately everything stays analytic and fairly simple, so we can find both
!*          d Lambda / d u , d^2 Lambda / du du     with u = {a,b,c,d}
!*      and then maximise by solving
!*
!*          d^2 Lambda / du du delta u + d Lambda / du = 0
!*      to give an update change in u={a,b,c,d}
!*
!*
        !use Lib_SafeExp        !   note: should change to safeExp
        use iso_fortran_env
        implicit none
        private

        external            ::      DSYSV           !   lapack symmetric solve

        real(kind=real64),private,parameter         ::      PI = 3.141592653590d0
        logical,public                              ::      LIB_HET_FIT_DBG = .false.

        public      ::      Heteroscedastic_Fit
        

    contains
!---^^^^^^^^

        subroutine Heteroscedastic_Fit( x,y , a,b,c,d , sig_min)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            real(kind=real64),dimension(:),intent(in)               ::      x
            real(kind=real64),dimension(:),intent(in)               ::      y
            real(Kind=real64),intent(out)                           ::      a,b     !   equation of best fit line
            real(Kind=real64),intent(out)                           ::      c,d     !   equation of error bars
            real(Kind=real64),intent(in)                            ::      sig_min    !    minimum error - ensure c x + d >= sig_min, the minimum permitted error
            
            real(kind=real64),dimension(size(x))        ::      ff
            integer                                     ::      step,loop
            real(kind=real64)                           ::      Lambda,oldLambda
            integer                                     ::      nPoints,ii
             
        !---    quick escape for a small number of points
            nPoints = size(x,dim=1)
            if (nPoints == 0) then
                !   problem - no points ???
                print *,"Lib_Heteroscedastic_Fit::Heteroscedastic_Fit error - trying to fit a line with no points"
                a = 0 ; b = 0 ; c = 0 ; d = sig_min
                return
            else if (nPoints == 1) then
                !   OK, I can fit an equation if I include the origin, but I can't fit the error 
                print *,"Lib_Heteroscedastic_Fit::Heteroscedastic_Fit warning - trying to fit a line with one point"
                if (x(1)/=0) then
                    a = y(1)/x(1) ; b = 0 ; c = 0 ; d = sig_min
                else
                    print *,"Lib_Heteroscedastic_Fit::Heteroscedastic_Fit warning - trying to fit a line with one point with x(1)=0"
                    a = 0 ; b = 0 ; c = 0 ; d = sig_min
                end if
                return
            else if (nPoints == 2) then
                !   OK, I can fit an equation, but I can't fit the error 
                print *,"Lib_Heteroscedastic_Fit::Heteroscedastic_Fit warning - trying to fit an error trend with two points"
                if (x(1)/=x(2)) then
                    a = (y(2) - y(1))/(x(2) - x(1)) 
                    b = (y(1)*x(2) - y(2)*x(1))/(x(2) - x(1)) 
                    c = 0 ; d = sig_min
                else
                    print *,"Lib_Heteroscedastic_Fit::Heteroscedastic_Fit warning - trying to fit a line with two points with x(1) = x(2)"
                    a = 0 ; b = 0 ; c = 0 ; d = sig_min
                end if
                return
            end if


                






        !   make an initial guess assuming all data points are good    
            ff = 1
            call Heteroscedastic_LinearFit( x,y , ff , a,b,c,d )
            oldLambda = computeLambda( x,y , ff , a,b,c,d , sig_min ) /nPoints
            
            !oldLambda = oldLambda/nPoints
            if (LIB_HET_FIT_DBG) &       
            write (*,fmt='(a,i8,a,f16.6,a,i8,a,g16.6,a,g16.6)') "Lib_Heteroscedastic_Fit::Heteroscedastic_Fit info - loop ",0," point density ",sum(ff),"/",nPoints," Lambda ",oldLambda," rss ",computeRss(x,y,ff,a,b)/sum(ff)
                
        !---    improve the fit by chucking out wild points
            if (npoints>3) then !need more than one data point to fit.
                do loop = 1,1
                    !d = d/2
                    do step = 1,1000
                        call findDodgyPoints( x,y,a,b,c,d, 4.0d0-loop, ff , sig_min )  
                        call fit( x,y , ff , a,b,c,d , Lambda , sig_min )   
                        Lambda = Lambda / sum(ff)
                        !print *,"a,b,c,d ",a,b,c,d
                        if (abs(Lambda - oldLambda) < abs(Lambda)*1.0d-8) then
                            !   converged
                            exit
                        end if
                        oldLambda = Lambda
                    end do
                    if (LIB_HET_FIT_DBG) &       
                    write (*,fmt='(a,i8,a,f16.6,a,i8,a,g16.6,a,g16.6)') "Lib_Heteroscedastic_Fit::Heteroscedastic_Fit info - loop ",loop," point density ",sum(ff),"/",nPoints," Lambda ",Lambda," rss ",computeRss(x,y,ff,a,b)/sum(ff)
                end do
            end if

            if (LIB_HET_FIT_DBG) then
                call findDodgyPoints( x,y,a,b,c,d, 1.0d0, ff , sig_min )  
                write (*,fmt='(a8,100a16)') "point ","x","y","a x + b","c x + d","f"
                do ii = 1,nPoints
                    write (*,fmt='(i8,100f16.8)') ii, x(ii),y(ii),a*x(ii)+b,c*x(ii)+d,ff(ii)
                end do
            end if


            return
        end subroutine Heteroscedastic_Fit 


        subroutine Heteroscedastic_LinearFit( x,y , f , a,b,c,d  )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*
            real(kind=real64),dimension(:),intent(in)               ::      x
            real(kind=real64),dimension(:),intent(in)               ::      y
            real(kind=real64),dimension(:),intent(in)               ::      f       !   logistic function goodness of fit = 1 for on line, =0 for not on line
            real(Kind=real64),intent(out)                           ::      a,b     !   equation of best fit line
            real(Kind=real64),intent(out)                           ::      c,d     !   equation of error bars
            !real(kind=real64),intent(out)                           ::      Lambda  !   log likelihood
            integer                 ::      nPoints
            real(kind=real64)       ::      sumx,sumxx,sumxy,sumy,sumf !,sumyy
            real(kind=real64)       ::      res,det
            integer                 ::      ii

            nPoints = size(x,dim=1)

            if (nPoints <= 1) then
                !   can't even solve for a straight line!
                a = 0 ; b = 0 ; c = 0 ; d = 0 
                return
            end if


        !---    make a guess for the equation of the line y = a x + b assuming c = 0, d = 1
        !       this is standard linear least squares. See eg wikipedia.
            c = 0 ; d = 1
            sumx = 0 ; sumxx = 0 ; sumxy = 0 ; sumy = 0  ; sumf = 0 !; sumyy = 0
            do ii = 1,nPoints
                sumx = sumx + f(ii)*x(ii)
                sumxx = sumxx + f(ii)*x(ii)*x(ii)
                sumy = sumy + f(ii)*y(ii)
                sumxy = sumxy + f(ii)*x(ii)*y(ii)
                !sumyy = sumyy + f(ii)*y(ii)*y(ii)
                sumf = sumf + f(ii)
            end do
            det = sumxx*sumf - sumx*sumx
            if (det <= abs(sumxx)*1.0d-8) then
                !   all the points are on top of each other!!
                a = 0 ; b = 0 ; c = 0 ; d = 0 
                return
            end if
            a = (sumxy*sumf - sumy*sumx)/det
            b = (sumy*sumxx - sumxy*sumx)/det

           ! print *,"Heteroscedastic_LinearFit ",sumx,sumxx,sumy,sumxy,sumf

        !---    make a guess for d given the equation of the line y=ax+b, still keeping c = 0?
        !   d = sqrt(max(0.0d0, ( sumyy - 2*a*sumxy - 2*b*sumy + a*a*sumxx + 2*a*b*sumx + b*b )/nPoints ))


        !---    make a guess for c,d given the distribution of the absolute residuals.
        !       note that I am saying abs(res) ~ c x + d as an estimate for sigma ~ c x + d
            sumxy = 0 ; sumy = 0  
            do ii = 1,nPoints
                res = abs( y(ii) - (a*x(ii)+b) )                
                sumy = sumy + f(ii)*res
                sumxy = sumxy + f(ii)*x(ii)*res
            end do
            c = (sumxy*sumf - sumy*sumx)/det
            d = (sumy*sumxx - sumxy*sumx)/det

            if (LIB_HET_FIT_DBG) &       
                write (*,fmt='(4(a,f16.6))') "Lib_Heteroscedastic_Fit::Heteroscedastic_Fit info - y = a x + b , ( a,b = ",a,",",b," )"
           ! print *,"Heteroscedastic_LinearFit ",a,b,c,d
        !---    improve the guess
            !call fit( x,y , f , a,b,c,d , Lambda , sig_min)

            return
        end subroutine Heteroscedastic_LinearFit

            
        subroutine fit( x,y , f , a,b,c,d , Lambda , sig_min)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*
            real(kind=real64),dimension(:),intent(in)               ::      x
            real(kind=real64),dimension(:),intent(in)               ::      y
            real(kind=real64),dimension(:),intent(in)               ::      f       !   logistic function goodness of fit = 1 for on line, =0 for not on line
            real(Kind=real64),intent(inout)                         ::      a,b     !   equation of best fit line
            real(Kind=real64),intent(inout)                         ::      c,d     !   equation of error bars
            real(Kind=real64),intent(in)                            ::      sig_min    !    ensure c x + d >= sig_min, the minimum permitted error
            real(kind=real64),intent(out)                           ::      Lambda  !   log likelihood
                                   
            real(kind=real64),dimension(4)          ::      dLambdadu
            real(kind=real64),dimension(4,4)        ::      d2Lambdadudu
            integer,dimension(4)                    ::      ipiv
            integer,dimension(66*4)                 ::      work
            real(kind=real64)       ::      oldLambda
            integer                 ::      step

            real(kind=real64),dimension(0:7),parameter    ::      delta =  (/ 0.0d0,0.1d0,0.5d0,0.9d0,1.0d0,2.0d0,3.0d0,5.0d0 /)
            real(kind=real64)       ::      bestLambda
            integer                 ::      bestDelta
            integer                 ::      jj

            !print *,"fit ",a,b,c,d
            oldLambda = computeLambda( x,y , f , a,b,c,d , sig_min )
            
            do step = 1,1000

                call computeDerivs( x,y , f , a,b,c,d , dLambdadu, d2Lambdadudu )

                call DSYSV("L",4,1,d2Lambdadudu,4,ipiv,dLambdadu,4,work,size(work),jj )

                bestLambda = oldLambda ; bestDelta = 0
                dLambdadu = -dLambdadu
                do jj = 1,ubound(delta,dim=1)
                    Lambda = computeLambda( x,y , f , a+delta(jj)*dLambdadu(1),b+delta(jj)*dLambdadu(2),c+delta(jj)*dLambdadu(3),d+delta(jj)*dLambdadu(4) , sig_min )
                    !print *,"test ",jj,delta(jj),Lambda
                    if (Lambda > bestLambda) then
                        bestLambda = Lambda
                        bestDelta = jj
                    end if
                end do

                a = a + delta(bestDelta)*dLambdadu(1)
                b = b + delta(bestDelta)*dLambdadu(2)   
                c = c + delta(bestDelta)*dLambdadu(3)   
                d = d + delta(bestDelta)*dLambdadu(4)   
                

                if (bestLambda - oldLambda < abs(bestLambda)*1.0d-8) then
                    !   converged
                    return
                end if
                oldLambda = bestLambda
                !stop
            end do

            Lambda = bestLambda
            return

        contains
    !---^^^^^^^^

            subroutine computeDerivs( x,y , f , a,b,c,d , dLambdadu, d2Lambdadudu )
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                real(kind=real64),dimension(:),intent(in)               ::      x
                real(kind=real64),dimension(:),intent(in)               ::      y
                real(kind=real64),dimension(:),intent(in)               ::      f       !   logistic function goodness of fit = 1 for on line, =0 for not on line
                real(Kind=real64),intent(in)                            ::      a,b     !   equation of best fit line
                real(Kind=real64),intent(in)                            ::      c,d     !   equation of error bars
                real(Kind=real64),dimension(4),intent(out)              ::      dLambdadu
                real(Kind=real64),dimension(4,4),intent(out)            ::      d2Lambdadudu
                integer                 ::      nPoints
                integer                 ::      ii
            
                real(Kind=real64)       ::      axpb,cxpd,ymaxpb,icxpd

                nPoints = size(x,dim=1)
                Lambda = 0
                dLambdadu = 0
                d2Lambdadudu = 0
                do ii = 1,nPoints
                    axpb = a*x(ii) + b
                    cxpd = c*x(ii) + d
                    ymaxpb = y(ii) - axpb
                    icxpd = 1/max(1.0d-16,cxpd)
 
                    dLambdadu(1) = dLambdadu(1) + f(ii)*x(ii)*ymaxpb*icxpd*icxpd
                    dLambdadu(2) = dLambdadu(2) + f(ii)*ymaxpb*icxpd*icxpd
                    dLambdadu(3) = dLambdadu(3) + f(ii)*x(ii)*(ymaxpb*ymaxpb - cxpd*cxpd)*icxpd*icxpd*icxpd
                    dLambdadu(4) = dLambdadu(4) + f(ii)*(ymaxpb*ymaxpb - cxpd*cxpd)*icxpd*icxpd*icxpd

                    d2Lambdadudu(1,1) = d2Lambdadudu(1,1) - f(ii)*x(ii)*x(ii)*icxpd*icxpd
                    d2Lambdadudu(2,1) = d2Lambdadudu(2,1) - f(ii)*x(ii)*icxpd*icxpd
                    d2Lambdadudu(3,1) = d2Lambdadudu(3,1) - f(ii)*2*x(ii)*x(ii)*ymaxpb*icxpd*icxpd*icxpd
                    d2Lambdadudu(4,1) = d2Lambdadudu(4,1) - f(ii)*2*x(ii)*ymaxpb*icxpd*icxpd*icxpd
                                        
                    d2Lambdadudu(2,2) = d2Lambdadudu(2,2) - f(ii)*icxpd*icxpd
                    d2Lambdadudu(3,2) = d2Lambdadudu(3,2) - f(ii)*2*x(ii)*ymaxpb*icxpd*icxpd*icxpd
                    d2Lambdadudu(4,2) = d2Lambdadudu(4,2) - f(ii)*2*ymaxpb*icxpd*icxpd*icxpd

                    d2Lambdadudu(3,3) = d2Lambdadudu(3,3) - f(ii)*x(ii)*x(ii)*(3*ymaxpb*ymaxpb - cxpd*cxpd)*icxpd*icxpd*icxpd*icxpd
                    d2Lambdadudu(4,3) = d2Lambdadudu(4,3) - f(ii)*x(ii)*(3*ymaxpb*ymaxpb - cxpd*cxpd)*icxpd*icxpd*icxpd*icxpd

                    d2Lambdadudu(4,4) = d2Lambdadudu(4,4) - f(ii)*(3*ymaxpb*ymaxpb - cxpd*cxpd)*icxpd*icxpd*icxpd*icxpd
                end do

                return
            end subroutine computeDerivs

 
        end subroutine fit

        
        pure real(kind=real64) function computeLambda( x,y , f , a,b,c,d , sig_min  )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute the log likelihood of the input data given the equation of best fit and the non-uniform error bars
    !*          y = a x + b  +/- eps(x)
    !*      with eps(x) = c x + d >= sig_min
            real(kind=real64),dimension(:),intent(in)               ::      x
            real(kind=real64),dimension(:),intent(in)               ::      y
            real(kind=real64),dimension(:),intent(in)               ::      f       !   logistic function goodness of fit = 1 for on line, =0 for not on line
            real(Kind=real64),intent(in)                            ::      a,b     !   equation of best fit line
            real(Kind=real64),intent(in)                            ::      c,d     !   equation of error bars
            real(Kind=real64),intent(in)                            ::      sig_min     !   minimum error bars
  
            integer                 ::      nPoints
            integer                 ::      ii
        
            real(Kind=real64)       ::      axpb,cxpd,ymaxpd,icxpd
            nPoints = size(x,dim=1)
            computeLambda = 0
            do ii = 1,nPoints
                axpb = a*x(ii) + b
                cxpd = max( sig_min , c*x(ii) + d )   !   ensure expected error c*x + d is positive, even if small.
                ymaxpd = y(ii) - axpb
                icxpd = 1/cxpd
                computeLambda = computeLambda - f(ii)*( log ( 2*PI*cxpd*cxpd ) + ymaxpd*ymaxpd*icxpd*icxpd )/2
            end do
            return
        end function computeLambda               

        pure real(kind=real64) function computeRss( x,y , f , a,b  )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),dimension(:),intent(in)               ::      x
            real(kind=real64),dimension(:),intent(in)               ::      y
            real(kind=real64),dimension(:),intent(in)               ::      f       !   logistic function goodness of fit = 1 for on line, =0 for not on line
            real(Kind=real64),intent(in)                            ::      a,b     !   equation of best fit line
  
            integer                 ::      nPoints
            integer                 ::      ii
        
            real(Kind=real64)       ::      res
            nPoints = size(x,dim=1)
            computeRss = 0
            do ii = 1,nPoints            
                res = y(ii) - (a*x(ii) + b)
                computeRss = computeRss + f(ii)*res*res
            end do
            return
        end function computeRss               


        subroutine findDodgyPoints( x,y,a,b,c,d, nsig, f, sig_min )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
    !*      given the points {xi,yi} and a guess for the straight line fit
    !*          y = a x + b     with errorbars  sigma = c x + d
    !*      find the points with residuals < multiple of sigma
    !*      give it a goodness fo fit 
    !*          f = 1/( 1 + exp( (res - nsig x sig )/(sig/3) ) )
    !*      eg nsig = 3 gives 1 out to about 2 sig, then quickly turning through 3 sig to zero at 4 sig.           
            real(kind=real64),dimension(:),intent(in)               ::      x
            real(kind=real64),dimension(:),intent(in)               ::      y
            real(kind=real64),intent(in)                            ::      a,b     !   equation of best fit line
            real(kind=real64),intent(in)                            ::      c,d     !   equation of error bars
            real(kind=real64),intent(in)                            ::      nsig
            real(kind=real64),intent(in)                            ::      sig_min !   minimum error
            real(kind=real64),dimension(:),intent(out)              ::      f
            integer                 ::      nPoints
            integer                 ::      ii

        
            real(Kind=real64)       ::      res,sig ,over_flow_check
            nPoints = size(x,dim=1)
            do ii = 1,nPoints
                res = y(ii) - (a*x(ii) + b)
                sig = max(sig_min,c*x(ii) + d)      
                !f(ii) = 1/( 1 + safeExp( 3*(res/sig - nsig) ))
                over_flow_check=3*(res/sig - nsig)
                if (over_flow_check > 50) then
                    f(ii)=0.0d0
                else
                    f(ii) = 1/( 1 + exp(over_flow_check)) 
                end if  
                            
                !f(ii) = 1 - erf( res*0.707106781d0/(nsig*sig) )   
            end do

            return
        end subroutine findDodgyPoints


    end module Lib_Heteroscedastic_Fit