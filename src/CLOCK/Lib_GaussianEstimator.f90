!*-----------------------------------------------------------------------------------------------------------------------------------
!*    Lib_GaussianEstimator from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
!*    Copyright (c) UKAEA Jan 2025  Daniel Mason

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
!*      A simple module to make an estimate of a gaussian function N(mu,sig)
!*      in 1d, given that we only have information about the gaussian over range (0:x) and not (-infty:infty)

    module Lib_GaussianEstimator
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      given data as a histogram 
        !*          h(i) = counts within range i dx:(i+1) dx
        !*      where 0<=i<N and N dx = x
        !*      find an estimator for the best fit single gaussian function to return these counts.
        !*
        !*      Version History
        !*          v.0.0.0     Jan 25      First working version
        !*
        
                use iso_fortran_env
                implicit none
                private
        
            !---    default values for tables
                integer,parameter,private           ::      NDIV0 = 2000
                real(kind=real64),parameter,private ::      X0 = 20.0d0
                character(len=256),private          ::      FILENAME0 = "GaussTables.dat"
        
            !---    current values of tables
                integer,private             ::      nDiv = 0
                real(kind=real64),private   ::      sigmax = X0
                real(kind=real64),private   ::      mumax = X0        
                real(kind=real64),dimension(:,:),private,allocatable      ::      fbar,f2bar
        
        
                public      ::      findMeanAndMeanSquare
                public      ::      integral
                public      ::      gaussianMeanAndMeanSquare
                public      ::      GaussianEstimate
                
                public      ::      writeTables
                public      ::      readTables
        
                interface findMeanAndMeanSquare
                    module procedure    findMeanAndMeanSquare0
                    module procedure    findMeanAndMeanSquare1
                end interface    
        
                interface integral
                    module procedure    integral0
                end interface    
                
        
            contains
        !---^^^^^^^^
        
                pure real(kind=real64) function integral0( mu,sig,x )
            !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !*      find the integral of the gaussian from 0:x
            !*      note that it is assumed that the gaussian is noramlised to integrate to 1 over real line.
                    real(kind=real64),intent(in)        ::      mu,sig          !   gaussian properties
                    real(kind=real64),intent(in)        ::      x               !   xmax
        
                    real(kind=real64)       ::      irt2s
                    if (sig < 1.0d-8*x) then
                        if (mu == 0) then
                            integral0 = 0.5d0
                        else if (mu == x) then
                            integral0 = 0.5d0
                        else if ( (mu<0).or.(mu>x) ) then
                            integral0 = 0.0d0
                        else
                            integral0 = 1.0d0
                        end if
                    else
                        if (mu + 10*sig < 0) then
                            integral0 = 0.0d0
                        else if (mu - 10*sig > x) then
                            integral0 = 0.0d0
                        else if ( (mu-10*sig > 0).and.(mu+10*sig < x) ) then
                            integral0 = 1.0d0
                        else
                            irt2s = 0.7071067811865d0/sig
                            integral0 = ( erf( (x-mu)*irt2s ) + erf( mu*irt2s ) )/2
                        end if
                    end if
                    return
                end function integral0
        
                pure subroutine findMeanAndMeanSquare0( h,x, hbar,h2bar )
            !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !*      given the histogram defined in range 0:x
            !*      with equally spaced bins
            !*      compute the mean and mean square
            !*               y2
            !*            _/|
            !*          _/  |
            !*      y1_/    |           f(x) = (y2-y1)/(x2-x1) x + (x2y1 - x1y2)/(x2-x1)
            !*        |     |       
            !*        |_____| 
            !*      x1      x2 = x1+dx
            !*                         
            !*              int f(x) dx     = (y1 + y2) dx/2
            !*              int x f(x) dx   = (x1 (y1+y2) dx)/2 + (y1+2y2) dx^2/6
            !*              int x^2 f(x) dx = (x1^2 (y1+y2) dx)/2 + x1 (y1+2y2) dx^2/3 + (y1+3y2) dx^3/12
            !*
        
                    integer,dimension(0:),intent(in)            ::      h
                    real(kind=real64),intent(in)                ::      x
                    real(kind=real64),intent(out)               ::      hbar,h2bar
                    integer                                     ::      N
        
                    integer             ::      ii,sumh
                    real(kind=real64)   ::      dx,xx
                    
        
                !---    set up problem
                    N = size(h)
                    hbar = 0
                    h2bar = 0
        
                    if (N == 0) then
                        return
                    else 
                        sumh = sum(h)
                        if (sumh == 0) return
                    end if
        
                     
                    dx = x/N
        
                    do ii = 0,N-1
                        xx = (ii+0.5d0)*dx
                        hbar = hbar + xx*h(ii)
                        h2bar = h2bar + xx*xx*h(ii)
                    end do
                    xx = 1.0d0/sumh
                    hbar = hbar*xx
                    h2bar = h2bar*xx
                    
         
                    return
                end subroutine findMeanAndMeanSquare0
        
                pure subroutine findMeanAndMeanSquare1( h,x, hbar,h2bar )
            !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !*      given the histogram defined in range 0:x
            !*      with equally spaced bins
            !*      compute the mean and mean square
            !*               y2
            !*            _/|
            !*          _/  |
            !*      y1_/    |           f(x) = (y2-y1)/(x2-x1) x + (x2y1 - x1y2)/(x2-x1)
            !*        |     |       
            !*        |_____| 
            !*      x1      x2 = x1+dx
            !*                         
            !*              int f(x) dx     = (y1 + y2) dx/2
            !*              int x f(x) dx   = (x1 (y1+y2) dx)/2 + (y1+2y2) dx^2/6
            !*              int x^2 f(x) dx = (x1^2 (y1+y2) dx)/2 + x1 (y1+2y2) dx^2/3 + (y1+3y2) dx^3/12
            !*
        
                    real(kind=real64),dimension(0:),intent(in)  ::      h
                    real(kind=real64),intent(in)                ::      x
                    real(kind=real64),intent(out)               ::      hbar,h2bar
                    integer                                     ::      N
        
                    integer             ::      ii
                    real(kind=real64)   ::      dx,xx,sumh
                    
        
                !---    set up problem
                    N = size(h)
                    hbar = 0
                    h2bar = 0
        
                    if (N == 0) then
                        return
                    else 
                        sumh = sum(h)
                        if (sumh <= 0) return
                    end if
        
                     
                    dx = x/N
        
                    do ii = 0,N-1
                        xx = (ii+0.5d0)*dx
                        hbar = hbar + xx*h(ii)
                        h2bar = h2bar + xx*xx*h(ii)
                    end do
                    xx = 1.0d0/sumh
                    hbar = hbar*xx
                    h2bar = h2bar*xx
                    
         
                    return
                end subroutine findMeanAndMeanSquare1
        
                
                subroutine gaussianMeanAndMeanSquare( mu,sig,x , fbar,f2bar )
            !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !*      given a normal distribution N(mu,sig)
            !*      over the interval 0:x with x>0
            !*      compute the mean and mean square
                    real(kind=real64),intent(in)                ::      mu,sig
                    real(kind=real64),intent(in)                ::      x
                    real(kind=real64),intent(out)               ::      fbar,f2bar
        
                    real(kind=real64),parameter                 ::      ISQRT2 = 0.707106781186548d0
                    real(kind=real64),parameter                 ::      SQRTTWOONPI = 0.797884560802865d0
                    real(kind=real64)               ::      aa,bb,dd,isig
        
                    !print *," gaussianMeanAndMeanSquare( mu,sig,x )",mu,sig,x
                    if (abs(sig)<x*1.0d-8) then
                        if ( (mu>=0).and.(mu<x) ) then
                            fbar = mu
                            f2bar = mu*mu
                        else
                            fbar = 0.0d0
                            f2bar = 0.0d0
                        end if
                        return                
                    end if
        
                    if (mu + 5*sig < 0) then
                        fbar = 0.0d0
                        f2bar = 0.0d0                
                        return                
                    else if (mu - 5*sig > x) then
                        fbar = 0.0d0
                        f2bar = 0.0d0                
                        return                
                    end if
        
                    isig = 1/abs(sig)
                    if (x < sig/5) then
                        !   when x is very small, we can take a Taylor series expansion
                        fbar = ( 1 + (x*isig*isig/6)*(mu - x/2) )*x/2
                        f2bar = ( 1 + mu*isig*isig/4 )*x*x/3
                        return
                    end if
        
                    aa = x*ISQRT2*isig
                    bb = mu*ISQRT2*isig
        
                    dd = 1/( erf( aa-bb ) + erf( bb ) )     !   note: x>0, so erf(aa-bb) > erf(-bb)
        
                    fbar = exp( -bb*bb ) - exp( -(aa-bb)*(aa-bb) )
                    fbar = mu + SQRTTWOONPI*sig*dd*fbar
         
        
                    f2bar = mu*exp( - bb*bb ) - (x + mu)*exp( -(aa-bb)*(aa-bb) ) 
                    f2bar = mu*mu + sig*sig + sig*SQRTTWOONPI*dd*f2bar
                    
                    if ((fbar /= fbar).or. (f2bar /= f2bar)) then
                        print *,"mu,sig,aa,bb,dd ",mu,sig,aa,bb,dd
                        print *,"fbar,f2bar ",fbar,f2bar
                        stop
                    end if
                 
        
                    
                    return
                end subroutine gaussianMeanAndMeanSquare
        
        
                subroutine writeTables( sigmax_in,mumax_in,nDiv_in , filename_in)
            !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !*      compute tables of fbar and f2bar 
            !*      for gaussian between 0 and 1
            !*      0<sig<sigmax
            !*      -mumax<mu<mumax
                    real(kind=real64),intent(in),optional   ::      sigmax_in,mumax_in
                    integer,intent(in),optional             ::      nDiv_in
                    character(len=*),intent(in),optional    ::      filename_in
                    character(len=256)      ::      filename
                    real(kind=real64)       ::      mu,sig,dmu,dsig
                    
                    integer                 ::      ii,jj
        
                    sigmax = X0
                    mumax = X0
                    nDiv = NDIV0
                    filename = FILENAME0
        
                    if (present(sigmax_in)) sigmax = sigmax_in
                    if (present(mumax_in)) mumax = mumax_in
                    if (present(nDiv_in)) nDiv = nDiv_in
                    if (present(filename_in)) filename = filename_in
        
                    dsig = sigmax/nDiv
                    dmu = mumax/nDiv
        
                    if (allocated(fbar)) then
                        deallocate(fbar)
                        deallocate(f2bar)
                    end if
                    allocate(fbar(-nDiv:nDiv,0:nDiv))
                    allocate(f2bar(-nDiv:nDiv,0:nDiv))
        
        
                    open(file=trim(filename),unit=500,action="write",form="unformatted")
                        write(unit=500) nDiv,sigmax,mumax
                        do ii = 0,nDiv
                            sig = dsig*ii
                            do jj = -nDiv,nDiv
                                mu = dmu*jj
                                call gaussianMeanAndMeanSquare( mu,sig,1.0d0 , fbar(jj,ii),f2bar(jj,ii) )
                                write(unit=500) fbar(jj,ii),f2bar(jj,ii)
                            end do
                        end do
                    close(unit=500)
        
                    return
                end subroutine writeTables
        
        
        
        
                subroutine readTables( filename_in )
            !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !*      read tables of fbar and f2bar 
            !*      for gaussian between 0 and 1
            !*      0<sig<sigmax
            !*      -mumax<mu<mumax
                    character(len=256),intent(in),optional           ::      filename_in
        
                    integer                 ::      ii,jj,ioerr
                    logical                 ::      ok
                    character(len=256)      ::      filename
        
                    filename = FILENAME0
                    if (present(filename_in)) filename = filename_in
                    if (allocated(fbar)) then
                        deallocate(fbar)
                        deallocate(f2bar)
                    end if
        
                    inquire(file=trim(filename),exist=ok)
                    if (ok) then
                        open(file=trim(filename),unit=600,action="read",form="unformatted")
                            read(unit=600,iostat=ioerr) nDiv,sigmax,mumax
                            ok = (ioerr==0) 
                            if (ok) then
                                allocate(fbar(-nDiv:nDiv,0:nDiv))
                                allocate(f2bar(-nDiv:nDiv,0:nDiv))
                    
                                do ii = 0,nDiv
                                    do jj = -nDiv,nDiv
                                        read(unit=600,iostat=ioerr) fbar(jj,ii),f2bar(jj,ii)
                                        ok = ok .and. (ioerr==0)
                                    end do
                                end do
                            end if
                        close(unit=600)
                    end if
        
                    if (.not. ok) then
                        print *,"Lib_GaussianEstimator::readTables() error - tables not found. Recomputing ..."
                        call writeTables( X0,X0,NDIV0 , filename)
                        print *,"Lib_GaussianEstimator::readTables() info - tables written to """//trim(filename)//""""
                    end if
        
                    return
                end subroutine readTables
        
                subroutine GaussianEstimate( fbar_in,f2bar_in,x, mubest,sigbest )
            !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !*      given the moments fbar_in,f2bar_in and the range 0:x
            !*      find a good estimate for the original gaussian mu,sig
                    real(kind=real64),intent(in)            ::          fbar_in,f2bar_in
                    real(kind=real64),intent(in)            ::          x
                    real(kind=real64),intent(out)           ::          mubest,sigbest
        
                    real(kind=real64)       ::      fbaronx,f2baronx2,dmu,dsig
                    integer                 ::      ii,jj 
                    real(kind=real64)       ::      ss,sa,sb,sbest , mu,sig
                    real(kind=real64)       ::      xa,ya,xb,yb , gg,hh
                    integer                 ::      indx
        
                    fbaronx = fbar_in/x
                    f2baronx2 = f2bar_in/(x*x)
        
                !---    do we have the tables?
                    if (nDiv==0) call readTables()
                    dmu = mumax/nDiv
                    dsig = sigmax/nDiv
                    
        
                !---    find the point in the table which crosses the answer, or the best ( least square )
                    mubest = fbaronx ; sigbest = sqrt(f2baronx2-fbaronx*fbaronx)  
                    call gaussianMeanAndMeanSquare( mubest,sigbest,1.0d0 , gg,hh ) 
                    sbest = ( gg - fbaronx )**2 + ( hh - f2baronx2 )**2 
                    do ii = 0,nDiv
                        sig = ii*dsig
        
                        do jj = -nDiv,nDiv
                            mu = jj*dmu 
                            ss = ( fbar(jj,ii) - fbaronx )**2 + ( f2bar(jj,ii) - f2baronx2 )**2 
                           ! print *,"test ",mu,sig,ss,sbest
                            indx = 0                !   bit pattern for whether we are above/below solution
                            if ( (ii<nDiv).and.(jj<nDiv)) then
                                if ( fbar(jj  ,ii  ) >fbaronx   ) indx = indx + 1
                                if ( fbar(jj+1,ii  ) >fbaronx   ) indx = indx + 2
                                if ( fbar(jj  ,ii+1) >fbaronx   ) indx = indx + 4
                                if ( fbar(jj+1,ii+1) >fbaronx   ) indx = indx + 8
                                if ( f2bar(jj  ,ii  )>f2baronx2 ) indx = indx + 16
                                if ( f2bar(jj+1,ii  )>f2baronx2 ) indx = indx + 32
                                if ( f2bar(jj  ,ii+1)>f2baronx2 ) indx = indx + 64
                                if ( f2bar(jj+1,ii+1)>f2baronx2 ) indx = indx + 128
                            end if
                            
        
                            if ( (iand(indx,15)==0).or.(iand(indx,15)==15).or.(iand(indx,240)==0).or.(iand(indx,240)==240) ) then
                                sa = huge(1.0)
                                sb = huge(1.0)
                            else
                                !   there is at least one high, one low for both fbar and f2bar in this square of 4.
                                !print *,"i,j,mu,sig,indx",ii,jj,mu,sig,indx
                                call bilinearFit( mu,sig,dmu,dsig,                                      &
                                        fbar(jj,ii)-fbaronx,fbar(jj+1,ii)-fbaronx,fbar(jj,ii+1)-fbaronx,fbar(jj+1,ii+1)-fbaronx,                &
                                        f2bar(jj,ii)-f2baronx2,f2bar(jj+1,ii)-f2baronx2,f2bar(jj,ii+1)-f2baronx2,f2bar(jj+1,ii+1)-f2baronx2,    &
                                        xa,ya,xb,yb )
        
                                !print *," xa,ya,xb,yb ",xa,ya,xb,yb
                                call gaussianMeanAndMeanSquare( xa,ya,1.0d0 , gg,hh ) 
                                sa = ( gg - fbaronx )**2 + ( hh - f2baronx2 )**2 
        
                                call gaussianMeanAndMeanSquare( xb,yb,1.0d0 , gg,hh ) 
                                sb = ( gg - fbaronx )**2 + ( hh - f2baronx2 )**2 
                                !print *,sa,sb
                                !if ( (sa/=sa).or.(sb/=sb) ) stop
                            end if
        
                            if (minval( (/ss,sa,sb/) )<sbest) then
                                if (ss < min(sa,sb)) then
                                    sbest = ss
                                    mubest = mu
                                    sigbest = sig
                                else if (sa < min(ss,sb)) then
                                    sbest = sa
                                    mubest = xa
                                    sigbest = ya
                                else 
                                    sbest = sb
                                    mubest = xb
                                    sigbest = yb
                                end if
                            end if
        
                        end do
                    end do
                    mubest = mubest*x
                    sigbest = sigbest*x 
        
                    return
                end subroutine GaussianEstimate
        
        
                subroutine bilinearFit( x1,y1,dx,dy, g1,g2,g3,g4 , h1,h2,h3,h4 , xa,ya,xb,yb )
            !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            !*      given the value of two functions on box corners,
            !*      x1 y1 ----- x2 y2           g1 ----- g2
            !*          |       |                |       !
            !*          |       |                !       !
            !*      x1 y2 ----- x2 y2           g3 ----- g4
            !*      find the points (xa,ya) (xb,yb) where the functions equal 0
                    real(kind=real64),intent(in)            ::      x1,y1,dx,dy
                    real(kind=real64),intent(in)            ::      g1,g2,g3,g4
                    real(kind=real64),intent(in)            ::      h1,h2,h3,h4
                    real(kind=real64),intent(out)           ::      xa,ya
                    real(kind=real64),intent(out)           ::      xb,yb
        
                    real(kind=real64)           ::      d1,d2
                    real(kind=real64)           ::      ss
        
        
                !---    check denominators for actual solution
                    d1 = (2*(-((g3 - g4)*(h1 - h2)) + (g1 - g2)*(h3 - h4)))
                    d2 = (2*(-((g2 - g4)*(h1 - h3)) + (g1 - g3)*(h2 - h4)))
                    if ( abs(d1*d2)<1.0d-16 ) then
                        xa = 0
                        xb = 0
                        ya = 0
                        yb = 0
                        return
                    end if
        
                !---    check square roots for solution
                    ss = (g4**2*h1**2 + g3**2*h2**2 + (g2*h3 - g1*h4)**2 - 2*g4*(g3*h1*h2 + g2*h1*h3 -   &
                        2*g1*h2*h3 + g1*h1*h4) - 2*g3*(g2*h2*h3 - 2*g2*h1*h4 + g1*h2*h4))
                    if ( ss<0 ) then
                        xa = 0
                        xb = 0
                        ya = 0
                        yb = 0
                        return
                    end if
         
                    ! print *,"bilinearFit g ",g1,g2,g3,g4
                    ! print *,"bilinearFit h ",h1,h2,h3,h4
                    ! print *,"denom ",(2*(-((g3 - g4)*(h1 - h2)) + (g1 - g2)*(h3 - h4)))     &
                    !                 ,(2*(-((g2 - g4)*(h1 - h3)) + (g1 - g3)*(h2 - h4)))     &
                    !                 ,(2*(-((g3 - g4)*(h1 - h2)) + (g1 - g2)*(h3 - h4)))     &
                    !                 ,(2*(-((g2 - g4)*(h1 - h3)) + (g1 - g3)*(h2 - h4)))
        
                    xa = -(dx*(2*g3*h1 - g4*h1 - g3*h2 - 2*g1*h3 + g2*h3 + g1*h4 +                          &
                        Sqrt(ss)) -                &
                        2*(-((g3 - g4)*(h1 - h2)) + (g1 - g2)*(h3 - h4))*x1)/                               &
                        d1
                    ya = -(-(dy*(g4*h1 + 2*g1*h2 - g3*h2 + g2*(-2*h1 + h3) - g1*h4 +                        &
                        Sqrt(ss))) -               &
                        2*(-((g2 - g4)*(h1 - h3)) + (g1 - g3)*(h2 - h4))*y1)/                               &
                        d2
                    xb = -(-(dx*(g4*h1 + g3*(-2*h1 + h2) + 2*g1*h3 - g2*h3 - g1*h4 +                        &
                        Sqrt(ss))) -               &
                        2*(-((g3 - g4)*(h1 - h2)) + (g1 - g2)*(h3 - h4))*x1)/                               &
                        d1
                    yb = -(dy*(2*g2*h1 - g4*h1 - 2*g1*h2 + g3*h2 - g2*h3 + g1*h4 +                          &
                        Sqrt(ss)) -                &
                        2*(-((g2 - g4)*(h1 - h3)) + (g1 - g3)*(h2 - h4))*y1)/                               &
                        d2
        
                    return
                end subroutine bilinearFit
        
        
            end module Lib_GaussianEstimator 