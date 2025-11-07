    program testLib_DrawEllipse
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*-----------------------------------------------------------------------------------------------------------------------------------
!*    testLib_DrawEllipse from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
!*    Copyright (C) 2024  James Heath

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
        use Lib_Colourscale
        use Lib_ColouredTerminal !needed
        use Lib_DrawEllipse
        use Lib_Png
        use Lib_UtilsForTests
        use iso_fortran_env
        implicit none

        real(kind=real64),parameter     ::      PI = 3.14159265390d0
        real(kind=real64),parameter     ::      DEGINRAD = PI/180.0d0
        real(kind=real64)       ::      s1 = 20, s2 = 10 , theta = 30*DEGINRAD

        real(kind=real64),dimension(2,2)        ::      D0,D1
        real(kind=real64),dimension(2)          ::      p0,p1

        real(kind=real64),dimension(:,:),allocatable            ::      img_grey
        real(kind=real64),dimension(:,:,:),allocatable          ::      img_rgb

        real(kind=real64)       ::      iou,f0,f1,expected,calculated,rssp,rssm
        real(kind=real64),dimension(6)          ::      ellipse,drss,spot_dat1,spot_dat2
        real(kind=real64),dimension(6,4)        ::      test_spot_list

        integer             ::      Nx 
        logical             ::      ok,overlap_check
        integer             ::      correctcount
        integer             ::      noofTests=15 !
        character(len=32)   ::      libName="testLib_DrawEllipse"
        logical             ::      tempCheck
        real(kind=real64)   ::      tolerance=1.0e-8 !floating ponit tolerance
        integer             ::      ii

        correctcount=0

        Nx = floor(10*max(s1,s2)) 

        D0 = getDmatrix( s1,s2,theta )
        p0 = Nx*0.5;

        D1 = getDmatrix( s2,s1,-theta )
        p1 = p0 + (/ s1,s2 /)
    
        print *,"Regression tests"

        !test 1, can we draw a greyscale ellipse?
        tempCheck=.false.
        allocate(img_grey(0:Nx-1,0:Nx-1))
        img_grey = 0
        call drawEllipse(D0,p0,img_grey,mode=LIB_DRAWELLIPSE_SHADE_RING)
        call drawEllipse(D1,p1,img_grey,mode=LIB_DRAWELLIPSE_SHADE_RING)
        
        call writePng("test_grey.png",img_grey)
        print*, "img grey sum is:", sum(img_grey)
        tempCheck = (abs(sum(img_grey)/295.40029830691043d0-1)<tolerance) 

        call announceSubTest(libName,"drawEllipse: LIB_DRAWELLIPSE_SHADE_RING",1,noofTests,tempCheck,correctcount)

        !test 2, can we draw an RGB 2D gaussian?
        allocate(img_rgb(3,0:Nx-1,0:Nx-1))
        img_rgb = 0
        call drawEllipse(D0,p0,img_rgb,COLOURSCALE_VIRIDIS,f=0.2d0,mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN)
        call drawEllipse(D1,p1,img_rgb,COLOURSCALE_VIRIDIS,f=0.8d0,mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN)
        call write_rgb_png("test_rgb.png",img_rgb)
        print*, "img rgb sum is:", sum(img_rgb)
        tempCheck = (abs(sum(img_rgb)/3249.1868163037811d0-1)<tolerance) 
        call announceSubTest(libName,"drawEllipse: LIB_DRAWELLIPSE_SHADE_GAUSSIAN",2,noofTests,tempCheck,correctcount)

        !test 3, Is the rectangular IoU of identical ellipses 1.0?
        iou = intersectionOverUnion(D0,p0,D0,p0,mode = LIB_DRAWELLIPSE_IOU_RECTANGLE)
        print *,"iou(rectangles) D0,D0 ",iou
        tempCheck = (abs(iou-1.0d0)<tolerance)  

        call announceSubTest(libName,"intersectionOverUnion: LIB_DRAWELLIPSE_IOU_RECTANGLE same",3,noofTests,tempCheck,correctcount)

        !test 4, Does the rectangular IoU of different ellipses equal a known value?
        iou = intersectionOverUnion(D0,p0,D1,p1,mode = LIB_DRAWELLIPSE_IOU_RECTANGLE)
        print *,"iou(rectangles) D0,D1 ",iou
        tempCheck = (abs(iou/0.16347731000546747d0-1)<tolerance) 
        call announceSubTest(libName,"intersectionOverUnion: LIB_DRAWELLIPSE_IOU_RECTANGLE different",4,noofTests,tempCheck,correctcount)
 
        !test 5, Does the circular IoU of identical ellipses equal 1.0? 
        iou = intersectionOverUnion(D0,p0,D0,p0,mode = LIB_DRAWELLIPSE_IOU_CIRCLE)
        print *,"iou(circles) D0,D0 ",iou
        tempCheck = (abs(iou-1.0d0)<tolerance) 
        call announceSubTest(libName,"intersectionOverUnion: LIB_DRAWELLIPSE_IOU_CIRCLE same",5,noofTests,tempCheck,correctcount)
        
        !test 6, Does the circualr IoU of different ellipses equal a known value?
        iou = intersectionOverUnion(D0,p0,D1,p1,mode = LIB_DRAWELLIPSE_IOU_CIRCLE)
        print *,"iou(circles) D0,D1 ",iou
        tempCheck = (abs(iou/0.19565460720676783d0-1)<tolerance) 
        call announceSubTest(libName,"intersectionOverUnion: LIB_DRAWELLIPSE_IOU_CIRCLE different",6,noofTests,tempCheck,correctcount)

 
        !test 7, does the numerical ellipse overlap IoU for indentical ellipses equal 1.0?
        iou = intersectionOverUnion(D0,p0,D0,p0,mode = LIB_DRAWELLIPSE_IOU_ELLIPSE)
        print *,"iou(ellipse) D0,D0 ",iou
        tempCheck = (abs(iou-1.0d0)<tolerance)  
        call announceSubTest(libName,"intersectionOverUnion: LIB_DRAWELLIPSE_IOU_ELLIPSE same",7,noofTests,tempCheck,correctcount)

        !test 8, does the numerical ellipse overlap IoU for differnet ellipses equal a known value?
        
        iou = intersectionOverUnion(D0,p0,D1,p1,mode = LIB_DRAWELLIPSE_IOU_ELLIPSE)
        print *,"iou(ellipse) D0,D1 ",iou
        tempCheck = (abs(iou/0.15290269828291086d0-1)<tolerance)  
        call announceSubTest(libName,"intersectionOverUnion: LIB_DRAWELLIPSE_IOU_ELLIPSE different",8,noofTests,tempCheck,correctcount)

        !test 9, can we check that two ellipses that are on top of each other overlap with IoU
        overlap_check=.false.
        spot_dat1=(/1.0d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0/)!x,y,dmaj,f0,dmin,theta
        spot_dat2=(/1.0d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0/)
        call do_spots_overlap(spot_dat1,spot_dat2,0.5d0,1.0d0,overlap_check)
        tempCheck=overlap_check
        call announceSubTest(libName,"do_spots_overlap: identical spots",9,noofTests,tempCheck,correctcount)

        !test 10, can we check that two ellipses that are 0.25 over the IoU threshold overlap
        overlap_check=.false.
        spot_dat1=(/1.0d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0/)!x,y,dmaj,f0,dmin,theta
        spot_dat2=(/0.5d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0/)
        call do_spots_overlap(spot_dat1,spot_dat2,0.5d0,1.0d0,overlap_check)
        tempCheck=overlap_check
        call announceSubTest(libName,"do_spots_overlap: 0.75 IoU overlap",10,noofTests,tempCheck,correctcount)

        !test 11, can we check that two ellipses that are exactly over the IoU threshold overlap.
        overlap_check=.false.
        spot_dat1=(/1.0d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0/)!x,y,dmaj,f0,dmin,theta
        spot_dat2=(/0.0d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0/)
        call do_spots_overlap(spot_dat1,spot_dat2,0.5d0,1.0d0,overlap_check)
        tempCheck=overlap_check
        call announceSubTest(libName,"do_spots_overlap: 0.5 IoU overlap",11,noofTests,tempCheck,correctcount)

        !test 12 ,can  we check that two ellipses that only just fail the IoU threshold do not overlap
        overlap_check=.true.
        spot_dat1=(/1.1d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0/)!x,y,dmaj,f0,dmin,theta
        spot_dat2=(/0.0d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0/)
        call do_spots_overlap(spot_dat1,spot_dat2,0.5d0,1.0d0,overlap_check)
        tempCheck=.not.overlap_check
        call announceSubTest(libName,"do_spots_overlap: <0.5 IoU overlap",12,noofTests,tempCheck,correctcount)

        !test 13, can we check that two ellipses that are 0.25 under the IoU threshold do not overlap
        overlap_check=.true.
        spot_dat1=(/1.5d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0/)!x,y,dmaj,f0,dmin,theta
        spot_dat2=(/0.0d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0/)
        call do_spots_overlap(spot_dat1,spot_dat2,0.5d0,1.0d0,overlap_check)
        tempCheck=.not.overlap_check
        call announceSubTest(libName,"do_spots_overlap: 0.25 IoU overlap",13,noofTests,tempCheck,correctcount)

        !test 14, can we check that two ellipses that are completely isolated do not overlap.
        overlap_check=.true.
        spot_dat1=(/3.0d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0/)!x,y,dmaj,f0,dmin,theta
        spot_dat2=(/0.0d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0/)
        call do_spots_overlap(spot_dat1,spot_dat2,0.5d0,1.0d0,overlap_check)
        tempCheck=.not.overlap_check
        call announceSubTest(libName,"do_spots_overlap: 0.0 IoU overlap",14,noofTests,tempCheck,correctcount)

        !test 15, can we check if a spot overlaps with an array of existing spots
        overlap_check=.false.
        spot_dat1=(/1.0d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0/)!x,y,dmaj,f0,dmin,theta
        test_spot_list=reshape((/5.0d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0, &
                                10.0d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0, &
                                0.5d0,1.0d0,2.0d0,1.0d0,1.0d0,0.0d0, &!!!!!!
                                1.0d0,10.0d0,2.0d0,1.0d0,1.0d0,0.0d0/),(/6,4/))

        call check_overlap_with_prexisting_spots(spot_dat1,test_spot_list,0.5d0,1.0d0,overlap_check)
        tempCheck=overlap_check
        call announceSubTest(libName,"check_overlap_with_prexisting_spots: list spots, 1 overlap",15,noofTests,tempCheck,correctcount)



    !---
        ok=haveAllSubTestsPassed(correctcount,noofTests)


    !---    unit tests DRM March 25


        print *,""
        print *,"unit tests for LIB_DRAWELLIPSE_SHADE_GAUSSIAN"

    !   Compute a gaussian with known intensity, position, orientation                 
        D0 = getDmatrix( s1,s2,theta )
        p0 = Nx*0.5d0
        f0 = 0.8d0
        
         

    !---    unit test 1: correct integral
    !   the integral should be
    !       int f0 Exp[ - (x-p).D(x-p) ] dA = f0 pi / sqrt( Dxx Dyy - Dxy^2 )
        img_grey = 0
        call drawEllipse(D0,p0,img_grey,mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN,f=f0)
        expected = PI * f0 / sqrt( D0(1,1)*D0(2,2)-D0(1,2)*D0(2,1) )
        calculated = sum(img_grey)
        tempCheck = (abs(calculated/expected-1)<0.001d0)               !   0.1% tolerance - we are comparing pixels in a finite range with analytic to infinity. Expect answer to be low, so no abs
        print *," unit test 1. gaussian integral calc ",calculated," expected ",expected,tempCheck
        ok = ok .and. tempCheck


    !---    unit test 1a: correct integral
    !   the integral should be
    !       int f0 Exp[ - (x-p).D(x-p) ] dA = f0 pi / sqrt( Dxx Dyy - Dxy^2 )        
        ellipse = (/ p0(1),p0(2),f0,D0(1,1),D0(1,2),D0(2,2) /)
        expected = PI * f0 / sqrt( D0(1,1)*D0(2,2)-D0(1,2)*D0(2,1) )
        calculated = getWeight( ellipse,mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN )
        tempCheck = (abs(calculated/expected-1)<0.001d0)               !   0.1% tolerance - we are comparing pixels in a finite range with analytic to infinity. Expect answer to be low, so no abs
        print *," unit test 1a. weight ",calculated," expected ",expected,tempCheck
        ok = ok .and. tempCheck

    !---    unit test 1b: correct integral
    !   the integral should be
    !       int f0 Exp[ - (x-p).D(x-p) ] dA = f0 pi / sqrt( Dxx Dyy - Dxy^2 )
        img_grey = 0
        call add( ellipse,img_grey,mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN )
        expected = PI * f0 / sqrt( D0(1,1)*D0(2,2)-D0(1,2)*D0(2,1) )
        calculated = sum(img_grey)
        tempCheck = (abs(calculated/expected-1)<0.001d0)               !   0.1% tolerance - we are comparing pixels in a finite range with analytic to infinity. Expect answer to be low, so no abs
        print *," unit test 1b. add ",calculated," expected ",expected,tempCheck
        ok = ok .and. tempCheck



    !---    unit test 2: correct rss easy mode - input = output. Should give zero.
        
        calculated = getRss( ellipse,img_grey,mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN )
        expected = 0.0d0
        tempCheck = (calculated<tolerance)                          !      tiny tolerance - basically checking get same answer twice.
        print *," unit test 2. rss self-test calc ",calculated," expected ",expected,tempCheck
        ok = ok .and. tempCheck


        
    !---    unit test 3: correct rss easyish mode - input = multiple of output
    !   should have int ( f0 Exp[ - (x-p).D(x-p) ] - f1 Exp[ - (x-p).D(x-p) ] )^2 dA = (f0-f1)2 pi / ( 2 sqrt( Dxx Dyy - Dxy^2 ) )
        f1 = 1.0d0     
        ellipse = (/ p0(1),p0(2),f1,D0(1,1),D0(1,2),D0(2,2) /)
        calculated = getRss( ellipse,img_grey,mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN )
        expected = (PI/2) * ((f0-f1)**2) / sqrt( D0(1,1)*D0(2,2)-D0(1,2)*D0(2,1) )
        tempCheck = (abs(calculated/expected - 1)<0.001d0)                         !     0.1% tolerance  
        print *," unit test 3. rss different intensity integral calc ",calculated," expected ",expected,tempCheck
        ok = ok .and. tempCheck


    !---    unit test 4: correct rss easyish mode - input = output, but now cut out half the input> should give zero.
        ellipse = (/ p0(1),p0(2),f0,D0(1,1),D0(1,2),D0(2,2) /)
        img_grey(Nx/2:,:) = LIB_DRAWELLIPSE_IGNORE
        calculated = getRss( ellipse,img_grey,mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN )
        expected = 0.0d0
        tempCheck = (calculated<tolerance)                              !      tiny tolerance - basically checking can ignore pixels
        print *," unit test 4. rss half image blanked ",calculated," expected ",expected,tempCheck
        ok = ok .and. tempCheck


    !---    unit test 5. derivatives. Draw image again
        call drawEllipse(D0,p0,img_grey,mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN,f=f0)
        call findDrss( ellipse,img_grey,mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN,rss=calculated,drss=drss)
    !   compute derivative numerically
        do ii = 1,6         !   ellipse parameters
            f1 = ellipse(ii)
            ellipse(ii) = f1 * 1.01d0
            rssp = getRss( ellipse,img_grey,mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN )
            ellipse(ii) = f1 * 0.99d0
            rssm = getRss( ellipse,img_grey,mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN )
            ellipse(ii) = f1
            calculated = drss(ii)
            expected = (rssp - rssm) / (f1*0.02d0)
            tempCheck = (max(abs(calculated),abs(expected))<1.0d-6).or.(abs(calculated/expected - 1)<0.01d0)                         !     1% tolerance. Or accept if there is derivative = 0
            ok = ok .and. tempCheck
            print *," unit test 5. drss(",ii,") ",calculated," expected ",expected,tempCheck
        end do



    !---    unit test 6. integral for circle
    !   t = int f0 Exp[ - (x-p).D(x-p) ] dA  )   where n is number of pixels within 1 std dev
    !     = ( f0 pi ( exp[-0.5] - 1 ) / Dxx )  
        s1 = 20.0d0
        D1 = getDmatrix( s1,s1,0.0d0 )
        f1 = 1.0d0
        ellipse = (/ p0(1),p0(2),f1,D1(1,1),D1(1,2),D1(2,2) /)

        
        calculated = getWeight( ellipse,mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN,mult=2.0d0 )
        expected = PI * f1 *( 1.0d0 - exp(-2.0d0) ) / D1(1,1)   !   integral inside 2 sigma
        tempCheck = (abs(calculated/expected - 1)<0.01d0)                         !     huge 1% tolerance, beacuse I've got a gauss circle issue here.
        print *," unit test 6. integral to 2 std dev ",calculated," expected ",expected,tempCheck
        ok = ok .and. tempCheck


    !---    unit test 7.T-value for circle
    !   t = int f0 Exp[ - (x-p).D(x-p) ] dA / (sqrt(n)*sigma)   where n is number of pixels within 1 std dev
    !     = ( f0 pi ( exp[-0.5] - 1 ) / Dxx ) / sqrt(Tr(D)^2/4 - det(D)) / sigma        
        
        
        calculated = getT( ellipse,mode=LIB_DRAWELLIPSE_SHADE_GAUSSIAN,sigma=1.0d0 )
        expected = PI * f1 *( 1.0d0 - exp(-0.5d0) ) / D1(1,1)   !   integral inside 1 sigma
        expected = expected / sqrt( PI/(2*D1(1,1)) )            !   divide by sqrt area inside 1 sigma
        tempCheck = (abs(calculated/expected - 1)<0.01d0)                         !     huge 1% tolerance, beacuse I've got a gauss circle issue here.
        print *," unit test 7. T value ",calculated," expected ",expected,tempCheck
        ok = ok .and. tempCheck

        print *,""
        







        !-----------------------------------------------
        

        call announcePassOrFail(ok)

    end program testLib_DrawEllipse