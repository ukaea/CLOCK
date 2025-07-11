
    program test_saltandpepperfilter
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*      A simple standalone code to launch flat
!*      with different input options
        use iso_fortran_env
        use Lib_ColouredTerminal
        use NBAX_StringTokenizers
        use Lib_CommandLineArguments
        use Lib_Png
                    use Lib_RandomSeed      !   for dev only
        implicit none

        
    !---    command line parameter input        
        type(CommandLineArguments)      ::      cla
        character(len=8)                ::      test = "sandp1"
        character(len=256)              ::      saltandpepper_exe = "saltAndPepperFilter"!flat_exe = "flat"


    !---    input image to test
        integer,parameter               ::      Nx = 8, Ny = 8
        real(kind=real64),dimension(Nx,Ny),target           ::      flat_grey = 0.5d0
        real(kind=real64),dimension(Nx,Ny),target           ::      flat_black = 0.0d0

        real(kind=real64),dimension(Nx*2,Ny*2),target       ::      noisy_grey = 0.0d0

    !---    images to pass to flat
        real(kind=real64),dimension(:,:),pointer                    ::      img_in
        real(kind=real64),dimension(:,:),allocatable                ::      img_out


    !---    command line to run 
        character(len=256)              ::      infile,outfile
        character(len=256)              ::      flat_cmd

    !---    test
        logical                         ::      ok        
        integer                         ::      ii


    !---    read command line arguments
        cla = CommandLineArguments_ctor(30)  
         
        call setProgramDescription( cla, "test_saltandpepperfilter" )
        call setProgramVersion( cla, "0.0.1" )   
          
        
        call get( cla,"t",test ,LIB_CLA_OPTIONAL," test to run (sandp1|sandp2)" )
        call get( cla,"f",saltandpepper_exe ,LIB_CLA_OPTIONAL," saltandpepper executable" )


        call report(cla)
        if (.not. allRequiredArgumentsSet(cla)) call errorExit("")
        if (hasHelpArgument(cla)) call errorExit("")
        call delete(cla)


        test = convertLowerCase(test)
        infile = "test_saltandpepperfilter."//trim(test)//".png"
        outfile = "test_saltandpepperfilter."//trim(test)//"_out.png"
        select case(trim(test))
            case ("sandp1")
                print *,"salt and pepper test on greyscale image"
                img_in => flat_grey
                flat_cmd = trim(saltandpepper_exe)//" -f "//trim(infile)//" -o "//trim(outfile)
            case ("sandp2")
                print *,"salt and pepper test on black image"
                img_in => flat_black
                flat_cmd = trim(saltandpepper_exe)//" -f "//trim(infile)//" -o "//trim(outfile)
            case ("sandp3")
                print *,"salt and pepper test on noisy grey image"

            !---   how to generate a new noisy image
                !   noisy_grey = reshape( gaussianVariate(size(noisy_grey)), (/size(noisy_grey,dim=1),size(noisy_grey,dim=2)/) )
                !   noisy_grey = 0.50d0 + 0.10d0*noisy_grey

            !---    use an existing noisy image
                noisy_grey = reshape( (/                                                                                                                                                                                                        &
                0.57319367d0, 0.54332448d0, 0.59894715d0, 0.43570930d0, 0.51441636d0, 0.49603114d0, 0.42990108d0, 0.38039276d0, 0.42150639d0, 0.62055195d0, 0.54890541d0, 0.55300602d0, 0.69775283d0, 0.50022892d0, 0.48666911d0, 0.53985581d0, &
                0.48786335d0, 0.53423348d0, 0.31448796d0, 0.56276431d0, 0.40198095d0, 0.54962231d0, 0.65278155d0, 0.43668399d0, 0.53946719d0, 0.51223295d0, 0.47658475d0, 0.57497283d0, 0.52887877d0, 0.36681710d0, 0.60472741d0, 0.45470399d0, &
                0.57112721d0, 0.48997634d0, 0.48070507d0, 0.46754431d0, 0.70395382d0, 0.57144379d0, 0.49387791d0, 0.36015124d0, 0.31878713d0, 0.40317313d0, 0.46527773d0, 0.54603442d0, 0.53075755d0, 0.58068363d0, 0.55430372d0, 0.49458930d0, &
                0.46793598d0, 0.36923011d0, 0.56888399d0, 0.49262617d0, 0.54309214d0, 0.51413543d0, 0.54712390d0, 0.33046700d0, 0.67458610d0, 0.37137507d0, 0.53110969d0, 0.48258209d0, 0.52101581d0, 0.52438574d0, 0.53314480d0, 0.54307609d0, &
                0.43283439d0, 0.58315524d0, 0.57037517d0, 0.42840374d0, 0.57296556d0, 0.49320688d0, 0.45301229d0, 0.40811290d0, 0.53126785d0, 0.57355494d0, 0.48758640d0, 0.70123919d0, 0.57317887d0, 0.55111354d0, 0.49645327d0, 0.38134442d0, &
                0.32397046d0, 0.43977879d0, 0.39252178d0, 0.32694397d0, 0.60301703d0, 0.51479819d0, 0.60951896d0, 0.39162957d0, 0.74364578d0, 0.50819566d0, 0.43206906d0, 0.45919923d0, 0.49675878d0, 0.35315086d0, 0.61262520d0, 0.39322737d0, &
                0.57881931d0, 0.68821410d0, 0.52939816d0, 0.55260909d0, 0.53509543d0, 0.34823733d0, 0.70422817d0, 0.55130615d0, 0.61460932d0, 0.52184115d0, 0.44241957d0, 0.59675744d0, 0.50981756d0, 0.61352569d0, 0.59021081d0, 0.62808629d0, &
                0.30464642d0, 0.45998630d0, 0.48026099d0, 0.56309014d0, 0.47554293d0, 0.49510087d0, 0.54009772d0, 0.53851526d0, 0.54463856d0, 0.43248184d0, 0.36470254d0, 0.47443187d0, 0.62579144d0, 0.45913266d0, 0.53215797d0, 0.42958882d0, &
                0.53965318d0, 0.44709256d0, 0.65146207d0, 0.54369736d0, 0.39649289d0, 0.53343182d0, 0.55846847d0, 0.80412864d0, 0.38571607d0, 0.57184867d0, 0.41207067d0, 0.46717176d0, 0.40714985d0, 0.60813259d0, 0.51543329d0, 0.44506363d0, &
                0.31495836d0, 0.56218390d0, 0.49947897d0, 0.38215978d0, 0.47292464d0, 0.24491331d0, 0.45527633d0, 0.44852286d0, 0.57040283d0, 0.46077525d0, 0.50044555d0, 0.41784802d0, 0.45702272d0, 0.55301217d0, 0.45604606d0, 0.45776533d0, &
                0.49178992d0, 0.55076508d0, 0.60453708d0, 0.35256675d0, 0.47454029d0, 0.47383282d0, 0.46041893d0, 0.67005617d0, 0.72104104d0, 0.63750909d0, 0.62546966d0, 0.43561766d0, 0.44032508d0, 0.68429083d0, 0.59208234d0, 0.53403225d0, &
                0.45658827d0, 0.37221292d0, 0.39688625d0, 0.44628985d0, 0.62070082d0, 0.48251805d0, 0.33408436d0, 0.50533195d0, 0.48355077d0, 0.53821533d0, 0.75835079d0, 0.69768776d0, 0.63168479d0, 0.63243647d0, 0.56713825d0, 0.38821490d0, &
                0.73321105d0, 0.47342826d0, 0.58346298d0, 0.69040755d0, 0.60615132d0, 0.54323669d0, 0.60319033d0, 0.42029692d0, 0.39963568d0, 0.38024893d0, 0.67769620d0, 0.53075107d0, 0.57226698d0, 0.41342846d0, 0.53662641d0, 0.45056158d0, &
                0.55088927d0, 0.37057010d0, 0.61547641d0, 0.52222805d0, 0.52539969d0, 0.62630703d0, 0.54254114d0, 0.63708432d0, 0.53485453d0, 0.46643178d0, 0.42206023d0, 0.46979406d0, 0.40025212d0, 0.29841381d0, 0.57092127d0, 0.52843569d0, &
                0.29304928d0, 0.57622211d0, 0.43821899d0, 0.50095764d0, 0.55820853d0, 0.39798033d0, 0.50941157d0, 0.50014316d0, 0.62930482d0, 0.44652830d0, 0.37846340d0, 0.55697809d0, 0.59046752d0, 0.41198573d0, 0.47806876d0, 0.64297111d0, &
                0.59788925d0, 0.46627967d0, 0.52105307d0, 0.48718376d0, 0.42878609d0, 0.63321384d0, 0.45131118d0, 0.53170626d0, 0.46459507d0, 0.33409122d0, 0.56790562d0, 0.64337340d0, 0.49515342d0, 0.33879201d0, 0.61687386d0, 0.46964333d0  &
                            /),(/Nx*2,Ny*2/) )


                img_in => noisy_grey
                do ii = 1,size(img_in,dim=2)
                    write(*,fmt='(100f16.8)') img_in(:,ii)
                end do
                flat_cmd = trim(saltandpepper_exe)//" -f "//trim(infile)//" -o "//trim(outfile)
            case default
                call errorExit( "error - test """//trim(test)//""" not recognised" )
        end select

        call writePng( infile,img_in )
        call system( trim(flat_cmd) ) 
        call readPng( outfile,img_out )
        
        !call system( "rm "//trim(infile) )
        !call system( "rm "//trim(outfile) )


        select case(trim(test))
            case ("sandp1")
                ok = exact_test( img_in, img_out )
            case ("sandp2")
                ok = exact_test( img_in, img_out )
            case ("sandp3")                 
                ok = dic_test( img_in, img_out, 0.93129463466161500d0 )
            case default
                call errorExit( "error - test """//trim(test)//""" not recognised" )
        end select


        if (ok) then
            print *,colour(GREEN,"PASS")
        else
            print *,colour(RED,"FAIL")
        end if



    contains
!---^^^^^^^^

        logical function exact_test( img_in, img_out ) result(ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if in = out
            real(kind=real64),dimension(:,:),intent(in)         ::      img_in
            real(kind=real64),dimension(:,:),intent(in)         ::      img_out
            
            
            ok = ( size(img_out) == size(img_in) )

            if (.not. ok) then
                print *,"test_flattenBackground::exact_test fail - flat returns different size output image"
            else
                print *,"test_flattenBackground::exact_test info - maxval( abs( img_out - img_in ) ) ",maxval( abs( img_out - img_in ) )
                ok = ( maxval( abs( img_out - img_in ) ) < 1.0/256 )
            end if

            return
        end function exact_test
        
        logical function dic_test( img_in, img_out , dic_target) result(ok)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if in = out
            real(kind=real64),dimension(:,:),intent(in)         ::      img_in
            real(kind=real64),dimension(:,:),intent(in)         ::      img_out
            real(kind=real64),intent(in)                        ::      dic_target
            real(kind=real64)           ::      dd
            
            ok = ( size(img_out) == size(img_in) )

            if (.not. ok) then
                print *,"test_flattenBackground::dic_test fail - flat returns different size output image"
            else
                dd = dic( img_in, img_out )
                print *,"test_flattenBackground::dic_test info - dic( img_in, img_out ) ",dd
                ok = ( abs(dd-dic_target) < 1.0d-4 )
            end if

            return
        end function dic_test
        


        pure real(kind=real64) function dic( f,g )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute dic = sum (f-<f>)(g-<g>) / sqrt( sum(f-<f>)^2 sum(g-<g>)^2 )
            real(kind=real64),dimension(:,:),intent(in)       ::      f,g
            integer             ::      nn,nx,ny,ii,jj
            real(kind=real64)   ::      fsum,gsum,f2sum,g2sum,fgsum
            real(kind=real64)   ::      fbar,gbar
            

            nx = size(f,dim=1)
            ny = size(f,dim=2)
            nn = nx*ny
            fsum = 0
            gsum = 0
            f2sum = 0
            g2sum = 0
            fgsum = 0
            do jj = 1,ny
                do ii = 1,nx
                    fsum    = fsum  + f(ii,jj)
                    f2sum   = f2sum + f(ii,jj)*f(ii,jj)
                    gsum    = gsum  + g(ii,jj)
                    g2sum   = g2sum + g(ii,jj)*g(ii,jj)
                    fgsum   = fgsum + f(ii,jj)*g(ii,jj)
                end do
            end do
            fbar = fsum/nn
            gbar = gsum/nn

           ! print *,nn,fsum,gsum,fbar,gbar,f2sum,g2sum,fgsum

            dic = (f2sum - nn*fbar*fbar)*(g2sum - nn*gbar*gbar)
            if (dic <= 0) then      !   <0 for -0 error
                dic = 0
                return
            end if

            dic = (fgsum - nn*fbar*gbar) / sqrt(dic)

            return
        end function dic

!        

        subroutine errorExit(message)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            character(len=*),intent(in),optional             ::      message
            if (present(message)) print *,"test_flattenBackground::"//trim(message)
            stop
        end subroutine errorExit
        


    end program test_saltandpepperfilter
