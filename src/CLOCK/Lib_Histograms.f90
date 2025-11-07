module Lib_Histograms
    !---^^^^^^^^^^^^^^^^^^
    !*-----------------------------------------------------------------------------------------------------------------------------------
    !*    AnalysePngHistogram from the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
    !*    Copyright (C) 2024  Daniel Mason

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

    !*      Generic histograms
        use Lib_Quicksort
        use iso_fortran_env
        implicit none
        private

        public              ::      automatic_xrange
        public              ::      automatic_binwidth
        public              ::      automatic_nBins

        contains
        
        pure real(kind=real64) function automatic_xrange(data) 
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      given the input set of data, choose a sensible histogram range
                real(kind=real64),dimension(:),intent(in)           ::      data !array of spot data in
        
        
            !---    list of possible -x settings
                integer,parameter                               ::      NAUTOX = 25 !number of logarithmic bins
                real(kind=real64),dimension(NAUTOX),parameter   ::      AUTOX = (/  1.0d4,8.0d3,5.0d3,4.0d3,3.0d3,2.0d3,        & !Logarithmic bins
                                                                                    1.0d3,8.0d2,5.0d2,4.0d2,3.0d2,2.0d2,        &
                                                                                    1.0d2,8.0d1,5.0d1,4.0d1,3.0d1,2.0d1,        &
                                                                                    1.0d1,8.0d0,5.0d0,4.0d0,3.0d0,2.0d0,1.0d0   /)             
                real(kind=real64),parameter                     ::      MIN_HIST_CAPTURE = 0.95d0           !   minimum proportion of diameters I want to capture
                integer,dimension(NAUTOX)                       ::      autox_bin_count !logaritmic histogram bin count                                            
                
                integer         ::      ii,jj,nData,minbinocc !diameter index, bin index, number of datapoints, minimum bin occupancy
                
            !---    check for pathological cases                    
                nData = size(data)
                if (nData == 0) then
                    automatic_xrange = 1.0d0
                    return
                else if (nData <= 4) then
                    automatic_xrange = maxval(data) 
                    return
                end if            
                
            !---    find a histogram of the data                        
                autox_bin_count = 0            
                do ii = 1,nData
                    !   find the setting just larger than the data
                    do jj = 2,NAUTOX
                        if (data(ii) > AUTOX(jj-1)) then
                            autox_bin_count(jj-1) = autox_bin_count(jj-1) + 1
                            exit
                        end if
                    end do
                end do 
                
                
            !---    now look at the largest settings found. I want to choose the largest setting 
            !       which has > 1 data ( or 0.1%, whichever is bigger ) in it, subject to constraint that I want to <5% of data points outside range
                minbinocc = max( 1, nint( nData*0.001d0 ) )
                ii = floor( (1-MIN_HIST_CAPTURE) * nData )            !    This is the largest number of data points I'm allowed to ignore                                      
                do jj = 2,NAUTOX
                
                    if (autox_bin_count(jj)>minbinocc) then
                        automatic_xrange = AUTOX(jj-1)
                        return
                    end if
                    
                    ii = ii - autox_bin_count(jj)
                    if (ii<0) then
                        !   If I look for a smaller setting, then I exclude too many data points
                        automatic_xrange = AUTOX(jj)
                        return
                    end if
                        
                end do                                                                    
            
            !--- I don't think its possible to get here, but I have to return something!                                                                                                
                automatic_xrange = 50.0d0
                return
            end function automatic_xrange                       
            
            
            
            
            real(kind=real64) function automatic_binwidth(data)
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      given the input data set , choose a sensible histogram bin width
        !*      using the Freedman-Diaconis' rule
        !*          https://en.wikipedia.org/wiki/Histogram#Number_of_bins_and_width
                real(kind=real64),dimension(:),intent(in)           ::      data !array of data
                
                real(kind=real64),dimension(size(data))         ::      data_sorted_list !array of data    
                integer             ::      nData !number of data points
                
                real(kind=real64)   ::      IQR         !   inter-quartile range
                
            !---    check for pathological cases            
                nData = size(data)            
                if (nData == 0) then
                    automatic_binwidth = 1.0d0
                    return
                else if (nData <= 10) then
                    automatic_binwidth = maxval(data)/2
                    return
                end if
                
            !---    sort the diameters into smallest to largest                            
                data_sorted_list = data 
                call quicksort( data_sorted_list )
                
            !---    find the interquartile range = q3 - q1
                IQR = data_sorted_list( nint( nData*0.75d0 ) ) - data_sorted_list( nint( nData*0.25d0 ) )
                
                automatic_binwidth = 2*IQR/(real(nData)**0.33333)
                
                
                return
            end function automatic_binwidth            
        
            
            pure subroutine automatic_nBins(x,dx,nBins)
        !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        !*      given the max bin size, and an ideal bin width 
        !*      return a nice-looking number of bins. Adjust the bin width as appropriate
                real(kind=real64),intent(in)                    ::      x !max bin value
                real(kind=real64),intent(inout)                 ::      dx !bin width
                integer,intent(out)                             ::      nBins !number of bins
                
                integer                 ::      ii,nn,besti !bin edge index, current number of bins, best number of bins, best bin width index.
                real(kind=real64)       ::      dn,bestdn,nBin_FD
            !---    list of possible dx settings
                integer,parameter                               ::      NAUTODX = 26 !number of logarithmic bin edges
                real(kind=real64),dimension(NAUTODX),parameter  ::      AUTODX = (/ 1.0d4,8.0d3,5.0d3,3.0d3,2.0d3,          & !logarithmic bin edges
                                                                                    1.0d3,8.0d2,5.0d2,3.0d2,2.0d2,          &
                                                                                    1.0d2,8.0d1,5.0d1,3.0d1,2.0d1,          &
                                                                                    1.0d1,8.0d0,5.0d0,3.0d0,2.0d0,          &
                                                                                    1.0d0,0.8d0,0.5d0,0.3d0,0.2d0,0.1d0     /)             
                nBin_FD = x/dx
                nBins = nint( nBin_FD )          !   this is the default
                
            !---    check pathological cases                               
                if (nBins<=5) then
                    dx = x/nBins
                    return
                end if
                
                if (nBins>20) then
                    !   pretty much never want > 20 bins
                    nBins = 20
                    dx = x/nBins
                    return
                end if
                
                        
                
            !---    OK, at this point what I want is a nice bin width.
                !print *,"unscaled x,nBins,dx ",x,nBins,x/nBins
                bestdn = huge(1.0)
                besti = NAUTODX
                do ii = 1,NAUTODX
                    
                    nn = nint( x/AUTODX(ii) )       !   how many bins would I have if I chose this bin with
                    dn = (nn*AUTODX(ii) - x) 
                    if (abs(dn) > AUTODX(ii)*1.0d-4) cycle      !   don't get a good integer number of bins with this choice
                    dn = (nn-nBin_FD)
                    dn = dn*dn                      !   this is the square difference between this number of bins and my stated ideal
                    if (dn <= bestdn) then
                        bestdn = dn
                        besti = ii
                    end if
                end do
                
                nBins = nint( x/AUTODX(besti) )   
                            
                 
                 
                                   
                            
                
                
                dx = x/nBins
                return
            end subroutine automatic_nBins   

        


end module Lib_Histograms