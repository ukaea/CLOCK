# The Culham LOop Counting Kit (CLOCK)
[![Build](https://github.com/ukaea/CLOCK/workflows/CI/badge.svg)](https://github.com/ukaea/CLOCK/actions)
[![codecov](https://codecov.io/gh/ukaea/CLOCK/branch/RSE_get_CI_working/graph/badge.svg?token=NR07LOGXI4)](https://codecov.io/gh/ukaea/CLOCK)

 A code for analysing TEM micrographs, it is run via command line. This is a preview build preceding a methods paper. There will be some additions/changes before submission. If you use CLOCK before the updated methods paper come out, please cite:
 - For background subtraction and black spot detection: Mason et al., Acta Materialia, Volume 144, February 2018, Pages 905-917
 - For denoising with the maximum likelihood filter: D. R. Mason & A. London, (2020) Ultramicroscopy 211:112940

## Codes included
    
    "bfdf" (built from: AnalysePngHistogram.f90)- decide whether an input .png file is dark-signal-on-bright-background or vice versa.
    "flat" (built from: flattenBackground.f90)- a collection of subroutines for filtering/cleaning up images.
    "count" (built from: analyseSpots.f90 countSpots.f90)- Detects/counts bright spot features.   

## Requirements

CLOCK runs on Linux (Ubuntu 22.04.1 LTS), Windows users can use Windows Subsystem for Linux 2 (WSL2) or Cygwin. Much of CLOCK's recent development has been performed on WSL2 but Cygwin has worked for older builds. In theory CLOCK should work on macOS but there may be some troubleshooting of paths required.

    You will need cmake
    You will need a fortran and a c compiler.
    You will need openmp.
    You will need libpng.
    You will also need LAPACK, BLAS, FFTW3
NB: pkgconf is optional but advised as it will automate locating FFTW3.

## To compile

    >   cd Clock_download_dir
    >   mkdir build
    >   cd build
    >   cmake ..
    >   make

Lapack and libpng should be detected automatically.
The fftw3 library should be located by pkgconf. Should that fail or pkgconf is unavailable to you, FFTW3 can added by hand. 
Note that the file "fftw3.f03" is included by src/Lib_FFTW3.f90. This is typically found in the /usr/include/ directory or similar.

NB: to compile for debugging (assumes you have compiled already, if not remove the "make clean"):

    >   cd build
    >   make clean
    >   cmake -DCMAKE_BUILD_TYPE=Debug ..
    >   make

There are debug config files in .vscode/launch.json for those who want to use vscode's interactive debugger.

## To test

    >   cd build
    >   ctest

- To run an individual test:

    >   cd build
    >   ctest -R NameOfTest

where NameOfTest is is defined by the NAME property in add_test called in the /test/testType CMakeLists.

- Some test executables contain multiple unit tests (e.g. test_LibGraphs), to see individual unit tests run in verbose mode:
    >   cd build
    >   ctest -V -R NameOfTest

- To run the short tests:

    >   cd build
    >   ctest -L SHORT

Note: for vscode users there is launch.json in /.vscode that contains launch configurations if you want to use vscode's debugger.

Some users have reported tests failing due to windows end of line character being inserted into shell scripts. Error messages may look like:
- $'\r': command not found
- ")syntax error: invalid arithmetic operator (error token is "
- syntax error near unexpected token `fi'

To fix, please run dos2unix test/\*/*.sh from the /CLOCK directory.

## Installation (optional)

If desired CLOCK can be installed from the project root directory with:

    >   cmake --install build --config Release

## Usage

CLOCK is designed to run on 16-bit greyscale pngs. 8-bit greyscale pngs will also work but the former is recommended. If you are working with 32-bit images (formats such as .dm3 or .dm4), these can be converted with imageJ. Images must be free of annotations. If you have a colour micrograpgh please insure you convert it to greyscale using equal weights for the red/gren/blue channels, e.g '(r+g+b)/3'.

NB:    
- Run programs from Clock directory.
- Run programs with -h, -help flags or without arguments for help.


### bfdf

Code to determine whether features are bright on a dark background or vice-versa.

This is done by analysing the histogram of pixel intensities as a sum of Gaussians.

To determine the probability that a fit is correct:
    $`P_{unnormalised} = exp[(AIC_{0} - AIC_{j})/2]`$
where $`AIC_{0}`$ is the best aic value we find and $`AIC_{j}`$ is the value for the jth trial set of Gaussians.


We imagine that the background is one Gaussian contribution to the histogram
and that the signal is the other Gaussians
(note: the background may not be > 50% of the signal, rather I will say that the probability that peak j is background is proportional to the area under j).

We then say that a secondary peak gives a weighting to dark signal given by its area left of the mean of peak j and a weighting to bright signal 
given by its area to the right of peak j.
This gives a probability for each fit to the data that each part of the signal is bright or dark,
multiplied together with the probability that the total fit is correct, and that peak j is the background. 
Sum these together to give a total probability that the signal is bright or dark.

src/AnalysePngHistogram.f90 is the main code
- checks the command line arguments
- reads in the .png file
- constructs the histogram of intensities
- calls src/Lib_FitGaussians1d to generate a whole bunch of trial fits
- finds the probability of bright-feature-on-dark

*Arguments*

         -f <char>                    input filename
        [ -n <int> ]                    number of histogram bins [ default 256 ]
        [ -m <int> ]                    number of gaussians to try [ default 10 ]

*example*:

    >   ./build/bin/bfdf -f data/fe6cr_n_irrad_1.8dpa_fib_dam_BF.png | grep "best"

*output*:

 AnalysePngHistogram.exe info - best multiple gaussian result
"data/fe6cr_n_irrad_1.8dpa_fib_dam_BF.png" Bright background  72.09% AIC (best)  -3275.831

*NB:* 

    >   ./build/bin/bfdf -f data/fe6cr_n_irrad_1.8dpa_fib_dam_BF.png > anameofyourchoosing.txt
    
for additional fitting info.


### flat

flattenBackground.f90 is a utility program designed to help feature detection in greyscale images by removing background fluctuations and performing a couple of simple cleaning tasks.

*Arguments*

         -f              filename 
         [ -o ]          output filename [ default ".flat.png" ] 
         [ -zero ]       interpret pixels with intensity zero as 'unset' rather than just very dark [ default T ] 
         [ -sandp ]      apply salt and pepper filter [ default F ] 
         [ -lambda ]     set range for max likelihood filter 
         [ -negative ]   invert image [ default F ] 
         [ -half ]       compute background using 1/2 size image (default T if Nx>2000) [ default F ] 
         [ -quarter ]    compute background using 1/4 size image (default T if Nx>4000) [ default F ] 
         [ -g ]          number of gaussian blurs to use for length scale estimate [ default 20 ] 
         [ -n ]          number of iterations [ default 2 ] 
         [ -p ]          power law for gradient contributions [ default 1 ] 
         [ -flat ]       flatten image using scale invariant background subtraction [ default T ] 
         [ -retone ]     retone image to preserve foreground/background intensity [ default F ] 
         [ -hpf ]        high pass filter - Fourier coefficients to remove [ default 0 ] 
         [ -vRC ]        minimum variable-Ridler-Calvard length scale

*example*:

Using a test fe6cr_n_irrad_1.8dpa_fib_dam_BF.png. Note that we want to preserve bright features on a dark background as this is a darkfield micrograph. For brightfield 
images use the -negative parameter to invert the image.

Filter the image with no background subtraction by removing salt-and-pepper noise and using a 4 px maximum likelihood filter and removing first two fourier coefficients.

    >   ./build/bin/flat -f data/fe6cr_n_irrad_1.8dpa_fib_dam_DF.png -o data/fe6cr_n_irrad_1.8dpa_fib_dam_DF.sandp.png -sandp -lambda 2.0 -noflat -hpf 2

Perform a background subtraction ( best to start from the filtered image if you have it, but you can do it all in one step ) assuming intensity zero pixels are intentional and using four loops.

    >   ./build/bin/flat -f data/fe6cr_n_irrad_1.8dpa_fib_dam_DF.sandp.png  -o data/fe6cr_n_irrad_1.8dpa_fib_dam_DF.flat.png -nozero -n 4

*output*:

Now compare fe6cr_n_irrad_1.8dpa_fib_dam_DF.png, fe6cr_n_irrad_1.8dpa_fib_dam_DF.sandp.png and fe6cr_n_irrad_1.8dpa_fib_dam_DF.flat.png.

### count    

Count detects bright spot features (such as an inverted brightfield image of "black spot" radiation damage), 
fits 2D gaussians to the spots and produces a size frequency histogram

*Arguments*

    file handling
        -f <char>                    input filename
        [-o <char> ]                 output filename
        [-png ]                      output .png of spot locations [ default F ]
        [-roi ]                      output roi map [ default F ]
        [-randc ]                    output Ridler&Calvard threshold map [ default F ]

    input image parameters
        [-s <float> ]                image resolution (nm per px) [ default 1.000 ]
        [-negative ]                 input negative image [ default F ]

    detection parameters
        [-t <float> ]                t* threshold [ default 0.000 ]
        [-d <float> ]                intensity dip between separable maxima [ default 0.8000 ]
        [-e <float> ]                max eccentricity [ default 3.000 ]
        [-b <float> ]                intensity threshold (f0/sigma) ( 0 for automatic ) [ default 0.000 ]
        [-tond <float> ]             t/d threshold [ default 0.000 ]
        [-q <float> ]                q threshold [ default 1.000 ]
        [-m <float> ]                maximum fraction of foreground pixels - ( 0 for automatic ) [ default 0.000 ]
        [-X <int> ]                  maximum roi dimension (px) [ default 200 ]
        [-R <float> ]                maximum roi padding (px) [ default 2.000 ]

    output options
        [-x <float> ]                max histogram bin size, if set to 0 CLOCK will find a good value [ default 50.00 ] 
        [-nbins <int> ]              number of histogram bins, note overridden if -x=0 [ default 20 ]
        [-scale <float> ]            diameter scaling [ default 2.000 ]
        [-a ]                        report avg diameter in histogram [ default F ]
        [-colour <char> ]            colour scale [ default "IBM" ]

    debug / misc
        [-dbg <float_array> ]        debug ROI in vicinity of point (x,y) [ default 0.0000,0.0000 ]
        [-count ]                  count spots [ default T ]

*example*:

NB: using the filtered image from the previous example.

    >   ./build/bin/count -f data/fe6cr_n_irrad_1.8dpa_fib_dam_DF.flat.png -png -x 0

*output*: 

* fe6cr_n_irrad_1.8dpa_fib_dam_DF.flat.count - size frequency histogram of detected spots.
* fe6cr_n_irrad_1.8dpa_fib_dam_DF.flat.spots - details of every spot detected.
* fe6cr_n_irrad_1.8dpa_fib_dam_DF.flat.spots.png - img9.flat.png with detected spots annotated. 

### Change log

4.2.0 --> 4.3.0 Default spot diameter colour scale set to "IBM", user options to change colour scale. start of implementation of new linear feaure code.

4.1.0 --> 4.2.0 Changed default behaviour of Count to automatically set the spot detection threshold to three.