#!/bin/bash
#   short script to test output of bfdf code.
#   returns "PASS" if the output of the command
#       ../bin/bfdf -f "$1/img1$i.png"
#   gives the same result Dark/Bright/Mixed as stored in the groundtruth file [img_files_directory]/bfdf.groundtruth.dat
#   usage:
#       testBfdf.sh [img_files_directory]
if [ "$#" -ne 1 ]; then
    echo "usage: testBfdf.sh [img_files_directory]"
    echo "FAIL"
else
    if [ -f bfdf.out ]; then
        rm bfdf.out
    fi
    declare -a test_images=("fe6cr_n_irrad_1.8dpa_fib_dam_BF.png" 
                            "fe6cr_n_irrad_1.8dpa_fib_dam_DF.png" 
                            "filtertest.png"
                            "flattest.png"
                            "test_cameraman.png"
                            "test_cat.png"
                            "test_fracture.png"
                            "test_grey.png"
                            "test_line.png"
                            "test_squares.png"
                            "test_synth_ms.png"
                            "W_30K_0.2dpa_DF.png")
    #for (( i = 0; i <= 12; i++ )); do echo -n "$1$test_images[i].png" " "  ; ../../src/BFDF/bfdf -f "$1/img1$i.png" | grep "AIC (best)" | awk '{print $2}' >> bfdf.out ; tail -n 1 bfdf.out ; done
    for i in "${test_images[@]}"; do echo -n $1$i; printf "\n"; ../../src/BFDF/bfdf -f "$1$i" | grep "AIC (best)" | awk '{print $2}' >> bfdf.out ; tail -n 1 bfdf.out ; done
    
    echo "NEW RESULT       CACHED RESULT"
    paste bfdf.out $1/bfdf.groundtruth.dat
    if [ -z `diff --strip-trailing-cr bfdf.out $1/bfdf.groundtruth.dat` ]; then
        echo "PASS"
        
    else
        echo "FAIL"        
    fi 
    rm bfdf.out 

fi

