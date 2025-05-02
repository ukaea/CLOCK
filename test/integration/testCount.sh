#!/bin/bash
#   short script to test output of count code.
#   returns "PASS" if the output of the command
#       ./bin/count.exe -f [img_files_directory]/img9.png 
#   gives the same histogram as stored in the groundtruth file [img_files_directory]/img9.groundtruth.count
#   usage:
#       testCount.sh [img_files_directory]

if [ "$#" -ne 1 ]; then
    echo "usage: testCount.sh [img_files_directory]"
    echo "FAIL"
else        
    ../../src/Count/count -quiet -f $1/W_30K_0.2dpa_DF.png | grep "headline result" > countresult
    nspots=`awk '{print $2}' countresult `
    avgdiam=`awk '{print $3}' countresult`
    # now find the error - either the absolute difference in number of spots, or the first place in decimals of the average diameter

    delta=`awk 'function max(a,b){return a>b?a:b} ; function abs(v) {return v < 0 ? -v : v}{print max( abs($2-1780), int(abs($3-12.269)*10) )}' countresult`
    echo "         NEW RESULT       CACHED RESULT"
    echo "n spots  " ${nspots}   " 1780 " 
    echo "<d>      " ${avgdiam}  " 12.269 "

    if [ "${delta}" -le 2 ]; then    
        echo "PASS"        
    else
        echo "FAIL"        
    fi 
    rm countresult   
fi 