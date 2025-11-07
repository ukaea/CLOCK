#!/bin/bash
#Integration test to check the F1score program is working.

counter=0

#test circular F1.
f1a="`../../src/Tools/f1score -f ../../../data/F1Dat/F1A.dat -g ../../../data/F1Dat/F1B.dat | grep -A 1 "nSpots" | tail -n1 | awk '{print $(NF-1)}'`"
f1b="`../../src/Tools/f1score -f ../../../data/F1Dat/F1C.dat -g ../../../data/F1Dat/F1B.dat | grep -A 1 "nSpots" | tail -n1 | awk '{print $(NF-1)}'`"
f1c="`../../src/Tools/f1score -f ../../../data/F1Dat/F1D.dat -g ../../../data/F1Dat/F1B.dat | grep -A 1 "nSpots" | tail -n1 | awk '{print $(NF-1)}'`"

if [[ "$f1a" == "1.00000000" ]]; then
    counter=$((counter+1))
    echo "testF1score: CORRECT F1 value of $f1a for 5 overalpping spots using circular overalap."
else
    echo "testF1score: INCORRECT F1 value of $f1a for 5 overalpping spots using circular overalap."
    fi
if [[ "$f1b" == "0.00000000" ]]; then
    counter=$((counter+1))
    echo "testF1score: CORRECT F1 value of $f1b for 0 overalpping spots using circular overalap."
else
    echo "testF1score: INCORRECT F1 value of $f1b for 0 overalpping spots using circular overalap."
    fi
if [[ "$f1c" == "0.75000000" ]]; then
    counter=$((counter+1))
    echo "testF1score: CORRECT F1 value of $f1c for 3 true +ve, 1 false -ve and 1 false +ve using circular overalap."
else
    echo "testF1score: INCORRECT F1 value of $f1c for 4/6 overalpping spots using circular overalap."
    fi

#test rectangular F1.
f2a="`../../src/Tools/f1score -f ../../../data/F1Dat/F1A.dat -g ../../../data/F1Dat/F1B.dat -m 1 | grep -A 1 "nSpots" | tail -n1 | awk '{print $(NF-1)}'`"
f2b="`../../src/Tools/f1score -f ../../../data/F1Dat/F1C.dat -g ../../../data/F1Dat/F1B.dat -m 1 | grep -A 1 "nSpots" | tail -n1 | awk '{print $(NF-1)}'`"
f2c="`../../src/Tools/f1score -f ../../../data/F1Dat/F1D.dat -g ../../../data/F1Dat/F1B.dat -m 1 | grep -A 1 "nSpots" | tail -n1 | awk '{print $(NF-1)}'`"

if [[ "$f2a" == "1.00000000" ]]; then
    counter=$((counter+1))
    echo "testF1score: CORRECT F1 value of $f2a for 5 overalpping spots using rectangular overalap."
else
    echo "testF1score: INCORRECT F1 value of $f2a for 5 overalpping spots using rectangular overalap."
    fi
if [[ "$f2b" == "0.00000000" ]]; then
    counter=$((counter+1))
    echo "testF1score: CORRECT F1 value of $f2b for 0 overalpping spots using rectangular overalap."
else
    echo "testF1score: INCORRECT F1 value of $f2b for 0 overalpping spots using rectangular overalap."
    fi
if [[ "$f2c" == "0.75000000" ]]; then
    counter=$((counter+1))
    echo "testF1score: CORRECT F1 value of $f2c for 3 true +ve, 1 false -ve and 1 false +ve using rectangular overalap."
else
    echo "testF1score: INCORRECT F1 value of $f2c for 4/6 overalpping spots using rectangular overalap."
    fi

#test numerical ellipse F1.
f3a="`../../src/Tools/f1score -f ../../../data/F1Dat/F1A.dat -g ../../../data/F1Dat/F1B.dat -m 2 | grep -A 1 "nSpots" | tail -n1 | awk '{print $(NF-1)}'`"
f3b="`../../src/Tools/f1score -f ../../../data/F1Dat/F1C.dat -g ../../../data/F1Dat/F1B.dat -m 2 | grep -A 1 "nSpots" | tail -n1 | awk '{print $(NF-1)}'`"
f3c="`../../src/Tools/f1score -f ../../../data/F1Dat/F1D.dat -g ../../../data/F1Dat/F1B.dat -m 2 | grep -A 1 "nSpots" | tail -n1 | awk '{print $(NF-1)}'`"

if [[ "$f3a" == "1.00000000" ]]; then
    counter=$((counter+1))
    echo "testF1score: CORRECT F1 value of $f3a for 5 overalpping spots using numerical ellipse overalap."
else
    echo "testF1score: INCORRECT F1 value of $f3a for 5 overalpping spots using numerical ellipse overalap."
    fi
if [[ "$f3b" == "0.00000000" ]]; then
    counter=$((counter+1))
    echo "testF1score: CORRECT F1 value of $f3b for 0 overalpping spots using numerical ellipse overalap."
else
    echo "testF1score: INCORRECT F1 value of $f3b for 0 overalpping spots using numerical ellipse overalap."
    fi
if [[ "$f3c" == "0.75000000" ]]; then
    counter=$((counter+1))
    echo "testF1score: CORRECT F1 value of $f3c for 3 true +ve, 1 false -ve and 1 false +ve using numerical ellipse overalap."
else
    echo "testF1score: INCORRECT F1 value of $f3c for 4/6 overalpping spots using numerical ellipse overalap."
    fi

 echo "$counter/9 programs ran successfully"
if [ $counter -eq 9 ]; then
        echo "PASS"
    fi