#!/bin/bash
#unit tests for chekcing all executables in bin run with just the version (-v) flag

counter=0

[[ ! -z `./../../src/Count/count -v | grep "version"` ]] && counter=$((counter+1))
[[ ! -z `./../../src/Flat/flat -v | grep "version"` ]] && counter=$((counter+1))
[[ ! -z `./../../src/BFDF/bfdf -v | grep "version"` ]] && counter=$((counter+1))
[[ ! -z `./../../src/Linear/linear -v | grep "version"` ]] && counter=$((counter+1))
 echo "$counter/4 programs ran successfully"
if [ $counter -eq 4 ]; then
        echo "PASS"
    fi