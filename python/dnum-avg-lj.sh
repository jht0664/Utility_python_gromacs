#!/bin/bash
# make average 1d-number profile from .dnums file
# $1 : begin frame
# $2 : tolerance for block avg
i1="a"
i2="b"

inputgen=".dnums.align"
avggen=".dnums"

input1=$i1$inputgen
avg1=$i1$avggen
output1=$i1$outputgen
input2=$i2$inputgen
output2=$i2$outputgen
avg2=$i2$avggen

if [ -f "$input1" ]
then
        echo "$input1 found."
else
        echo "$input1 not found."
        exit 1
fi

if [ -f "$input2" ]
then
        echo "$input2 found."
else
        echo "$input2 not found."
        exit 1
fi


## extract trajectory
python3 ~/Utility/python/dnum-avg.py -i $input1 -avg $avg1 -b $1 -tol $2 | tee $il1.dnum-avg.log
python3 ~/Utility/python/dnum-avg.py -i $input2 -avg $avg2 -b $1 -tol $2 | tee $il2.dnum-avg.log
