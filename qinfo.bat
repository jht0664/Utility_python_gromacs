#!/bin/bash -f

A=$(qstat -u htjung | grep htjung > qstat.imsi ; wc qstat.imsi | awk '{print$1}' ; rm qstat.imsi)
B=$(qstat -u htjung | grep ' R ' > qstat.imsi ; wc qstat.imsi | awk '{print$1}' ; rm qstat.imsi)
C=$(qstat -u htjung | grep ' Q ' > qstat.imsi ; wc qstat.imsi | awk '{print$1}' ; rm qstat.imsi)
D=$(qstat -u htjung | grep ' E ' > qstat.imsi ; wc qstat.imsi | awk '{print$1}' ; rm qstat.imsi)
echo
echo "**********************************"
echo "  Total number of Process: $A"
echo "  Number of Wokrs in R:    $B"
echo "  Number of Wokrs in Q:    $C"
echo "  Number of Wokrs in E:    $D"
echo "**********************************"
echo
