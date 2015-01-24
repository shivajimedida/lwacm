#!/bin/bash

for i in $(seq 5 100)
do
   t_max=$((10000000/($i*$i*$i)))
   echo "N_X = $i, T_MAX = $t_max"
   #./lwacm s $i t $t_max 
done

