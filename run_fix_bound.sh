#!/bin/bash

for i in {5..100..1}
do
   t_max=$((20000000/($i*$i*$i)))
   echo "N_X = $i, T_MAX = $t_max"
   ./lwacm s $i t $t_max 
done
