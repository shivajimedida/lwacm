#!/bin/bash

for i in {5..100..1}
do
   t_max=$((20000000/($i*$i*$i)))
   echo "N_X = $i, T_MAX = $t_max"
   gcc -DT_MAX=$t_max -DN_X=$i -DN_Y=$i -DN_Z=$i -O3 lwacm_raw.c -o lwacm
   ./lwacm
done


