#!/bin/bash

for i in {2..200..1}
do
   t_max=$((24000000/($i*$i*$i)))
   echo "N_X = $i, T_MAX = $t_max"
   ./lwacm s $i t $t_max 
done

