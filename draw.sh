#!/bin/bash

for i in log_*; do
    octave draw.m $i
done


