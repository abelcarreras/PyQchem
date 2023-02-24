#!/bin/bash

module load Python/3.7.4-Anaconda3-2019.10
molecule="cucl4_2-"
#reference="d10"
#basis="6-31Gd def2-TZVP"

for file in $molecule/$molecule*/*.out; do
 echo $file | cut -f3 -d "/"
 python g_main.py $file
done

