#!/bin/bash
#Joseph Salazar
#salazjos@oregonstate.edu
#proj1bash.sh
#bash script for running openMPproject1.cpp 

#number of threads
for t in 1 2 4 6 8
do
    echo NUMT = $t
    #number of trails
    for s in 200000 400000 600000 800000 1000000
    do
         echo NUMTRIALS = $t
         g++ -DNUMTRIALS=$s -DNUMT=$t openMPproject1.cpp -o prog1 -lm -fopenmp
         ./prog1
    done
done