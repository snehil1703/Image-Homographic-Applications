#!/bin/bash

make clean
make

# part3

time ./a2 part3 a2-images/lincoln.png

path="a2-images/part3_images/seq1/"
rfile=($path/*)
curr="${rfile[RANDOM % ${#rfile[@]}]}"
declare -a temp
i=0
for n in $path/*
do
    temp[$i]=$n
    ((i++))
done
time ./a2 part3 $curr $( echo ${temp[*]} )

path="a2-images/part3_images/seq2/"
rfile=($path/*)
curr="${rfile[RANDOM % ${#rfile[@]}]}"
declare -a temp
i=0
for n in $path/*
do
    temp[$i]=$n
    ((i++))
done
time ./a2 part3 $curr $( echo ${temp[*]} )
