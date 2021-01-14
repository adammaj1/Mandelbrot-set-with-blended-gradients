#!/bin/bash
 
# script file for BASH 
# which bash
# save this file as g.sh
# chmod +x g.sh
# ./g.sh
# checked in https://www.shellcheck.net/
# shellcheck ./g.sh



echo "make ppm files"
gcc d.c -lm -Wall -march=native -fopenmp
time ./a.out > a.txt

echo "convert all ppm files to png using Image Magic convert"
# for all ppm files in this directory
for file in *.ppm ; do
  # b is name of file without extension
  b=$(basename "$file" .ppm)
  # convert  using ImageMagic
  convert "${b}".ppm -resize 2000x2000 "${b}".png
  echo "$file"
done



echo  "delete all ppm files"
rm ./*.ppm


# display OpenMP info
export  OMP_DISPLAY_ENV="TRUE"	

echo OK
# end
