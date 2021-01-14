#!/bin/bash
 
# script file for BASH 
# which bash
# save this file as g.sh
# chmod +x g.sh
# ./g.sh
##  checked in https://www.shellcheck.net/
## shellcheck ./g.sh

CURRENTDATE=$(date +"%Y-%m-%d-%T")

echo "make ppm files"
gcc d.c -lm -Wall -march=native -fopenmp
time ./a.out > "${CURRENTDATE}".txt

echo "convert all ppm files to png using Image Magic convert"
# for all ppm files in this directory
for file in *.ppm ; do
  # read comment 
  c=$(identify -verbose "$file"|grep "comment")	
  # b is name of file without extension
  b=$(basename "$file" .ppm)
  # convert  using ImageMagic
  convert "${b}".ppm -resize 2000x2000  -set comment "${c}" "${b}".png
  echo "$file"
done

# display OpenMP info
export  OMP_DISPLAY_ENV="TRUE"	

echo  "delete all ppm files"
rm ./*.ppm




echo OK
# end
