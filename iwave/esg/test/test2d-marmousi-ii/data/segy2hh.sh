#!/bin/bash
#   This script is used to get 2D elastic Marmousi model from the website,
#   and convert it to rsf format.

#   download file from http://www.agl.uh.edu/downloads/downloads.htm
curl -O http://www.agl.uh.edu/downloads/modeling_grids/MODEL_S-WAVE_VELOCITY_1.25m.segy.zip
curl -O http://www.agl.uh.edu/downloads/modeling_grids/MODEL_P-WAVE_VELOCITY_1.25m.segy.zip
curl -O http://www.agl.uh.edu/downloads/modeling_grids/MODEL_DENSITY_1.25m.segy.zip

#   unzip files
for file in *.zip
do
unzip $file
rm $file
done

#   convert to rsf format
for file in *.segy
do
sfsegyread < $file tfile=tfile.rsf hfile=hfile.rsf bfile=bfile.rsf > "${file%segy}hh" --out=stdout
rm $file
done

#   rm t(/h/b)file.rsf
rm *file*
