#!/bin/bash

mesh=$1
max=50000
min=1000

# uniform material
../../build/bin/gen_vox_mtrl mesh=$mesh type=hexs mode=000 maxv=5000 minv=5000 mtr_file=mtr-0.txt

# # inhomo
# i=1
# for mode in 010 100 001 110 111; do
#     echo 'mtr '$i
#     ../../build/bin/gen_vox_mtrl mesh=$mesh type=hexs mode=$mode maxv=$max minv=$min mtr_file=cube.sub3-mtr-$i.txt outmesh=temp$i.vtk
#     let i=i+1
# done
