#!/bin/bash
# specific conversion script for my html files to php
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig
i=0
j=0
dest_dir=/net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/feature/perdoch_hesaff_sift
mkdir ${dest_dir}
for frame_fold in /net/per610a/export/das11f/ledduy/plsang/nvtiep/INS/INS2013/query/frames_png/*; do
        i=$(($i+1))
        echo ${i}
        query_id=${frame_fold##/*/}
		dest=${dest_dir}/${query_id}
		test -d ${dest} && continue
        mkdir ${dest}
        for frame in ${frame_fold}/*.png; do
            j=$(($j+1))
            default_feature_fullname=${frame%.png}.png.hesaff.sift
            frame_filename=${frame##/*/}
            feature_fullname=${dest}/${frame_filename%.png}.txt
            /net/per610a/export/das11f/ledduy/plsang/nvtiep/tool/hesaff ${frame} 
            mv ${default_feature_fullname} ${feature_fullname}
        done
done
echo ${i}
echo ${j}
