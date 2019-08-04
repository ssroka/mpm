!/bin/bash

for d in */; do
	grep -il "ppm" $d*.m >> files_w_ppm.txt
	grep -il "ppt" $d*.m >> files_w_ppt.txt
done

for d in */*/; do
	grep -il "ppm" $d*.m >> files_w_ppm.txt
	grep -il "ppt" $d*.m >> files_w_ppt.txt
done




