gcc -o3 time-series.c -lm -o time-series.ce
gcc -o3 param-scan.c -lm -o param-scan.ce
gcc -o3 zoom-scan.c -lm -o zoom-scan.ce
gcc -o3 stoch.c -lm -o stoch.ce
gcc -o3 hill.c -lm -o hill.ce

./param-scan.ce > tmp-p &
./zoom-scan.ce > tmp-z &
./time-series.ce > tmp-t &
./stoch.ce > tmp-s &
./hill.ce > tmp-h &
