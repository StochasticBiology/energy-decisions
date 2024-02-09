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

gcc -o3 stoch-scan.c -lm -o stoch-scan.ce
./stoch-scan.ce 0 > tmp-ss-0 &
./stoch-scan.ce 1 > tmp-ss-1 &
./stoch-scan.ce 2 > tmp-ss-2 &
./stoch-scan.ce 3 > tmp-ss-3 &
./stoch-scan.ce 4 > tmp-ss-4 &
./stoch-scan.ce 5 > tmp-ss-5 &
./stoch-scan.ce 6 > tmp-ss-6 &
./stoch-scan.ce 7 > tmp-ss-7 &
./stoch-scan.ce 8 > tmp-ss-8 &

gcc -o3 param-matrix.c -lm -o param-matrix.ce
./param-matrix.ce 2 0 > tmp-pm-0 &
./param-matrix.ce 2 1 > tmp-pm-1 &
./param-matrix.ce 2 2 > tmp-pm-2 &
./param-matrix.ce 2 3 > tmp-pm-3 &
