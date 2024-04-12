gcc -o3 grn-sim.c -lm -o grn-sim.ce

./grn-sim.ce 0 > tmp-0 &
./grn-sim.ce 1 1 > tmp-1.1 &
./grn-sim.ce 1 2 > tmp-1.2 &
./grn-sim.ce 1 3 > tmp-1.3 &
./grn-sim.ce 1 4 > tmp-1.4 &
# ./grn-sim.ce 2 1 > tmp-2.1 &
./grn-sim.ce 2 2 > tmp-2.2 &
# ./grn-sim.ce 2 3 > tmp-2.3 &
# ./grn-sim.ce 2 4 > tmp-2.4 &
./grn-sim.ce 3 > tmp-3 &
./grn-sim.ce 4 > tmp-4 &
./grn-sim.ce 5 > tmp-5 &
./grn-sim.ce 6 > tmp-6 &
./grn-sim.ce 7 > tmp-7 &

rm grn-sim-8.csv
./grn-sim.ce 8 0 > tmp-8.0 &
./grn-sim.ce 8 1 > tmp-8.1 &
./grn-sim.ce 8 2 > tmp-8.2 &
./grn-sim.ce 8 3 > tmp-8.3 &
./grn-sim.ce 8 4 > tmp-8.4 &
./grn-sim.ce 8 5 > tmp-8.5 &
./grn-sim.ce 8 6 > tmp-8.6 &
./grn-sim.ce 8 7 > tmp-8.7 &
