#/bin/bash
#
# Warning: You must adapt this script to your runtime environment.

ulimit -s unlimited

gcc -L$GSL_ROOT/lib -I$GSL_ROOT/include SGBM_ArithBasket.c -lgsl -lgslcblas -lm -O3
# gcc SGBM_ArithBasket.c -lgsl -lgslcblas -lm -O3

#for((i=0;i<30;i++))
#do
#./a.out
#done
