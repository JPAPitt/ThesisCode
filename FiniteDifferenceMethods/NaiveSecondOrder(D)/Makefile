CC = gcc
CFLAGS = -I/home/jp/anaconda2/include/python2.7 -Wall -fPIC -O3


NaiveSecondOrderSolver_wrap.c:
	swig -python NaiveSecondOrderSolver.i

NaiveSecondOrderSolver: NaiveSecondOrderSolver_wrap.o
	ld -shared NaiveSecondOrderSolver.o NaiveSecondOrderSolver_wrap.o -o _NaiveSecondOrderSolver.so 

%.o : %.c
	$(CC) -c NaiveSecondOrderSolver.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./NaiveSecondOrderSolver.py ./NaiveSecondOrderSolver_wrap.c ./*.pyc ./*.so

