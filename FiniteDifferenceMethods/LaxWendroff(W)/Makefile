CC = gcc
CFLAGS = -I/home/jp/anaconda2/include/python2.7 -Wall -fPIC -O3


LaxWendroffSolver_wrap.c:
	swig -python LaxWendroffSolver.i

LaxWendroffSolver: LaxWendroffSolver_wrap.o
	ld -shared LaxWendroffSolver.o LaxWendroffSolver_wrap.o -o _LaxWendroffSolver.so 

%.o : %.c
	$(CC) -c LaxWendroffSolver.c $^ $(CFLAGS)


clean:
	rm -f ./*.o *~ ./LaxWendroffSolver.py ./LaxWendroffSolver_wrap.c ./*.pyc ./*.so

