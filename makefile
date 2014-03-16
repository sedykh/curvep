CC=g++
CFLAGS=-c -m32 -O3 -Wall -static
all: curvep

curvep: curvep.o apstring.o set.o core.o stack.o
	$(CC) curvep.o apstring.o set.o core.o stack.o -m32 -o curvep

apstring.o: apstring.cpp
	$(CC) $(CFLAGS) apstring.cpp
set.o:
	$(CC) $(CFLAGS) set.cpp

core.o: core.cpp
	$(CC) $(CFLAGS) core.cpp

stack.o: stack.cpp
	$(CC) $(CFLAGS) stack.cpp

curvep.o: curvep.cpp
	$(CC) $(CFLAGS) curvep.cpp

clean:
	rm -rf *.o
