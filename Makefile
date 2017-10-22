CXX=g++

all: szum.o
	${CXX} -o szuuuumy szum.o

szum.o:
	${CXX} -c szum.cpp

clean:
	rm -f szum.o szuuuumy gauss.dat szuum.dat
