LFLAGS=-pthread -lgsl -lgslcblas -lm
all: simControlledEscapeTime.cpp parSolver.cpp
	${CXX} -g -std=c++11 simControlledEscapeTime.cpp parSolver.cpp ${LFLAGS}
