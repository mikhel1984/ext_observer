INCDIR = -I /usr/include/eigen3

.PHONY: test kalman rnea dyn

# check observers with M, C, G matrices
TEST_SRC = test.cpp  ../lib/kalman_filter_continous.cpp
test: ${TEST_SRC}
	g++ -Wall ${INCDIR} ${TEST_SRC} -o test -lm
	./test
	gnuplot -p force.gnuplot

# check observers with RNEA
RNEA_SRC = test_rnea.cpp  ../lib/robot_dynamics_rnea.cpp
rnea: $(RNEA_SRC)
	g++ -Wall $(RNEA_SRC) -o rnea -lm
	./rnea
	gnuplot -p force.gnuplot	

# check Kalman filter
kalman: test_kalman.cpp
	g++ -Wall test_kalman.cpp -o kalman -lm
	./kalman
	gnuplot -p kalman.gnuplot	

# build dynamic library
dyn: observers.cpp
	g++ -shared -o libobservers.so -fPIC observers.cpp
