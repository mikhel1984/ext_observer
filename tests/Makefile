
.PHONY: test kalman rnea dyn

# check observers with M, C, G matrices
test: test.cpp
	g++ -Wall test.cpp -o test -lm
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
	
#dyn: try_dyn.cpp
#	g++ -shared -o libtrydin.so -fPIC try_dyn.cpp
