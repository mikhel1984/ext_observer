
.PHONY: test kalman rnea

# check observers with M, C, G matrices
test: test.cpp
	g++ -Wall test.cpp -o test -lm
	./test
	gnuplot -p force.gnuplot

# check observers with RNEA
RNEA_SRC = test_rnea.cpp  ../src/robot_dynamics_rnea.cpp
rnea: $(RNEA_SRC)
	g++ -Wall $(RNEA_SRC) -o rnea -lm
	./rnea
	gnuplot -p force.gnuplot	

# check Kalman filter
kalman: test_kalman.cpp
	g++ -Wall test_kalman.cpp -o kalman -lm
	./kalman
	gnuplot -p kalman.gnuplot	
	
