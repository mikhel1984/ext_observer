
.PHONY: test kalman rnea

# check observers with M, C, G matrices
test: test.cpp
	g++ -Wall test.cpp -o test -lm 
	./test
	gnuplot -p force.gnuplot

# check observers with RNEA
rnea: test_rnea.cpp
	g++ -Wall test_rnea.cpp -o rnea -lm
	./rnea
	gnuplot -p force.gnuplot	

# check Kalman filter
kalman: test_kalman.cpp
	g++ -Wall test_kalman.cpp -o kalman -lm
	./kalman
	gnuplot -p kalman.gnuplot	
	
