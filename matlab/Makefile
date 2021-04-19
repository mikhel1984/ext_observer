CC = g++
CFLAGS = -Wall
INCDIR = -I /usr/include/eigen3
      
SRC = bin/observers.o bin/fcont.o 
      
.PHONY: all clear test

# get information from csv file
all: ${SRC}
	${CC} ${SRC} -shared -o observers.so 

# remove object files
clear:
	rm -f bin/*.o

# make object files

bin/observers.o: observers.cpp
	${CC} ${CFLAGS} observers.cpp -c -fPIC -o $@
	
bin/fcont.o: ../lib/kalman_filter_continous.cpp
	${CC} ${CFLAGS} ../lib/kalman_filter_continous.cpp -c -fPIC -o $@	
	
test: test_lib.cpp
	${CC} ${CFLAGS} -Wl,-rpath=$(shell pwd) test_lib.cpp -o test_lib  -L. -l:observers.so
	./test_lib
	

	



