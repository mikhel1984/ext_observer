CC = g++
CFLAGS = -Wall 

MDIR = M_mtrx_fcn/M_mtrx_fcn
CDIR = C_mtrx_fcn/C_mtrx_fcn
GDIR = G_vctr_fcn/G_vctr_fcn
FDIR = F_vctr_fcn/F_vctr_fcn
       
SRC = bin/main.o bin/robot.o bin/mfnc.o bin/minit.o bin/mterm.o \
      bin/cfnc.o bin/cinit.o bin/cterm.o bin/gfnc.o bin/ginit.o \
      bin/gterm.o bin/ffnc.o bin/finit.o bin/fterm.o 

SRC_TCP = bin/listener.o bin/robot.o bin/mfnc.o bin/minit.o bin/mterm.o \
      bin/cfnc.o bin/cinit.o bin/cterm.o bin/gfnc.o bin/ginit.o \
      bin/gterm.o bin/ffnc.o bin/finit.o bin/fterm.o 


.PHONY: all test clear listen

# get information from csv file
all: ${SRC}
	${CC} ${SRC} -o test -lm 

# get information from TCP connection
listen: ${SRC_TCP}
	${CC} ${SRC_TCP} -o listen -lm 

# remove object files
clear:
	rm -f bin/*.o

# make object files
bin/main.o: test.cpp
	${CC} ${CFLAGS} test.cpp -c -o $@
	
bin/listener.o: test_listen.cpp
	${CC} ${CFLAGS} test_listen.cpp -c -o $@

bin/robot.o: ur_robot.cpp 
	${CC} ${CFLAGS} ur_robot.cpp -c -o $@
	
bin/mfnc.o: ${MDIR}.c
	${CC} ${CFLAGS} ${MDIR}.c -c -o $@
	
bin/minit.o: ${MDIR}_initialize.c
	${CC} ${CFLAGS} ${MDIR}_initialize.c -c -o $@

bin/mterm.o: ${MDIR}.c
	${CC} ${CFLAGS} ${MDIR}_terminate.c -c -o $@
	
bin/cfnc.o: ${CDIR}.c
	${CC} ${CFLAGS} ${CDIR}.c -c -o $@
	
bin/cinit.o: ${CDIR}_initialize.c
	${CC} ${CFLAGS} ${CDIR}_initialize.c -c -o $@

bin/cterm.o: ${CDIR}.c
	${CC} ${CFLAGS} ${CDIR}_terminate.c -c -o $@
	
bin/gfnc.o: ${GDIR}.c
	${CC} ${CFLAGS} ${GDIR}.c -c -o $@
	
bin/ginit.o: ${GDIR}_initialize.c
	${CC} ${CFLAGS} ${GDIR}_initialize.c -c -o $@

bin/gterm.o: ${GDIR}.c
	${CC} ${CFLAGS} ${GDIR}_terminate.c -c -o $@
	
bin/ffnc.o: ${FDIR}.c
	${CC} ${CFLAGS} ${FDIR}.c -c -o $@
	
bin/finit.o: ${FDIR}_initialize.c
	${CC} ${CFLAGS} ${FDIR}_initialize.c -c -o $@

bin/fterm.o: ${FDIR}.c
	${CC} ${CFLAGS} ${FDIR}_terminate.c -c -o $@



