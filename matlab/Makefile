CC = g++
CFLAGS = -Wall 

URDIR = ../ur
MDIR = ${URDIR}/M_mtrx_fcn/M_mtrx_fcn
CDIR = ${URDIR}/C_mtrx_fcn/C_mtrx_fcn
GDIR = ${URDIR}/G_vctr_fcn/G_vctr_fcn
FDIR = ${URDIR}/F_vctr_fcn/F_vctr_fcn
       
SRC = bin/observers.o bin/robot.o bin/mfnc.o bin/minit.o bin/mterm.o \
      bin/cfnc.o bin/cinit.o bin/cterm.o bin/gfnc.o bin/ginit.o \
      bin/gterm.o bin/ffnc.o bin/finit.o bin/fterm.o 

.PHONY: all clear

# get information from csv file
all: ${SRC}
	${CC} ${SRC} -shared -o observers.so 

# remove object files
clear:
	rm -f bin/*.o

# make object files

bin/observers.o: observers.cpp
	${CC} ${CFLAGS} observers.cpp -c -fPIC -o $@
	
bin/robot.o: ${URDIR}/ur_robot.cpp 
	${CC} ${CFLAGS} ${URDIR}/ur_robot.cpp -c -fPIC -o $@
	
bin/mfnc.o: ${MDIR}.c
	${CC} ${CFLAGS} ${MDIR}.c -c -fPIC -o $@
	
bin/minit.o: ${MDIR}_initialize.c
	${CC} ${CFLAGS} ${MDIR}_initialize.c -c -fPIC -o $@

bin/mterm.o: ${MDIR}.c
	${CC} ${CFLAGS} ${MDIR}_terminate.c -c -fPIC -o $@
	
bin/cfnc.o: ${CDIR}.c
	${CC} ${CFLAGS} ${CDIR}.c -c -fPIC -o $@
	
bin/cinit.o: ${CDIR}_initialize.c
	${CC} ${CFLAGS} ${CDIR}_initialize.c -c -fPIC -o $@

bin/cterm.o: ${CDIR}.c
	${CC} ${CFLAGS} ${CDIR}_terminate.c -c -fPIC -o $@
	
bin/gfnc.o: ${GDIR}.c
	${CC} ${CFLAGS} ${GDIR}.c -c -fPIC -o $@
	
bin/ginit.o: ${GDIR}_initialize.c
	${CC} ${CFLAGS} ${GDIR}_initialize.c -c -fPIC -o $@

bin/gterm.o: ${GDIR}.c
	${CC} ${CFLAGS} ${GDIR}_terminate.c -c -fPIC -o $@
	
bin/ffnc.o: ${FDIR}.c
	${CC} ${CFLAGS} ${FDIR}.c -c -fPIC -o $@
	
bin/finit.o: ${FDIR}_initialize.c
	${CC} ${CFLAGS} ${FDIR}_initialize.c -c -fPIC -o $@

bin/fterm.o: ${FDIR}.c
	${CC} ${CFLAGS} ${FDIR}_terminate.c -c -fPIC -o $@



