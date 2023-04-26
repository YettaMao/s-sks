########################################################################
#                                                                      #
# makefile script for generating rays                                  #
# this version is for rf                                               #
# Written by: Liang Zhao, IGG, CAS                                     #
# Last modified: 20131113                                              #
#                                                                      #
########################################################################

CC   = gcc  
CC   = cc 

F77  = gfortran 
BIN  = ./

all: ray_zl ray_sh_zl ray_p ray_new ray_ref clean

ray_p: \
	ray_p.o 
	$(F77) -o $(BIN)ray_p \
	ray_p.o 

ray_p.o: \
	ray_p.f
	$(F77) -c ray_p.f

ray_zl: \
	ray_zl.o 
	$(F77) -o $(BIN)ray_zl \
	ray_zl.o 

ray_zl.o: \
	ray_zl.f
	$(F77) -c ray_zl.f

ray_sh_zl: \
	ray_sh_zl.o 
	$(F77) -o $(BIN)ray_sh_zl \
	ray_sh_zl.o 

ray_sh_zl.o: \
	ray_sh_zl.f
	$(F77) -c ray_sh_zl.f


ray_ref: \
	ray_ref.o 
	$(F77) -o $(BIN)ray_ref \
	ray_ref.o 

ray_ref.o: \
	ray_ref.f
	$(F77) -c ray_ref.f

ray_new: \
	ray_new.o 
	$(F77) -o $(BIN)ray_new \
	ray_new.o 

ray_new.o: \
	ray_new.f
	$(F77) -c ray_new.f
clean: 
	rm -f \
	*.o

########################################################################
