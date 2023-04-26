########################################################################
#                                                                      #
# Spectral element method code makefile script                         #
#                                                                      #
# Written by: Liang Zhao, IGG CAS                                      #
# Begin time: Nov 19, 2007.                                            #
# Last modified: 20140321: using openMP for sem2d.c                    #
# 20140928: adding SH system                                           #
########################################################################

RECIPES = Recipes_C/
CC   = gcc  
CC   = cc 
FLAG = -fopenmp

LIBS = -lm -Lgetpar/lib -lget

GETPARINC  = -Igetpar/include
F77  = f77
F95  = f95  
BIN  = ../bin/

all: sem2d sem2d_psvsh sem2dMP asersh_sem aserpsv_sem demult demult_sh clean

sem2d:  \
	sem2d.o readmodel.o interpl.o lagrange.o 
	$(CC) -o $(BIN)sem2d \
	sem2d.o readmodel.o interpl.o lagrange.o $(LIBS) 

sem2d.o: \
	sem2d.c sem2d.h
	$(CC) -c sem2d.c

sem2d_psvsh:  \
	sem2d_psvsh.o readmodel.o interpl.o lagrange.o 
	$(CC) -o $(BIN)sem2d_psvsh \
	sem2d_psvsh.o readmodel.o interpl.o lagrange.o $(LIBS) 

sem2d_psvsh.o: \
	sem2d_psvsh.c sem2d.h
	$(CC) -c sem2d_psvsh.c

sem2dMP:  \
	sem2dMP.o readmodel.o interpl.o lagrange.o 
	$(CC) $(FLAG) -o $(BIN)sem2dMP \
	sem2dMP.o readmodel.o interpl.o lagrange.o $(LIBS) 

sem2dMP.o: \
	sem2dMP.c sem2d.h
	$(CC) -c $(FLAG) sem2dMP.c


aserpsv_sem: \
	aserpsv_sem.o aser_psv.o source.o readmodel.o aser.rout_noplot.o dcsource.o nrutil.o lagrange.o
	$(F77) -o $(BIN)aserpsv_sem aserpsv_sem.o source.o readmodel.o aser_psv.o \
	aser.rout_noplot.o dcsource.o nrutil.o lagrange.o $(LIBS)                   

asersh_sem: \
	asersh_sem.o aser_sh.o source.o readmodel.o aser.rout_noplot.o dcsource.o nrutil.o lagrange.o
	$(F77) -o $(BIN)asersh_sem asersh_sem.o source.o readmodel.o aser_sh.o \
	aser.rout_noplot.o dcsource.o nrutil.o lagrange.o $(LIBS)                   

aser_psv.o: \
	aser_psv.f psv.p
	$(F77) -c aser_psv.f
# aser_sh.f is copied from SHHY/aser_psv.f, 20141025
aser_sh.o: \
	aser_sh.f psv.p
	$(F77) -c aser_sh.f

aser.rout_noplot.o: \
	aser.rout_noplot.f psv.p
	$(F77) -c aser.rout_noplot.f

dcsource.o: \
	dcsource.f
	$(F77) -c dcsource.f

interpl.o: \
	interpl.c 
	$(CC) -c interpl.c	

nrutil.o: \
	$(RECIPES)nrutil.c 
	$(CC) -c  -I$(RECIPES) $(RECIPES)nrutil.c -o nrutil.o

demult: \
	demult.o readmodel.o
	$(CC) -o $(BIN)demult \
	demult.o readmodel.o $(LIBS) 

demult.o: \
	demult.c
	$(CC) -c demult.c

demult_sh: \
	demult_sh.o readmodel.o
	$(CC) -o $(BIN)demult_sh \
	demult_sh.o readmodel.o $(LIBS) 

demult_sh.o: \
	demult_sh.c
	$(CC) -c demult_sh.c

readmodel.o: \
	readmodel.c 
	$(CC) -c readmodel.c

lagrange.o: \
	lagrange.c 
	$(CC) -c lagrange.c	

source.o: \
	source.f
	$(F77) -c source.f		

clean: 
	rm -f \
	*.o

########################################################################
