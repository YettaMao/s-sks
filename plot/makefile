########################################################################
#                                                                      #
# Plot code  makefile script                                           #
#      for plotting the results from aserpsv_sem                       #
# Written by: Liang Zhao, iGG, CAS                                     #
# Last modified: Feb 1, 2008.                                          #
#                                                                      #
########################################################################

RECIPES = Recipes_C/
CC   = gcc  
CC   = cc 
#-xdepend  -xchip=ultra -xarch=generic  
CC   = cc 
LIBS = -lm -Lgetpar/lib -lget
GETPARINC  = -Igetpar/include
F77  = g77
F77  = f77  

all: read_grt read_grt_sh read_uvw read_sac readsac_seg sac_plot sac_plot2 sac_plot_f uvw2sac clean


read_uvw:  \
	read_uvw.o readmodel.o
	$(CC) -o $(BIN)read_uvw \
	read_uvw.o readmodel.o $(LIBS) 

read_uvw.o: \
	read_uvw.c 
	$(CC) -c read_uvw.c	

uvw2sac:  \
	uvw2sac.o sac_rw.o source.o readmodel.o
	$(F77) -o $(BIN)uvw2sac \
	uvw2sac.o sac_rw.o source.o readmodel.o $(LIBS) 

uvw2sac.o: \
	uvw2sac.c sachead.h
	$(CC) -c uvw2sac.c

read_sac:  \
	read_sac.o sac_rw.o
	$(F77) -o $(BIN) read_sac \
	read_sac.o sac_rw.o $(LIBS) 

read_sac.o: \
	read_sac.c sachead.h
	$(CC) -c read_sac.c

readsac_seg:  \
	readsac_seg.o sac_rw.o
	$(F77) -o $(BIN) readsac_seg \
	readsac_seg.o sac_rw.o $(LIBS) 

readsac_seg.o: \
	readsac_seg.c sachead.h
	$(CC) -c readsac_seg.c

misfit:  \
	misfit.o sac_rw.o
	$(F77) -o $(BIN) misfit \
	misfit.o sac_rw.o $(LIBS) 

misfit.o: \
	misfit.c sachead.h
	$(CC) -c misfit.c

sac_plot:  \
	sac_plot.o sac_rw.o
	$(CC) -o $(BIN)sac_plot \
	sac_plot.o sac_rw.o $(LIBS) 

sac_plot.o: \
	sac_plot.c sachead.h
	$(CC) -c sac_plot.c		

sac_plot2:  \
	sac_plot2.o sac_rw.o
	$(CC) -o $(BIN)sac_plot2 \
	sac_plot2.o sac_rw.o $(LIBS) 

sac_plot2.o: \
	sac_plot2.c sachead.h
	$(CC) -c sac_plot2.c

sac_plot_f:  \
	sac_plot_f.o sac_rw.o
	$(CC) -o $(BIN)sac_plot_f \
	sac_plot_f.o sac_rw.o $(LIBS) 

sac_plot_f.o: \
	sac_plot_f.c sachead.h
	$(CC) -c sac_plot_f.c	

read_grt:  \
	read_grt.o readmodel.o
	$(CC) -o $(BIN)read_grt \
	read_grt.o readmodel.o $(LIBS) 

read_grt.o: \
	read_grt.c 
	$(CC) -c read_grt.c	

read_grt_sh:  \
	read_grt_sh.o readmodel.o
	$(CC) -o $(BIN)read_grt_sh \
	read_grt_sh.o readmodel.o $(LIBS) 

read_grt_sh.o: \
	read_grt_sh.c 
	$(CC) -c read_grt_sh.c	

interpl.o: \
	interpl.c 
	$(CC) -c interpl.c	
	
nrutil.o: \
	nrutil.c 
	$(CC) -c nrutil.c	

readmodel.o: \
	readmodel.c 
	$(CC) -c readmodel.c	

source.o: \
	source.f
	$(F77) -c source.f	

sac_rw.o: \
	sac_rw.c sachead.h 
	$(CC) -c sac_rw.c

clean: 
	rm -f \
	*.o

########################################################################
