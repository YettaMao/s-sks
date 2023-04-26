########################################################################
#                                                                      #
# Model code  makefile script                                          #
#      for generating elastic constant matrix                          #
# Written by: Liang Zhao, iGG, CAS                                     #
# Last modified: July. 4, 2006.                                        #
#                                                                      #
########################################################################

CC   = gcc  
CC   = cc 

LIBS = -lm -L/home/zhl/work/aniso/code/plot/getpar/lib
GETPARINC  = -I/home/zhl/work/aniso/code/plot/getpar/include
F77  = g77 

all: elastic_rotate elastic_dip cjcbi0 clean

elastic_rotate:  \
	elastic_rotate.o
	$(CC) -o $(BIN)elastic_rotate \
	elastic_rotate.o $(LIBS) 

elastic_rotate.o: \
	elastic_rotate.c 
	$(CC) -c elastic_rotate.c

elastic_dip:  \
	elastic_dip.o
	$(CC) -o $(BIN)elastic_dip \
	elastic_dip.o $(LIBS) 

elastic_dip.o: \
	elastic_dip.c 
	$(CC) -c elastic_dip.c

cjcbi0:  \
	cjcbi0.o
	$(CC) -o $(BIN)cjcbi0 \
	cjcbi0.o  $(LIBS) 

cjcbi0.o: \
	cjcbi0.c 
	$(CC) -c cjcbi0.c		

clean: 
	rm -f \
	*.o

########################################################################
