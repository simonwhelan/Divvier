CC=gcc
CPP=g++
OPTIMISER = -O3
CPPFLAGS = $(OPTIMISER) -std=c++14 -Wall -Wmissing-prototypes -Wshadow -fmessage-length=0 -msse2 -mfpmath=sse
CFLAGS = $(OPTIMISER) -Wall -Wmissing-prototypes -Wshadow 
#CFLAGS = -g 

# C Code
LIB = -lm 
OBJS = trim.o utils.o hmm.o matrices.o
HDR = trim.h utils.h hmm.h matrices.h
INC = -I/usr/local/include
# CPP code
CPPOBJS = Cluster.o Tree.o Random.o Divvier.o bionj.o
CPPHDR = Cluster.h Tree.h Random.h Divvier.h

all : fastZorro

# Zorro
trim.o : trim.c trim.h
	$(CC) $(CFLAGS) $(INC) -c trim.c

utils.o : utils.c utils.h
	$(CC) $(CFLAGS) $(INC) -c utils.c

hmm.o : hmm.c hmm.h
	$(CC) $(CFLAGS) $(INC) -c hmm.c

matrices.o : matrices.c matrices.h
	$(CC) $(CFLAGS) $(INC) -c matrices.c

fastZorro : $(CPPOBJS) $(OBJS) $(HDR)
	$(CPP) $(CPPFLAGS) $(INC) $(LIB) -o divvier $(OBJS) $(CPPOBJS)

#Cluster
Cluster.o : Cluster.cpp Cluster.h
	$(CPP) $(CPPFLAGS) $(INC) -c Cluster.cpp

Tree.o : Tree.cpp Tree.h
	$(CPP) $(CPPFLAGS) $(INC) -c Tree.cpp

Random.o : Random.cpp Random.h
	$(CPP) $(CPPFLAGS) $(INC) -c Random.cpp

#bionj
bionj.o : bionj.cxx
	$(CPP) $(CPPFLAGS) $(INC) -c bionj.cxx

# Divvier
Divvier.o : Divvier.cpp Divvier.h
	$(CPP) $(CPPFLAGS) $(INC) -c Divvier.cpp

clean:
	rm -f $(OBJS) $(CPPOBJS)
	rm -f divvier
	rm -f core
	rm -f *~
	rm -f a.out

