CC=gcc

CFLAGS = -O3 -Wall -Wmissing-prototypes -Wshadow 
#CFLAGS = -g 


#LIB = -lm -lmysqlclient -L/panfs/panfs.gm.berkeley.edu/download/souravc/local/lib/mysql/
LIB = -lm 
OBJS = trim.o utils.o hmm.o matrices.o
HDR = trim.h utils.h hmm.h matrices.h
INC = -I/usr/local/include


all : probmask

.c.o: $(HDR)
	$(CC) $(CFLAGS) $(INC) -c $*.c

probmask : $(OBJS) $(HDR)
	$(CC) $(CFLAGS) $(INC) $(LIB) -o probmask $(OBJS)

clean:
	rm -f $(OBJS)
	rm -f probmask
	rm -f core
	rm -f *~
	rm -f a.out

