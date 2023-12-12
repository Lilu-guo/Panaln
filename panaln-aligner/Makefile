CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O3 -std=c99 -fopenmp -no-pie
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DHAVE_PTHREAD
DBGFLAGS= 	#-DDEBUG_ENABLED
OBJS=		align.o exact_match.o inexact_match.o main.o io.o FMapi.o
PROG=		clean panaln
INCLUDES=	
LIBS=		-lm -lz -lpthread
DFLAGS +=   -DFM
# DFLAGS +=   -DmyPrint

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $(DBGFLAGS) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)
		
panaln:$(OBJS)
		cp ./FM/fm.a .
		$(CXX) $(CFLAGS) $(DFLAGS) $(OBJS) fm.a -o $@ $(LIBS)
fm.a:
		cp ./FM/fm.a .
clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a
