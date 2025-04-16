CC=			gcc
CXX=		g++
CFLAGS=		-g -w -O3 -std=c99 -fopenmp -no-pie -D_GNU_SOURCE  
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DHAVE_PTHREAD
DBGFLAGS= 	#-DDEBUG_ENABLED
OBJS=		bwt.o align.o exact_match.o inexact_match.o is.o main.o io.o FMapi.o comb.o data_prep.o
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
panaln:$(OBJS) fm.a
		cp ./FM/fm.a .
		$(CXX) $(CFLAGS) $(DFLAGS) $(OBJS) fm.a -o $@ $(LIBS) -L ./semiWFA/lib -lwfa
fm.a:
		cp ./FM/fm.a .
clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a
