CC=gcc  
CFLAGS=  -O3 -Wall -no-pie -fno-pie -fPIC -DHAVE_CONFIG_H 
libdivsufsort.a:divsufsort.o sssort.o trsort.o utils.o
	ar rc libdivsufsort.a divsufsort.o sssort.o trsort.o utils.o

%.o:%.cpp *.h
	$(CC) -c $(CFLAGS) $< -o $@
clean:
	rm *.a *.o
