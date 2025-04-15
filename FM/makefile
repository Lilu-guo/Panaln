CXX=g++
# CFLAGS=  -O3 -no-pie -g -Wall -mpopcnt #配合popcnt1中的__builtin_popcountll使用,但是并没有明显加速
CFLAGS=  -O3 -no-pie -g -Wall # 正式发布用
# CFLAGS=  -O -no-pie -g -Wall # 调试用

fm:div main.o fm.a
	$(CXX) -o fm main.o fm.a

all: fm clean
#adaptive: div man.o fm.a
#	g++ test.o fm.a $(CFLAGS) -o adaptive -lrt

# 打包目标文件，静态库
fm.a:ABS_WT.o Balance_WT.o Huffman_WT.o Hutacker_WT.o FM.o newDmap.o BitMap.o UseCount.o WT_Handle.o InArray.o loadkit.o savekit.o divsufsort.o sssort.o trsort.o utils.o
	ar rc fm.a ABS_WT.o Balance_WT.o Huffman_WT.o Hutacker_WT.o FM.o newDmap.o BitMap.o UseCount.o WT_Handle.o  InArray.o loadkit.o savekit.o divsufsort.o sssort.o trsort.o utils.o

%.o:%.cpp *.h
	$(CXX) -c  $(CFLAGS) $< -o $@

main.o:main.cpp  FM.h
	$(CXX) -c  main.cpp

test.o: test.cpp FM.h
	$(CXX) -c $(CFLAGS) test.cpp -lrt

clean:
	rm *.a *.o ./divsufsort/*.a ./divsufsort/*.o

div:
	make -C ./divsufsort/; cp divsufsort/libdivsufsort.a .; ar x libdivsufsort.a; rm libdivsufsort.a 
