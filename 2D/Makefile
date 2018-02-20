#CPP = g++ 
#FLAGS = -std=c++11

.PHONY: clean

asm: main.cpp main.asm mat2.h simulator.h io.h #assembly
	nasm -f elf64 -g -F dwarf main.asm -o mainasm.o
	g++ -DUSE_ASM -c -m64 -std=c++11 main.cpp -o main.o
	g++ -o main_asm -m64 mainasm.o main.o

cpp: main.cpp mat2.h io.h
	g++ -DUSE_CPP main.cpp -std=c++11 -c
	g++ -o main_cpp main.o



cpp_o1: main.cpp mat2.h io.h
	g++ -DUSE_CPP main.cpp -O1 -std=c++11 -c
	g++ -o main_cpp_o1 main.o


cpp_o2: main.cpp mat2.h io.h
	g++ -DUSE_CPP main.cpp -O2 -std=c++11 -c
	g++ -o main_cpp_o2 main.o


cpp_o3: main.cpp mat2.h io.h
	g++ -DUSE_CPP main.cpp -O3 -std=c++11 -c
	g++ -o main_cpp_o3 main.o


cpp_ofast: main.cpp mat2.h io.h
	g++ -DUSE_CPP main.cpp -Ofast -std=c++11 -c
	g++ -o main_cpp_ofast main.o

asm_omp: main.cpp main.asm mat2.h simulator.h io.h #assembly
	nasm -f elf64 -g -F dwarf main.asm -o mainasm.o
	g++ -DUSE_ASM -c -m64 -std=c++11 main.cpp -o main.o
	g++ -fopenmp -o main_asm_omp -m64 mainasm.o main.o

cpp_omp:  main.cpp mat2.h io.h
	g++ -DUSE_CPP -fopenmp main.cpp -std=c++11  -o main_cpp_omp

icc: main.cpp mat2.h io.h
	/opt/intel/compilers_and_libraries/linux/bin/intel64/icpc -DUSE_CPP -DUSE_ICC  -fp_speculation=safe main.cpp -std=c++11 -c
	/opt/intel/compilers_and_libraries/linux/bin/intel64/icpc -o main_icc main.o

clean:
	rm -f main_asm
	rm -f main_cpp
	rm -f main_cpp_o1
	rm -f main_cpp_o2
	rm -f main_cpp_o3
	rm -f main_cpp_ofast
	rm -f main_asm_omp
	rm -f main_cpp_omp
	rm -f *.o
	rm -f ./out/*
