#CPP = g++ 
#FLAGS = -std=c++11

.PHONY: clean
all: asm asm_o1 asm_o2 asm_o3 asm_ofast cpp cpp_o1 cpp_o2 cpp_o3 cpp_ofast asm_omp cpp_omp cpp_omp_o1 cpp_omp_o2 cpp_omp_o3 cpp_omp_ofast icc icc_o1 icc_o2 icc_o3 icc_ofast


asm: main.cpp main.asm mat2.h simulator.h io.h #assembly
	nasm -f elf64 -g -F dwarf main.asm -o mainasm.o
	g++ -DUSE_ASM -c -m64 -std=c++11 main.cpp -o main.o
	g++ -o main_asm -m64 mainasm.o main.o


asm_o1: main.cpp main.asm mat2.h simulator.h io.h #assembly
	nasm -f elf64 -g -F dwarf main.asm -o mainasm.o
	g++ -DUSE_ASM -c -m64 -std=c++11 -O1 main.cpp -o main.o
	g++ -o main_asm_o1 -m64 mainasm.o main.o

asm_o2: main.cpp main.asm mat2.h simulator.h io.h #assembly
	nasm -f elf64 -g -F dwarf main.asm -o mainasm.o
	g++ -DUSE_ASM -c -m64 -std=c++11 -O2 main.cpp -o main.o
	g++ -o main_asm_o2 -m64 mainasm.o main.o


asm_o3: main.cpp main.asm mat2.h simulator.h io.h #assembly
	nasm -f elf64 -g -F dwarf main.asm -o mainasm.o
	g++ -DUSE_ASM -c -m64 -std=c++11 -O3 main.cpp -o main.o
	g++ -o main_asm_o3 -m64 mainasm.o main.o


asm_ofast: main.cpp main.asm mat2.h simulator.h io.h #assembly
	nasm -f elf64 -g -F dwarf main.asm -o mainasm.o
	g++ -DUSE_ASM -c -m64 -std=c++11 -Ofast main.cpp -o main.o
	g++ -o main_asm_ofast -m64 mainasm.o main.o


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

cpp_omp_o1:  main.cpp mat2.h io.h
	g++ -DUSE_CPP -fopenmp main.cpp -std=c++11 -O1 -o main_cpp_omp_o1

cpp_omp_o2:  main.cpp mat2.h io.h
	g++ -DUSE_CPP -fopenmp main.cpp -std=c++11 -O2 -o main_cpp_omp_o2

cpp_omp_o3:  main.cpp mat2.h io.h
	g++ -DUSE_CPP -fopenmp main.cpp -std=c++11 -O3 -o main_cpp_omp_o3

cpp_omp_ofast:  main.cpp mat2.h io.h
	g++ -DUSE_CPP -fopenmp main.cpp -std=c++11 -Ofast -o main_cpp_omp_ofast

icc: main.cpp mat2.h io.h
	/opt/intel/compilers_and_libraries/linux/bin/intel64/icpc -DUSE_CPP -DUSE_ICC main.cpp -std=c++11 -c
	/opt/intel/compilers_and_libraries/linux/bin/intel64/icpc -o main_icc main.o

icc_o1: main.cpp mat2.h io.h
	/opt/intel/compilers_and_libraries/linux/bin/intel64/icpc -DUSE_CPP -DUSE_ICC main.cpp -std=c++11 -O1 -c
	/opt/intel/compilers_and_libraries/linux/bin/intel64/icpc -o main_icc_o1 main.o

icc_o2: main.cpp mat2.h io.h
	/opt/intel/compilers_and_libraries/linux/bin/intel64/icpc -DUSE_CPP -DUSE_ICC main.cpp -std=c++11 -O2 -c
	/opt/intel/compilers_and_libraries/linux/bin/intel64/icpc -o main_icc_o2 main.o

icc_o3: main.cpp mat2.h io.h
	/opt/intel/compilers_and_libraries/linux/bin/intel64/icpc -DUSE_CPP -DUSE_ICC main.cpp -std=c++11 -O3 -c
	/opt/intel/compilers_and_libraries/linux/bin/intel64/icpc -o main_icc_o3 main.o

icc_ofast: main.cpp mat2.h io.h
	/opt/intel/compilers_and_libraries/linux/bin/intel64/icpc -DUSE_CPP -DUSE_ICC main.cpp -std=c++11 -Ofast -c
	/opt/intel/compilers_and_libraries/linux/bin/intel64/icpc -o main_icc_ofast main.o


objdump: main.o
	objdump -D -M intel -S main.o>dumpFile.txt

help:
	@echo 'Opciones para el Makefile:'
	@echo ' Make cpp, para compilación normal con g++'
	@echo ' Make cpp_o1, para compilación con optimizaciones -O1'
	@echo ' Make cpp_o2, para compilación con optimizaciones -O2'
	@echo ' Make cpp_o3, para compilación con optimizaciones -O3'
	@echo ' Make cpp_ofast, para compilación con optimizaciones -Ofast'
	@echo ' Make cpp_omp, para utilizar soporte de OpenMP'
	@echo ' Make asm, para compilar la versión assembler'
	@echo ' Make objdump, para generar código assembler de la versión cpp compilada.'

clean:
	rm -f main_asm
	rm -f main_asm_o1
	rm -f main_asm_o2
	rm -f main_asm_o3
	rm -f main_asm_ofast
	rm -f main_cpp
	rm -f main_cpp_o1
	rm -f main_cpp_o2
	rm -f main_cpp_o3
	rm -f main_cpp_ofast
	rm -f main_asm_omp
	rm -f main_cpp_omp
	rm -f main_cpp_omp_o1
	rm -f main_cpp_omp_o2
	rm -f main_cpp_omp_o3
	rm -f main_cpp_omp_ofast
	rm -f main_icc
	rm -f main_icc_o1
	rm -f main_icc_o2
	rm -f main_icc_o3
	rm -f main_icc_ofast
	rm -f *.o
	rm -f ./out/*
