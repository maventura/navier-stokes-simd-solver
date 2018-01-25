#CPP = g++ 
#FLAGS = -std=c++11

.PHONY: clean

asm: main.cpp main.asm mat2.h simulator.h io.h #assembly
	nasm -f elf64 -g -F dwarf main.asm -o mainasm.o
	g++ -DUSE_ASM -c -m64 -std=c++11 main.cpp -o main.o
	g++ -o mainASM -m64 mainasm.o main.o

cpp:  main.cpp mat2.h io.h
	g++ -DUSE_CPP main.cpp -std=c++11  -o mainCPP

asm_omp: main.cpp main.asm mat2.h simulator.h io.h #assembly
	nasm -f elf64 -g -F dwarf main.asm -o mainasm.o
	g++ -DUSE_ASM -c -m64 -std=c++11 main.cpp -o main.o
	g++ -fopenmp -o mainASM_omp -m64 mainasm.o main.o

cpp_omp:  main.cpp mat2.h io.h
	g++ -DUSE_CPP -fopenmp main.cpp -std=c++11  -o mainCPP_omp

clean:
	rm -f mainASM
	rm -f mainCPP
	rm -f *.o
	rm ./out/*
