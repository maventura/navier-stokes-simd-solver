#CPP = g++ 
#FLAGS = -std=c++11

.PHONY: clean

asm: main.cpp vvp_asm.asm mat3.h simulator.h io.h #assembly
	nasm -f elf64 -g -F dwarf vvp_asm.asm -o vvp_asm.o
	g++ -DUSE_ASM -c -m64 -std=c++11 main.cpp -o main.o
	g++ -o main -m64 vvp_asm.o main.o

cpp:  main.cpp mat3.h simulator.h io.h #cpp
	g++ -DUSE_CPP -std=c++11 main.cpp -o main

clean:
	rm -f main
	rm -f *.o
	rm ./out/*
