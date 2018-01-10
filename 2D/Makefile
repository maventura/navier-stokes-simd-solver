#CPP = g++ 
#FLAGS = -std=c++11

.PHONY: clean

asm: main.cpp main.asm mat2.h simulator.h io.h #assembly
	nasm -f elf64 -g -F dwarf main.asm -o mainasm.o
	g++ -DUSE_ASM -c -m64 -std=c++11 main.cpp -o main.o
	g++ -o main -m64 mainasm.o main.o

cpp:  main.cpp mat2.h io.h
	g++ main.cpp -std=c++11  -o main

clean:
	rm -f main
	rm -f *.o
	rm ./out/*
