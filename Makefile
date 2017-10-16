#CPP = g++ 
#FLAGS = -std=c++11

.PHONY: clean

main: main.cpp mat3.h simulator.h io.h
	g++ -std=c++11 main.cpp -o main

clean:
	rm -f main
	rm -f asm
	rm ./out/*

help: 
	@echo 'For compiling the c++ version, run make.'
	@echo 'For compiling the assembly version, run make assembly.'
	@echo 'To clean the executable file, run clean.'
	@echo 'That is all for now...'

assembly: main.cpp vvp_asm.asm mat3.h simulator.h io.h
	nasm -f elf64 -g -F dwarf vvp_asm.asm -o vvp_asm.o
	g++ -c -m64 -std=c++11 main.cpp -o main.o
	g++ -o main -m64 vvp_asm.o main.o