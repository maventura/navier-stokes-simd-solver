#CPP = g++ 
#FLAGS = -std=c++11

.PHONY: clean

main: main.cpp mat3.h simulator.h io.h
	g++ -std=c++11 main.cpp -o main

clean:
	rm -f main
	rm -f assmain
	rm ./out/*

help: 
	@echo 'For principal cpp executable, run make.'
	@echo 'For compiling the assembly version, run make assembly.'
	@echo 'To clean the executable file, run clean.'
	@echo 'That is all for now...'

assembly: main.cpp vvp_asm.asm mat3.h simulator.h io.h
	nasm -f elf64 vvp_asm.asm -o vvp_asm.o
	g++ -o assmain -std=c++11 main.cpp vvp_asm.o