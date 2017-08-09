#CPP = g++ 
#FLAGS = -std=c++11

.PHONY: clean

main : main.cpp mat3.h simulator.h io.h
	g++ -std=c++11 main.cpp -o main
	mkdir out

clean:
	rm -f main
	rm -r -f out

help: 
	@echo 'For principal cpp executable, run make.'
	@echo 'That is all for now'
