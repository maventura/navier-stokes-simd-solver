#include <iostream>
#include <stdio.h>
#include <string.h>

using namespace std;

int main(){
	bool stress_ram = false;
	bool stress_cpu = true;

	if(stress_ram){
		while(1){
			char str[10000];
			memset (str,'-',10000);
		}
	}

	if(stress_cpu){
		float a = 500.0;
		while(1){
			a = a*1.05 + 0.1;
			if(a > 10000.0){
				a = 500.0;
			}
		}
	}

}