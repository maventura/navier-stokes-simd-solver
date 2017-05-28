#include <iostream>
#include <vector>

using namespace std;
using mat3 = vector<vector<vector<long double > > >;
using mat2 = vector<vector<long double > >;
using mat1 = vector<long double>;

int main(){
	
	//1D matrix, a line.
	mat1 a1;
	a1.push_back(1);
	a1.push_back(2);
	a1.push_back(3);

	//2D matrix, a square.
	mat2 a2;
	a2.push_back(a1);
	a2.push_back(a1);
	a2.push_back(a1);

	//3D matrix, a cube.
	mat3 a3;
	a3.push_back(a2);
	a3.push_back(a2);
	a3.push_back(a2);

	cout << a3[2][1][2] << endl;
}
