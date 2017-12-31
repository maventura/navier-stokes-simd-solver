#include <stdio.h>
#include <string.h>
#include <chrono>       // for high_resolution_clock.
#include <math.h>
#include <vector>
#include <utility>      // std::pair
#include <sstream>
#include "mat2.h"
#include "io.h"
#include "simulator.h"

using namespace std;


float distance(mat2 &A, mat2 &B) {
    float sum = 0;
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            sum += pow(B.at(i, j) - A.at(i, j), 2);
        }
    }
    return sqrt(sum);
}


int main() {

    simulator sim;
    sim.process();
}
