#include <iostream>
#include <vector>
#include <chrono>  // for high_resolution_clock.
#include <sstream>
#include <cmath>
#include "simulator.h"

using namespace std;


int main(int argc, char const *argv[]) {
  simulator sim;
  sim.process();
  return 0;
}
