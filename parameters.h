#ifndef PARAMETERS_H
#define PARAMETERS_H


#include <math.h>
#include "io.h"

using namespace std;

double xMax;
double yMax;
double zMax;
double tMax;
double nu;
double rho;
double dx;
double dy;
double dz;
double dt;
int nX = round(xMax / dx) + 1;
int nY = round(yMax / dy) + 1;
int nZ = round(yMax / dy) + 1;
int nT = round(tMax / dt) + 1;double al;
double fixedPointError;
double minFixedPointIters;
int printPercentageSteps;
double pi = atan(1)*4;
double C_d;
//TODO: area should depende upon intersection
// of the element with the fan.
int maxSteps;
int stepsUntilPrint;
bool printPercentage;

void readParameters(string file_name){
  io file(file_name, io::type_read);
  string tag;
  while(file.readWord(tag)){

    if(tag == "xMax") file.readDouble(xMax);
    if(tag == "yMax") file.readDouble(yMax);
    if(tag == "zMax") file.readDouble(zMax);
    if(tag == "tMax") file.readDouble(tMax);
    if(tag == "nu") file.readDouble(nu);
    if(tag == "rho") file.readDouble(rho);
    if(tag == "dx") file.readDouble(dx);
    if(tag == "dy") file.readDouble(dy);
    if(tag == "dz") file.readDouble(dy);
    if(tag == "dt") file.readDouble(dt);
    if(tag == "al") file.readDouble(al);
    if(tag == "fixedPointError") file.readDouble(fixedPointError);
    if(tag == "minFixedPointIters") file.readDouble(minFixedPointIters);
    if(tag == "printPercentageSteps") file.readInt(printPercentageSteps);
    if(tag == "C_d") file.readDouble(C_d);
    if(tag == "maxSteps") file.readInt(maxSteps);
    if(tag == "stepsUntilPrint") file.readInt(stepsUntilPrint);
    if(tag == "printPercentage") file.readBool(printPercentage);
  }
  file.close();
}

#endif
