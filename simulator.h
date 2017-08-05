#ifndef simulator_H
#define simulator_H

#include <math.h>
#include "mat3.h"
#include "io.h"

class simulator{

public:
  simulator();
  void process();

private:
  double xMax,  yMax,  zMax,  tMax;
  double nu, rho, C_d;
  double dx,  dy,  dz,  dt;
  int nX, nY, nZ, nT;

  double al, fixedPointError, minFixedPointIters;
  double pi = atan(1)*4;
  int stepsUntilPrint, printPercentageSteps, maxSteps;
  bool printPercentage;
  double t;

  //Calculation variables,
  double U1x, U2x, U1y, U2y, U1z, U2z, U1xx, U2xx, U1yy, U2yy, U1zz, U2zz;
  double P1x, P2x, P1y, P2y, P1z, P2z, V1x, V2x, V1y, V2y, V1z, V2z;
  double V1xx, V2xx, V1yy, V2yy, V1zz, V2zz, P1xx, P1yy, P1zz, P2xx, P2yy, P2zz;
  mat3 U0, V0, W0, P0;
  mat3 U1, V1, W1, P1;
  mat3 U2, V2, W2, P2;

  void calcTerms(int i, int j, int k);
  void saveVtk(mat3 state, string file_name);
  void setBorderConditions();
  void readParameters(string file_name);
};

simulator::simulator(){

  readParameters("parameters.txt");

  nX = round(xMax / dx) + 1;
  nY = round(yMax / dy) + 1;
  nZ = round(yMax / dy) + 1;
  nT = round(tMax / dt) + 1;double al;

  U0.confAndInit(nX, nY, nZ, 0);
  V0.confAndInit(nX, nY, nZ, 0);
  W0.confAndInit(nX, nY, nZ, 0);
  P0.confAndInit(nX, nY, nZ, 0);

  U1.confAndInit(nX, nY, nZ, 0);
  V1.confAndInit(nX, nY, nZ, 0);
  W1.confAndInit(nX, nY, nZ, 0);
  P1.confAndInit(nX, nY, nZ, 0);

  U2.confAndInit(nX, nY, nZ, 0);
  V2.confAndInit(nX, nY, nZ, 0);
  W2.confAndInit(nX, nY, nZ, 0);
  P2.confAndInit(nX, nY, nZ, 0);

  setBorderConditions();
}

void simulator::readParameters(string file_name){
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

void simulator::setBorderConditions(){
      //TODO: Add if(cavityFlow) here and in the parameters.
      for (int i = 0; i < nX; ++i) {
          for (int k = 0; k < nZ; ++k) {
              U0.set(i, nY - 1, k, 1.0);
              U1.set(i, nY - 1, k, 1.0);
              U2.set(i, nY - 1, k, 1.0);
          }
      }
}

void simulator::process(){
  for (t = 0; t < tMax; t = t + dt) {
    U0.print();
    V0.print();
    W0.print();
    for (int i = 1; i < nX - 1; ++i) {
      for (int j = 1; j < nY - 1; ++j) {
        for (int k = 0; k < nZ; ++k) {
          for (int iter = 0; iter < 10; ++iter) { //TODO: Iter termination condition

            calcTerms(i,j,k);

            U2.set(i,j,k,U2.at(i,j,k)+1);
            V2.set(i,j,k,V2.at(i,j,k)+10);
            W2.set(i,j,k,W2.at(i,j,k)+100);
            P2.set(i,j,k,P2.at(i,j,k)+1000);

            // double diff = sqrt( pow(U3.at(i,j,k) - oldU, 2) + pow(V3.at(i,j,k) - oldV, 2) + pow(W3.at(i,j,k) - oldW, 2 ) );
            // if (diff < 0.01) {
            //     break;
            // } else if (iter > 50) {
            //     cerr << "WARNING: unstable." << endl;
            //     break;
            // }
          }
        }
      }
    }
    //TODO: the 3 matrices exist so that we work on them without changing old values of the 2 matrices.
    // But they represent the new 2 matrices. Change name.

    U0.setAll(U1);
    V0.setAll(V1);
    W0.setAll(W1);
    P0.setAll(P1);

    U1.setAll(U2);
    V1.setAll(V2);
    W1.setAll(W2);
    P1.setAll(P2);
  }
  cerr << "Message: Processing finished correctly" << endl << flush;
}


void simulator::saveVtk(mat3 m, string file_name){
// Save a 3-D scalar array in VTK format.

  io out(file_name, io::type_write);
  out.write("# vtk DataFile Version 2.0");
  out.newLine();
  out.write("ASCII");
  out.newLine();
  out.write("DATASET STRUCTURED_POINTS");
  out.newLine();
  out.write("DIMENSIONS    " + to_string(nX) + " " +  to_string(nY) + " " + to_string(nZ));
  out.newLine();
  out.write("ORIGIN    0.000   0.000   0.000");
  out.newLine();
  out.write("SPACING    1.000   1.000   1.000");
  out.newLine();
  out.write("POINT_DATA   " + to_string(nX*nY));
  out.newLine();
  out.write("SCALARS scalars float");
  out.newLine();
  out.write("LOOKUP_TABLE default");
  out.newLine();

  for(int a = 1; a < nX; a++){
    for(int b = 1; b < nY; b++){
      for(int c = 1; c < nZ; c++){
          out.writeDouble(m.at(a,b,c));
      }
      out.newLine();
    }
  }
  out.close();
  return;
}


void simulator::calcTerms(int i, int j, int k){
  U1x = (U1.at(i + 1, j, k) - U1.at(i - 1, j, k)) / (2 * dx);
  U2x = (U2.at(i + 1, j, k) - U2.at(i - 1, j, k)) / (2 * dx);
  U1y = (U1.at(i,j + 1,k) - U1.at(i,j - 1,k)) / (2 * dy);
  U2y = (U2.at(i,j + 1,k) - U2.at(i,j - 1,k)) / (2 * dy);
  U1z = (U1.at(i,j,k + 1) - U1.at(i,j,k - 1)) / (2 * dz);
  U2z = (U2.at(i,j,k + 1) - U2.at(i,j,k - 1)) / (2 * dz);

  U1xx = (U1.at(i + 1,j,k) - 2 * U1.at(i,j,k) + U1.at(i - 1,j,k)) / (dx * dx);
  U2xx = (U2.at(i + 1,j,k) - 2 * U2.at(i,j,k) + U2.at(i - 1,j,k)) / (dx * dx);
  U1yy = (U1.at(i,j + 1,k) - 2 * U1.at(i,j,k) + U1.at(i,j - 1,k)) / (dy * dy);
  U2yy = (U2.at(i,j + 1,k) - 2 * U2.at(i,j,k) + U2.at(i,j - 1,k)) / (dy * dy);
  U1zz = (U1.at(i,j,k + 1) - 2 * U1.at(i,j,k) + U1.at(i,j,k - 1)) / (dz * dz);
  U2zz = (U2.at(i,j,k + 1) - 2 * U2.at(i,j,k) + U2.at(i,j,k - 1)) / (dz * dz);

  P1x = (P1.at(i + 1,j,k) - P1.at(i - 1,j,k)) / (2 * dx);
  P2x = (P2.at(i + 1,j,k) - P2.at(i - 1,j,k)) / (2 * dx);
  P1y = (P1.at(i,j + 1,k) - P1.at(i,j - 1,k)) / (2 * dy);
  P2y = (P2.at(i,j + 1,k) - P2.at(i,j - 1,k)) / (2 * dy);
  P1z = (P1.at(i,j,k + 1) - P1.at(i,j,k - 1)) / (2 * dz);
  P2z = (P2.at(i,j,k + 1) - P2.at(i,j,k - 1)) / (2 * dz);

  V1x = (V1.at(i + 1,j,k) - V1.at(i - 1,j,k)) / (2 * dx);
  V2x = (V2.at(i + 1,j,k) - V2.at(i - 1,j,k)) / (2 * dx);
  V1y = (V1.at(i,j + 1,k) - V1.at(i,j - 1,k)) / (2 * dy);
  V2y = (V2.at(i,j + 1,k) - V2.at(i,j - 1,k)) / (2 * dy);
  V1z = (V1.at(i,j,k + 1) - V1.at(i,j,k - 1)) / (2 * dz);
  V2z = (V2.at(i,j,k + 1) - V2.at(i,j,k - 1)) / (2 * dz);

  V1xx = (V1.at(i + 1,j,k) - 2 * V1.at(i,j,k) + V1.at(i - 1,j,k)) / (dx * dx);
  V2xx = (V2.at(i + 1,j,k) - 2 * V2.at(i,j,k) + V2.at(i - 1,j,k)) / (dx * dx);
  V1yy = (V1.at(i,j + 1,k) - 2 * V1.at(i,j,k) + V1.at(i,j - 1,k)) / (dy * dy);
  V2yy = (V2.at(i,j + 1,k) - 2 * V2.at(i,j,k) + V2.at(i,j - 1,k)) / (dy * dy);
  V1zz = (V1.at(i,j,k + 1) - 2 * V1.at(i,j,k) + V1.at(i,j,k - 1)) / (dz * dz);
  V2zz = (V2.at(i,j,k + 1) - 2 * V2.at(i,j,k) + V2.at(i,j,k - 1)) / (dz * dz);

  P1xx = (P1.at(i + 1,j,k) - 2 * P1.at(i,j,k) + P1.at(i - 1,j,k)) / (dx * dx);
  P1yy = (P1.at(i,j + 1,k) - 2 * P1.at(i,j,k) + P1.at(i,j - 1,k)) / (dy * dy);
  P1zz = (P1.at(i,j,k + 1) - 2 * P1.at(i,j,k) + P1.at(i,j,k - 1)) / (dz * dz);
  P2xx = (P2.at(i + 1,j,k) - 2 * P2.at(i,j,k) + P2.at(i - 1,j,k)) / (dx * dx);
  P2yy = (P2.at(i,j + 1,k) - 2 * P2.at(i,j,k) + P2.at(i,j - 1,k)) / (dy * dy);
  P2zz = (P2.at(i,j,k + 1) - 2 * P2.at(i,j,k) + P2.at(i,j,k - 1)) / (dz * dz);
}
#endif
