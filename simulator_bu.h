#ifndef simulator_H
#define simulator_H

#include <math.h>
#include "mat3.h"
//#include "mat3.h"
#include "io.h"

class simulator{

public:
  simulator();
  void process();

private:
  double xMax,  yMax,  zMax,  tMax;
  double nu, rho, C_d;
  double dx,  dy,  dz,  dt;
  int nX, nY, nZ, nT, iter;

  double al, fixedPointError, minFixedPointIters;
  double pi = atan(1)*4;
  int stepsUntilPrint, printPercentageSteps, maxSteps;
  bool printPercentage;
  double t;

  //Calculation variables,
  double U1x, U2x, U1y, U2y, U1z, U2z, U1xx, U2xx, U1yy, U2yy, U1zz, U2zz;
  double V1xx, V2xx, V1yy, V2yy, V1zz, V2zz, P1xx, P1yy, P1zz, P2xx, P2yy, P2zz;
  double W1x, W2x, W1y, W2y, W1z, W2z, W1xx, W2xx, W1yy, W2yy, W1zz, W2zz;
  double P1x, P2x, P1y, P2y, P1z, P2z, V1x, V2x, V1y, V2y, V1z, V2z;

  double U1xb, U2xb, U1yb, U2yb, U1zb, U2zb;
  double V1xb, V2xb, V1yb, V2yb, V1zb, V2zb;
  double W1xb, W2xb, W1yb, W2yb, W1zb, W2zb;

  double U1xf, U2xf, U1yf, U2yf, U1zf, U2zf;
  double V1xf, V2xf, V1yf, V2yf, V1zf, V2zf;
  double W1xf, W2xf, W1yf, W2yf, W1zf, W2zf;

  mat3 U0, V0, W0, P0;
  mat3 U1, V1, W1, P1;
  mat3 U2, V2, W2, P2;

  void calcTerms(int i, int j, int k);
  void saveVtk(mat3 &state, string file_name);
  void saveCSV(mat3 &state, string file_name);

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
    if(tag == "dz") file.readDouble(dz);
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
        // U0.set(i, nY - 1, k, 1.0);
        // U1.set(i, nY - 1, k, 1.0);
        // U2.set(i, nY - 1, k, 1.0);

        U1.set(2,2,2,0.1);
        U2.set(2,2,2,0.1);
        U0.set(2,2,2,0.1);


        U1.set(8,2,2,-0.1);
        U2.set(8,2,2,-0.1);
        U0.set(8,2,2,-0.1);


        P0.set(i,nY-2,k, P0.at(i,nY-1,k));
        P1.set(i,nY-2,k, P0.at(i,nY-1,k));
        P2.set(i,nY-2,k, P0.at(i,nY-1,k));

        P0.set(i,1,k, P0.at(i,0,k));
        P1.set(i,1,k, P1.at(i,0,k));
        P2.set(i,1,k, P2.at(i,0,k));

        P0.set(nX-1,i,k, P0.at(nX-2,i,k));
        P1.set(nX-1,i,k, P1.at(nX-2,i,k));
        P2.set(nX-1,i,k, P2.at(nX-2,i,k));

        P0.set(0,i,k, P0.at(1,i,k));
        P1.set(0,i,k, P1.at(1,i,k));
        P2.set(0,i,k, P2.at(1,i,k));

        P0.set(i,k,nZ-1, P0.at(i,k,nZ-2));
        P1.set(i,k,nZ-1, P1.at(i,k,nZ-2));
        P2.set(i,k,nZ-1, P2.at(i,k,nZ-2));

        P0.set(i,k,0, P0.at(i,k,1));
        P1.set(i,k,0, P1.at(i,k,1));
        P2.set(i,k,0, P2.at(i,k,1));
    }
}
//some point has to be always equal to zero to act
//as reference for the potential.
P0.set(5,5,0, 0);
P1.set(5,5,0, 0);
P2.set(5,5,0, 0);
}


void simulator::process(){
  iter = 0; //TODO: Iter has two uses, fix.
  for (t = 0; t < tMax; t = t + dt) {
    cerr << 100*t/tMax << "%" << endl;
    iter++;
    ostringstream U_name;
    ostringstream V_name;
    ostringstream W_name;
    ostringstream P_name;

    U_name << "./out/U_" << iter << ".vtk";
    ostringstream V0_name;
    V_name << "./out/V_" << iter << ".vtk";
    ostringstream W0_name;
    W_name << "./out/W_" << iter << ".vtk";
    ostringstream P0_name;
    P_name << "./out/P_" << iter << ".vtk";
    saveVtk(U0,U_name.str());
    saveVtk(V0,V_name.str());
    saveVtk(W0,W_name.str());
    saveVtk(P0,P_name.str());

    for (int iter = 0; iter < 10; ++iter) { //TODO: Iter termination condition
      for (int i = 1; i < nX - 1; ++i) {
        for (int j = 1; j < nY - 1; ++j) {
          for (int k = 1; k < nZ - 1; ++k) {

            calcTerms(i,j,k);

            //spatial backward:
            // double newU = U1.at(i,j,k) + dt * (-al * U1.at(i,j,k) * U1xb - (1-al) * U1.at(i,j,k) * U2xb - al * V1.at(i,j,k) * U1yb - (1-al) * V1.at(i,j,k) * U2yb - al * W1.at(i,j,k) * U1zb - (1-al) * W1.at(i,j,k) * U2zb - (1/rho) * (al * P1x + (1-al) * P2x) + nu * (al * U1xx + (1-al) * U2xx + al * U1yy + (1-al) * U2yy + al * U1zz + (1-al) * U2zz) );
            // double newV = V1.at(i,j,k) + dt * (-al * U1.at(i,j,k) * V1xb - (1-al) * U1.at(i,j,k) * V2xb - al * V1.at(i,j,k) * V1yb - (1-al) * V1.at(i,j,k) * V2yb - al * W1.at(i,j,k) * V1zb - (1-al) * W1.at(i,j,k) * V2zb - (1/rho) * (al * P1y + (1-al) * P2y) + nu * (al * V1xx + (1-al) * V2xx + al * V1yy + (1-al) * V2yy + al * V1zz + (1-al) * V2zz) );
            // double newW = W1.at(i,j,k) + dt * (-al * U1.at(i,j,k) * W1xb - (1-al) * U1.at(i,j,k) * W2xb - al * V1.at(i,j,k) * W1yb - (1-al) * V1.at(i,j,k) * W2yb - al * W1.at(i,j,k) * W1zb - (1-al) * W1.at(i,j,k) * W2zb - (1/rho) * (al * P1z + (1-al) * P2z) + nu * (al * W1xx + (1-al) * W2xx + al * W1yy + (1-al) * W2yy + al * W1zz + (1-al) * W2zz) );

            //Spatial centered.
            double newU = U1.at(i,j,k) + dt * (-al * U1.at(i,j,k) * U1x - (1-al) * U1.at(i,j,k) * U2x - al * V1.at(i,j,k) * U1y - (1-al) * V1.at(i,j,k) * U2y - al * W1.at(i,j,k) * U1z - (1-al) * W1.at(i,j,k) * U2z - (1/rho) * (al * P1x + (1-al) * P2x) + nu * (al * U1xx + (1-al) * U2xx + al * U1yy + (1-al) * U2yy + al * U1zz + (1-al) * U2zz) );
            double newV = V1.at(i,j,k) + dt * (-al * U1.at(i,j,k) * V1x - (1-al) * U1.at(i,j,k) * V2x - al * V1.at(i,j,k) * V1y - (1-al) * V1.at(i,j,k) * V2y - al * W1.at(i,j,k) * V1z - (1-al) * W1.at(i,j,k) * V2z - (1/rho) * (al * P1y + (1-al) * P2y) + nu * (al * V1xx + (1-al) * V2xx + al * V1yy + (1-al) * V2yy + al * V1zz + (1-al) * V2zz) );
            double newW = W1.at(i,j,k) + dt * (-al * U1.at(i,j,k) * W1x - (1-al) * U1.at(i,j,k) * W2x - al * V1.at(i,j,k) * W1y - (1-al) * V1.at(i,j,k) * W2y - al * W1.at(i,j,k) * W1z - (1-al) * W1.at(i,j,k) * W2z - (1/rho) * (al * P1z + (1-al) * P2z) + nu * (al * W1xx + (1-al) * W2xx + al * W1yy + (1-al) * W2yy + al * W1zz + (1-al) * W2zz) );


            //original:
            //double A = (1/dt)*(U1x + V1y + W1z) - (U1x*U1x + V1y*V1y + W1z*W1z + 2*U1y*V1x + 2*U1z*W1x + 2*W1y*V1z);
            //double newP = P2.at(i+1,j,k)+P2.at(i-1,j,k) - dx*dx*rho*A;

            //P1:
            //double A = (1/dt)*(U1x + V1y + W1z) - (U1x*U1x + V1y*V1y + W1z*W1z + 2*U1y*V1x + 2*U1z*W1x + 2*W1y*V1z);
            //double newP = P1.at(i+1,j,k)+P1.at(i-1,j,k) - dx*dx*rho*A;

            //literal:
            //double newP = P1.at(i,j,k) -P1xx  -P1yy -P1zz -rho*(U1x*U1x + V1y*V1y + W1z*W1z + 2*U1y*V1x + 2*U1z*W1x + 2*W1y*V1z);

            //discretized literal:
             double newP  = -0.5*(-(P1.at(i + 1,j,k) + P1.at(i - 1,j,k)) + dx*dx*(-( P1yy + P1zz)+rho*((1/dt)*(U1xx + V1yy) - (U1x*U1x + V1y*V1y + W1z*W1z + 2*U1y*V1x + 2*U1z*W1x + 2*W1y*V1z))));

            if(isnan(newP)){
              cerr << "Error: Nan found. Returning." << endl;
              return;
            }

            U2.set(i,j,k, newU);
            V2.set(i,j,k, newV);
            W2.set(i,j,k, newW);
            P2.set(i,j,k, newP);

          //simple diffusion for testing.
          //  U2.set(i,j,k,(U1.at(i,j,k)+U1.at(i+1,j,k)+U1.at(i,j+1,k)+U1.at(i,j,k+1)+U1.at(i-1,j,k)+U1.at(i,j-1,k) +U1.at(i,j,k-1))/7.0 );
          //  if(i == 10 && j == 10 && k == 10) U2.set(i,j,k, 50);

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


void simulator::saveVtk(mat3 &m, string file_name){
// Save a 3-D scalar array in VTK format.
  io out(file_name, io::type_write);
  out.write("# vtk DataFile Version 2.0");
  out.newLine();
  out.write("Comment goes here");
  out.newLine();
  out.write("ASCII");
  out.newLine();
  out.newLine();

  out.write("DATASET STRUCTURED_POINTS");
  out.newLine();
  out.write("DIMENSIONS    " + to_string(nX) + " " +  to_string(nY) + " " + to_string(nZ));
  out.newLine();
  out.newLine();

  out.write("ORIGIN    0.000   0.000   0.000");
  out.newLine();
  out.write("SPACING    1.000   1.000   1.000");
  out.newLine();
  out.newLine();
  out.write("POINT_DATA   " + to_string(nX*nY*nZ));
  out.newLine();
  out.write("SCALARS scalars float");
  out.newLine();
  out.write("LOOKUP_TABLE default");
  out.newLine();
  out.newLine();

  for (int i = 0; i < nX; ++i) {
    for (int j = 0; j < nY; ++j) {
      for (int k = 0; k < nZ; ++k) {
          out.writeDouble(m.at(i,j,k));
          out.write(" ");
      }
      out.newLine();
    }
  }
  out.close();
}



void simulator::saveCSV(mat3 &m, string file_name){
// Save a 3-D scalar array in VTK format.
  io out(file_name, io::type_write);

  for (int i = 0; i < nX - 1; ++i) {
    for (int j = 0; j < nY - 1; ++j) {
      for (int k = 0; k < nZ; ++k) {
          out.writeDouble(m.at(i,j,k));
          out.write(", ");

      }
      out.newLine();
    }
  }
  out.close();
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

  W1x = (W1.at(i + 1,j,k) - W1.at(i - 1,j,k)) / (2 * dx);
  W2x = (W2.at(i + 1,j,k) - W2.at(i - 1,j,k)) / (2 * dx);
  W1y = (W1.at(i,j + 1,k) - W1.at(i,j - 1,k)) / (2 * dy);
  W2y = (W2.at(i,j + 1,k) - W2.at(i,j - 1,k)) / (2 * dy);
  W1z = (W1.at(i,j,k + 1) - W1.at(i,j,k - 1)) / (2 * dz);
  W2z = (W2.at(i,j,k + 1) - W2.at(i,j,k - 1)) / (2 * dz);

  W1xx = (W1.at(i + 1,j,k) - 2 * W1.at(i,j,k) + W1.at(i - 1,j,k)) / (dx * dx);
  W2xx = (W2.at(i + 1,j,k) - 2 * W2.at(i,j,k) + W2.at(i - 1,j,k)) / (dx * dx);
  W1yy = (W1.at(i,j + 1,k) - 2 * W1.at(i,j,k) + W1.at(i,j - 1,k)) / (dy * dy);
  W2yy = (W2.at(i,j + 1,k) - 2 * W2.at(i,j,k) + W2.at(i,j - 1,k)) / (dy * dy);
  W1zz = (W1.at(i,j,k + 1) - 2 * W1.at(i,j,k) + W1.at(i,j,k - 1)) / (dz * dz);
  W2zz = (W2.at(i,j,k + 1) - 2 * W2.at(i,j,k) + W2.at(i,j,k - 1)) / (dz * dz);

  P1x = (P1.at(i + 1,j,k) - P1.at(i - 1,j,k)) / (2 * dx);
  P2x = (P2.at(i + 1,j,k) - P2.at(i - 1,j,k)) / (2 * dx);
  P1y = (P1.at(i,j + 1,k) - P1.at(i,j - 1,k)) / (2 * dy);
  P2y = (P2.at(i,j + 1,k) - P2.at(i,j - 1,k)) / (2 * dy);
  P1z = (P1.at(i,j,k + 1) - P1.at(i,j,k - 1)) / (2 * dz);
  P2z = (P2.at(i,j,k + 1) - P2.at(i,j,k - 1)) / (2 * dz);

  P1xx = (P1.at(i + 1,j,k) - 2 * P1.at(i,j,k) + P1.at(i - 1,j,k)) / (dx * dx);
  P1yy = (P1.at(i,j + 1,k) - 2 * P1.at(i,j,k) + P1.at(i,j - 1,k)) / (dy * dy);
  P1zz = (P1.at(i,j,k + 1) - 2 * P1.at(i,j,k) + P1.at(i,j,k - 1)) / (dz * dz);
  P2xx = (P2.at(i + 1,j,k) - 2 * P2.at(i,j,k) + P2.at(i - 1,j,k)) / (dx * dx);
  P2yy = (P2.at(i,j + 1,k) - 2 * P2.at(i,j,k) + P2.at(i,j - 1,k)) / (dy * dy);
  P2zz = (P2.at(i,j,k + 1) - 2 * P2.at(i,j,k) + P2.at(i,j,k - 1)) / (dz * dz);

  //Backward.
    U1xb = (U1.at(i + 1, j, k) - U1.at(i , j, k)) / dx;
    U2xb = (U2.at(i + 1, j, k) - U2.at(i, j, k)) / dx;
    U1yb = (U1.at(i,j + 1, k) - U1.at(i,j,k)) / dy;
    U2yb = (U2.at(i,j + 1, k) - U2.at(i,j,k)) / dy;
    U1zb = (U1.at(i,j,k + 1) - U1.at(i,j,k)) / dz;
    U2zb = (U2.at(i,j,k + 1) - U2.at(i,j,k)) / dz;

    V1xb = (V1.at(i + 1,j,k) - V1.at(i,j,k)) / dx;
    V2xb = (V2.at(i + 1,j,k) - V2.at(i,j,k)) / dx;
    V1yb = (V1.at(i,j + 1,k) - V1.at(i,j,k)) / dy;
    V2yb = (V2.at(i,j + 1,k) - V2.at(i,j,k)) / dy;
    V1zb = (V1.at(i,j,k + 1) - V1.at(i,j,k)) / dz;
    V2zb = (V2.at(i,j,k + 1) - V2.at(i,j,k)) / dz;

    W1xb = (W1.at(i + 1,j,k) - W1.at(i,j,k)) / dx;
    W2xb = (W2.at(i + 1,j,k) - W2.at(i,j,k)) / dx;
    W1yb = (W1.at(i,j + 1,k) - W1.at(i,j,k)) / dy;
    W2yb = (W2.at(i,j + 1,k) - W2.at(i,j,k)) / dy;
    W1zb = (W1.at(i,j,k + 1) - W1.at(i,j,k)) / dz;
    W2zb = (W2.at(i,j,k + 1) - W2.at(i,j,k)) / dz;

    //Forward.
      U1xf = (U1.at(i, j, k) - U1.at(i - 1, j, k)) / dx;
      U2xf = (U2.at(i, j, k) - U2.at(i - 1, j, k)) / dx;
      U1yf = (U1.at(i,j,k) - U1.at(i,j - 1,k)) / dy;
      U2yf = (U2.at(i,j,k) - U2.at(i,j - 1,k)) / dy;
      U1zf = (U1.at(i,j,k) - U1.at(i,j,k - 1)) / dz;
      U2zf = (U2.at(i,j,k) - U2.at(i,j,k - 1)) / dz;

      V1xf = (V1.at(i,j,k) - V1.at(i - 1,j,k)) / dx;
      V2xf = (V2.at(i,j,k) - V2.at(i - 1,j,k)) / dx;
      V1yf = (V1.at(i,j,k) - V1.at(i,j - 1,k)) / dy;
      V2yf = (V2.at(i,j,k) - V2.at(i,j - 1,k)) / dy;
      V1zf = (V1.at(i,j,k) - V1.at(i,j,k - 1)) / dz;
      V2zf = (V2.at(i,j,k) - V2.at(i,j,k - 1)) / dz;

      W1xf = (W1.at(i,j,k) - W1.at(i - 1,j,k)) / dx;
      W2xf = (W2.at(i,j,k) - W2.at(i - 1,j,k)) / dx;
      W1yf = (W1.at(i,j,k) - W1.at(i,j - 1,k)) / dy;
      W2yf = (W2.at(i,j,k) - W2.at(i,j - 1,k)) / dy;
      W1zf = (W1.at(i,j,k) - W1.at(i,j,k - 1)) / dz;
      W2zf = (W2.at(i,j,k) - W2.at(i,j,k - 1)) / dz;


}
#endif