#include <iostream>
#include <vector>
#include <chrono>  // for high_resolution_clock.
#include <sstream>
#include <cmath>
#include "mat3.h"
using namespace std;


int main() {
    long double xMax = 2.0;
    long double yMax = 2.0;
    long double zMax = 2.0;
    long double tMax = 1.0;
    long double nu = 0.1;       //viscosity
    long double rho = 1;        //density
    long double dx = 0.1;
    long double dy = 0.1;
    long double dz = 0.1;
    long double dt = 0.1;
    int nX = round(xMax / dx) + 1;
    int nY = round(yMax / dy) + 1;
    int nZ = round(yMax / dy) + 1;
    int nT = round(tMax / dt) + 1;
    long double al = 0.5;

    mat3 U0(nX, nY, nZ, 0);
    mat3 V0(nX, nY, nZ, 0);
    mat3 W0(nX, nY, nZ, 0);
    mat3 P0(nX, nY, nZ, 0);

    mat3 U1(nX, nY, nZ, 0);
    mat3 V1(nX, nY, nZ, 0);
    mat3 W1(nX, nY, nZ, 0);
    mat3 P1(nX, nY, nZ, 0);

    mat3 U2(nX, nY, nZ, 0);
    mat3 V2(nX, nY, nZ, 0);
    mat3 W2(nX, nY, nZ, 0);
    mat3 P2(nX, nY, nZ, 0);

    mat3 U3(nX, nY, nZ, 0);
    mat3 V3(nX, nY, nZ, 0);
    mat3 W3(nX, nY, nZ, 0);
    mat3 P3(nX, nY, nZ, 0);

    for (int i = 0; i < nX; ++i) {
        for (int k = 0; k < nZ; ++k) {
            U0.set(i, nY - 1, k, 1.0);
            U1.set(i, nY - 1, k, 1.0);
            U2.set(i, nY - 1, k, 1.0);
            U3.set(i, nY - 1, k, 1.0);
        }
    }

    for (long double t = 0.0; t < tMax; t = t + dt) {
        U0.print();
        V0.print();
        W0.print();
        for (int i = 1; i < nX - 1; ++i) {
            for (int j = 1; j < nY - 1; ++j) {
                for (int k = 0; k < nZ; ++k) {
                    for (int iter = 0; iter < 10; ++iter) { //TODO: Iter termination condition

                        /*long double oldU = U2[i][j][k];
                        long double oldV = V2[i][j][k];
                        long double oldW = W2[i][j][k];*/

                        long double U1x = (U1.at(i + 1, j, k) - U1.at(i - 1, j, k)) / (2 * dx);
                        long double U2x = (U2.at(i + 1, j, k) - U2.at(i - 1, j, k)) / (2 * dx);
                        long double U1y = (U1[i][j + 1][k] - U1[i][j - 1][k]) / (2 * dy);
                        long double U2y = (U2[i][j + 1][k] - U2[i][j - 1][k]) / (2 * dy);
                        long double U1z = (U1[i][j][k + 1] - U1[i][j][k - 1]) / (2 * dz);
                        long double U2z = (U2[i][j][k + 1] - U2[i][j][k - 1]) / (2 * dz);

                        long double U1xx = (U1[i + 1][j][k] - 2 * U1[i][j][k] + U1[i - 1][j][k]) / (dx * dx);
                        long double U2xx = (U2[i + 1][j][k] - 2 * U2[i][j][k] + U2[i - 1][j][k]) / (dx * dx);
                        long double U1yy = (U1[i][j + 1][k] - 2 * U1[i][j][k] + U1[i][j - 1][k]) / (dy * dy);
                        long double U2yy = (U2[i][j + 1][k] - 2 * U2[i][j][k] + U2[i][j - 1][k]) / (dy * dy);
                        long double U1zz = (U1[i][j][k + 1] - 2 * U1[i][j][k] + U1[i][j][k - 1]) / (dz * dz);
                        long double U2zz = (U2[i][j][k + 1] - 2 * U2[i][j][k] + U2[i][j][k - 1]) / (dz * dz);

                        long double P1x = (P1[i + 1][j][k] - P1[i - 1][j][k]) / (2 * dx);
                        long double P2x = (P2[i + 1][j][k] - P2[i - 1][j][k]) / (2 * dx);
                        long double P1y = (P1[i][j + 1][k] - P1[i][j - 1][k]) / (2 * dy);
                        long double P2y = (P2[i][j + 1][k] - P2[i][j - 1][k]) / (2 * dy);
                        long double P1z = (P1[i][j][k + 1] - P1[i][j][k - 1]) / (2 * dz);
                        long double P2z = (P2[i][j][k + 1] - P2[i][j][k - 1]) / (2 * dz);

                        long double V1x = (V1[i + 1][j][k] - V1[i - 1][j][k]) / (2 * dx);
                        long double V2x = (V2[i + 1][j][k] - V2[i - 1][j][k]) / (2 * dx);
                        long double V1y = (V1[i][j + 1][k] - V1[i][j - 1][k]) / (2 * dy);
                        long double V2y = (V2[i][j + 1][k] - V2[i][j - 1][k]) / (2 * dy);
                        long double V1z = (V1[i][j][k + 1] - V1[i][j][k - 1]) / (2 * dz);
                        long double V2z = (V2[i][j][k + 1] - V2[i][j][k - 1]) / (2 * dz);

                        long double V1xx = (V1[i + 1][j][k] - 2 * V1[i][j][k] + V1[i - 1][j][k]) / (dx * dx);
                        long double V2xx = (V2[i + 1][j][k] - 2 * V2[i][j][k] + V2[i - 1][j][k]) / (dx * dx);
                        long double V1yy = (V1[i][j + 1][k] - 2 * V1[i][j][k] + V1[i][j - 1][k]) / (dy * dy);
                        long double V2yy = (V2[i][j + 1][k] - 2 * V2[i][j][k] + V2[i][j - 1][k]) / (dy * dy);
                        long double V1zz = (V1[i][j][k + 1] - 2 * V1[i][j][k] + V1[i][j][k - 1]) / (dz * dz);
                        long double V2zz = (V2[i][j][k + 1] - 2 * V2[i][j][k] + V2[i][j][k - 1]) / (dz * dz);

                        long double P1xx = (P1[i + 1][j][k] - 2 * P1[i][j][k] + P1[i - 1][j][k]) / (dx * dx);
                        long double P1yy = (P1[i][j + 1][k] - 2 * P1[i][j][k] + P1[i][j - 1][k]) / (dy * dy);
                        long double P1zz = (P1[i][j][k + 1] - 2 * P1[i][j][k] + P1[i][j][k - 1]) / (dz * dz);
                        long double P2xx = (P2[i + 1][j][k] - 2 * P2[i][j][k] + P2[i - 1][j][k]) / (dx * dx);
                        long double P2yy = (P2[i][j + 1][k] - 2 * P2[i][j][k] + P2[i][j - 1][k]) / (dy * dy);
                        long double P2zz = (P2[i][j][k + 1] - 2 * P2[i][j][k] + P2[i][j][k - 1]) / (dz * dz);






                        /* POR QUÉ ESTO?
                        U3[i][j][k] =  U2[i][j][k] + 0.001;
                        V3[i][j][k] =  V2[i][j][k] + 0.002;
                        P3[i][j][k] =  P2[i][j][k] + 0.003;*/

                        // long double diff = sqrt( pow(U3[i][j][k] - oldU, 2) + pow(V3[i][j][k] - oldV, 2) + pow(W3[i][j][k] - oldW, 2 ) );
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
        /*U2 = U3;
        V2 = V3;
        W2 = W3;
        P2 = P3;

        U0 = U1;
        V0 = V1;
        W0 = W1;
        P0 = P1;

        U1 = U2;
        V1 = V2;
        W1 = W2;
        P1 = P2;*/
    }
    cerr << "DBG: process returned without errors" << endl;
}
