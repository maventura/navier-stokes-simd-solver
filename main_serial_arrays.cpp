#include <stdio.h>
#include <string.h>
#include <iostream>
#include <chrono>       // for high_resolution_clock.
#include <math.h>
#include <vector>
#include <utility>      // std::pair
#include <sstream>
#include "mat2.h"
#include "parameters.h"

using namespace std;

void clearScreen() {
    cerr << string( 100, '\n' );
}

long double diff(mat2 &A, mat2 &B) {
    long double sum = 0;
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            sum += pow(B.at(i, j) - A.at(i, j), 2);
        }
    }
    return sqrt(sum);
}


void setPBorders(mat2 &P0, mat2 &P1, mat2 &P2, int nX, int nY) {
    for (int c = 0; c < nX; ++c) {
        P0.set(c, 0, P0.at(c, 1));
        P1.set(c, 0, P1.at(c, 1));
        P2.set(c, 0, P2.at(c, 1));

        P0.set(c, nY - 1, P0.at(c, nY - 2));
        P1.set(c, nY - 1, P1.at(c, nY - 2));
        P2.set(c, nY - 1, P2.at(c, nY - 2));
    }
    for (int c = 0; c < nY; ++c) {
        P0.set(nX - 1, c , P0.at(nX - 2, c));
        P1.set(nX - 1, c , P1.at(nX - 2, c));
        P2.set(nX - 1, c , P2.at(nX - 2, c));
        P0.set(0, c, P0.at(1, c));
        P1.set(0, c, P1.at(1, c));
        P2.set(0, c, P2.at(1, c));
    }
}

int main() {



    //0 es tiempo n-1, 1 es tiempo n, 2 es tiempo n+1
    //3 es simplemente las matrices auxiliares para no perder datos de el tiempo 2
    mat2 U0(nX, nY, 0);
    mat2 V0(nX, nY, 0);
    mat2 P0(nX, nY, 0);
    mat2 U1(nX, nY, 0);
    mat2 V1(nX, nY, 0);
    mat2 P1(nX, nY, 0);
    mat2 U2(nX, nY, 0);
    mat2 V2(nX, nY, 0);
    mat2 P2(nX, nY, 0);
    mat2 U2_old(nX, nY, 0);
    mat2 V2_old(nX, nY, 0);
    mat2 P2_old(nX, nY, 0);

    unsigned long int step = 0;

    for (long double t = 0.0; t < tMax; t = t + dt) {
        if (step % printPercentageSteps == 0) cerr << 100 * t / tMax << "%" << endl;
        step++;
        if(100 * t / tMax > percentageStop) return 0;

        long double dFanAngle = fanTurns * 2 * pi / nT; //
        fanAngle += dFanAngle; //TODO: solo funciona en el 1er y tercer cuadrante.
        if (fanAngle > 2 * pi) fanAngle = 0;

        //cerr << 100 * t / tMax << "%" << endl ;
        if (isnan(U1.at(3, 3))) {
            cerr << "ERROR: nan found" << endl;
            exit(EXIT_FAILURE);
        }
       
        setPBorders(P0, P1, P2, nX, nY);

        for (int i = 1; i < nX - 1; ++i) {
            for (int j = 1; j < nY - 1; ++j) {

                long double U1x = (U1.at(i + 1, j) - U1.at(i - 1, j)) / (2.0 * dx);
                long double U2x = (U2.at(i + 1, j) - U2.at(i - 1, j)) / (2.0 * dx);
                long double U1y = (U1.at(i, j + 1) - U1.at(i, j - 1)) / (2.0 * dy);
                long double U2y = (U2.at(i, j + 1) - U2.at(i, j - 1)) / (2.0 * dy);
                long double U1xx = (U1.at(i + 1, j) - 2.0 * U1.at(i, j) + U1.at(i - 1, j)) / (dx * dx);
                long double U2xx = (U2.at(i + 1, j) - 2.0 * U2.at(i, j) + U2.at(i - 1, j)) / (dx * dx);
                long double U1yy = (U1.at(i, j + 1) - 2.0 * U1.at(i, j) + U1.at(i, j - 1)) / (dy * dy);
                long double U2yy = (U2.at(i, j + 1) - 2.0 * U2.at(i, j) + U2.at(i, j - 1)) / (dy * dy);
                long double P1x = (P1.at(i + 1, j) - P1.at(i - 1, j)) / (2.0 * dx);
                long double P2x = (P2.at(i + 1, j) - P2.at(i - 1, j)) / (2.0 * dx);
                long double P1y = (P1.at(i, j + 1) - P1.at(i, j - 1)) / (2.0 * dy);
                long double P2y = (P2.at(i, j + 1) - P2.at(i, j - 1)) / (2.0 * dy);
                long double V1x = (V1.at(i + 1, j) - V1.at(i - 1, j)) / (2.0 * dx);
                long double V2x = (V2.at(i + 1, j) - V2.at(i - 1, j)) / (2.0 * dx);
                long double V1y = (V1.at(i, j + 1) - V1.at(i, j - 1)) / (2.0 * dy);
                long double V2y = (V2.at(i, j + 1) - V2.at(i, j - 1)) / (2.0 * dy);
                long double V1xx = (V1.at(i + 1, j) - 2.0 * V1.at(i, j) + V1.at(i - 1, j)) / (dx * dx);
                long double V2xx = (V2.at(i + 1, j) - 2.0 * V2.at(i, j) + V2.at(i - 1, j)) / (dx * dx);
                long double V1yy = (V1.at(i, j + 1) - 2.0 * V1.at(i, j) + V1.at(i, j - 1)) / (dy * dy);
                long double V2yy = (V2.at(i, j + 1) - 2.0 * V2.at(i, j) + V2.at(i, j - 1)) / (dy * dy);
                long double P1xx = (P1.at(i + 1, j) - 2.0 * P1.at(i, j) + P1.at(i - 1, j)) / (dx * dx);
                long double P1yy = (P1.at(i, j + 1) - 2.0 * P1.at(i, j) + P1.at(i, j - 1)) / (dy * dy);
                long double P2yy = (P2.at(i, j + 1) - 2.0 * P2.at(i, j) + P2.at(i, j - 1)) / (dy * dy);



                long double u2val = U1.at(i, j) - U1.at(i, j) * (dt / dx) * (U1.at(i, j) - U1.at(i - 1, j)) - V1.at(i, j) * (dt / dy) * (U1.at(i, j) - U1.at(i, j - 1))
                                    - (dt / (rho * 2 * dx)) * (P1.at(i + 1, j) - P1.at(i - 1, j))
                                    + nu * ((dt / (dx * dx)) * (U1.at(i + 1, j) - 2 * U1.at(i, j) + U1.at(i - 1, j)) + (dt / (dy * dy)) * (U1.at(i, j + 1) - 2 * U1.at(i, j) + U1.at(i, j - 1)));
                U2.set(i, j, u2val);

                long double v2val = V1.at(i, j) - U1.at(i, j) * (dt / dx) * (V1.at(i, j) - V1.at(i - 1, j)) - V1.at(i, j) * (dt / dy) * (V1.at(i, j) - V1.at(i, j - 1))
                                    - (dt / (rho * 2 * dy)) * (P1.at(i, j + 1) - P1.at(i, j - 1))
                                    + nu * ((dt / (dx * dx)) * (V1.at(i + 1, j) - 2 * V1.at(i, j) + V1.at(i - 1, j)) + (dt / (dy * dy)) * (V1.at(i, j + 1) - 2 * V1.at(i, j) + V1.at(i, j - 1)));
                V2.set(i, j, v2val);

                long double finalTerm = (1 / dt) * (U1x + V1y) - U1x * U1x - 2 * U1y * V1x - V1y * V1y;
                long double p2val = ((P1.at(i + 1, j) + P1.at(i - 1, j)) * dy * dy + (P1.at(i, j + 1) + P1.at(i, j - 1)) * dx * dx) * (1.0 / (2 * (dx * dx + dy * dy)));
                p2val -= (rho * dx * dx * dy * dy / (2 * dx * dx + 2 * dy * dy)) * finalTerm;
                P2.set(i, j, p2val);



                long double x = i * dx - xc;
                long double y = j * dy - yc;


                 if (y > nY / 2) y -= nY / 2;
                    long double theta = atan2(y , x);
                    theta += pi;
                    long double beta = theta + (pi / 2.0);

                    long double r = sqrt(x * x + y * y);
                    long double tanVel = dFanAngle * r;
                    long double tanVelX = tanVel*cos(beta);
                    long double tanVelY = tanVel*sin(beta);

                    long double relVelX = tanVelX - U2.at(i, j);
                    long double relVelY = tanVelY - V2.at(i, j);
                    long double Fx = 0.5 * fanArea * C_d*relVelX / dt;
                    long double Fy = 0.5 * fanArea * C_d*relVelY / dt; //Mass(dt*dx since water density is 1)
                    //times tangent expected speed
                    //divided by time


                    long double Fu = Fx;//F * cos(beta);
                    long double Fv = Fy;//F * sin(beta);

                    long double angleDif = fabs(theta - fanAngle);
                    if ( angleDif < fanWidth  && r > rMin && r < rMax) {
                        long double strMult = 1 - angleDif / fanWidth;
                        U2.add(i, j, -Fu * dt * strMult);
                        V2.add(i, j, -Fv * dt * strMult);
                    }

            }
        }

        U0.setAll(U1);
        U1.setAll(U2);
        V0.setAll(V1);
        V1.setAll(V2);
        P0.setAll(P1);
        P1.setAll(P2);

         if (step % stepsUntilPrint == 0) {
            U0.print();
            V0.print();
        }

    }
}