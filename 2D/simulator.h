#ifndef simulator_H
#define simulator_H

#include <math.h>
#include "mat2.h"
#include <cmath>
#include "io.h"
#include <algorithm>    // std::min


class simulator {

  public:
    simulator();
    void process();
  private:
    float xMax, yMax, tMax;
    float nu, rho, C_d;
    float dx, dy, dt, h;
    int nX, nY, nT, step, maxIters;

    float xc, yc;
    float rMax, rMin, fanTurns, startinAngle;

    float al, fixedPointError, minFixedPointIters;
    float pi = atan(1) * 4;

    float U1x, U2x, U2y, U1y;
    float U1xx ,U2xx ,U1yy,U2yy;
    float P1x, P2x, P1y, P2y;
    float V1x, V2x, V1y, V2y;
    float V1xx, V2xx, V1yy, V2yy;
    float P1xx, P1yy, P2yy;

    int printPercentageSteps, maxSteps;
    bool printPercentage;
    float t;

    float percentageStop, fanArea, fanWidth;
    //TODO: Add option to import text.
    //0 es tiempo n-1, 1 es tiempo n, 2 es tiempo n+1
    mat2 U0, V0, P0;
    mat2 U1, V1, P1;
    mat2 U2, V2, P2;

    float Re, Fx, Fy, Fz;

    void centralSpeed();
    void centralPressure(int i, int j, int k);

    void calcTerms(int i, int j);
    void calcVelocities(int i, int j);

    void saveVtk(mat2 &state, string file_name);
    void saveStreamVtk(mat2 &m1, mat2 &m2, string file_name);

    void setCavityFlowSpeeds();
    void setPBorders();
    void readParameters(string file_name);

    void simpleDiffusion(int i, int j, int k);

    void saveVelocitiesToFile();
    void saveVelocitiesToTxt();
    void appendTxt(string file_name, mat2 &M);

    bool isBorder(int i, int j, int k);


};

simulator::simulator() {
    readParameters("parameters.txt");
    h = dx; //for square grids.
    nX = round(xMax / dx) + 1;
    nY = round(yMax / dy) + 1;
    nT = round(tMax / dt) + 1;
    Re = (xMax * 0.01) / (nu); //TODO: 0.01 is the inicial velocity at y = yMax
    U0.confAndInit(nX, nY, 0);
    V0.confAndInit(nX, nY, 0);
    P0.confAndInit(nX, nY, 1);
    U1.confAndInit(nX, nY, 0);
    V1.confAndInit(nX, nY, 0);
    P1.confAndInit(nX, nY, 1);
    U2.confAndInit(nX, nY, 0);
    V2.confAndInit(nX, nY, 0);
    P2.confAndInit(nX, nY, 1);

    fanArea = dx * dy;
    xc = nX / 2;
    yc = nY / 2;
}

void simulator::readParameters(string file_name) {
    io file(file_name, io::type_read);
    string tag;
    while (file.readWord(tag)) {
        if (tag == "xMax") file.readFloat(xMax);
        if (tag == "yMax") file.readFloat(yMax);
        if (tag == "tMax") file.readFloat(tMax);
        if (tag == "nu") file.readFloat(nu);
        if (tag == "rho") file.readFloat(rho);
        if (tag == "dx") file.readFloat(dx);
        if (tag == "dy") file.readFloat(dy);
        if (tag == "dt") file.readFloat(dt);
        if (tag == "al") file.readFloat(al);
        if (tag == "fixedPointError") file.readFloat(fixedPointError);
        if (tag == "minFixedPointIters") file.readFloat(minFixedPointIters);
        if(tag == "rMax")file.readFloat(rMax);
        if(tag == "rMin")file.readFloat(rMin);
        if(tag == "fanTurns")file.readFloat(fanTurns);
        if(tag == "startinAngle")file.readFloat(startinAngle);
        if (tag == "printPercentageSteps") file.readInt(printPercentageSteps);
        if (tag == "C_d") file.readFloat(C_d);
        if (tag == "maxSteps") file.readInt(maxSteps);
        if (tag == "printPercentage") file.readBool(printPercentage);
        if (tag == "maxIters") file.readInt(maxIters);
        if (tag == "percentageStop") file.readFloat(percentageStop);
        if(tag == "fanWidth") file.readFloat(fanWidth);
    }
    file.close();
}


void simulator::setCavityFlowSpeeds() {
    for (int i = 0; i < nX; ++i) {
        //Set Cavity flow conditions.
        U0.set(i, nY - 1, 0.01);
        U1.set(i, nY - 1, 0.01);
        U2.set(i, nY - 1, 0.01);
    }
}
void simulator::setPBorders() {
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


void simulator::process() {
    cout << "cols: " << U2.cols() << ", rows: " << U2.rows() << endl;
    for (t = 0.0; t < tMax; t = t + dt) {
        if (step % printPercentageSteps == 0) 
            cerr << 100 * t / tMax << "% \r";
        step++;

        //saveVelocitiesToFile();
        saveVelocitiesToTxt();
        if(100 * t / tMax > percentageStop) return ;

        float dFanAngle = fanTurns * 2 * pi / nT; //
        startinAngle += dFanAngle; //TODO: solo funciona en el 1er y tercer cuadrante.
        if (startinAngle > 2 * pi) startinAngle = 0;

        //cerr << 100 * t / tMax << "%" << endl ;
        if (isnan(U1.at(3, 3))) {
            cerr << "ERROR: nan found" << endl;
            exit(EXIT_FAILURE);
        }
       
        setPBorders();

        for (int i = 1; i < nX - 1; ++i) {
            for (int j = 1; j < nY - 1; ++j) {

                calcTerms(i,j);
                centralSpeed();
                calcVelocities(i,j);




                //from fan scenario. 
               // float x = i * dx - xc;
               // float y = j * dy - yc;
                //TODO: Fan scenario needed this fix, why?
               // if (y > nY / 2) y -= nY / 2;
               // float theta = atan2(y , x);
               // theta += pi;
               // float beta = theta + (pi / 2.0);

               // float r = sqrt(x * x + y * y);
               // float tanVel = dFanAngle * r;
               // float tanVelX = tanVel*cos(beta);
               // float tanVelY = tanVel*sin(beta);

               // float relVelX = tanVelX - U2.at(i, j);
               // float relVelY = tanVelY - V2.at(i, j);
              //  float Fx = 0.5 * fanArea * C_d*relVelX / dt;
              //  float Fy = 0.5 * fanArea * C_d*relVelY / dt; //Mass(dt*dx since water density is 1)
                //times tangent expected speed
                //divided by time


              //  float Fu = Fx;//F * cos(beta);
               // float Fv = Fy;//F * sin(beta);

               // float angleDif = fabs(theta - startinAngle);
               // if ( angleDif < fanWidth  && r > rMin && r < rMax) {
               //     float strMult = 1 - angleDif / fanWidth;
               //     U2.add(i, j, -Fu * dt * strMult);
               //     V2.add(i, j, -Fv * dt * strMult);
               // }

            }
        }

        U0.copyAll(U1);
        U1.copyAll(U2);
        V0.copyAll(V1);
        V1.copyAll(V2);
        P0.copyAll(P1);
        P1.copyAll(P2);


    }
}

void simulator::calcVelocities(int i, int j){
    float u2val = U1.at(i, j) - U1.at(i, j) * (dt / dx) * (U1.at(i, j) - U1.at(i - 1, j)) - V1.at(i, j) * (dt / dy) * (U1.at(i, j) - U1.at(i, j - 1)) - (dt / (rho * 2 * dx)) * (P1.at(i + 1, j) - P1.at(i - 1, j))
                                    + nu * ((dt / (dx * dx)) * (U1.at(i + 1, j) - 2 * U1.at(i, j) + U1.at(i - 1, j)) + (dt / (dy * dy)) * (U1.at(i, j + 1) - 2 * U1.at(i, j) + U1.at(i, j - 1)));
    U2.set(i, j, u2val);

    float v2val = V1.at(i, j) - U1.at(i, j) * (dt / dx) * (V1.at(i, j) - V1.at(i - 1, j)) - V1.at(i, j) * (dt / dy) * (V1.at(i, j) - V1.at(i, j - 1))
                                    - (dt / (rho * 2 * dy)) * (P1.at(i, j + 1) - P1.at(i, j - 1))
                                    + nu * ((dt / (dx * dx)) * (V1.at(i + 1, j) - 2 * V1.at(i, j) + V1.at(i - 1, j)) + (dt / (dy * dy)) * (V1.at(i, j + 1) - 2 * V1.at(i, j) + V1.at(i, j - 1)));
    V2.set(i, j, v2val);

    float finalTerm = (1 / dt) * (U1x + V1y) - U1x * U1x - 2 * U1y * V1x - V1y * V1y;
    float p2val = ((P1.at(i + 1, j) + P1.at(i - 1, j)) * dy * dy + (P1.at(i, j + 1) + P1.at(i, j - 1)) * dx * dx) * (1.0 / (2 * (dx * dx + dy * dy)));
    p2val -= (rho * dx * dx * dy * dy / (2 * dx * dx + 2 * dy * dy)) * finalTerm;
    P2.set(i, j, p2val);
}

void simulator::calcTerms(int i, int j) {
    U1x = (U1.at(i + 1, j) - U1.at(i - 1, j)) / (2.0 * dx);
    U2x = (U2.at(i + 1, j) - U2.at(i - 1, j)) / (2.0 * dx);
    U1y = (U1.at(i, j + 1) - U1.at(i, j - 1)) / (2.0 * dy);
    U2y = (U2.at(i, j + 1) - U2.at(i, j - 1)) / (2.0 * dy);
    U1xx = (U1.at(i + 1, j) - 2.0 * U1.at(i, j) + U1.at(i - 1, j)) / (dx * dx);
    U2xx = (U2.at(i + 1, j) - 2.0 * U2.at(i, j) + U2.at(i - 1, j)) / (dx * dx);
    U1yy = (U1.at(i, j + 1) - 2.0 * U1.at(i, j) + U1.at(i, j - 1)) / (dy * dy);
    U2yy = (U2.at(i, j + 1) - 2.0 * U2.at(i, j) + U2.at(i, j - 1)) / (dy * dy);
    P1x = (P1.at(i + 1, j) - P1.at(i - 1, j)) / (2.0 * dx);
    P2x = (P2.at(i + 1, j) - P2.at(i - 1, j)) / (2.0 * dx);
    P1y = (P1.at(i, j + 1) - P1.at(i, j - 1)) / (2.0 * dy);
    P2y = (P2.at(i, j + 1) - P2.at(i, j - 1)) / (2.0 * dy);
    V1x = (V1.at(i + 1, j) - V1.at(i - 1, j)) / (2.0 * dx);
    V2x = (V2.at(i + 1, j) - V2.at(i - 1, j)) / (2.0 * dx);
    V1y = (V1.at(i, j + 1) - V1.at(i, j - 1)) / (2.0 * dy);
    V2y = (V2.at(i, j + 1) - V2.at(i, j - 1)) / (2.0 * dy);
    V1xx = (V1.at(i + 1, j) - 2.0 * V1.at(i, j) + V1.at(i - 1, j)) / (dx * dx);
    V2xx = (V2.at(i + 1, j) - 2.0 * V2.at(i, j) + V2.at(i - 1, j)) / (dx * dx);
    V1yy = (V1.at(i, j + 1) - 2.0 * V1.at(i, j) + V1.at(i, j - 1)) / (dy * dy);
    V2yy = (V2.at(i, j + 1) - 2.0 * V2.at(i, j) + V2.at(i, j - 1)) / (dy * dy);
    P1xx = (P1.at(i + 1, j) - 2.0 * P1.at(i, j) + P1.at(i - 1, j)) / (dx * dx);
    P1yy = (P1.at(i, j + 1) - 2.0 * P1.at(i, j) + P1.at(i, j - 1)) / (dy * dy);
    P2yy = (P2.at(i, j + 1) - 2.0 * P2.at(i, j) + P2.at(i, j - 1)) / (dy * dy);
}


bool simulator::isBorder(int i, int j, int k){
    return (j == nY - 1  || i == nX - 1 || i * j  == 0);
}

// ----------Functions for file management ----------

void simulator::saveVelocitiesToFile(){
    ostringstream U_name;
    ostringstream V_name;
    U_name << "./out/U_" << step << ".vtk";
    ostringstream V0_name;
    V_name << "./out/V_" << step << ".vtk";
    
    ostringstream stream_name;
    stream_name << "./out/stream_" << step << ".vtk";

    saveVtk(U2, U_name.str());
    saveVtk(V2, V_name.str());
}


void simulator::saveVelocitiesToTxt(){
    ostringstream name;
    name << "./out/UV.txt";

    appendTxt(name.str(),U2);
    appendTxt(name.str(),V2);
}




void simulator::centralSpeed() {
    for (int i = 0; i < nX; ++i) {
        for (int j = 0; j < nY; ++j) {
            float r = sqrt((i - xc) * (i - xc) + (j - yc) * (j - yc));
            if (r < 2.5) {
                float max = 0.01;
                U2.set(i, j, min(max / r, max));
                U1.set(i, j, min(max / r, max));
                U0.set(i, j, min(max / r, max));
                //TODO: we were using a custom min function, dont know why, test this one.
            }
        }
    }
}


#endif


void simulator::saveVtk(mat2 &m, string file_name) {
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
    out.write("DIMENSIONS    " + to_string(nX) + " " +  to_string(nY));
    out.newLine();
    out.newLine();
    out.write("ORIGIN    0.000   0.000   0.000");
    out.newLine();
    out.write("SPACING    1.000   1.000   1.000");
    out.newLine();
    out.newLine();
    out.write("POINT_DATA   " + to_string(nX * nY));
    out.newLine();
    out.write("SCALARS scalars float");
    out.newLine();
    out.write("LOOKUP_TABLE default");
    out.newLine();
    out.newLine();
    for (int i = 0; i < nX; ++i) {
        for (int j = 0; j < nY; ++j) {
            out.writeFloat(m.at(i, j));
            out.write(" ");
        }
    }
    out.close();
}


void simulator::appendTxt(string file_name, mat2 &M){
    io out(file_name, io::type_read_write_appending);
    for (int i = 0; i < M.rows(); ++i)
    {
        for (int j = 0; j < M.cols(); ++j)
        {
            out.writeFloat(M.at(i,j));
            out.write(" ");

        }
        out.newLine();
    }
    out.newLine();
    out.close();

}
