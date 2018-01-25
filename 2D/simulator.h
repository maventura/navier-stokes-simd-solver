#ifndef simulator_H
#define simulator_H

#include <math.h>
#include "mat2.h"
#include <cmath>
#include "io.h"
#include <algorithm>    // std::min
#include <omp.h>

#ifdef USE_ASM
extern "C" {
    void vvp_asm(float *mats[], int pos, int offsetI, 
        float dt, float dx, float dy, float rho, float nu);}
#endif

class simulator {

  public:
    simulator();
    void process();
  private:
    float xMax, yMax, tMax;
    float nu, rho;
    float dx, dy, dt, h;
    int nX, nY, nT, step;

    float xc, yc;

    float al;
    float pi = atan(1) * 4;

    float U1x, U2x, U2y, U1y;
    float U1xx ,U2xx ,U1yy,U2yy;
    float P1x, P2x, P1y, P2y;
    float V1x, V2x, V1y, V2y;
    float V1xx, V2xx, V1yy, V2yy;
    float P1xx, P1yy, P2yy;

    int printPercentageSteps;
    bool printPercentage;
    float t;

    float percentageStop;
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

    void saveVelocitiesToFile();
    void saveVtk(mat2 &state, string file_name);

    void setCavityFlowSpeeds();
    void setPBorders();
    void readParameters(string file_name);

    void simpleDiffusion(int i, int j);
    void calcVelocitiesAsm(int i, int j);

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
    step = 0;
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
        if (tag == "printPercentageSteps") file.readInt(printPercentageSteps);
        if (tag == "printPercentage") file.readBool(printPercentage);
        if (tag == "percentageStop") file.readFloat(percentageStop);
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
    for (t = 0.0; t < tMax; t = t + dt) {
        if (step % printPercentageSteps == 0) 
            cerr << static_cast<int>(100 * t / tMax) << "% \r" << std::flush;
        step++;

        //saveVelocitiesToFile();
        if(100 * t / tMax > percentageStop) return ;

        setPBorders();
        setCavityFlowSpeeds();
        omp_set_dynamic(0);     // Explicitly disable dynamic teams
        omp_set_num_threads(4); // Use 4 threads for all consecutive parallel regions
        #pragma omp parallel for
        for (int i = 1; i < nX - 1; ++i) {
            for (int j = 1; j < nY - 1; ++j) {
                #ifdef USE_ASM
                    if(j > nY-5){
                        calcTerms(i,j);
                        calcVelocities(i,j);
                    }else{
                        calcVelocitiesAsm(i, j);
                        j += 3;
                    }
                    

                #endif
                #ifdef USE_CPP
                    calcTerms(i,j);
                    calcVelocities(i,j);
                #endif
            }
        }

        U0.copyAll(U1);
        U1.copyAll(U2);
        V0.copyAll(V1);
        V1.copyAll(V2);
        P0.copyAll(P1);
        P1.copyAll(P2);
    }
    cerr << 100 << "% \n";
}

void simulator::calcVelocities(int i, int j){
    float u2val = U1.at(i, j) 
                - U1.at(i, j) * (dt / dx) * (U1.at(i, j) - U1.at(i - 1, j)) 
                - V1.at(i, j) * (dt / dy) * (U1.at(i, j) - U1.at(i, j - 1)) 
                - (dt / (rho * 2 * dx)) * (P1.at(i + 1, j) - P1.at(i - 1, j))
                + nu * ((dt / (dx * dx)) * (U1.at(i + 1, j) - 2 * U1.at(i, j) + U1.at(i - 1, j)) 
                    + (dt / (dy * dy)) * (U1.at(i, j + 1) - 2 * U1.at(i, j) + U1.at(i, j - 1)) );
    U2.set(i, j, u2val);

    float v2val = V1.at(i, j) 
                - U1.at(i, j) * (dt / dx) * (V1.at(i, j) - V1.at(i - 1, j)) 
                - V1.at(i, j) * (dt / dy) * (V1.at(i, j) - V1.at(i, j - 1))
                - (dt / (rho * 2 * dy)) * (P1.at(i, j + 1) - P1.at(i, j - 1))
                + nu * ((dt / (dx * dx)) * (V1.at(i + 1, j) - 2 * V1.at(i, j) + V1.at(i - 1, j)) 
                    + (dt / (dy * dy)) * (V1.at(i, j + 1) - 2 * V1.at(i, j) + V1.at(i, j - 1)) );
    V2.set(i, j, v2val);

    float finalTerm = (1 / dt) * (U1x + V1y) - U1x * U1x - 2 * U1y * V1x - V1y * V1y;
    float p2val =   (1.0 / (2 * (dx * dx + dy * dy)))
                    * ((P1.at(i + 1, j) + P1.at(i - 1, j)) * dy * dy 
                        + (P1.at(i, j + 1) + P1.at(i, j - 1)) * dx * dx)
                    - (rho * dx * dx * dy * dy / (2 * dx * dx + 2 * dy * dy)) * finalTerm;
    P2.set(i, j, p2val);
}


#ifdef USE_ASM
void simulator::calcVelocitiesAsm(int i, int j) {

    float *mats[] = {U2.data, V2.data, P2.data, U1.data, V1.data, P1.data};
    int pos = i*nX+j;
    vvp_asm(mats, 4*pos, nX*4, dt, dx, dy, rho, nu);
}
#endif


void simulator::calcTerms(int i, int j) {
    U1x  = (U1.at(i + 1, j) - U1.at(i - 1, j)) / (2.0 * dx);
    U2x  = (U2.at(i + 1, j) - U2.at(i - 1, j)) / (2.0 * dx);
    U1y  = (U1.at(i, j + 1) - U1.at(i, j - 1)) / (2.0 * dy);
    U2y  = (U2.at(i, j + 1) - U2.at(i, j - 1)) / (2.0 * dy);
    U1xx = (U1.at(i + 1, j) - 2.0 * U1.at(i, j) + U1.at(i - 1, j)) / (dx * dx);
    U2xx = (U2.at(i + 1, j) - 2.0 * U2.at(i, j) + U2.at(i - 1, j)) / (dx * dx);
    U1yy = (U1.at(i, j + 1) - 2.0 * U1.at(i, j) + U1.at(i, j - 1)) / (dy * dy);
    U2yy = (U2.at(i, j + 1) - 2.0 * U2.at(i, j) + U2.at(i, j - 1)) / (dy * dy);
    P1x  = (P1.at(i + 1, j) - P1.at(i - 1, j)) / (2.0 * dx);
    P2x  = (P2.at(i + 1, j) - P2.at(i - 1, j)) / (2.0 * dx);
    P1y  = (P1.at(i, j + 1) - P1.at(i, j - 1)) / (2.0 * dy);
    P2y  = (P2.at(i, j + 1) - P2.at(i, j - 1)) / (2.0 * dy);
    V1x  = (V1.at(i + 1, j) - V1.at(i - 1, j)) / (2.0 * dx);
    V2x  = (V2.at(i + 1, j) - V2.at(i - 1, j)) / (2.0 * dx);
    V1y  = (V1.at(i, j + 1) - V1.at(i, j - 1)) / (2.0 * dy);
    V2y  = (V2.at(i, j + 1) - V2.at(i, j - 1)) / (2.0 * dy);
    V1xx = (V1.at(i + 1, j) - 2.0 * V1.at(i, j) + V1.at(i - 1, j)) / (dx * dx);
    V2xx = (V2.at(i + 1, j) - 2.0 * V2.at(i, j) + V2.at(i - 1, j)) / (dx * dx);
    V1yy = (V1.at(i, j + 1) - 2.0 * V1.at(i, j) + V1.at(i, j - 1)) / (dy * dy);
    V2yy = (V2.at(i, j + 1) - 2.0 * V2.at(i, j) + V2.at(i, j - 1)) / (dy * dy);
    P1xx = (P1.at(i + 1, j) - 2.0 * P1.at(i, j) + P1.at(i - 1, j)) / (dx * dx);
    P1yy = (P1.at(i, j + 1) - 2.0 * P1.at(i, j) + P1.at(i, j - 1)) / (dy * dy);
    P2yy = (P2.at(i, j + 1) - 2.0 * P2.at(i, j) + P2.at(i, j - 1)) / (dy * dy);
}


// ---------- Functions for file management ----------

void simulator::saveVtk(mat2 &m, string file_name) {
    // Save a 2-D scalar array in VTK format.
    io out(file_name, io::type_write);
    out.write("# vtk DataFile Version 2.0");
    out.newLine();
    out.write("Comment goes here");
    out.newLine();
    out.write("ASCII");
    out.newLine();
    out.write("DATASET STRUCTURED_POINTS");
    out.newLine();
    out.write("DIMENSIONS " + to_string(nX) + " " +  to_string(nY) + " 1");
    out.newLine();
    out.write("ORIGIN 0 0 0");
    out.newLine();
    out.write("SPACING " + to_string(dx) + " " + to_string(dy) + " 1");
    out.newLine();
    out.write("POINT_DATA " + to_string(nX * nY));
    out.newLine();
    out.write("SCALARS volume_scalars float 1");
    out.newLine();
    out.write("LOOKUP_TABLE default");
    out.newLine();
    for (int i = 0; i < nX; ++i) {
        for (int j = 0; j < nY; ++j) {
            out.writeFloat(m.at(i, j));
            out.write(" ");
        }
    }
    out.close();
}


void simulator::saveVelocitiesToFile(){
    ostringstream U_name;
    ostringstream V_name;
    ostringstream P_name;
    ostringstream norm_name;

    norm_name << "./out/norm_" << step << ".vtk";
    U_name << "./out/U_" << step << ".vtk";
    V_name << "./out/V_" << step << ".vtk";
    P_name << "./out/P_" << step << ".vtk";

    mat2 norm(U2.rows(), U2.cols());
    for (int i = 0; i < norm.rows(); ++i)
    {
        for (int j = 0; j < norm.cols(); ++j)
        {
            norm.set(i,j,sqrt(U2.at(i,j)*U2.at(i,j) + V2.at(i,j)*V2.at(i,j)));
        }
    }

    saveVtk(U2, U_name.str());
    saveVtk(V2, V_name.str());
    saveVtk(P2, P_name.str());
    saveVtk(norm, norm_name.str());

}


// ---------- Functions for testing ----------

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