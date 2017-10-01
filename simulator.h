#ifndef simulator_H
#define simulator_H

#include <math.h>
#include "mat3.h"
#include <cmath>
#include "io.h"

extern "C" {int vvp_asm(double *data, int pos);}

class simulator {

  public:
    simulator();
    void process();

  private:
    double xMax,  yMax,  zMax,  tMax;
    double nu, rho, C_d;
    double dx,  dy,  dz,  dt, h;
    int nX, nY, nZ, nT, step, maxIters;

    double al, fixedPointError, minFixedPointIters;
    double pi = atan(1) * 4;

    double W1y, V1z, W2y, V2z, U1z,
           W1x, U2z, W2x, V1x, U1y, V2x, U2y;

    int stepsUntilPrint, printPercentageSteps, maxSteps;
    bool printPercentage;
    double t;

    mat3 U0, V0, W0, P0;
    mat3 U1, V1, W1, P1;
    mat3 U2, V2, W2, P2;

    mat3 omx0, omy0, omz0;
    mat3 omx1, omy1, omz1;
    mat3 omx2, omy2, omz2;

    mat3 psix0, psiy0, psiz0;
    mat3 psix1, psiy1, psiz1;
    mat3 psix2, psiy2, psiz2;

    mat3 U_aux_0, V_aux_0, W_aux_0;
    mat3 U_aux_1, V_aux_1, W_aux_1;
    mat3 U_aux_2, V_aux_2, W_aux_2;

    mat3 color;


    double Re, Fx, Fy, Fz;

    void centralSpeed();
    void centralPressure(int i, int j, int k);

    void calcTerms(int i, int j, int k);
    void saveVtk(mat3 &state, string file_name);
    void saveStreamVtk(mat3 &m1, mat3 &m2, mat3 &m3, string file_name);


    void setBorderConditions();
    void readParameters(string file_name);

    void simpleDiffusion(int i, int j, int k);
    void vorticityVectorPotencial(int i, int j, int k);
    void setVorticityVectorPotencialBorders(int i, int j, int k);
    void calcular_V(int i, int j, int k);

    void centralColor();
    void gridColor();
    void runColorTest(int i, int j, int k);

};

simulator::simulator() {

    readParameters("parameters.txt");
    h = dx; //for square grids.
    nX = round(xMax / dx) + 1;
    nY = round(yMax / dy) + 1;
    nZ = round(yMax / dy) + 1;
    nT = round(tMax / dt) + 1;
    Re = (xMax * 0.01) / (nu); //TODO: 0.01 is the inicial velocity at y = yMax

    U0.confAndInit(nX, nY, nZ, 0);
    V0.confAndInit(nX, nY, nZ, 0);
    W0.confAndInit(nX, nY, nZ, 0);
    P0.confAndInit(nX, nY, nZ, 1);

    U1.confAndInit(nX, nY, nZ, 0);
    V1.confAndInit(nX, nY, nZ, 0);
    W1.confAndInit(nX, nY, nZ, 0);
    P1.confAndInit(nX, nY, nZ, 1);

    U2.confAndInit(nX, nY, nZ, 0);
    V2.confAndInit(nX, nY, nZ, 0);
    W2.confAndInit(nX, nY, nZ, 0);
    P2.confAndInit(nX, nY, nZ, 1);

    psix0.confAndInit(nX, nY, nZ, 0);
    psix1.confAndInit(nX, nY, nZ, 0);
    psix2.confAndInit(nX, nY, nZ, 0);
    psiy0.confAndInit(nX, nY, nZ, 0);
    psiy1.confAndInit(nX, nY, nZ, 0);
    psiy2.confAndInit(nX, nY, nZ, 0);
    psiz0.confAndInit(nX, nY, nZ, 0);
    psiz1.confAndInit(nX, nY, nZ, 0);
    psiz2.confAndInit(nX, nY, nZ, 0);

    omx0.confAndInit(nX, nY, nZ, 0);
    omx1.confAndInit(nX, nY, nZ, 0);
    omx2.confAndInit(nX, nY, nZ, 0);
    omy0.confAndInit(nX, nY, nZ, 0);
    omy1.confAndInit(nX, nY, nZ, 0);
    omy2.confAndInit(nX, nY, nZ, 0);
    omz0.confAndInit(nX, nY, nZ, 0);
    omz1.confAndInit(nX, nY, nZ, 0);
    omz2.confAndInit(nX, nY, nZ, 0);

    U_aux_0.confAndInit(nX, nY, nZ, 0);
    V_aux_0.confAndInit(nX, nY, nZ, 0);
    W_aux_0.confAndInit(nX, nY, nZ, 0);
    U_aux_1.confAndInit(nX, nY, nZ, 0);
    V_aux_1.confAndInit(nX, nY, nZ, 0);
    W_aux_1.confAndInit(nX, nY, nZ, 0);
    U_aux_2.confAndInit(nX, nY, nZ, 0);
    V_aux_2.confAndInit(nX, nY, nZ, 0);
    W_aux_2.confAndInit(nX, nY, nZ, 0);


    color.confAndInit(nX, nY, nZ, 0);

    setBorderConditions();
}

void simulator::readParameters(string file_name) {
    io file(file_name, io::type_read);
    string tag;
    while (file.readWord(tag)) {

        if (tag == "xMax") file.readDouble(xMax);
        if (tag == "yMax") file.readDouble(yMax);
        if (tag == "zMax") file.readDouble(zMax);
        if (tag == "tMax") file.readDouble(tMax);
        if (tag == "nu") file.readDouble(nu);
        if (tag == "rho") file.readDouble(rho);
        if (tag == "dx") file.readDouble(dx);
        if (tag == "dy") file.readDouble(dy);
        if (tag == "dz") file.readDouble(dz);
        if (tag == "dt") file.readDouble(dt);
        if (tag == "al") file.readDouble(al);
        if (tag == "fixedPointError") file.readDouble(fixedPointError);
        if (tag == "minFixedPointIters") file.readDouble(minFixedPointIters);
        if (tag == "printPercentageSteps") file.readInt(printPercentageSteps);
        if (tag == "C_d") file.readDouble(C_d);
        if (tag == "maxSteps") file.readInt(maxSteps);
        if (tag == "stepsUntilPrint") file.readInt(stepsUntilPrint);
        if (tag == "printPercentage") file.readBool(printPercentage);
        if (tag == "maxIters") file.readInt(maxIters);

    }
    file.close();
}


void simulator::setBorderConditions() {
    //TODO: Add if(cavityFlow) here and in the parameters.
    for (int i = 0; i < nX; ++i) {
        for (int k = 0; k < nZ; ++k) {

            //Set Cavity flow conditions.
            U0.set(i, nY - 1, k, 0.01);
            U1.set(i, nY - 1, k, 0.01);
            U2.set(i, nY - 1, k, 0.01);
        }
    }

    //Vorticity conditions:
    for (int i = 1; i < nX - 1; ++i) {
        for (int j = 1; j < nY - 1; ++j) {
            for (int k = 1; k < nZ - 1; ++k) {
                calcTerms(i, j, k);
                omz1.set(i, j, k, W1y - V1z);
                omz2.set(i, j, k, W2y - V2z);

                omy1.set(i, j, k, U1z - W1x);
                omy2.set(i, j, k, U2z - W2x);

                omx1.set(i, j, k, V1x - U1y);
                omx2.set(i, j, k, V2x - U2y);
            }
        }
    }
}


void simulator::process() {
    step = 0;
    //centralColor();
    gridColor();
    for (t = 0; t < tMax; t = t + dt) {
        cerr << 100 * t / tMax << "%" << endl;
        step++;
        ostringstream U_name;
        ostringstream V_name;
        ostringstream W_name;
        ostringstream P_name;
        U_name << "./out/U_" << step << ".vtk";
        ostringstream V0_name;
        V_name << "./out/V_" << step << ".vtk";
        ostringstream W0_name;
        W_name << "./out/W_" << step << ".vtk";


        ostringstream psix0_name;
        psix0_name << "./out/psix0_" << step << ".vtk";
        ostringstream omx0_name;
        omx0_name << "./out/omx0_" << step << ".vtk";

        ostringstream psiy0_name;
        psiy0_name << "./out/psiy0_" << step << ".vtk";
        ostringstream omy0_name;
        omy0_name << "./out/omy0_" << step << ".vtk";

        ostringstream psiz0_name;
        psiz0_name << "./out/psiz0_" << step << ".vtk";
        ostringstream omz0_name;
        omz0_name << "./out/omz0_" << step << ".vtk";


        ostringstream stream_name;
        stream_name << "./out/stream_" << step << ".vtk";

      //  if(step%5 == 0){
        ostringstream color_name;
        color_name << "./out/color_" << step << ".vtk";
        saveVtk(color, color_name.str());
        //}


        //  saveVtk(U2, U_name.str());
        // saveVtk(V2, V_name.str());
        // saveVtk(W2, W_name.str());

        // saveVtk(psix0, psix0_name.str());
        // saveVtk(psiy0, psiy0_name.str());
         //saveVtk(psiz0, psiz0_name.str());

        // saveVtk(omx0, omx0_name.str());
        // saveVtk(omy0, omy0_name.str());
        //saveVtk(omz0, omz0_name.str());


        //saveStreamVtk(U2, V2, W2, stream_name.str());

        for (int iter = 0; iter < maxIters; ++iter) {
            for (int i = 0; i < nX; ++i) {
                for (int j = 0; j < nY; ++j) {
                    //if(j == nY-1) continue;
                    for (int k = 0; k < nZ; ++k) {
                        runColorTest(i, j, k);
                        setVorticityVectorPotencialBorders(i, j, k);
                        bool inside = !(j == nY - 1  || i == nX - 1  || k == nZ - 1 || i * j * k == 0);
                        U1.setAll(1);
                        vvp_asm(U1.data, nX * nY * k + i * nY + j);


                        U1.print();
                        cout << U1.at(0,0,0) << endl;


                    }
                }
            }
        }

        //TODO: actualizar las velocidades acá también genera transiciones bruscas entre t y t+1
        U1.copyAll(U2);
        V1.copyAll(V2);
        W1.copyAll(W2);

        psix1.copyAll(psix2);
        psiy1.copyAll(psiy2);
        psiz1.copyAll(psiz2);

        omx1.copyAll(omx2);
        omy1.copyAll(omy2);
        omz1.copyAll(omz2);

        //TODO:mas bien reseña, establecer a las px_2 en cero no esta bueno, genera resultados con transiciones menos suaves.
    }
    cerr << "Message: Processing finished correctly" << endl << flush;
}

void simulator::vorticityVectorPotencial(int i, int j, int k) {
    /*
        for (int f = 0; f < nX; ++f) {
            for (int q = 0; q < nZ; ++q) {

                //Set Cavfty flow condftfons.
                U0.set(f, nY - 1, q, -0.01);
                U1.set(f, nY - 1, q, -0.01);
                U2.set(f, nY - 1, q, -0.01);
            }
        }

    */
    double h = dx;
    double q = dt / (2 * h);
    double r = dt / (Re * h * h);
    double wt = 0.8;
    double aux_omx2 = omx2.at(i, j, k);
    double aux_omy2 = omy2.at(i, j, k);
    double aux_omz2 = omz2.at(i, j, k);
    double aux_psix2 = psix2.at(i, j, k);
    double aux_psiy2 = psiy2.at(i, j, k);
    double aux_psiz2 = psiz2.at(i, j, k);

    //Eje x.
    calcular_V(i, j, k);
    double delta = (1 - q * U2.at(i + 1, j, k) + q * U2.at(i - 1, j, k) + 6 * r);
    double p1 = (-q * U2.at(i + 1, j, k) + r) / delta;
    double p2 = (-q * V2.at(i, j + 1, k) + r) / delta;
    double p3 = (q * U2.at(i - 1, j, k) + r) / delta;
    double p4 = (q * V2.at(i, j - 1, k) + r) / delta;
    double p5 = (q * W2.at(i, j, k + 1) + r) / delta;
    double p6 = (q * W2.at(i, j, k - 1) + r) / delta;

    double ax1 = q * omy2.at(i, j, k) * U2.at(i, j + 1, k);
    double ax2 = -q * omy2.at(i, j, k) * U2.at(i, j - 1, k);
    double ax3 = q * omz2.at(i, j, k) * U2.at(i, j, k + 1);
    double ax4 = -q * omz2.at(i, j, k) * U2.at(i, j, k - 1);

    omx2.set(i, j, k, p1 * omx2.at(i + 1, j, k) + p2 * omx2.at(i, j + 1, k)
             + p3 * omx2.at(i - 1, j, k) + p4 * omx2.at(i, j - 1, k)
             + p5 * omx2.at(i, j, k + 1) + p6 * omx2.at(i, j, k - 1)
             + omx1.at(i, j, k) / delta
             + ax1 / delta + ax2 / delta + ax3 / delta + ax4 / delta);

    omx2.set(i, j, k, (1.0 - wt)*aux_omx2 + wt * omx2.at(i, j, k));

    //Eliptica eje x.
    //p1=p2=p3=p4=p5=p6=1/6.0;
    psix2.set(i, j, k, (psix2.at(i + 1, j, k) + psix2.at(i, j + 1, k)
                        + psix2.at(i - 1, j, k) + psix2.at(i, j - 1, k)
                        + psix2.at(i, j, k + 1) + psix2.at(i, j, k - 1)
                        + h * h * omx2.at(i, j, k)) / 6.0);
    psix2.set(i, j, k, (1.0 - wt)*aux_psix2 + wt * psix2.at(i, j, k));

    //Eje y.
    calcular_V(i, j, k);
    delta = (1 - q * V2.at(i, j + 1, k) + q * V2.at(i, j - 1, k) + 6 * r);
    p1 = -q * U2.at(i + 1, j, k) + r;
    p2 = -q * V2.at(i, j + 1, k) + r;
    p3 = q * U2.at(i - 1, j, k) + r;
    p4 = q * V2.at(i, j - 1, k) + r;
    p5 = q * W2.at(i, j, k + 1) + r;
    p6 = q * W2.at(i, j, k - 1) + r;
    p1 = p1 / delta;
    p2 = p2 / delta;
    p3 = p3 / delta;
    p4 = p4 / delta;
    p5 = p5 / delta;
    p6 = p6 / delta;

    ax1 = q * omx2.at(i, j, k) * V2.at(i + 1, j, k);
    ax2 = -q * omx2.at(i, j, k) * V2.at(i - 1, j, k);
    ax3 = q * omz2.at(i, j, k) * V2.at(i, j, k + 1);
    ax4 = -q * omz2.at(i, j, k) * V2.at(i, j, k - 1);

    omy2.set(i, j, k, p1 * omy2.at(i + 1, j, k) + p2 * omy2.at(i, j + 1, k)
             + p3 * omy2.at(i - 1, j, k) + p4 * omy2.at(i, j - 1, k)
             + p5 * omy2.at(i, j, k + 1) + p6 * omy2.at(i, j, k - 1)
             + omy1.at(i, j, k) / delta
             + ax1 / delta + ax2 / delta + ax3 / delta + ax4 / delta);

    omy2.set(i, j, k, (1.0 - wt)*aux_omy2 + wt * omy2.at(i, j, k));

    //Eliptica eje y.
    psiy2.set(i, j, k, (psiy2.at(i + 1, j, k) + psiy2.at(i, j + 1, k)
                        + psiy2.at(i - 1, j, k) + psiy2.at(i, j - 1, k)
                        + psiy2.at(i, j, k + 1) + psiy2.at(i, j, k - 1)
                        + h * h * omy2.at(i, j, k)) / 6.0);
    psiy2.set(i, j, k, (1.0 - wt)*aux_psiy2 + wt * psiy2.at(i, j, k));

    //Eje z.
    calcular_V(i, j, k);
    delta = (1 - q * W2.at(i, j, k + 1) + q * W2.at(i, j, k - 1) + 6 * r);
    p1 = -q * U2.at(i + 1, j, k) + r;
    p2 = -q * V2.at(i, j + 1, k) + r;
    p3 = q * U2.at(i - 1, j, k) + r;
    p4 = q * V2.at(i, j - 1, k) + r;
    p5 = q * W2.at(i, j, k + 1) + r;
    p6 = q * W2.at(i, j, k - 1) + r;
    p1 = p1 / delta;
    p2 = p2 / delta;
    p3 = p3 / delta;
    p4 = p4 / delta;
    p5 = p5 / delta;
    p6 = p6 / delta;

    ax1 = q * omx2.at(i, j, k) * W2.at(i + 1, j, k);
    ax2 = -q * omx2.at(i, j, k) * W2.at(i - 1, j, k);
    ax3 = q * omy2.at(i, j, k) * W2.at(i, j + 1, k);
    ax4 = -q * omy2.at(i, j, k) * W2.at(i, j - 1, k);

    omz2.set(i, j, k, p1 * omz2.at(i + 1, j, k) + p2 * omz2.at(i, j + 1, k)
             + p3 * omz2.at(i - 1, j, k) + p4 * omz2.at(i, j - 1, k)
             + p5 * omz2.at(i, j, k + 1) + p6 * omz2.at(i, j, k - 1)
             + omz1.at(i, j, k) / delta
             + ax1 / delta + ax2 / delta + ax3 / delta + ax4 / delta);

    omz2.set(i, j, k, (1.0 - wt)*aux_omz2 + wt * omz2.at(i, j, k));

    //Eliptica eje z.
    //p1=p2=p3=p4=p5=p6=1/6.0;
    psiz2.set(i, j, k, (psiz2.at(i + 1, j, k) + psiz2.at(i, j + 1, k)
                        + psiz2.at(i - 1, j, k) + psiz2.at(i, j - 1, k)
                        + psiz2.at(i, j, k + 1) + psiz2.at(i, j, k - 1)
                        + h * h * omz2.at(i, j, k)) / 6.0);
    psiz2.set(i, j, k, (1.0 - wt)*aux_psiz2 + wt * psiz2.at(i, j, k));

    calcular_V(i, j, k);

}

void simulator::simpleDiffusion(int i, int j, int k) {
    //simple diffusion for testing.
    U2.set(i, j, k, (U1.at(i, j, k) + U1.at(i + 1, j, k) + U1.at(i, j + 1, k) + U1.at(i, j, k + 1) + U1.at(i - 1, j, k) + U1.at(i, j - 1, k) + U1.at(i, j, k - 1)) / 7.0 );
    V2.set(i, j, k, (V1.at(i, j, k) + V1.at(i + 1, j, k) + V1.at(i, j + 1, k) + V1.at(i, j, k + 1) + V1.at(i - 1, j, k) + V1.at(i, j - 1, k) + V1.at(i, j, k - 1)) / 7.0 );
    W2.set(i, j, k, (W1.at(i, j, k) + W1.at(i + 1, j, k) + W1.at(i, j + 1, k) + W1.at(i, j, k + 1) + W1.at(i - 1, j, k) + W1.at(i, j - 1, k) + W1.at(i, j, k - 1)) / 7.0 );

    if (i == 10 && j == 10 && k == 10) U2.set(i, j, k, 0.1);
    if (i == 10 && j == 10 && k == 10) V2.set(i, j, k, 0.1);
    if (i == 10 && j == 10 && k == 10) W2.set(i, j, k, 0.1);

}

void simulator::saveStreamVtk(mat3 &m1, mat3 &m2, mat3 &m3, string file_name) {
    // Save a 3-D scalar array in VTK format.
    io out(file_name, io::type_write);
    out.write("# vtk DataFile Version 2.0");
    out.newLine();
    out.write("Comment goes here");
    out.newLine();
    out.write("ASCII");
    out.newLine();
    out.newLine();

    out.write("DATASET RECTILINEAR_GRID");
    out.newLine();
    out.write("DIMENSIONS    " + to_string(nX) + " " +  to_string(nY) + " " + to_string(nZ));
    out.newLine();
    out.newLine();

    out.write("X_COORDINATES " + to_string(nX) + " float");
    out.newLine();
//for (int i = 0; i < nX; ++i) out.write(to_string(i) + " ");
//    out.newLine();

    out.write("Y_COORDINATES " + to_string(nY) + " float");
    out.newLine();
//for (int i = 0; i < nY; ++i) out.write(to_string(i) + " ");
//    out.newLine();

    out.write("Z_COORDINATES " + to_string(nZ) + " float");
    out.newLine();
//for (int i = 0; i < nZ; ++i) out.write(to_string(i) + " ");
//    out.newLine();

    out.write("POINT_DATA " + to_string(nX * nY * nZ));
    out.newLine();

    out.write("SCALARS Xvelocity float 1");
    out.newLine();
    out.write("LOOKUP_TABLE default");
    out.newLine();
    for (int i = 0; i < nX * nY * nZ; ++i) {
        out.write(to_string(i) + ".0");
        out.newLine();
    }


    out.write("SCALARS Yvelocity float 1");
    out.newLine();
    out.write("LOOKUP_TABLE default");
    out.newLine();
    for (int i = 0; i < nX * nY * nZ; ++i) {
        out.write(to_string(i) + ".0");
        out.newLine();
    }


    out.write("SCALARS Zvelocity float 1");
    out.newLine();
    out.write("LOOKUP_TABLE default");
    out.newLine();
    for (int i = 0; i < nX * nY * nZ; ++i) {
        out.write(to_string(i) + ".0");
        out.newLine();
    }

    out.write("VECTORS VecVelocity float");


    for (int i = 0; i < nX; ++i) {
        for (int j = 0; j < nY; ++j) {
            for (int k = 0; k < nZ; ++k) {
                out.write(to_string(m1.at(i, j, k)) + " " + to_string(m2.at(i, j, k))  + " " + to_string(m3.at(i, j, k)));
                out.newLine();
            }
        }
    }
    out.close();
}



void simulator::saveVtk(mat3 &m, string file_name) {
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
    out.write("POINT_DATA   " + to_string(nX * nY * nZ));
    out.newLine();
    out.write("SCALARS scalars float");
    out.newLine();
    out.write("LOOKUP_TABLE default");
    out.newLine();
    out.newLine();

    for (int i = 0; i < nX; ++i) {
        for (int j = 0; j < nY; ++j) {
            for (int k = 0; k < nZ; ++k) {
                out.writeDouble(m.at(i, j, k));
                out.write(" ");
            }
            out.newLine();
        }
    }
    out.close();
}


void simulator::setVorticityVectorPotencialBorders(int i, int j, int k) {
//Condiciones de borde.
    if (j == 0) {
        psix2.set(i, j, k, 0);
        //Se agrega nueva condicion con la derivada normal en dos puntos.
        psiy2.set(i, j, k,  (4 * psiy2.at(i, j + 1, k) - psiy2.at(i, j + 2, k)) / 3);
        psiz2.set(i, j, k, 0);

        U2.set(i, j, k, 0);
        V2.set(i, j, k, 0);
        W2.set(i, j, k, 0);
        //Se abusa de que la velocidad en el punto es 0.
        omx2.set(i, j, k,  (W2.at(i, j + 1, k)) / h);
        omy2.set(i, j, k, 0);
        omz2.set(i, j, k,  -(U2.at(i, j + 1, k)) / h);

    } else if (j == nY - 1) {
        //Anodo.
        psix2.set(i, j, k, 0);
        //Se agrega nueva condicion con la derivada normal en dos puntos.
        psiy2.set(i, j, k,  (4 * psiy2.at(i, j - 1, k) - psiy2.at(i, j - 2, k)) / 3);
        psiz2.set(i, j, k, 0);

        U2.set(i, j, k, 0);
        V2.set(i, j, k, 0);
        W2.set(i, j, k, 0);
        //Se abusa de que la velocidad en el punto es 0.
        omx2.set(i, j, k,  (-W2.at(i, j - 1, k)) / h);
        omy2.set(i, j, k, 0);
        omz2.set(i, j, k,  (U2.at(i, j - 1, k)) / h);

    } else if (i == 0) {
        //Pared Lateral.
        //Se agrega nueva condici�n con la derivada normal en dos puntos.
        psix2.set(i, j, k,  (4 * psix2.at(i + 1, j, k) - psix2.at(i + 2, j, k)) / 3);
        psiy2.set(i, j, k, 0);
        psiz2.set(i, j, k, 0);
        //10/06/2002: Se agrega fijar las velocidades de
        //este punto en 0.
        U2.set(i, j, k, 0);
        V2.set(i, j, k, 0);
        W2.set(i, j, k, 0);
        //Se abusa de que la velocidad en el punto es 0.
        omx2.set(i, j, k, 0);
        omy2.set(i, j, k, (-W2.at(i + 1, j, k)) / h);
        omz2.set(i, j, k, (V2.at(i + 1, j, k)) / h);

    } else if (i == nX - 1) {
        //Pared Lateral.
        //Se agrega nueva condicion con la derivada normal en dos puntos.
        psix2.set(i, j, k,  (4 * psix2.at(i - 1, j, k) - psix2.at(i - 2, j, k)) / 3);
        psiy2.set(i, j, k, 0);
        psiz2.set(i, j, k, 0);
        //10/06/2002: Se agrega fijar las velocidades de
        //este punto en 0.
        U2.set(i, j, k, 0);
        V2.set(i, j, k, 0);
        W2.set(i, j, k, 0);
        //Se abusa de que la velocidad en el punto es 0.
        omx2.set(i, j, k, 0);
        omy2.set(i, j, k, ( W2.at(i - 1, j, k)) / h);
        omz2.set(i, j, k, (-V2.at(i - 1, j, k)) / h);

    } else if (k == 0) {
        //Piso.
        psix2.set(i, j, k, 0);
        psiy2.set(i, j, k, 0);
        //Se agrega nueva condicion con la derivada normal en dos puntos.
        psiz2.set(i, j, k,  (4 * psiz2.at(i, j, k + 1) - psiz2.at(i, j, k + 2)) / 3);

        U2.set(i, j, k, 0);
        V2.set(i, j, k, 0);
        W2.set(i, j, k, 0);
        //Se abusa de que la velocidad en el punto es 0.
        omx2.set(i, j, k, (-V2.at(i, j, k + 1)) / h);
        omy2.set(i, j, k, (U2.at(i, j, k + 1)) / h);
        omz2.set(i, j, k, 0);

    } else if (k == nZ - 1) {
        //Techo.
        psix2.set(i, j, k, 0);
        psiy2.set(i, j, k, 0);
        //Se agrega nueva condicion con la derivada normal en dos puntos.
        psiz2.set(i, j, k,  (4 * psiz2.at(i, j, k - 1) - psiz2.at(i, j, k - 2)) / 3);

        U2.set(i, j, k, 0);
        V2.set(i, j, k, 0);
        W2.set(i, j, k, 0);
        //Se abusa de que la velocidad en el punto es 0.
        omx2.set(i, j, k, (V2.at(i, j, k - 1)) / h);
        omy2.set(i, j, k, (-U2.at(i, j, k - 1)) / h);
        omz2.set(i, j, k, 0);

    }

}

void simulator::calcular_V(int i, int j, int k) {
    U2.set(i, j, k, (psiz2.at(i, j + 1, k) - psiz2.at(i, j - 1, k)) / 2 / h -
           (psiy2.at(i, j, k + 1) - psiy2.at(i, j, k - 1)) / 2 / h);
    V2.set(i, j, k, (psix2.at(i, j, k + 1) - psix2.at(i, j, k - 1)) / 2 / h -
           (psiz2.at(i + 1, j, k) - psiz2.at(i - 1, j, k)) / 2 / h);
    W2.set(i, j, k, (psiy2.at(i + 1, j, k) - psiy2.at(i - 1, j, k)) / 2 / h -
           (psix2.at(i, j + 1, k) - psix2.at(i, j - 1, k)) / 2 / h);
}

void simulator::calcTerms(int i, int j, int k) {

    U1y = (U1.at(i, j + 1, k) - U1.at(i, j - 1, k)) / (2 * dy);
    U2y = (U2.at(i, j + 1, k) - U2.at(i, j - 1, k)) / (2 * dy);
    U1z = (U1.at(i, j, k + 1) - U1.at(i, j, k - 1)) / (2 * dz);
    U2z = (U2.at(i, j, k + 1) - U2.at(i, j, k - 1)) / (2 * dz);

    V1x = (V1.at(i + 1, j, k) - V1.at(i - 1, j, k)) / (2 * dx);
    V2x = (V2.at(i + 1, j, k) - V2.at(i - 1, j, k)) / (2 * dx);

    V1z = (V1.at(i, j, k + 1) - V1.at(i, j, k - 1)) / (2 * dz);
    V2z = (V2.at(i, j, k + 1) - V2.at(i, j, k - 1)) / (2 * dz);

    W1x = (W1.at(i + 1, j, k) - W1.at(i - 1, j, k)) / (2 * dx);
    W2x = (W2.at(i + 1, j, k) - W2.at(i - 1, j, k)) / (2 * dx);
    W1y = (W1.at(i, j + 1, k) - W1.at(i, j - 1, k)) / (2 * dy);
    W2y = (W2.at(i, j + 1, k) - W2.at(i, j - 1, k)) / (2 * dy);

}

double minDouble(double a, double b){
  if(a > b)return b;
  return a;
}

void simulator::centralColor() {

    for (int i = 0; i < nX; ++i) {
        for (int j = 0; j < nY; ++j) {
            for (int k = 0; k < nZ; ++k) {
                double xc = nX / 2;
                double yc = nY / 2;
                double zc = nZ / 2;
                double r = sqrt((i - xc) * (i - xc) + (j - yc) * (j - yc) + (k - zc) * (k - zc));


                if (r < 3.5) {
                  double max = 5;
                    color.set(i, j, k, minDouble(5/r,5));

                }
            }
        }
    }

}



void simulator::centralSpeed() {

    for (int i = 0; i < nX; ++i) {
        for (int j = 0; j < nY; ++j) {
            for (int k = 0; k < nZ; ++k) {
                double xc = nX / 2;
                double yc = nY / 2;
                double zc = nZ / 2;
                double r = sqrt((i - xc) * (i - xc) + (j - yc) * (j - yc) + (k - zc) * (k - zc));


                if (r < 2.5) {
                  double max = 0.01;
                    U2.set(i, j, k, minDouble(max/r,max));
                    U1.set(i, j, k, minDouble(max/r,max));
                    U0.set(i, j, k, minDouble(max/r,max));

                }
            }
        }
    }

}




void simulator::gridColor() {

    for (int i = 3; i < nX-2; ++i) {
        for (int j = 3; j < nY-2; ++j) {

            for (int k = 3; k < nZ-2; ++k) {
              if(k%20==0 && j%20==0 && i%20==0) {
                  color.set(i, j, k, 5);
                  color.set(i+1, j, k, 5);
                  color.set(i-1, j, k, 5);
                  color.set(i, j+1, k, 5);
                  color.set(i, j-1, k, 5);
                  color.set(i, j, k+1, 5);
                  color.set(i, j, k-1, 5);

                  color.set(i, j+1, k, 5);
                  color.set(i+1, j+1, k, 5);
                  color.set(i-1, j+1, k, 5);
                  color.set(i, j+1, k+1, 5);
                  color.set(i, j+1, k-1, 5);

                  color.set(i, j-1, k, 5);
                  color.set(i+1, j-1, k, 5);
                  color.set(i-1, j-1, k, 5);
                  color.set(i, j-1, k+1, 5);
                  color.set(i, j-1, k-1, 5);


                  color.set(i+1, j, k+1, 5);
                  color.set(i-1, j, k+1, 5);
                  color.set(i, j+1, k+1, 5);
                  color.set(i, j-1, k+1, 5);

                  color.set(i+1, j, k-1, 5);
                  color.set(i-1, j, k-1, 5);
                  color.set(i, j+1, k-1, 5);
                  color.set(i, j-1, k-1, 5);

                  color.set(i+1, j+1, k, 5);
                  color.set(i+1, j-1, k, 5);
                  color.set(i+1, j, k+1, 5);
                  color.set(i+1, j, k-1, 5);

                  color.set(i-1, j+1, k, 5);
                  color.set(i-1, j-1, k, 5);
                  color.set(i-1, j, k+1, 5);
                  color.set(i-1, j, k-1, 5);

                }
            }
        }
    }

}

void simulator::runColorTest(int i, int j, int k) {
    bool border = (j >= nY - 1  || i >= nX - 1  || k >= nZ - 1 || i * j * k == 0);
    if (border) return;

    double c = 0;
    c += color.at(i, j, k);
    c += color.at(i - 1, j, k) * U2.at(i - 1, j, k)*2000;
    c += -color.at(i + 1, j, k) * U2.at(i + 1, j, k)*2000;
    c += color.at(i, j - 1, k) * V2.at(i, j - 1, k)*2000;
    c += -color.at(i, j + 1, k) * V2.at(i, j + 1, k)*2000;
    c += color.at(i, j, k - 1) * W2.at(i, j, k - 1)*2000;
    c += -color.at(i, j, k + 1) * W2.at(i, j, k + 1)*2000;
    color.set(i, j, k, c);
}


#endif
