#ifndef simulator_H
#define simulator_H

#include <math.h>
#include "mat3.h"
#include <cmath>
#include "io.h"

extern "C" {void vvp_asm(float *mats[], int pos, float r, int h, float q, int offsetI, int offsetJ);}
                        
class simulator {

  public:
    simulator();
    void process();

  private:
    float xMax,  yMax,  zMax,  tMax;
    float nu, rho, C_d;
    float dx,  dy,  dz,  dt, h;
    int nX, nY, nZ, nT, step, maxIters;

    float al, fixedPointError, minFixedPointIters;
    float pi = atan(1) * 4;

    float W1y, V1z, W2y, V2z, U1z,
           W1x, U2z, W2x, V1x, U1y, V2x, U2y;

    int stepsUntilPrint, printPercentageSteps, maxSteps;
    bool printPercentage;
    float t;

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


    float Re, Fx, Fy, Fz;

    void centralSpeed();
    void centralPressure(int i, int j, int k);

    void calcTerms(int i, int j, int k);
    void saveVtk(mat3 &state, string file_name);
    void saveStreamVtk(mat3 &m1, mat3 &m2, mat3 &m3, string file_name);

    void setCavityFlowSpeeds();
    void setBorderConditions();
    void readParameters(string file_name);

    void simpleDiffusion(int i, int j, int k);
    void vorticityVectorPotencial(int i, int j, int k);
    void vorticityVectorPotencialAsm(int i, int j, int k);

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

        if (tag == "xMax") file.readFloat(xMax);
        if (tag == "yMax") file.readFloat(yMax);
        if (tag == "zMax") file.readFloat(zMax);
        if (tag == "tMax") file.readFloat(tMax);
        if (tag == "nu") file.readFloat(nu);
        if (tag == "rho") file.readFloat(rho);
        if (tag == "dx") file.readFloat(dx);
        if (tag == "dy") file.readFloat(dy);
        if (tag == "dz") file.readFloat(dz);
        if (tag == "dt") file.readFloat(dt);
        if (tag == "al") file.readFloat(al);
        if (tag == "fixedPointError") file.readFloat(fixedPointError);
        if (tag == "minFixedPointIters") file.readFloat(minFixedPointIters);
        if (tag == "printPercentageSteps") file.readInt(printPercentageSteps);
        if (tag == "C_d") file.readFloat(C_d);
        if (tag == "maxSteps") file.readInt(maxSteps);
        if (tag == "stepsUntilPrint") file.readInt(stepsUntilPrint);
        if (tag == "printPercentage") file.readBool(printPercentage);
        if (tag == "maxIters") file.readInt(maxIters);

    }
    file.close();
}


void simulator::setCavityFlowSpeeds(){
        for (int i = 0; i < nX; ++i) {
        for (int k = 0; k < nZ; ++k) {

            //Set Cavity flow conditions.
            U0.set(i, nY - 1, k, 0.01);
            U1.set(i, nY - 1, k, 0.01);
            U2.set(i, nY - 1, k, 0.01);
        }
    }

}
void simulator::setBorderConditions() {
    //TODO: Add if(cavityFlow) here and in the parameters.
    setCavityFlowSpeeds();

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
        setCavityFlowSpeeds();
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


        saveVtk(U2, U_name.str());
        saveVtk(V2, V_name.str());
        saveVtk(W2, W_name.str());

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
                        //runColorTest(i, j, k);
                        setVorticityVectorPotencialBorders(i, j, k);
                        bool inside = !(j == nY - 1  || i == nX - 1  || k == nZ - 1 || i * j * k == 0);
                        if(inside){ //TODO!! watafac paso aca
                            //vorticityVectorPotencialAsm(i,j,k);
                            //k += 3;
                            vorticityVectorPotencial(i,j,k);
                        }

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


void simulator::vorticityVectorPotencialAsm(int i, int j, int k) {
    /*
    mat_arr *mats;
    mats->U2 = U2.data;
    mats->V2 = V2.data;
    mats->W2 = W2.data;
    mats->omx2 = omx2.data;
    mats->omy2 = omy2.data;
    mats->omz2 = omz2.data;
    mats->psix2 = psix2.data;
    mats->psiy2 = psiy2.data;
    mats->psiz2 = psiz2.data;
    */
    float h = dx;
    float r = dt / (Re * h * h);
    float q = dt / (2 * h);
    int index = nZ * nY * i + j * nZ + k;
    float *mats[] = {U2.data, V2.data, W2.data, omx2.data, omy2.data, omz2.data, psix2.data, psiy2.data, psiz2.data, omx1.data, omy1.data, omz1.data};
    vvp_asm(mats, nZ * nY * i + j * nZ + k, r, h, q, nY * nZ, nZ);
}    


void simulator::vorticityVectorPotencial(int i, int j, int k) {
    /*for (int f = 0; f < nX; ++f) {
            for (int q = 0; q < nZ; ++q) {

                //Set Cavity flow conditions.
                U0.set(f, nY - 1, q, -0.01);
                U1.set(f, nY - 1, q, -0.01);
                U2.set(f, nY - 1, q, -0.01);
            }
        }
    */

    float h = dx; //TODO: ver que pasa con h si no es cuadrada la matriz o que
    float q = dt / (2 * h);
    float r = dt / (Re * h * h);
    float wt = 0.8;

    float aux_omx2 = omx2.at(i, j, k);
    float aux_omy2 = omy2.at(i, j, k);
    float aux_omz2 = omz2.at(i, j, k);
    float aux_psix2 = psix2.at(i, j, k);
    float aux_psiy2 = psiy2.at(i, j, k);
    float aux_psiz2 = psiz2.at(i, j, k);

    //Eje x.
    //calcular_V(i, j, k);
    float delta = (1 - q * U2.at(i + 1, j, k) + q * U2.at(i - 1, j, k) + 6 * r);
    float p1 = (-U2.at(i + 1, j, k) + r); 
    float p2 = (-V2.at(i, j + 1, k) + r);  
    float p3 = (U2.at(i - 1, j, k) + r);  
    float p4 = (V2.at(i, j - 1, k) + r);  
    float p5 = (W2.at(i, j, k + 1) + r);  
    float p6 = (W2.at(i, j, k - 1) + r);  
                                                        
    float ax1 = omy2.at(i, j, k) * U2.at(i, j + 1, k); 
    float ax2 = -omy2.at(i, j, k) * U2.at(i, j - 1, k);
    float ax3 = omz2.at(i, j, k) * U2.at(i, j, k + 1);
    float ax4 = -omz2.at(i, j, k) * U2.at(i, j, k - 1);

    float axs = ax1 + ax2 + ax3 + ax4;

    omx2.set(i, j, k, p1 * omx2.at(i + 1, j, k) + p2 * omx2.at(i, j + 1, k)
             + p3 * omx2.at(i - 1, j, k) + p4 * omx2.at(i, j - 1, k)
             + p5 * omx2.at(i, j, k + 1) + p6 * omx2.at(i, j, k - 1)
             + (1.0 / q)*omx1.at(i, j, k)                           
             + axs);
    omx2.set(i, j, k, omx2.at(i, j, k) * (q / delta)); 
    omx2.set(i, j, k, (1.0 - wt)*aux_omx2 + wt * omx2.at(i, j, k));


    //Eliptica eje x.
    //p1=p2=p3=p4=p5=p6=1/6.0;                                    
    psix2.set(i, j, k, (psix2.at(i + 1, j, k) + psix2.at(i, j + 1, k)
                        + psix2.at(i - 1, j, k) + psix2.at(i, j - 1, k)
                        + psix2.at(i, j, k + 1) + psix2.at(i, j, k - 1) 
                        + h * h * omx2.at(i, j, k)) / 6.0);
    psix2.set(i, j, k, (1.0 - wt)*aux_psix2 + wt * psix2.at(i, j, k));
 
    //Eje y.
    //calcular_V(i, j, k);
    p1 = -U2.at(i + 1, j, k) + r;
    p2 = -V2.at(i, j + 1, k) + r;
    p3 = U2.at(i - 1, j, k) + r;
    p4 = V2.at(i, j - 1, k) + r;
    p5 = W2.at(i, j, k + 1) + r;
    p6 = W2.at(i, j, k - 1) + r;

    delta = (1 - q * V2.at(i, j + 1, k) + q * V2.at(i, j - 1, k) + 6 * r);

    ax1 = omx2.at(i, j, k) * V2.at(i + 1, j, k);
    ax2 = -omx2.at(i, j, k) * V2.at(i - 1, j, k);
    ax3 = omz2.at(i, j, k) * V2.at(i, j, k + 1);
    ax4 = -omz2.at(i, j, k) * V2.at(i, j, k - 1);

    omy2.set(i, j, k, (p1) * omy2.at(i + 1, j, k) + (p2) * omy2.at(i, j + 1, k)
             + (p3) * omy2.at(i - 1, j, k) + (p4) * omy2.at(i, j - 1, k)
             + (p5) * omy2.at(i, j, k + 1) + (p6) * omy2.at(i, j, k - 1)
             + (1.0 / q)*omy1.at(i, j, k)
             + ax1  + ax2  + ax3  + ax4 );
    omy2.set(i, j, k, omy2.at(i, j, k) * (q / delta));
    omy2.set(i, j, k, (1.0 - wt)*aux_omy2 + wt * omy2.at(i, j, k));

    //Eliptica eje y.
    psiy2.set(i, j, k, (psiy2.at(i + 1, j, k) + psiy2.at(i, j + 1, k)
                        + psiy2.at(i - 1, j, k) + psiy2.at(i, j - 1, k)
                        + psiy2.at(i, j, k + 1) + psiy2.at(i, j, k - 1)
                        + h * h * omy2.at(i, j, k)) / 6.0);
    psiy2.set(i, j, k, (1.0 - wt)*aux_psiy2 + wt * psiy2.at(i, j, k));

    //Eje z.
    //calcular_V(i, j, k);
    delta = (1 - q * W2.at(i, j, k + 1) + q * W2.at(i, j, k - 1) + 6 * r);
    p1 = -U2.at(i + 1, j, k) + r;
    p2 = -V2.at(i, j + 1, k) + r;
    p3 = U2.at(i - 1, j, k) + r;
    p4 = V2.at(i, j - 1, k) + r;
    p5 = W2.at(i, j, k + 1) + r;
    p6 = W2.at(i, j, k - 1) + r;

    ax1 = omx2.at(i, j, k) * W2.at(i + 1, j, k);
    ax2 = -omx2.at(i, j, k) * W2.at(i - 1, j, k);
    ax3 = omy2.at(i, j, k) * W2.at(i, j + 1, k);
    ax4 = -omy2.at(i, j, k) * W2.at(i, j - 1, k);

    omz2.set(i, j, k, p1 * omz2.at(i + 1, j, k) + p2 * omz2.at(i, j + 1, k)
             + p3 * omz2.at(i - 1, j, k) + p4 * omz2.at(i, j - 1, k)
             + p5 * omz2.at(i, j, k + 1) + p6 * omz2.at(i, j, k - 1)
             + (1.0 / q)*omz1.at(i, j, k)
             + ax1 + ax2 + ax3 + ax4);
    omz2.set(i, j, k, omz2.at(i, j, k) * (q / delta));

    omz2.set(i, j, k, (1.0 - wt)*aux_omz2 + wt * omz2.at(i, j, k));

    //Eliptica eje z.
    //p1=p2=p3=p4=p5=p6=1/6.0;
    psiz2.set(i, j, k, (psiz2.at(i + 1, j, k) + psiz2.at(i, j + 1, k)
                        + psiz2.at(i - 1, j, k) + psiz2.at(i, j - 1, k)
                        + psiz2.at(i, j, k + 1) + psiz2.at(i, j, k - 1)
                        + h * h * omz2.at(i, j, k)) / 6.0);
    psiz2.set(i, j, k, (1.0 - wt)*aux_psiz2 + wt * psiz2.at(i, j, k));

    //calcular_V(i, j, k);

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
                out.writeFloat(m.at(i, j, k));
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
    U2.set(i, j, k,(psiz2.at(i, j + 1, k) - psiz2.at(i, j - 1, k)
            - psiy2.at(i, j, k + 1) + psiy2.at(i, j, k - 1) ) / 2 / h );
    V2.set(i, j, k, (psix2.at(i, j, k + 1) - psix2.at(i, j, k - 1) -
           psiz2.at(i + 1, j, k) + psiz2.at(i - 1, j, k)) / 2 / h);
    W2.set(i, j, k, (psiy2.at(i + 1, j, k) - psiy2.at(i - 1, j, k) -
           psix2.at(i, j + 1, k) + psix2.at(i, j - 1, k)) / 2 / h);
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


float minFloat(float a, float b) {
    if (a > b)return b;
    return a;
}


void simulator::centralColor() {

    for (int i = 0; i < nX; ++i) {
        for (int j = 0; j < nY; ++j) {
            for (int k = 0; k < nZ; ++k) {
                float xc = nX / 2;
                float yc = nY / 2;
                float zc = nZ / 2;
                float r = sqrt((i - xc) * (i - xc) + (j - yc) * (j - yc) + (k - zc) * (k - zc));


                if (r < 3.5) {
                    float max = 5;
                    color.set(i, j, k, minFloat(5 / r, 5));

                }
            }
        }
    }

}


void simulator::centralSpeed() {

    for (int i = 0; i < nX; ++i) {
        for (int j = 0; j < nY; ++j) {
            for (int k = 0; k < nZ; ++k) {
                float xc = nX / 2;
                float yc = nY / 2;
                float zc = nZ / 2;
                float r = sqrt((i - xc) * (i - xc) + (j - yc) * (j - yc) + (k - zc) * (k - zc));


                if (r < 2.5) {
                    float max = 0.01;
                    U2.set(i, j, k, minFloat(max / r, max));
                    U1.set(i, j, k, minFloat(max / r, max));
                    U0.set(i, j, k, minFloat(max / r, max));

                }
            }
        }
    }

}


void simulator::gridColor() {

    for (int i = 2; i < nX - 2; ++i) {
        for (int j = 2; j < nY - 2; ++j) {

            for (int k = 2; k < nZ - 2; ++k) {
                if (k % 10 == 3 && j % 10 == 3 && i % 10 == 3) {
                    //we want the first element non on the border to be painted
                    //so that edge behabiour is observable from the first step.
                    float color_intensity = 100;
                    color.set(i, j, k, color_intensity);
                    color.set(i + 1, j, k, color_intensity);
                    color.set(i - 1, j, k, color_intensity);
                    color.set(i, j + 1, k, color_intensity);
                    color.set(i, j - 1, k, color_intensity);
                    color.set(i, j, k + 1, color_intensity);
                    color.set(i, j, k - 1, color_intensity);

                    color.set(i, j + 1, k, color_intensity);
                    color.set(i + 1, j + 1, k, color_intensity);
                    color.set(i - 1, j + 1, k, color_intensity);
                    color.set(i, j + 1, k + 1, color_intensity);
                    color.set(i, j + 1, k - 1, color_intensity);

                    color.set(i, j - 1, k, color_intensity);
                    color.set(i + 1, j - 1, k, color_intensity);
                    color.set(i - 1, j - 1, k, color_intensity);
                    color.set(i, j - 1, k + 1, color_intensity);
                    color.set(i, j - 1, k - 1, color_intensity);


                    color.set(i + 1, j, k + 1, color_intensity);
                    color.set(i - 1, j, k + 1, color_intensity);
                    color.set(i, j + 1, k + 1, color_intensity);
                    color.set(i, j - 1, k + 1, color_intensity);

                    color.set(i + 1, j, k - 1, color_intensity);
                    color.set(i - 1, j, k - 1, color_intensity);
                    color.set(i, j + 1, k - 1, color_intensity);
                    color.set(i, j - 1, k - 1, color_intensity);

                    color.set(i + 1, j + 1, k, color_intensity);
                    color.set(i + 1, j - 1, k, color_intensity);
                    color.set(i + 1, j, k + 1, color_intensity);
                    color.set(i + 1, j, k - 1, color_intensity);

                    color.set(i - 1, j + 1, k, color_intensity);
                    color.set(i - 1, j - 1, k, color_intensity);
                    color.set(i - 1, j, k + 1, color_intensity);
                    color.set(i - 1, j, k - 1, color_intensity);

                }
            }
        }
    }

}


void simulator::runColorTest(int i, int j, int k) {
    bool border = (j >= nY - 1  || i >= nX - 1  || k >= nZ - 1 || i * j * k == 0);
    if (border) return;
    float multiplier = 20000;


    float c = 0;
    
    c += color.at(i, j, k);


    c += color.at(i - 1, j, k) * U2.at(i - 1, j, k) * multiplier;
    c += -color.at(i + 1, j, k) * U2.at(i + 1, j, k) * multiplier;
    c += color.at(i, j - 1, k) * V2.at(i, j - 1, k) * multiplier;
    c += -color.at(i, j + 1, k) * V2.at(i, j + 1, k) * multiplier;
    c += color.at(i, j, k - 1) * W2.at(i, j, k - 1) * multiplier;
    c += -color.at(i, j, k + 1) * W2.at(i, j, k + 1) * multiplier;
    color.set(i, j, k, c);
}


#endif
