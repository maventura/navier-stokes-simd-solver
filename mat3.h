#include <iostream>
using namespace std;

class mat3 {
    
    public:
    long double * data;
    int iMax;
    int jMax;
    int kMax;

    mat3() {
        iMax = 0;
        jMax = 0;
        kMax = 0;
        data = NULL;
    }

    mat3(int iMax, int jMax, int kMax) {
        this->iMax = iMax;
        this->jMax = jMax;
        this->kMax = kMax;
        inic(0);
    }


    mat3(int iMax, int jMax, int kMax, long double val) {
        this->iMax = iMax;
        this->jMax = jMax;
        this->kMax = kMax;
        inic(val);
    }

    void inic(long double val) {
        int n = iMax * jMax * kMax;
        data = new long double[n];

        for (int k = 0; k < n; k++)
            data[k] = val;
    }

    void set(int i, int j, int k, long double val) {
        data[iMax * jMax * k + i * jMax + j] = val;
    }

  //en este caso, setea una sola row de la k-esima matriz
  //TODO: ver si es esta la funcionalidad buscada
    void setRow(int i, int k, mat3 &B) {
        for (int j = 0; j < cols(); ++j)
            set(i, j, k, B.at(0, j, 0));
    }

    void setAll(mat3 &B) {
        for (int i = 0; i < rows(); ++i) {
            for (int j = 0; j < cols(); ++j) {
                for (int k = 0; k < dims(); ++k) {
                    set(i, j, k, B.at(i, j, k));
                }
            }
        }
    }


    void add(int i, int j, int k, long double val) {
        data[iMax * jMax * k + i * jMax + j] += val;
    }


    long double at(int i, int j, int k) {
        return data[iMax * jMax * k + i * jMax + j];
    }


    long double cols() {
        return jMax;
    }


    long double rows() {
        return iMax;
    }
    
    
    long double dims() {
        return kMax;
    }


    void print() {
        for(int k = 0; k < kMax; k++) {
            for (int i = 0; i < iMax; i++) {
                for (int j = 0; j < jMax; j++) {
                    cout << data[iMax * jMax * k + i * jMax + j] << " ";
                }
                cout << endl;
            }
            cout << endl;
            cout << endl;
        }
        cout << endl;
    }





    ~mat3() {
        if (data != NULL)
            delete [] data;
    }
};
