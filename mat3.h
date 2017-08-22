#ifndef MAT3_H
#define MAT3_H

#include <iostream>
using namespace std;

class mat3 {

public:
    mat3();
    ~mat3();
    mat3(int i_max_, int j_max_, int k_max_);
    mat3(int i_max_, int j_max_, int k_max_, double val);
    void confAndInit(int i_max_, int j_max_, int k_max_, double val);
    //TODO: No control over multiple allocations using confAndInit.

    void inic(double val);
    void set(int i, int j, int k, double val);
    void copyAll(mat3 &B);
    void setAll(double val);

    void add(int i, int j, int k, double val);
    double at(int i, int j, int k);
    double cols();
    double rows();
    double dims();
    void print();


    private:
    double * data;
    int i_max_;
    int j_max_;
    int k_max_;


};

    mat3::mat3() {
        i_max_ = 0;
        j_max_ = 0;
        k_max_ = 0;
        data = NULL;
    }

    mat3::mat3(int i_max_, int j_max_, int k_max_) {
        this->i_max_ = i_max_;
        this->j_max_ = j_max_;
        this->k_max_ = k_max_;
        inic(0);
    }


    mat3::mat3(int i_max_, int j_max_, int k_max_, double val) {
        this->i_max_ = i_max_;
        this->j_max_ = j_max_;
        this->k_max_ = k_max_;
        inic(val);
    }

    void mat3::inic(double val) {
        int size = i_max_ * j_max_ * k_max_;
        data = new double[size];

        for (int k = 0; k < size; k++)
            data[k] = val;
    }


    void mat3::setAll(double val) {
        int size = i_max_ * j_max_ * k_max_;
        for (int k = 0; k < size; k++)
            data[k] = val;
    }



    void mat3::set(int i, int j, int k, double val) {
      if(std::isnan(val)){
      cerr << "Error: Setting Nan value. Terminating." << endl;
      std::exit(1);
      }
      data[k_max_ * j_max_ * i + j*k_max_ + k] = val;
    }

    void mat3::copyAll(mat3 &B) {
      double val;
        for (int i = 0; i < rows(); ++i) {
            for (int j = 0; j < cols(); ++j) {
                for (int k = 0; k < dims(); ++k) {
                    set(i, j, k, B.at(i, j, k));
                }
            }
        }
    }

    void mat3::add(int i, int j, int k, double val) {
        data[i_max_ * j_max_ * k + i * j_max_ + j] += val;
    }


    double mat3::at(int i, int j, int k) {
        return data[k_max_ * j_max_ * i + j*k_max_ + k];
    }


    double mat3::cols() {
        return j_max_;
    }


    double mat3::rows() {
        return i_max_;
    }


    double mat3::dims() {
        return k_max_;
    }


    void mat3::print() {
        for(int k = 0; k < k_max_; k++) {
            for (int i = 0; i < i_max_; i++) {
                for (int j = 0; j < j_max_; j++) {
                    cout << data[i_max_ * j_max_ * k + i * j_max_ + j] << " ";
                }
                cout << endl;
            }
            cout << endl;
            cout << endl;
        }
        cout << endl;
    }

    mat3::~mat3() {
        if (data != NULL)
            delete [] data;
    }

    double diff(mat3 &A, mat3 &B) {
        double sum = 0;
        for (int i = 0; i < A.rows(); ++i) {
            for (int j = 0; j < A.cols(); ++j) {
                        for (int k = 0; k < A.cols(); ++k) {

                sum += pow(B.at(i, j, k) - A.at(i, j, k), 2);
            }
        }
    }
        return sqrt(sum);
    }

    void mat3::confAndInit(int i_max_, int j_max_, int k_max_, double val){
      this->i_max_ = i_max_;
      this->j_max_ = j_max_;
      this->k_max_ = k_max_;
      inic(val);
    }

#endif
