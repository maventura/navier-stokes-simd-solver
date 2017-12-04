#ifndef MAT3_H
#define MAT3_H

#include <iostream>


using namespace std;

class mat2 {

  public:
    float *data;

    mat2();
    ~mat2();
    mat2(int i_max_, int j_max_);
    mat2(int i_max_, int j_max_, float val);
    void confAndInit(int i_max_, int j_max_, float val);
    //TODO: No control over multiple allocations using confAndInit.

    void inic(float val);
    void set(int i, int j, float val);
    void copyAll(mat2 &B);
    void setAll(float val);

    void add(int i, int j, float val);
    float at(int i, int j);
    float cols();
    float rows();
    void print();
    bool checkForNan();


  private:
    int i_max_;
    int j_max_;

};



bool mat2::checkForNan(){
    /*for (int i = 0; i < rows(); ++i) {
        for (int j = 0; j < cols(); ++j) {
            if (std::isnan(at(i,j))) {
                cerr << "NaN fount at (" << i << ", " << j <<")" << endl; 
                return true;
            }
        }
    }*/
    return false;
}

mat2::mat2() {
    i_max_ = 0;
    j_max_ = 0;
    data = NULL;
}

mat2::mat2(int i_max_, int j_max_) {
    this->i_max_ = i_max_;
    this->j_max_ = j_max_;
    inic(0);
}


mat2::mat2(int i_max_, int j_max_, float val) {
    this->i_max_ = i_max_;
    this->j_max_ = j_max_;
    inic(val);
}

void mat2::inic(float val) {
    int size = i_max_ * j_max_;
    data = new float[size];

    for (int idx = 0; idx < size; idx++)
        data[idx] = val;
}


void mat2::setAll(float val) {
    int size = i_max_ * j_max_;
    for (int idx = 0; idx < size; idx++)
        data[idx] = val;
}



void mat2::set(int i, int j, float val) {
    /*if (std::isnan(val)) {
        cerr << "Error: Setting Nan value at (" << i << "," << j << "). Terminating." << endl;
        std::exit(1);
    }*/
    data[j_max_ * i + j] = val;
}

void mat2::copyAll(mat2 &B) {
    float val;
    for (int i = 0; i < rows(); ++i) {
        for (int j = 0; j < cols(); ++j) {
            set(i, j, B.at(i, j));
        }
    }
}

void mat2::add(int i, int j, float val) {
    data[j_max_ * i + j ] += val;
}


float mat2::at(int i, int j) {
    return data[j_max_ * i + j];
}


float mat2::cols() {
    return j_max_;
}


float mat2::rows() {
    return i_max_;
}


void mat2::print() {
    for (int i = 0; i < i_max_; i++) {
        for (int j = 0; j < j_max_; j++) {
            float num = data[j_max_ * i + j];
            cout << num << " ";
        }
        cout << endl;
        cout << endl;
    }
    cout << endl;
}

mat2::~mat2() {
    if (data != NULL)
        delete [] data;
}

float diff(mat2 &A, mat2 &B) {
    float sum = 0;
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            sum += pow(B.at(i, j) - A.at(i, j), 2);
        }
    }
    return sqrt(sum);
}

void mat2::confAndInit(int i_max_, int j_max_, float val) {
    this->i_max_ = i_max_;
    this->j_max_ = j_max_;
    inic(val);
}

#endif
