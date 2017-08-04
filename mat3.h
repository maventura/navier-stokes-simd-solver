#include <iostream>
using namespace std;

class mat3 {

public:
    mat3();
    ~mat3();
    mat3(int i_max_, int j_max_, int k_max_);
    mat3(int i_max_, int j_max_, int k_max_, long double val);
    void inic(long double val);
    void set(int i, int j, int k, long double val);
    void setAll(mat3 &B);
    void add(int i, int j, int k, long double val);
    long double at(int i, int j, int k);
    long double cols();
    long double rows();
    long double dims();
    void print();



    private:
    long double * data;
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


    mat3::mat3(int i_max_, int j_max_, int k_max_, long double val) {
        this->i_max_ = i_max_;
        this->j_max_ = j_max_;
        this->k_max_ = k_max_;
        inic(val);
    }

    void mat3::inic(long double val) {
        int n = i_max_ * j_max_ * k_max_;
        data = new long double[n];

        for (int k = 0; k < n; k++)
            data[k] = val;
    }

    void mat3::set(int i, int j, int k, long double val) {
        data[k_max_ * j_max_ * i + j*k_max_ + k] = val;
    }

    void mat3::setAll(mat3 &B) {
        for (int i = 0; i < rows(); ++i) {
            for (int j = 0; j < cols(); ++j) {
                for (int k = 0; k < dims(); ++k) {
                    set(i, j, k, B.at(i, j, k));
                }
            }
        }
    }


    void mat3::add(int i, int j, int k, long double val) {
        data[i_max_ * j_max_ * k + i * j_max_ + j] += val;
    }


    long double mat3::at(int i, int j, int k) {
        return data[k_max_ * j_max_ * i + j*k_max_ + k];
    }


    long double mat3::cols() {
        return j_max_;
    }


    long double mat3::rows() {
        return i_max_;
    }


    long double mat3::dims() {
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
