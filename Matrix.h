#ifndef Matrix_H
#define Matrix_H

#include "LoadData.h"
using namespace std;

template<typename T>
class Matrix {
public:
    Data<T> data;
    unsigned int row;
    unsigned int col;
public:
    Matrix();

    Matrix(unsigned int row, unsigned int col, T lamd, const string &type); // overloading constructor

    void initMatrix(unsigned int row, unsigned int col,T lamd, const string &type);

    void LoadData(const char* filename);

    void print();

    Matrix<T> copyMatrix(Matrix<T> &m);

    Matrix<T> getRow(unsigned int iRow);

    Matrix<T> getCol(unsigned int jCol);

    void deleteRow(unsigned int iRow);

    void deleteCol(unsigned int jCol);

    Matrix<T> transpose();

    Matrix<T> addMatrix(const Matrix<T> &m1, const Matrix<T> &m2);  // Matrix addition

    Matrix<T> subMatrix(const Matrix<T> &m1, const Matrix<T> &m2); // Matrix subtraction

    Matrix<T> multsMatrix(const Matrix<T> &m1, const Matrix<T> &m2); // Matrix multiplication

    Matrix<T> dotMultMatrix(const Matrix<T> &m1, const Matrix<T> &m2); // dot Matrix multiplication

    Matrix<float> multsScalar(const Matrix<T> &m2, double alpha);

    // double detMatrix();
    // template <typename T>
    Matrix<T> operator+(const Matrix<T> &m2){
        return this->addMatrix(*this, m2);
    };

    Matrix<T> operator*(const Matrix<T> &m2){
        return this->multsMatrix(*this, m2);
    }

    Matrix<float> operator*(double alpha){
        return this->multsScalar(*this, alpha);
    }

    Matrix<T> operator-(const Matrix<T> &m2){
        return this->subMatrix(*this, m2);
    }

    Matrix<T> operator+=(const Matrix<T> &m2) {
        Matrix<T>tmp =  this->addMatrix(*this, m2);
        this->data = tmp.data; 
    }
    
    Matrix<float> operator/(double alpha){
        return this->multsScalar(*this, 1/alpha);
    }


    // Matrix<T>
};

#include "Matrix.cpp"

#endif