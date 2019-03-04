#include "Matrix.h"
#include "MatrixOpe.h"
#include <stdlib.h>

template<typename T> 
Matrix<T> operator+(Matrix<T> &m1, Matrix<T> &m2) {
    return m1.addMatrix(m1, m2);
};