#include <iostream>
#include "Num.h"
#include <typeinfo>
#include <vector>
#include "Matrix.h"
// #include "MatrixOpe.h"
#include <string>

using namespace std;

// Matrix<int> operator+(Matrix<int> m1, Matrix<int> m2) {
//     return m1.addMatrix(m1, m2);  
// };

//template<typename T>
//Matrix<T> operator +(const Matrix<T>&m1, const Matrix<T> &m2) {
//    return m1.addMatrix(m2);
//}

int main() {
    // Num<int> num(35);
    // cout << num.getNum() << endl;
    vector<int> a{1,2,3,4};
    for(vector<int>::iterator i = a.begin(); i <= a.end(); i++) {
        *i += 1;
    }
    vector<vector<int>> m = {{1,2}, {1,2}};
    int i = 0, j = 0;

    for(i = 0; i < m.size(); i++){
        j = 0;
        for(vector<int>::iterator it = m[i].begin(); it < m[i].end(); it++, j++)
            if(j == 1) m[i].erase(it);
    }

    for(auto &i : m){
        for(auto &j : i) cout << j << " ";
        cout << "\n";
    }
    const string name("not");
    // test Matrix
    unsigned int c = 7000, b = 7000;
    float n = 1.1;  
    Matrix<float> m1(c, b, n, name);
    
    cout << "size of m1 : " << " row, col : " << 
            m1.row << " , " << m1.col << endl;

    Matrix<float> m2 = m1.copyMatrix(m1);
    Matrix<float> m3 = m1.multsMatrix(m1, m2);
    Matrix<float> m4 = m1 / 3.1;

    cout << "size of m3 : " << " row, col : " << 
            m3.row << " , " << m3.col << endl;

    cout << "size of m4 : " << " row, col : " << 
            m4.row << " , " << m4.col << endl;
    return 0;
}