#include "Matrix.h"
#include <iostream>

using namespace std;

template <typename T>
Matrix<T>::Matrix() {
};

template <typename T>
Matrix<T>::Matrix(unsigned int row, unsigned int col, T lamd, const string &type){
    initMatrix(row, col, lamd, type);
}

template <typename T>
void Matrix<T>::initMatrix(unsigned int row, unsigned int col, T lamd, const string &type){
    if(row == 0|| col == 0) {
        cout << "Not a valid matrix" << endl;
        return;
    }

    // initialize 
    RowData<T> colData(col);
    Data<T> da(row, colData);
    unsigned int i = 0,j = 0;

    if(type!="diag") {
        colData = RowData<T>(col, lamd);
        da = Data<T>(row, colData);
    } else {
        unsigned int n = min(row, col);
        for(int i=0; i < n; i++) da[i][i] = lamd;
    }

    // assign values to objects
    this->data = da;
    this->row = row;
    this->col = col;
}

template <typename T>
void Matrix<T>::LoadData(const char* filename) {
    loadData<T>(this->data,filename);
    this->row=this->data.size();
    this->col=this->data[0].size();
}

template <typename T>
void Matrix<T>::print() {
    unsigned int i,j;
    for(i=0; i<row; i++) {
        for(j=0; j<col; j++)
            cout<<data[i][j]<<"  ";
        cout<<endl;
    }
}

template <typename T>
Matrix<T> Matrix<T>::copyMatrix(Matrix<T> &m)
{
    // replicated matrix to be stored
    Matrix<T> cp;
    unsigned int i=0,j=0;
    RowData<T> cda(this->col);
    Data<T> da(this->row,cda);

    for(i=0; i < this->row; i++){
        for(j=0; j < this->col; j++) 
            da[i][j] = this->data[i][j];
    }

    cp.row = this->row;
    cp.col = this->col;
    cp.data = da;
    return cp;
}

template <typename T>
Matrix<T> Matrix<T>::getRow(unsigned int iRow)
{
    if(iRow >= this->row) {
        throw ("row index out of range");
    }
    // return an row vector that is represented in 1*col
    RowData<T> theRow = this->data[iRow];
    Matrix<T> matrix;
    matrix->data = Data<T>(1, this->data[iRow]);
    matrix->row = 1;
    matrix->col = this->col;
    return matrix;
}

template <typename T>
Matrix<T> Matrix<T>::getCol(unsigned int jCol) {
    if(jCol >= this->col) {
        throw("row index out of range");
    }
    // return a column taht is represented in row * 1
    RowData<T> oneColRow(1);
    Matrix<T> matrix;
    matrix->data = Data<T>(this->row, oneColRow);
    matrix.col=1;
    matrix.row=row;
    for(int i = 0; i < this->data.size(); i++) {
        //cout<<i<<"="<<this->data[i][jCol]<<endl;
        matrix.data[i][0]=this->data[i][jCol];
    }
    return matrix;
}

template <typename T>
void Matrix<T>::deleteRow(unsigned int iRow){
    unsigned int i=0;
    for(typename Data<T>::iterator it = data.begin(); it <= data.end(); it++,i++) {
        if(i==iRow) this->data.erase(it);
    }
    this->row--;
}

template <typename T>
void Matrix<T>::deleteCol(unsigned int iCol) {

    unsigned int i=0,j=0;
    Matrix cp=this->copyMatrix();

    for(typename Data<T>::iterator it = data.begin(); it < data.end(); it++, i++){
        j = 0;
        for(typename RowData<T>::iterator jt = *it.begin(); jt < *it.end(); jt++, j++)
            if(j == iCol) this->data[i].erase(jt);
    }
    this->col--;
}

template<typename T>
Matrix<T> Matrix<T>::transpose()//矩阵形式的转置
{
    unsigned int i=0, j=0;
    Matrix<T> matrixT;
    ColData<T> coldata(row);
    Data<T> da(col,coldata);
    matrixT->data=this->da;
    matrixT->row=this->col;
    matrixT->col=this->row;
    for(i=0; i<col; i++)
        for(j=0; j<row; j++)
            matrixT.data[i][j]=data[j][i];
    return matrixT;
}

template <typename T>
Matrix<T> Matrix<T>::addMatrix(const Matrix<T> &m1,const Matrix<T> &m2) {
    if(m1.col!=m2.col || m1.row!=m2.row){
        cout<<"shape of two matrices unmatch"<<endl;
        exit(-1);
    }

    // inialize new matrix
    RowData<T> rowdata(m1.col);
    Data<T> da(m1.row, rowdata);
    Matrix<T> addition;

    addition.data = da;
    addition.row = row;
    addition.col = col;

    // add up each element
    unsigned int i,j;
    for(i=0; i<m1.row; i++)
        for(j=0; j<m1.col; j++)
            addition.data[i][j]=m1.data[i][j] + m2.data[i][j];
    return addition;
}

template <typename T>
Matrix<T> Matrix<T>::subMatrix(const Matrix<T> &m1,const Matrix<T> &m2) {   
    if(m1.col != m2.col || m1.row != m2.row) {
        cout << "shape of two matrices mismatch" << endl;
        exit(-1);
    }

    // initialize new matrix
    RowData<T> rowdata(m1.col);
    Data<T> da(m1.row, rowdata);
    Matrix<T> subtraction;
    subtraction.data=da;
    subtraction.row=m1.row;
    subtraction.col=m1.col;

    // subtract value in m2 from m1, element-wise
    unsigned int i,j;
    for(i=0; i<m1.row; i++) 
        for(j=0; j < m1.col; j++)
            subtraction.data[i][j]=m1.data[i][j]-m2.data[i][j];
    return subtraction;
}

template <typename T>
Matrix<T> Matrix<T>::multsMatrix(const Matrix<T> &m1, const Matrix<T> &m2) {
    if(m1.col!=m2.row) {
        cout<<"multsData error"<<endl;
        exit(-1);
    }
    unsigned int i,j,k;
    Matrix<T> mults;
    ColData<T> cda(m2.col, 0);
    Data<T> da(m1.row, cda);
    mults.data = da;
    mults.row=m1.row;
    mults.col=m2.col;
    for(i=0; i<m1.row; i++)
        for(j=0; j<m2.col; j++)
            mults.data[i][j]=0;
    
    for(i = 0; i < m1.row; i++)
        for(j = 0; j < m2.col; j++)
            for(k = 0; k < m1.col; k++)
                mults.data[i][j]+=m1.data[i][k]*m2.data[k][j];
    return mults;
}

template <typename T>
Matrix<T> Matrix<T>::dotMultMatrix(const Matrix<T> &m1, const Matrix<T> &m2)//矩阵形式的相乘
{
    if(m1->row!=m2->row || m1->col!=m2->col) {
        cout<<"shape not match"<<endl;
        exit(-1);
    }
    unsigned int i,j;
    Matrix<T> dotmults;
    ColData<T> cda(m1->col);
    Data<T> da(m1->row, cda);
    for(i=0; i < m1->row; i++)
        for(j = 0; j < m2->col; j++)
            da[i][j]=m1->data[i][j] * m2->data[i][j];
    dotmults->data = da;
    return dotmults;
}

template <typename T>
Matrix<float> Matrix<T>::multsScalar(const Matrix<T> &m2, double alpha){

    unsigned int i,j,k;
    Matrix<T> mults;
    ColData<T> cda(m2.col, 0);
    Data<T> da(m2.row, cda);
    mults.data = da;
    mults.row=m2.row;
    mults.col=m2.col;

    for(i = 0; i < m2.row; i++)
        for(j = 0; j < m2.col; j++) 
            mults.data[i][j] = this->data[i][j] * alpha;
    return mults;
}
/*
//行列式
double Matrix::detMatrix()
{
    if(row!=col)
    {
        cout<<"Data det is no"<<endl;
        exit(-1);
    }
    Matrix mCopy=*this;
    double det=1;
    unsigned int i=0,j=0,k=0;
    double max=-9999999;
    int swap=-1;
    double temp;
    ColData cda(col);
    Data aij(row,cda);
    for(k=0; k<mCopy.col-1; k++)//k表示第k次消元，一共需要n-1次
    {
        for(i=0; i<mCopy.row; i++)
        {
            if(mCopy.data[i][k]>max)//每一次消元都是比较第k列的元素，选出第k列中最大的一行
            {
                swap=i;
            }
        }//找到第k次列主元消去的最大行的下标
        if(swap==-1||mCopy.data[swap][k]==0)
            return -1;//最大主元为0
        for(j=0; j<mCopy.col; j++)
        {
            temp=mCopy.data[k][j];
            mCopy.data[k][j]=mCopy.data[swap][j];
            mCopy.data[swap][j]=temp;
        }//第k次消元，选出最大的一行是swap行，与第k行交换
        for(i=k+1; i<mCopy.row; i++)
        {
            aij[i][k]=mCopy.data[i][k]/mCopy.data[k][k];// 第k次消元，主元素为第k行第k列，把第k行以下的行都进行消元
            for(j=k; j<mCopy.col; j++)//对于k行以下的每一行的每一列元素都减去主行与消元因子的乘积
            {
                mCopy.data[i][j]-=aij[i][k]*mCopy.data[k][j];
            }
        }
    }
    for(i=0; i<mCopy.row; i++)
    {
        det*=mCopy.data[i][i];
    }
    //cout<<"det="<<det<<endl;
    return det;
}
//高斯消元矩阵求逆,特别注意，LU分解不能进行行列式变换
/*Matrix Matrix::niMatrix()
{
    if(row!=col)
    {
        cout<<"Data ni is no "<<endl;
        exit(-1);
    }
    if(this->detMatrix()==0)//这里调用求行列式进行了列主元消去改变了参数矩阵，如何传递不改变是一个问题
    {
        cout<<"Data det is no so ni is no "<<endl;
        exit(-1);
    }
    unsigned int i=0,j=0,k=0;//这里存在-1的情况，务必定义为int型
    double temp;
    Matrix mCopy=*this;
    Matrix UMatrix=*this;
    Matrix LMatrix=*this;
    Matrix UniMatrix=*this;
    Matrix LniMatrix=*this;
    ColData cda(col);
    Data aij(row,cda);
    for(k=0; k<col-1; k++)//k表示第k次消元，一共需要n-1次
    {
        for(i=k+1; i<row; i++)
        {
            aij[i][k]=data[i][k]/data[k][k];// 第k次消元，主元素为第k行第k列，把第k行以下的行都进行消元
            for(j=k; j<col; j++)//对于k行以下的每一行的每一列元素都减去主行与消元因子的乘积
            {
                data[i][j]-=aij[i][k]*data[k][j];
            }
        }
    }
    UMatrix=*this;
    for(j=0; j<col; j++)
    {
        for(i=j+1; i<row; i++)
        {
            temp=0;
            for(k=0; k<j; k++)
            {
                temp=LMatrix.data[i][k]*UMatrix.data[k][j];
            }
            LMatrix.data[i][j]=1/UMatrix.data[j][j]*(mCopy.data[i][j]-temp);
        }
    }
    for(i=0; i<row; i++)
    {
        for(j=0; j<col; j++)
        {
            if(i==j)
                LMatrix.data[i][j]=1;
            if(j>i)
                LMatrix.data[i][j]=0;
        }
    }
    Matrix mults;
    mults=*this;
    mults=mults.multsMatrix(LMatrix,UMatrix);
    Matrix LU=mults;
    //cout<<"lu"<<endl;
    //mults.print();

    //计算u逆
    for(j=0; j<col; j++)
    {
        for(i=j; (int)i>=0; i--)
        {
            if(i==j)
                UniMatrix.data[i][j]=1/UMatrix.data[i][j];
            else
            {
                temp=0;
                for(k=j; k>i; k--)
                {
                    temp+=UMatrix.data[i][k]*UniMatrix.data[k][j];
                }
                UniMatrix.data[i][j]=-1/UMatrix.data[i][i]*temp;
            }
        }
        ///关键，将下三角清零
        for(i=j+1; i<row; i++)
            UniMatrix.data[i][j]=0;
    }
    //计算l逆
    for(j=0; j<col; j++)
    {
        for(i=0; i<row; i++)
        {
            if(j==i)
                LniMatrix.data[i][j]=1;
            else
            {
                temp=0;
                for(k=j; k<i; k++)
                {
                    temp+=(LMatrix.data[i][k]*LniMatrix.data[k][j]);
                }
                LniMatrix.data[i][j]=-temp;
            }
        }
    }

    mults=mults.multsMatrix(UniMatrix,LniMatrix);
    *this=mCopy;
    Matrix I=*this;
    I=I.multsMatrix(LU,mults);
    //cout<<"LU"<<"*"<<"LUni"<<endl;
    //I.print();
    return mults;
}*/

//template <typename T>
//Matrix<T>::operator+(const Matrix<T> &m2){
//    return addMatrix(*this, m2);
//}